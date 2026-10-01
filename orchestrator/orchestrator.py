"""
orchestrator/orchestrator.py — SEA Orchestrator
================================================
Single entry-point module that accepts a viral + host protein sequence pair,
routes them through the three existing SEA modules, and returns a unified,
ranked output describing autoimmune risk via the Super-epitope Architecture
(SEA) mechanism.

Pipeline
--------
1.  Profile viral sequence  → find_all_degradation_motifs()
2.  Profile host sequence   → find_all_degradation_motifs()
3.  Extract host KFERQ motifs → provided_motifs (for SEAModule)
4.  Check host CMA membership → check_cma_network_membership()
5a. Source A: find_homologous_pairs()     → seq_aligner pairs
5b. Source B: mclachlan_to_pairs()        → mclachlan pairs   [if mclachlan_hits supplied]
5c. Source C: mclachlan_to_pairs(sw_hits) → sw pairs          [if sw_hits supplied]
6.  Merge + _dedup_pairs(tolerance=4)    → merged pairs, re-ranked
7.  SEAModule.run(merged_pairs, ...)      → sea_results
7b. apply_epitope_proximity_bonus()       → bonus scores for hits in known epitopes
8.  Optional calculate_cma_score()        → cma_score
9.  Build OrchestratorResult + risk_summary

Usage
-----
>>> from orchestrator.orchestrator import orchestrate
>>> result = orchestrate(viral_seq, host_seq, "HCV", "polyprotein", "CYP2E1",
...                      host_accession="P05181")
>>> print(result.risk_summary)
>>> result.report("/home/sandbox/reports/hcv_cyp2e1_sea.md")

The Orchestrator is a thin integration layer.
It does NOT reimplement hinge/jammer/KFERQ logic — it calls the existing modules.
"""

from __future__ import annotations

import os
import sys
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

# ── sys.path: resolve sibling modules under /home/sandbox/ ───────────────────
_REPO_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Import modules (fail fast with a clear message) ──────────────────────────
try:
    from motif_finder import find_all_degradation_motifs, find_kferq_motifs
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        f"motif_finder not found on sys.path={sys.path!r}\n"
        f"Original error: {exc}"
    )

try:
    from cma_network import check_cma_network_membership, calculate_cma_score
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        f"cma_network not found on sys.path={sys.path!r}\n"
        f"Original error: {exc}"
    )

try:
    from sea.sea_module import (
        SEAModule, SEAConfig, find_homologous_pairs,
        HingeScanner, JammerScanner, HingeTier,
    )
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        f"sea.sea_module not found on sys.path={sys.path!r}\n"
        f"Original error: {exc}"
    )


# ════════════════════════════════════════════════════════════════════════════
#  RESULT DATACLASS
# ════════════════════════════════════════════════════════════════════════════

@dataclass
class OrchestratorResult:
    """Unified output from the Orchestrator pipeline."""

    # ── Identity ─────────────────────────────────────────────────────────────
    virus_name:         str
    viral_protein_name: str
    host_protein_name:  str

    # ── Core outputs ─────────────────────────────────────────────────────────
    sea_results:          List[Any]          # List[SEAResult], typed as Any to avoid circular
    viral_motif_profile:  Dict
    host_motif_profile:   Dict
    host_cma_membership:  Optional[Dict]
    cma_score:            Optional[float]

    # ── Summary ──────────────────────────────────────────────────────────────
    risk_summary:        Dict               = field(default_factory=dict)
    pair_source_counts:  Dict[str, int]     = field(default_factory=dict)

    # ── Epitope proximity bonus ───────────────────────────────────────────────
    host_accession:                  Optional[str]  = None
    epitope_proximity_hits:          List[Dict]     = field(default_factory=list)
    epitope_proximity_bonus_applied: float          = 0.0

    # ── Hinge / jammer full-sequence profiles ─────────────────────────────────
    viral_hinge_jammer_profile:  Dict  = field(default_factory=dict)
    host_hinge_jammer_profile:   Dict  = field(default_factory=dict)

    # ── Phase 2 — Degradation pathway convergence ─────────────────────────────
    phase2_enriched_pairs: List[Dict]  = field(default_factory=list)
    calibration_result:    Dict        = field(default_factory=dict)

    # ── Phase 3 — TRUE super-epitope pairs (viral bridged pairs) ──────────────
    super_epitope_pairs:   List[Dict]  = field(default_factory=list)

    # ── Meta ─────────────────────────────────────────────────────────────────
    timestamp:           str                = field(
        default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    )

    # ── Report ───────────────────────────────────────────────────────────────
    def report(self, output_path: str) -> str:
        """
        Write a structured Markdown report to *output_path*.

        Returns the absolute path of the written file.
        The report uses consistent H1/H2/H3 headings so it drops cleanly
        into an Obsidian vault.
        """
        lines = _build_report_lines(self)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        return os.path.abspath(output_path)


# ════════════════════════════════════════════════════════════════════════════
#  KNOWLEDGE BASE
# ════════════════════════════════════════════════════════════════════════════

_DEFAULT_KB_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "protein_knowledge_base.json",
)


def load_protein_knowledge_base(path: Optional[str] = None) -> Dict:
    """
    Load ``protein_knowledge_base.json`` from *path*.

    Defaults to ``orchestrator/protein_knowledge_base.json`` (the file
    co-located with this module).

    A missing or malformed file is **non-fatal**: the function returns ``{}``
    so the pipeline continues without any epitope proximity bonuses.

    Parameters
    ----------
    path : explicit path to the JSON file, or None to use the default.

    Returns
    -------
    Parsed dict from the JSON file, or {} on any failure.
    """
    resolved = path if path is not None else _DEFAULT_KB_PATH
    if not os.path.exists(resolved):
        return {}
    try:
        with open(resolved, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError):
        return {}


def get_epitope_ranges(
    kb_data:   Dict,
    accession: str,
) -> List[Tuple[int, int, str]]:
    """
    Extract mapped autoantibody epitope ranges for *accession* from *kb_data*.

    Only entries with non-null ``residue_start`` **and** ``residue_end`` are
    included (conformational epitopes with null positions are skipped).

    Parameters
    ----------
    kb_data   : loaded knowledge base dict (from load_protein_knowledge_base)
    accession : UniProt primary accession (e.g. "P05181")

    Returns
    -------
    List of ``(residue_start, residue_end, label)`` tuples, 1-based inclusive.
    Empty list if the accession is absent or has no mapped linear epitopes.
    """
    entry    = kb_data.get("proteins", {}).get(accession, {})
    epitopes = entry.get("known_autoantibody_epitopes", [])
    ranges: List[Tuple[int, int, str]] = []
    for ep in epitopes:
        start = ep.get("residue_start")
        end   = ep.get("residue_end")
        label = ep.get("label", "unknown")
        if start is not None and end is not None:
            ranges.append((int(start), int(end), label))
    return ranges


def apply_epitope_proximity_bonus(
    sea_results:    List[Any],   # List[SEAResult]
    host_accession: str,
    kb_data:        Dict,
    bonus:          float,
) -> Tuple[List[Any], List[Dict]]:
    """
    Add *bonus* to the ``final_sea_score`` of every SEAResult whose host
    anchor lands inside a known autoantibody epitope for *host_accession*.

    "Host anchor" is ``result.position2`` (0-based).  It is converted to a
    1-based residue position and checked against each ``(start, end)`` range
    from the knowledge base.  The bonus is applied **at most once per
    result** even if the anchor overlaps multiple epitope ranges (first
    matching range wins).  Each boosted result also gets an explanatory
    entry appended to its ``notes`` list.

    After all bonuses are applied, *sea_results* is re-sorted by
    ``final_sea_score`` descending (in place).

    Parameters
    ----------
    sea_results    : list of SEAResult objects — **mutated in place**
    host_accession : UniProt accession for the host protein
    kb_data        : knowledge base dict (from load_protein_knowledge_base)
    bonus          : score delta (e.g. 2.0).  Pass 0.0 to disable.

    Returns
    -------
    (modified_sea_results, proximity_hits)
        modified_sea_results : input list, re-sorted by final_sea_score DESC
        proximity_hits       : list of dicts — one per boosted result —
                               with keys: rank, seq1, seq2, host_pos_1based,
                               epitope_label, epitope_range, bonus_applied,
                               new_score
    """
    if bonus == 0.0 or not sea_results:
        return sea_results, []

    epitope_ranges = get_epitope_ranges(kb_data, host_accession)
    if not epitope_ranges:
        return sea_results, []

    proximity_hits: List[Dict] = []

    for result in sea_results:
        host_pos_1based = result.position2 + 1   # 0-based → 1-based
        for (ep_start, ep_end, ep_label) in epitope_ranges:
            if ep_start <= host_pos_1based <= ep_end:
                result.final_sea_score = round(result.final_sea_score + bonus, 4)
                result.notes.append(
                    f"Epitope proximity bonus +{bonus:.2f}: host pos {host_pos_1based} "
                    f"in known epitope '{ep_label}' ({ep_start}–{ep_end})"
                )
                proximity_hits.append({
                    "rank":            result.rank,
                    "seq1":            result.seq1,
                    "seq2":            result.seq2,
                    "host_pos_1based": host_pos_1based,
                    "epitope_label":   ep_label,
                    "epitope_range":   f"{ep_start}–{ep_end}",
                    "bonus_applied":   bonus,
                    "new_score":       result.final_sea_score,
                })
                break   # at most one bonus per result

    sea_results.sort(key=lambda r: r.final_sea_score, reverse=True)
    return sea_results, proximity_hits


def annotate_epitope_overlap(
    sea_results:    List[Any],   # List[SEAResult]
    host_accession: str,
    kb_data:        Dict,
) -> List[Dict]:
    """
    Annotate each SEAResult whose host anchor overlaps a known autoantibody
    epitope for *host_accession*.

    Does **not** modify ``final_sea_score``.  Appends a note to
    ``result.notes`` for traceability only.

    Parameters
    ----------
    sea_results    : list of SEAResult objects — notes are mutated in place
    host_accession : UniProt accession for the host protein
    kb_data        : knowledge base dict (from load_protein_knowledge_base)

    Returns
    -------
    List[Dict]
        One entry per matching result, with keys:
        seq1, seq2, position1, position2, host_pos_1based,
        epitope_label, epitope_range.
    """
    if not sea_results:
        return []

    epitope_ranges = get_epitope_ranges(kb_data, host_accession)
    if not epitope_ranges:
        return []

    overlap_hits: List[Dict] = []
    for result in sea_results:
        host_pos_1based = result.position2 + 1   # 0-based → 1-based
        for (ep_start, ep_end, ep_label) in epitope_ranges:
            if ep_start <= host_pos_1based <= ep_end:
                result.notes.append(
                    f"Known epitope overlap: host pos {host_pos_1based} "
                    f"in '{ep_label}' ({ep_start}–{ep_end}) [annotation only]"
                )
                overlap_hits.append({
                    "seq1":            result.seq1,
                    "seq2":            result.seq2,
                    "position1":       result.position1,
                    "position2":       result.position2,
                    "host_pos_1based": host_pos_1based,
                    "epitope_label":   ep_label,
                    "epitope_range":   f"{ep_start}–{ep_end}",
                })
                break   # at most one annotation per result
    return overlap_hits


# ════════════════════════════════════════════════════════════════════════════
#  McLACHLAN CONVERSION
# ════════════════════════════════════════════════════════════════════════════

def _scan_hinge_jammer_profile(seq: str, cfg: Optional[SEAConfig] = None) -> Dict:
    """
    Full-sequence hinge and jammer scan.

    Slides a 7-residue target window across *seq* with step=7, scanning
    the flanking PROXIMITY_WINDOW on each side.  Unique matches are
    de-duplicated by (position, sequence) so each site is counted once.

    Returns a dict with keys:
        hinges_t1  : List[HingeMatch]  — Tier 1 (phosphoserine-type)
        hinges_t2  : List[HingeMatch]  — Tier 2 (Ser/Thr flexible)
        jammers    : List[JammerMatch]
        n_hinges_t1, n_hinges_t2, n_jammers  : int counts
    """
    if cfg is None:
        cfg = SEAConfig()

    hs  = HingeScanner()
    js  = JammerScanner()
    win = cfg.PROXIMITY_WINDOW
    step = 7

    seen_h: Dict[tuple, Any] = {}
    seen_j: Dict[tuple, Any] = {}

    for start in range(0, len(seq), step):
        end = start + step
        for h in hs.scan(seq, start, end, win, cfg):
            key = (h.position, h.sequence, h.tier)
            if key not in seen_h:
                seen_h[key] = h
        for j in js.scan(seq, start, end, win, cfg):
            key = (j.position, j.sequence)
            if key not in seen_j:
                seen_j[key] = j

    hinges_t1 = [h for h in seen_h.values() if h.tier == HingeTier.TIER1]
    hinges_t2 = [h for h in seen_h.values() if h.tier == HingeTier.TIER2]
    jammers   = list(seen_j.values())

    return {
        "hinges_t1":   hinges_t1,
        "hinges_t2":   hinges_t2,
        "jammers":     jammers,
        "n_hinges_t1": len(hinges_t1),
        "n_hinges_t2": len(hinges_t2),
        "n_jammers":   len(jammers),
    }


def mclachlan_to_pairs(
    mclachlan_hits:    List[Dict],
    viral_seq:         str,
    min_composite:     float = 15.0,
    max_pairs:         int   = 200,
    source:            str   = "mclachlan",
) -> List[Dict]:
    """
    Convert McLachlan v3 hit dicts (from the cross-protein comparator) to the
    pair dict schema expected by SEAModule.run().

    McLachlan positions are 1-based (cross_protein_comparator.py line 180).
    SEAModule expects 0-based positions — we subtract 1 from the start.

    Parameters
    ----------
    mclachlan_hits  : list of hit dicts from Q9WMX2_vs_P05181_scored_v3.json
                      ('hits' key) or any compatible McLachlan output
    viral_seq       : full viral protein sequence (stored in protein1 field)
    min_composite   : composite_primary threshold — hits below this are skipped
    max_pairs       : maximum number of pairs to return (sorted by composite_primary DESC)
    source          : value written into the 'source' field of every pair dict.
                      Use "mclachlan" for gapless McLachlan v3 hits (default)
                      and "sw" for gap-tolerant Smith-Waterman hits.

    Returns
    -------
    List of pair dicts with keys:
        seq1             str    viral fragment
        seq2             str    host fragment
        position1        int    0-based start in viral_seq
        position2        int    0-based start in host_seq
        protein1         str    full viral_seq
        similarity_score float  min(layer2_cross_mean / 4.0, 1.0)
        rank             int    0  (re-assigned by Orchestrator after dedup)
        source           str    value of the source parameter
    """
    # Filter
    filtered = [
        h for h in mclachlan_hits
        if h.get("composite_primary", 0.0) >= min_composite
    ]

    # Sort best first
    filtered.sort(key=lambda h: h.get("composite_primary", 0.0), reverse=True)

    pairs: List[Dict] = []
    for hit in filtered[:max_pairs]:
        try:
            p1_str = hit["p1_pos1"]   # e.g. "384-390" (1-based)
            p2_str = hit["p2_pos"]    # e.g. "48-54"  (1-based)

            p1_start = int(p1_str.split("-")[0]) - 1   # → 0-based
            p2_start = int(p2_str.split("-")[0]) - 1   # → 0-based

            raw_sim = hit.get("layer2_cross_mean", 0.0)
            similarity_score = min(raw_sim / 4.0, 1.0)

            pairs.append({
                "seq1":             hit["motif"],
                "seq2":             hit.get("p2_window", ""),
                "position1":        p1_start,
                "position2":        p2_start,
                "protein1":         viral_seq,
                "similarity_score": round(similarity_score, 4),
                "rank":             0,
                "source":           source,
                # carry forward for diagnostics
                "_composite_primary": hit.get("composite_primary", 0.0),
                "_layer2_cross_mean": raw_sim,
            })
        except (KeyError, ValueError, IndexError):
            # Malformed hit — skip silently
            continue

    return pairs


# ════════════════════════════════════════════════════════════════════════════
#  DEDUPLICATION
# ════════════════════════════════════════════════════════════════════════════

def _dedup_pairs(pairs: List[Dict], tolerance: int = 4) -> List[Dict]:
    """
    Remove duplicate pairs where both positions are within *tolerance* residues
    of an already-accepted pair.  Keep the higher similarity_score.

    Strategy: sort descending by similarity_score, then greedily accept pairs
    that are not within *tolerance* of any already-accepted position block.

    Parameters
    ----------
    pairs     : list of pair dicts (each must have position1, position2,
                similarity_score)
    tolerance : minimum position distance in BOTH coordinates to be kept

    Returns
    -------
    Deduplicated list, sorted by similarity_score descending.
    """
    sorted_pairs = sorted(pairs, key=lambda p: p.get("similarity_score", 0.0), reverse=True)
    accepted:  List[Dict]        = []
    seen_pos:  List[tuple]       = []   # (position1, position2) of accepted pairs

    for pair in sorted_pairs:
        p1 = pair.get("position1", 0)
        p2 = pair.get("position2", 0)

        is_dup = any(
            abs(p1 - ap1) < tolerance and abs(p2 - ap2) < tolerance
            for (ap1, ap2) in seen_pos
        )
        if not is_dup:
            accepted.append(pair)
            seen_pos.append((p1, p2))

    return accepted


# ════════════════════════════════════════════════════════════════════════════
#  RISK SUMMARY
# ════════════════════════════════════════════════════════════════════════════

def _build_risk_summary(
    sea_results:           List[Any],
    host_cma:              Optional[Dict],
    viral_profile:         Dict,
    host_profile:          Dict,
    cma_score_result:      Optional[Dict],
    viral_hj_profile:      Optional[Dict] = None,
    host_hj_profile:       Optional[Dict] = None,
) -> Dict:
    """Build a concise risk summary dict from Orchestrator outputs."""

    vhj = viral_hj_profile or {}
    hhj = host_hj_profile  or {}

    if not sea_results:
        return {
            "top_architecture":        "NONE",
            "top_final_score":         0.0,
            "n_complete_sea":          0,
            "n_sandwiched":            0,
            "n_pairs_scored":          0,
            "host_is_cma_member":      False,
            "host_cma_category":       None,
            "viral_n_kferq":           0,
            "host_n_kferq":            0,
            "viral_n_hinges_t1":       vhj.get("n_hinges_t1", 0),
            "viral_n_hinges_t2":       vhj.get("n_hinges_t2", 0),
            "viral_n_jammers":         vhj.get("n_jammers", 0),
            "host_n_hinges_t1":        hhj.get("n_hinges_t1", 0),
            "host_n_hinges_t2":        hhj.get("n_hinges_t2", 0),
            "host_n_jammers":          hhj.get("n_jammers", 0),
            "overall_risk":            "LOW",
            "n_epitope_proximity_hits": 0,
            "epitope_proximity_bonus": 0.0,
        }

    top     = sea_results[0]
    supers  = sum(1 for r in sea_results if r.is_complete_sea)
    sands   = sum(1 for r in sea_results if r.is_sandwiched and not r.is_complete_sea)
    top_arc = top.architecture_class.name if hasattr(top.architecture_class, "name") else str(top.architecture_class)

    host_cma_member   = host_cma is not None
    host_cma_category = host_cma.get("category") if host_cma else None

    v_kferq = viral_profile.get("summary", {}).get("n_kferq", 0)
    h_kferq = host_profile.get("summary",  {}).get("n_kferq", 0)

    # Simple ordinal risk tier
    score = top.final_sea_score
    if supers > 0 and score >= 15.0:
        overall_risk = "CRITICAL"
    elif supers > 0 or score >= 10.0:
        overall_risk = "HIGH"
    elif sands > 0 or score >= 5.0:
        overall_risk = "MODERATE"
    else:
        overall_risk = "LOW"

    # Elevate if host is a known CMA substrate
    if host_cma_member and host_cma_category in ("substrate", "regulator") and overall_risk in ("LOW", "MODERATE"):
        overall_risk = "HIGH"

    return {
        "top_architecture":        top_arc,
        "top_final_score":         round(score, 4),
        "n_complete_sea":          supers,
        "n_sandwiched":            sands,
        "n_pairs_scored":          len(sea_results),
        "host_is_cma_member":      host_cma_member,
        "host_cma_category":       host_cma_category,
        "viral_n_kferq":           v_kferq,
        "host_n_kferq":            h_kferq,
        "viral_n_hinges_t1":       vhj.get("n_hinges_t1", 0),
        "viral_n_hinges_t2":       vhj.get("n_hinges_t2", 0),
        "viral_n_jammers":         vhj.get("n_jammers", 0),
        "host_n_hinges_t1":        hhj.get("n_hinges_t1", 0),
        "host_n_hinges_t2":        hhj.get("n_hinges_t2", 0),
        "host_n_jammers":          hhj.get("n_jammers", 0),
        "overall_risk":            overall_risk,
        "cma_score":               (
            cma_score_result.get("cma_score") if cma_score_result else None
        ),
        # Populated in orchestrate() after apply_epitope_proximity_bonus()
        "n_epitope_proximity_hits": 0,
        "epitope_proximity_bonus":  0.0,
    }


# ════════════════════════════════════════════════════════════════════════════
#  REPORT BUILDER  (5-Phase Architecture)
# ════════════════════════════════════════════════════════════════════════════

def _arch_distribution(sea_results: List[Any]) -> Dict[str, int]:
    """Count SEA results by architecture class name."""
    dist: Dict[str, int] = {}
    for r in sea_results:
        name = r.architecture_class.name if hasattr(r.architecture_class, "name") else str(r.architecture_class)
        dist[name] = dist.get(name, 0) + 1
    return dict(sorted(dist.items(), key=lambda x: -x[1]))


def _build_report_lines(result: OrchestratorResult) -> List[str]:
    """Return a list of Markdown lines for the phased Orchestrator report."""
    rs  = result.risk_summary
    psc = result.pair_source_counts
    vhj = result.viral_hinge_jammer_profile  or {}
    hhj = result.host_hinge_jammer_profile   or {}

    lines: List[str] = []
    a = lines.append

    # ── Header ────────────────────────────────────────────────────────────
    a(f"# SEA Orchestrator Report: {result.virus_name} vs {result.host_protein_name}")
    a(f"")
    a(f"**Generated:** {result.timestamp}")
    host_id = f" (`{result.host_accession}`)" if result.host_accession else ""
    a(f"**Viral protein:** {result.viral_protein_name}  |  "
      f"**Host protein:** {result.host_protein_name}{host_id}")
    a(f"")
    a("---")
    a("")

    # ════════════════════════════════════════════════════════════
    # RISK SUMMARY  (top-of-report — synthesised from all phases)
    # ════════════════════════════════════════════════════════════
    a("## Risk Summary")
    a("")
    n_ep = rs.get("n_epitope_proximity_hits", 0)
    ep_b = rs.get("epitope_proximity_bonus",  0.0)
    ep_cell = f"+{ep_b:.2f} × {n_ep} hit(s)" if result.host_accession else "n/a"
    a("| Field | Value |")
    a("|-------|-------|")
    a(f"| Overall Risk | **{rs.get('overall_risk', 'N/A')}** |")
    a(f"| Top Architecture | {rs.get('top_architecture', 'N/A')} |")
    a(f"| Top Final Score | {rs.get('top_final_score', 0.0):.4f} |")
    a(f"| Complete-SEA Hits | {rs.get('n_complete_sea', 0)} |")
    a(f"| Sandwiched Hits | {rs.get('n_sandwiched', 0)} |")
    a(f"| Pairs Scored | {rs.get('n_pairs_scored', 0)} |")
    a(f"| Host CMA Member | {rs.get('host_is_cma_member', False)} ({rs.get('host_cma_category', 'N/A')}) |")
    a(f"| Viral KFERQ Motifs | {rs.get('viral_n_kferq', 0)} |")
    a(f"| Host KFERQ Motifs | {rs.get('host_n_kferq', 0)} |")
    a(f"| Viral Hinges (T1 / T2) | {vhj.get('n_hinges_t1', 0)} / {vhj.get('n_hinges_t2', 0)} |")
    a(f"| Viral Jammers | {vhj.get('n_jammers', 0)} |")
    a(f"| Host Hinges (T1 / T2) | {hhj.get('n_hinges_t1', 0)} / {hhj.get('n_hinges_t2', 0)} |")
    a(f"| Host Jammers | {hhj.get('n_jammers', 0)} |")
    if rs.get("cma_score") is not None:
        a(f"| CMA Score | {rs['cma_score']:.4f} |")
    a(f"| Epitope Proximity Bonus | {ep_cell} |")
    a("")
    a("---")
    a("")

    # ════════════════════════════════════════════════════════════
    # PHASE 1 — Amino Acid Homology
    # ════════════════════════════════════════════════════════════
    a("## Phase 1 — Amino Acid Homology")
    a("")
    a("*Source: sequence-aligner sliding-window + McLachlan physicochemical pair scoring.*")
    a("")
    a("### Pair Source Counts")
    a("")
    a("| Source | Count |")
    a("|--------|-------|")
    for k, v in psc.items():
        a(f"| {k} | {v} |")
    a("")
    a("### Top Homologous Pairs (pre-SEA, by composite_primary)")
    a("")
    # Collect pairs from sea_results; sort by base_score as proxy for composite_primary
    top_pairs = sorted(result.sea_results, key=lambda r: r.base_score, reverse=True)[:15]
    if top_pairs:
        a("| # | Viral Seq | Host Seq | Viral Pos | Host Pos | Base Score |")
        a("|---|-----------|----------|-----------|----------|-----------|")
        for i, r in enumerate(top_pairs, 1):
            a(f"| {i} | `{r.seq1}` | `{r.seq2}` | {r.position1} | {r.position2} | {r.base_score:.3f} |")
    else:
        a("*No pairs available.*")
    a("")
    a("---")
    a("")

    # ════════════════════════════════════════════════════════════
    # PHASE 2 — Degradation Motif Analysis
    # ════════════════════════════════════════════════════════════
    a("## Phase 2 — Degradation Motif Analysis")
    a("")
    a("*Source: protein_degradation module — KFERQ, LIR, D-box, hinge, jammer scans.*")
    a("")
    a("### Viral Protein Motif Profile")
    a("")
    vs = result.viral_motif_profile.get("summary", {})
    a(f"| Motif Type | Value |")
    a(f"|-----------|-------|")
    a(f"| CMA Category | `{vs.get('cma_category', 'N/A')}` |")
    a(f"| KFERQ Motifs | {vs.get('n_kferq', 0)} |")
    a(f"| LIR Motifs | {vs.get('n_lir', 0)} |")
    a(f"| D-box Degrons | {vs.get('n_dbox', 0)} |")
    a(f"| Hinges — Tier 1 | {vhj.get('n_hinges_t1', 0)} |")
    a(f"| Hinges — Tier 2 | {vhj.get('n_hinges_t2', 0)} |")
    a(f"| Jammers | {vhj.get('n_jammers', 0)} |")
    a(f"| Destabilising N-term | {vs.get('has_destabilising_n_term', False)} |")
    a(f"| Destabilising C-term | {vs.get('has_destabilising_c_term', False)} |")
    a("")
    a("### Host Protein Motif Profile")
    a("")
    hs = result.host_motif_profile.get("summary", {})
    a(f"| Motif Type | Value |")
    a(f"|-----------|-------|")
    a(f"| CMA Category | `{hs.get('cma_category', 'N/A')}` |")
    a(f"| KFERQ Motifs | {hs.get('n_kferq', 0)} |")
    a(f"| LIR Motifs | {hs.get('n_lir', 0)} |")
    a(f"| D-box Degrons | {hs.get('n_dbox', 0)} |")
    a(f"| Hinges — Tier 1 | {hhj.get('n_hinges_t1', 0)} |")
    a(f"| Hinges — Tier 2 | {hhj.get('n_hinges_t2', 0)} |")
    a(f"| Jammers | {hhj.get('n_jammers', 0)} |")
    a(f"| Destabilising N-term | {hs.get('has_destabilising_n_term', False)} |")
    a(f"| Destabilising C-term | {hs.get('has_destabilising_c_term', False)} |")
    a("")
    a("### Host CMA Network Membership")
    a("")
    if result.host_cma_membership:
        hc = result.host_cma_membership
        a(f"| Field | Value |")
        a(f"|-------|-------|")
        a(f"| Symbol | `{hc.get('symbol', 'N/A')}` |")
        a(f"| Category | {hc.get('category', 'N/A')} |")
        a(f"| Direction | {'+1 (CMA activating)' if hc.get('direction', 0) == 1 else '-1 (CMA inhibiting)'} |")
        a(f"| Weight | {hc.get('weight', 'N/A')} |")
    else:
        a(f"*{result.host_protein_name} is not a known CMA network member.*")
    a("")

    # ── Phase 2a — KFERQ Motif Proximity ────────────────────────────────
    a("### Phase 2a — Degradation Motif Proximity (Additive Event Scoring)")
    a("")
    a("*Distance-weighted proximity to LIR and KFERQ degradation motifs on both viral and host*")
    a("*sequences is computed as an additive event score (0–1) and multiplied by the*")
    a("*proximity_event_weight to boost final_sea_score.  No tier floors are applied;*")
    a("*tier classification is output-only, derived from total_event_signal.*")
    a("")

    p2 = result.phase2_enriched_pairs
    proximal_kferq_pairs = [p for p in p2 if p.get("proximal_kferq_motif")]
    proximal_lir_pairs   = [p for p in p2 if p.get("proximal_lir_motif")]
    proximal_both_pairs  = [p for p in p2 if p.get("proximal_kferq_motif") and p.get("proximal_lir_motif")]
    proximal_pairs = [p for p in p2 if p.get("proximal_deg_motif")]

    a(f"**Proximity summary:** KFERQ-proximal={len(proximal_kferq_pairs)}  |  "
      f"LIR-proximal={len(proximal_lir_pairs)}  |  "
      f"Both KFERQ+LIR={len(proximal_both_pairs)}")
    a("")
    if proximal_pairs:
        proximal_sorted = sorted(proximal_pairs, key=lambda x: x["final_sea_score"], reverse=True)
        a("| # | Viral Seq | V-Pos | Host Seq | H-Pos | V-KFdist | H-KFdist | V-LIRdist | H-LIRdist | Prox-Score | Prox-Boost | Tier | Score | KE |")
        a("|---|-----------|-------|----------|-------|----------|----------|-----------|-----------|------------|------------|------|-------|----|")
        for i, p in enumerate(proximal_sorted, 1):
            v_pos  = p["position1"] + 1
            h_pos  = p["position2"] + 1
            vkf    = str(p["viral_kferq_dist"]) if p["viral_kferq_dist"] is not None else "—"
            hkf    = str(p["host_kferq_dist"])  if p["host_kferq_dist"]  is not None else "—"
            vlir   = str(p["viral_lir_dist"])   if p.get("viral_lir_dist") is not None else "—"
            hlir   = str(p["host_lir_dist"])    if p.get("host_lir_dist")  is not None else "—"
            pscr   = f"{p.get('proximity_event_score', 0.0):.4f}"
            pboost = f"{p.get('proximity_boost', 0.0):.4f}"
            ke     = "★" if p.get("overlaps_known_epitope") else "—"
            a(f"| {i} | `{p['seq1']}` | {v_pos} | `{p['seq2']}` | {h_pos} "
              f"| {vkf} | {hkf} | {vlir} | {hlir} | {pscr} | {pboost} | {p['phase2_tier']} | {p['final_sea_score']:.4f} | {ke} |")
    else:
        a("*No homologous pairs found within the KFERQ or LIR proximity threshold.*")
        a("*Interpretation: degradation motif co-localisation is absent at this threshold;*")
        a("*Phase 2b convergence scoring remains the primary ranking signal.*")
    a("")
    a("---")
    a("")

    # ── Phase 2 Convergence Scoring ─────────────────────────────────────
    a("### Phase 2b — Convergence Scoring — Proximity to Degradation Motifs")
    a("")
    a("*Degradation pathway convergence: both viral and host fragments must approach*")
    a("*the same CMA machinery (hinges + jammers) to enable co-presentation by MHC.*")
    a("")

    p2 = result.phase2_enriched_pairs
    if p2:
        # Tier breakdown
        t1_cnt = sum(1 for p in p2 if p["phase2_tier"] == "T1")
        t2_cnt = sum(1 for p in p2 if p["phase2_tier"] == "T2")
        t3_cnt = sum(1 for p in p2 if p["phase2_tier"] == "T3")
        a("#### Tier Distribution")
        a("")
        a("| Tier | Criteria | Count |")
        a("|------|----------|-------|")
        a(f"| **T1** | total_event_signal ≥ 0.45 OR dual-sandwich OR (≥ 0.30 + single sandwich) | {t1_cnt} |")
        a(f"| **T2** | total_event_signal ≥ 0.15 | {t2_cnt} |")
        a(f"| **T3** | total_event_signal < 0.15 | {t3_cnt} |")
        a("")

        # Top 20 pairs sorted by convergence
        top_conv = sorted(p2, key=lambda x: x["convergence_score"], reverse=True)[:20]
        a("#### Top Pairs by Convergence Score")
        a("")
        a("*Positions are 1-based.  Score = final SEA score (Phase 2 + Phase 3 boosts, no epitope bias).  KE = ★ if host sequence overlaps a known autoantibody epitope (annotation only).*")
        a("")
        a("| # | Viral Seq | V-Pos | Host Seq | H-Pos | Score | Tier | Conv | V-ctx | H-ctx | V-dT1 | V-dJM | H-dT1 | H-dJM | Sandwich | KE |")
        a("|---|-----------|-------|----------|-------|-------|------|------|-------|-------|-------|-------|-------|-------|----------|----|")
        for i, p in enumerate(top_conv, 1):
            v_sw = "V" if p["viral_sandwich"] else "—"
            h_sw = "H" if p["host_sandwich"]  else "—"
            sw   = "+".join(filter(lambda x: x != "—", [v_sw, h_sw])) or "—"
            vdt1 = str(p["viral_dist_t1"]) if p["viral_dist_t1"] is not None else "—"
            vdjm = str(p["viral_dist_jammer"]) if p["viral_dist_jammer"] is not None else "—"
            hdt1 = str(p["host_dist_t1"]) if p["host_dist_t1"] is not None else "—"
            hdjm = str(p["host_dist_jammer"]) if p["host_dist_jammer"] is not None else "—"
            v_pos = p["position1"] + 1   # 0-based → 1-based
            h_pos = p["position2"] + 1
            score = p["final_sea_score"]
            ke    = "★" if p.get("overlaps_known_epitope") else "—"
            a(f"| {i} | `{p['seq1']}` | {v_pos} | `{p['seq2']}` | {h_pos} | {score:.4f} | {p['phase2_tier']} "
              f"| {p['convergence_score']:.3f} | {p['viral_context_score']:.3f} "
              f"| {p['host_context_score']:.3f} | {vdt1} | {vdjm} | {hdt1} | {hdjm} | {sw} | {ke} |")
    else:
        a("*Phase 2 convergence data not available.*")
    a("")
    a("---")
    a("")

    # ════════════════════════════════════════════════════════════
    # PHASE 3 — SEA Finder Analysis
    # ════════════════════════════════════════════════════════════
    a("## Phase 3 — SEA Finder Analysis")
    a("")
    a("*Source: SEA module — architecture scoring, hinge/jammer context, CMA weighting.*")
    a("")
    a("### Architecture Class Distribution")
    a("")
    arch_dist = _arch_distribution(result.sea_results)
    if arch_dist:
        a("| Architecture | Count |")
        a("|-------------|-------|")
        for arc, cnt in arch_dist.items():
            a(f"| {arc} | {cnt} |")
    else:
        a("*No SEA results.*")
    a("")
    a("### Known Epitope Overlap")
    a("")
    a("*Annotation only — these sequences are flagged because their host anchor overlaps a known*")
    a("*autoantibody epitope in the knowledge base.  No score boost is applied; this label serves*")
    a("*as ground-truth for calibration only.*")
    a("")
    if result.epitope_proximity_hits:
        a(f"**Host accession:** `{result.host_accession}`  |  "
          f"**Overlapping hits:** {len(result.epitope_proximity_hits)}")
        a("")
        a("| Viral Seq | Host Seq | Host Pos | Epitope | Range |")
        a("|-----------|----------|----------|---------|-------|")
        for h in result.epitope_proximity_hits:
            a(f"| `{h['seq1']}` | `{h['seq2']}` | "
              f"{h['host_pos_1based']} | {h['epitope_label']} | {h['epitope_range']} |")
    else:
        a("*No host sequences overlap a known autoantibody epitope range.*")
        if result.host_accession:
            a(f"*(Host accession: `{result.host_accession}`)*")
    a("")
    a("### Top SEA Hits (ranked by final_sea_score)")
    a("")
    # Build lookups from position pair → phase2_tier and overlaps_known_epitope
    _p2_tier_lookup: Dict[tuple, str] = {
        (p["position1"], p["position2"]): p["phase2_tier"]
        for p in result.phase2_enriched_pairs
    }
    _p2_ke_lookup: Dict[tuple, bool] = {
        (p["position1"], p["position2"]): p.get("overlaps_known_epitope", False)
        for p in result.phase2_enriched_pairs
    }
    top_n = result.sea_results[:20]
    if top_n:
        # Build phospho_exposure lookup from phase2_pairs
        _p2_pexp_lookup: Dict[tuple, float] = {
            (p["position1"], p["position2"]): p.get("phospho_exposure_score", 0.0)
            for p in result.phase2_enriched_pairs
        }
        a("| Rank | Viral Seq | Host Seq | V-Pos | H-Pos | Architecture | P2 Tier | Phospho-Exp | KE | Final Score |")
        a("|------|-----------|----------|-------|-------|-------------|---------|-------------|-----|-------------|")
        for i, r in enumerate(top_n, 1):
            arc   = r.architecture_class.name if hasattr(r.architecture_class, "name") else str(r.architecture_class)
            tier  = _p2_tier_lookup.get((r.position1, r.position2), "—")
            ke    = "★" if _p2_ke_lookup.get((r.position1, r.position2), False) else "—"
            pexp  = _p2_pexp_lookup.get((r.position1, r.position2), 0.0)
            v_pos = r.position1 + 1
            h_pos = r.position2 + 1
            a(f"| {i} | `{r.seq1}` | `{r.seq2}` | {v_pos} | {h_pos} | {arc} | {tier} | {pexp:.3f} | {ke} | {r.final_sea_score:.4f} |")
    else:
        a("*No SEA hits found.*")
    a("")
    a("### Complete-SEA Details")
    a("")
    super_hits = [r for r in result.sea_results if r.is_complete_sea]
    if super_hits:
        for i, r in enumerate(super_hits, 1):
            a(f"#### Super-Epitope #{i} — viral pos {r.position1} / host pos {r.position2}")
            a("")
            a(f"- **Viral fragment:** `{r.seq1}`")
            a(f"- **Host fragment:** `{r.seq2}`")
            a(f"- **Final score:** {r.final_sea_score:.4f}")
            a(f"- **Score components:** base={r.base_score:.3f}, hinge={r.hinge_score:.3f}, "
              f"jammer_density={r.jammer_density_score:.3f}, degradation={r.degradation_score:.3f}")
            a(f"- **Architecture bonus:** ×{r.architecture_bonus:.1f}  "
              f"Cell-type weight: ×{r.cell_type_weight:.1f}")
            if r.notes:
                a("- **Notes:**")
                for note in r.notes:
                    a(f"  - {note}")
            a("")
    else:
        a("*No super-epitope hits in this run.*")
    a("")
    a("---")
    a("")

    # ── TRUE SUPER-EPITOPE PAIRS ─────────────────────────────────────────
    a("### TRUE SUPER-EPITOPE Pairs")
    a("")
    a("*A TRUE SUPER-EPITOPE requires TWO viral hit fragments bridged by a T1 phospho-hinge.*")
    a("*When the hinge is phosphorylated, both fragments are simultaneously exposed,*")
    a("*acting as a combined large epitope that can be recognised by autoantibodies.*")
    a("*Viral span ≤ 150 aa | Host span ≤ 150 aa | ≥1 T1 hinge position between viral fragments.*")
    a("")
    sep_list = result.super_epitope_pairs
    if sep_list:
        a(f"**{len(sep_list)} TRUE SUPER-EPITOPE pair(s) detected.**")
        a("")
        a("| # | Hit-A Rank | Hit-B Rank | Hit-A Viral | Hit-B Viral | Hit-A Host | Hit-B Host | V-Span | H-Span | Bridging T1 Hinges | A-Score | B-Score |")
        a("|---|-----------|-----------|------------|------------|-----------|-----------|--------|--------|-------------------|---------|---------|")
        for i, sp in enumerate(sep_list[:10], 1):
            bridging_str = ",".join(str(p) for p in sp["bridging_t1_hinges"][:5])
            a(f"| {i} | #{sp['hit_a_rank']} | #{sp['hit_b_rank']} "
              f"| `{sp['hit_a_seq1']}` @{sp['hit_a_pos1']+1} "
              f"| `{sp['hit_b_seq1']}` @{sp['hit_b_pos1']+1} "
              f"| `{sp['hit_a_seq2']}` @{sp['hit_a_pos2']+1} "
              f"| `{sp['hit_b_seq2']}` @{sp['hit_b_pos2']+1} "
              f"| {sp['viral_span']} | {sp['host_span']} "
              f"| {bridging_str} "
              f"| {sp['hit_a_score']:.4f} | {sp['hit_b_score']:.4f} |")
    else:
        a("*No TRUE SUPER-EPITOPE pairs detected at this threshold (viral ≤ 150 aa, host ≤ 150 aa, T1-bridged).*")
        a("*Interpretation: no two hit fragments are brought into proximity by a single T1 hinge event.*")
    a("")
    a("---")
    a("")

    # ── RANKING TALLY ────────────────────────────────────────────────────
    a("### Ranking Tally — Known Epitope Sequences")
    a("")
    a("*Tracks current rank of sequences that overlap the known autoantibody epitope (★ KE annotation).*")
    a("*Use this table across tuning iterations to monitor convergence toward placing known immunogenic*")
    a("*sequences at the top of the prioritised list.*")
    a("")
    # Build a fresh position→rank lookup from sorted sea_results
    _rank_lookup: Dict[tuple, int] = {
        (r.position1, r.position2): r.rank
        for r in result.sea_results
    }
    ke_tally = [
        p for p in result.phase2_enriched_pairs
        if p.get("overlaps_known_epitope")
    ]
    if ke_tally:
        ke_tally_sorted = sorted(ke_tally, key=lambda p: _rank_lookup.get((p["position1"], p["position2"]), 9999))
        a("| # | Viral Seq → Host Seq | V-Pos | H-Pos | Rank | Score | Phospho-Exp | KFERQ-Prox | LIR-Prox | P2 Tier |")
        a("|---|---------------------|-------|-------|------|-------|-------------|-----------|---------|---------|")
        for i, p in enumerate(ke_tally_sorted, 1):
            cur_rank = _rank_lookup.get((p["position1"], p["position2"]), "?")
            v_pos    = p["position1"] + 1
            h_pos    = p["position2"] + 1
            kf_flag  = "YES" if p.get("proximal_kferq_motif") else "—"
            lir_flag = "YES" if p.get("proximal_lir_motif")   else "—"
            pexp     = p.get("phospho_exposure_score", 0.0)
            a(f"| {i} | `{p['seq1']}` → `{p['seq2']}` | {v_pos} | {h_pos} "
              f"| **#{cur_rank}** | {p['final_sea_score']:.4f} | {pexp:.3f} "
              f"| {kf_flag} | {lir_flag} | {p['phase2_tier']} |")
    else:
        a("*No sequences overlap a known autoantibody epitope — supply host_accession and")
        a(" protein_knowledge_base.json to enable this tally.*")
    a("")
    a("---")
    a("")

    # ════════════════════════════════════════════════════════════
    # PHASE 4 — 3D Structural Analysis  (Pending)
    # ════════════════════════════════════════════════════════════
    a("## Phase 4 — 3D Structural Analysis")
    a("")
    a("> **Status: Pending** — methodology not yet implemented.")
    a("")
    a("Phase 4 will assess whether a SEA-identified mimicry pair organises two")
    a("discontinuous fragments of the host protein into a functional")
    a("**super-epitope architecture** after hinge-mediated phosphorylation and")
    a("unfolding of the intervening sequence.  Candidate pairs for Phase 4 review")
    a("are listed below.")
    a("")
    a("### Candidate Hits for Phase 4 (COMPLETE_SEA, score ≥ 10.0)")
    a("")
    ph4_candidates = [
        r for r in result.sea_results
        if r.is_complete_sea and r.final_sea_score >= 10.0
    ]
    if ph4_candidates:
        a("| # | Viral Seq | Host Seq | Viral Pos | Host Pos | Score |")
        a("|---|-----------|----------|-----------|----------|-------|")
        for i, r in enumerate(ph4_candidates, 1):
            a(f"| {i} | `{r.seq1}` | `{r.seq2}` | {r.position1} | {r.position2} | {r.final_sea_score:.4f} |")
    else:
        a("*No candidates meeting threshold (COMPLETE_SEA + score ≥ 10.0) found in this run.*")
    a("")
    a("### Planned Methodology")
    a("")
    a("1. Fetch PDB structure(s) for host protein via RCSB API.")
    a("2. Extract Cα coordinates for each Phase 3 super-epitope fragment.")
    a("3. Identify intervening hinge residues (from Phase 2 Tier-1 scan).")
    a("4. Model phospho-hinge unfolding and estimate fragment proximity (Å).")
    a("5. Score geometric feasibility of a **conformational super-epitope**.")
    a("6. Pass candidate structures to Phase 6 (Data Visualization) for PDB rendering.")
    a("")
    a("---")
    a("")

    # ════════════════════════════════════════════════════════════
    # PHASE 5 — Orchestrator Synthesis
    # ════════════════════════════════════════════════════════════
    a("## Phase 5 — Orchestrator Synthesis")
    a("")
    a("*Cross-phase key findings, risk interpretation, and recommended follow-on actions.*")
    a("")
    a("### Key Findings")
    a("")
    overall_risk = rs.get("overall_risk", "UNKNOWN")
    top_score    = rs.get("top_final_score", 0.0)
    n_super      = rs.get("n_complete_sea", 0)
    n_pairs      = rs.get("n_pairs_scored", 0)
    n_ep_hits    = rs.get("n_epitope_proximity_hits", 0)
    top_arc      = rs.get("top_architecture", "N/A")

    a(f"1. **Risk level:** {overall_risk} — top SEA score {top_score:.4f} across {n_pairs} scored pairs.")
    if n_super > 0:
        a(f"2. **Complete-SEA architecture detected:** {n_super} hit(s) classified as COMPLETE_SEA, "
          f"indicating a single viral fragment is sandwiched by structural features AND flanked by CMA motifs.")
    else:
        a("2. **No COMPLETE_SEA architecture detected** in this scoring run.")
    if n_ep_hits > 0:
        a(f"3. **Known epitope overlap (annotation):** {n_ep_hits} hit(s) whose host anchor "
          f"overlaps a mapped autoantibody epitope for `{result.host_accession}`. "
          f"These sequences were NOT score-boosted; their ranking is driven entirely by "
          f"Phase 2 convergence + Phase 3 architecture scoring.")
    else:
        a("3. **No known epitope overlap** detected with mapped autoantibody epitope ranges.")

    # Phase 2a finding — KFERQ proximity
    n_proximal = sum(1 for p in result.phase2_enriched_pairs if p.get("proximal_deg_motif"))
    if n_proximal > 0:
        a(f"6. **KFERQ-proximal pairs (Phase 2a):** {n_proximal} homologous pair(s) lie within "
          f"30 residues of a KFERQ-like CMA-targeting motif. These are the highest-priority "
          f"mechanistic candidates — they represent direct intersection of molecular mimicry "
          f"with the host chaperone-mediated autophagy pathway. Review Phase 2a table above.")
    else:
        a("6. **No KFERQ-proximal pairs detected** at the 30-residue threshold. "
          f"SEA priority is driven by Phase 2b convergence and Phase 3 architecture alone.")
    v_t1 = vhj.get("n_hinges_t1", 0)
    v_t2 = vhj.get("n_hinges_t2", 0)
    v_jm = vhj.get("n_jammers", 0)
    h_t1 = hhj.get("n_hinges_t1", 0)
    h_t2 = hhj.get("n_hinges_t2", 0)
    h_jm = hhj.get("n_jammers", 0)
    a(f"4. **Viral architectural complexity:** {v_t1} Tier-1 hinges, {v_t2} Tier-2 hinges, "
      f"{v_jm} jammers — a high-jammer density increases immune evasion probability.")
    a(f"5. **Host architectural context:** {h_t1} Tier-1 hinges, {h_t2} Tier-2 hinges, "
      f"{h_jm} jammers — hinge-flanked regions are priority Phase 4 candidates.")
    a("")
    a("### Risk Interpretation")
    a("")
    if overall_risk == "CRITICAL":
        a("The CRITICAL risk designation indicates the viral protein carries sequence motifs")
        a("capable of mimicking host epitopes at multiple positions, with structural features")
        a("(COMPLETE_SEA architecture, high jammer density) that suggest active immune evasion.")
        a("Priority action: Phase 4 structural validation of top super-epitope candidates.")
    elif overall_risk == "HIGH":
        a("HIGH risk: significant mimicry potential detected. Functional validation of top")
        a("SEA hits is warranted before clinical or experimental extrapolation.")
    elif overall_risk == "MODERATE":
        a("MODERATE risk: mimicry signals present but below threshold for high-confidence")
        a("clinical interpretation.  Additional sequence coverage (seq_aligner with BioPython)")
        a("may clarify borderline hits.")
    else:
        a("LOW / UNKNOWN risk: insufficient evidence for actionable mimicry signal in this run.")
    a("")
    a("### Pre-Phase 4 Calibration Check")
    a("")
    cal = result.calibration_result
    if cal:
        cal_status = cal.get("status", "N/A")
        cal_msg    = cal.get("message", "")
        cal_top    = cal.get("top_n", 10)
        cal_hits   = cal.get("known_epitope_hits", [])
        cal_best   = cal.get("best_known_rank")
        status_fmt = f"**{cal_status}**" if cal_status in ("PASS", "FAIL") else cal_status
        a(f"**Status:** {status_fmt}  |  "
          f"**Message:** {cal_msg}")
        a("")
        a("*Scores shown are unbiased — no epitope proximity bonus was applied.*")
        a("*A PASS here means Phase 2 convergence + Phase 3 architecture scoring*")
        a("*surfaced known-epitope sequences on their own merits.*")
        a("")
        if cal_hits:
            a(f"*Known-epitope hits in top {cal_top} results (unbiased ranking):*")
            a("")
            a("| Rank | Host Pos | Epitope | Range | Score |")
            a("|------|----------|---------|-------|-------|")
            for h in cal_hits:
                a(f"| {h['rank']} | {h['host_pos']} | {h['epitope_label']} | "
                  f"{h['epitope_range']} | {h['score']:.4f} |")
            a("")
        if cal_status == "FAIL":
            a("> :warning: **Calibration FAIL** — unbiased scoring did not surface")
            a("> known autoantibody epitope positions in the top results.")
            a("> Consider: relaxing mclachlan_min_composite, broadening epitope ranges")
            a("> in protein_knowledge_base.json, or tuning SEAConfig thresholds.")
            a("")
    else:
        a("*Calibration check not run (no host_accession supplied).*")
    a("")
    a("### Recommended Follow-On Actions")
    a("")
    a("- [ ] **Phase 4:** Fetch PDB structures and run geometric super-epitope assessment")
    a(f"      on {len(ph4_candidates)} candidate hit(s) listed above.")
    a("- [ ] **Install BioPython** in the analysis environment to enable seq_aligner pairs,")
    a("      recovering the ~12 pairs currently missing from the scoring pool.")
    a("- [ ] **Phase 6:** Generate PDB-mapped visualization of hinges, jammers, and top")
    a("      SEA hits once Phase 4 structural coordinates are available.")
    a("- [ ] **Update protein_knowledge_base.json** with any newly validated autoantibody")
    a("      epitope positions to sharpen future epitope_proximity_bonus scoring.")
    a("")
    a("---")
    a("")

    # ── Footer ────────────────────────────────────────────────────────────
    a("> *Generated by SEA Orchestrator v6 — BiologicExplorer/Tools*  ")
    a(f"> *Phases 1–3 complete (Phase 3 v2: phospho-exposure + TRUE SUPER-EPITOPE pairs + ranking tally) | Phase 4 pending | Phase 6 visualization queued*")

    return lines

# ════════════════════════════════════════════════════════════════════════════
#  PHASE 2 — DEGRADATION PATHWAY CONVERGENCE HELPERS
# ════════════════════════════════════════════════════════════════════════════

def _min_dist(pos: int, elements: List[Any], seq_window: int = 7) -> Optional[int]:
    """
    Return the minimum residue distance from *pos* (0-based, fragment start)
    to the nearest element in *elements* (each must have a ``.position`` attr).

    Distance is measured centre-to-centre: the fragment centre is
    ``pos + seq_window // 2``.  Returns None when *elements* is empty.
    """
    if not elements:
        return None
    center = pos + seq_window // 2
    return min(abs(center - e.position) for e in elements)


def _phospho_exposure_score(
    sea_result: Any,
    window:     int = 15,
) -> float:
    """
    Convert SEAResult.phospho_t1_hinge_dist (computed by sea_module inside
    its per-pair proximity window) to a 0-1 phospho-exposure score.

    The orchestrator does NOT re-scan hinge positions — that work is owned by
    sea_module.  This function is a thin normalisation layer only.

    Score = max(0.0, 1.0 - dist / window), rounded to 4 dp.
    Returns 0.0 when sea_module found no T1 hinge within its proximity window.
    """
    d = getattr(sea_result, "phospho_t1_hinge_dist", None)
    if d is None:
        return 0.0
    return round(max(0.0, 1.0 - d / window), 4)



def _compute_pair_degradation_convergence(
    viral_pos:       int,
    host_pos:        int,
    viral_hj:        Dict,
    host_hj:         Dict,
    seq_window:      int   = 7,
    sandwich_window: int   = 20,
    decay_factor:    float = 10.0,
) -> Dict:
    """
    Compute the degradation pathway convergence score for one pair.

    Context score formula (per side):
        context = min(1.0,
            0.70 × max(0, 1 − dist_t1  / decay_factor)
          + 0.45 × max(0, 1 − dist_t2  / decay_factor)
          + 0.30 × max(0, 1 − dist_jmr / decay_factor))

    Convergence = viral_context × host_context   (multiplicative)

    Sandwich flag: pair centre lies between a jammer and a hinge that are
    both within *sandwich_window* residues of each other on the same side.

    Tier assignment:
        T1 : convergence ≥ 0.45  OR  (v_sandwich AND h_sandwich)
             OR  (convergence ≥ 0.30 AND (v_sandwich OR h_sandwich))
        T2 : convergence ≥ 0.15
        T3 : below 0.15

    Returns a flat dict with all distances, scores, flags, and tier.
    """

    def _context(pos: int, hj: Dict, win: int):
        t1d = _min_dist(pos, hj.get("hinges_t1", []), win)
        t2d = _min_dist(pos, hj.get("hinges_t2", []), win)
        jmd = _min_dist(pos, hj.get("jammers",   []), win)

        def contrib(d, w):
            if d is None:
                return 0.0
            return w * max(0.0, 1.0 - d / decay_factor)

        score = contrib(t1d, 0.70) + contrib(t2d, 0.45) + contrib(jmd, 0.30)
        return min(1.0, score), t1d, t2d, jmd

    def _sandwich(pos: int, hj: Dict, win: int) -> bool:
        """True if pos centre lies between a (jammer, hinge) pair ≤ sandwich_window apart."""
        jammers = hj.get("jammers",   [])
        hinges  = hj.get("hinges_t1", []) + hj.get("hinges_t2", [])
        center  = pos + win // 2
        for jm in jammers:
            for hg in hinges:
                lo   = min(jm.position, hg.position)
                hi   = max(jm.position, hg.position)
                span = hi - lo
                if span <= sandwich_window and lo <= center <= hi:
                    return True
        return False

    v_ctx, v_t1d, v_t2d, v_jmd = _context(viral_pos, viral_hj, seq_window)
    h_ctx, h_t1d, h_t2d, h_jmd = _context(host_pos,  host_hj,  seq_window)

    convergence = round(v_ctx * h_ctx, 4)

    v_sand = _sandwich(viral_pos, viral_hj, seq_window)
    h_sand = _sandwich(host_pos,  host_hj,  seq_window)

    if (convergence >= 0.45
            or (v_sand and h_sand)
            or (convergence >= 0.30 and (v_sand or h_sand))):
        tier = "T1"
    elif convergence >= 0.15:
        tier = "T2"
    else:
        tier = "T3"

    return {
        "viral_dist_t1":       v_t1d,
        "viral_dist_t2":       v_t2d,
        "viral_dist_jammer":   v_jmd,
        "host_dist_t1":        h_t1d,
        "host_dist_t2":        h_t2d,
        "host_dist_jammer":    h_jmd,
        "viral_context_score": round(v_ctx, 4),
        "host_context_score":  round(h_ctx, 4),
        "convergence_score":   convergence,
        "viral_sandwich":      v_sand,
        "host_sandwich":       h_sand,
        "phase2_tier":         tier,
    }


def _check_degradation_proximity(
    viral_pos:        int,
    host_pos:         int,
    viral_kferq:      List[Dict],
    host_kferq:       List[Dict],
    seq_window:       int        = 7,
    proximity_window: int        = 30,
    viral_lir:        List[Dict] = None,   # find_all_degradation_motifs() lir key
    host_lir:         List[Dict] = None,
) -> Dict:
    """
    Phase 2a — KFERQ motif proximity check for a homologous pair.

    For each side (viral, host), compute the centre-to-centre distance from
    the fragment centre to the nearest KFERQ-like motif.

    Fragment centre  = pos + seq_window // 2  (0-based).
    KFERQ motif centre = (start - 1) + KFERQ_LEN // 2  (start is 1-based
    from find_kferq_motifs(); converted here to 0-based).

    Parameters
    ----------
    viral_pos, host_pos  : 0-based fragment start positions.
    viral_kferq          : find_kferq_motifs() output for the viral sequence.
    host_kferq           : find_kferq_motifs() output for the host sequence.
    seq_window           : fragment length (default 7, matching SEAModule).
    proximity_window     : residues; pair is flagged proximal if either
                           distance ≤ this threshold (default 30).

    Returns
    -------
    {
      "viral_kferq_dist":   int | None,
      "host_kferq_dist":    int | None,
      "proximal_deg_motif": bool,
    }
    """
    _KFERQ_LEN = 5  # canonical/phospho/acetyl pentapeptide

    def _nearest_kferq(center: int, motifs: List[Dict]) -> Optional[int]:
        if not motifs:
            return None
        return min(
            abs(center - ((m["start"] - 1) + _KFERQ_LEN // 2))
            for m in motifs
        )

    v_center = viral_pos + seq_window // 2
    h_center = host_pos  + seq_window // 2

    v_dist = _nearest_kferq(v_center, viral_kferq)
    h_dist = _nearest_kferq(h_center, host_kferq)

    proximal = (
        (v_dist is not None and v_dist <= proximity_window) or
        (h_dist is not None and h_dist <= proximity_window)
    )

    # ── LIR proximity ────────────────────────────────────────────────────
    _LIR_LEN = 4   # canonical LIR core: [WFY]xx[ILV]

    def _nearest_lir(center: int, motifs: List[Dict]) -> Optional[int]:
        if not motifs:
            return None
        return min(
            abs(center - ((m["start"] - 1) + _LIR_LEN // 2))
            for m in motifs
        )

    vl_dist = _nearest_lir(v_center, viral_lir or [])
    hl_dist = _nearest_lir(h_center, host_lir  or [])

    proximal_kferq = (
        (v_dist  is not None and v_dist  <= proximity_window) or
        (h_dist  is not None and h_dist  <= proximity_window)
    )
    proximal_lir = (
        (vl_dist is not None and vl_dist <= proximity_window) or
        (hl_dist is not None and hl_dist <= proximity_window)
    )
    proximal_combined = proximal_kferq or proximal_lir

    return {
        "viral_kferq_dist":    v_dist,
        "host_kferq_dist":     h_dist,
        "viral_lir_dist":      vl_dist,
        "host_lir_dist":       hl_dist,
        "proximal_kferq_motif": proximal_kferq,
        "proximal_lir_motif":   proximal_lir,
        "proximal_deg_motif":   proximal_combined,   # combined for backwards compat
    }


def _compute_proximity_event_score(
    dp: Dict,
    window: int  = 30,
    w_viral_lir:   float = 1.0,
    w_host_lir:    float = 1.5,
    w_viral_kferq: float = 1.0,
    w_host_kferq:  float = 2.0,
) -> Dict:
    """
    Distance-weighted additive proximity event score.

    For each of the four degradation-event signals (viral LIR, host LIR,
    viral KFERQ, host KFERQ) a component score is computed as:

        component = max(0, 1 - dist / window) × weight

    where *dist* is the distance in residues to the nearest motif of that
    type (``None`` → treated as out-of-window).  The four components are
    summed and normalised to [0, 1] by dividing by the maximum possible
    weight sum.

    Parameters
    ----------
    dp     : dict returned by ``_check_degradation_proximity()``
    window : decay window in residues (default 30)
    w_*    : per-signal weights

    Returns
    -------
    dict with keys:
        viral_lir_score, host_lir_score,
        viral_kferq_score, host_kferq_score,
        proximity_event_score   (0–1, normalised)
    """
    max_possible = w_viral_lir + w_host_lir + w_viral_kferq + w_host_kferq

    def _component(dist_key: str, weight: float) -> float:
        dist = dp.get(dist_key)
        if dist is None:
            return 0.0
        return max(0.0, 1.0 - dist / window) * weight

    vl = _component("viral_lir_dist",   w_viral_lir)
    hl = _component("host_lir_dist",    w_host_lir)
    vk = _component("viral_kferq_dist", w_viral_kferq)
    hk = _component("host_kferq_dist",  w_host_kferq)

    raw_sum   = vl + hl + vk + hk
    normalised = raw_sum / max_possible if max_possible > 0 else 0.0

    return {
        "viral_lir_score":       round(vl / w_viral_lir   if w_viral_lir   else 0.0, 4),
        "host_lir_score":        round(hl / w_host_lir    if w_host_lir    else 0.0, 4),
        "viral_kferq_score":     round(vk / w_viral_kferq if w_viral_kferq else 0.0, 4),
        "host_kferq_score":      round(hk / w_host_kferq  if w_host_kferq  else 0.0, 4),
        "proximity_event_score": round(min(1.0, normalised), 4),
    }


def _detect_super_epitope_pairs(
    sea_results:   List[Any],
    viral_hj:      Dict,
    viral_window:  int = 150,
    host_window:   int = 150,
) -> List[Dict]:
    """
    Detect TRUE SUPER-EPITOPE pairs: two ranked viral hit fragments where:
      1. The viral fragment start positions are within *viral_window* residues.
      2. The host  fragment start positions are within *host_window*  residues.
      3. At least one T1 hinge position lies BETWEEN the two viral fragments
         in the linear sequence (strict interior).

    Returns a list of dicts (one per detected pair), sorted by combined score.
    """
    t1_positions = [h.position for h in viral_hj.get("hinges_t1", [])]
    pairs: List[Dict] = []

    n = len(sea_results)
    for i in range(n):
        for j in range(i + 1, n):
            r1 = sea_results[i]
            r2 = sea_results[j]
            v_lo = min(r1.position1, r2.position1)
            v_hi = max(r1.position1, r2.position1)
            h_lo = min(r1.position2, r2.position2)
            h_hi = max(r1.position2, r2.position2)

            viral_span = v_hi - v_lo
            host_span  = h_hi - h_lo

            if viral_span > viral_window or host_span > host_window:
                continue

            bridging = [p for p in t1_positions if v_lo < p < v_hi]
            if not bridging:
                continue

            pairs.append({
                "hit_a_rank":         r1.rank,
                "hit_b_rank":         r2.rank,
                "hit_a_seq1":         r1.seq1,
                "hit_b_seq1":         r2.seq1,
                "hit_a_seq2":         r1.seq2,
                "hit_b_seq2":         r2.seq2,
                "hit_a_pos1":         r1.position1,
                "hit_b_pos1":         r2.position1,
                "hit_a_pos2":         r1.position2,
                "hit_b_pos2":         r2.position2,
                "viral_span":         viral_span,
                "host_span":          host_span,
                "bridging_t1_hinges": bridging,
                "hit_a_score":        r1.final_sea_score,
                "hit_b_score":        r2.final_sea_score,
            })

    # Sort by combined score of the two hits
    pairs.sort(key=lambda p: p["hit_a_score"] + p["hit_b_score"], reverse=True)
    return pairs


def _calibration_check(
    sea_results:    List[Any],   # List[SEAResult]
    host_accession: str,
    kb_data:        Dict,
    top_n:          int = 10,
) -> Dict:
    """
    Pre-Phase 4 calibration check.

    Verifies that at least one known autoantibody epitope position for
    *host_accession* appears in the top *top_n* SEA results (after all
    score boosts have been applied).

    Returns
    -------
    Dict with keys:
        status            : "PASS" | "FAIL" | "NO_EPITOPES"
        top_n             : int (as supplied)
        best_known_rank   : Optional[int] — 1-based rank of best known-epitope hit
        known_epitope_hits: List[Dict]    — hits in top_n that overlap known epitope
        epitopes_checked  : int           — number of linear epitopes loaded
        message           : str           — human-readable summary
    """
    epitope_ranges = get_epitope_ranges(kb_data, host_accession)
    if not epitope_ranges:
        return {
            "status":              "NO_EPITOPES",
            "top_n":               top_n,
            "best_known_rank":     None,
            "known_epitope_hits":  [],
            "epitopes_checked":    0,
            "message": (
                f"No linear epitopes found in knowledge base for {host_accession}."
            ),
        }

    known_hits: List[Dict] = []
    for rank_1b, r in enumerate(sea_results[:top_n], 1):
        host_pos = r.position2 + 1   # 0-based → 1-based
        for ep_start, ep_end, ep_label in epitope_ranges:
            if ep_start <= host_pos <= ep_end:
                known_hits.append({
                    "rank":          rank_1b,
                    "host_pos":      host_pos,
                    "epitope_label": ep_label,
                    "epitope_range": f"{ep_start}–{ep_end}",
                    "score":         r.final_sea_score,
                })
                break   # one entry per result

    best_rank = min((h["rank"] for h in known_hits), default=None)

    if known_hits:
        status  = "PASS"
        message = (
            f"PASS — {len(known_hits)} known-epitope hit(s) in top {top_n}; "
            f"best at rank #{best_rank}."
        )
    else:
        status  = "FAIL"
        message = (
            f"FAIL — no known-epitope host positions found in top {top_n} "
            f"SEA results.  Methodology refinement recommended."
        )

    return {
        "status":              status,
        "top_n":               top_n,
        "best_known_rank":     best_rank,
        "known_epitope_hits":  known_hits,
        "epitopes_checked":    len(epitope_ranges),
        "message":             message,
    }


# ════════════════════════════════════════════════════════════════════════════
#  MAIN ENTRY POINT
# ════════════════════════════════════════════════════════════════════════════

def orchestrate(
    viral_seq:                 str,
    host_seq:                  str,
    virus_name:                str,
    viral_protein_name:        str,
    host_protein_name:         str,
    # McLachlan pair source
    mclachlan_hits:            Optional[List[Dict]] = None,
    mclachlan_min_composite:   float = 15.0,
    mclachlan_max_pairs:       int   = 200,
    # Smith-Waterman pair source
    sw_hits:                   Optional[List[Dict]] = None,
    sw_min_composite:          float = 14.5,
    sw_max_pairs:              int   = 200,
    # Sequence aligner pair source
    seq_aligner_min_identity:  float = 0.33,
    seq_aligner_window:        int   = 12,
    seq_aligner_step:          int   = 4,
    seq_aligner_max_pairs:     int   = 40,
    # Optional expression data for CMA scoring
    expression_dict:           Optional[Dict[str, float]] = None,
    reference_dict:            Optional[Dict[str, float]] = None,
    # SEA config override
    sea_config:                Optional[SEAConfig] = None,
    # Knowledge base / epitope proximity bonus
    host_accession:            Optional[str]   = None,
    kb_path:                   Optional[str]   = None,
    epitope_proximity_bonus:   float           = 2.0,
    # Phase 2 — convergence scoring
    convergence_bonus_weight:  float           = 3.0,
    # Phase 2a — KFERQ proximity window (residues)
    deg_proximity_window:      int             = 30,
    # Phase 2a — additive proximity event weight (multiplier applied to
    # the normalised proximity_event_score before adding to final_sea_score)
    proximity_event_weight:    float           = 5.0,
) -> OrchestratorResult:
    """
    Run the full SEA Orchestrator pipeline.

    Parameters
    ----------
    viral_seq               : full viral protein sequence
    host_seq                : full host protein sequence
    virus_name              : short virus name (e.g. 'HCV') — sets APC tropism weight
    viral_protein_name      : descriptive name (e.g. 'polyprotein Q9WMX2')
    host_protein_name       : descriptive name (e.g. 'CYP2E1 P05181')
    mclachlan_hits          : pre-computed McLachlan hit dicts (from v3 JSON 'hits' key)
                              Pass None to skip this pair source
    mclachlan_min_composite : minimum composite_primary score to include a McLachlan hit
    mclachlan_max_pairs     : maximum McLachlan pairs to use
    sw_hits                 : pre-computed Smith-Waterman hit dicts (compatible schema)
                              Pass None to skip this pair source
    sw_min_composite        : minimum composite_primary score to include an SW hit
    sw_max_pairs            : maximum SW pairs to use  
    seq_aligner_min_identity: minimum identity fraction for find_homologous_pairs()
    seq_aligner_window      : sliding window length for find_homologous_pairs()
    seq_aligner_step        : stride for find_homologous_pairs()
    seq_aligner_max_pairs   : max pairs from find_homologous_pairs()
    expression_dict         : gene expression values for CMA score (optional)
    reference_dict          : reference expression values (optional)
    sea_config              : custom SEAConfig (uses defaults if None)
    host_accession          : UniProt primary accession of the host protein
                              (e.g. 'P05181').  When supplied, the pipeline
                              loads protein_knowledge_base.json and adds
                              epitope_proximity_bonus to any SEAResult whose
                              host anchor falls inside a known autoantibody
                              epitope range.  Pass None to skip this step.
    kb_path                 : path to protein_knowledge_base.json.  Defaults
                              to orchestrator/protein_knowledge_base.json.
    epitope_proximity_bonus : score bonus added when a hit lands in a known
                              epitope (default 2.0).  Set 0.0 to disable while
                              still supplying host_accession.

    Returns
    -------
    OrchestratorResult
    """

    # ── Step 1 & 2: Motif profiles ────────────────────────────────────────
    viral_profile = find_all_degradation_motifs(viral_seq)
    host_profile  = find_all_degradation_motifs(host_seq)

    # ── Step 1c: Full-sequence hinge/jammer profiles ──────────────────────
    viral_hj = _scan_hinge_jammer_profile(viral_seq, sea_config)
    host_hj  = _scan_hinge_jammer_profile(host_seq,  sea_config)

    # ── Step 3: Convert host KFERQ motifs → SEA degradation_motifs format ─
    # find_kferq_motifs returns: {start (1-based), end, motif, type, notes}
    # SEAModule.run() expects:   {motif, position (0-based), protein}
    host_kferq_raw = find_kferq_motifs(host_seq)
    provided_motifs: List[Dict] = [
        {
            "motif":    m["motif"],
            "position": m["start"] - 1,   # 1-based → 0-based
            "protein":  "protein2",
        }
        for m in host_kferq_raw
    ]

    # ── Step 4: CMA membership ────────────────────────────────────────────
    host_cma = check_cma_network_membership(host_protein_name)

    # ── Step 5a: Sequence aligner pairs ───────────────────────────────────
    try:
        aligner_pairs = find_homologous_pairs(
            viral_seq,
            host_seq,
            min_identity = seq_aligner_min_identity,
            window       = seq_aligner_window,
            step         = seq_aligner_step,
            max_pairs    = seq_aligner_max_pairs,
        )
        # Tag source
        for p in aligner_pairs:
            p["source"] = "seq_aligner"
    except ImportError:
        aligner_pairs = []

    n_seq_aligner = len(aligner_pairs)

    # ── Step 5b: McLachlan pairs ──────────────────────────────────────────
    if mclachlan_hits is not None:
        mc_pairs = mclachlan_to_pairs(
            mclachlan_hits,
            viral_seq,
            min_composite = mclachlan_min_composite,
            max_pairs     = mclachlan_max_pairs,
        )
    else:
        mc_pairs = []

    n_mclachlan = len(mc_pairs)

    # ── Step 5c: Smith-Waterman pairs (optional third source) ─────────────
    if sw_hits is not None:
        sw_pairs = mclachlan_to_pairs(
            sw_hits,
            viral_seq,
            min_composite = sw_min_composite,
            max_pairs     = sw_max_pairs,
            source        = "sw",
        )
    else:
        sw_pairs = []

    n_sw = len(sw_pairs)

    # ── Step 6: Merge, dedup, re-rank ─────────────────────────────────────
    merged = aligner_pairs + mc_pairs + sw_pairs
    merged_after_dedup = _dedup_pairs(merged, tolerance=4)

    # Assign consecutive 1-based ranks (sorted by similarity_score DESC)
    for i, p in enumerate(merged_after_dedup, 1):
        p["rank"] = i

    pair_source_counts = {
        "seq_aligner": n_seq_aligner,
        "mclachlan":   n_mclachlan,
        "sw":           n_sw,
        "merged":      len(merged),
        "after_dedup": len(merged_after_dedup),
    }

    # ── Step 7: SEA scoring ───────────────────────────────────────────────
    sea_module  = SEAModule(virus_name=virus_name, config=sea_config)
    sea_results = sea_module.run(
        homologous_pairs    = merged_after_dedup,
        autoimmune_diseases = [],
        degradation_motifs  = provided_motifs,
    )
    # sea_results is already sorted by final_sea_score DESC by SEAModule.run()

    # ── Knowledge base (loaded once; used by Steps 7b, 7d) ───────────────
    kb_data: Dict = {}
    if host_accession is not None:
        kb_data = load_protein_knowledge_base(kb_path)

    # ── Step 7b: Known epitope overlap annotation (no score boost) ───────
    # Sequences that overlap a known autoantibody epitope are flagged for
    # traceability.  The annotation is ONLY used as a post-hoc label — it
    # does not alter final_sea_score.  Ranking is driven entirely by Phase 2
    # convergence scoring (Step 7c), keeping calibration unbiased.
    proximity_hits: List[Dict] = []
    if host_accession is not None:
        proximity_hits = annotate_epitope_overlap(sea_results, host_accession, kb_data)

    # Build position lookup for phase2_pairs stamping (Step 7c)
    _overlap_pos_set: set = {
        (h["position1"], h["position2"]) for h in proximity_hits
    }

    # ── Step 7b-prime: Phase 2a — KFERQ motif proximity check ────────────
    # For each homologous pair from Phase 1, compute the distance to the
    # nearest KFERQ-like degradation motif on the viral AND host side.
    # This is the primary mechanistic filter: a viral sequence that mimics
    # a host sequence near a KFERQ motif is the highest-priority SEA event.
    # Annotation only — no score mutation.  A proximal pair is floored at
    # Tier 2 in Step 7c to ensure it is never buried below T3.
    viral_kferq_raw   = viral_profile.get("kferq", [])
    host_kferq_raw_p2 = find_kferq_motifs(host_seq)   # reuse raw 1-based list
    viral_lir_raw     = viral_profile.get("lir", [])
    host_lir_raw      = host_profile.get("lir", [])
    _deg_proximity_map: Dict[tuple, Dict] = {}
    for r in sea_results:
        dp = _check_degradation_proximity(
            viral_pos        = r.position1,
            host_pos         = r.position2,
            viral_kferq      = viral_kferq_raw,
            host_kferq       = host_kferq_raw_p2,
            seq_window       = 7,
            proximity_window = deg_proximity_window,
            viral_lir        = viral_lir_raw,
            host_lir         = host_lir_raw,
        )
        _deg_proximity_map[(r.position1, r.position2)] = dp
        if dp["proximal_kferq_motif"]:
            r.notes.append(
                f"Phase2a KFERQ-proximal "
                f"(v_dist={dp['viral_kferq_dist']}, "
                f"h_dist={dp['host_kferq_dist']})"
            )
        if dp["proximal_lir_motif"]:
            r.notes.append(
                f"Phase2a LIR-proximal "
                f"(v_dist={dp['viral_lir_dist']}, "
                f"h_dist={dp['host_lir_dist']})"
            )

    # ── Step 7c: Phase 2 — degradation pathway convergence scoring ────────
    phase2_pairs: List[Dict] = []
    for r in sea_results:
        conv = _compute_pair_degradation_convergence(
            viral_pos       = r.position1,
            host_pos        = r.position2,
            viral_hj        = viral_hj,
            host_hj         = host_hj,
            seq_window      = 7,
            sandwich_window = 20,
            decay_factor    = 10.0,
        )
        boost = round(conv["convergence_score"] * convergence_bonus_weight, 4)
        if boost > 0.0:
            r.final_sea_score = round(r.final_sea_score + boost, 4)
            r.notes.append(
                f"Phase2 convergence boost +{boost:.4f} "
                f"(v_ctx={conv['viral_context_score']:.3f} × "
                f"h_ctx={conv['host_context_score']:.3f}, "
                f"tier={conv['phase2_tier']})"
            )

        # ── Phase 2a: Additive proximity event scoring ──────────────────
        dp = _deg_proximity_map.get((r.position1, r.position2), {})
        prox = _compute_proximity_event_score(
            dp, window=deg_proximity_window
        )
        prox_boost = round(prox["proximity_event_score"] * proximity_event_weight, 4)
        if prox_boost > 0:
            r.final_sea_score = round(r.final_sea_score + prox_boost, 4)
            r.notes.append(
                f"Phase2a proximity event boost +{prox_boost:.4f} "
                f"(prox_score={prox['proximity_event_score']:.4f}, "
                f"weight={proximity_event_weight})"
            )

        # ── Tier assignment (output-only, no floors) ─────────────────────
        # total_event_signal combines convergence + proximity on [0, 1].
        # Tier 1 = highest priority (alphanumerically first).
        v_sand = conv.get("viral_sandwich", False)
        h_sand = conv.get("host_sandwich",  False)
        total_event_signal = min(
            1.0,
            conv["convergence_score"] + prox["proximity_event_score"]
        )
        if (
            total_event_signal >= 0.45
            or (v_sand and h_sand)
            or (total_event_signal >= 0.30 and (v_sand or h_sand))
        ):
            effective_tier = "T1"
        elif total_event_signal >= 0.15:
            effective_tier = "T2"
        else:
            effective_tier = "T3"


        phospho_exp = _phospho_exposure_score(r)

        phase2_pairs.append({
            "rank":                   r.rank,
            "seq1":                   r.seq1,
            "seq2":                   r.seq2,
            "position1":              r.position1,
            "position2":              r.position2,
            "convergence_score":      conv["convergence_score"],
            "phase2_tier":            effective_tier,
            "viral_context_score":    conv["viral_context_score"],
            "host_context_score":     conv["host_context_score"],
            "viral_dist_t1":          conv["viral_dist_t1"],
            "viral_dist_t2":          conv["viral_dist_t2"],
            "viral_dist_jammer":      conv["viral_dist_jammer"],
            "host_dist_t1":           conv["host_dist_t1"],
            "host_dist_t2":           conv["host_dist_t2"],
            "host_dist_jammer":       conv["host_dist_jammer"],
            "viral_sandwich":         conv["viral_sandwich"],
            "host_sandwich":          conv["host_sandwich"],
            "convergence_boost":      boost,
            "final_sea_score":        r.final_sea_score,
            "overlaps_known_epitope": (r.position1, r.position2) in _overlap_pos_set,
            "viral_kferq_dist":       dp.get("viral_kferq_dist"),
            "host_kferq_dist":        dp.get("host_kferq_dist"),
            "viral_lir_dist":         dp.get("viral_lir_dist"),
            "host_lir_dist":          dp.get("host_lir_dist"),
            "proximal_kferq_motif":   dp.get("proximal_kferq_motif", False),
            "proximal_lir_motif":     dp.get("proximal_lir_motif",   False),
            "proximal_deg_motif":     dp.get("proximal_deg_motif",   False),
            "viral_lir_score":        prox["viral_lir_score"],
            "host_lir_score":         prox["host_lir_score"],
            "viral_kferq_score":      prox["viral_kferq_score"],
            "host_kferq_score":       prox["host_kferq_score"],
            "proximity_event_score":  prox["proximity_event_score"],
            "proximity_boost":        prox_boost,
            "total_event_signal":     round(total_event_signal, 4),
            "phospho_exposure_score": phospho_exp,
        })

    # Re-sort and re-rank after convergence boosts
    sea_results.sort(key=lambda r: r.final_sea_score, reverse=True)
    for i, r in enumerate(sea_results, 1):
        r.rank = i

    # ── Step 7e: TRUE SUPER-EPITOPE pair detection ───────────────────────
    super_ep_pairs = _detect_super_epitope_pairs(sea_results, viral_hj)

    # ── Step 7d: Pre-Phase 4 calibration check ────────────────────────────
    calibration: Dict = {}
    if host_accession is not None:
        calibration = _calibration_check(
            sea_results, host_accession, kb_data, top_n=10
        )

    # ── Step 8: Optional CMA score ────────────────────────────────────────
    cma_score_result: Optional[Dict] = None
    cma_score_val:    Optional[float] = None
    if expression_dict is not None:
        try:
            cma_score_result = calculate_cma_score(expression_dict, reference_dict)
            cma_score_val    = cma_score_result.get("cma_score")
        except Exception:
            pass   # expression data failures do not block the pipeline

    # ── Step 9: Risk summary ──────────────────────────────────────────────
    risk_summary = _build_risk_summary(
        sea_results,
        host_cma,
        viral_profile,
        host_profile,
        cma_score_result,
        viral_hj_profile = viral_hj,
        host_hj_profile  = host_hj,
    )
    # Patch in epitope proximity stats (not available inside _build_risk_summary)
    risk_summary["n_epitope_proximity_hits"] = len(proximity_hits)
    risk_summary["epitope_proximity_bonus"]  = (
        epitope_proximity_bonus if host_accession is not None else 0.0
    )

    return OrchestratorResult(
        virus_name                       = virus_name,
        viral_protein_name               = viral_protein_name,
        host_protein_name                = host_protein_name,
        sea_results                      = sea_results,
        viral_motif_profile              = viral_profile,
        host_motif_profile               = host_profile,
        host_cma_membership              = host_cma,
        cma_score                        = cma_score_val,
        risk_summary                     = risk_summary,
        pair_source_counts               = pair_source_counts,
        host_accession                   = host_accession,
        epitope_proximity_hits           = proximity_hits,
        epitope_proximity_bonus_applied  = 0.0,   # annotation only — no score boost applied
        viral_hinge_jammer_profile       = viral_hj,
        host_hinge_jammer_profile        = host_hj,
        phase2_enriched_pairs            = phase2_pairs,
        calibration_result               = calibration,
        super_epitope_pairs              = super_ep_pairs,
    )
