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
        source           str    'mclachlan'
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
                "source":           "mclachlan",
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
            "n_super_epitope":         0,
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
    supers  = sum(1 for r in sea_results if r.is_super_epitope)
    sands   = sum(1 for r in sea_results if r.is_sandwiched and not r.is_super_epitope)
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
        "n_super_epitope":         supers,
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
    a(f"| Super-Epitope Hits | {rs.get('n_super_epitope', 0)} |")
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
    a("### Epitope Proximity Bonus")
    a("")
    if result.epitope_proximity_hits:
        a(f"**Host accession:** `{result.host_accession}`  |  "
          f"**Bonus per hit:** +{result.epitope_proximity_bonus_applied:.2f}  |  "
          f"**Hits boosted:** {len(result.epitope_proximity_hits)}")
        a("")
        a("| Rank | Viral Seq | Host Seq | Host Pos | Epitope | Range | New Score |")
        a("|------|-----------|----------|----------|---------|-------|-----------|")
        for h in result.epitope_proximity_hits:
            a(f"| {h['rank']} | `{h['seq1']}` | `{h['seq2']}` | "
              f"{h['host_pos_1based']} | {h['epitope_label']} | "
              f"{h['epitope_range']} | {h['new_score']:.4f} |")
    else:
        a("*No hits received an epitope proximity bonus.*")
        if result.host_accession:
            a(f"*(Host accession: `{result.host_accession}`, "
              f"bonus configured: {result.epitope_proximity_bonus_applied:.2f})*")
    a("")
    a("### Top SEA Hits (ranked by final_sea_score)")
    a("")
    top_n = result.sea_results[:20]
    if top_n:
        a("| Rank | Viral Seq | Host Seq | Viral Pos | Host Pos | Architecture | Final Score |")
        a("|------|-----------|----------|-----------|----------|-------------|-------------|")
        for i, r in enumerate(top_n, 1):
            arc = r.architecture_class.name if hasattr(r.architecture_class, "name") else str(r.architecture_class)
            a(f"| {i} | `{r.seq1}` | `{r.seq2}` | {r.position1} | {r.position2} | {arc} | {r.final_sea_score:.4f} |")
    else:
        a("*No SEA hits found.*")
    a("")
    a("### Super-Epitope Details")
    a("")
    super_hits = [r for r in result.sea_results if r.is_super_epitope]
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
    a("### Candidate Hits for Phase 4 (SUPER_EPITOPE, score ≥ 10.0)")
    a("")
    ph4_candidates = [
        r for r in result.sea_results
        if r.is_super_epitope and r.final_sea_score >= 10.0
    ]
    if ph4_candidates:
        a("| # | Viral Seq | Host Seq | Viral Pos | Host Pos | Score |")
        a("|---|-----------|----------|-----------|----------|-------|")
        for i, r in enumerate(ph4_candidates, 1):
            a(f"| {i} | `{r.seq1}` | `{r.seq2}` | {r.position1} | {r.position2} | {r.final_sea_score:.4f} |")
    else:
        a("*No candidates meeting threshold (SUPER_EPITOPE + score ≥ 10.0) found in this run.*")
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
    n_super      = rs.get("n_super_epitope", 0)
    n_pairs      = rs.get("n_pairs_scored", 0)
    n_ep_hits    = rs.get("n_epitope_proximity_hits", 0)
    top_arc      = rs.get("top_architecture", "N/A")

    a(f"1. **Risk level:** {overall_risk} — top SEA score {top_score:.4f} across {n_pairs} scored pairs.")
    if n_super > 0:
        a(f"2. **Super-epitope architecture detected:** {n_super} hit(s) classified as SUPER_EPITOPE, "
          f"indicating the viral sequence may scaffold a multi-component mimicry epitope.")
    else:
        a("2. **No super-epitope architecture detected** in this scoring run.")
    if n_ep_hits > 0:
        a(f"3. **Epitope proximity overlap:** {n_ep_hits} hit(s) land inside a known autoantibody "
          f"epitope for host protein `{result.host_accession}` (+{result.epitope_proximity_bonus_applied:.2f} bonus each).")
    else:
        a("3. **No epitope proximity overlap** detected with known autoantibody epitope ranges.")
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
        a("(SUPER_EPITOPE architecture, high jammer density) that suggest active immune evasion.")
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
    a("> *Generated by SEA Orchestrator v5 — BiologicExplorer/Tools*  ")
    a(f"> *Phases 1–3 complete | Phase 4 pending | Phase 6 visualization queued*")

    return lines

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

    # ── Step 6: Merge, dedup, re-rank ─────────────────────────────────────
    merged = aligner_pairs + mc_pairs
    merged_after_dedup = _dedup_pairs(merged, tolerance=4)

    # Assign consecutive 1-based ranks (sorted by similarity_score DESC)
    for i, p in enumerate(merged_after_dedup, 1):
        p["rank"] = i

    pair_source_counts = {
        "seq_aligner": n_seq_aligner,
        "mclachlan":   n_mclachlan,
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

    # ── Step 7b: Epitope proximity bonus ──────────────────────────────────
    proximity_hits: List[Dict] = []
    if host_accession is not None and epitope_proximity_bonus != 0.0:
        kb_data = load_protein_knowledge_base(kb_path)
        sea_results, proximity_hits = apply_epitope_proximity_bonus(
            sea_results,
            host_accession,
            kb_data,
            epitope_proximity_bonus,
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
        epitope_proximity_bonus_applied  = (
            epitope_proximity_bonus if host_accession is not None else 0.0
        ),
        viral_hinge_jammer_profile       = viral_hj,
        host_hinge_jammer_profile        = host_hj,
    )
