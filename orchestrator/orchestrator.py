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
8.  Optional calculate_cma_score()        → cma_score
9.  Build OrchestratorResult + risk_summary

Usage
-----
>>> from orchestrator.orchestrator import orchestrate
>>> result = orchestrate(viral_seq, host_seq, "HCV", "polyprotein", "CYP2E1")
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
from typing import Any, Dict, List, Optional

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
    from sea.sea_module import SEAModule, SEAConfig, find_homologous_pairs
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
#  McLACHLAN CONVERSION
# ════════════════════════════════════════════════════════════════════════════

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
    sea_results:       List[Any],
    host_cma:          Optional[Dict],
    viral_profile:     Dict,
    host_profile:      Dict,
    cma_score_result:  Optional[Dict],
) -> Dict:
    """Build a concise risk summary dict from Orchestrator outputs."""

    if not sea_results:
        return {
            "top_architecture":   "NONE",
            "top_final_score":    0.0,
            "n_super_epitope":    0,
            "n_sandwiched":       0,
            "n_pairs_scored":     0,
            "host_is_cma_member": False,
            "host_cma_category":  None,
            "viral_n_kferq":      0,
            "host_n_kferq":       0,
            "overall_risk":       "LOW",
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
        "top_architecture":   top_arc,
        "top_final_score":    round(score, 4),
        "n_super_epitope":    supers,
        "n_sandwiched":       sands,
        "n_pairs_scored":     len(sea_results),
        "host_is_cma_member": host_cma_member,
        "host_cma_category":  host_cma_category,
        "viral_n_kferq":      v_kferq,
        "host_n_kferq":       h_kferq,
        "overall_risk":       overall_risk,
        "cma_score":          (
            cma_score_result.get("cma_score") if cma_score_result else None
        ),
    }


# ════════════════════════════════════════════════════════════════════════════
#  REPORT BUILDER
# ════════════════════════════════════════════════════════════════════════════

def _build_report_lines(result: OrchestratorResult) -> List[str]:
    """Return a list of Markdown lines for the Orchestrator report."""
    rs = result.risk_summary
    psc = result.pair_source_counts

    lines: List[str] = []
    a = lines.append   # shorthand

    a(f"# SEA Orchestrator Report: {result.virus_name} vs {result.host_protein_name}")
    a(f"")
    a(f"**Generated:** {result.timestamp}")
    a(f"**Viral protein:** {result.viral_protein_name}  |  **Host protein:** {result.host_protein_name}")
    a(f"")

    # ── Risk Summary ──────────────────────────────────────────────────────
    a("## Risk Summary")
    a("")
    a(f"| Field | Value |")
    a(f"|-------|-------|")
    a(f"| Overall Risk | **{rs.get('overall_risk', 'N/A')}** |")
    a(f"| Top Architecture | {rs.get('top_architecture', 'N/A')} |")
    a(f"| Top Final Score | {rs.get('top_final_score', 0.0):.4f} |")
    a(f"| Super-Epitope Hits | {rs.get('n_super_epitope', 0)} |")
    a(f"| Sandwiched Hits | {rs.get('n_sandwiched', 0)} |")
    a(f"| Pairs Scored | {rs.get('n_pairs_scored', 0)} |")
    a(f"| Host CMA Member | {rs.get('host_is_cma_member', False)} ({rs.get('host_cma_category', 'N/A')}) |")
    a(f"| Viral KFERQ motifs | {rs.get('viral_n_kferq', 0)} |")
    a(f"| Host KFERQ motifs | {rs.get('host_n_kferq', 0)} |")
    if rs.get("cma_score") is not None:
        a(f"| CMA Score | {rs['cma_score']:.4f} |")
    a("")

    # ── Pair Source Counts ────────────────────────────────────────────────
    a("## Pair Source Counts")
    a("")
    a(f"| Source | Count |")
    a(f"|--------|-------|")
    for k, v in psc.items():
        a(f"| {k} | {v} |")
    a("")

    # ── Top SEA Hits ──────────────────────────────────────────────────────
    a("## Top SEA Hits (ranked by final_sea_score)")
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

    # ── Viral Motif Profile ───────────────────────────────────────────────
    a("## Viral Degradation Motif Profile")
    a("")
    vs = result.viral_motif_profile.get("summary", {})
    a(f"- CMA category: `{vs.get('cma_category', 'N/A')}`")
    a(f"- KFERQ motifs: {vs.get('n_kferq', 0)}")
    a(f"- LIR motifs: {vs.get('n_lir', 0)}")
    a(f"- D-box degrons: {vs.get('n_dbox', 0)}")
    a(f"- Destabilising N-terminus: {vs.get('has_destabilising_n_term', False)}")
    a(f"- Destabilising C-terminus: {vs.get('has_destabilising_c_term', False)}")
    a("")

    # ── Host Motif Profile ────────────────────────────────────────────────
    a("## Host Degradation Motif Profile")
    a("")
    hs = result.host_motif_profile.get("summary", {})
    a(f"- CMA category: `{hs.get('cma_category', 'N/A')}`")
    a(f"- KFERQ motifs: {hs.get('n_kferq', 0)}")
    a(f"- LIR motifs: {hs.get('n_lir', 0)}")
    a(f"- D-box degrons: {hs.get('n_dbox', 0)}")
    a(f"- Destabilising N-terminus: {hs.get('has_destabilising_n_term', False)}")
    a(f"- Destabilising C-terminus: {hs.get('has_destabilising_c_term', False)}")
    a("")

    # ── CMA Membership ────────────────────────────────────────────────────
    a("## Host CMA Network Membership")
    a("")
    if result.host_cma_membership:
        hc = result.host_cma_membership
        a(f"- Symbol: `{hc.get('symbol', 'N/A')}`")
        a(f"- Category: {hc.get('category', 'N/A')}")
        a(f"- Direction: {'+1 (CMA activating)' if hc.get('direction', 0) == 1 else '-1 (CMA inhibiting)'}")
        a(f"- Weight: {hc.get('weight', 'N/A')}")
    else:
        a(f"*{result.host_protein_name} is not a known CMA network member.*")
    a("")

    # ── Notes on all Super-Epitopes ───────────────────────────────────────
    super_hits = [r for r in result.sea_results if r.is_super_epitope]
    if super_hits:
        a("## Super-Epitope Details")
        a("")
        for i, r in enumerate(super_hits, 1):
            a(f"### Super-Epitope #{i} — pos {r.position1} (viral) / pos {r.position2} (host)")
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

    # ── Footer ────────────────────────────────────────────────────────────
    a("---")
    a("*Generated by SEA Orchestrator — BiologicExplorer/Tools*")

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

    Returns
    -------
    OrchestratorResult
    """

    # ── Step 1 & 2: Motif profiles ────────────────────────────────────────
    viral_profile = find_all_degradation_motifs(viral_seq)
    host_profile  = find_all_degradation_motifs(host_seq)

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
    )

    return OrchestratorResult(
        virus_name          = virus_name,
        viral_protein_name  = viral_protein_name,
        host_protein_name   = host_protein_name,
        sea_results         = sea_results,
        viral_motif_profile = viral_profile,
        host_motif_profile  = host_profile,
        host_cma_membership = host_cma,
        cma_score           = cma_score_val,
        risk_summary        = risk_summary,
        pair_source_counts  = pair_source_counts,
    )
