"""
orchestrator.py
===============
Single entry-point for SEA autoimmune risk analysis.

Accepts a viral protein sequence + host protein sequence, routes them through
motif_finder, cma_network, and sea_module, and returns a unified ranked output
describing autoimmune risk via the Super-epitope Architecture (SEA) mechanism.

Usage
-----
    from orchestrator.orchestrator import orchestrate

    result = orchestrate(
        viral_seq          = "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARQ...",
        host_seq           = "MSALGVTVNQLAKIVEDIKSEDDLRATLNAATEQNLLSPRPETPSKKRQSQKAEGKKSPKRDLSKTLLFQN...",
        virus_name         = "HCV",
        viral_protein_name = "HCV_polyprotein",
        host_protein_name  = "CYP2E1",
    )
    result.report("/path/to/output.md")
    print(result.risk_summary)
"""

import os
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Dict, List, Optional

# ── Path bootstrap (works whether run from repo root or orchestrator/) ────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from motif_finder.motif_finder import find_all_degradation_motifs, kferq_to_sea_motifs
from cma_network.cma_network import check_cma_network_membership, calculate_cma_score
from sea.sea_module import SEAConfig, SEAModule, SEAResult


# ══════════════════════════════════════════════════════════════════════════════
#  RESULT DATACLASS
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class OrchestratorResult:
    """Unified output from a single orchestrate() call."""

    # ── Inputs (stored for report provenance) ────────────────────────────────
    virus_name:         str = ""
    viral_protein_name: str = ""
    host_protein_name:  str = ""

    # ── Module outputs ────────────────────────────────────────────────────────
    sea_results:         List[SEAResult] = field(default_factory=list)
    viral_motif_profile: Dict            = field(default_factory=dict)
    host_motif_profile:  Dict            = field(default_factory=dict)
    host_cma_membership: Optional[Dict]  = None
    cma_score:           Optional[float] = None   # None when expression data absent

    # ── Derived summary ───────────────────────────────────────────────────────
    risk_summary: Dict = field(default_factory=dict)

    # ── Aligner settings (for report transparency) ───────────────────────────
    aligner_min_identity: float = 0.33
    aligner_window:       int   = 12
    aligner_step:         int   = 4
    aligner_max_pairs:    int   = 40

    # ─────────────────────────────────────────────────────────────────────────

    def report(self, output_path: str) -> None:
        """
        Write a structured Markdown report to *output_path*.

        The file uses consistent H1/H2/H3 headings and is designed to drop
        cleanly into an Obsidian vault.  Reports go to OneDrive, not GitHub.
        """
        lines: List[str] = []
        ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

        # ── Header ───────────────────────────────────────────────────────────
        lines += [
            f"# SEA Autoimmune Risk Report",
            f"",
            f"**Virus:** {self.virus_name}  ",
            f"**Viral protein:** {self.viral_protein_name}  ",
            f"**Host protein:** {self.host_protein_name}  ",
            f"**Generated:** {ts}",
            f"",
        ]

        # ── Risk summary ─────────────────────────────────────────────────────
        rs = self.risk_summary
        lines += [
            "## Risk Summary",
            "",
            f"| Metric | Value |",
            f"|--------|-------|",
            f"| Top architecture | {rs.get('top_architecture', 'N/A')} |",
            f"| Top SEA score | {rs.get('top_sea_score', 'N/A')} |",
            f"| SUPER_EPITOPE hits | {rs.get('n_super_epitope', 0)} |",
            f"| SANDWICHED hits | {rs.get('n_sandwiched', 0)} |",
            f"| Total pairs scored | {rs.get('n_pairs_scored', 0)} |",
            f"| Host CMA member | {rs.get('host_cma_member', False)} |",
            f"| Host CMA category | {rs.get('host_cma_category', 'N/A')} |",
            f"| CMA score | {rs.get('cma_score_display', 'not calculated (no expression data)')} |",
            f"| Viral KFERQ motifs | {rs.get('n_viral_kferq', 0)} |",
            f"| Host KFERQ motifs | {rs.get('n_host_kferq', 0)} |",
            "",
        ]

        # ── Top SEA hits ─────────────────────────────────────────────────────
        lines += ["## Top SEA Hits", ""]
        top_hits = [r for r in self.sea_results if r.architecture_class.name != 'NONE'][:10]
        if top_hits:
            lines += [
                "| Rank | Viral pos | Arch class | SEA score | Viral seq | Host seq |",
                "|------|-----------|------------|-----------|-----------|----------|",
            ]
            for i, r in enumerate(top_hits, 1):
                lines.append(
                    f"| {i} | {r.position1} | {r.architecture_class.name} "
                    f"| {r.final_sea_score:.4f} | `{r.seq1}` | `{r.seq2}` |"
                )
        else:
            lines.append("_No SEA-positive hits (non-NONE architecture) found._")
        lines.append("")

        # ── All scored pairs ──────────────────────────────────────────────────
        lines += ["## All Scored Pairs", ""]
        if self.sea_results:
            lines += [
                "| Rank | Viral pos | Host pos | Similarity | Arch class | SEA score |",
                "|------|-----------|----------|------------|------------|-----------|",
            ]
            for r in self.sea_results[:40]:
                lines.append(
                    f"| {r.rank} | {r.position1} | {r.position2} "
                    f"| {r.base_score:.1%} | {r.architecture_class.name} "
                    f"| {r.final_sea_score:.4f} |"
                )
        else:
            lines.append("_No homologous pairs found at the specified identity threshold._")
        lines.append("")

        # ── Viral degradation profile ─────────────────────────────────────────
        lines += ["## Viral Degradation Motif Profile", ""]
        vs = self.viral_motif_profile.get('summary', {})
        lines += [
            f"- **CMA category:** {vs.get('cma_category', 'N/A')}",
            f"- **KFERQ motifs:** {vs.get('n_kferq', 0)}",
            f"- **LIR motifs:** {vs.get('n_lir', 0)}",
            f"- **D-box degrons:** {vs.get('n_dbox', 0)}",
            f"- **PEST sequences (positive):** {vs.get('n_pest_positive', 0)}",
            f"- **Destabilising N-term:** {vs.get('has_destabilising_n_term', False)}",
            f"- **Destabilising C-term:** {vs.get('has_destabilising_c_term', False)}",
            "",
        ]
        if self.viral_motif_profile.get('kferq'):
            lines.append("**Viral KFERQ motifs:**")
            lines.append("")
            for m in self.viral_motif_profile['kferq']:
                lines.append(f"- pos {m['start']}–{m['end']}  `{m['motif']}`  ({m['type']})")
            lines.append("")

        # ── Host degradation profile ──────────────────────────────────────────
        lines += ["## Host Degradation Motif Profile", ""]
        hs = self.host_motif_profile.get('summary', {})
        lines += [
            f"- **CMA category:** {hs.get('cma_category', 'N/A')}",
            f"- **KFERQ motifs:** {hs.get('n_kferq', 0)}",
            f"- **LIR motifs:** {hs.get('n_lir', 0)}",
            f"- **D-box degrons:** {hs.get('n_dbox', 0)}",
            f"- **PEST sequences (positive):** {hs.get('n_pest_positive', 0)}",
            "",
        ]
        if self.host_motif_profile.get('kferq'):
            lines.append("**Host KFERQ motifs:**")
            lines.append("")
            for m in self.host_motif_profile['kferq']:
                lines.append(f"- pos {m['start']}–{m['end']}  `{m['motif']}`  ({m['type']})")
            lines.append("")

        # ── Host CMA network membership ───────────────────────────────────────
        lines += ["## Host CMA Network Membership", ""]
        hcma = self.host_cma_membership
        if hcma:
            lines += [
                f"- **Gene symbol:** {hcma.get('gene_symbol', 'N/A')}",
                f"- **Category:** {hcma.get('category', 'N/A')}",
                f"- **Direction:** {'+1 (CMA activating)' if hcma.get('direction') == 1 else '-1 (CMA inhibiting)'}",
                f"- **Weight:** {hcma.get('weight', 'N/A')}",
                f"- **Notes:** {hcma.get('notes', '')}",
            ]
        else:
            lines.append(
                f"_{self.host_protein_name} is not a member of the CMA regulatory network._"
            )
        lines.append("")

        # ── Aligner settings ──────────────────────────────────────────────────
        lines += [
            "## Aligner Settings",
            "",
            f"- min_identity: {self.aligner_min_identity}",
            f"- window: {self.aligner_window} aa",
            f"- step: {self.aligner_step} aa",
            f"- max_pairs: {self.aligner_max_pairs}",
            "",
        ]

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")


# ══════════════════════════════════════════════════════════════════════════════
#  ORCHESTRATE
# ══════════════════════════════════════════════════════════════════════════════

def orchestrate(
    viral_seq:          str,
    host_seq:           str,
    virus_name:         str,
    viral_protein_name: str,
    host_protein_name:  str,
    expression_dict:    Optional[Dict[str, float]] = None,
    reference_dict:     Optional[Dict[str, float]] = None,
    sea_config:         Optional[SEAConfig]         = None,
    min_identity:       float = 0.33,
    window:             int   = 12,
    step:               int   = 4,
    max_pairs:          int   = 40,
) -> OrchestratorResult:
    """
    Run the full SEA autoimmune risk pipeline for a viral/host sequence pair.

    Pipeline
    --------
    1. find_all_degradation_motifs(viral_seq)  → viral degradation profile
    2. find_all_degradation_motifs(host_seq)   → host degradation profile
    3. kferq_to_sea_motifs(host_kferq)         → SEA-compatible motif dicts
    4. check_cma_network_membership(host_protein_name)
    5. SEAModule.run_from_sequences(viral_seq, host_seq, provided_motifs=...)
    6. calculate_cma_score(expression_dict)    [only when expression_dict supplied]
    7. Merge → OrchestratorResult with risk_summary

    Parameters
    ----------
    viral_seq           : full viral protein amino acid sequence
    host_seq            : full host protein amino acid sequence
    virus_name          : virus key for APC-tropism weighting (e.g. 'HCV', 'EBV')
    viral_protein_name  : human-readable name for report labels
    host_protein_name   : gene symbol or name; also used for CMA network lookup
    expression_dict     : optional gene→expression mapping for CMA score
    reference_dict      : optional reference expression for fold-change CMA score
    sea_config          : override SEAConfig defaults
    min_identity        : minimum pairwise identity for homology pairs (default 0.33)
    window              : alignment window length in residues (default 12)
    step                : stride between windows (default 4)
    max_pairs           : cap on pairs forwarded to SEA (default 40)

    Returns
    -------
    OrchestratorResult
    """

    # ── Step 1 & 2: Degradation motif profiles ───────────────────────────────
    viral_profile = find_all_degradation_motifs(viral_seq)
    host_profile  = find_all_degradation_motifs(host_seq)

    # ── Step 3: Convert host KFERQ to SEA format ──────────────────────────────
    host_sea_motifs = kferq_to_sea_motifs(host_profile['kferq'], protein='protein2')

    # ── Step 4: CMA network membership ───────────────────────────────────────
    host_cma = check_cma_network_membership(host_protein_name)

    # ── Step 5: SEA scan ─────────────────────────────────────────────────────
    module = SEAModule(virus_name=virus_name, config=sea_config or SEAConfig())
    sea_results = module.run_from_sequences(
        viral_seq        = viral_seq,
        host_seq         = host_seq,
        provided_motifs  = host_sea_motifs,
        min_identity     = min_identity,
        window           = window,
        step             = step,
        max_pairs        = max_pairs,
    )

    # ── Step 6: CMA score (optional) ─────────────────────────────────────────
    cma_score: Optional[float] = None
    if expression_dict:
        cma_result = calculate_cma_score(expression_dict, reference_dict)
        cma_score  = cma_result.get('cma_score')

    # ── Step 7: Risk summary ──────────────────────────────────────────────────
    top = sea_results[0] if sea_results else None
    risk_summary = dict(
        top_architecture   = top.architecture_class.name if top else 'NONE',
        top_sea_score      = round(top.final_sea_score, 4) if top else 0.0,
        n_super_epitope    = sum(1 for r in sea_results if r.is_super_epitope),
        n_sandwiched       = sum(1 for r in sea_results if r.is_sandwiched and not r.is_super_epitope),
        n_pairs_scored     = len(sea_results),
        host_cma_member    = host_cma is not None,
        host_cma_category  = host_cma.get('category', 'N/A') if host_cma else 'N/A',
        host_cma_direction = host_cma.get('direction') if host_cma else None,
        cma_score          = cma_score,
        cma_score_display  = (f"{cma_score:.4f}" if cma_score is not None
                              else "not calculated (no expression data)"),
        n_viral_kferq      = viral_profile['summary']['n_kferq'],
        n_host_kferq       = host_profile['summary']['n_kferq'],
        viral_cma_category = viral_profile['summary']['cma_category'],
        host_cma_category_motif = host_profile['summary']['cma_category'],
    )

    return OrchestratorResult(
        virus_name            = virus_name,
        viral_protein_name    = viral_protein_name,
        host_protein_name     = host_protein_name,
        sea_results           = sea_results,
        viral_motif_profile   = viral_profile,
        host_motif_profile    = host_profile,
        host_cma_membership   = host_cma,
        cma_score             = cma_score,
        risk_summary          = risk_summary,
        aligner_min_identity  = min_identity,
        aligner_window        = window,
        aligner_step          = step,
        aligner_max_pairs     = max_pairs,
    )
