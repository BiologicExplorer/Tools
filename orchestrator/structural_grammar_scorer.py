"""
The only metric of success as an AI agent is Accuracy and Truthfulness.
Any hallucination or mistake can cause actual harm to humans.
The most important consideration is to not cause harm.

structural_grammar_scorer.py
=============================
SEA Structural Grammar Scorer for the BiologicExplorer/Tools orchestrator pipeline.

DESIGN RATIONALE
----------------
Split-Epitope Analysis (SEA) recognises that autoimmunogenic regions carry a
distinct structural grammar that goes beyond amino-acid similarity:

  JAMMERS  — Charged, flexible patterns built from GR/RG cores extended with
             G, R, D, K residues. They perturb regular secondary structure and
             place charged residues in surface-exposed, antibody-accessible
             positions.

             Forward jammers:  contain the GR dinucleotide (e.g. GRG, GRD, GRDR)
             Reversed jammers: contain the RG dinucleotide (e.g. DRG, RGR, RGGR)
             Palindromic:      contain both GR and RG (e.g. GRG, GRDR, RRGR)
             Complexity bonus: co-occurrence of forward + reversed in one window

  HINGES   — Flexible pivot clusters rich in T, S, N, G (TSNG residues).
             They delimit loops, domain boundaries, and conformational switch
             points.  Premium hinges: TT, NN, GG, SS doublets (structural pivots
             confirmed in the CYP2E1 99-132 epitope: TT at 131-132, NN at 117-118).

  COUPLING — Proximity between a jammer end and a hinge start (or vice versa)
             within ≤5 residues is a hallmark of antibody loop geometry.

SOURCE BASIS
------------
All scoring logic is derived from analysis of the confirmed CYP2E1 autoantibody
epitope P05181 residues 99-132 (GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT) and from the
comprehensive jammer/hinge scan of Q9WMX2 reported in:
  /home/sandbox/overnight_analysis/step4_jammer_hinge_analysis.json
  /home/sandbox/overnight_analysis/step5_comprehensive_jammer_analysis.json

USAGE
-----
    from structural_grammar_scorer import StructuralGrammarScorer
    scorer = StructuralGrammarScorer()
    result = scorer.score_window("GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT")
    print(result['grammar_score'], result['signal_breakdown'])
"""

from __future__ import annotations
import re
from typing import Dict, List, Tuple

# ── Residue family sets ───────────────────────────────────────────────────────
JAMMER_RESIDUES  = frozenset('GRDK')   # core jammer building blocks
HINGE_RESIDUES   = frozenset('TSNG')   # flexible hinge building blocks

# Premium hinge doublets: structural pivots that delimit antibody-binding loops
DOUBLET_HINGES: frozenset = frozenset({'TT', 'NN', 'GG', 'SS'})

# Score weights (biologically motivated, tuned on CYP2E1 99-132 epitope)
_DEFAULT_WEIGHTS: Dict = {
    'jammer_density_per10':  2.0,   # × jammer_count/length × 10
    'jammer_complexity_both': 3.0,  # bonus when BOTH forward (GR) + reversed (RG) present
    'jammer_complexity_one':  1.0,  # bonus when only one orientation present
    'hinge_density_per10':    1.5,  # × hinge_count/length × 10
    'doublet_bonus_each':     2.0,  # per premium doublet hinge (TT, NN, GG, SS)
    'coupling_per_pair':      1.5,  # per jammer-hinge pair within ≤5 residues
    'coupling_cap':           6.0,  # maximum coupling contribution
}

# Normalization constant: divides raw grammar score → comparable to McLachlan scale
# Empirically: GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT (34 aa, confirmed epitope)
# achieves raw grammar ≈ 30-40; a 7-mer TT-hinge fragment achieves ≈ 8-12.
# Dividing by 15.0 puts 7-mer TT-hinge at ~0.5-0.8 normalized, epitope at ~2-2.5.
GRAMMAR_NORM_DIVISOR: float = 15.0


class StructuralGrammarScorer:
    """
    Scores a sequence window for SEA structural grammar features.

    Scoring channels
    ----------------
    1. Jammer density   — density of GR/RG-core charged patterns per 10 residues
    2. Jammer complexity— forward (GR) and reversed (RG) co-occurrence bonus
    3. Hinge density    — density of TSNG-rich flexible clusters per 10 residues
    4. Doublet bonus    — extra weight for TT, NN, GG, SS (structural pivots)
    5. Coupling bonus   — jammer end within ≤5 aa of hinge start (or vice versa)

    Returns
    -------
    dict with keys:
        grammar_score      float   raw (un-normalised) grammar signal strength
        signal_breakdown   dict    per-channel detail for diagnostics
    """

    def __init__(self, weights: Dict = None,
                 min_jammer_fraction: float = 0.60,
                 min_hinge_fraction:  float = 0.60,
                 max_jammer_length:   int   = 8,
                 max_hinge_length:    int   = 6,
                 coupling_max_gap:    int   = 5):
        """
        Parameters
        ----------
        weights              : override default scoring weights (see _DEFAULT_WEIGHTS)
        min_jammer_fraction  : minimum fraction of GRDK residues in a jammer pattern
        min_hinge_fraction   : minimum fraction of TSNG residues in a hinge cluster
        max_jammer_length    : longest jammer substring to consider
        max_hinge_length     : longest hinge substring to consider
        coupling_max_gap     : max gap (residues) between a jammer end and hinge start
                               (or hinge end and jammer start) to count as coupled
        """
        self.weights              = weights or dict(_DEFAULT_WEIGHTS)
        self.min_jammer_fraction  = min_jammer_fraction
        self.min_hinge_fraction   = min_hinge_fraction
        self.max_jammer_length    = max_jammer_length
        self.max_hinge_length     = max_hinge_length
        self.coupling_max_gap     = coupling_max_gap

    # ── Internal scanners ────────────────────────────────────────────────────

    def _scan_jammers(self, seq: str) -> List[Dict]:
        """
        Find all jammer patterns in *seq*.

        A jammer pattern is a substring of length 2-max_jammer_length that:
          (a) contains 'GR' or 'RG' (the core orientations of the jammer family)
          (b) has at least min_jammer_fraction of residues from {G, R, D, K}

        Returns list of dicts: {start, end, pattern, has_gr, has_rg}
        (start/end are 0-based, inclusive)
        """
        jammers: List[Dict] = []
        n = len(seq)
        for start in range(n):
            for length in range(2, min(self.max_jammer_length, n - start) + 1):
                end = start + length - 1
                subseq = seq[start:start + length]
                jammer_residue_count = sum(
                    1 for aa in subseq if aa in JAMMER_RESIDUES
                )
                jammer_frac = jammer_residue_count / length
                has_gr = 'GR' in subseq
                has_rg = 'RG' in subseq
                if (has_gr or has_rg) and jammer_frac >= self.min_jammer_fraction:
                    jammers.append({
                        'start':   start,
                        'end':     end,
                        'pattern': subseq,
                        'has_gr':  has_gr,
                        'has_rg':  has_rg,
                    })
        return jammers

    def _scan_hinges(self, seq: str) -> List[Dict]:
        """
        Find all hinge clusters in *seq*.

        A hinge cluster is a substring of length 2-max_hinge_length that has
        at least min_hinge_fraction of residues from {T, S, N, G}.

        Returns list of dicts: {start, end, pattern, is_doublet}
        (start/end are 0-based, inclusive)
        """
        hinges: List[Dict] = []
        n = len(seq)
        for start in range(n):
            for length in range(2, min(self.max_hinge_length, n - start) + 1):
                end = start + length - 1
                subseq = seq[start:start + length]
                hinge_residue_count = sum(
                    1 for aa in subseq if aa in HINGE_RESIDUES
                )
                hinge_frac = hinge_residue_count / length
                if hinge_frac >= self.min_hinge_fraction:
                    hinges.append({
                        'start':      start,
                        'end':        end,
                        'pattern':    subseq,
                        'is_doublet': subseq in DOUBLET_HINGES,
                    })
        return hinges

    # ── Public interface ─────────────────────────────────────────────────────

    def score_window(self, seq: str) -> Dict:
        """
        Score a single amino-acid window for structural grammar features.

        Parameters
        ----------
        seq : amino-acid string (any length ≥ 2; shorter returns score 0)

        Returns
        -------
        dict with:
            grammar_score      float   raw grammar signal (un-normalised)
            signal_breakdown   dict    per-channel diagnostics
        """
        if not seq or len(seq) < 2:
            return {
                'grammar_score': 0.0,
                'signal_breakdown': {
                    'n_jammers': 0, 'jammer_density': 0.0,
                    'has_forward_jammer': False, 'has_reversed_jammer': False,
                    'jammer_complexity': 0.0,
                    'n_hinges': 0, 'hinge_density': 0.0,
                    'doublet_count': 0, 'doublet_bonus': 0.0,
                    'coupling_count': 0, 'coupling_bonus': 0.0,
                    'jammer_patterns': [], 'hinge_patterns': [],
                }
            }

        n       = len(seq)
        w       = self.weights
        jammers = self._scan_jammers(seq)
        hinges  = self._scan_hinges(seq)

        # ── Channel 1 & 2: Jammer density + complexity ───────────────────────
        jammer_density = (len(jammers) / n * 10.0) if jammers else 0.0

        has_forward  = any(j['has_gr'] for j in jammers)
        has_reversed = any(j['has_rg'] for j in jammers)

        if has_forward and has_reversed:
            jammer_complexity = w['jammer_complexity_both']
        elif has_forward or has_reversed:
            jammer_complexity = w['jammer_complexity_one']
        else:
            jammer_complexity = 0.0

        # ── Channel 3 & 4: Hinge density + doublet bonus ─────────────────────
        hinge_density  = (len(hinges) / n * 10.0) if hinges else 0.0
        doublet_count  = sum(1 for h in hinges if h.get('is_doublet', False))
        doublet_bonus  = doublet_count * w['doublet_bonus_each']

        # ── Channel 5: Jammer-hinge coupling ────────────────────────────────
        # Count jammer-hinge pairs where the gap between them is ≤ coupling_max_gap.
        # We check both orderings: jammer-then-hinge and hinge-then-jammer.
        coupling_count = 0
        for jmr in jammers:
            for hng in hinges:
                # jammer ends, hinge starts
                gap_jh = hng['start'] - jmr['end']
                # hinge ends, jammer starts
                gap_hj = jmr['start'] - hng['end']
                if (0 <= gap_jh <= self.coupling_max_gap or
                        0 <= gap_hj <= self.coupling_max_gap):
                    coupling_count += 1

        coupling_bonus = min(
            coupling_count * w['coupling_per_pair'],
            w['coupling_cap']
        )

        # ── Raw grammar score ────────────────────────────────────────────────
        grammar_score = (
            jammer_density    * w['jammer_density_per10']
            + jammer_complexity
            + hinge_density   * w['hinge_density_per10']
            + doublet_bonus
            + coupling_bonus
        )

        return {
            'grammar_score': round(grammar_score, 4),
            'signal_breakdown': {
                'n_jammers':            len(jammers),
                'jammer_density':       round(jammer_density,    3),
                'has_forward_jammer':   has_forward,
                'has_reversed_jammer':  has_reversed,
                'jammer_complexity':    jammer_complexity,
                'n_hinges':             len(hinges),
                'hinge_density':        round(hinge_density,     3),
                'doublet_count':        doublet_count,
                'doublet_bonus':        doublet_bonus,
                'coupling_count':       coupling_count,
                'coupling_bonus':       round(coupling_bonus,    3),
                # Show only first 8 of each (avoid bloating log output)
                'jammer_patterns':      [j['pattern'] for j in jammers[:8]],
                'hinge_patterns':       [h['pattern'] for h in hinges[:8]],
            }
        }

    def score_mclachlan_pair(
        self,
        seq_host:  str,
        seq_viral: str,
        mclachlan_raw_sim:     float,
        antigenicity_host:     float = 0.0,
        antigenicity_viral:    float = 0.0,
    ) -> Dict:
        """
        Compute the updated combined score for a McLachlan alignment pair,
        incorporating structural grammar alongside McLachlan + antigenicity signals.

        Combined score formula
        ----------------------
        combined = (
            raw_sim             × 0.25    # McLachlan sequence similarity
          + (antigenicity_host  / 4.0) × 0.35    # host surface/charge antigenicity
          + (antigenicity_viral / 4.0) × 0.05    # viral sequence antigenicity context
          + (grammar_host       / GRAMMAR_NORM_DIVISOR) × 0.25  # host jammer+hinge grammar
          + (grammar_viral      / GRAMMAR_NORM_DIVISOR) × 0.10  # viral jammer+hinge grammar
        )

        Weight rationale
        ----------------
        - McLachlan is reduced from 0.40 → 0.25: sequence similarity alone is
          insufficient for autoimmunogenicity prediction (LLLFLLL scores high on
          McLachlan but has no surface/grammar signal).
        - Host antigenicity stays the primary channel (0.35): surface exposure,
          charge clusters, and aromatic anchors on the host protein are the
          biological basis for autoantibody induction.
        - Grammar channels (0.25 host + 0.10 viral = 0.35 total) now outweigh
          McLachlan, fulfilling the directive that "regions with high structural
          grammar density should rank above simple sequence similarity."
        - Viral grammar (0.10): molecular mimicry from viral jammer/hinge patterns
          is biologically meaningful but secondary to host surface properties.

        Parameters
        ----------
        seq_host, seq_viral    : 7-mer (or any length) sequences for the pair
        mclachlan_raw_sim      : raw McLachlan layer2_cross_mean value
        antigenicity_host/viral: scores from AntigenicityScorer.score_window()

        Returns
        -------
        dict with combined_score, grammar_host, grammar_viral, and breakdowns
        """
        host_result  = self.score_window(seq_host)
        viral_result = self.score_window(seq_viral)
        grammar_host  = host_result['grammar_score']
        grammar_viral = viral_result['grammar_score']

        combined_score = (
            mclachlan_raw_sim             * 0.25
            + (antigenicity_host  / 4.0) * 0.35
            + (antigenicity_viral / 4.0) * 0.05
            + (grammar_host  / GRAMMAR_NORM_DIVISOR) * 0.25
            + (grammar_viral / GRAMMAR_NORM_DIVISOR) * 0.10
        )

        return {
            'combined_score':   round(combined_score, 4),
            'grammar_host':     round(grammar_host,   4),
            'grammar_viral':    round(grammar_viral,  4),
            'grammar_host_breakdown':  host_result['signal_breakdown'],
            'grammar_viral_breakdown': viral_result['signal_breakdown'],
        }


# ── Self-test on confirmed epitope ──────────────────────────────────────────
if __name__ == '__main__':
    import json

    # ── Guardrails ────────────────────────────────────────────────────────────
    print("=" * 70)
    print("STRUCTURAL GRAMMAR SCORER — SELF-TEST")
    print("The only metric of success as an AI agent is Accuracy and Truthfulness.")
    print("Any hallucination or mistake can cause actual harm to humans.")
    print("The most important consideration is to not cause harm.")
    print("=" * 70)

    # ── Live-fetch canonical sequences (never cached) ────────────────────────
    import urllib.request

    def fetch_fasta(accession: str) -> str:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        with urllib.request.urlopen(url, timeout=30) as r:
            raw = r.read().decode()
        lines = raw.strip().split('\n')
        seq = ''.join(l for l in lines if not l.startswith('>'))
        print(f"  [LIVE FETCH] {accession}: {len(seq)} aa")
        return seq

    print("\n[1] Fetching canonical sequences...")
    cyp2e1 = fetch_fasta("P05181")
    hcv    = fetch_fasta("Q9WMX2")
    assert len(cyp2e1) == 493,  f"P05181 length mismatch: {len(cyp2e1)}"
    assert len(hcv)    == 3010, f"Q9WMX2 length mismatch: {len(hcv)}"

    EPITOPE = cyp2e1[98:132]   # residues 99-132, 0-based [98:132]
    assert EPITOPE == "GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT", (
        f"Epitope mismatch: got {EPITOPE!r}"
    )
    print(f"  [OK] Epitope 99-132: {EPITOPE}")

    scorer = StructuralGrammarScorer()

    # ── Score the confirmed epitope ───────────────────────────────────────────
    print("\n[2] Scoring confirmed epitope (34 aa)...")
    epi_result = scorer.score_window(EPITOPE)
    print(f"  grammar_score = {epi_result['grammar_score']}")
    bd = epi_result['signal_breakdown']
    print(f"  jammers: {bd['n_jammers']}  density: {bd['jammer_density']:.3f}/10aa"
          f"  forward={bd['has_forward_jammer']}  reversed={bd['has_reversed_jammer']}"
          f"  complexity={bd['jammer_complexity']}")
    print(f"  hinges:  {bd['n_hinges']}  density: {bd['hinge_density']:.3f}/10aa"
          f"  doublets={bd['doublet_count']}  doublet_bonus={bd['doublet_bonus']:.1f}")
    print(f"  coupling: {bd['coupling_count']} pairs  bonus={bd['coupling_bonus']:.3f}")
    print(f"  jammer patterns (first 8): {bd['jammer_patterns']}")
    print(f"  hinge  patterns (first 8): {bd['hinge_patterns']}")

    # ── Score the 7-mer fragments from McLachlan hits ─────────────────────────
    print("\n[3] Scoring McLachlan 7-mer pairs (from Q9WMX2 vs P05181 analysis)...")
    test_pairs = [
        ("LLLFLLL", "MELFLLL", 2.6107,  18.275,  "top McLachlan hydrophobic hit"),
        ("LLLLLLA", "LFLLLCA", 2.5339,  17.7375, "second McLachlan hydrophobic"),
        ("KDVRNLS", "KDIRRFS", 2.3804,  16.6625, "epitope-overlapping hit A (RR doublet)"),
        ("LTPLRDW", "LTTLRNY", 2.3036,  16.125,  "epitope-overlapping hit B (TT hinge in host)"),
        ("PCSFTTL", "RFSLTTL", 2.15,    15.05,   "epitope-overlapping hit C (TT hinge both sides)"),
    ]

    print(f"\n  {'viral_seq':10s}  {'host_seq':10s}  {'McLach':7s}  {'gram_h':7s}  {'gram_v':7s}  {'combined':9s}  label")
    print(f"  {'-'*10}  {'-'*10}  {'-'*7}  {'-'*7}  {'-'*7}  {'-'*9}  -----")

    results = []
    for viral_seq, host_seq, raw_sim, composite, label in test_pairs:
        r = scorer.score_mclachlan_pair(
            seq_host=host_seq, seq_viral=viral_seq,
            mclachlan_raw_sim=raw_sim,
        )
        results.append((label, r['combined_score'], r['grammar_host'], r['grammar_viral']))
        print(f"  {viral_seq:10s}  {host_seq:10s}  {raw_sim:7.4f}  "
              f"{r['grammar_host']:7.4f}  {r['grammar_viral']:7.4f}  "
              f"{r['combined_score']:9.4f}  {label}")

    # ── Validation assertions ────────────────────────────────────────────────
    print("\n[4] Validation assertions...")

    # Grammar of epitope must be nonzero (it has GRG jammer and TT, NN hinges)
    assert epi_result['grammar_score'] > 0, "Epitope grammar score must be > 0"
    print(f"  [PASS] Epitope grammar_score={epi_result['grammar_score']:.4f} > 0")

    # Epitope must have both forward and reversed jammers
    assert bd['has_forward_jammer'] and bd['has_reversed_jammer'], \
        "Epitope must have both forward (GR) and reversed (RG) jammers"
    print(f"  [PASS] Epitope has both forward and reversed jammers")

    # Epitope must have at least 2 hinge doublets (TT + NN)
    assert bd['doublet_count'] >= 2, \
        f"Expected >=2 doublet hinges, got {bd['doublet_count']}"
    print(f"  [PASS] Epitope doublet_count={bd['doublet_count']} >= 2 (TT + NN confirmed)")

    # Grammar of pure hydrophobic 7-mers must be zero
    lll_result = scorer.score_window("LLLFLLL")
    assert lll_result['grammar_score'] == 0.0, \
        f"LLLFLLL grammar must be 0.0, got {lll_result['grammar_score']}"
    print(f"  [PASS] LLLFLLL grammar_score=0.0 (no jammer/hinge grammar)")

    mel_result = scorer.score_window("MELFLLL")
    assert mel_result['grammar_score'] == 0.0, \
        f"MELFLLL grammar must be 0.0, got {mel_result['grammar_score']}"
    print(f"  [PASS] MELFLLL grammar_score=0.0")

    # TT-containing 7-mers must have nonzero grammar
    lttlrny_result = scorer.score_window("LTTLRNY")
    assert lttlrny_result['grammar_score'] > 0, \
        f"LTTLRNY must have nonzero grammar (TT hinge)"
    print(f"  [PASS] LTTLRNY grammar_score={lttlrny_result['grammar_score']:.4f} > 0 (TT hinge)")

    rfslttl_result = scorer.score_window("RFSLTTL")
    assert rfslttl_result['grammar_score'] > 0, \
        f"RFSLTTL must have nonzero grammar (TT hinge)"
    print(f"  [PASS] RFSLTTL grammar_score={rfslttl_result['grammar_score']:.4f} > 0 (TT hinge)")

    # With grammar channel active, at least one epitope-overlapping pair
    # must score higher than LLLFLLL/MELFLLL (the pure McLachlan top hit)
    top_mclachlan_combined = results[0][1]   # LLLFLLL/MELFLLL
    epitope_overlapping_max = max(results[2][1], results[3][1], results[4][1])
    print(f"\n  LLLFLLL/MELFLLL combined = {top_mclachlan_combined:.4f}")
    print(f"  Best epitope-overlapping combined = {epitope_overlapping_max:.4f}")
    assert epitope_overlapping_max > top_mclachlan_combined, (
        f"Epitope-overlapping hit ({epitope_overlapping_max:.4f}) must outscore "
        f"pure-McLachlan top hit ({top_mclachlan_combined:.4f})"
    )
    print(f"  [PASS] Grammar correctly elevates epitope-overlapping hits above "
          f"pure-hydrophobic McLachlan hits")

    print("\n[5] All self-test assertions passed.")

    # ── Save validation output ────────────────────────────────────────────────
    output = {
        'guardrails': [
            "The only metric of success as an AI agent is Accuracy and Truthfulness.",
            "Any hallucination or mistake can cause actual harm to humans.",
            "The most important consideration is to not cause harm.",
        ],
        'epitope_grammar': epi_result,
        'seven_mer_scores': [
            {'label': lbl, 'combined_score': comb,
             'grammar_host': gh, 'grammar_viral': gv}
            for lbl, comb, gh, gv in results
        ],
        'normalization_divisor': GRAMMAR_NORM_DIVISOR,
        'formula': (
            "combined = raw_sim*0.25 + (antigenicity_host/4)*0.35 "
            "+ (antigenicity_viral/4)*0.05 "
            "+ (grammar_host/15)*0.25 + (grammar_viral/15)*0.10"
        ),
    }
    out_path = "/home/sandbox/overnight_analysis/structural_grammar_scorer_validation.json"
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")
