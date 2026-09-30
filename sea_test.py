"""
sea_test.py
===========
Synthetic regression tests for the SEA Scanner module.

Each test case is hand-crafted so that its expected architecture class is
unambiguous given the SEA scoring rules.  The homologous region is always
'FLKEKFKEQ' (contains a KFERQ-like 5-mer: KFKEQ → Q·F·K·K·E — wait, let's
pick a clean one).

KFERQ biochemical rule (from SEAConfig):
  – starts with Q or N (glutamine or asparagine)
  – contains exactly one F or Y (aromatic)
  – contains at least one K or R (basic)
  – contains at least one L, I, or V (hydrophobic)
  – 5 residues

Clean KFERQ-like 5-mer used here: 'QFVLK'
  Q=start, F=aromatic, V=hydrophobic, L=hydrophobic, K=basic  ✓

Sequence layout convention (window = 15 residues each side):
  [N-wing 15aa] [homologous 9aa] [C-wing 15aa]

We embed hinges/jammers/motifs inside the wings to test each scenario.
"""

import sys
import os
sys.path.insert(0, '/home/sandbox')

from sea_module import SEAModule, SEAConfig
from dataclasses import dataclass
from typing import List

# ---------------------------------------------------------------------------
# Helper: build a 15-residue wing with an optional embedded feature
# ---------------------------------------------------------------------------
def _pad(seq: str, length: int = 15, fill: str = 'A') -> str:
    """Pad/trim seq to exactly `length` with fill character."""
    if len(seq) >= length:
        return seq[:length]
    return seq + fill * (length - len(seq))

HOMO = 'ALVMTLRPQ'   # 9-residue homologous core (neutral, no features)

# KFERQ-like motif — 5 residues: Q(start) + F(aromatic) + L(hydrophobic) + V(hydrophobic) + K(basic)
KFERQ = 'QFLVK'

# Tier-1 hinge: pure SS cluster  →  'SS' somewhere in the wing
T1_HINGE = 'AAASSAAAAAAAAAA'   # 'SS' at pos 3-4, padded to 15

# Tier-2 hinge: ST/TS/TT  →  'ST' in wing
T2_HINGE = 'AAASTAAAAAAAAAA'   # 'ST' at pos 3-4

# Jammer: GR repeat × 3  →  'GRGRGR'
JAMMER_3 = 'AAGRGRGR' + 'A' * 7  # 8+7=15

# Neutral wing (no features)
NEUTRAL = 'A' * 15


# ---------------------------------------------------------------------------
# Build sequence pairs
# ---------------------------------------------------------------------------
# Protein 2 (host) is kept neutral throughout — we don't want it to produce
# features that bleed into protein 1 test logic.
P2_SEQ = NEUTRAL + HOMO + NEUTRAL


def make_pair(n_wing: str, c_wing: str) -> str:
    """Return viral protein1 sequence with embedded wings."""
    return _pad(n_wing) + HOMO + _pad(c_wing)


# ── Case definitions ────────────────────────────────────────────────────────
# Each entry: (label, p1_seq, expected_arch_class)
# homologous_start = 15  (0-based index, length of n_wing)

HOMO_START = 15
HOMO_LEN   = len(HOMO)   # 9

# Provided degradation motifs: we supply the KFERQ motif pre-embedded in the
# N-wing for cases that need it, but also pass it in the provided list so the
# DegradationMotifScanner always finds it even if inline scan misses edge cases.

# For cases that need a KFERQ in wing, we embed it at the very start of the
# wing so it's within the 15-residue window.
# N-wing with KFERQ at positions 0-4:  QFLVK + 10 neutrals
N_KFERQ = _pad('QFLVK', 15)   # 'QFLVKAAAAAAAAA'
C_KFERQ = _pad('QFLVK', 15)


CASES = [
    # ── 1. SUPER_EPITOPE ──────────────────────────────────────────────────
    # N-wing: T1 hinge (SS) + KFERQ embedded;  C-wing: jammer (GRGRGRGR)
    # Sandwiched (hinge N-side, jammer C-side) + degradation → SUPER_EPITOPE
    {
        'label':    'SUPER_EPITOPE',
        'p1_seq':   _pad('AASSQFLVK', 15) + HOMO + _pad('GRGRGRGR', 15),
        'p2_seq':   P2_SEQ,
        # provided_motifs uses keys: motif, position, protein
        'provided_motifs': [{'motif': 'QFLVK', 'position': 4, 'protein': 'protein1'}],
        'expected': 'SUPER_EPITOPE',
    },

    # ── 2. SANDWICHED (no degradation in window) ─────────────────────────
    # N-wing: T1 hinge only;  C-wing: jammer only;  no KFERQ anywhere
    {
        'label':    'SANDWICHED',
        'p1_seq':   make_pair(T1_HINGE, JAMMER_3),
        'p2_seq':   P2_SEQ,
        'provided_motifs': [],
        'expected': 'SANDWICHED',
    },

    # ── 3. HINGE_AND_DEGRADATION ─────────────────────────────────────────
    # N-wing: T1 hinge + KFERQ;  C-wing: neutral;  no C-side feature
    {
        'label':    'HINGE_AND_DEGRADATION',
        'p1_seq':   _pad('AASSQFLVK', 15) + HOMO + NEUTRAL,
        'p2_seq':   P2_SEQ,
        'provided_motifs': [{'motif': 'QFLVK', 'position': 4, 'protein': 'protein1'}],
        'expected': 'HINGE_AND_DEGRADATION',
    },

    # ── 4. JAMMER_AND_DEGRADATION ────────────────────────────────────────
    # N-wing: neutral;  C-wing: jammer + KFERQ;  no N-side feature
    {
        'label':    'JAMMER_AND_DEGRADATION',
        'p1_seq':   NEUTRAL + HOMO + _pad('GRGRQFLVK', 15),
        'p2_seq':   P2_SEQ,
        'provided_motifs': [{'motif': 'QFLVK', 'position': 15 + len(HOMO) + 4,
                              'protein': 'protein1'}],
        'expected': 'JAMMER_AND_DEGRADATION',
    },

    # ── 5. HINGE_ONLY ────────────────────────────────────────────────────
    # N-wing: T2 hinge (ST);  C-wing: neutral;  no KFERQ
    {
        'label':    'HINGE_ONLY',
        'p1_seq':   make_pair(T2_HINGE, NEUTRAL),
        'p2_seq':   P2_SEQ,
        'provided_motifs': [],
        'expected': 'HINGE_ONLY',
    },

    # ── 6. JAMMER_ONLY ───────────────────────────────────────────────────
    # N-wing: neutral;  C-wing: jammer (GRGRGR);  no KFERQ
    {
        'label':    'JAMMER_ONLY',
        'p1_seq':   make_pair(NEUTRAL, JAMMER_3),
        'p2_seq':   P2_SEQ,
        'provided_motifs': [],
        'expected': 'JAMMER_ONLY',
    },

    # ── 7. DEGRADATION_ONLY ──────────────────────────────────────────────
    # N-wing: KFERQ at start;  C-wing: neutral;  no hinge or jammer
    {
        'label':    'DEGRADATION_ONLY',
        'p1_seq':   make_pair(N_KFERQ, NEUTRAL),
        'p2_seq':   P2_SEQ,
        'provided_motifs': [{'motif': 'QFLVK', 'position': 0, 'protein': 'protein1'}],
        'expected': 'DEGRADATION_ONLY',
    },

    # ── 8. NONE ───────────────────────────────────────────────────────────
    # All neutral — no hinge, no jammer, no KFERQ
    {
        'label':    'NONE',
        'p1_seq':   make_pair(NEUTRAL, NEUTRAL),
        'p2_seq':   P2_SEQ,
        'provided_motifs': [],
        'expected': 'NONE',
    },
]


# ---------------------------------------------------------------------------
# Build the homologous_pairs input format expected by SEAModule
# ---------------------------------------------------------------------------
def build_pairs(case: dict, rank: int = 1) -> List[dict]:
    return [
        {
            'seq1':             HOMO,            # homologous region in viral protein
            'seq2':             HOMO,            # homologous region in host protein
            'position1':        HOMO_START,      # absolute start in protein1
            'position2':        HOMO_START,      # absolute start in protein2
            'protein1':         case['p1_seq'],  # FULL viral protein (scanned for features)
            'similarity_score': 0.85,
            'rank':             rank,
        }
    ]


# ---------------------------------------------------------------------------
# Run all cases
# ---------------------------------------------------------------------------
def run_tests():
    print("=" * 70)
    print("SEA Scanner — Architecture Class Regression Tests")
    print("=" * 70)

    passed = 0
    failed = 0
    results_summary = []

    for idx, case in enumerate(CASES):
        pairs  = build_pairs(case, rank=idx+1)
        module = SEAModule(
            virus_name = 'EBV',      # use EBV for max APC weight → easier to see score differences
            config     = SEAConfig(),
        )
        results = module.run(
            homologous_pairs    = pairs,
            autoimmune_diseases = ['SLE', 'MS', 'RA'],
            degradation_motifs  = case['provided_motifs'],
        )

        assert len(results) == 1, f"Expected 1 result, got {len(results)}"
        r = results[0]
        got      = r.architecture_class.name
        expected = case['expected']
        ok       = (got == expected)

        status = 'PASS' if ok else 'FAIL'
        if ok:
            passed += 1
        else:
            failed += 1

        line = (
            f"[{status}]  {case['label']:<28}  "
            f"expected={expected:<28}  got={got:<28}  "
            f"score={r.final_sea_score:.4f}"
        )
        print(line)
        results_summary.append({
            'label':    case['label'],
            'expected': expected,
            'got':      got,
            'score':    r.final_sea_score,
            'ok':       ok,
            'result':   r,
        })

    print("-" * 70)
    print(f"Results: {passed} passed, {failed} failed out of {len(CASES)} cases")
    print()

    # ── Score ordering check ─────────────────────────────────────────────────
    # SUPER_EPITOPE score must exceed SANDWICHED, which must exceed HINGE_AND_DEG
    # etc.  If architecture classes matched but scores are inverted, flag it.
    score_map = {r['label']: r['score'] for r in results_summary}

    ordering_checks = [
        ('SUPER_EPITOPE',          'SANDWICHED'),
        ('SANDWICHED',             'HINGE_AND_DEGRADATION'),
        ('HINGE_AND_DEGRADATION',  'HINGE_ONLY'),
        ('HINGE_AND_DEGRADATION',  'DEGRADATION_ONLY'),
        ('HINGE_ONLY',             'NONE'),
        ('JAMMER_ONLY',            'NONE'),
        ('DEGRADATION_ONLY',       'NONE'),
    ]
    print("Score ordering assertions:")
    order_ok = True
    for higher, lower in ordering_checks:
        h, l = score_map.get(higher, 0), score_map.get(lower, 0)
        ok = h > l
        sym = '✓' if ok else '✗'
        print(f"  {sym}  {higher} ({h:.4f}) > {lower} ({l:.4f})")
        if not ok:
            order_ok = False

    print()
    if order_ok:
        print("All score ordering assertions passed.")
    else:
        print("WARNING: one or more score ordering assertions FAILED.")

    # ── Detailed diagnostics for SUPER_EPITOPE ──────────────────────────────
    print()
    print("── SUPER_EPITOPE detail ─────────────────────────────────────────")
    sea = [r for r in results_summary if r['label'] == 'SUPER_EPITOPE'][0]['result']
    print(f"  hinges detected    : {len(sea.hinges)}")
    for h in sea.hinges:
        print(f"    {h}")
    print(f"  jammers detected   : {len(sea.jammers)}")
    for j in sea.jammers:
        print(f"    {j}")
    print(f"  deg motifs detected: {len(sea.degradation_motifs)}")
    for d in sea.degradation_motifs:
        print(f"    {d}")
    print(f"  sandwiched         : {sea.is_sandwiched}")
    print(f"  arch class         : {sea.architecture_class.name}")
    print(f"  base_score         : {sea.base_score:.4f}")
    print(f"  hinge_score        : {sea.hinge_score:.4f}")
    print(f"  jammer_score       : {sea.jammer_score:.4f}")
    print(f"  jammer_density_scr : {sea.jammer_density_score:.4f}")
    print(f"  deg_score          : {sea.degradation_score:.4f}")
    print(f"  arch_bonus         : {sea.architecture_bonus:.4f}")
    print(f"  cell_weight        : {sea.cell_type_weight:.4f}")
    print(f"  FINAL sea score    : {sea.final_sea_score:.4f}")

    print()
    print("── APC tropism check ────────────────────────────────────────────")
    # Run the NONE case again with a default (unknown) virus → weight should be 1.0
    case_none = [c for c in CASES if c['label'] == 'NONE'][0]
    pairs_none = build_pairs(case_none, rank=99)
    module_default = SEAModule(
        virus_name = 'unknown_virus_xyz',
        config     = SEAConfig(),
    )
    r_default = module_default.run(
        homologous_pairs    = pairs_none,
        autoimmune_diseases = [],
        degradation_motifs  = [],
    )[0]
    ebv_none = [r for r in results_summary if r['label'] == 'NONE'][0]
    print(f"  NONE with EBV    (weight=2.0): score={ebv_none['score']:.4f}")
    print(f"  NONE with unknown(weight=1.0): score={r_default.final_sea_score:.4f}")
    ratio = ebv_none['score'] / r_default.final_sea_score if r_default.final_sea_score else float('nan')
    print(f"  Ratio (should be ~2.0)       : {ratio:.4f}")

    print()
    print("=" * 70)
    return failed == 0 and order_ok


# ---------------------------------------------------------------------------
# EBNA density amplification test
# ---------------------------------------------------------------------------
def run_ebna_density_test():
    """
    Verify that a protein with an EBNA-style high-density jammer region
    produces a substantially higher jammer_density_score than one with a
    single isolated GA dipeptide — demonstrating the action-potential
    amplitude analogy: nearby copies amplify the signal.

    Test layout (both cases):
      N-wing (15aa): neutral
      Homologous (9aa): ALVMTLRPQ  (neutral HOMO)
      C-wing (15aa): varies

    Case A — single GA:
      C-wing: 'GA' + 13 neutral 'A's → exactly one GA dipeptide in range

    Case B — EBNA-style 30-residue run (15 in C-wing + overflow):
      C-wing: 'GAGAGAGAGAGAGAG' (15aa of alternating GA) — densely packed;
      additional GA copies continue just past the window boundary so the
      density scan (window=30) picks up many hits.

    Assertion:
      ebna_density  > single_density
      amplification > 5.0×  (expected ~20× based on design math)
    """
    print()
    print("── EBNA Density Amplification Test ─────────────────────────────")

    cfg = SEAConfig()

    # ── Case A: single GA ────────────────────────────────────────────────
    p1_single = NEUTRAL + HOMO + _pad('GA', 15)   # 'GA' + 13 'A's in C-wing
    pairs_single = [{
        'seq1': HOMO, 'seq2': HOMO,
        'position1': HOMO_START, 'position2': HOMO_START,
        'protein1': p1_single,
        'similarity_score': 0.85, 'rank': 1,
    }]
    module_single = SEAModule(virus_name='unknown_virus_xyz', config=cfg)
    r_single = module_single.run(
        homologous_pairs=pairs_single,
        autoimmune_diseases=[],
        degradation_motifs=[],
    )[0]
    single_density = r_single.jammer_density_score

    # ── Case B: EBNA-style dense GA repeat ──────────────────────────────
    # Fill the entire C-wing (15aa) with alternating GA and extend beyond
    # by making the full protein tail a 30-residue GAGAGA... run so the
    # density scan window (30 each side) is saturated.
    ga_repeat_30 = 'GA' * 15   # 30 residues of GAGAGA...
    p1_ebna = NEUTRAL + HOMO + ga_repeat_30   # C-wing extends 30aa
    pairs_ebna = [{
        'seq1': HOMO, 'seq2': HOMO,
        'position1': HOMO_START, 'position2': HOMO_START,
        'protein1': p1_ebna,
        'similarity_score': 0.85, 'rank': 1,
    }]
    module_ebna = SEAModule(virus_name='unknown_virus_xyz', config=cfg)
    r_ebna = module_ebna.run(
        homologous_pairs=pairs_ebna,
        autoimmune_diseases=[],
        degradation_motifs=[],
    )[0]
    ebna_density = r_ebna.jammer_density_score

    # ── Report ───────────────────────────────────────────────────────────
    amplification = ebna_density / single_density if single_density > 0 else float('inf')
    print(f"  Single GA density score  : {single_density:.4f}")
    print(f"  EBNA 30-aa repeat score  : {ebna_density:.4f}")
    print(f"  Amplification ratio      : {amplification:.2f}×")

    assert ebna_density > single_density, (
        f"EBNA density ({ebna_density:.4f}) should exceed single GA ({single_density:.4f})"
    )
    assert amplification > 5.0, (
        f"Amplification ratio {amplification:.2f}× should exceed 5× minimum"
    )
    print("  [PASS]  EBNA density > single GA  ✓")
    print(f"  [PASS]  Amplification > 5×  ({amplification:.2f}×)  ✓")
    print()
    return True


if __name__ == '__main__':
    success = run_tests()
    ebna_ok = run_ebna_density_test()
    sys.exit(0 if (success and ebna_ok) else 1)
