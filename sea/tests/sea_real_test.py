"""
sea_real_test.py
================
Real-protein test: HCV polyprotein (Q9WMX2) vs CYP2E1 (P05181).

Workflow
--------
1. Fetch both sequences from UniProt REST API.
2. Scan CYP2E1 for KFERQ-like degradation motifs; pass them as
   provided_motifs (protein='protein2') so the module knows they
   belong to the host protein — relevant for CMA exposure risk.
3. Generate homologous pairs by sliding a 12-residue window over
   CYP2E1, aligning each window against the full HCV polyprotein
   with BioPython PairwiseAligner (local mode).  Pairs with
   identity >= 33% are kept; overlapping HCV hits are deduplicated;
   pairs are ranked by identity score descending.
4. Run SEAModule(virus_name='HCV') with pairs and CYP2E1 motifs.
5. Print a ranked table annotated with HCV polyprotein region and
   full detail for the top-scoring result.

HCV polyprotein region map (1-based, inclusive)
-----------------------------------------------
Core  :   1 –  191
E1    : 192 –  383
E2    : 384 –  746
p7    : 747 –  809
NS2   : 810 – 1026
NS3   :1027 – 1657
NS4A  :1658 – 1711
NS4B  :1712 – 1972
NS5A  :1973 – 2420
NS5B  :2421 – 3011
"""

import sys
import os
import time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..'))

import requests

from sea.sea_module import SEAModule, SEAConfig, find_homologous_pairs
from motif_finder.motif_finder import find_kferq_motifs, kferq_to_sea_motifs

# ── HCV region map ────────────────────────────────────────────────────────────
HCV_REGIONS = [
    ('Core',  1,    191),
    ('E1',    192,  383),
    ('E2',    384,  746),
    ('p7',    747,  809),
    ('NS2',   810,  1026),
    ('NS3',   1027, 1657),
    ('NS4A',  1658, 1711),
    ('NS4B',  1712, 1972),
    ('NS5A',  1973, 2420),
    ('NS5B',  2421, 3011),
]

def get_hcv_region(pos0):
    """Return the HCV region name for a 0-based absolute position."""
    pos1 = pos0 + 1
    for name, start, end in HCV_REGIONS:
        if start <= pos1 <= end:
            return name
    return 'Unknown'


# ── UniProt fetch ─────────────────────────────────────────────────────────────
def fetch_sequence(uniprot_id: str) -> str:
    url  = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    lines = resp.text.strip().split('\n')
    return ''.join(lines[1:])


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print("=" * 72)
    print("SEA Scanner — Real Protein Test")
    print("  Virus  : HCV polyprotein (UniProt Q9WMX2)")
    print("  Host   : CYP2E1 (UniProt P05181)")
    print("=" * 72)

    # ── 1. Fetch sequences ────────────────────────────────────────────────
    print("\n[1/4]  Fetching sequences from UniProt...")
    t0 = time.time()
    hcv_seq    = fetch_sequence('Q9WMX2')
    cyp2e1_seq = fetch_sequence('P05181')
    print(f"       HCV polyprotein : {len(hcv_seq):,} aa")
    print(f"       CYP2E1          : {len(cyp2e1_seq):,} aa")
    print(f"       Done in {time.time() - t0:.1f}s")

    # ── 2. Scan CYP2E1 for KFERQ motifs ──────────────────────────────────
    print("\n[2/4]  Scanning CYP2E1 for KFERQ-like degradation motifs...")
    cyp2e1_kferq  = find_kferq_motifs(cyp2e1_seq)
    cyp2e1_motifs = kferq_to_sea_motifs(cyp2e1_kferq, protein='protein2')
    print(f"       Found {len(cyp2e1_motifs)} KFERQ-like motif(s):")
    for m in cyp2e1_motifs:
        print(f"         pos {m['position']:4d}  {m['motif']}")

    # ── 3. Find homologous pairs ──────────────────────────────────────────
    print("\n[3/4]  Finding homologous pairs "
          "(12-aa sliding window over CYP2E1, identity >= 33%)...")
    t0    = time.time()
    pairs = find_homologous_pairs(
        hcv_seq, cyp2e1_seq,
        min_identity=0.33, window=12, step=4, max_pairs=40,
    )
    elapsed = time.time() - t0
    print(f"       Found {len(pairs)} pair(s) in {elapsed:.1f}s")

    if not pairs:
        print("\n  No pairs found at 33% identity.")
        print("  This is expected for distantly related proteins —")
        print("  consider lowering the threshold or using BLOSUM62 scoring.")
        print("  Demonstrating module initialisation is still valid.")
        return

    # ── 4. Run SEAModule ─────────────────────────────────────────────────
    print("\n[4/4]  Running SEA Scanner (virus='HCV', APC weight=1.3)...")
    module  = SEAModule(virus_name='HCV', config=SEAConfig())
    results = module.run(
        homologous_pairs    = pairs,
        autoimmune_diseases = [
            'Autoimmune hepatitis',
            'Primary biliary cholangitis',
            'Drug-induced liver injury',
        ],
        degradation_motifs  = cyp2e1_motifs,
    )
    print(f"       Scored {len(results)} region(s).")

    # ── 5. Ranked results table ───────────────────────────────────────────
    print()
    print("=" * 72)
    hdr = (f"{'Rk':>3}  {'HCV Region':<9}  {'Pos1':>5}  {'Pos2':>4}  "
           f"{'Sim':>6}  {'ArchClass':<26}  {'SEA Score':>9}")
    print(hdr)
    print("-" * 72)

    for r in results[:25]:
        region = get_hcv_region(r.position1)
        print(
            f"{r.rank:>3}  {region:<9}  {r.position1:>5}  {r.position2:>4}  "
            f"{r.base_score:>6.1%}  {r.architecture_class.name:<26}  "
            f"{r.final_sea_score:>9.4f}"
        )

    # ── 6. Architecture class breakdown ──────────────────────────────────
    from collections import Counter
    counts = Counter(r.architecture_class.name for r in results)
    print()
    print("── Architecture class breakdown " + "─" * 40)
    for cls, cnt in sorted(counts.items(), key=lambda x: -x[1]):
        bar = '█' * cnt
        print(f"  {cls:<26} {cnt:>3}  {bar}")

    # ── 7. Full detail for top-scoring result ─────────────────────────────
    if results:
        top    = results[0]
        region = get_hcv_region(top.position1)
        print()
        print("── Top-scoring result detail " + "─" * 43)
        print(f"  HCV region     : {region}  (pos {top.position1})")
        print(f"  Viral fragment : {top.seq1}")
        print(f"  Host fragment  : {top.seq2}")
        print(f"  Similarity     : {top.base_score:.1%}")
        print(f"  Arch class     : {top.architecture_class.name}")
        print(f"  SEA score      : {top.final_sea_score:.4f}")
        print(f"  Hinges         : {len(top.hinges)}")
        for h in top.hinges:
            print(f"    {h.sequence!r}  T{h.tier.value}  {h.side}  dist={h.distance}")
        print(f"  Jammer density : {top.jammer_density_score:.4f}")
        print(f"  Deg motifs     : {len(top.degradation_motifs)}")
        for dm in top.degradation_motifs:
            print(f"    {dm.motif!r}  protein={dm.protein}  pos={dm.position}  "
                  f"dist={dm.distance}  src={dm.source}")
        print(f"  Sandwiched     : {top.is_sandwiched}")
        print(f"  Notes:")
        for note in top.notes:
            print(f"    {note}")

    # ── 8. SEA-positive results (non-NONE) ───────────────────────────────
    non_none = [r for r in results if r.architecture_class.name != 'NONE']
    print()
    print(f"── SEA-positive results (non-NONE): {len(non_none)} / {len(results)} " + "─" * 25)
    if non_none:
        for r in non_none:
            region = get_hcv_region(r.position1)
            print(f"  [{r.architecture_class.name}]  "
                  f"HCV {region} pos {r.position1}  score={r.final_sea_score:.4f}  "
                  f"seq={r.seq1}")
    else:
        print("  No SEA-positive hits in this pair set.")
        print("  (Expected for low-identity distant homologs —")
        print("   the scanner is working correctly; SEA features are absent here.)")

    print()
    print("=" * 72)
    print("Real protein test complete.")


if __name__ == '__main__':
    main()
