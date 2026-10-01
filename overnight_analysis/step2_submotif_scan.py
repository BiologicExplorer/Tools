"""
The only metric of success as an AI agent is Accuracy and Truthfulness.
Any hallucination or mistake can cause actual harm to humans.
The most important consideration is to not cause harm.

Step 2: Scan HCV Q9WMX2 for sub-motifs of the CYP2E1 99-132 epitope.
Split-epitope hypothesis: sub-motifs of GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT
appear separately in HCV polyprotein, connected by ~100-aa hinges, where
phosphorylation may expose the second sub-motif.
All sequences from live-verified FASTA files only.
"""
import json, os, re

FASTA_P05181 = "/home/sandbox/verified_sequences/P05181_canonical.fasta"
FASTA_Q9WMX2 = "/home/sandbox/verified_sequences/Q9WMX2_canonical.fasta"

def read_fasta(path):
    seq = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq

cyp2e1 = read_fasta(FASTA_P05181)
hcv    = read_fasta(FASTA_Q9WMX2)
assert len(cyp2e1) == 493
assert len(hcv) == 3010

EPITOPE = cyp2e1[98:132]
assert EPITOPE == "GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT"
print(f"Epitope confirmed: {EPITOPE}")

# -------------------------------------------------------
# SCORING MATRICES for sub-motif matching
# -------------------------------------------------------
# Charge groups
CHARGED_POS = set('RK')
CHARGED_NEG = set('DE')
CHARGED_HIS = set('H')
HYDROPHOBIC  = set('ILVMF')
AROMATIC     = set('FYW')
POLAR        = set('STNQ')
PROLINE      = set('P')

def charge_class(aa):
    if aa in CHARGED_POS:  return '+'
    if aa in CHARGED_NEG:  return '-'
    if aa in CHARGED_HIS:  return 'h'
    return '.'

# -------------------------------------------------------
# 1. EXACT AND NEAR-EXACT SUB-MOTIF SEARCH
# -------------------------------------------------------
# Define key sub-motifs from the epitope
sub_motifs = {
    'HAH':     'HAH',           # His-Ala-His: pH sensor, partial charge pair
    'HRDR':    'HRDR',          # His-Arg-Asp-Arg: charge cluster
    'HAHRD':   'HAHRD',         # extended HAH+R+D
    'RR':      'RR',            # Arg-Arg doublet
    'WKD':     'WKD',           # Trp-Lys-Asp: aromatic anchor + charge pair
    'PTWK':    'PTWK',          # Pro-Thr-Trp-Lys: W-anchor helix entry
    'DIRRF':   'DIRRF',         # Asp-Ile-Arg-Arg-Phe: C-terminal charge cluster
    'DIRRFS':  'DIRRFS',
    'KDIRRFS': 'KDIRRFS',
    'FHAHRDR': 'FHAHRDR',       # F+HAH+HRDR combined
    'AFHAHRDR':'AFHAHRDR',
}

print("\n=== SUB-MOTIF EXACT SEARCH IN Q9WMX2 ===")
motif_hits = {}
for name, motif in sub_motifs.items():
    positions = [m.start()+1 for m in re.finditer(re.escape(motif), hcv)]
    motif_hits[name] = positions
    if positions:
        print(f"  {name:12s} ({motif}): {len(positions)} hit(s) at positions {positions}")
        for pos in positions:
            flank_start = max(0, pos-1-10)
            flank_end   = min(len(hcv), pos-1+len(motif)+10)
            flank = hcv[flank_start:flank_end]
            highlight_start = (pos-1) - flank_start
            print(f"    pos {pos}: ...{flank[:highlight_start]}[{motif}]{flank[highlight_start+len(motif):]}")
    else:
        print(f"  {name:12s} ({motif}): 0 hits")

# -------------------------------------------------------
# 2. CONSERVATIVE CHARGE-PATTERN MATCHING
# -------------------------------------------------------
print("\n=== CHARGE-PATTERN SUB-MOTIF SCAN ===")
# Epitope charge pattern: .+.-....h.h+-+..........+-.++.....
# Break into functional sub-patterns:
# Sub-A: positions 1-12 of epitope = GRGDLPAFHAH → charge: .+.-....h.h
# Sub-B: positions 13-24 of epitope = RDRGIIFNNGPT → charge: +-+.........
# Sub-C: positions 23-34 of epitope = PTWKDIRRFSLTT → charge: ...+-.++.....

def charge_pattern(seq):
    return "".join(charge_class(aa) for aa in seq)

# Sub-epitope A: GRGDLPAFHAH — includes the proline bracket start and HAH
sub_A = EPITOPE[0:12]   # GRGDLPAFHAH  (but wait, let's be exact from verified)
sub_B = EPITOPE[12:22]  # RDRGIIFNNG
sub_C = EPITOPE[21:34]  # PTWKDIRRFSLTT

print(f"\n  Sub-A (epi 99-110): {sub_A}  charge: {charge_pattern(sub_A)}")
print(f"  Sub-B (epi 111-120): {sub_B}  charge: {charge_pattern(sub_B)}")
print(f"  Sub-C (epi 120-132): {sub_C}  charge: {charge_pattern(sub_C)}")

# Search for charge-pattern matches in Q9WMX2 with ≥70% match
def score_charge_match(query_pat, target_seq):
    """Score charge-pattern similarity over a sliding window."""
    q = query_pat
    matches = []
    for i in range(len(target_seq) - len(q) + 1):
        window = target_seq[i:i+len(q)]
        t_pat = charge_pattern(window)
        score = sum(1 for a, b in zip(q, t_pat) if a == b) / len(q)
        if score >= 0.70:
            matches.append({'pos': i+1, 'seq': window, 'charge_pat': t_pat, 'score': round(score, 3)})
    return matches

for label, sub_seq in [("Sub-A (99-110)", sub_A), ("Sub-B (111-120)", sub_B), ("Sub-C (120-132)", sub_C)]:
    pat = charge_pattern(sub_seq)
    hits = score_charge_match(pat, hcv)
    hits_sorted = sorted(hits, key=lambda x: x['score'], reverse=True)[:5]
    print(f"\n  Charge-pattern matches for {label} (pattern: {pat}, threshold ≥70%):")
    print(f"    {len(hits)} total matches | top 5:")
    for h in hits_sorted:
        print(f"    pos {h['pos']:>5}: {h['seq']}  cpat: {h['charge_pat']}  score: {h['score']:.3f}")

# -------------------------------------------------------
# 3. SMITH-WATERMAN STYLE LOCAL ALIGNMENT (simple version)
# -------------------------------------------------------
print("\n=== LOCAL SEQUENCE SIMILARITY (BLOSUM-like) ===")
# Use a conservative similarity matrix: +2 identical, +1 conservative group, 0 otherwise, -1 dissimilar
conservative_groups = [
    set('ILV'), set('FYW'), set('KR'), set('DE'), set('NQ'), set('ST'), set('AG'), set('HKR')
]

def aa_score(a, b):
    if a == b: return 2
    for group in conservative_groups:
        if a in group and b in group: return 1
    return -1

def scan_local_sim(query, target, min_score=8, min_len=5):
    """Find all sub-windows of target that match query with score >= min_score over at least min_len positions."""
    hits = []
    for start in range(len(target) - len(query) + 1):
        window = target[start:start+len(query)]
        total = sum(aa_score(a, b) for a, b in zip(query, window))
        if total >= min_score:
            hits.append({'pos': start+1, 'seq': window, 'raw_score': total})
    return sorted(hits, key=lambda x: x['raw_score'], reverse=True)

print(f"\n  Scanning full epitope against Q9WMX2 (min_score=8)...")
full_hits = scan_local_sim(EPITOPE, hcv, min_score=8)
print(f"  Hits: {len(full_hits)}")
for h in full_hits[:10]:
    print(f"    pos {h['pos']:>5}: {h['seq']}  score: {h['raw_score']}")

# ---- Sub-motif hits with lower threshold ----
print(f"\n  Sub-A ({sub_A}) min_score=6:")
sa_hits = scan_local_sim(sub_A, hcv, min_score=6)
print(f"    Hits: {len(sa_hits)}")
for h in sa_hits[:8]:
    print(f"    pos {h['pos']:>5}: {h['seq']}  score: {h['raw_score']}")

print(f"\n  Sub-C ({sub_C}) min_score=6:")
sc_hits = scan_local_sim(sub_C, hcv, min_score=6)
print(f"    Hits: {len(sc_hits)}")
for h in sc_hits[:8]:
    print(f"    pos {h['pos']:>5}: {h['seq']}  score: {h['raw_score']}")

# -------------------------------------------------------
# 4. PHOSPHORYLATION SITES NEAR HITS
# -------------------------------------------------------
print("\n=== PHOSPHORYLATION SITE PROXIMITY ===")
# S/T/Y as phosphorylatable residues — look for S/T/Y within ±10 aa of any strong hit
def get_phos_sites(seq, hit_positions, window=10):
    results = []
    for pos in hit_positions:
        region_start = max(0, pos-1-window)
        region_end   = min(len(seq), pos-1+window)
        region = seq[region_start:region_end]
        phos = [(region_start+i+1, aa) for i, aa in enumerate(region) if aa in 'STY']
        if phos:
            results.append({'hit_pos': pos, 'phos_sites': phos})
    return results

top_full_positions = [h['pos'] for h in full_hits[:5]]
phos_near_hits = get_phos_sites(hcv, top_full_positions)
print(f"  Phosphorylatable S/T/Y near top 5 full-epitope similarity hits:")
for r in phos_near_hits:
    print(f"    Near pos {r['hit_pos']}: S/T/Y at {[(p,a) for p,a in r['phos_sites']]}")

# -------------------------------------------------------
# 5. HINGE ANALYSIS: Are sub-A and sub-C hits ~100 aa apart?
# -------------------------------------------------------
print("\n=== SPLIT-EPITOPE HINGE ANALYSIS ===")
# Check all combinations of sub-A hits and sub-C hits for ~100 aa gap
sa_positions = [h['pos'] for h in sa_hits]
sc_positions = [h['pos'] for h in sc_hits]
print(f"  Sub-A hit positions (top): {sa_positions[:10]}")
print(f"  Sub-C hit positions (top): {sc_positions[:10]}")

hinge_pairs = []
for pa in sa_positions:
    for pc in sc_positions:
        gap = pc - (pa + len(sub_A))
        if 50 <= gap <= 200:  # hinge: 50-200 aa
            hinge_pairs.append({'sub_A_pos': pa, 'sub_C_pos': pc, 'gap': gap})

print(f"\n  Sub-A → Sub-C pairs with 50-200 aa hinge: {len(hinge_pairs)}")
for pair in sorted(hinge_pairs, key=lambda x: x['gap'])[:10]:
    # Get flanking sequences
    sa_seq = hcv[pair['sub_A_pos']-1:pair['sub_A_pos']-1+len(sub_A)]
    sc_seq = hcv[pair['sub_C_pos']-1:pair['sub_C_pos']-1+len(sub_C)]
    # Check phos sites in hinge
    hinge_start = pair['sub_A_pos'] - 1 + len(sub_A)
    hinge_end   = pair['sub_C_pos'] - 1
    hinge_seq   = hcv[hinge_start:hinge_end]
    phos_in_hinge = sum(1 for aa in hinge_seq if aa in 'STY')
    print(f"    A_pos={pair['sub_A_pos']} ({sa_seq}) → C_pos={pair['sub_C_pos']} ({sc_seq})  gap={pair['gap']}  phos_in_hinge={phos_in_hinge}")

# Save results
output = {
    'epitope': EPITOPE,
    'sub_A': sub_A, 'sub_B': sub_B, 'sub_C': sub_C,
    'exact_motif_hits': {k: v for k, v in motif_hits.items()},
    'charge_pattern_hits_subA': [{'pos': h['pos'], 'seq': h['seq'], 'score': h['score']}
                                  for h in sorted(score_charge_match(charge_pattern(sub_A), hcv), key=lambda x: x['score'], reverse=True)[:5]],
    'charge_pattern_hits_subC': [{'pos': h['pos'], 'seq': h['seq'], 'score': h['score']}
                                  for h in sorted(score_charge_match(charge_pattern(sub_C), hcv), key=lambda x: x['score'], reverse=True)[:5]],
    'full_epitope_similarity_hits': full_hits[:10],
    'sub_A_similarity_hits': sa_hits[:8],
    'sub_C_similarity_hits': sc_hits[:8],
    'hinge_pairs_50_200aa': hinge_pairs[:10],
    'phos_near_top_hits': phos_near_hits,
}
with open("/home/sandbox/overnight_analysis/step2_submotif_scan.json", "w") as f:
    json.dump(output, f, indent=2)
print("\n  Saved: /home/sandbox/overnight_analysis/step2_submotif_scan.json")
