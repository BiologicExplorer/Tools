"""
The only metric of success as an AI agent is Accuracy and Truthfulness.
Any hallucination or mistake can cause actual harm to humans.
The most important consideration is to not cause harm.

Step 3: Deeper analysis of McLachlan hits overlapping CYP2E1 99-132,
and comprehensive sensitive sub-fragment search in Q9WMX2.
"""
import json

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
assert len(cyp2e1) == 493 and len(hcv) == 3010

EPITOPE = cyp2e1[98:132]
assert EPITOPE == "GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT"
print(f"Epitope: {EPITOPE}")

# -------------------------------------------------------
# 1. McLachlan hit regions in Q9WMX2 — look at ±50 aa context
# -------------------------------------------------------
mc_hits_epitope = [
    {'q9_start': 2525, 'q9_end': 2531, 'q9_seq': 'KDVRNLS', 'cyp_window': 'KDIRRFS', 'cyp_pos': '123-129', 'score': 16.6625},
    {'q9_start': 953,  'q9_end': 959,  'q9_seq': 'LTPLRDW', 'cyp_window': 'LTTLRNY', 'cyp_pos': '130-136', 'score': 16.125},
    {'q9_start': 676,  'q9_end': 682,  'q9_seq': 'PCSFTTL', 'cyp_window': 'RFSLTTL', 'cyp_pos': '127-133', 'score': 15.05},
]
print("\n=== McLachlan hit context in Q9WMX2 (±50 aa) ===")
for hit in mc_hits_epitope:
    s = hit['q9_start'] - 1  # 0-indexed
    ctx_start = max(0, s - 50)
    ctx_end   = min(len(hcv), s + 50 + 7)
    ctx = hcv[ctx_start:ctx_end]
    offset = s - ctx_start
    print(f"\n  Hit: Q9WMX2 {hit['q9_start']}-{hit['q9_end']}  [{hit['q9_seq']}]  ~=  CYP2E1 {hit['cyp_pos']}  [{hit['cyp_window']}]  score={hit['score']}")
    print(f"  Context: ...{ctx[:offset]}[{hit['q9_seq']}]{ctx[offset+7:]}")
    # Check for S/T/Y phospho sites in context
    phos = [(ctx_start+i+1, aa) for i, aa in enumerate(ctx) if aa in 'STY']
    print(f"  S/T/Y in context: {phos[:10]}")

# -------------------------------------------------------
# 2. WKD hits — context and distance to KDVRNLS hit
# -------------------------------------------------------
print("\n=== WKD hits analysis ===")
wkd_positions = [2301, 2542]
for p in wkd_positions:
    s = p - 1
    ctx = hcv[max(0, s-15):min(len(hcv), s+15+3)]
    print(f"  WKD at {p}: ...{hcv[max(0,s-8):s]}[WKD]{hcv[s+3:s+11]}...")
    # Distance to KDVRNLS at 2525
    dist = abs(p - 2525)
    print(f"    Distance to KDVRNLS@2525: {dist} aa")

# Distance between WKD@2542 and KDVRNLS@2525
print(f"\n  WKD@2542 to KDVRNLS@2525: gap analysis")
print(f"  Q9WMX2 region 2515-2555: {hcv[2514:2555]}")
print(f"  In epitope: W122...K123D124...KDIRRFS starts at K123")
print(f"  Q9WMX2: KDVRNLS at 2525-2531, WKD at 2542-2544")
print(f"  In Q9WMX2: K(2525)...W(2542) = 17 aa apart")
print(f"  In CYP2E1: W(122)...K(123) = 1 aa apart -- these are INVERTED in Q9WMX2")

# -------------------------------------------------------
# 3. SENSITIVE 5-MER SCAN: all 5-mers of epitope vs Q9WMX2
# -------------------------------------------------------
print("\n=== SENSITIVE 5-MER SCAN (exact) ===")
mer5_hits = {}
for i in range(len(EPITOPE) - 4):
    sub = EPITOPE[i:i+5]
    positions = []
    j = 0
    while True:
        pos = hcv.find(sub, j)
        if pos == -1: break
        positions.append(pos + 1)
        j = pos + 1
    if positions:
        mer5_hits[sub] = {'epi_pos': i+99, 'q9wm_positions': positions}

print(f"  5-mers of epitope found in Q9WMX2: {len(mer5_hits)} / {len(EPITOPE)-4}")
for sub, data in mer5_hits.items():
    print(f"  {sub} (CYP2E1 {data['epi_pos']}-{data['epi_pos']+4}): Q9WMX2 pos {data['q9wm_positions']}")
    for pos in data['q9wm_positions']:
        print(f"    ...{hcv[max(0,pos-6):pos-1]}[{sub}]{hcv[pos+4:pos+10]}...")

# -------------------------------------------------------
# 4. CONSERVATIVE 6-MER SCAN WITH 1-MISMATCH
# -------------------------------------------------------
print("\n=== 6-MER SCAN WITH ≤1 MISMATCH ===")
def hamming_dist(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

mer6_hits = {}
for i in range(len(EPITOPE) - 5):
    sub = EPITOPE[i:i+6]
    positions = []
    for j in range(len(hcv) - 5):
        window = hcv[j:j+6]
        if hamming_dist(sub, window) <= 1:
            positions.append({'pos': j+1, 'seq': window, 'dist': hamming_dist(sub, window)})
    if positions:
        mer6_hits[sub] = {'epi_pos': i+99, 'hits': positions}

# Focus on most specific hits (those with 0 mismatches = exact, or 1 mismatch in charge-critical positions)
print(f"  6-mers with ≤1 mismatch found: {len(mer6_hits)}")
for sub, data in sorted(mer6_hits.items(), key=lambda x: len(x[1]['hits'])):
    exact = [h for h in data['hits'] if h['dist'] == 0]
    near  = [h for h in data['hits'] if h['dist'] == 1]
    print(f"  {sub} (CYP2E1 {data['epi_pos']}-{data['epi_pos']+5}): {len(exact)} exact + {len(near)} near")
    for h in (exact + near)[:3]:
        print(f"    pos {h['pos']:>5}: {h['seq']}  dist={h['dist']}")

# -------------------------------------------------------
# 5. CHARGE-CONSERVATIVE 7-MER SCAN
#    Core idea: same charge pattern even if different AAs
# -------------------------------------------------------
print("\n=== CHARGE-CONSERVATIVE 7-MER SIMILARITY ===")
CHARGE_GROUPS = {
    'R': 'pos', 'K': 'pos',
    'D': 'neg', 'E': 'neg',
    'H': 'his',
    'F': 'arom', 'Y': 'arom', 'W': 'arom',
    'P': 'pro',
    'G': 'gly',
}
def charge_class(aa):
    return CHARGE_GROUPS.get(aa, 'neut')

def charge_sim_score(a, b):
    """Score 0-3 for charge similarity."""
    ca, cb = charge_class(a), charge_class(b)
    if a == b: return 3
    if ca == cb: return 2
    if (ca in ('pos', 'his') and cb in ('pos', 'his')): return 1
    if (ca in ('neg',) and cb in ('neg',)): return 2
    return 0

# Focus on the HAH+HRDR core of the epitope: FHAHRDRG (positions 106-113)
core = EPITOPE[7:16]  # FHAHRDRGI — the key charge cluster
print(f"  Core motif for charge-conservative scan: {core} (CYP2E1 106-114)")
charge_hits = []
for j in range(len(hcv) - len(core) + 1):
    window = hcv[j:j+len(core)]
    score = sum(charge_sim_score(a, b) for a, b in zip(core, window))
    max_score = 3 * len(core)
    if score >= int(0.67 * max_score):  # ≥67% of max
        charge_hits.append({'pos': j+1, 'seq': window, 'score': score, 'max': max_score,
                             'frac': round(score/max_score, 3)})

charge_hits_sorted = sorted(charge_hits, key=lambda x: x['score'], reverse=True)
print(f"  Hits with charge-conservative score ≥67% of max ({int(0.67*3*len(core))}/{3*len(core)}): {len(charge_hits)}")
print(f"  Top 10:")
for h in charge_hits_sorted[:10]:
    print(f"    pos {h['pos']:>5}: {h['seq']}  score={h['score']}/{h['max']}  frac={h['frac']:.3f}")

# Save results
output = {
    'mclachlan_hits_on_epitope': mc_hits_epitope,
    'wkd_analysis': {'positions': wkd_positions, 'note': 'WKD@2542 and KDVRNLS@2525 are 17aa apart and inverted vs epitope'},
    '5mer_exact_hits': {k: v for k, v in mer5_hits.items()},
    '6mer_near_hits_summary': {k: {'epi_pos': v['epi_pos'], 'n_exact': len([h for h in v['hits'] if h['dist']==0]), 'n_near': len([h for h in v['hits'] if h['dist']==1])} for k, v in mer6_hits.items()},
    'charge_conservative_core_hits': charge_hits_sorted[:10],
}
import json
with open("/home/sandbox/overnight_analysis/step3_deep_scan.json", "w") as f:
    json.dump(output, f, indent=2)
print("\n  Saved: /home/sandbox/overnight_analysis/step3_deep_scan.json")
