#!/usr/bin/env python3
"""
The only metric of success as an AI agent is Accuracy and Truthfulness.
Any hallucination or mistake can cause actual harm to humans.
The most important consideration is to not cause harm.

Phase 1 Two-Pass McLachlan Search
==================================
Pass 1 : 11-mer gapless sliding window  Q9WMX2 x P05181, raw sum >= 22.0
Pass 2 : best 5-mer sub-window per P1 hit, raw sum >= 17.5 (avg >= 3.5/pos)

Sequences loaded LIVE from canonical FASTA files (no cached strings).
"""

import sys, json, time
sys.path.insert(0, '/home/sandbox/orchestrator')
from mclachlan_aligner import _raw

def load_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith('>'):
                seq.append(line.strip())
    return ''.join(seq).upper()

viral_seq = load_fasta('/home/sandbox/verified_sequences/Q9WMX2_canonical.fasta')
host_seq  = load_fasta('/home/sandbox/verified_sequences/P05181_canonical.fasta')

print(f"Q9WMX2 loaded: {len(viral_seq)} aa")
print(f"P05181 loaded: {len(host_seq)} aa")
assert len(viral_seq) > 3000
assert len(host_seq)  > 490

print("Building score table ...")
t0 = time.time()
nv = len(viral_seq)
nh = len(host_seq)
score_table = [
    [_raw(viral_seq[i], host_seq[j]) for j in range(nh)]
    for i in range(nv)
]
print(f"  done in {time.time()-t0:.1f}s")

# ---- Pass 1 ----
WIN1    = 11
THRESH1 = 22.0
print(f"\nPass 1: WIN={WIN1}, threshold={THRESH1}")
t0 = time.time()

raw_hits = []
for vi in range(nv - WIN1 + 1):
    for hi in range(nh - WIN1 + 1):
        s = sum(score_table[vi + k][hi + k] for k in range(WIN1))
        if s >= THRESH1:
            raw_hits.append((s, vi, hi))

print(f"  Raw hits: {len(raw_hits)},  time: {time.time()-t0:.1f}s")

raw_hits.sort(key=lambda x: -x[0])
used_viral = set()
pass1_hits = []

for score, vi, hi in raw_hits:
    ve = vi + WIN1 - 1
    if any(p in used_viral for p in range(vi, ve + 1)):
        continue
    used_viral.update(range(vi, ve + 1))
    pass1_hits.append({
        'viral_seq':   viral_seq[vi : vi + WIN1],
        'host_seq':    host_seq[hi  : hi + WIN1],
        'viral_start': vi + 1,
        'viral_end':   vi + WIN1,
        'host_start':  hi + 1,
        'host_end':    hi + WIN1,
        'score_11mer': round(score, 4),
    })

print(f"  After dedup: {len(pass1_hits)} hits")

# ---- Pass 2 ----
WIN2    = 5
THRESH2 = 17.5
print(f"\nPass 2: WIN={WIN2}, threshold={THRESH2}")

pass2_hits = []
for h in pass1_hits:
    vs0 = h['viral_start'] - 1
    hs0 = h['host_start']  - 1

    best_s5     = -1
    best_offset = -1

    for d in range(WIN1 - WIN2 + 1):
        s5 = sum(score_table[vs0 + d + k][hs0 + d + k] for k in range(WIN2))
        if s5 > best_s5:
            best_s5     = s5
            best_offset = d

    if best_s5 < THRESH2:
        continue

    d   = best_offset
    v5  = viral_seq[vs0 + d : vs0 + d + WIN2]
    h5  = host_seq [hs0 + d : hs0 + d + WIN2]
    vs5 = vs0 + d + 1
    ve5 = vs5 + WIN2 - 1
    hs5 = hs0 + d + 1
    he5 = hs5 + WIN2 - 1

    aa_scores_11 = [score_table[vs0 + k][hs0 + k] for k in range(WIN1)]
    aa_scores_5  = [score_table[vs0 + d + k][hs0 + d + k] for k in range(WIN2)]

    pass2_hits.append({
        'viral_seq':        h['viral_seq'],
        'host_seq':         h['host_seq'],
        'viral_start':      h['viral_start'],
        'viral_end':        h['viral_end'],
        'host_start':       h['host_start'],
        'host_end':         h['host_end'],
        'score_11mer':      h['score_11mer'],
        'aa_scores_11mer':  aa_scores_11,
        'core_viral_seq':   v5,
        'core_host_seq':    h5,
        'core_viral_start': vs5,
        'core_viral_end':   ve5,
        'core_host_start':  hs5,
        'core_host_end':    he5,
        'core_score_5mer':  round(best_s5, 4),
        'core_aa_scores':   aa_scores_5,
        'core_offset':      d,
    })

print(f"  Pass 2 hits: {len(pass2_hits)}")

# ---- Annotate ----
NS3_RANGE  = (1215, 1650)
NS5A_RANGE = (2008, 2170)
NS5B_RANGE = (2421, 2989)

def assign_protein(vs, ve):
    mid = (vs + ve) / 2
    if NS3_RANGE[0]  <= mid <= NS3_RANGE[1]:  return 'NS3'
    if NS5A_RANGE[0] <= mid <= NS5A_RANGE[1]: return 'NS5A'
    if NS5B_RANGE[0] <= mid <= NS5B_RANGE[1]: return 'NS5B'
    return 'other'

PDB_VIRAL_RANGES = {
    '3KQL': (1027, 1660),
    '1ZH1': (1973, 2180),
    '3FQL': (2420, 2995),
}

def viral_pdbs(vs, ve):
    return [pdb for pdb, (ps, pe) in PDB_VIRAL_RANGES.items() if vs <= pe and ve >= ps]

EPITOPE_START = 99
EPITOPE_END   = 132

def aa_class(score):
    if score >= 4: return 'S'
    if score >= 2: return 'P'
    return 'N'

for h in pass2_hits:
    h['protein']         = assign_protein(h['viral_start'], h['viral_end'])
    h['viral_pdbs']      = viral_pdbs(h['viral_start'], h['viral_end'])
    h['host_pdb']        = '3KOH'
    hs = h['core_host_start']
    he = h['core_host_end']
    h['in_epitope']      = (hs <= EPITOPE_END and he >= EPITOPE_START)
    h['aa_classes_5mer'] = [aa_class(s) for s in h['core_aa_scores']]
    h['n_strong']        = sum(1 for c in h['aa_classes_5mer'] if c == 'S')
    h['n_partial']       = sum(1 for c in h['aa_classes_5mer'] if c == 'P')

pass2_hits.sort(key=lambda h: (-h['core_score_5mer'], -h['score_11mer']))

output = {
    'method':            'two_pass_mclachlan_v1',
    'pass1_window':      WIN1,
    'pass1_threshold':   THRESH1,
    'pass2_window':      WIN2,
    'pass2_threshold':   THRESH2,
    'pass1_raw_count':   len(raw_hits),
    'pass1_dedup_count': len(pass1_hits),
    'pass2_count':       len(pass2_hits),
    'viral_seq_len':     len(viral_seq),
    'host_seq_len':      len(host_seq),
    'hits':              pass2_hits,
}

outfile = '/home/sandbox/phase4_viz/phase1_twopass_results.json'
with open(outfile, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved {len(pass2_hits)} hits -> {outfile}")

print(f"\n{'#':>3}  {'Prot':>5}  {'V-11mer':>13}  {'V-range':>11}  "
      f"{'H-11mer':>13}  {'H-range':>9}  {'Cv':>6}  {'Ch':>6}  "
      f"{'5mer':>5}  {'11mer':>6}  {'Cls':>5}  Ep  PDBs")
for i, h in enumerate(pass2_hits[:25]):
    vr   = f"{h['viral_start']}-{h['viral_end']}"
    hr   = f"{h['host_start']}-{h['host_end']}"
    cls  = ''.join(h['aa_classes_5mer'])
    pdbs = ','.join(h['viral_pdbs']) or '-'
    print(f"{i+1:>3}  {h['protein']:>5}  {h['viral_seq']:>13}  {vr:>11}  "
          f"{h['host_seq']:>13}  {hr:>9}  {h['core_viral_seq']:>6}  "
          f"{h['core_host_seq']:>6}  {h['core_score_5mer']:>5.1f}  "
          f"{h['score_11mer']:>6.1f}  {cls:>5}  "
          f"{'Y' if h['in_epitope'] else 'N'}  {pdbs}")

from collections import Counter
prot_counts = Counter(h['protein'] for h in pass2_hits)
ep_count    = sum(1 for h in pass2_hits if h['in_epitope'])
pdb_covered = sum(1 for h in pass2_hits if h['viral_pdbs'])
print(f"\nProtein breakdown : {dict(prot_counts)}")
print(f"Epitope overlaps  : {ep_count}")
print(f"PDB-covered hits  : {pdb_covered}")
print(f"Total pass2 hits  : {len(pass2_hits)}")
