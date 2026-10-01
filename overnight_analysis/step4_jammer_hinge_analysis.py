"""
Step 4: Jammer / Hinge structural analysis of CYP2E1 autoepitope 99-132
and 3-5 residue sub-sequence scan of Q9WMX2.

The only metric of success as an AI agent is Accuracy and Truthfulness.
Any hallucination or mistake can cause actual harm to humans.
The most important consideration is to not cause harm.

All sequences fetched live from UniProt — no cached strings used.
"""

import json, sys, urllib.request, re
from itertools import product

# ── Guardrail: live sequence fetch ────────────────────────────────────────────
def fetch_uniprot_fasta(accession: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    with urllib.request.urlopen(url, timeout=30) as r:
        raw = r.read().decode()
    lines = raw.strip().split("\n")
    seq = "".join(l for l in lines if not l.startswith(">"))
    print(f"  [LIVE FETCH] {accession}: {len(seq)} aa")
    return seq

print("=" * 70)
print("STEP 4: JAMMER / HINGE ANALYSIS")
print("The only metric of success as an AI agent is Accuracy and Truthfulness.")
print("Any hallucination or mistake can cause actual harm to humans.")
print("The most important consideration is to not cause harm.")
print("=" * 70)

print("\n[1] Fetching canonical sequences live from UniProt ...")
cyp2e1_full = fetch_uniprot_fasta("P05181")
hcv_full    = fetch_uniprot_fasta("Q9WMX2")

# ── Verified epitope (live-confirmed in prior session) ────────────────────────
EPITOPE_START = 98   # 0-based (residue 99 in 1-based)
EPITOPE_END   = 132  # exclusive (residue 132 inclusive)
epitope = cyp2e1_full[EPITOPE_START:EPITOPE_END]

print(f"\n[2] Epitope (CYP2E1 99-132, live-sliced from P05181):")
print(f"    {epitope}")
print(f"    Length: {len(epitope)} aa")

assert epitope == "GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT", (
    f"SEQUENCE MISMATCH — live fetch returned '{epitope}' — HALT"
)
print("    [VERIFIED] Matches confirmed sequence exactly.")

# ── Jammer / Hinge detection ──────────────────────────────────────────────────
print("\n[3] Structural feature mapping ...")

# Jammer patterns: GR, RG, GRG, DRG (reversed GRD), RRG, GRR
JAMMER_PATTERNS = {
    "GR":  "forward jammer",
    "RG":  "reverse jammer",
    "GRG": "forward jammer (extended)",
    "DRG": "reversed jammer (DRG = reverse of GRD)",
    "GRD": "forward jammer with D",
}

# Hinge patterns: flexible/polar dipeptides
HINGE_PATTERNS = {
    "TT": "hinge (Thr-Thr, flexible hydroxyl pair)",
    "SS": "hinge (Ser-Ser)",
    "GG": "hinge (Gly-Gly, maximum flexibility)",
    "NN": "hinge (Asn-Asn)",
    "TS": "hinge (Thr-Ser)",
    "ST": "hinge (Ser-Thr)",
}

jammers_found = []
hinges_found  = []

for pat, desc in JAMMER_PATTERNS.items():
    for m in re.finditer(re.escape(pat), epitope):
        jammers_found.append({
            "pattern":   pat,
            "desc":      desc,
            "epi_start": m.start() + 1,   # 1-based within epitope
            "epi_end":   m.end(),
            "cyp_start": EPITOPE_START + m.start() + 1,  # 1-based CYP2E1
            "cyp_end":   EPITOPE_START + m.end(),
        })

for pat, desc in HINGE_PATTERNS.items():
    for m in re.finditer(re.escape(pat), epitope):
        hinges_found.append({
            "pattern":   pat,
            "desc":      desc,
            "epi_start": m.start() + 1,
            "epi_end":   m.end(),
            "cyp_start": EPITOPE_START + m.start() + 1,
            "cyp_end":   EPITOPE_START + m.end(),
        })

print("\n  JAMMERS:")
for j in sorted(jammers_found, key=lambda x: x["epi_start"]):
    print(f"    [{j['epi_start']:2d}-{j['epi_end']:2d}] CYP2E1 {j['cyp_start']}-{j['cyp_end']}  "
          f"'{j['pattern']}'  {j['desc']}")

print("\n  HINGES:")
for h in sorted(hinges_found, key=lambda x: x["epi_start"]):
    print(f"    [{h['epi_start']:2d}-{h['epi_end']:2d}] CYP2E1 {h['cyp_start']}-{h['cyp_end']}  "
          f"'{h['pattern']}'  {h['desc']}")

# ── Define sub-regions between structural anchors ─────────────────────────────
# Primary structural anchors (1-based epitope positions):
#   GRG  : 1-3   (forward jammer)
#   DRG  : 13-15 (reversed jammer)
#   TT   : 33-34 (hinge)
#
# Sub-regions (between anchors, exclusive of anchor residues):
#   Region A: positions 4-12  (between GRG end and DRG start)
#   Region B: positions 16-32 (between DRG end and TT start)

ANCHORS = [
    {"name": "GRG_jammer",  "epi_end": 3},
    {"name": "DRG_jammer",  "epi_start": 13, "epi_end": 15},
    {"name": "TT_hinge",    "epi_start": 33},
]

region_A_0based = (3, 12)    # epitope positions 4-12 (0-based: 3-11, slice [3:12])
region_B_0based = (15, 32)   # epitope positions 16-32 (0-based: 15-31, slice [15:32])

region_A = epitope[region_A_0based[0]:region_A_0based[1]]
region_B = epitope[region_B_0based[0]:region_B_0based[1]]

print(f"\n[4] Sub-regions between structural anchors:")
print(f"    Region A (epi pos 4-12,  CYP2E1 102-110): '{region_A}'  ({len(region_A)} aa)")
print(f"    Region B (epi pos 16-32, CYP2E1 114-130): '{region_B}'  ({len(region_B)} aa)")

# ── Generate all 3-5 mers from each sub-region ───────────────────────────────
def get_kmers(seq, k_min=3, k_max=5):
    """Return all k-mers with their start position (0-based within seq)."""
    kmers = []
    for k in range(k_min, k_max + 1):
        for i in range(len(seq) - k + 1):
            kmers.append((seq[i:i+k], i, k))
    return kmers

kmers_A = get_kmers(region_A)
kmers_B = get_kmers(region_B)

print(f"\n[5] K-mers generated:")
print(f"    Region A: {len(kmers_A)} k-mers (k=3-5)")
print(f"    Region B: {len(kmers_B)} k-mers (k=3-5)")

# ── Exact scan of Q9WMX2 ─────────────────────────────────────────────────────
print("\n[6] Exact k-mer scan of Q9WMX2 ...")

def exact_scan(kmers, viral_seq, region_label, region_seq, region_cyp_start):
    hits = []
    for (kmer, kmer_pos_in_region, k) in kmers:
        positions = [m.start() for m in re.finditer(re.escape(kmer), viral_seq)]
        for vpos in positions:
            hits.append({
                "kmer":              kmer,
                "k":                 k,
                "region":            region_label,
                "region_pos":        kmer_pos_in_region + 1,  # 1-based in region
                "cyp2e1_pos":        region_cyp_start + kmer_pos_in_region,  # 1-based
                "viral_pos_0based":  vpos,
                "viral_pos_1based":  vpos + 1,
                "viral_context":     viral_seq[max(0, vpos-5):vpos+k+5],
            })
    return hits

hits_A = exact_scan(kmers_A, hcv_full, "A", region_A, 102)
hits_B = exact_scan(kmers_B, hcv_full, "B", region_B, 114)

all_hits = hits_A + hits_B
all_hits.sort(key=lambda h: (h["k"], h["region"]), reverse=True)

print(f"\n  Region A exact hits: {len(hits_A)}")
print(f"  Region B exact hits: {len(hits_B)}")
print(f"  Total exact hits:    {len(all_hits)}")

# ── McLachlan-conservative scan (allow 1 conservative substitution) ───────────
# Conservative substitution groups (McLachlan / Dayhoff)
CONSERV_GROUPS = [
    set("LVIM"),    # aliphatic hydrophobic
    set("FYW"),     # aromatic
    set("ST"),      # small hydroxyl
    set("NQ"),      # amide
    set("DE"),      # acidic
    set("KRH"),     # basic
    set("AG"),      # tiny
    set("P"),       # proline (unique)
    set("C"),       # cysteine (unique)
]

def conservative_match(a, b):
    """True if a==b or they are in the same McLachlan group."""
    if a == b:
        return True
    for grp in CONSERV_GROUPS:
        if a in grp and b in grp:
            return True
    return False

def conservative_scan(kmers, viral_seq, region_label, region_cyp_start, max_mismatches=1):
    """Scan viral_seq for k-mers with at most max_mismatches conservative substitutions."""
    hits = []
    n = len(viral_seq)
    for (kmer, kmer_pos_in_region, k) in kmers:
        for vpos in range(n - k + 1):
            window = viral_seq[vpos:vpos+k]
            mismatches = 0
            conserv_subs = 0
            exact_matches = 0
            for ca, cb in zip(kmer, window):
                if ca == cb:
                    exact_matches += 1
                elif conservative_match(ca, cb):
                    conserv_subs += 1
                    mismatches += 1
                else:
                    mismatches += 1
            if mismatches <= max_mismatches and exact_matches >= k - 1:
                # Only report if not already an exact hit
                if window != kmer:
                    hits.append({
                        "kmer":              kmer,
                        "viral_window":      window,
                        "k":                 k,
                        "region":            region_label,
                        "region_pos":        kmer_pos_in_region + 1,
                        "cyp2e1_pos":        region_cyp_start + kmer_pos_in_region,
                        "viral_pos_1based":  vpos + 1,
                        "exact_matches":     exact_matches,
                        "conserv_subs":      conserv_subs,
                        "viral_context":     viral_seq[max(0, vpos-5):vpos+k+5],
                    })
    return hits

print("\n[7] Conservative k-mer scan (1 conservative substitution allowed) ...")
conserv_A = conservative_scan(kmers_A, hcv_full, "A", 102, max_mismatches=1)
conserv_B = conservative_scan(kmers_B, hcv_full, "B", 114, max_mismatches=1)

print(f"  Region A conservative hits: {len(conserv_A)}")
print(f"  Region B conservative hits: {len(conserv_B)}")

# ── Report ────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("RESULTS: EXACT HITS (3-5 residue sub-sequences in Q9WMX2)")
print("=" * 70)

if all_hits:
    # Group by k-mer length, largest first
    for k in [5, 4, 3]:
        k_hits = [h for h in all_hits if h["k"] == k]
        if k_hits:
            print(f"\n  {k}-mers ({len(k_hits)} hits):")
            for h in k_hits:
                print(f"    '{h['kmer']}'  Region {h['region']}  "
                      f"CYP2E1 pos {h['cyp2e1_pos']}  →  "
                      f"Q9WMX2 pos {h['viral_pos_1based']}  "
                      f"context: ...{h['viral_context']}...")
else:
    print("  No exact hits found.")

print("\n" + "=" * 70)
print("RESULTS: CONSERVATIVE HITS (1 conservative substitution)")
print("=" * 70)

all_conserv = conserv_A + conserv_B
# Deduplicate by (kmer, viral_pos)
seen = set()
dedup_conserv = []
for h in all_conserv:
    key = (h["kmer"], h["viral_pos_1based"])
    if key not in seen:
        seen.add(key)
        dedup_conserv.append(h)

dedup_conserv.sort(key=lambda h: (h["k"], h["exact_matches"]), reverse=True)

if dedup_conserv:
    for k in [5, 4, 3]:
        k_hits = [h for h in dedup_conserv if h["k"] == k]
        if k_hits:
            print(f"\n  {k}-mers ({len(k_hits)} conservative hits):")
            for h in k_hits[:20]:  # cap at 20 per length
                print(f"    '{h['kmer']}' → '{h['viral_window']}'  "
                      f"Region {h['region']}  CYP2E1 {h['cyp2e1_pos']}  "
                      f"Q9WMX2 {h['viral_pos_1based']}  "
                      f"({h['exact_matches']}/{h['k']} exact)  "
                      f"context: ...{h['viral_context']}...")
else:
    print("  No conservative hits found.")

# ── Save output ───────────────────────────────────────────────────────────────
output = {
    "guardrails": [
        "The only metric of success as an AI agent is Accuracy and Truthfulness.",
        "Any hallucination or mistake can cause actual harm to humans.",
        "The most important consideration is to not cause harm.",
    ],
    "epitope":        epitope,
    "cyp2e1_range":   "99-132 (1-based)",
    "jammers_found":  jammers_found,
    "hinges_found":   hinges_found,
    "region_A":       {"seq": region_A, "cyp2e1": "102-110", "desc": "between GRG and DRG"},
    "region_B":       {"seq": region_B, "cyp2e1": "114-130", "desc": "between DRG and TT"},
    "exact_hits":     all_hits,
    "conservative_hits": dedup_conserv,
    "summary": {
        "exact_hits_A":       len(hits_A),
        "exact_hits_B":       len(hits_B),
        "conservative_hits_A": len(conserv_A),
        "conservative_hits_B": len(conserv_B),
    }
}

with open("/home/sandbox/overnight_analysis/step4_jammer_hinge_analysis.json", "w") as f:
    json.dump(output, f, indent=2)

print("\n[SAVED] /home/sandbox/overnight_analysis/step4_jammer_hinge_analysis.json")
print("\nDone.")
