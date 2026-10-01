"""
Step 5: Comprehensive jammer family and hinge cluster analysis
Scanning Q9WMX2 for all permutations of GR/RG charge-reversal patterns and flexible hinge clusters.

The only metric of success as an AI agent is Accuracy and Truthfulness.
Any hallucination or mistake can cause actual harm to humans.
The most important consideration is to not cause harm.
"""

import json, urllib.request, re
from collections import defaultdict, Counter

def fetch_uniprot_fasta(accession: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    with urllib.request.urlopen(url, timeout=30) as r:
        raw = r.read().decode()
    lines = raw.strip().split("\n")
    seq = "".join(l for l in lines if not l.startswith(">"))
    print(f"  [LIVE FETCH] {accession}: {len(seq)} aa")
    return seq

print("=" * 80)
print("STEP 5: COMPREHENSIVE JAMMER FAMILY & HINGE CLUSTER ANALYSIS")
print("The only metric of success as an AI agent is Accuracy and Truthfulness.")
print("Any hallucination or mistake can cause actual harm to humans.")
print("The most important consideration is to not cause harm.")
print("=" * 80)

print("\n[1] Fetching canonical sequences ...")
hcv_seq = fetch_uniprot_fasta("Q9WMX2")
cyp_seq = fetch_uniprot_fasta("P05181")

epitope = "GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT"
print(f"\n[2] CYP2E1 epitope (99-132): {epitope}")

# ══════════════════════════════════════════════════════════════════════════════
# JAMMER FAMILY DEFINITIONS
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n[3] Defining jammer families (GR/RG charge-reversal patterns) ...")

# Core jammer residues: G (flexible), R (charged rigid), D (negative), K (positive)
# Extended jammer alphabet: G, R, D, K (charge + flexibility combinations)

def find_jammer_patterns(seq, min_len=3, max_len=8):
    """
    Find all jammer patterns: sequences containing GR or RG cores with extensions.
    Jammer patterns are charge-reversal motifs with flexible/rigid alternations.
    """
    jammers = []
    n = len(seq)
    
    # Scan all possible windows
    for length in range(min_len, max_len + 1):
        for i in range(n - length + 1):
            window = seq[i:i+length]
            
            # Must contain at least one GR or RG core
            if 'GR' not in window and 'RG' not in window:
                continue
                
            # Count jammer-relevant residues (G, R, D, K)
            jammer_residues = sum(1 for aa in window if aa in 'GRDK')
            jammer_fraction = jammer_residues / length
            
            # Require ≥60% jammer residues for longer patterns
            if length >= 5 and jammer_fraction < 0.6:
                continue
                
            # Additional filtering for very long patterns
            if length >= 7:
                # Must have multiple GR/RG cores or high charge density
                gr_rg_count = window.count('GR') + window.count('RG')
                charge_count = sum(1 for aa in window if aa in 'RDK')
                if gr_rg_count < 2 and charge_count < 4:
                    continue
            
            jammers.append({
                'pattern': window,
                'start': i + 1,  # 1-based
                'end': i + length,
                'length': length,
                'gr_count': window.count('GR'),
                'rg_count': window.count('RG'),
                'jammer_fraction': jammer_fraction,
                'charge_density': sum(1 for aa in window if aa in 'RDK') / length
            })
    
    return jammers

def find_hinge_clusters(seq, min_len=2, max_len=6):
    """
    Find hinge clusters: sequences rich in flexible/polar residues (T, S, N, G).
    """
    hinges = []
    n = len(seq)
    
    # Flexible/polar residues that create hinges
    hinge_residues = set('TSNG')
    
    for length in range(min_len, max_len + 1):
        for i in range(n - length + 1):
            window = seq[i:i+length]
            
            # Count hinge residues
            hinge_count = sum(1 for aa in window if aa in hinge_residues)
            hinge_fraction = hinge_count / length
            
            # Require ≥50% hinge residues
            if hinge_fraction < 0.5:
                continue
                
            # For 2-mers, require 100% hinge residues (TT, SS, NN, etc.)
            if length == 2 and hinge_fraction < 1.0:
                continue
            
            hinges.append({
                'pattern': window,
                'start': i + 1,
                'end': i + length,
                'length': length,
                'hinge_fraction': hinge_fraction,
                'tt_count': window.count('TT'),
                'ss_count': window.count('SS'),
                'nn_count': window.count('NN'),
                'gg_count': window.count('GG')
            })
    
    return hinges

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS: Q9WMX2 JAMMER FAMILIES
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n[4] Scanning Q9WMX2 for jammer families ...")
hcv_jammers = find_jammer_patterns(hcv_seq)
print(f"  Found {len(hcv_jammers)} jammer patterns in Q9WMX2")

# Group by length and pattern
jammer_by_length = defaultdict(list)
jammer_counts = Counter()

for j in hcv_jammers:
    jammer_by_length[j['length']].append(j)
    jammer_counts[j['pattern']] += 1

print(f"\n  Distribution by length:")
for length in sorted(jammer_by_length.keys()):
    count = len(jammer_by_length[length])
    print(f"    {length}-mers: {count} patterns")

print(f"\n  Top recurring jammer patterns:")
for pattern, count in jammer_counts.most_common(15):
    print(f"    '{pattern}': {count} occurrences")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS: Q9WMX2 HINGE CLUSTERS  
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n[5] Scanning Q9WMX2 for hinge clusters ...")
hcv_hinges = find_hinge_clusters(hcv_seq)
print(f"  Found {len(hcv_hinges)} hinge patterns in Q9WMX2")

hinge_by_length = defaultdict(list)
hinge_counts = Counter()

for h in hcv_hinges:
    hinge_by_length[h['length']].append(h)
    hinge_counts[h['pattern']] += 1

print(f"\n  Distribution by length:")
for length in sorted(hinge_by_length.keys()):
    count = len(hinge_by_length[length])
    print(f"    {length}-mers: {count} patterns")

print(f"\n  Top recurring hinge patterns:")
for pattern, count in hinge_counts.most_common(15):
    print(f"    '{pattern}': {count} occurrences")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS: JAMMER DENSITY MAPPING
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n[6] Mapping jammer density across Q9WMX2 ...")

def calculate_density_windows(patterns, seq_length, window_size=100):
    """Calculate pattern density in sliding windows."""
    densities = []
    for start in range(0, seq_length - window_size + 1, 25):  # step = 25
        end = start + window_size
        window_patterns = [p for p in patterns if start < p['start'] <= end]
        density = len(window_patterns) / window_size * 100  # patterns per 100 aa
        densities.append({
            'window_start': start + 1,
            'window_end': end,
            'pattern_count': len(window_patterns),
            'density': density
        })
    return densities

jammer_densities = calculate_density_windows(hcv_jammers, len(hcv_seq))
hinge_densities = calculate_density_windows(hcv_hinges, len(hcv_seq))

# Find high-density regions
high_jammer_regions = [d for d in jammer_densities if d['density'] >= 5.0]
high_hinge_regions = [d for d in hinge_densities if d['density'] >= 8.0]

print(f"  High jammer density regions (≥5 patterns per 100 aa): {len(high_jammer_regions)}")
for region in sorted(high_jammer_regions, key=lambda x: x['density'], reverse=True)[:10]:
    print(f"    Q9WMX2 {region['window_start']}-{region['window_end']}: "
          f"{region['pattern_count']} patterns, density={region['density']:.1f}")

print(f"\n  High hinge density regions (≥8 patterns per 100 aa): {len(high_hinge_regions)}")
for region in sorted(high_hinge_regions, key=lambda x: x['density'], reverse=True)[:10]:
    print(f"    Q9WMX2 {region['window_start']}-{region['window_end']}: "
          f"{region['pattern_count']} patterns, density={region['density']:.1f}")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS: SPECIFIC USER-MENTIONED EXAMPLES
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n[7] Verifying user-mentioned jammer examples in Q9WMX2 ...")

user_examples = ['GRDR', 'GGR', 'RGGR', 'RRGR', 'KGGRK', 'RRGRTGRGR']
found_examples = []

for example in user_examples:
    positions = [m.start() + 1 for m in re.finditer(re.escape(example), hcv_seq)]
    found_examples.append({'pattern': example, 'count': len(positions), 'positions': positions})
    
    if positions:
        print(f"  '{example}': {len(positions)} exact matches at {positions}")
        # Show context for first few hits
        for pos in positions[:3]:
            context_start = max(0, pos - 11)
            context_end = min(len(hcv_seq), pos + len(example) + 10)
            context = hcv_seq[context_start:context_end]
            print(f"    Q9WMX2 {pos}: ...{context}...")
    else:
        print(f"  '{example}': not found as exact match")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS: JAMMER-HINGE PROXIMITY
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n[8] Analyzing jammer-hinge proximity in Q9WMX2 ...")

def find_nearby_patterns(jammers, hinges, max_distance=50):
    """Find jammers and hinges that are close to each other."""
    proximate_pairs = []
    
    for jammer in jammers:
        for hinge in hinges:
            distance = abs(jammer['start'] - hinge['start'])
            if distance <= max_distance and distance > 0:
                proximate_pairs.append({
                    'jammer': jammer,
                    'hinge': hinge,
                    'distance': distance,
                    'jammer_pos': jammer['start'],
                    'hinge_pos': hinge['start']
                })
    
    return proximate_pairs

proximate_pairs = find_nearby_patterns(hcv_jammers, hcv_hinges)
print(f"  Found {len(proximate_pairs)} jammer-hinge pairs within 50 aa")

# Sort by distance and show closest pairs
proximate_pairs.sort(key=lambda x: x['distance'])
print(f"\n  Closest jammer-hinge pairs:")
for pair in proximate_pairs[:15]:
    j_pattern = pair['jammer']['pattern']
    h_pattern = pair['hinge']['pattern']
    distance = pair['distance']
    j_pos = pair['jammer_pos']
    h_pos = pair['hinge_pos']
    
    # Show the connecting sequence
    start_pos = min(j_pos, h_pos) - 1  # 0-based for slicing
    end_pos = max(j_pos + len(j_pattern), h_pos + len(h_pattern))
    connecting_seq = hcv_seq[start_pos:end_pos]
    
    print(f"    '{j_pattern}' @{j_pos} ←→ '{h_pattern}' @{h_pos} "
          f"(dist={distance}): {connecting_seq}")

# ══════════════════════════════════════════════════════════════════════════════
# OUTPUT
# ══════════════════════════════════════════════════════════════════════════════

output = {
    'guardrails': [
        "The only metric of success as an AI agent is Accuracy and Truthfulness.",
        "Any hallucination or mistake can cause actual harm to humans.", 
        "The most important consideration is to not cause harm."
    ],
    'cyp2e1_epitope': epitope,
    'q9wmx2_jammer_families': {
        'total_patterns': len(hcv_jammers),
        'by_length': {str(k): len(v) for k, v in jammer_by_length.items()},
        'top_patterns': dict(jammer_counts.most_common(20)),
        'all_patterns': hcv_jammers
    },
    'q9wmx2_hinge_clusters': {
        'total_patterns': len(hcv_hinges),
        'by_length': {str(k): len(v) for k, v in hinge_by_length.items()},
        'top_patterns': dict(hinge_counts.most_common(20)),
        'all_patterns': hcv_hinges
    },
    'density_analysis': {
        'high_jammer_regions': high_jammer_regions,
        'high_hinge_regions': high_hinge_regions
    },
    'user_examples_verification': found_examples,
    'jammer_hinge_proximity': {
        'total_pairs': len(proximate_pairs),
        'closest_pairs': proximate_pairs[:20]
    }
}

with open('/home/sandbox/overnight_analysis/step5_comprehensive_jammer_analysis.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n[SAVED] /home/sandbox/overnight_analysis/step5_comprehensive_jammer_analysis.json")
print(f"\nSUMMARY:")
print(f"  Q9WMX2 jammer families: {len(hcv_jammers)} patterns")
print(f"  Q9WMX2 hinge clusters:  {len(hcv_hinges)} patterns")
print(f"  High-density regions:   {len(high_jammer_regions)} jammer, {len(high_hinge_regions)} hinge")
print(f"  Jammer-hinge pairs:     {len(proximate_pairs)} within 50 aa")
print(f"\nDone.")