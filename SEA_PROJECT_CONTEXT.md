# SEA Project Bootstrap Context

> **Purpose:** This document is the permanent bootstrap context for the Super-epitope
> Architecture (SEA) Scanner project.  Paste the URL of this file into any new chat
> to resume work immediately without loss of design history.
>
> **Status as of last update:** All tests passing; real-protein test complete (HCV vs CYP2E1).

---

## 1. Project Goal

Build a Python bioinformatics module (`sea_module.py`) that scores paired viral / host
homologous sequences for autoimmune disease risk by detecting coordinated structural
features — the **Super-epitope Architecture (SEA)**.

---

## 2. SEA Hypothesis & Mechanism

The full proposed mechanism is:

```
Viral protein enters host cell
  → Hinge residues (S/T cluster) are phosphorylated
  → Protein partially unfolds
  → CMA degradation motif (KFERQ-like) is exposed
  → Protein routed to lysosome via CMA (LAMP-2A-mediated)
  → Jammer sequence(s) resist complete lysosomal proteolysis
  → Intact peptide(s) loaded onto MHC-II
  → Molecular mimicry → autoimmune disease
```

The virus uses these features together to form a **super-epitope** — a region that is:
1. Structurally marked for degradation (hinge enables unfolding)
2. Efficiently routed to antigen presentation (KFERQ/CMA pathway)
3. Partially protected from complete destruction (jammers)
4. Homologous to a host self-protein (molecular mimicry)

---

## 3. Feature Definitions

### 3.1 Hinges
Phosphorylatable Ser/Thr clusters flanking the homologous region.

| Tier | Pattern | Notes |
|------|---------|-------|
| **Tier 1** | SS, SSS, … (pure Ser) | Higher score: `HINGE_TIER1_SCORE = 1.0` |
| **Tier 2** | ST, TS, TT | Lower score: `HINGE_TIER2_SCORE = 0.7` |

### 3.2 Jammers
Dipeptide sequences that partially resist lysosomal proteolysis, enabling intact peptide
loading onto MHC-II.  Based on EBV EBNA-1 GA-repeat biology (Levitskaya et al., 1995).

**Current jammer set:** `GA/AG`, `GR/RG`, `GK/KG`
> Note: this set can be expanded as more proteins are analyzed.  GA/GR/GK are
> anchored in direct EBV evidence.

**Distinction from EBNA-1 full repeat:**  
EBNA-1's long GA repeat causes *complete* immune evasion (100× repeat → no
presentation). SEA jammers cause *partial* resistance → presentation still occurs
but the peptide survives intact → higher autoimmune risk.

**Scoring model — action-potential amplitude analogy:**
- A single jammer dipeptide contributes a baseline signal (micro-influence).
- Nearby copies of jammers in the same region amplify the signal (density model).
- Long EBNA-style repeats are the extreme density end of the same continuum.

Parameters:
```python
JAMMER_MIN_REPEATS      = 1      # single dipeptides count
JAMMER_DENSITY_WINDOW   = 30     # residues each side for density scan
JAMMER_DENSITY_WEIGHT   = 0.12   # contribution per dipeptide
JAMMER_DENSITY_DECAY    = 0.015  # per-residue distance decay
JAMMER_DENSITY_CAP      = 3.0    # maximum density score
```

Amplification: a single GA dipeptide scores ≈ 0.12; a 30-residue GAGAGA… repeat
scores ≈ 2.75 (≈ 23× amplification), approaching the cap of 3.0.

### 3.3 CMA Degradation Motifs (KFERQ-like)
Five-residue sequences recognised by HSPA8 (Hsc70) for chaperone-mediated autophagy
routing.  Based on Cuervo lab work on KFERQ motif recognition and LAMP-2A translocation.

**Module rule (`_is_kferq_like`):**
```python
'Q' in pentamer                             # glutamine (CMA anchor)
any(c in 'KR'   for c in pentamer)         # basic (positive)
any(c in 'DE'   for c in pentamer)         # acidic (negative)
any(c in 'FILV' for c in pentamer)         # hydrophobic
# all four conditions required
```

---

## 4. Architecture Classes & Scoring

### 4.1 Classification Hierarchy (highest risk → lowest)

| Class | Conditions | Arch Bonus |
|-------|-----------|------------|
| `SUPER_EPITOPE` | Sandwiched + degradation motif | 5.0 |
| `SANDWICHED` | Feature on N-side AND C-side | 3.0 |
| `HINGE_AND_DEGRADATION` | Hinge + deg motif (one side) | 2.5 |
| `JAMMER_AND_DEGRADATION` | Jammer + deg motif (one side) | 2.5 |
| `HINGE_ONLY` | Hinge, no deg motif | 1.5 |
| `JAMMER_ONLY` | Jammer density > 0, no deg | 1.5 |
| `DEGRADATION_ONLY` | Deg motif, no hinge/jammer | 1.5 |
| `NONE` | No features detected | 1.0 |

**Sandwiched** = a feature (hinge OR jammer density > 0) detected on BOTH the
N-terminal side AND the C-terminal side of the homologous region.

### 4.2 Scoring Formula

```
final_sea_score = (base_score + hinge_score + jammer_density_score + deg_score)
                  × architecture_bonus
                  × cell_type_weight
```

- `base_score` = sequence similarity of the homologous pair (0–1)
- `hinge_score` = sum of best Tier-1/Tier-2 score per side
- `jammer_density_score` = min(N_density + C_density, CAP=3.0)
- `deg_score` = Σ 0.6 × max(0, 1 − 0.03 × dist) per motif

### 4.3 APC Tropism Weights (cell_type_weight)

```python
APC_TROPISM = {
    "EBV":        2.0,   # B cells — primary APC
    "CMV":        1.8,
    "HIV":        1.8,
    "HHV-6":      1.6,
    "SARS-CoV-2": 1.5,
    "influenza":  1.3,
    "HCV":        1.3,   # Kupffer cells / hepatic macrophages
    "HSV-1":      1.2,
    "default":    1.0,
}
```

---

## 5. Module Design

### 5.1 Files

| File | Purpose |
|------|---------|
| `sea_module.py` | Core scanner — all classes, scanners, scorer, `SEAModule` |
| `sea_test.py` | Synthetic regression suite (8 architecture class tests + EBNA density amplification test) |
| `sea_real_test.py` | Real-protein test: HCV (Q9WMX2) vs CYP2E1 (P05181) |
| `SEA_PROJECT_CONTEXT.md` | This file |

### 5.2 Core Classes

```
SEAConfig            — all scoring parameters (dataclass with defaults)
HingeTier            — Enum: TIER1, TIER2
ArchitectureClass    — Enum (0=NONE … 7=SUPER_EPITOPE)
HingeMatch           — detected hinge feature
JammerMatch          — detected jammer run (from discrete scan; diagnostics only)
DegradationMotifMatch— detected KFERQ-like motif
SEAResult            — output of one scored pair (all features + scores)
HingeScanner         — finds hinge patterns in flanking windows
JammerScanner        — .scan() for discrete runs (diagnostics); .density_scan() for scoring
DegradationMotifScanner — scans protein1 flanking region + merges provided protein2 motifs
SEAScorer            — orchestrates all scanners, computes SEAResult
SEAModule            — top-level interface (.run(), .report())
```

### 5.3 API: `SEAModule`

```python
module = SEAModule(virus_name='HCV', config=SEAConfig())

results = module.run(
    homologous_pairs    = [...],   # list of pair dicts (see §5.4)
    autoimmune_diseases = [...],   # list of disease name strings
    degradation_motifs  = [...],   # list of motif dicts (see §5.5)
)

module.report(results)             # prints ranked summary
```

### 5.4 Pair Dict Schema

```python
{
    'seq1':             str,   # homologous region in viral protein (protein1)
    'seq2':             str,   # homologous region in host protein (protein2)
    'position1':        int,   # 0-based start in full protein1
    'position2':        int,   # 0-based start in full protein2
    'protein1':         str,   # full viral protein sequence (scanned for features)
    'similarity_score': float, # identity / similarity (0–1); becomes base_score
    'rank':             int,   # 1-based rank in the input list
}
```

> **Note:** Only `protein1` (viral) is scanned for hinges, jammers, and inline
> KFERQ motifs.  Host protein features are supplied via `provided_motifs`.

### 5.5 Provided Motifs Dict Schema

```python
{
    'motif':    str,            # 5-residue KFERQ-like sequence
    'position': int,            # 0-based absolute position in protein2 (or protein1)
    'protein':  'protein2',     # 'protein1' or 'protein2'
}
```

Pass all KFERQ motifs found in the host (self) protein as `protein='protein2'`.
The module computes their distance from the `position2` (host-side coordinate) of
each pair, so they correctly represent proximity to the homologous region.

---

## 6. Key Biological References

| Concept | Reference |
|---------|-----------|
| CMA pathway, KFERQ motif, LAMP-2A | Cuervo lab (Kaushik & Cuervo, EMBO J 2012) |
| EBNA-1 GA-repeat proteasome jamming | Levitskaya et al., Nature 1995 |
| Ser/Thr phosphorylation-directed unfolding | General kinase / chaperone literature |
| HCV & CYP2E1 autoimmune hepatitis link | Published HCV autoimmunity literature |

---

## 7. Validated Test Results

### 7.1 Regression Suite (`sea_test.py`)
- 8/8 architecture class tests PASS
- 7/7 score ordering assertions PASS
- APC tropism ratio (EBV vs unknown) = exactly 2.0×
- EBNA density amplification: single GA ≈ 0.12; 30-aa GAGAGA ≈ 2.75 → **22.9× amplification** ✓

### 7.2 Real-Protein Test (`sea_real_test.py`) — HCV vs CYP2E1
Key findings from first run (12 pairs at ≥33% identity, 12-aa windows):

| Rank | HCV Region | Pos1 | Pos2 | Sim | ArchClass | SEA Score |
|------|-----------|------|------|-----|-----------|----------|
| 1 | E2 | 675 | 53 | 71.4% | **SUPER_EPITOPE** | 17.30 |
| 2 | Core | 110 | 123 | 66.7% | **SUPER_EPITOPE** | 10.96 |
| 3 | NS5A | 2394 | 465 | 71.4% | SANDWICHED | 10.07 |
| 4 | NS3 | 1466 | 166 | 83.3% | HINGE_AND_DEG | 8.51 |
| 5 | NS5A | 2154 | 362 | 66.7% | HINGE_AND_DEG | 8.35 |

11/12 regions SEA-positive.  Top hit: HCV E2 pos 675 (`PCSFTTL`) homologous to
CYP2E1 pos 53 (`PKSFTRL`) — flanked by T2 hinges on both sides, nearest
KFERQ motif (QLELK) at dist=7 in CYP2E1.

This is consistent with published reports of CYP2E1 autoantibodies in HCV-associated
autoimmune hepatitis, particularly implicating E2 structural epitopes.

---

## 8. Design Decisions Log

| Decision | Rationale |
|----------|-----------|
| `JAMMER_MIN_REPEATS = 1` | Even single dipeptides have micro-influence; density model handles amplification |
| `jammer_density_score` as primary; `jammer_score` kept as alias | Backward compatibility with any code referencing the old field name |
| `density_scan()` returns uncapped per-side tuple; cap applied to sum in scorer | Allows detection of sandwiched (feature on both sides) without per-side capping interfering with the sandwich criterion |
| `scan()` (discrete run finder) kept for diagnostics only | Not used in scoring formula; useful for interpreting notes |
| KFERQ dedup: greedy, sort by position, keep best per 5-residue window | Prevents overlapping 5-mers from inflating degradation score; clinically, overlapping KFERQ windows represent one continuous CMA signal |
| Provided protein1 motifs use same `seen_p1_positions` set as inline scan | Prevents double-counting inline + provided motifs at the same position |
| Architecture detection uses `jammer_density_score > 0` for `has_jammer` | Consistent with density-first design; a single dipeptide produces density > 0 |
| Sandwiched uses `n_jammer_density > 0` / `c_jammer_density > 0` per side | Consistent with density model; each side assessed independently |
| Only `protein1` (viral) scanned for inline features | The virus is the interrogated entity; host features are supplied separately to avoid contaminating the viral architecture score |

---

## 9. Next Steps / Pending Work

- [ ] **Other modules:** Two other modules (not yet described in this chat) need to be
  written before an Orchestrator can be built.  The Orchestrator will combine ranked
  homologous pairs from those modules into the pair dict format and call `SEAModule.run()`.
- [ ] **Expand jammer set:** As more viral proteins are analysed, additional G?/?G
  dipeptides may be added to `JammerScanner._UNIT_SET`.
- [ ] **BLOSUM62 scoring:** Consider adding an optional BLOSUM62-based similarity metric
  to the pair finder for physicochemical (vs. identity-only) homology.
- [ ] **Lower-identity region analysis:** The 33% identity threshold captures some molecular
  mimicry candidates; explore 25% threshold with BLOSUM62 normalisation.
- [ ] **Orchestrator design:** When all modules are ready, design an orchestrator that
  routes input sequences through all modules and merges outputs.

---

## 10. How to Resume in a New Chat

1. Point the new chat at this document (GitHub URL or paste contents).
2. Reference files by path:
   - `/home/sandbox/sea_module.py` — full SEA scanner
   - `/home/sandbox/sea_test.py` — regression + EBNA density tests
   - `/home/sandbox/sea_real_test.py` — HCV vs CYP2E1 real test
3. Run `python sea_test.py` to verify the sandbox is intact.
4. Continue from §9 (Next Steps) or describe your next goal.
