# Protein Degradation Motif Finder — Agent Bootstrap Guide v2.0

## Purpose
This document gives a new agent everything it needs to find seven classes of
protein degradation motifs **without reading any PDFs or making external calls**.
Read this file, run the self-test, and you are ready.

---

## Script Location (GitHub)

| File | Raw URL |
|---|---|
| `motif_finder.py` | `https://raw.githubusercontent.com/BiologicExplorer/Tools/main/motif_finder.py` |
| `KFERQ_LIR_KNOWLEDGE.md` | `https://raw.githubusercontent.com/BiologicExplorer/Tools/main/KFERQ_LIR_KNOWLEDGE.md` |

**Bootstrap instruction for a fresh agent:**
> "Using ONLY the GitHub integration: get raw file content of `motif_finder.py`
> from repository `BiologicExplorer/Tools`, save to `/home/sandbox/motif_finder.py`,
> then run `python3 /home/sandbox/motif_finder.py --self-test`"

---

## Related Tool — CMA Network Lookup & Scoring

For pathway-level analysis (gene membership in the CMA regulatory network and
CMA activation scoring from expression data), see the **separate** companion tool:

| File | Raw URL |
|---|---|
| `cma_network.py` | `https://raw.githubusercontent.com/BiologicExplorer/Tools/main/cma_network.py` |
| `CMA_NETWORK_KNOWLEDGE.md` | `https://raw.githubusercontent.com/BiologicExplorer/Tools/main/CMA_NETWORK_KNOWLEDGE.md` |

`cma_network.py` and `motif_finder.py` are intentionally separate: `motif_finder.py`
operates on protein sequences (residue level); `cma_network.py` operates on gene
names and expression datasets (pathway level).

---

## The Script

**Location:** `BiologicExplorer/Tools` → `motif_finder.py`  
**Version:** 2.0  
**Quick self-test:**

```bash
python3 /home/sandbox/motif_finder.py --self-test
# Expected output: "ALL 31 TESTS PASSED — motif_finder.py v2.0 READY"
```

---

## Quick-Start for a Fresh Agent

```bash
# 1. Verify the script and knowledge file are present
ls -lh /home/sandbox/motif_finder.py /home/sandbox/KFERQ_LIR_KNOWLEDGE.md

# 2. Run the full self-test — all 31 lines must read PASS
python3 /home/sandbox/motif_finder.py --self-test

# 3. Read this file (you are doing that now)

# 4. Scan a sequence for all degradation motifs
python3 /home/sandbox/motif_finder.py --sequence "MSQKFERQSTPESTDDD"
```

---

## Unified Python API

```python
import sys; sys.path.insert(0, '/home/sandbox')
from motif_finder import find_all_degradation_motifs

result = find_all_degradation_motifs("MSQKFERQSTPESTDDD")
# Returns a dict with 8 keys:
# {
#   'kferq':    [ {start, end, motif, type, notes}, ... ],
#   'lir':      [ {start, end, motif, type, notes}, ... ],
#   'dbox':     [ {start, end, motif, type, notes}, ... ],
#   'ken':      [ {start, end, motif, type, notes}, ... ],
#   'n_degron': {start, end, motif, type, notes}   # single dict or None
#   'c_degron': [ {start, end, motif, type, notes}, ... ],
#   'pest':     [ {start, end, motif, type, notes}, ... ],
#   'summary':  { 'kferq': int, 'lir': int, ... }  # counts per class
# }
```

**Coordinate convention:** `start` and `end` are **1-based, inclusive**.

Individual engine functions are also importable:

```python
from motif_finder import (
    find_kferq_motifs,       # → list of dicts
    find_lir_motifs,         # → list of dicts
    find_dbox_motifs,        # → list of dicts
    find_ken_motifs,         # → list of dicts
    classify_n_degron,       # → single dict or None
    find_c_degron_motifs,    # → list of dicts
    find_pest_sequences,     # → list of dicts
    classify_protein_cma,    # → str (category label)
)
```

---

## CLI Usage

```bash
# All motifs, single sequence
python3 /home/sandbox/motif_finder.py --sequence MSQIRLKEFERQSLP...

# FASTA file (one or many proteins)
python3 /home/sandbox/motif_finder.py proteins.fasta

# Include N-bearing KFERQ motifs; suppress LIR output
python3 /home/sandbox/motif_finder.py proteins.fasta --n-bearing --no-lir

# Canonical KFERQ only (no phospho, no acetyl)
python3 /home/sandbox/motif_finder.py proteins.fasta --no-phospho --no-acetyl
```

---

## Motif Engine 1 — KFERQ-Like Motifs (CMA)

**Source:** Kirchner P et al. (2019). *Proteome-wide analysis of chaperone-mediated
autophagy targeting motifs.* PLoS Biology 17(6): e3000301.
DOI: 10.1371/journal.pbio.3000301

### Biological context
KFERQ-like motifs are pentapeptides recognized by the cytosolic chaperone HSC70.
Binding leads to protein translocation into the lysosome via LAMP2A —
chaperone-mediated autophagy (CMA). The motif is **necessary and sufficient** for
HSC70 recognition (Q-bearing) or **necessary but not sufficient alone** (N-bearing).

---

### Canonical motif

A sliding 5-residue window is **canonical** when:

| Position | Rule |
|---|---|
| **pos 1 OR pos 5** (flanking) | Must be **Q** (glutamine) |
| **Remaining 4 positions** (core) | Satisfy ALL three constraints below |

**Core constraints** (each residue must fall into exactly one class):

| Class | Residues | Count |
|---|---|---|
| Positively charged (basic) | **K, R** | 1 or 2 |
| Hydrophobic | **I, L, V, F** | 1 or 2 |
| Negatively charged (acidic) | **D, E** | exactly 1 |

Valid compositions: 1B+2H+1A or 2B+1H+1A (all four core positions must be
assigned; any unclassified residue invalidates the motif).

**Validated examples:**

| Motif | Q pos | Core | Protein |
|---|---|---|---|
| KFERQ | 5 | K(B) F(H) E(A) R(B) | Ribonuclease A |
| ILKEQ | 5 | I(H) L(H) K(B) E(A) | DJ-1 / PARK7 |
| VKKDQ | 5 | V(H) K(B) K(B) D(A) | α-Synuclein |
| LDRLQ | 5 | L(H) D(A) R(B) L(H) | PLIN3 |

---

### Phosphorylation-activated motif

Same as canonical but **S, T, or Y** occupies the acidic slot (negative charge
acquired upon phosphorylation). Non-functional until pS/pT/pY.
Enable with `allow_phospho=True` (default).

### Acetylation-activated motif

Same as canonical but the flanking position (1 or 5) is **K** rather than Q.
Acetylation neutralises K's charge, mimicking Q.
Enable with `allow_acetyl=True` (default).

### N-bearing motif (advanced)

Flanking Q replaced by **N** (asparagine). Validated in GAPDH and HIF1α but
**necessary not sufficient** — requires additional determinants.
Enable with `allow_n_bearing=True` (off by default in systematic scans).

### Protein-level classification

`classify_protein_cma(motifs)` returns one of:
`canonical` > `phospho_activated` > `acetyl_activated` > `n_bearing` > `no_motif`

---

## Motif Engine 2 — LIR / AIM Motifs (Selective Autophagy)

**Core rule sources:**
- Cheng et al. (2023) preprint; DOI: 10.1101/2022.09.25.509395
- Karin et al. (2014); DOI: 10.6084/m9.figshare.1601912.v1
- McLaughlin et al. (2020) *Biochem J* 477(16):3115–3130; DOI: 10.1042/bcj20200714
- Rozenknop A et al. (2011) *J Mol Biol* 410:622–631; DOI: 10.1016/j.jmb.2011.05.003

**xLIR / phospho-flanking sources:**
- Klionsky et al. (2014, iLIR); DOI: 10.4161/AUTO.28260
- Kohler et al. (2022, TEX264/CK2); DOI: 10.1101/2022.02.11.480038

**Non-canonical source:**
- AF2 SLiM screen (2025); DOI: 10.6084/m9.figshare.28929437.v1

### Biological context
LIR (LC3-Interacting Region) / AIM (Atg8-Interacting Motif) motifs dock into
the LIR-docking site (LDS) of ATG8-family proteins (LC3A/B/C, GABARAP,
GABARAPL1/L2). They mediate selective autophagy by linking cargo receptors,
autophagy adaptors, and ATG proteins to the autophagosome membrane.

### Core grammar

```
[W / F / Y] - x₁ - x₂ - [I / L / V]
```

Notation: Θ₀ (aromatic) — x₁ — x₂ — Γ₃ (aliphatic).
Affinity order of Θ₀: **W > F > Y** (Rozenknop 2011).
GABARAPL1 shows strong preference for Trp at Θ₀.

### xLIR / flanking classification

The three upstream positions relative to Θ₀ are labelled −3, −2, −1.

| Type tag | Upstream −3..−1 | Meaning |
|---|---|---|
| `lir_xlir_acidic` | ≥1 D or E | acidic context |
| `lir_xlir_phospho` | ≥1 S or T (no D/E) | phosphorylatable context; CK2 phosphorylates S/T at −1 or −2 (TEX264) |
| `lir_xlir_optimal` | ≥1 D/E **and** ≥1 S/T | optimal context |
| `lir_canonical` | none of the above | canonical, no xLIR context |

The upstream window is the up-to-3 residues immediately before Θ₀ (fewer
residues are used when Θ₀ is near the N-terminus).

### Non-canonical flag

If the match overlaps with the pattern `[DE]W[DE]` (from the AF2 SLiM screen),
the `notes` field gains `+DEWDE`.

### Output type strings

`lir_canonical` | `lir_xlir_acidic` | `lir_xlir_phospho` | `lir_xlir_optimal`
(plus optional `+DEWDE` suffix in `notes`)

---

## Motif Engine 3 — D-Box (Destruction Box)

**Sources:**
- Glotzer M et al. (1991) *Nature* 349:132–138 (minimal RxxL definition)
- Andrews PD et al. (2024) preprint; DOI: 10.1101/2024.04.30.590460 (extended RxxLxxxxN)

### Biological context
D-boxes (Destruction boxes) are recognized by the Cdh1/Cdc20 co-activators of
the APC/C ubiquitin ligase complex during mitotic exit and the G1/S transition.
They drive cell-cycle-regulated proteasomal degradation of cyclins, securin,
and other substrates.

### Grammar

| Pattern | Regex | Type tag | Note |
|---|---|---|---|
| Extended | `R.{1}[LI].{4}N` (9 aa) | `dbox_extended` | Higher specificity |
| Minimal | `R.{1}[LI]` (4 aa) | `dbox_minimal` | High false-positive risk |

**Default behaviour:** extended hits are reported; overlapping minimal hits within
an extended match are suppressed.

**Caution:** Minimal `RxxL` alone is a high false-positive pattern. Prefer the
extended motif for confident predictions; use minimal only for exploratory scans.

---

## Motif Engine 4 — KEN Box

**Source:** Pfleger CM & Kirschner MW (2000) *Genes Dev* 14:655–665
(DOI: 10.1101/gad.14.6.655)

### Biological context
KEN boxes are recognized by the APC/C co-activator Cdh1 (not Cdc20) during
mitotic exit and G1. They mediate ubiquitination and proteasomal degradation of
mitotic regulators alongside or independently of D-boxes.

### Grammar

```
K - E - N    (exact tripeptide)
```

Type tag: `ken_box`.
Note: KEN box and D-box can cooperate on the same substrate; both engines should
be run together for APC/C substrates.

---

## Motif Engine 5 — N-Degron

**Sources:**
- Timms RT & Koren I (2020) *Biochem Soc Trans* 48(4):1355–1364;
  DOI: 10.1042/BST20191094
- Dissmeyer N et al. (2016) *Curr Protein Pept Sci* 17(1):4–34;
  DOI: 10.2174/0929866523666160108115809

### Biological context
N-degrons are degradation signals determined by the identity of the N-terminal
residue after co- or post-translational processing. They feed into N-degron
pathways (formerly N-end rule) via sequential enzymatic steps and are ultimately
recognised by N-recognin E3 ubiquitin ligases (e.g., UBR1–UBR5) for
proteasomal degradation.

### N-terminal residue classification

| Class | Residues | Fate |
|---|---|---|
| Type I primary destabilising | R, H, K | Direct N-recognin recognition |
| Type II primary destabilising | F, Y, W, L, I | Direct N-recognin recognition |
| Secondary destabilising | D, E | Arginylated by ATE1 → type I |
| Tertiary destabilising | N, Q | Deamidated (NTAN1/NTAQ1) → D/E → secondary |
| Stabilising | M, A, G, V, C, P, S, T | Long half-life |

### MAP cleavage prediction (Methionine Aminopeptidase)

When residue 1 is **Met (M)** and residue 2 is in {A, C, G, P, S, T, V}, the
initiator Met is predicted to be cotranslationally cleaved, exposing residue 2
as the effective N-terminus. The `type` field in the returned dict reflects the
**post-cleavage** classification.

The engine returns a **single dict** (not a list), describing residue 1 (or the
post-cleavage residue if MAP cleavage is predicted):

```python
{
  'start': 1, 'end': 1,
  'motif': 'R',
  'type': 'type_I_primary_destabilising',
  'notes': 'N-terminal residue; N-degron class'
}
```

Returns `None` for sequences of length 0.

---

## Motif Engine 6 — C-Degron

**Source:** Timms RT & Koren I (2020) *Biochem Soc Trans* 48(4):1355–1364;
DOI: 10.1042/BST20191094

### Biological context
C-degrons are degradation signals at or near the protein C-terminus, recognised
by dedicated CRL2 and CRL4 E3 ubiquitin ligase complexes. Two well-characterised
classes are implemented:

### Grammar

| Type tag | Pattern | Recognising complex |
|---|---|---|
| `cdegron_appbp2` | `R.{2,4}G.{0,3}$` | CRL2-APPBP2 (also known as KLHDC2) |
| `cdegron_dcaf12` | `EE$` | CRL4-DCAF12 |

Both patterns are anchored to the protein C-terminus (`$`).

---

## Motif Engine 7 — PEST Sequences

**Sources:**
- Rogers S et al. (1986) *Science* 234:364–368
- Rechsteiner M & Rogers SW (1996) *Trends Biochem Sci* 21:267–271

### Biological context
PEST sequences are regions enriched in Pro (P), Glu (E), Ser (S), and Thr (T),
flanked by positively charged residues (K, R, H). They act as cis-acting
instability elements that promote rapid proteasomal or calpain-mediated
degradation. Classic substrates include c-Fos, c-Myc, IκBα, and HIF-2α.

### Structural criteria (all must be satisfied)

1. Flanked by K, R, or H on both sides (these flanking residues are NOT part of
   the PEST segment itself)
2. Internal segment length ≥ 12 amino acids
3. Contains ≥ 1 Pro (P) within the internal segment
4. Contains ≥ 1 Asp (D) or Glu (E) within the internal segment
5. No K, R, or H inside the internal segment

### PEST fraction scoring

PEST fraction = (count of P + E + S + T residues) / segment length

| Score | Type tag |
|---|---|
| ≥ 0.30 | `pest_positive` |
| < 0.30 | `pest_putative` |

**Implementation note:** This uses a simplified PEST fraction (not the original
Rogers hydrophilicity formula). This mirrors the approach of EMBOSS epestfind.
The `max_length` parameter (default 100 aa) prevents spuriously long segments
spanning entire disordered proteins.

---

## Output Format — All Engines

Every engine returns dicts with this structure:

```python
{
    'start': int,   # 1-based, inclusive
    'end':   int,   # 1-based, inclusive
    'motif': str,   # matched sequence substring
    'type':  str,   # motif type label (see each engine above)
    'notes': str    # additional information; empty string if none
}
```

N-degron engine returns a **single dict or None** (not a list).

---

## Reference Sequences for Manual Verification

| Sequence | Expected hits |
|---|---|
| `KFERQAAA` | KFERQ (canonical, pos 1–5) |
| `AAAAILKEQ` | ILKEQ (canonical, pos 5–9) |
| `AAAVKKDQ` | VKKDQ (canonical, pos 4–8) |
| `QKFSR` | QKFSR (phospho_activated) |
| `KFERRXXX` | KFERR (acetyl_activated) |
| `DDDWTHLSS` | WTHLS → lir_xlir_acidic |
| `SSWTILAA` | WTIL → lir_xlir_phospho |
| `DSWTILAA` | WTIL → lir_xlir_optimal |
| `RTALGSSSN` | dbox_extended |
| `AAKENAAB` | ken_box at pos 3–5 |
| `RAAAGGG` | N-degron type_I_primary_destabilising |
| `MAAAGGG` | MAP cleavage → A at pos 2 → stabilising |
| `DAAAGGG` | N-degron secondary_destabilising |
| `AAAAEE` | cdegron_dcaf12 |
| `AAAARAAGG` | cdegron_appbp2 |
| `RSDEPSTSSTTPPEK` | pest_positive |

---

## Full Reference List

| Citation | DOI | Engine |
|---|---|---|
| Kirchner P et al. (2019) PLoS Biol 17(6):e3000301 | 10.1371/journal.pbio.3000301 | KFERQ |
| Rozenknop A et al. (2011) J Mol Biol 410:622–631 | 10.1016/j.jmb.2011.05.003 | LIR affinity |
| Cheng et al. (2023) preprint | 10.1101/2022.09.25.509395 | LIR core |
| Karin et al. (2014) figshare | 10.6084/m9.figshare.1601912.v1 | LIR core |
| McLaughlin et al. (2020) Biochem J 477:3115 | 10.1042/bcj20200714 | LIR core |
| Klionsky et al. (2014) iLIR/Autophagy | 10.4161/AUTO.28260 | xLIR/iLIR |
| Kohler et al. (2022, TEX264) preprint | 10.1101/2022.02.11.480038 | xLIR phospho |
| AF2 SLiM screen (2025) figshare | 10.6084/m9.figshare.28929437.v1 | LIR non-canonical |
| Glotzer M et al. (1991) Nature 349:132 | — | D-box minimal |
| Andrews PD et al. (2024) preprint | 10.1101/2024.04.30.590460 | D-box extended |
| Pfleger CM & Kirschner MW (2000) Genes Dev 14:655 | 10.1101/gad.14.6.655 | KEN box |
| Timms RT & Koren I (2020) Biochem Soc Trans 48:1355 | 10.1042/BST20191094 | N-degron, C-degron |
| Dissmeyer N et al. (2016) Curr Protein Pept Sci 17:4 | 10.2174/0929866523666160108115809 | N-degron |
| Rogers S et al. (1986) Science 234:364 | — | PEST |
| Rechsteiner M & Rogers SW (1996) TIBS 21:267 | — | PEST |

---

## Important Caveats

1. **Accessibility:** The scan is purely sequence-based. Buried motifs require
   structural data to assess functional relevance; this script does not assess
   solvent accessibility.
2. **D-box false positives:** `RxxL` is a very short, frequent pattern. Use
   `dbox_extended` (default) for confident predictions.
3. **PEST scoring:** The PEST fraction here differs from Rogers' original
   hydrophilicity-weighted score; results may differ from published PEST predictors.
4. **LIR context-dependence:** xLIR classification is based on upstream sequence
   context only; phosphorylation status is not directly assayed.
5. **N-degron MAP cleavage:** The MAP cleavage prediction is a rule-based heuristic
   (Met + small/uncharged penultimate residue). Actual cleavage depends on the
   MAP isoform and cellular context.
6. **Runtime dependencies:** This script runs fully offline. No external API calls
   are made at any point. All rules are encoded directly in the script.