# CMA Network Gene Lookup & Scoring — Agent Bootstrap Guide v1.0

## Purpose
This document gives a new agent everything it needs to use `cma_network.py`
**without reading any PDFs or making external calls**. Read this file, run the
self-test, and you are ready.

---

## Script Location (GitHub)

| File | Raw URL |
|---|---|
| `cma_network.py` | `https://raw.githubusercontent.com/BiologicExplorer/Tools/main/cma_network.py` |
| `CMA_NETWORK_KNOWLEDGE.md` | `https://raw.githubusercontent.com/BiologicExplorer/Tools/main/CMA_NETWORK_KNOWLEDGE.md` |

**See also:** `KFERQ_LIR_KNOWLEDGE.md` / `motif_finder.py` — companion tool for
sequence-level motif finding (KFERQ, LIR, D-box, KEN, N/C-degron, PEST).

**Bootstrap instruction for a fresh agent:**
> "Using ONLY the GitHub integration: get raw file content of `cma_network.py`
> from repository `BiologicExplorer/Tools`, save to `/home/sandbox/cma_network.py`,
> then run `python3 /home/sandbox/cma_network.py --self-test`"

---

## The Script

**Location:** `BiologicExplorer/Tools` → `cma_network.py`
**Version:** 1.0
**Self-test:**

```bash
python3 /home/sandbox/cma_network.py --self-test
# Expected: "ALL 35 TESTS PASSED — cma_network.py v1.0 READY"
```

---

## Scope — What This Tool Does

`cma_network.py` operates at the **gene / pathway level**. It answers two questions:

1. **Is this gene part of the CMA regulatory network?** (`check_cma_network_membership`)
2. **What is the net CMA activity level in this expression dataset?** (`calculate_cma_score`)

It does **not** scan protein sequences for motifs — use `motif_finder.py` for that.

| Tool | Input | What it finds | Level |
|---|---|---|---|
| `motif_finder.py` | Protein sequence | KFERQ, LIR, D-box, KEN, N/C-degron, PEST | Residue / motif |
| `cma_network.py` | Gene name or expression data | CMA network membership, CMA activation score | Gene / pathway |

---

## Biological Background

**Chaperone-Mediated Autophagy (CMA)** is a selective lysosomal degradation pathway
in which cytosolic chaperone HSC70 (HSPA8) recognises KFERQ-like motifs on substrate
proteins and escorts them to the lysosomal membrane receptor LAMP-2A for translocation.

The CMA regulatory network (Kirchner et al. 2019) is a set of genes whose expression
levels predict and modulate CMA activity. Genes with direction +1 promote CMA;
genes with direction -1 suppress it. LAMP-2A (LAMP2) is the rate-limiting step and
carries double weight in the activation score (Bourdenx et al. 2021).

---

## Python API

```python
import sys; sys.path.insert(0, '/home/sandbox')
from cma_network import (
    check_cma_network_membership,   # gene/UniProt → dict or None
    list_cma_network,               # → List[dict] of all 16 genes
    calculate_cma_score,            # expression_dict → score dict
    CMA_NETWORK,                    # the raw dict (16 genes)
)
```

---

## Function Reference

### `check_cma_network_membership(query)` → dict or None

Looks up a gene by symbol, alias, or human UniProt ID (case-insensitive).

**Returns** a dict if found, `None` if not in network:

```python
{
    "gene_symbol":    "LAMP2",
    "protein_name":   "Lysosome-associated membrane protein 2 (LAMP-2A isoform)",
    "aliases":        ["LAMP2", "LAMP2A", "LAMP-2A", "CD107B"],
    "uniprot_human":  "P13473",
    "category":       "effector",           # effector | lysosomal_modulator | extra_lysosomal_modulator
    "direction":      +1,                   # +1 activates CMA | -1 suppresses CMA
    "weight":         2,                    # LAMP2 = 2; all others = 1
    "notes":          "Rate-limiting receptor for CMA..."
}
```

**Examples:**

```python
check_cma_network_membership("LAMP2")   # → dict  (gene symbol)
check_cma_network_membership("LAMP2A")  # → dict  (alias)
check_cma_network_membership("P13473")  # → dict  (UniProt ID)
check_cma_network_membership("HSC70")   # → dict  (alias → resolves to HSPA8)
check_cma_network_membership("NRF2")    # → dict  (alias → resolves to NFE2L2)
check_cma_network_membership("PKB")     # → dict  (alias → resolves to AKT1)
check_cma_network_membership("TP53")    # → None  (not in network)
```

---

### `list_cma_network()` → List[dict]

Returns all 16 CMA network genes as a list of dicts (same structure as above).

```python
genes = list_cma_network()
len(genes)   # → 16
```

---

### `calculate_cma_score(expression_dict, reference_dict=None)` → dict

Implements the **Bourdenx 2021 CMA activation score**:

```
score = Σ(weight_i × direction_i × expr_i) / Σ(weight_i)
```

**Parameters:**
- `expression_dict` — `{"GENE": value, ...}` where values are **log-normalised**
  expression counts (e.g. log1p-normalised scRNA-seq or bulk RNA-seq counts).
  Accepts gene symbols, aliases, or UniProt IDs as keys.
  Non-network genes are silently ignored.
- `reference_dict` (optional) — healthy-control reference values with the same
  keys. When provided, each gene's contribution becomes
  `weight × direction × (expr - ref)`, i.e. fold-change from baseline.

**Returns:**

```python
{
    "cma_score":      float | None,   # None if no network genes were found
    "n_genes_scored": int,
    "genes_scored":   [str, ...],     # canonical gene symbols that contributed
    "warnings":       [str, ...]      # e.g. reference gene missing from expression
}
```

**Interpretation:**
| Score | Meaning |
|---|---|
| > 0 | Net CMA activation |
| < 0 | Net CMA suppression |
| ≈ 0 | Balanced / no signal |
| None | No CMA network genes in input |

LAMP2 (weight=2) dominates because it is the rate-limiting receptor.

**Examples:**

```python
# Basic score from absolute expression values
result = calculate_cma_score({
    "LAMP2": 5.2,    # activator, weight 2
    "HSPA8": 3.1,    # activator, weight 1
    "AKT1":  1.8,    # inhibitor, weight 1
})
# result["cma_score"] > 0 because LAMP2 (weight=2) dominates

# Fold-change score relative to healthy control
result = calculate_cma_score(
    {"LAMP2": 6.0, "HSPA8": 4.0, "AKT1": 2.0},
    reference_dict={"LAMP2": 3.0, "HSPA8": 3.0, "AKT1": 2.0}
)
# Each gene contributes (expr - ref), so only up/down-regulated genes add signal
```

---

## CLI Usage

```bash
# Check membership of a single gene
python3 /home/sandbox/cma_network.py --check LAMP2
python3 /home/sandbox/cma_network.py --check P13473     # UniProt ID
python3 /home/sandbox/cma_network.py --check HSC70      # alias

# List all 16 network genes as a table
python3 /home/sandbox/cma_network.py --list-network

# Calculate CMA score from JSON expression data
python3 /home/sandbox/cma_network.py --score '{"LAMP2": 5.2, "HSPA8": 3.1, "AKT1": 1.8}'

# With reference (fold-change scoring)
python3 /home/sandbox/cma_network.py \
    --score '{"LAMP2": 6.0, "HSPA8": 4.0, "AKT1": 2.0}' \
    --reference '{"LAMP2": 3.0, "HSPA8": 3.0, "AKT1": 2.0}'

# Run self-test suite
python3 /home/sandbox/cma_network.py --self-test
```

---

## CMA Network Gene List (Kirchner 2019, 16 genes)

Source: Kirchner P et al. (2019), Figure 4 + legend.
Directions cross-validated: Bourdenx M et al. (2021), Figure S10L immunoblot.

### Effectors — Core CMA Machinery

| Gene | Protein | Dir | Wt | UniProt | Notes |
|---|---|---|---|---|---|
| LAMP2 | LAMP-2A | +1 | **2** | P13473 | Rate-limiting receptor; LAMP-2A splice isoform only |
| HSPA8 | HSC70 | +1 | 1 | P11142 | Cytosolic chaperone; recognises KFERQ motifs |
| HSP90AA1 | HSP90α | +1 | 1 | P07900 | Co-chaperone; assists substrate delivery |
| DNAJB1 | HSP40 / DNAJB1 | +1 | 1 | P25685 | J-domain co-chaperone |

### Lysosomal Modulators

| Gene | Protein | Dir | Wt | UniProt | Notes |
|---|---|---|---|---|---|
| MT-RNR2 | Humanin | +1 | 1 | None* | Mitochondria-derived peptide; stabilises LAMP2A |
| GFAP | GFAP | +1 | 1 | P14136 | Glial fibrillary acidic protein |
| RAC1 | Rac1 | +1 | 1 | P63000 | GTPase; positive regulator |
| PHLPP1 | PHLPP1 | +1 | 1 | O60346[*] | PP2C phosphatase; dephosphorylates AKT |
| EEF1A1 | EF-1α | +1 | 1 | P68104 | Elongation factor; moonlights at lysosome |
| RICTOR | Rictor | -1 | 1 | Q6IA86[*] | mTORC2 component; phosphorylates AKT → inhibits CMA |
| AKT1 | AKT1/PKB | -1 | 1 | P31749 | Serine/threonine kinase; inhibits LAMP2A |
| CTSA | Cathepsin A | -1 | 1 | P10619 | Lysosomal protease; degrades LAMP2A |

### Extra-Lysosomal Modulators

| Gene | Protein | Dir | Wt | UniProt | Notes |
|---|---|---|---|---|---|
| NFATC1 | NFATc1 | +1 | 1 | O95644 | Transcription factor; upregulates CMA genes |
| NFE2L2 | NRF2 | +1 | 1 | Q16236 | Antioxidant TF; induces CMA under oxidative stress |
| RAB11A | Rab11A | -1 | 1 | P62491 | GTPase; promotes LAMP2A recycling away from surface |
| RARA | RARα | -1 | 1 | P10276 | Retinoic acid receptor; transcriptional repressor of CMA |

*MT-RNR2 (Humanin): mitochondrially encoded, no standard UniProt ID.
[*] Moderately confident UniProt assignments — verify at https://www.uniprot.org

---

## CMA Score Formula Detail

**Source:** Bourdenx M et al. (2021), Cell 184(10):2696–2714, PMC8152331.
Supplementary Methods + Figure S10L.

```
CMA_score = Σ(weight_i × direction_i × log_expr_i) / Σ(weight_i)
```

- Inputs must be **log-normalised** (e.g. log1p of TPM, CPM, or normalised scRNA-seq counts)
- LAMP2 weight = 2 (all others = 1) because LAMP-2A is the rate-limiting step
- Direction: +1 for activators, -1 for suppressors
- When `reference_dict` is provided:
  `expr_i` is replaced by `(sample_expr_i - reference_expr_i)` — fold-change from baseline

**Arithmetic example (3 genes):**
```
Input: LAMP2=5.2 (w=2,d=+1), HSPA8=3.1 (w=1,d=+1), AKT1=1.8 (w=1,d=-1)
score = (2×(+1)×5.2 + 1×(+1)×3.1 + 1×(-1)×1.8) / (2+1+1)
      = (10.4 + 3.1 - 1.8) / 4
      = 11.7 / 4
      = 2.925   → net CMA activation
```

---

## Reference List

| Citation | DOI | Used for |
|---|---|---|
| Kirchner P et al. (2019) PLoS Biol 17(6):e3000301 | 10.1371/journal.pbio.3000301 | Gene list, categories, directions |
| Bourdenx M et al. (2021) Cell 184(10):2696–2714 | 10.1016/j.cell.2021.03.048 | CMA score formula; LAMP2 weight=2 |

---

## Important Caveats

1. **Gene list size:** This implementation encodes 16 genes from Kirchner 2019 Fig 4.
   The full published CMA network (available as an Excel file from the Cuervo lab)
   may contain additional genes not yet included.
2. **Expression input units:** The score formula requires log-normalised values.
   Raw counts will give incorrect results.
3. **LAMP2 isoforms:** Only the LAMP-2A splice isoform is the CMA receptor. If your
   expression data distinguishes isoforms, use LAMP-2A-specific counts.
4. **Tissue specificity:** CMA network gene expression and CMA activity vary by cell
   type and stress context. The score is a relative indicator, not an absolute rate.
5. **No runtime calls:** This script is fully offline. No external API calls are made.