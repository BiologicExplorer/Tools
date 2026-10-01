# Orchestrator — Bootstrap Context

## What This Document Is

Paste the URL of this file into a new chat to bootstrap it for Orchestrator development:

```
Read https://raw.githubusercontent.com/BiologicExplorer/Tools/main/orchestrator/ORCHESTRATOR_CONTEXT.md
and use it to bootstrap this chat for Orchestrator development.
```

---

## Project Goal

`orchestrator/orchestrator.py` is built — a single entry-point module that accepts a
**viral protein sequence + host protein sequence**, routes them through all three existing
modules, and returns a **unified, ranked output** describing autoimmune risk via the
Super-epitope Architecture (SEA) mechanism.

The Orchestrator is a thin integration layer. It does **not** reimplement logic from the
three modules — it calls them and merges their results.

---

## Mechanism Summary (SEA)

Viral protein → hinge phosphorylated → protein unfolds → CMA motif exposed → CMA routing
→ jammer resists full proteolysis → intact peptide loaded onto MHC-II → molecular mimicry
→ autoimmune disease.

The Orchestrator scores how well a viral/host sequence pair satisfies this pathway
end-to-end.

---

## Repository Structure

```
BiologicExplorer/Tools        (GitHub, branch: main)
├── README.md
├── orchestrator/
│   ├── orchestrator.py
│   ├── mclachlan_aligner.py     ← gap-tolerant SW aligner (McLachlan 1972 matrix)
│   ├── protein_knowledge_base.json  ← CYP2E1 / Q9WMX2 epitope & CMA knowledge
│   ├── PROTEIN_BIOLOGY_CONTEXT.md
│   ├── ORCHESTRATOR_CONTEXT.md  (this file)
│   └── tests/
│       └── orchestrator_test.py
├── sea/
│   ├── sea_module.py
│   ├── SEA_PROJECT_CONTEXT.md
│   └── tests/
│       ├── sea_test.py
│       └── sea_real_test.py
├── motif_finder/
│   ├── motif_finder.py
│   └── KFERQ_LIR_KNOWLEDGE.md
└── cma_network/
    ├── cma_network.py
    └── CMA_NETWORK_KNOWLEDGE.md
```

Local workspace: `/home/sandbox/` mirrors this structure.

---

## The Three Modules — Verified APIs

### 1. `motif_finder/motif_finder.py`

**Primary entry point:**
```python
find_all_degradation_motifs(sequence: str, ...) -> Dict
```
Returns:
```python
{
  'kferq'    : List[Dict],   # {start, end, motif, type, notes}  — 1-based inclusive
  'lir'      : List[Dict],
  'dbox'     : List[Dict],
  'ken'      : List[Dict],
  'n_degron' : Dict,
  'c_degron' : List[Dict],
  'pest'     : List[Dict],
  'summary'  : Dict,         # per-engine counts and boolean flags
}
```
`summary` flags:
- `cma_category` — `"canonical"`, `"phospho"`, `"acetyl"`, or `"none"`
- `n_kferq`, `n_lir`, `n_dbox` — motif counts
- `has_destabilising_n_term`, `has_destabilising_c_term`

**SEA integration helper (added):**
```python
kferq_to_sea_motifs(
    kferq_list: List[Dict],   # output of find_kferq_motifs() or profile['kferq']
    protein:    str = 'protein2',
) -> List[Dict]               # [{motif, position (0-based), protein}, ...]
```
Converts motif_finder's 1-based positions to the 0-based format expected by SEAModule.
Pass `protein='protein2'` for host motifs (signals host CMA exposure risk to SEA).

---

### 2. `cma_network/cma_network.py`

**Primary entry points:**
```python
check_cma_network_membership(query: str) -> Optional[Dict]
```
Query by gene symbol, alias, or UniProt accession (case-insensitive).
Returns `None` if not found, or:
```python
{
  'gene_symbol'   : str,
  'protein_name'  : str,
  'aliases'       : List[str],
  'uniprot_human' : str or None,
  'category'      : str,   # 'effector' | 'lysosomal_modulator' | 'extralysosomal_modulator'
  'direction'     : int,   # +1 = CMA activating, -1 = CMA inhibiting
  'weight'        : int,   # 2 for LAMP2A; 1 for all others
  'notes'         : str,
  'in_network'    : True
}
```

```python
calculate_cma_score(expression_dict: Dict[str, float],
                    reference_dict: Optional[Dict[str, float]] = None) -> Dict
```
Requires gene expression data. The Orchestrator calls this **only** when
`expression_dict` is supplied; otherwise `cma_score` is `None`.

---

### 3. `sea/sea_module.py`

#### `find_homologous_pairs()` — module-level function (added)

```python
find_homologous_pairs(
    viral_seq:    str,
    host_seq:     str,
    min_identity: float = 0.33,
    window:       int   = 12,
    step:         int   = 4,
    max_pairs:    int   = 40,
) -> List[Dict]
```
Slides a *window*-residue window over *host_seq* (stride *step*), locally aligns each
window against the full *viral_seq* with BioPython PairwiseAligner, and keeps hits with
ungapped identity ≥ *min_identity* and aligned block ≥ 6 residues. Deduplicates by
viral position (≤50% overlap). Returns pairs sorted by similarity descending, 1-based rank.

Each pair dict:
```python
{
  'seq1'            : str,    # viral fragment
  'seq2'            : str,    # host fragment
  'position1'       : int,    # 0-based start in viral_seq
  'position2'       : int,    # 0-based start in host_seq
  'protein1'        : str,    # full viral_seq (needed by SEAScorer)
  'similarity_score': float,  # ungapped identity 0–1
  'rank'            : int,    # 1-based
}
```
Requires BioPython (`pip install biopython`).

#### `SEAModule` — primary class

```python
SEAModule(virus_name: str = "default", config: Optional[SEAConfig] = None)
```

**`run()` — low-level entry point (takes pre-built pairs):**
```python
module.run(
    homologous_pairs:    List[Dict],   # output of find_homologous_pairs()
    autoimmune_diseases: List[str],    # candidate disease names (stored for downstream)
    degradation_motifs:  List[Dict],   # [{motif, position (0-based), protein}, ...]
) -> List[SEAResult]
```

**`run_from_sequences()` — convenience entry point (added):**
```python
module.run_from_sequences(
    viral_seq:           str,
    host_seq:            str,
    provided_motifs:     Optional[List[Dict]] = None,   # SEA-format motif dicts
    autoimmune_diseases: Optional[List[str]]  = None,
    min_identity:        float = 0.33,
    window:              int   = 12,
    step:                int   = 4,
    max_pairs:           int   = 40,
) -> List[SEAResult]
```
Calls `find_homologous_pairs()` then `run()` — the standard Orchestrator entry point.

#### `SEAResult` — verified field names

```python
SEAResult:
  .rank               : int
  .seq1               : str          # viral fragment
  .seq2               : str          # host fragment
  .position1          : int          # 0-based start in viral protein
  .position2          : int          # 0-based start in host protein
  .hinges             : List[HingeMatch]
  .jammers            : List[JammerMatch]
  .degradation_motifs : List[DegradationMotifMatch]
  .base_score         : float
  .hinge_score        : float
  .jammer_density_score: float
  .degradation_score  : float
  .architecture_bonus : float
  .cell_type_weight   : float        # APC-tropism multiplier
  .final_sea_score    : float        # primary sort key
  .architecture_class : ArchitectureClass
  .is_sandwiched      : bool
  .is_super_epitope   : bool
  .notes              : List[str]
```

#### `ArchitectureClass` — actual enum values (ascending risk)
```
NONE(0) → DEGRADATION_ONLY(1) → HINGE_ONLY(2) → JAMMER_ONLY(3)
→ HINGE_AND_DEGRADATION(4) → JAMMER_AND_DEGRADATION(5) → SANDWICHED(6) → SUPER_EPITOPE(7)
```

**Key config values** (defaults in `SEAConfig`):
- `PROXIMITY_WINDOW = 15`, `JAMMER_DENSITY_WINDOW = 30`
- Architecture bonuses: `SUPER_EPITOPE = 5.0`, `SANDWICHED = 3.0`
- `APC_TROPISM`: EBV=2.0, CMV=1.8, HIV=1.8, HCV=1.3, default=1.0

**`provided_motifs` format** (critical — differs from old context doc):
```python
# Correct format (0-based position):
[{'motif': 'QLELK', 'position': 46, 'protein': 'protein2'}, ...]

# Use kferq_to_sea_motifs() from motif_finder to produce this from motif_finder output.
```

---

## Orchestrator API

```python
orchestrate(
    viral_seq:              str,
    host_seq:               str,
    virus_name:             str,
    viral_protein_name:     str,
    host_protein_name:      str,
    host_accession:         Optional[str]           = None,
    expression_dict:        Optional[Dict[str, float]] = None,
    reference_dict:         Optional[Dict[str, float]] = None,
    sea_config:             Optional[SEAConfig]     = None,
    # ── pair source A: BioPython seq aligner ─────────────────────────────
    min_identity:           float = 0.33,
    window:                 int   = 12,
    step:                   int   = 4,
    max_pairs:              int   = 40,
    # ── pair source B: McLachlan v3 gapless hits ──────────────────────────
    mclachlan_hits:         Optional[List[Dict]]    = None,
    mclachlan_min_composite: float = 15.0,
    mclachlan_max_pairs:    int   = 200,
    # ── pair source C: Smith-Waterman gap-tolerant hits ───────────────────
    sw_hits:                Optional[List[Dict]]    = None,
    sw_min_composite:       float = 14.5,
    sw_max_pairs:           int   = 200,
) -> OrchestratorResult
```

### Pair Sources
The orchestrator accepts up to three independent pair sources and merges + deduplicates them before SEA scoring:

| Source | Parameter | Description |
|--------|-----------|-------------|
| A: seq_aligner | *(none — always runs)* | BioPython pairwise2 sliding-window aligner (requires `biopython`). Returns 0 if BioPython not installed. |
| B: mclachlan | `mclachlan_hits` | Pre-computed McLachlan v3 gapless 7-mer hit dicts from `cross_protein_comparator.py`. Pass `None` to skip. |
| C: sw | `sw_hits` | Pre-computed Smith-Waterman gap-tolerant hit dicts from `mclachlan_aligner.py`. Same core schema as McLachlan hits; extra `_has_gaps`, `_sw_score`, `_v_span`, `_h_span` keys are safely ignored. Pass `None` to skip. |

**SW hit schema** (produced by `mclachlan_aligner.hits_from_sw()`):
```python
{
  "motif":             str,   # viral fragment (aligned, no gaps in motif)
  "p1_pos1":           str,   # "start-end" 1-based, viral
  "p2_window":         str,   # host fragment (may contain gap chars)
  "p2_pos":            str,   # "start-end" 1-based, host
  "composite_primary": float, # McLachlan-equivalent score (threshold 14.5)
  "layer2_cross_mean": float, # mean cross-residue score (→ similarity_score)
  "_aln1":             str,   # gapped alignment string, viral strand
  "_aln2":             str,   # gapped alignment string, host strand
  "_n_matched":        int,   # number of matched residue pairs
  "_sw_score":         float, # raw SW matrix score
  "_has_gaps":         bool,  # True if alignment contains indels
  "_v_span":           int,   # residues spanned in viral sequence
  "_h_span":           int,   # residues spanned in host sequence
}
```

### Processing Pipeline
1. `find_all_degradation_motifs(viral_seq)` → viral degradation profile
2. `find_all_degradation_motifs(host_seq)` → host degradation profile
3. `kferq_to_sea_motifs(host_profile['kferq'], protein='protein2')` → SEA motif dicts
4. `check_cma_network_membership(host_protein_name)` → host CMA status
5a. Source A: `find_homologous_pairs()` → seq_aligner pairs
5b. Source B: `mclachlan_to_pairs(mclachlan_hits, source="mclachlan")` → gapless pairs
5c. Source C: `mclachlan_to_pairs(sw_hits, source="sw")` → gap-tolerant pairs
6. `_dedup_pairs(tolerance=4)` → merged, deduplicated, re-ranked pairs
7. `SEAModule.run(merged_pairs, ...)` → SEA results
7b. `apply_epitope_proximity_bonus()` → knowledge-base bonus for known-epitope hits
8. *(optional)* `calculate_cma_score(expression_dict)` if expression data supplied
9. Build `OrchestratorResult` + `risk_summary`

### `OrchestratorResult`
```python
OrchestratorResult:
  .sea_results          : List[SEAResult]    # sorted by final_sea_score desc
  .phase2_enriched_pairs: List[Dict]         # pairs with proximity event fields
  .viral_motif_profile  : Dict               # find_all_degradation_motifs on viral seq
  .host_motif_profile   : Dict               # find_all_degradation_motifs on host seq
  .host_cma_membership  : Optional[Dict]     # check_cma_network_membership result
  .cma_score            : Optional[float]    # None if no expression data
  .pair_source_counts   : Dict               # see keys below
  .risk_summary         : Dict               # see keys below
  .report(output_path)  : None               # writes .md file

pair_source_counts keys:
  seq_aligner   int   pairs from BioPython aligner (0 if BioPython absent)
  mclachlan     int   pairs from McLachlan gapless hits
  sw            int   pairs from Smith-Waterman gap-tolerant hits (0 if sw_hits=None)
  merged        int   total before dedup
  after_dedup   int   total after _dedup_pairs(tolerance=4)

risk_summary keys:
  top_architecture, top_sea_score, n_super_epitope, n_sandwiched,
  n_pairs_scored, host_cma_member, host_cma_category, host_cma_direction,
  cma_score, cma_score_display, n_viral_kferq, n_host_kferq,
  viral_cma_category, host_cma_category_motif
```

### Output Policy
Code lives on GitHub. Reports are written as `.md` files via `result.report(path)` and
pushed to OneDrive (not committed to GitHub). The report uses consistent H1/H2/H3
headings suitable for an Obsidian vault.

---

## Verified Regression Tests (HCV NS3/NS4A vs CYP2E1)

**Sequences:** HCV polyprotein Q9WMX2 (3,010 aa) vs CYP2E1 P05181 (493 aa)

**Test suite:** 102 tests, 1 known skip (BioPython not installed), 0 failures.

### McLachlan-only pipeline (no SW source)
| Metric | Expected |
|--------|----------|
| mclachlan pairs (threshold 15.0) | 134 |
| after_dedup | 76 |
| Top final_sea_score | ≥ 17.0 |
| Overall risk | CRITICAL |

### Full pipeline with SW source
| Metric | Expected |
|--------|----------|
| SW hits (threshold 14.5) | 128 (66 gapped-new vs v3 gapless scan) |
| SW hit score range | 30–42 composite_primary |
| sw key in pair_source_counts | present, > 0 when sw_hits supplied |

### Known-epitope calibration
| Metric | Expected |
|--------|----------|
| JHDN-5 range (CYP2E1 Gly113–Leu135) | 3 McLachlan hits in range |
| Host KFERQ motifs (canonical) | 5 (pos 47, 123, 158, 341, 358) |

> **Note on KFERQ count:** The original context doc cited 7 motifs from a simplified
> inline scanner. `motif_finder.find_kferq_motifs()` applies canonical KFERQ biochemistry
> (Q-flanked with strict core pattern) and returns 5 validated motifs — the more
> accurate count for CMA risk assessment.

Run regression:
```bash
cd /home/sandbox
python -m unittest orchestrator.tests.orchestrator_test -v
```

---

## Notes for the Agent

- Do **not** reimplement hinge/jammer/KFERQ logic — import and call the existing modules
- `sys.path` pattern for test files: `sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..'))`
- GitHub PAT is available in the user's environment — ask the user to provide it for push operations
- If expression data is unavailable, `cma_score` is `None` and noted in the report; does not block the pipeline
- Mouse experiment is unrelated — do not mix with this work
- BioPython is required for `find_homologous_pairs()`: `pip install biopython`
