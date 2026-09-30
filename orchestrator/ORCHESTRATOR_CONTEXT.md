# Orchestrator — Bootstrap Context

## What This Document Is

Paste the URL of this file into a new chat to bootstrap it for Orchestrator
development or extension:

```
Read https://raw.githubusercontent.com/BiologicExplorer/Tools/main/orchestrator/ORCHESTRATOR_CONTEXT.md
and use it to bootstrap this chat for Orchestrator development.
```

---

## Status

**COMPLETE and tested.** `orchestrator/orchestrator.py` is built, passing 30/30
regression tests, and live on GitHub (commits `51484d03` / `07de516f`).

This document describes the system **as-built**.

---

## Project Goal

`orchestrator/orchestrator.py` is a single entry-point module that accepts a
**viral protein sequence + host protein sequence**, routes them through all three
existing SEA modules, and returns a unified, ranked output describing autoimmune
risk via the Super-epitope Architecture (SEA) mechanism.

The Orchestrator is a thin integration layer.  It does **not** reimplement logic
from the three core modules — it calls them and merges their results.

---

## Mechanism Summary (SEA)

Viral protein → hinge phosphorylated → protein unfolds → CMA motif exposed →
CMA routing → jammer resists full proteolysis → intact peptide loaded onto
MHC-II → molecular mimicry → autoimmune disease.

The Orchestrator scores how well a viral/host sequence pair satisfies this
pathway end-to-end.

---

## Repository Structure

```
BiologicExplorer/Tools        (GitHub, branch: main)
├── README.md
├── orchestrator/
│   ├── orchestrator.py              ← BUILT
│   ├── ORCHESTRATOR_CONTEXT.md      (this file)
│   └── tests/
│       ├── __init__.py
│       └── orchestrator_test.py     ← BUILT (30/30 PASS)
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

## The Three Modules — APIs the Orchestrator Calls

### 1. `motif_finder/motif_finder.py`

**Primary entry point:**
```python
find_all_degradation_motifs(sequence: str, ...) -> Dict
```
Returns a unified dict:
```python
{
  'kferq'    : List[Dict],   # KFERQ-like CMA motifs  {start, end, motif, type, notes}
  'lir'      : List[Dict],   # LIR/AIM autophagy motifs
  'dbox'     : List[Dict],   # D-box degrons
  'ken'      : List[Dict],   # KEN box degrons
  'n_degron' : Dict,         # N-terminal residue classification
  'c_degron' : List[Dict],   # C-terminal degrons
  'pest'     : List[Dict],   # PEST rapid-turnover sequences
  'summary'  : Dict,         # per-engine counts and boolean flags
}
```
`summary` flags used by the Orchestrator:
- `cma_category` — `"canonical"`, `"phospho"`, `"acetyl"`, or `"none"`
- `n_kferq`, `n_lir`, `n_dbox` — motif counts
- `has_destabilising_n_term`, `has_destabilising_c_term`

All start/end positions are **1-based inclusive**.

Also used individually: `find_kferq_motifs(sequence)` → `List[Dict]`

---

### 2. `cma_network/cma_network.py`

**Primary entry points:**
```python
check_cma_network_membership(query: str) -> Optional[Dict]
```
Returns `None` if not found, or:
```python
{
  'symbol'    : str,
  'category'  : str,   # "substrate", "regulator", "chaperone"
  'direction' : int,   # +1 = CMA activating, -1 = CMA inhibiting
  'weight'    : float,
  'aliases'   : List[str],
}
```

```python
calculate_cma_score(expression_dict: Dict[str, float],
                    reference_dict: Optional[Dict[str, float]] = None) -> Dict
```
Called **optionally** when expression data is supplied.
Returns `{'cma_score': float or None, 'n_genes_scored': int, ...}`.

---

### 3. `sea/sea_module.py`

**Three functions used by the Orchestrator:**

```python
# 1. Sequence-aligner pair source
find_homologous_pairs(
    viral_seq:    str,
    host_seq:     str,
    min_identity: float = 0.33,
    window:       int   = 12,
    step:         int   = 4,
    max_pairs:    int   = 40,
) -> List[Dict]
```
Returns pair dicts tagged `source='seq_aligner'` (added by Orchestrator).

```python
# 2. SEA scorer
SEAModule(virus_name: str, config: Optional[SEAConfig] = None)
module.run(
    homologous_pairs:    List[Dict],   # merged pairs from both sources
    autoimmune_diseases: List,         # pass [] if unknown
    degradation_motifs:  List[Dict],   # host KFERQ motifs in SEA schema
) -> List[SEAResult]
```

**CRITICAL — SEAModule.run() pair dict schema.**
`SEAModule.run()` reads only these keys from each pair dict:
`protein1`, `seq1`, `position1`, `position2`, `rank`, `similarity_score`.
The `protein2` field is NOT used by the scorer.

`SEAResult` fields:
```python
SEAResult:
  .seq1                : str
  .seq2                : str
  .position1           : int     # 0-based, viral
  .position2           : int     # 0-based, host
  .hinges              : list
  .jammers             : list
  .degradation_motifs  : list
  .base_score          : float
  .hinge_score         : float
  .jammer_score        : float
  .jammer_density_score: float
  .degradation_score   : float
  .architecture_bonus  : float
  .cell_type_weight    : float
  .final_sea_score     : float
  .architecture_class  : ArchitectureClass   (enum; .name → str)
  .is_sandwiched       : bool
  .is_super_epitope    : bool
  .notes               : List[str]
```

`ArchitectureClass` values (ascending risk):
`NONE → HINGE_ONLY → JAMMER_ONLY → DEGRADATION_ONLY →
 HINGE_AND_DEGRADATION → JAMMER_AND_DEGRADATION → SANDWICHED → SUPER_EPITOPE`

**Key SEAConfig defaults:**
- `PROXIMITY_WINDOW = 15`, `JAMMER_DENSITY_WINDOW = 30`
- Architecture score bonuses: `SUPER_EPITOPE = 5.0`, `SANDWICHED = 3.0`
- `APC_TROPISM` multipliers: `EBV=2.0`, `CMV=1.8`, `HIV=1.8`, `HCV=1.3`, default=1.0

```python
# 3. Host KFERQ motif bridge
# find_kferq_motifs returns {start (1-BASED), end, motif, type, notes}
# Convert to SEA schema before passing to module.run():
provided_motifs = [
    {"motif": m["motif"], "position": m["start"] - 1, "protein": "protein2"}
    for m in find_kferq_motifs(host_seq)
]
```

---

## The Two Pair Sources

### Source A — Sequence Aligner (`find_homologous_pairs`)

Sliding-window identity scan over the two full sequences.
Controlled by `seq_aligner_min_identity`, `seq_aligner_window`,
`seq_aligner_step`, `seq_aligner_max_pairs`.

Pairs are tagged `source='seq_aligner'` (by the Orchestrator, not the function).

### Source B — McLachlan (`mclachlan_to_pairs`)

Pre-computed cross-protein comparator output (`Q9WMX2_vs_P05181_scored_v3.json`)
filtered to `composite_primary ≥ min_composite` and converted to pair dict schema.

Pass the `"hits"` list from the v3 JSON as `mclachlan_hits`.  Pass `None` to
skip this source entirely.

**Normalization:**
```
similarity_score = min(layer2_cross_mean / 4.0, 1.0)
```
McLachlan positions (`p1_pos1`, `p2_pos`) are **1-based** in the JSON.
`mclachlan_to_pairs()` subtracts 1 to produce **0-based** positions for SEA.

---

## Orchestrator API — `orchestrate()`

```python
from orchestrator.orchestrator import orchestrate

result = orchestrate(
    viral_seq                = str,          # full viral protein sequence
    host_seq                 = str,          # full host protein sequence
    virus_name               = str,          # e.g. "HCV" — sets APC tropism weight
    viral_protein_name       = str,          # e.g. "polyprotein Q9WMX2"
    host_protein_name        = str,          # e.g. "CYP2E1 P05181"
    # McLachlan pair source
    mclachlan_hits           = Optional[List[Dict]],   # hits list from v3 JSON; None = skip
    mclachlan_min_composite  = float,        # default 15.0
    mclachlan_max_pairs      = int,          # default 200
    # Sequence aligner pair source
    seq_aligner_min_identity = float,        # default 0.33
    seq_aligner_window       = int,          # default 12
    seq_aligner_step         = int,          # default 4
    seq_aligner_max_pairs    = int,          # default 40
    # Optional expression data
    expression_dict          = Optional[Dict[str, float]],
    reference_dict           = Optional[Dict[str, float]],
    # SEA config override
    sea_config               = Optional[SEAConfig],
) -> OrchestratorResult
```

---

## Pipeline (9 steps)

```
Step 1   find_all_degradation_motifs(viral_seq)        → viral_profile
Step 2   find_all_degradation_motifs(host_seq)         → host_profile
Step 3   find_kferq_motifs(host_seq) → convert 1-based → 0-based SEA schema
Step 4   check_cma_network_membership(host_protein_name)
Step 5a  find_homologous_pairs(...)                    → aligner_pairs (Source A)
Step 5b  mclachlan_to_pairs(mclachlan_hits, ...)       → mc_pairs      (Source B)
Step 6   aligner_pairs + mc_pairs → _dedup_pairs(tolerance=4) → merged_after_dedup
          Re-assign consecutive 1-based ranks sorted by similarity_score DESC
Step 7   SEAModule.run(merged_after_dedup, [], provided_motifs) → sea_results
Step 8   calculate_cma_score() — optional, skipped if expression_dict is None
Step 9   _build_risk_summary() → risk_summary
         Return OrchestratorResult
```

---

## `OrchestratorResult` Dataclass

```python
@dataclass
class OrchestratorResult:
    # Identity
    virus_name:          str
    viral_protein_name:  str
    host_protein_name:   str

    # Core outputs
    sea_results:         List[SEAResult]     # sorted by final_sea_score DESC
    viral_motif_profile: Dict
    host_motif_profile:  Dict
    host_cma_membership: Optional[Dict]
    cma_score:           Optional[float]

    # Summary
    risk_summary:        Dict
    pair_source_counts:  Dict[str, int]      # keys: seq_aligner, mclachlan, merged, after_dedup

    # Meta
    timestamp:           str                 # ISO-8601 UTC

    # Method
    def report(self, output_path: str) -> str
    # Writes structured .md to output_path; returns absolute path.
    # H1/H2/H3 headings — drops cleanly into Obsidian vault.
```

**`risk_summary` keys:**
`top_architecture`, `top_final_score`, `n_super_epitope`, `n_sandwiched`,
`n_pairs_scored`, `host_is_cma_member`, `host_cma_category`,
`viral_n_kferq`, `host_n_kferq`, `overall_risk`, `cma_score`.

**`overall_risk` tiers:**
| Condition | Risk |
|-----------|------|
| ≥1 SUPER_EPITOPE **and** top score ≥ 15.0 | `CRITICAL` |
| ≥1 SUPER_EPITOPE **or** top score ≥ 10.0 | `HIGH` |
| ≥1 SANDWICHED **or** top score ≥ 5.0 | `MODERATE` |
| Otherwise | `LOW` |
| Host is known CMA substrate/regulator and tier ≤ MODERATE | elevated to `HIGH` |

---

## `mclachlan_to_pairs()` — Public Function

```python
mclachlan_to_pairs(
    mclachlan_hits:  List[Dict],   # 'hits' list from scored_v3.json
    viral_seq:       str,          # full viral sequence (stored in protein1)
    min_composite:   float = 15.0,
    max_pairs:       int   = 200,
) -> List[Dict]
```

Each output pair dict:
```python
{
    "seq1":              str,    # viral fragment (hit["motif"])
    "seq2":              str,    # host fragment  (hit["p2_window"])
    "position1":         int,    # 0-based viral start
    "position2":         int,    # 0-based host start
    "protein1":          str,    # full viral_seq
    "similarity_score":  float,  # min(layer2_cross_mean / 4.0, 1.0)
    "rank":              int,    # 0 (re-assigned by Orchestrator after dedup)
    "source":            str,    # "mclachlan"
    "_composite_primary": float, # diagnostic carry-through
    "_layer2_cross_mean": float, # diagnostic carry-through
}
```

---

## `_dedup_pairs()` — Internal Function

```python
_dedup_pairs(pairs: List[Dict], tolerance: int = 4) -> List[Dict]
```

Removes pairs where BOTH `position1` and `position2` are within `tolerance`
residues of an already-accepted pair.  Resolves duplicates by keeping the
higher `similarity_score`.

Strategy: sort descending by `similarity_score`, then greedily accept.

**Locked-in default:** `tolerance = 4` (4 residues in both coordinates).

---

## Design Decisions (Locked)

| Decision | Value | Rationale |
|----------|-------|-----------|
| McLachlan position encoding | 1-based in JSON (`p1_pos1`) | `cross_protein_comparator.py` line 180 |
| SEA position encoding | 0-based | `SEAModule.run()` requirement |
| Conversion | subtract 1 from McLachlan start | in `mclachlan_to_pairs()` |
| Similarity normalization | `min(layer2_cross_mean / 4.0, 1.0)` | McLachlan scale is 0–4 |
| Dedup tolerance | 4 residues in both positions | avoids near-duplicate inflation |
| KFERQ position conversion | `start - 1` (1-based → 0-based) | `find_kferq_motifs()` is 1-based |
| `protein2` in pair dict | NOT used by `SEAModule.run()` | scorer reads only `protein1`/`seq1`/`position1`/`position2`/`rank`/`similarity_score` |
| CMA score | optional; `None` if no expression data | does not block pipeline |
| `seq_aligner` failure | caught, returns `[]` | `ImportError` fallback |

---

## Regression Test Suite — `orchestrator/tests/orchestrator_test.py`

**Result: 30/30 PASS (1.866 s)**

Four test classes:

| Class | Tests | What it validates |
|-------|-------|-------------------|
| `TestMcLachlanToPairs` | 6 | Unit: filter, positions, score range, cap, count=134, `protein1` field |
| `TestDedupPairs` | 5 | Unit: identical, within-tolerance, outside-tolerance, cross-source, empty |
| `TestOrchestratorRegression` | 14 | Integration: HCV/CYP2E1 full run with both sources |
| `TestOrchestratorMcLachlanOnly` | 3 | Integration: seq_aligner blocked (min_identity=0.999) |
| `TestOrchestratorSeqAlignerOnly` | 2 | Integration: mclachlan_hits=None |

**Run:**
```bash
cd /home/sandbox
python orchestrator/tests/orchestrator_test.py
# or
pytest orchestrator/tests/orchestrator_test.py -v
```

**Path setup pattern** (for all test files under `orchestrator/tests/`):
```python
_REPO_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
```

**Data loading pattern** (avoid repeated UniProt API hits):
```python
# Use module-level cache; sequences fetched from UniProt once per session
_VIRAL_SEQ = None
_HOST_SEQ  = None

def _get_test_data():
    global _VIRAL_SEQ, _HOST_SEQ
    if _VIRAL_SEQ is None:
        from protein_degradation.repeat_finder.sequence_loader import load_from_uniprot
        _, _VIRAL_SEQ = load_from_uniprot("Q9WMX2")   # [1] for sequence, not header
        _, _HOST_SEQ  = load_from_uniprot("P05181")
    return _VIRAL_SEQ, _HOST_SEQ
```

**NOTE:** `load_from_uniprot(id)` returns a `(header, sequence)` **tuple**.
Always take `[1]` for the sequence.  The `p1_header`/`p2_header` keys in the
v3 JSON are plain description strings, NOT FASTA — do not pass them as sequences.

---

## Verified Regression Targets (HCV Q9WMX2 vs CYP2E1 P05181)

| Metric | Value |
|--------|-------|
| `seq_aligner` pairs | 12 |
| `mclachlan` pairs (threshold 15.0) | 134 |
| `merged` pairs | 146 |
| `after_dedup` pairs | 72 |
| Top `final_sea_score` | **17.2198** |
| Top architecture | **SUPER_EPITOPE** |
| Total SUPER_EPITOPE hits | **12** |
| Total SANDWICHED hits | **14** |
| `overall_risk` | **CRITICAL** |
| HCV APC tropism weight | 1.3 |

Test assertions (must hold in regression):
- `top_final_sea_score ≥ 17.0`
- `top_architecture == "SUPER_EPITOPE"`
- `n_super_epitope ≥ 1`
- `mclachlan pair count ≥ 100`
- `after_dedup ≤ merged`
- `host KFERQ motifs ≥ 5`
- `overall_risk != "LOW"`

---

## Output Policy

Code lives on GitHub.  Reports and analysis outputs are written as `.md` files
and may be pushed to OneDrive — **not committed to GitHub**.  The `report()`
method produces a single well-structured `.md` file with consistent H1/H2/H3
headings so it drops cleanly into an Obsidian vault.

---

## Notes for the Agent

- Do **not** reimplement hinge/jammer/KFERQ logic — import and call the existing modules.
- `mclachlan_hits` must be the `"hits"` list extracted from the v3 JSON, not the whole file.
- `SEAModule.run()` does NOT use `protein2` from the pair dict — only `protein1`, `seq1`,
  `position1`, `position2`, `rank`, `similarity_score` are read.
- `find_kferq_motifs()` positions are **1-based**; subtract 1 when building `provided_motifs`.
- McLachlan `p1_pos1` / `p2_pos` are **1-based** strings like `"384-390"`;
  `mclachlan_to_pairs()` converts them to 0-based integers.
- If expression data is unavailable, `cma_score` is `None` — this does not block the pipeline.
- `_dedup_pairs()` is an internal helper (prefixed `_`); do not expose in new public APIs.
- GitHub PAT is available in the user's environment — ask the user to provide it for push operations.
- Mouse experiment is unrelated — do not mix with this work.
