# Orchestrator — Bootstrap Context

## What This Document Is

Paste the URL of this file into a new chat to bootstrap it for Orchestrator development:

```
Read https://raw.githubusercontent.com/BiologicExplorer/Tools/main/orchestrator/ORCHESTRATOR_CONTEXT.md
and use it to bootstrap this chat for Orchestrator development.
```

---

## Project Goal

Build `orchestrator/orchestrator.py` — a single entry-point module that accepts a **viral protein sequence + host protein sequence**, routes them through all three existing modules, and returns a **unified, ranked output** describing autoimmune risk via the Super-epitope Architecture (SEA) mechanism.

The Orchestrator is a thin integration layer. It does **not** reimplement logic from the three modules — it calls them and merges their results.

---

## Mechanism Summary (SEA)

Viral protein → hinge phosphorylated → protein unfolds → CMA motif exposed → CMA routing → jammer resists full proteolysis → intact peptide loaded onto MHC-II → molecular mimicry → autoimmune disease.

The Orchestrator scores how well a viral/host sequence pair satisfies this pathway end-to-end.

---

## Repository Structure

```
BiologicExplorer/Tools        (GitHub, branch: main)
├── README.md
├── orchestrator/
│   ├── orchestrator.py       ← TO BUILD
│   └── ORCHESTRATOR_CONTEXT.md  (this file)
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
`summary` flags useful to Orchestrator:
- `cma_category` — `"canonical"`, `"phospho"`, `"acetyl"`, or `"none"`
- `n_kferq`, `n_lir`, `n_dbox` — motif counts
- `has_destabilising_n_term`, `has_destabilising_c_term`

All start/end positions are **1-based inclusive**.

Also available as individual functions: `find_kferq_motifs()`, `find_lir_motifs()`, `find_dbox_motifs()`, etc.

---

### 2. `cma_network/cma_network.py`

**Primary entry points:**
```python
check_cma_network_membership(query: str) -> Optional[Dict]
```
Looks up whether a gene/protein is a CMA network member. Returns `None` if not found, or:
```python
{
  'symbol'    : str,
  'category'  : str,   # e.g. "substrate", "regulator", "chaperone"
  'direction' : int,   # +1 = CMA activating, -1 = CMA inhibiting
  'weight'    : float,
  'aliases'   : List[str],
}
```

```python
calculate_cma_score(expression_dict: Dict[str, float],
                    reference_dict: Optional[Dict[str, float]] = None) -> Dict
```
Requires gene expression data. Not always available in a sequence-only workflow — design the Orchestrator to call it **optionally** when expression data is supplied.

Returns `{'cma_score': float or None, 'n_genes_scored': int, ...}`.

```python
list_cma_network() -> List[Dict]
```
Returns all CMA network members — useful for intersection checks.

---

### 3. `sea/sea_module.py`

**Primary entry point:**
```python
SEAModule(virus_name: str, config: Optional[SEAConfig] = None)
module.run(
    protein1: str,              # viral sequence (scanned inline)
    protein2: str,              # host sequence
    protein1_name: str,
    protein2_name: str,
    provided_motifs: Optional[List[str]] = None,  # host KFERQ motifs as strings
) -> List[SEAResult]
```
Returns a list of `SEAResult` objects (one per hinge site found in viral protein), each with:
```python
SEAResult:
  .position          : int
  .hinge_match       : HingeMatch
  .jammer_matches    : List[JammerMatch]
  .degradation_matches: List[DegradationMotifMatch]
  .architecture      : ArchitectureClass    # NONE, HINGE_ONLY, SUPER_EPITOPE, etc.
  .sea_score         : float
  .apc_multiplier    : float
  .final_score       : float
  .notes             : List[str]
```

`ArchitectureClass` values (ascending risk):
`NONE → HINGE_ONLY → JAMMER_ONLY → DEGRADATION_ONLY → HINGE_AND_DEGRADATION → JAMMER_AND_DEGRADATION → SANDWICHED → SUPER_EPITOPE`

**Key config values** (defaults in `SEAConfig`):
- `PROXIMITY_WINDOW = 15`, `JAMMER_DENSITY_WINDOW = 30`
- `JAMMER_BASE_SCORE = 0.8`, `JAMMER_DENSITY_WEIGHT = 0.12`
- Architecture score bonuses: SUPER_EPITOPE = 5.0, SANDWICHED = 3.0
- `APC_TROPISM` multipliers: EBV=2.0, CMV=1.8, HIV=1.8, HCV=1.3, default=1.0

**Providing host KFERQ motifs to SEAModule:**
Run `find_kferq_motifs(host_seq)` from `motif_finder` first, then pass the motif strings as `provided_motifs` to `module.run()`. This is the integration bridge between the two modules.

---

## Orchestrator Design Intent

### Inputs
```python
orchestrate(
    viral_seq:    str,
    host_seq:     str,
    virus_name:   str,
    viral_protein_name: str,
    host_protein_name:  str,
    expression_dict: Optional[Dict[str, float]] = None,  # for CMA score
    reference_dict:  Optional[Dict[str, float]] = None,
    sea_config:   Optional[SEAConfig] = None,
) -> OrchestratorResult
```

### Processing Pipeline
1. Run `find_all_degradation_motifs(viral_seq)` → viral degradation profile
2. Run `find_all_degradation_motifs(host_seq)` → host degradation profile
3. Extract `kferq` motif strings from host result → `provided_motifs`
4. Run `check_cma_network_membership(host_protein_name)` → host CMA status
5. Run `SEAModule.run(viral_seq, host_seq, ..., provided_motifs=...)` → SEA results
6. If `expression_dict` provided → run `calculate_cma_score()`
7. Merge into `OrchestratorResult` with ranked hits and a unified risk summary

### Output
`OrchestratorResult` should contain at minimum:
- `sea_results` — ranked `List[SEAResult]` (top hits)
- `viral_motif_profile` — output of `find_all_degradation_motifs` on viral sequence
- `host_motif_profile` — output of `find_all_degradation_motifs` on host sequence
- `host_cma_membership` — `check_cma_network_membership` result
- `cma_score` — optional float if expression data was supplied
- `risk_summary` — Dict with top architecture class, top final_score, n_super_epitope_hits, etc.
- A `report()` method that writes a structured `.md` file suitable for OneDrive output

### Output Policy
Code lives on GitHub. Reports and analysis outputs are written as `.md` files and pushed to OneDrive (not committed to GitHub). The `report()` method should produce a single well-structured `.md` file with consistent H1/H2/H3 headings so it drops cleanly into an Obsidian vault later.

---

## Verified Test Case (for regression)

**HCV NS3/NS4A (Q9WMX2) vs CYP2E1 (P05181)**
- 7 host KFERQ motifs found: pos 46, 143, 145, 146, 340, 354, 355
- 12 SEA pairs scored; top hit: E2 pos 675, SUPER_EPITOPE, final_score = 17.30
- Second hit: Core pos 110, SUPER_EPITOPE, final_score = 10.96

These values should be reproducible from the Orchestrator's output.

---

## What to Build First

1. Define `OrchestratorResult` dataclass
2. Implement `orchestrate()` following the pipeline above
3. Write `orchestrator/tests/orchestrator_test.py` — regression against the HCV/CYP2E1 case
4. Confirm top score ≥ 17.0, architecture = SUPER_EPITOPE

---

## Notes for the Agent

- Do **not** reimplement hinge/jammer/KFERQ logic — import and call the existing modules
- `sys.path` pattern for test files: `sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..'))`
  (tests live two levels below repo root at `orchestrator/tests/`)
- GitHub PAT is available in the user's environment — ask the user to provide it for push operations
- Mouse experiment is unrelated — do not mix with this work
- If expression data is unavailable, `cma_score` should be `None` and noted in the report; do not block the pipeline
