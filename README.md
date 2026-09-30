# BiologicExplorer / Tools

A collection of independent bioinformatics modules and experimental tools.
Each project lives in its own directory. Tests are under `<project>/tests/` and
can be run or deleted independently.

---

## Projects

### `sea/` — Super-epitope Architecture (SEA) Scanner
Scores homologous viral/host sequence pairs for autoimmune disease risk by
detecting coordinated structural features: phosphorylatable hinges (SS/ST/TT),
lysosomal-resistance jammers (GA/GR/GK), and CMA degradation motifs (KFERQ-like).

- **Module:** `sea/sea_module.py`
- **Bootstrap context:** `sea/SEA_PROJECT_CONTEXT.md`
- **Tests:** `sea/tests/sea_test.py` (regression), `sea/tests/sea_real_test.py` (HCV vs CYP2E1)
- **Run tests:** `python sea/tests/sea_test.py`

### `motif_finder/` — KFERQ / LIR Motif Finder
Scans protein sequences for CMA targeting motifs (KFERQ-like) and
LC3-interaction region (LIR) motifs involved in selective autophagy.

- **Module:** `motif_finder/motif_finder.py`
- **Knowledge base:** `motif_finder/KFERQ_LIR_KNOWLEDGE.md`

### `cma_network/` — CMA Network Module
Models the chaperone-mediated autophagy (CMA) pathway network, including
LAMP-2A substrate interactions and pathway activity scoring.

- **Module:** `cma_network/cma_network.py`
- **Knowledge base:** `cma_network/CMA_NETWORK_KNOWLEDGE.md`

---

## Architecture

These three modules feed a planned **Orchestrator** that will combine ranked
homologous pairs from all upstream tools and call `SEAModule.run()` for final
autoimmune risk scoring.

```
motif_finder  ──┐
                ├──► Orchestrator ──► SEAModule ──► ranked risk output
cma_network   ──┘
```

---

## Output

Analysis outputs, reports, and documents written about these systems are stored
in OneDrive (not in this repository).

---

## Adding a New Project

1. Create a top-level directory: `mkdir <project_name>/`
2. Add a `tests/` subdirectory for any test files
3. Add a brief entry to this README
