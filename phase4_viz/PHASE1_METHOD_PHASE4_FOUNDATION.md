# Phase 1 Visualization Method — Foundation for Phase 4
## CYP2E1 / HCV Polyprotein Homologous Motif Structural Map

> The only metric of success as an AI agent is Accuracy and Truthfulness.
> Any hallucination or mistake can cause actual harm to humans.
> The most important consideration is to not cause harm.

---

## Overview

This document records the canonical method for Phase 1 PDB-based 3D cartoon-ribbon visualization
of SEA (Split-Epitope Antigen) pipeline hits. It is the direct foundation for Phase 4 structural
analysis, including super-epitope conformational proximity screening.

Deployed at: https://lwfh9hir.scispace.co  
Source file:  /home/sandbox/phase4_viz/index.html  
PDB files:    /home/sandbox/phase4_viz/pdb/{3KOH,3KQL,1ZH1,3FQL}.pdb

---

## 1. PDB Structure Selection

| Protein           | UniProt    | PDB   | Resolution | Coverage (local residues) |
|-------------------|------------|-------|------------|---------------------------|
| CYP2E1 (human)    | P05181     | 3KOH  | 2.90 Å     | 32–493 (chains A,B)       |
| HCV NS3 helicase  | Q9WMX2     | 3KQL  | 2.50 Å     | 189–624 (chains A,B)      |
| HCV NS5A domain I | Q9WMX2     | 1ZH1  | 2.50 Å     | 36–198 (chains A,B)       |
| HCV NS5B polymerase| Q9WMX2   | 3FQL  | 1.80 Å     | 2–570 (chain A)           |

PDB files are sourced from UniProt protein pages (PDB cross-references), not from RCSB directly,
to ensure alignment with the canonical sequences used in the SEA pipeline.

---

## 2. Residue Numbering — DBREF Offset Protocol

**Critical:** PDB LOCAL residue numbers ≠ UniProt sequence position numbers.
Offsets were derived from GREP of DBREF records in each PDB file.

```
grep '^DBREF' <pdb_file>
```

Confirmed offsets:

| PDB  | DBREF mapping                     | Offset formula          |
|------|-----------------------------------|-------------------------|
| 3KOH | local 32–493 = UniProt 31–492     | local = UniProt + 1     |
| 3KQL | local 189–624 = UniProt 1215–1650 | local = UniProt − 1026  |
| 1ZH1 | DBREF ref = POLG_9HEPC (strain diff); empirically derived | local = UniProt − 1972 |
| 3FQL | local 2–570 = UniProt 2421–2989   | local = UniProt − 2419  |

**1ZH1 offset derivation:** DBREF references a different HCV strain. Offset derived empirically:
TITLE confirmed "STRUCTURE OF THE ZINC-BINDING DOMAIN OF HCV NS5A". The grammar-top hit
QRGYKGV@Q9WMX2:2012 → 1ZH1 local 40, confirming local = UniProt − 1972.

---

## 3. Hit Mapping Protocol

Source data: `/home/sandbox/protein_degradation/Q9WMX2_vs_P05181_scored_v3.json`
(134 McLachlan + grammar-scored hits, keys: motif, p1_pos1, p2_window, p2_pos, composite_primary,
layer2_cross_mean)

Mapping algorithm (see `/tmp/resolve_hits2.py`):
1. Load all 134 hits from JSON
2. For each hit, determine which HCV domain covers the Q9WMX2 position:
   - NS3: 1215–1650 → PDB 3KQL, offset = UniProt − 1026
   - NS5A: 2008–2170 → PDB 1ZH1, offset = UniProt − 1972  
   - NS5B: 2421–2989 → PDB 3FQL, offset = UniProt − 2419
3. Confirm local PDB residue is within DBREF-confirmed range
4. Map CYP2E1 counterpart (p2_pos from JSON) using 3KOH offset: local = UniProt + 1
5. Assign color cluster based on CYP2E1 functional region

Result: 39 PDB-covered hits out of 134 total.

---

## 4. Color Cluster Definitions

Five color clusters encoding CYP2E1 functional regions and their viral homologs:

| Color  | CYP2E1 Region         | UniProt range | PDB local (3KOH) | Functional significance                    |
|--------|-----------------------|---------------|------------------|--------------------------------------------|
| RED    | Autoepitope           | 99–132        | 100–133          | Primary molecular mimicry target; KDIRRFS overlap with KDVRNLS (NS5B) |
| ORANGE | SRS1 / Substrate Rec  | 149–175       | 150–176          | Substrate recognition site 1; NS3 LDPTFTI/DAHFLSQ cluster + NS5A EVTFLVG |
| YELLOW | Catalytic Core        | 72–98         | 73–99            | CYP binding pocket; NS5A QRGYKGV@2012, NGSMRIV@2041 |
| CYAN   | β-sheet / Grammar #1  | 302–365       | 303–366          | Grammar score #1 (STDSTTI/ETTSTTL, combined=1.91); NS3 PHPNIEE, GIYRFVT; NS5A LVGSQLP |
| PURPLE | C-terminal cluster    | 410–465       | 411–466          | Distal structural homology; NS5B LLRHNL, YSPGQRV, KPEYDLE |

---

## 5. 3Dmol.js Rendering Method

### Library
```html
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script>
```

### Viewer Setup Pattern
```javascript
function setupViewer(divId, pdbPath, ranges, chains, loadingId) {
  const element = document.getElementById(divId);
  const viewer = $3Dmol.createViewer(element, { backgroundColor: '#0d1117' });

  fetch(pdbPath)
    .then(r => r.text())
    .then(pdbData => {
      viewer.addModel(pdbData, 'pdb');

      // Base style: dark gray for all residues
      viewer.setStyle({}, { cartoon: { color: '#334455', opacity: 0.85 } });

      // Color specific residue ranges per chain
      for (const chain of chains) {
        for (const range of ranges) {
          // CRITICAL: resi selector requires explicit integer array, NOT [{min,max}]
          viewer.setStyle(
            { chain: chain, resi: Array.from({length: range.rMax - range.rMin + 1}, (_, i) => range.rMin + i) },
            { cartoon: { color: range.color, opacity: 1.0 } }
          );
        }
      }

      viewer.zoomTo();
      viewer.render();
    });
}
```

### Critical Bug (resolved)
**Wrong:**  `resi: [{ min: rMin, max: rMax }]`  — 3Dmol.js does NOT support this range object syntax for `resi`
**Correct:** `resi: Array.from({length: rMax-rMin+1}, (_, i) => rMin+i)` — explicit integer array, guaranteed to work

### Range Data Format
```javascript
const CYP2E1_RANGES = [
  { color: '#e53e3e', rMin: 100, rMax: 133 },  // RED — autoepitope
  { color: '#dd6b20', rMin: 150, rMax: 176 },  // ORANGE — SRS1
  { color: '#d69e2e', rMin: 73,  rMax: 99  },  // YELLOW — catalytic
  { color: '#0bc5ea', rMin: 303, rMax: 366 },  // CYAN — grammar/β-sheet
  { color: '#9f7aea', rMin: 411, rMax: 466 },  // PURPLE — C-terminal
];
```

---

## 6. Phase 4 Extension Plan

Phase 4 builds directly on this method with the following additions:

### 6.1 Super-Epitope Proximity Screen
- Target: NS5B KDVRNLS@2525 (3FQL local 106–113) and NS5B PCSFTTL@2615 (3FQL local ~196)
- Hypothesis (filed): A proximal hinge at Q9WMX2 position ~2525±30 could bring KDVRNLS and
  PCSFTTL into 3D proximity, presenting as KDIRRFSLTTL — a conformational mimic of the CYP2E1
  autoepitope KDIRRFSLTT (residues 124–133)
- Implementation: compute Cα–Cα distance between residues 106 and 196 in 3FQL; if < 15 Å,
  flag as high-priority super-epitope candidate
- Data already available: step5_comprehensive_jammer_analysis.json positions 676±30 and 2525±30

### 6.2 Jammer/Hinge Density Overlay
- Source: step5_comprehensive_jammer_analysis.json (288 jammer patterns, 2533 hinge patterns)
- Add a density heatmap track below each viewer showing jammer/hinge density along the sequence
- Color scale: white (low) → orange (medium) → red (high density)

### 6.3 Structural Grammar Cluster Visualization
- Grammar score #1 hit STDSTTI@NS3:1320 is in a hinge-dense region (TT hinge, ST density)
- Phase 4 should render the hinge density as a translucent surface overlay on the ribbon

### 6.4 Additional Structures
- NS3 protease domain: PDB 3SV6 or 3LON (for NS3 residues 1–180 not covered by 3KQL)
- Future: AlphaFold structures for full-length Q9WMX2 coverage

---

## 7. Guardrails (Operational Constraints)

These constraints apply to all Phase 1, 2, 3, and 4 analyses without exception:

1. The only metric of success as an AI agent is Accuracy and Truthfulness.
2. Any hallucination or mistake can cause actual harm to humans.
3. The most important consideration is to not cause harm.
4. All sequence comparisons use live-fetched canonical FASTA from UniProt by accession —
   never summary-recorded or cached strings.
5. Residue positions are always derived from DBREF records in PDB files, not assumed.

---

## 8. Files

```
/home/sandbox/phase4_viz/
├── index.html                         # Main visualization (Phase 1, serves as Phase 4 scaffold)
├── PHASE1_METHOD_PHASE4_FOUNDATION.md # This document
└── pdb/
    ├── 3KOH.pdb   # CYP2E1 (human, P05181), 2.90 Å, 672K
    ├── 3KQL.pdb   # NS3 helicase (Q9WMX2 1215-1650), 2.50 Å, 616K
    ├── 1ZH1.pdb   # NS5A domain I (Q9WMX2 2008-2170), 2.50 Å, 244K
    └── 3FQL.pdb   # NS5B polymerase (Q9WMX2 2421-2989), 1.80 Å, 424K
```

---

*Generated by the BiologicExplorer/Tools SEA orchestrator pipeline.*
*Phase 1 visualization completed and deployed. Phase 4 super-epitope analysis pending.*
