# Meaningful Findings — CYP2E1 Autoimmunogenicity SEA Analysis
**Project:** BiologicExplorer / Tools — SEA Orchestrator Pipeline  
**Host protein:** CYP2E1 · UniProt P05181 · 493 aa  
**Viral proteome:** HCV polyprotein · UniProt Q9WMX2 · 3010 aa  
**Method:** McLachlan 1972 matrix · Two-Pass Homology Search (11-mer broad seed → 5-mer core)  
**Last updated:** 2026-10-01

> This document accumulates only findings with direct biological or immunological significance.
> Implementation decisions and pipeline parameters are tracked separately in PHASE1_METHOD_PHASE4_FOUNDATION.md.

---

## F-01 — CYP2E1 Autoantibody Epitope Confirmed (Residues 99–132)

**Finding:** The canonical autoantibody-reactive epitope on CYP2E1 is the 34-aa linear sequence  
`GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT` (P05181 residues 99–132).  
**Verification:** Live-fetched from UniProt P05181 canonical FASTA — no cached string.  
**Significance:** Establishes the ground-truth immunological target against which all molecular mimicry candidates are evaluated. All homology pairs that overlap this window are elevated to priority status.

---

## F-02 — 17 Statistically Distinct McLachlan Homology Pairs Identified

**Finding:** Two-pass McLachlan search across the full HCV polyprotein × CYP2E1 sequence space yielded **17 confirmed homology pairs** (from 184 broad-neighborhood candidates).  
**Pass 1:** 11-mer broad seed, threshold ≥ 22.0 → 184 hits  
**Pass 2:** Best 5-mer sub-window per hit, sum ≥ 17.5 · avg ≥ 3.5/position → 17 final pairs  
**Pair distribution by viral protein:**

| Viral Protein | Pairs | Score Range |
|---------------|-------|-------------|
| NS3 Helicase  | 5 (P01–P05) | 24–27 |
| NS5A Domain I | 4 (P06–P09) | 21–24 |
| NS5B Polymerase | 8 (P10–P17) | 20–26 |

**Significance:** NS5B accounts for nearly half (8/17) of all homology pairs, making it the most sequence-similar HCV protein to CYP2E1 in this analysis.

---

## F-03 — All Three Epitope-Overlapping Pairs Originate Exclusively from NS5B

**Finding:** Of the 17 pairs, exactly 3 overlap the confirmed autoantibody epitope (residues 99–132). All three are NS5B polymerase hits:

| Pair | NS5B locus (UniProt) | CYP2E1 window | Core | Score |
|------|---------------------|---------------|------|-------|
| P15 ★ | 2521–2531 (local 102–112) | 120–130 | AKDVR / WKDIR | 23 |
| P16 ★ | 2541–2551 (local 122–132) | 122–132 | VWKDL / TWKDI | 23 |
| P17 ★ | 2643–2653 (local 224–234) | 117–127 | TENDI / TWKDI | 20 |

**Significance:** The complete absence of NS3 or NS5A epitope-overlapping pairs, combined with NS5B's disproportionate total pair count, positions **NS5B as the primary molecular mimicry candidate** for HCV-driven CYP2E1 autoimmunogenicity. This is consistent with the established role of NS5B in chronic HCV replication and sustained antigen presentation.

---

## F-04 — Three Distinct NS5B Loci Converge on a Single CYP2E1 Epitope Window

**Finding:** P15, P16, and P17 originate from three non-overlapping loci within NS5B (local positions 102–112, 122–132, and 224–234 respectively), yet all map host-side onto the narrow 34-aa CYP2E1 epitope window (residues 99–132).  
**Significance:** Multiple independent NS5B sequence fragments target the same host epitope region. This convergence suggests the epitope window may present a physicochemical surface that is specifically mimicked by the NS5B fold, rather than being an artifact of a single local homology event.

---

## F-05 — Super-Epitope / Assembled Epitope Hypothesis (P15 · P16 · P17)

**Finding:** In the CYP2E1 3D structure (PDB 3KOH), the host-side residues of P15 (host 120–130), P16 (host 122–132), and P17 (host 117–127) **cluster spatially near the autoantibody epitope**, as visually confirmed in the interactive 3D panel.  
**Implication:** Even though P15/P16/P17 derive from different regions of NS5B (separated by ~110–120 aa in the polyprotein), their host-side footprints converge in 3D space on CYP2E1. If the corresponding NS5B residues are brought into proximity by the NS5B fold (proximal hinge), the viral surface could **present an assembled mimic** of the CYP2E1 epitope.  
**Coined term:** "Assembled epitope" or "super-epitope."  
**Status:** Visually confirmed in 3D panel. Quantitative hinge-proximity screen (Cα–Cα distance between NS5B residues 106 and 196) is filed for Phase 4 structural analysis.

---

## F-06 — KDVRNLS + PCSFTTL Candidate Super-Epitope Fragment Pair

**Finding:** Two short sequences on the HCV polyprotein — `KDVRNLS` and `PCSFTTL` — match CYP2E1 residues `KDIRRFS` and `RFSLTTL` respectively.  
**Hypothesis:** If a proximal hinge in the viral protein brings KDVRNLS and PCSFTTL into 3D proximity, the assembled surface could present as `KDIRRFSLTTL` — a contiguous 11-mer spanning the core of the confirmed autoantibody epitope (KDIRRFSLTT is residues 122–131 of the epitope).  
**Status:** Filed as Phase 4 super-epitope hypothesis. No structural data analyzed yet.

---

## F-07 — Epitope-Overlapping Pairs Score Below the Top-Tier Linear Matches

**Finding:** The three epitope-overlapping pairs (P15=23, P16=23, P17=20) rank 12th, 13th, and 17th respectively out of 17 by core score. The top linear scorers (P01=27 NS3, P02=26 NS3, P10=26 NS5B, P11=26 NS5B) have no epitope overlap.  
**Implication:** Linear sequence homology score alone does not predict biological priority. The epitope-overlapping pairs are immunologically most relevant despite having intermediate McLachlan scores. This motivates the Phase 3 SEA proximity scoring framework (jammer density, hinge clustering) and Phase 4 combined score as necessary layers beyond raw McLachlan homology.

---

## F-08 — NS3 Produces the Strongest Raw Linear Homology Scores

**Finding:** NS3 helicase yields the top-scoring pair (P01, score=27) and three of the top-seven pairs by core score. NS3 pairs P01–P04 all carry strong McLachlan cores of ≥25.  
**Significance:** Despite not overlapping the autoimmune epitope, the high NS3 homology could contribute to autoimmunogenicity via epitope spreading or cross-reactive T-cell priming. These pairs warrant independent investigation in Phase 4 structural analysis.

---

## F-09 — Two-Pass Search Recovers Biologically Relevant Low-Scoring Hits

**Finding:** P17 (score=20, the lowest of all 17 pairs) is one of the three epitope-overlapping pairs. A single-pass search with a higher threshold would have excluded it.  
**Significance:** Validates the two-pass design: using a broad 11-mer neighborhood in Pass 1 to anchor context, then scoring the best internal 5-mer in Pass 2, preserves biologically relevant hits that a strict linear threshold would discard. The assembled epitope signature (F-05) depends on recovering P17.

---

*Document is updated as new findings are confirmed. Each entry requires either: (a) direct computational verification against live UniProt/PDB data, or (b) explicit user confirmation from visual/experimental inspection.*
