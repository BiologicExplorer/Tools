# PROTEIN_BIOLOGY_CONTEXT.md

**BiologicExplorer / SEA Orchestrator Pipeline**
**Living document — update when new proteins or literature are added**

---

## Purpose

This document explains the `protein_knowledge_base.json` schema, describes how to add new entries, and provides the biological and literature context for the proteins currently encoded. It is intended to be consumed by the orchestrator pipeline at runtime to inform scoring decisions, and by analysts adding or reviewing entries.

The knowledge base is **protein-agnostic**: every field in the schema applies equally to any host or viral protein. CYP2E1 (P05181) and the HCV polyprotein (Q9WMX2) are the first entries and serve as worked examples of a general pattern, not as special cases.

---

## Schema Overview

`protein_knowledge_base.json` is a JSON object with four top-level keys:

| Key | Description |
|-----|-------------|
| `_schema_version` | Semantic version string for the schema itself |
| `_description` | Short human-readable description of the file |
| `_template` | A blank entry showing every field and its type (copy this to add a new protein) |
| `proteins` | Dict of `{ "UNIPROT_ACCESSION": { ...entry... } }` |

---

## Per-Protein Entry Fields

### Identification

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `uniprot_accession` | string | yes | UniProt primary accession (e.g. `"P05181"`) |
| `gene_name` | string | recommended | HGNC gene symbol or viral gene name |
| `protein_name` | string | recommended | Full protein name from UniProt |
| `organism` | string | recommended | Latin binomial (e.g. `"Homo sapiens"`) |
| `taxon_id` | int or null | optional | NCBI taxonomy ID |

### CMA Substrate Status

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `cma_substrate` | bool or null | recommended | `true` if evidence supports CMA substrate status; `null` if unknown |
| `cma_evidence` | array of citation objects | required if true | List of supporting citations (see Citation Object below). Internal pipeline evidence (KFERQ motif count) is valid — use `"journal": "internal"` and describe the finding in `"note"` |

CMA substrate status is determined by one or more of:
1. Presence of ≥1 KFERQ-like pentapeptide motif (detected by `motif_finder` module)
2. Published biochemical evidence of LAMP-2A–mediated lysosomal delivery
3. Published proteomic CMA substrate catalogues

### Known Autoantibody Epitopes

`known_autoantibody_epitopes` is an array. Each element represents one experimentally validated or computationally predicted immunodominant region.

| Field | Type | Description |
|-------|------|-------------|
| `label` | string | Short identifier for the epitope (e.g. `"JHDN-5"`) |
| `residue_start` / `residue_end` | int or null | 1-based inclusive residue range |
| `sequence` | string | Sequence or shorthand (e.g. `"Gly113...Leu135"`) |
| `mhc_restriction` | string | One of: `"MHC_I"`, `"MHC_II"`, `"conformational"`, `"unknown"` |
| `antibody_types` | array of strings | Detected Ig isotypes (e.g. `["IgG4"]`) |
| `critical_residues` | array | Residues with special functional roles (see below) |
| `disease_associations` | array of strings | Diseases where this epitope is documented |
| `notes` | string | Free-text summary of functional significance |
| `evidence` | array of citation objects | Supporting literature |

**Critical residue object:**
```json
{ "position": 123, "amino_acid": "Lys", "note": "modification required for immunogenicity" }
```

### Known Mimicry Partners

`known_mimicry_partners` is an array. Each element describes one molecular mimicry relationship between this protein and a partner protein.

| Field | Type | Description |
|-------|------|-------------|
| `partner_accession` | string | UniProt accession of the mimicry partner |
| `mimicry_type` | string | e.g. `"molecular_mimicry_MHC_II"`, `"CTL_mimicry"`, `"conformational_mimicry"` |
| `host_region_start` / `host_region_end` | int or null | Residue range on this protein involved in mimicry (null if not mapped) |
| `partner_region_start` / `partner_region_end` | int or null | Residue range on the partner involved in mimicry |
| `notes` | string | Mechanism, evidence quality, caveats |
| `evidence` | array of citation objects | Supporting literature |

### Citation Object

Used uniformly inside `cma_evidence`, `known_autoantibody_epitopes[*].evidence`, and `known_mimicry_partners[*].evidence`:

```json
{
  "citation_key": "AuthorYear",
  "authors": "Last F, ...",
  "year": 2024,
  "title": "Full title",
  "journal": "Journal Name",
  "doi": "10.xxxx/...",
  "pmid": "12345678",
  "note": "(optional) clarifying note about scope or relevance"
}
```

For internal pipeline evidence use `"journal": "internal"` and leave `doi` and `pmid` empty.

### Other Fields

| Field | Type | Description |
|-------|------|-------------|
| `disease_associations` | array of strings | All diseases the protein is associated with |
| `notes` | string | General biology notes not covered by the structured fields |

---

## How to Add a New Entry

1. Copy the `_template` object from `protein_knowledge_base.json`.
2. Set `uniprot_accession` to the primary UniProt ID of the new protein.
3. Fill in identification fields from the UniProt record.
4. For `cma_substrate`: run the pipeline (`motif_finder` module) on the protein sequence. If ≥1 KFERQ-like motif is found, set to `true` and record the motif count in `cma_evidence` with `"journal": "internal"`. If published evidence exists, add that citation as well. Leave `null` if unknown.
5. For `known_autoantibody_epitopes`: add only experimentally validated or SYFPEITHI/NetMHC-predicted epitopes supported by ≥1 peer-reviewed citation.
6. For `known_mimicry_partners`: add entries where published literature or pipeline McLachlan scoring (composite_primary ≥ 15.0 at ≥1 SUPER_EPITOPE hit in the known epitope region) supports the relationship.
7. Add the protein to both directions: if A mimics B, add `A.known_mimicry_partners → B` AND `B.known_mimicry_partners → A`.
8. Commit both `protein_knowledge_base.json` and (if substantive biology was added) an update to this file.

---

## CYP2E1 (P05181) — Biology and Literature Summary

### Overview

Cytochrome P450 2E1 (CYP2E1, UniProt P05181) is a hepatic microsomal enzyme expressed primarily in hepatocytes. Its physiological role is oxidative biotransformation of low-molecular-weight xenobiotics — ethanol, halogenated anesthetics (halothane, isoflurane), acetaminophen — and certain endogenous substrates. CYP2E1 expression is induced by ethanol consumption, obesity, and diabetes. It generates reactive oxygen species (ROS) as a by-product, contributing to oxidative stress in the diseased liver.

CYP2E1 is transported from the Golgi apparatus via secretory vesicles to the hepatocyte plasma membrane, where surface-expressed enzyme can be recognized by circulating autoantibodies, enabling complement-dependent antibody-mediated cytotoxicity (Sutti et al. 2014).

### Autoimmune Significance

Anti-CYP2E1 autoantibodies are detected in three clinical settings:

1. **Halogenated anesthetic-induced hepatitis** — CYP2E1 metabolizes the anesthetic to a reactive intermediate that forms protein adducts, modifying Lys123 and breaking self-tolerance.
2. **Alcoholic liver disease (ALD)** — Anti-CYP2E1 IgG detected in ~30% of patients with advanced ALD. Presence is an independent predictor of necro-inflammation and fibrosis.
3. **Chronic hepatitis C (CHC)** — Anti-CYP2E1 IgG detectable in CHC patients; anti-CYP2E1 IgG4 targeting the JHDN-5 epitope independently associated with severe hepatic fibrosis (P=0.0142; McCarthy et al. 2018).

In ALD, the breaking of tolerance is driven by **reactive metabolites** (hydroxyethyl radicals, HER) that modify CYP2E1 structurally. In CHC, **molecular mimicry** with HCV polyprotein sequences is the proposed primary mechanism (Sutti et al. 2014).

### Immunodominant Epitope: JHDN-5 (Gly113–Leu135)

McCarthy et al. (2018, mSphere, DOI: 10.1128/msphere.00453-18) identified **Gly113–Leu135** as a novel, MHC II-restricted immunodominant epitope shared between anesthetic-induced and viral hepatitis. Key findings:

- SYFPEITHI epitope prediction identified Gly113–Leu135 as a candidate MHC II-presented peptide.
- Immunization with Lys123-modified JHDN-5 (but not unmodified) induced hepatitis and CYP2E1 autoantibodies in mice.
- JHDN-5 antiserum recognized mitochondria and ER; upregulated HSP27; induced mitochondrial oxidative stress via complex I inhibition (P<0.001); and inhibited CYP2E1 enzymatic activity.
- JHDN-5 IgG4 was elevated in 200 HCV-positive patients from the ALIVE cohort versus controls (P<0.001).
- After covariate adjustment, elevated JHDN-5 IgG4 was independently associated with severe hepatic fibrosis (P=0.0142).

**Lys123 is the critical residue**: modification (adduct formation) is required to generate immunogenicity. This is directly relevant to pipeline scoring — the McLachlan hit at viral positions mirroring host Lys123 is the highest-priority finding from the Q9WMX2 vs P05181 run.

### Conformational Epitopes in ALD and CHC

In addition to the MHC II linear epitope JHDN-5, Sutti et al. (2014) established that anti-CYP2E1 autoantibodies in ALD and CHC predominantly target **conformational epitopes** that are co-located on the molecular surface of CYP2E1. These are not linear peptide sequences and are not readily decomposable into start/end residues from current literature. Their surface-exposed location is consistent with the membrane-targeting mechanism.

### CMA Substrate Evidence

The BiologicExplorer pipeline (`motif_finder` module) detected **5 KFERQ-like pentapeptide motifs** in the CYP2E1 (P05181) sequence. The presence of KFERQ motifs is the primary biochemical basis for CMA substrate classification; HSC70 recognition of these motifs initiates lysosomal delivery via LAMP-2A. No published biochemical study directly demonstrating CYP2E1 CMA degradation was identified in the current literature review; the CMA substrate flag is therefore based on motif analysis (internal evidence) and should be validated experimentally if used as a key decision variable.

---

## HCV Polyprotein (Q9WMX2) — Biology and Literature Summary

### Overview

The Hepatitis C virus (HCV) polyprotein (UniProt Q9WMX2) is a ~3000-residue precursor cleaved by viral and host proteases into ten mature proteins: Core, E1, E2, p7, NS2, NS3, NS4A, NS4B, NS5A, NS5B. The Core protein occupies approximately the first 191 residues of the polyprotein. The N-terminal core region (approx. residues 40–150) is hydrophobic and RNA-binding, and is implicated in molecular mimicry of hepatic CYPs.

### CYP2E1 Mimicry

Sutti et al. (2014) identified molecular mimicry as the mechanism by which HCV infection breaks self-tolerance to CYP2E1, distinct from the reactive-metabolite mechanism in ALD. This is supported by the pipeline:

The McLachlan SEA run (Q9WMX2 vs P05181, composite_primary ≥ 15.0 cutoff) yielded **18 SUPER_EPITOPE tier hits** with overall risk **CRITICAL** and top score 20.2898. Three scored hits land directly in the CYP2E1 JHDN-5 epitope range (Gly113–Leu135):

| Viral position | Host position | Viral window | Host window | Score |
|---------------|---------------|-------------|------------|-------|
| 2525–2531 | 123–129 | KDVRNLS | KDIRRFS | 16.66 |
| 953–959 | 130–136 | LTPLRDW | LTTLRNY | 16.13 |
| 676–682 | 127–133 | PCSFTTL | RFSLTTL | 15.05 |

Super-Epitope #10 in the SEA report places viral pos 110 against host pos 123 — precisely Lys123, the residue whose modification triggers immunogenicity. The pipeline found this independently, without prior knowledge of the epitope.

### Related but Distinct: CYP2A6 Mimicry

Kammer et al. (1999, J Exp Med, DOI: 10.1084/jem.190.2.169) reported that HCV core peptide 178–187 (LLALLSCLTV) mimics **CYP2A6** residues 8–17 at the level of CD8+ CTL recognition (MHC I). This is a separate finding from the CYP2E1 / MHC II mimicry described above. Both CYP2A6 and CYP2E1 are members of the CYP2 subfamily, and HCV may use cross-reactive mimicry across multiple CYP2 members. The CYP2A6 finding is cited in the Q9WMX2 entry as related context, with an explicit note that it is not CYP2E1.

---

## Design Rationale

### Why a Generalized Schema?

The pipeline is designed to compare any viral protein against any host protein. Hardcoding logic or data structures for specific proteins (CYP2E1, Q9WMX2) would break generalizability and require code changes every time a new protein pair is studied. The knowledge base externalizes protein-specific facts as data, which the orchestrator reads at runtime. The orchestrator itself remains a thin integration layer.

### What the Orchestrator Can Do With This Data (Future Integration)

The following scoring enhancements are identified for a future session (do **not** implement until explicitly scoped):

1. **Epitope proximity bonus** — if a McLachlan hit's host anchor position falls within a `known_autoantibody_epitopes` range for the host protein, add a configurable bonus to the composite score.
2. **CMA override** — if both viral and host proteins have `cma_substrate: true`, the degradation pathway analysis gets a higher confidence weight.
3. **Known Epitope Coverage metric** — report what fraction of known epitope residues are covered by SUPER_EPITOPE-tier hits, as a separate output field.
4. **Mimicry partner flag** — if the queried viral protein is already in `known_mimicry_partners` for the host protein, flag this in the report header.

### Version History

| Version | Date | Change |
|---------|------|--------|
| 1.0.0 | 2026-09-30 | Initial schema. First entries: P05181 (CYP2E1) and Q9WMX2 (HCV polyprotein). Literature: McCarthy 2018, Sutti 2014, Kammer 1999. |

---

## References

1. **McCarthy LR et al. (2018).** Identification of a Novel MHC-Restricted CYP2E1 Epitope Associated with Anesthetic and Viral Hepatitis. *mSphere*. DOI: [10.1128/msphere.00453-18](https://doi.org/10.1128/msphere.00453-18)

2. **Sutti S, Vidali M, Albano E. (2014).** CYP2E1 autoantibodies in liver diseases. *Redox Biology* 3:72–78. DOI: [10.1016/j.redox.2014.11.004](https://doi.org/10.1016/j.redox.2014.11.004) PMID: 25460731

3. **Kammer AR et al. (1999).** Molecular mimicry of human cytochrome P450 by hepatitis C virus at the level of cytotoxic T cell recognition. *J Exp Med* 190(2):169–176. DOI: [10.1084/jem.190.2.169](https://doi.org/10.1084/jem.190.2.169) PMID: 10432280
