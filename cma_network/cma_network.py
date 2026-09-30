#!/usr/bin/env python3
"""
cma_network.py  v1.0
--------------------
Chaperone-Mediated Autophagy (CMA) Network Gene Lookup and Scoring Tool.
Self-contained: no runtime web or API calls.

Sources
-------
  Kirchner P et al. (2019) "Proteome-wide analysis of chaperone-mediated
  autophagy targeting motifs." PLoS Biology 17(6):e3000301.
  DOI: 10.1371/journal.pbio.3000301
  -> Defines CMA network gene list, categories, and +/- directions.

  Bourdenx M et al. (2021) "Chaperone-mediated autophagy prevents collapse
  of the neuronal metastable proteome." Cell 184(10):2696-2714.
  DOI: 10.1016/j.cell.2021.03.048   PMC: PMC8152331
  -> Defines CMA activation score formula:
       score = sum(weight_i * direction_i * expr_i) / sum(weight_i)
     LAMP-2A weight=2 (rate-limiting step); all other genes weight=1.

Functions
---------
  check_cma_network_membership(query)        -> dict or None
  list_cma_network()                         -> List[dict]
  calculate_cma_score(expression_dict, ...)  -> dict

Import usage:
  from cma_network import (check_cma_network_membership,
                            calculate_cma_score,
                            list_cma_network,
                            CMA_NETWORK)

CLI usage:
  python3 cma_network.py --check LAMP2
  python3 cma_network.py --check P13473
  python3 cma_network.py --list-network
  python3 cma_network.py --score '{"LAMP2": 5.2, "HSPA8": 3.1, "AKT1": 1.8}'
  python3 cma_network.py --self-test
"""

import sys
import json
import argparse
from typing import Dict, List, Optional


# ===============================================================================
#  CMA NETWORK DATABASE
#
#  Source: Kirchner et al. 2019, PLoS Biol 17(6):e3000301 (Fig 4 + legend).
#  Direction (+1/-1) cross-validated against:
#    Bourdenx et al. 2021, Cell 184:2696 (Figure S10L immunoblot panel).
#
#  UniProt IDs: human canonical isoforms.
#  High-confidence IDs are well-established reference proteins.
#  IDs marked [*] are moderately confident -- verify at https://www.uniprot.org
# ===============================================================================

CMA_NETWORK: Dict[str, Dict] = {

    # --------------------------------------------------------------------------
    # EFFECTORS  (Core Machinery)
    # Essential structural/functional components; loss abolishes CMA.
    # LAMP-2A is the rate-limiting receptor -> weight=2 per Bourdenx 2021.
    # --------------------------------------------------------------------------

    "LAMP2": {
        "protein_name":  "Lysosome-associated membrane protein 2 (LAMP-2A isoform)",
        "aliases":       ["LAMP2", "LAMP2A", "LAMP-2A", "CD107B"],
        "uniprot_human": "P13473",
        "category":      "effector",
        "direction":     +1,
        "weight":        2,
        "notes": (
            "Rate-limiting receptor for CMA; only the LAMP-2A splice isoform "
            "participates. Substrate-Hsc70 complexes bind its cytosolic tail to "
            "initiate translocation. Weight=2 per Bourdenx 2021 (Cell 184:2696)."
        ),
    },

    "HSPA8": {
        "protein_name":  "Heat shock cognate 71 kDa protein (HSC70)",
        "aliases":       ["HSPA8", "HSC70", "HSP73", "HSPA10"],
        "uniprot_human": "P11142",
        "category":      "effector",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "Cytosolic chaperone that recognises KFERQ-like motifs on CMA substrates "
            "and escorts them to the lysosomal membrane."
        ),
    },

    "HSP90AA1": {
        "protein_name":  "Heat shock protein HSP 90-alpha (HSP90)",
        "aliases":       ["HSP90AA1", "HSP90", "HSP90A", "HSP86", "HSPCA"],
        "uniprot_human": "P07900",
        "category":      "effector",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "Stabilises LAMP-2A at the lysosomal membrane; required for assembly "
            "of the translocation complex."
        ),
    },

    "DNAJB1": {
        "protein_name":  "DnaJ homolog subfamily B member 1 (HSP40)",
        "aliases":       ["DNAJB1", "HSP40", "HSP40B1", "HSPF1"],
        "uniprot_human": "P25685",
        "category":      "effector",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "J-domain co-chaperone that stabilises substrate-Hsc70 complexes "
            "and stimulates Hsc70 ATPase activity."
        ),
    },

    # --------------------------------------------------------------------------
    # LYSOSOMAL MODULATORS
    # Act at or within the lysosomal membrane to tune CMA activity.
    # --------------------------------------------------------------------------

    "MT-RNR2": {
        "protein_name":  "Humanin (mitochondria-derived peptide)",
        "aliases":       ["MT-RNR2", "Humanin", "MTRNR2", "HN"],
        "uniprot_human": None,   # mitochondria-encoded peptide; no standard UniProt
        "category":      "lysosomal_modulator",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "Mitochondria-derived peptide encoded within the 16S rRNA locus. "
            "Activates CMA at the lysosomal membrane. "
            "Partially conserved across species (Kirchner 2019)."
        ),
    },

    "GFAP": {
        "protein_name":  "Glial fibrillary acidic protein",
        "aliases":       ["GFAP"],
        "uniprot_human": "P14136",
        "category":      "lysosomal_modulator",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "Organises lipid microdomains at the lysosomal membrane; "
            "facilitates LAMP-2A multimerisation into the translocation complex. "
            "Partially conserved across species (Kirchner 2019)."
        ),
    },

    "RAC1": {
        "protein_name":  "Ras-related C3 botulinum toxin substrate 1",
        "aliases":       ["RAC1", "RAC-1", "TC25", "p21-Rac1"],
        "uniprot_human": "P63000",
        "category":      "lysosomal_modulator",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "GTPase that organises the lysosomal membrane microtubule network "
            "and promotes LAMP-2A multimerisation."
        ),
    },

    "PHLPP1": {
        "protein_name":  "PH domain and leucine-rich repeat protein phosphatase 1",
        "aliases":       ["PHLPP1", "PHLPP", "SCOP"],
        "uniprot_human": "O60346",   # [*] moderately confident
        "category":      "lysosomal_modulator",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "Phosphatase that dephosphorylates AKT1; "
            "indirectly activates CMA by relieving AKT1-mediated inhibition of LAMP-2A."
        ),
    },

    "EEF1A1": {
        "protein_name":  "Elongation factor 1-alpha 1 (eF1-alpha)",
        "aliases":       ["EEF1A1", "EF1A", "EEF1A", "eEF1A1"],
        "uniprot_human": "P68104",
        "category":      "lysosomal_modulator",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "Translation elongation factor with a non-canonical role at the "
            "lysosomal membrane supporting LAMP-2A function."
        ),
    },

    "RICTOR": {
        "protein_name":  "Rapamycin-insensitive companion of mTOR (RICTOR)",
        "aliases":       ["RICTOR"],
        "uniprot_human": "Q6IA86",   # [*] moderately confident
        "category":      "lysosomal_modulator",
        "direction":     -1,
        "weight":        1,
        "notes": (
            "mTORC2 scaffold subunit; the mTORC2 complex phosphorylates and "
            "inhibits LAMP-2A at the lysosomal membrane."
        ),
    },

    "AKT1": {
        "protein_name":  "RAC-alpha serine/threonine-protein kinase (AKT1 / PKB)",
        "aliases":       ["AKT1", "AKT", "PKB", "PKBalpha", "RAC-PK-alpha"],
        "uniprot_human": "P31749",
        "category":      "lysosomal_modulator",
        "direction":     -1,
        "weight":        1,
        "notes": (
            "Serine/threonine kinase that phosphorylates LAMP-2A, "
            "reducing its stability and inhibiting CMA activity."
        ),
    },

    "CTSA": {
        "protein_name":  "Lysosomal protective protein / Cathepsin A (CathA)",
        "aliases":       ["CTSA", "CathA", "PPCA", "GLB2", "PPGB"],
        "uniprot_human": "P10619",
        "category":      "lysosomal_modulator",
        "direction":     -1,
        "weight":        1,
        "notes": (
            "Lysosomal serine protease that cleaves LAMP-2A at the lysosomal "
            "membrane, promoting its degradation and inhibiting CMA."
        ),
    },

    # --------------------------------------------------------------------------
    # EXTRA-LYSOSOMAL MODULATORS
    # Act upstream of the lysosome (nucleus, cytoplasm, endosomes).
    # --------------------------------------------------------------------------

    "NFATC1": {
        "protein_name":  "Nuclear factor of activated T-cells, cytoplasmic 1 (NFAT)",
        "aliases":       ["NFATC1", "NFAT", "NFAT2", "NF-ATc"],
        "uniprot_human": "O95644",   # [*] Kirchner 2019 refers to 'NFAT' generically
        "category":      "extralysosomal_modulator",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "Transcription factor that upregulates LAMP-2A gene expression. "
            "Partially conserved across species (Kirchner 2019)."
        ),
    },

    "NFE2L2": {
        "protein_name":  "Nuclear factor erythroid 2-related factor 2 (NRF2)",
        "aliases":       ["NFE2L2", "NRF2", "NRF-2", "HEBP1"],
        "uniprot_human": "Q16236",
        "category":      "extralysosomal_modulator",
        "direction":     +1,
        "weight":        1,
        "notes": (
            "Transcription factor upregulating LAMP-2A and HSP90 expression. "
            "Partially conserved across species (Kirchner 2019)."
        ),
    },

    "RAB11A": {
        "protein_name":  "Ras-related protein Rab-11A",
        "aliases":       ["RAB11A", "RAB11", "RAB-11"],
        "uniprot_human": "P62491",
        "category":      "extralysosomal_modulator",
        "direction":     -1,
        "weight":        1,
        "notes": (
            "GTPase that promotes recycling of LAMP-2A from the lysosomal lumen "
            "back into internal vesicles, reducing membrane LAMP-2A availability."
        ),
    },

    "RARA": {
        "protein_name":  "Retinoic acid receptor alpha (RAR-alpha)",
        "aliases":       ["RARA", "RAR-alpha", "RARa", "NR1B1"],
        "uniprot_human": "P10276",
        "category":      "extralysosomal_modulator",
        "direction":     -1,
        "weight":        1,
        "notes": (
            "Sequesters LAMP-2A within lipid rafts at the lysosomal membrane, "
            "preventing its translocation to functional CMA complexes."
        ),
    },
}


# ===============================================================================
#  LOOKUP INDEX  (built once at import time)
# ===============================================================================

def _build_lookup_index() -> Dict[str, str]:
    """Return a case-insensitive alias/UniProt -> canonical gene_symbol mapping."""
    index: Dict[str, str] = {}
    for gene, info in CMA_NETWORK.items():
        index[gene.upper()] = gene
        for alias in info["aliases"]:
            index[alias.upper()] = gene
        if info["uniprot_human"]:
            index[info["uniprot_human"].upper()] = gene
    return index


_LOOKUP_INDEX: Dict[str, str] = _build_lookup_index()


# ===============================================================================
#  PUBLIC API
# ===============================================================================

def check_cma_network_membership(query: str) -> Optional[Dict]:
    """
    Check whether a protein is a member of the CMA regulatory network.

    Parameters
    ----------
    query : str
        Gene symbol, common alias, or UniProt accession (case-insensitive).
        Examples: "LAMP2", "LAMP-2A", "HSC70", "P13473", "AKT1", "PKB"

    Returns
    -------
    dict  if the protein is in the CMA network:
        {
          'gene_symbol'   : str,   canonical HGNC symbol
          'protein_name'  : str,
          'aliases'       : list[str],
          'uniprot_human' : str or None,
          'category'      : str,   effector | lysosomal_modulator |
                                   extralysosomal_modulator
          'direction'     : int,   +1 (activates CMA) or -1 (inhibits CMA)
          'weight'        : int,   2 for LAMP2A; 1 for all others
          'notes'         : str,
          'in_network'    : True
        }
    None  if the protein is not in the CMA network.

    Examples
    --------
    >>> check_cma_network_membership("LAMP2")['direction']
    1
    >>> check_cma_network_membership("P13473")['gene_symbol']
    'LAMP2'
    >>> check_cma_network_membership("TP53") is None
    True
    """
    canonical = _LOOKUP_INDEX.get(query.strip().upper())
    if canonical is None:
        return None
    result = dict(CMA_NETWORK[canonical])
    result["gene_symbol"] = canonical
    result["in_network"] = True
    return result


def list_cma_network() -> List[Dict]:
    """
    Return all 16 CMA network genes as a list of membership dicts,
    sorted by category (effectors first) then gene symbol.

    Returns
    -------
    List[dict]  — each dict has the same structure as check_cma_network_membership().
    """
    category_order = {
        "effector": 0,
        "lysosomal_modulator": 1,
        "extralysosomal_modulator": 2,
    }
    results = []
    for gene, info in CMA_NETWORK.items():
        d = dict(info)
        d["gene_symbol"] = gene
        d["in_network"] = True
        results.append(d)
    results.sort(key=lambda x: (category_order.get(x["category"], 9), x["gene_symbol"]))
    return results


def calculate_cma_score(
    expression_dict: Dict[str, float],
    reference_dict: Optional[Dict[str, float]] = None,
) -> Dict:
    """
    Calculate the CMA activation score for a sample.

    The score is the weighted, direction-signed average of (log-normalised)
    expression values for CMA network genes, following the formula in
    Bourdenx et al. 2021 (Cell 184:2696, DOI: 10.1016/j.cell.2021.03.048).

    Formula
    -------
      score = sum(weight_i * direction_i * expr_i) / sum(weight_i)

    where the sum runs over all CMA network genes present in expression_dict.
    A positive score indicates net CMA activation; negative indicates inhibition.

    Parameters
    ----------
    expression_dict : dict {str: float}
        Log-normalised expression values keyed by gene name or alias
        (case-insensitive; aliases accepted).
        Recommended normalisation: log1p(counts / total_counts * 10_000),
        consistent with Bourdenx 2021 / SCANPY log-normalise defaults.
        Non-network genes are silently ignored.

    reference_dict : dict {str: float}, optional
        Reference (e.g. healthy-control) expression values for the same genes.
        When supplied, each gene's expression is first expressed as
        fold-change (expr / ref) before scoring, consistent with the
        Bourdenx 2021 instruction to express the score as "fold of healthy
        individuals within a given cell type."
        If a gene is missing from reference_dict, it is excluded from scoring
        and a warning is recorded.

    Returns
    -------
    dict:
        {
          'cma_score'       : float or None,   None if no genes matched
          'n_genes_scored'  : int,
          'n_genes_missing' : int,
          'scored_genes'    : list[dict],      per-gene breakdown
          'missing_genes'   : list[str],       canonical symbols not in input
          'warnings'        : list[str],
          'formula_ref'     : str,
        }

    Each entry in 'scored_genes':
        {
          'gene_symbol'  : str,
          'category'     : str,
          'direction'    : int,
          'weight'       : int,
          'expression'   : float,   raw or fold-change, as used in scoring
          'contribution' : float,   weight * direction * expression
        }
    """
    warnings_out: List[str] = []
    scored_genes: List[Dict] = []
    missing_genes: List[str] = []

    total_weighted = 0.0
    total_weight   = 0.0

    # Pre-resolve reference dict to canonical symbols
    ref_canonical: Dict[str, float] = {}
    if reference_dict:
        for k, v in reference_dict.items():
            canon = _LOOKUP_INDEX.get(k.strip().upper())
            if canon:
                ref_canonical[canon] = float(v)

    # Pre-resolve expression dict to canonical symbols (last value wins for dups)
    expr_canonical: Dict[str, float] = {}
    for k, v in expression_dict.items():
        canon = _LOOKUP_INDEX.get(k.strip().upper())
        if canon:
            expr_canonical[canon] = float(v)

    for gene, info in CMA_NETWORK.items():
        expr_val = expr_canonical.get(gene)

        if expr_val is None:
            missing_genes.append(gene)
            continue

        # Fold-change relative to reference
        if reference_dict is not None:
            ref_val = ref_canonical.get(gene)
            if ref_val is None:
                warnings_out.append(
                    f"{gene}: present in expression_dict but missing from "
                    f"reference_dict; excluded from fold-change score."
                )
                missing_genes.append(gene)
                continue
            if ref_val == 0.0:
                warnings_out.append(
                    f"{gene}: reference value is 0; skipped to avoid division by zero."
                )
                missing_genes.append(gene)
                continue
            expr_val = expr_val / ref_val

        w = info["weight"]
        d = info["direction"]
        contribution = w * d * expr_val
        total_weighted += contribution
        total_weight   += w

        scored_genes.append({
            "gene_symbol":  gene,
            "category":     info["category"],
            "direction":    d,
            "weight":       w,
            "expression":   expr_val,
            "contribution": contribution,
        })

    if total_weight == 0.0:
        cma_score = None
        warnings_out.append(
            "No CMA network genes were found in expression_dict; score is None."
        )
    else:
        cma_score = total_weighted / total_weight

    return {
        "cma_score":       cma_score,
        "n_genes_scored":  len(scored_genes),
        "n_genes_missing": len(missing_genes),
        "scored_genes":    scored_genes,
        "missing_genes":   missing_genes,
        "warnings":        warnings_out,
        "formula_ref": (
            "Bourdenx et al. 2021, Cell 184:2696. "
            "DOI: 10.1016/j.cell.2021.03.048. "
            "score = sum(weight * direction * expr) / sum(weight). "
            "LAMP-2A weight=2; all others weight=1."
        ),
    }


# ===============================================================================
#  CLI HELPERS
# ===============================================================================

def _fmt_dir(direction: int) -> str:
    return "ACTIVATES CMA (+1)" if direction == 1 else "INHIBITS CMA (-1)"


def _print_membership(result: Optional[Dict]) -> None:
    if result is None:
        print("  NOT FOUND -- not a member of the CMA network.")
        return
    print(f"  Gene symbol  : {result['gene_symbol']}")
    print(f"  Protein name : {result['protein_name']}")
    print(f"  UniProt      : {result['uniprot_human'] or 'N/A (mitochondria-encoded)'}")
    print(f"  Category     : {result['category']}")
    print(f"  Direction    : {_fmt_dir(result['direction'])}")
    print(f"  Weight       : {result['weight']}")
    print(f"  Aliases      : {', '.join(result['aliases'])}")
    print(f"  Notes        : {result['notes']}")


def _print_network_table() -> None:
    genes = list_cma_network()
    header = f"  {'Gene':<12} {'Category':<28} {'Dir':>4} {'Wt':>3}  {'UniProt':<10}  Protein"
    print(header)
    print("  " + "-" * 90)
    current_cat = None
    for g in genes:
        if g["category"] != current_cat:
            current_cat = g["category"]
            print(f"\n  [{current_cat.upper().replace('_', ' ')}]")
        dir_str = "+1" if g["direction"] == 1 else "-1"
        uniprot = g["uniprot_human"] or "N/A"
        print(
            f"  {g['gene_symbol']:<12} {g['category']:<28} {dir_str:>4} "
            f"{g['weight']:>3}  {uniprot:<10}  {g['protein_name']}"
        )


def _print_score_result(result: Dict) -> None:
    if result["cma_score"] is not None:
        print(f"\n  CMA Activation Score : {result['cma_score']:+.4f}")
    else:
        print("\n  CMA Activation Score : N/A")
    print(f"  Genes scored         : {result['n_genes_scored']} / {len(CMA_NETWORK)}")
    print(f"  Missing genes        : {', '.join(result['missing_genes']) or 'none'}")
    if result["warnings"]:
        print("\n  Warnings:")
        for w in result["warnings"]:
            print(f"    - {w}")
    print("\n  Per-gene contributions:")
    print(f"    {'Gene':<12} {'Dir':>4} {'Wt':>3}  {'Expr':>9}  {'Contribution':>13}")
    print("    " + "-" * 50)
    for g in result["scored_genes"]:
        dir_str = "+1" if g["direction"] == 1 else "-1"
        print(
            f"    {g['gene_symbol']:<12} {dir_str:>4} {g['weight']:>3}"
            f"  {g['expression']:>9.4f}  {g['contribution']:>13.4f}"
        )
    print(f"\n  Formula: {result['formula_ref']}")


# ===============================================================================
#  SELF-TEST
#  Run: python3 cma_network.py --self-test
# ===============================================================================

def _self_test() -> None:
    """
    Validate check_cma_network_membership(), list_cma_network(), and
    calculate_cma_score() against known values.
    Raises AssertionError on failure.
    """
    passed = failed = 0

    def _check(label: str, condition: bool, detail: str = "") -> None:
        nonlocal passed, failed
        if condition:
            print(f"  PASS  {label}")
            passed += 1
        else:
            print(f"  FAIL  {label}  {detail}")
            failed += 1

    # ── Membership lookup ────────────────────────────────────────────────────
    print("-" * 60)
    print("  MEMBERSHIP LOOKUP TESTS")
    print("-" * 60)

    r = check_cma_network_membership("LAMP2")
    _check("LAMP2 -> in network",            r is not None)
    _check("LAMP2 -> category effector",     r is not None and r["category"] == "effector")
    _check("LAMP2 -> direction +1",          r is not None and r["direction"] == +1)
    _check("LAMP2 -> weight 2",              r is not None and r["weight"] == 2)
    _check("LAMP2 -> UniProt P13473",        r is not None and r["uniprot_human"] == "P13473")

    r = check_cma_network_membership("LAMP2A")          # alias
    _check("LAMP2A alias -> resolves LAMP2", r is not None and r["gene_symbol"] == "LAMP2")

    r = check_cma_network_membership("P13473")          # UniProt
    _check("P13473 UniProt -> resolves LAMP2", r is not None and r["gene_symbol"] == "LAMP2")

    r = check_cma_network_membership("HSC70")           # alias for HSPA8
    _check("HSC70 alias -> resolves HSPA8", r is not None and r["gene_symbol"] == "HSPA8")

    r = check_cma_network_membership("akt1")            # lowercase
    _check("akt1 lowercase -> in network",   r is not None)
    _check("AKT1 -> direction -1",           r is not None and r["direction"] == -1)
    _check("AKT1 -> lysosomal_modulator",    r is not None and r["category"] == "lysosomal_modulator")

    r = check_cma_network_membership("PKB")             # alias for AKT1
    _check("PKB alias -> resolves AKT1",     r is not None and r["gene_symbol"] == "AKT1")

    r = check_cma_network_membership("RICTOR")
    _check("RICTOR -> inhibits CMA (-1)",    r is not None and r["direction"] == -1)

    r = check_cma_network_membership("NRF2")            # alias for NFE2L2
    _check("NRF2 alias -> resolves NFE2L2",  r is not None and r["gene_symbol"] == "NFE2L2")

    r = check_cma_network_membership("TP53")
    _check("TP53 -> not in network",         r is None)

    r = check_cma_network_membership("GAPDH")
    _check("GAPDH -> not in network",        r is None)

    # ── list_cma_network ─────────────────────────────────────────────────────
    print()
    print("-" * 60)
    print("  list_cma_network() TESTS")
    print("-" * 60)

    network = list_cma_network()
    _check(f"Returns all {len(CMA_NETWORK)} genes",
           len(network) == len(CMA_NETWORK))
    _check("All entries have gene_symbol",
           all("gene_symbol" in g for g in network))
    _check("LAMP2 present",
           any(g["gene_symbol"] == "LAMP2" for g in network))
    effectors = [g for g in network if g["category"] == "effector"]
    _check("Exactly 4 effectors",  len(effectors) == 4)
    lys_mods = [g for g in network if g["category"] == "lysosomal_modulator"]
    _check("Exactly 8 lysosomal modulators", len(lys_mods) == 8)
    extra_mods = [g for g in network if g["category"] == "extralysosomal_modulator"]
    _check("Exactly 4 extra-lysosomal modulators", len(extra_mods) == 4)

    # ── calculate_cma_score ──────────────────────────────────────────────────
    print()
    print("-" * 60)
    print("  CMA SCORE CALCULATION TESTS")
    print("-" * 60)

    # All activators -> positive score
    r = calculate_cma_score({"LAMP2": 5.0, "HSPA8": 3.0, "GFAP": 2.0})
    _check("All-activator input -> positive score",
           r["cma_score"] is not None and r["cma_score"] > 0)
    _check("All-activator: 3 genes scored", r["n_genes_scored"] == 3)

    # All inhibitors -> negative score
    r = calculate_cma_score({"AKT1": 4.0, "RICTOR": 3.0, "CTSA": 2.0})
    _check("All-inhibitor input -> negative score",
           r["cma_score"] is not None and r["cma_score"] < 0)

    # Alias keys resolved
    r = calculate_cma_score({"HSC70": 3.0, "PKB": 2.0})   # HSPA8 +1, AKT1 -1
    _check("Alias keys (HSC70, PKB) resolved correctly",
           r["n_genes_scored"] == 2)

    # LAMP2 weight=2 exact arithmetic
    # score = (2*1*1.0 + 1*(-1)*1.0) / (2+1) = 1/3
    r = calculate_cma_score({"LAMP2": 1.0, "AKT1": 1.0})
    expected = (2 * 1 * 1.0 + 1 * (-1) * 1.0) / (2 + 1)
    _check("LAMP2 weight=2 exact arithmetic",
           r["cma_score"] is not None and abs(r["cma_score"] - expected) < 1e-9,
           f"got {r['cma_score']}, expected {expected:.9f}")
    _check("LAMP2+AKT1: LAMP2 weight dominates -> score > 0",
           r["cma_score"] is not None and r["cma_score"] > 0)

    # Empty input
    r = calculate_cma_score({})
    _check("Empty expression_dict -> cma_score is None",  r["cma_score"] is None)
    _check("Empty input -> 0 genes scored",               r["n_genes_scored"] == 0)

    # Non-network genes silently ignored
    r = calculate_cma_score({"LAMP2": 5.0, "TP53": 9.9, "BRCA1": 8.8})
    _check("Non-network genes silently ignored",           r["n_genes_scored"] == 1)
    _check("Non-network ignored: score uses only LAMP2",   r["cma_score"] is not None and r["cma_score"] > 0)

    # Fold-change relative to reference
    r = calculate_cma_score({"LAMP2": 6.0, "HSPA8": 3.0},
                             reference_dict={"LAMP2": 3.0, "HSPA8": 3.0})
    # LAMP2 fold=2.0 (weight=2, dir=+1), HSPA8 fold=1.0 (weight=1, dir=+1)
    # score = (2*1*2.0 + 1*1*1.0) / (2+1) = 5/3
    expected_fc = (2 * 1 * 2.0 + 1 * 1 * 1.0) / (2 + 1)
    _check(f"Fold-change score LAMP2x2 + HSPA8x1 = {expected_fc:.4f}",
           r["cma_score"] is not None and abs(r["cma_score"] - expected_fc) < 1e-9,
           f"got {r['cma_score']}")

    # Missing reference gene -> warning + gene excluded
    r = calculate_cma_score({"LAMP2": 6.0, "HSPA8": 3.0},
                             reference_dict={"LAMP2": 3.0})   # HSPA8 missing from ref
    _check("Missing ref gene generates warning",   len(r["warnings"]) > 0)
    _check("Missing ref gene excluded from score", r["n_genes_scored"] == 1)

    # ── Summary ──────────────────────────────────────────────────────────────
    print()
    bar = "=" * 60
    total = passed + failed
    print(bar)
    if failed == 0:
        print(f"  ALL {total} TESTS PASSED -- cma_network.py v1.0 READY")
    else:
        print(f"  {passed}/{total} passed -- {failed} FAILED")
    print(bar)
    if failed:
        raise AssertionError(f"{failed} self-test(s) failed.")


# ===============================================================================
#  CLI ENTRY POINT
# ===============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="CMA Network Gene Lookup and Scoring  (cma_network.py v1.0)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--check", metavar="GENE",
        help="Check CMA network membership of a gene symbol or UniProt ID",
    )
    parser.add_argument(
        "--list-network", action="store_true",
        help="Print all 16 CMA network genes as a table",
    )
    parser.add_argument(
        "--score", metavar="JSON",
        help=(
            "Calculate CMA activation score from a JSON dict of "
            "log-normalised expression values, e.g. "
            "'{\"LAMP2\": 5.2, \"HSPA8\": 3.1, \"AKT1\": 1.8}'"
        ),
    )
    parser.add_argument(
        "--reference", metavar="JSON",
        help=(
            "Optional reference (control) expression dict for fold-change scoring, "
            "e.g. '{\"LAMP2\": 3.0, \"HSPA8\": 2.5}'"
        ),
    )
    parser.add_argument(
        "--self-test", action="store_true",
        help="Run the built-in validation suite",
    )
    args = parser.parse_args()

    if args.self_test:
        _self_test()
        return

    if args.check:
        print(f"\nCMA Network Membership: {args.check}")
        print("-" * 55)
        _print_membership(check_cma_network_membership(args.check))
        print()
        return

    if args.list_network:
        print(
            "\nCMA Network Genes  (Kirchner et al. 2019, PLoS Biol 17:e3000301)"
        )
        print("=" * 95)
        _print_network_table()
        print(f"\n  Total: {len(CMA_NETWORK)} genes\n")
        return

    if args.score:
        try:
            expr = json.loads(args.score)
        except json.JSONDecodeError as e:
            print(f"ERROR: invalid JSON in --score: {e}", file=sys.stderr)
            sys.exit(1)

        ref = None
        if args.reference:
            try:
                ref = json.loads(args.reference)
            except json.JSONDecodeError as e:
                print(f"ERROR: invalid JSON in --reference: {e}", file=sys.stderr)
                sys.exit(1)

        result = calculate_cma_score(expr, reference_dict=ref)
        _print_score_result(result)
        return

    parser.print_help()


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        _self_test()
    else:
        main()