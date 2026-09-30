"""
motif_finder/__init__.py
Re-exports the public API from motif_finder/motif_finder.py so the package can
be imported flat:  from motif_finder import find_all_degradation_motifs, ...
"""
from motif_finder.motif_finder import (
    find_all_degradation_motifs,
    find_kferq_motifs,
    kferq_to_sea_motifs,
    find_lir_motifs,
    find_dbox_motifs,
    find_ken_motifs,
    find_pest_sequences,
    classify_protein_cma,
    classify_n_degron,
    find_cdegron_motifs,
    parse_fasta,
)

__all__ = [
    "find_all_degradation_motifs",
    "find_kferq_motifs",
    "kferq_to_sea_motifs",
    "find_lir_motifs",
    "find_dbox_motifs",
    "find_ken_motifs",
    "find_pest_sequences",
    "classify_protein_cma",
    "classify_n_degron",
    "find_cdegron_motifs",
    "parse_fasta",
]
