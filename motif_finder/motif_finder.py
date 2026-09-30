#!/usr/bin/env python3
"""
motif_finder.py  v2.0
---------------------
Identifies protein degradation motifs in amino acid sequences.
Self-contained: no runtime web or API calls.

Engines
-------
  1. KFERQ-like (CMA / chaperone-mediated autophagy)
  2. LIR / AIM  (LC3-interacting region / selective autophagy)
  3. D-box       (RxxL — APC/C degron)
  4. KEN box     (KEN — APC/C degron)
  5. N-degron    (N-terminal residue classification; N-end rule)
  6. C-degron    (APPBP2 and DCAF12 C-terminal degrons)
  7. PEST        (hydrophilic rapid-turnover signal)

Unified entry point:
  find_all_degradation_motifs(seq) → dict with all engine results

All motif dicts share the standard format:
  {start: int, end: int, motif: str, type: str, notes: str}
  (start and end are 1-based inclusive)

CLI usage:
  python3 motif_finder.py sequences.fasta
  python3 motif_finder.py --sequence MSQIRLKEFERQSLP...
  python3 motif_finder.py sequences.fasta --no-phospho --n-bearing

Import usage:
  from motif_finder import (find_kferq_motifs, find_lir_motifs,
                             find_dbox_motifs, find_ken_motifs,
                             classify_n_degron, find_cdegron_motifs,
                             find_pest_sequences,
                             find_all_degradation_motifs)
"""

import re
import sys
import argparse
from typing import Dict, List, Optional, Tuple


# ═══════════════════════════════════════════════════════════════════════════════
#  KFERQ-LIKE MOTIF ENGINE
#
#  Source
#  ------
#  Kirchner P et al. (2019) "Proteome-wide analysis of chaperone-mediated
#  autophagy targeting motifs." PLoS Biology 17(6): e3000301.
#  DOI: 10.1371/journal.pbio.3000301
# ═══════════════════════════════════════════════════════════════════════════════

BASIC       = frozenset("KR")
HYDROPHOBIC = frozenset("ILVF")
ACIDIC      = frozenset("DE")
PHOS_STY    = frozenset("STY")


def _validate_core(core4: str, allow_phospho: bool = False) -> Tuple[bool, Optional[str]]:
    """
    Validate the 4 non-flanking positions of a candidate KFERQ-like pentapeptide.

    Rules (Kirchner et al. 2019):
      • Every residue must be basic | hydrophobic | acidic (no neutral filler)
      • Basic (K, R):           count = 1 or 2
      • Hydrophobic (I,L,V,F):  count = 1 or 2
      • Acidic (D, E):          count = exactly 1
        — OR — S/T/Y when allow_phospho=True

    Returns (is_valid, acid_subtype)
    """
    n_basic = n_hydro = n_acid = n_phos = 0
    for aa in core4:
        if   aa in BASIC:                         n_basic += 1
        elif aa in HYDROPHOBIC:                   n_hydro += 1
        elif aa in ACIDIC:                        n_acid  += 1
        elif allow_phospho and aa in PHOS_STY:    n_phos  += 1
        else:
            return False, None

    total_acid = n_acid + n_phos
    if not (1 <= n_basic <= 2): return False, None
    if not (1 <= n_hydro <= 2): return False, None
    if total_acid != 1:         return False, None
    return True, ("phospho" if n_phos else "canonical")


def find_kferq_motifs(
    sequence: str,
    allow_phospho:   bool = True,
    allow_acetyl:    bool = True,
    allow_n_bearing: bool = False,
) -> List[Dict]:
    """
    Scan for KFERQ-like (CMA targeting) motifs.

    Q-flanked canonical  : [QKILVF][KR]{1,2}[ILVF]{1,2}[DE][ILVF]{1,2}[QKILVF]
    Q-flanked phospho    : same but S/T/Y occupies the acidic slot
    K-flanked acetyl     : K replaces Q at flank; requires K acetylation in vivo
    N-flanked (advanced) : N replaces Q; necessary but NOT sufficient alone

    Returns list of standard motif dicts sorted by start position.
    """
    seq  = sequence.upper().replace(" ", "").replace("\n", "")
    hits = []
    seen = set()

    for i in range(len(seq) - 4):
        pep = seq[i: i + 5]

        for flank_pos in (0, 4):
            flanking = pep[flank_pos]
            core4    = pep[1:] if flank_pos == 0 else pep[:4]

            # Q-flanked -------------------------------------------------------
            if flanking == "Q":
                ok, sub = _validate_core(core4, allow_phospho=False)
                if ok:
                    key = (i + 1, "canonical")
                    if key not in seen:
                        seen.add(key)
                        hits.append(dict(
                            start=i + 1, end=i + 5, motif=pep,
                            type="canonical",
                            notes="Q-flanked | D/E acidic | constitutive CMA signal"
                        ))
                if allow_phospho:
                    ok, sub = _validate_core(core4, allow_phospho=True)
                    if ok and sub == "phospho":
                        key = (i + 1, "phospho_activated")
                        if key not in seen:
                            seen.add(key)
                            hits.append(dict(
                                start=i + 1, end=i + 5, motif=pep,
                                type="phospho_activated",
                                notes="Q-flanked | S/T/Y in acidic slot | requires phosphorylation"
                            ))

            # K-flanked (acetylation-activated) --------------------------------
            elif flanking == "K" and allow_acetyl:
                ok, sub = _validate_core(core4, allow_phospho=False)
                if ok and sub == "canonical":
                    key = (i + 1, "acetyl_activated")
                    if key not in seen:
                        seen.add(key)
                        hits.append(dict(
                            start=i + 1, end=i + 5, motif=pep,
                            type="acetyl_activated",
                            notes=(
                                f"K at position {i + flank_pos + 1} replaces Q | "
                                "requires acetylation of that K"
                            )
                        ))

            # N-flanked (advanced) --------------------------------------------
            elif flanking == "N" and allow_n_bearing:
                ok, sub = _validate_core(core4, allow_phospho=False)
                if ok and sub == "canonical":
                    key = (i + 1, "n_bearing")
                    if key not in seen:
                        seen.add(key)
                        hits.append(dict(
                            start=i + 1, end=i + 5, motif=pep,
                            type="n_bearing",
                            notes=(
                                "N replaces Q (advanced) | necessary but NOT sufficient "
                                "for HSC70 binding alone"
                            )
                        ))

    return sorted(hits, key=lambda x: x["start"])


def classify_protein_cma(motifs: List[Dict]) -> str:
    """
    Assign the highest-priority CMA substrate category (Kirchner et al. 2019):
    canonical > phospho_activated > acetyl_activated > n_bearing > no_motif
    """
    types = {m["type"] for m in motifs}
    for t in ("canonical", "phospho_activated", "acetyl_activated", "n_bearing"):
        if t in types:
            return t
    return "no_motif"


# ═══════════════════════════════════════════════════════════════════════════════
#  LIR / AIM MOTIF ENGINE
#
#  Sources
#  -------
#  Core [WFY]-x-x-[ILV] consensus:
#    Cheng X et al. (2023) bioRxiv DOI: 10.1101/2022.09.25.509395
#    Karin M et al. (2014) figshare  DOI: 10.6084/m9.figshare.1601912.v1
#    McLaughlin M et al. (2020) Biochem J DOI: 10.1042/bcj20200714
#  xLIR: upstream acidic (D/E) and phospho-activatable (S/T) flanking:
#    Kohler A et al. (2022) bioRxiv DOI: 10.1101/2022.02.11.480038
#    McLaughlin M et al. (2020) Biochem J DOI: 10.1042/bcj20200714
#  iLIR tool (xLIR + PSSM):
#    Klionsky DJ et al. (2014) Autophagy DOI: 10.4161/AUTO.28260
#  W > F > Y affinity hierarchy / GABARAPL1 preference:
#    Rozenknop A et al. (2011) J Mol Biol DOI: 10.1016/j.jmb.2011.05.003
#  Non-canonical variants (2HP-LIR, [DE]W[DE]-LIR, HP0-LIR):
#    AlphaFold2 SLiM screen (2025) figshare DOI: 10.6084/m9.figshare.28929437.v1
# ═══════════════════════════════════════════════════════════════════════════════

_LIR_AROMATIC  = frozenset("WFY")
_LIR_ALIPHATIC = frozenset("ILV")
_LIR_ACID_UP   = frozenset("DE")    # upstream acidic (xLIR-acidic)
_LIR_PHOS_UP   = frozenset("ST")    # upstream phospho-activatable (xLIR-phospho)

_LIR_CORE_RE = re.compile(r"(?=([WFY]..[ILV]))")   # overlapping 4-mer scan


def find_lir_motifs(sequence: str) -> List[Dict]:
    """
    Scan for LIR/AIM (LC3-interacting region) motifs.

    Core: [W/F/Y]-x-x-[I/L/V]  (positions Θ₀-x₁-x₂-Γ₃)

    Each hit is classified by its upstream window (positions -3, -2, -1 relative
    to the aromatic position):

    lir_canonical      No acidic or phospho-activatable residues upstream.
    lir_xlir_acidic    ≥1 D/E upstream; stronger constitutive binding.
    lir_xlir_phospho   ≥1 S/T upstream (no D/E); phospho-activatable binding.
    lir_xlir_optimal   Both D/E and S/T upstream; highest predicted affinity.

    The [DE]W[DE] non-canonical variant flag is appended as "+DEWDE" where
    applicable (AF2 SLiM screen; DOI: 10.6084/m9.figshare.28929437.v1).

    Note: W at position 0 gives highest affinity (W > F > Y; Rozenknop 2011).
    For full PSSM-weighted xLIR scoring use iLIR (DOI: 10.4161/AUTO.28260).
    """
    seq  = sequence.upper().replace(" ", "").replace("\n", "")
    hits = []

    for m in _LIR_CORE_RE.finditer(seq):
        i   = m.start()         # 0-based position of W/F/Y
        pep = m.group(1)        # 4-residue core [WFY]xx[ILV]

        aromatic  = pep[0]      # W, F, or Y
        aliphatic = pep[3]      # I, L, or V

        # Upstream residues: positions -3, -2, -1 relative to aromatic position
        upstream = seq[max(0, i - 3): i]
        n_acidic  = sum(1 for aa in upstream if aa in _LIR_ACID_UP)
        n_phos    = sum(1 for aa in upstream if aa in _LIR_PHOS_UP)

        # [DE]W[DE] non-canonical check (only when aromatic is W)
        is_dewde = (
            aromatic == "W"
            and len(upstream) >= 1
            and upstream[-1] in _LIR_ACID_UP          # -1 position D/E
            and i + 1 < len(seq)
            and seq[i + 1] in _LIR_ACID_UP            # +1 position D/E
        )

        # Classification
        if n_acidic >= 1 and n_phos >= 1:
            lir_type = "lir_xlir_optimal"
            flanking = f"upstream D/E ({n_acidic}) + S/T ({n_phos}); high-affinity xLIR"
        elif n_acidic >= 1:
            lir_type = "lir_xlir_acidic"
            flanking = f"upstream D/E ({n_acidic}); constitutive xLIR, stronger binding"
        elif n_phos >= 1:
            lir_type = "lir_xlir_phospho"
            flanking = f"upstream S/T ({n_phos}); phospho-activatable xLIR"
        else:
            lir_type = "lir_canonical"
            flanking = "no acidic/phospho flanking in -3 to -1"

        if is_dewde:
            lir_type += "+DEWDE"

        affinity = {"W": "high (W)", "F": "intermediate (F)", "Y": "lower (Y)"}[aromatic]

        hits.append(dict(
            start = i + 1,
            end   = i + 4,
            motif = pep,
            type  = lir_type,
            notes = (
                f"core {aromatic}xx{aliphatic} | aromatic affinity: {affinity} | {flanking}"
            ),
        ))

    return hits


# ═══════════════════════════════════════════════════════════════════════════════
#  D-BOX ENGINE
#
#  Sources
#  -------
#  Minimal RxxL (Destruction box):
#    Glotzer M et al. (1991) Nature 349:132–138  [original D-box]
#    King RW et al. (1996) Science 274:1652–1659
#  Extended RxxLxxxxN:
#    Andrews PD et al. (2024) bioRxiv DOI: 10.1101/2024.04.30.590460
#
#  WARNING: RxxL occurs very frequently by chance. Biological interpretation
#  requires APC/C substrate context (cell-cycle proteins, disordered regions).
# ═══════════════════════════════════════════════════════════════════════════════

_DBOX_EXT_RE = re.compile(r"(?=(R.{2}L.{4}N))")   # RxxLxxxxN  (9 residues)
_DBOX_MIN_RE = re.compile(r"(?=(R.{2}L))")         # RxxL       (4 residues)


def find_dbox_motifs(sequence: str, extended: bool = True) -> List[Dict]:
    """
    Scan for D-box (Destruction box) degrons.

    Parameters
    ----------
    extended : If True (default), RxxLxxxxN hits (higher specificity) are
               reported as 'dbox_extended'; overlapping RxxL are suppressed.
               Remaining RxxL-only hits are reported as 'dbox_minimal'.

    Returns list of standard motif dicts sorted by start position.
    """
    seq  = sequence.upper().replace(" ", "").replace("\n", "")
    hits = []
    seen_ext: set = set()

    # Extended RxxLxxxxN first (preferred, higher specificity)
    if extended:
        for m in _DBOX_EXT_RE.finditer(seq):
            i = m.start()
            seen_ext.add(i)
            hits.append(dict(
                start = i + 1,
                end   = i + 9,
                motif = m.group(1),
                type  = "dbox_extended",
                notes = "RxxLxxxxN | APC/C-Cdc20 degron (Andrews et al. 2024)"
            ))

    # Minimal RxxL (skip positions already covered by extended)
    for m in _DBOX_MIN_RE.finditer(seq):
        i = m.start()
        if i in seen_ext:
            continue
        hits.append(dict(
            start = i + 1,
            end   = i + 4,
            motif = m.group(1),
            type  = "dbox_minimal",
            notes = "RxxL | minimal D-box; APC/C-Cdc20 (Glotzer et al. 1991)"
        ))

    return sorted(hits, key=lambda x: x["start"])


# ═══════════════════════════════════════════════════════════════════════════════
#  KEN BOX ENGINE
#
#  Sources
#  -------
#  Pfleger CM & Kirschner MW (2000) Genes Dev 14:655–665  [original KEN box]
#  Andrews PD et al. (2024) bioRxiv DOI: 10.1101/2024.04.30.590460
#    "KEN motif docks to the top surface of the WD40 propeller of the co-activator"
# ═══════════════════════════════════════════════════════════════════════════════

_KEN_RE = re.compile(r"(?=(KEN))")


def find_ken_motifs(sequence: str) -> List[Dict]:
    """
    Scan for KEN box motifs (exact tripeptide K-E-N).
    Recognised by APC/C-Cdh1 (and in some substrates by APC/C-Cdc20).
    """
    seq  = sequence.upper().replace(" ", "").replace("\n", "")
    hits = []
    for m in _KEN_RE.finditer(seq):
        i = m.start()
        hits.append(dict(
            start = i + 1,
            end   = i + 3,
            motif = m.group(1),
            type  = "ken_box",
            notes = "KEN | APC/C-Cdh1 degron (Pfleger & Kirschner 2000)"
        ))
    return hits


# ═══════════════════════════════════════════════════════════════════════════════
#  N-DEGRON ENGINE
#
#  Sources
#  -------
#  Timms RT & Koren I (2020) Biochem Soc Trans DOI: 10.1042/BST20191094
#  Dissmeyer N et al. (2016) Curr Protein Pept Sci
#    DOI: 10.2174/0929866523666160108115809
#
#  N-end rule hierarchy (mammals)
#  ───────────────────────────────
#  Stabilising:                 M, A, G, V, C, P, S, T
#  Type I primary destabilising (basic):         R, H, K
#    → directly recognised by UBR1/UBR2 type I binding site
#  Type II primary destabilising (aromatic/bulky): F, Y, W, L, I
#    → directly recognised by UBR1/UBR2 type II binding site
#  Secondary destabilising (acidic): D, E
#    → arginylated by ATE1 arginyltransferase → R at N-term (type I)
#  Tertiary destabilising:           N, Q
#    → deamidated by NTAN1/NTAQ1 → D/E → arginylated → R (type I)
#
#  Methionine aminopeptidase (MAP) cleavage: MAP removes initiator Met when
#  the penultimate residue is A, C, G, P, S, T, or V (small/neutral).
# ═══════════════════════════════════════════════════════════════════════════════

_NDEG_STABLE   = frozenset("MAGVCPST")
_NDEG_PRIM_I   = frozenset("RHK")    # basic
_NDEG_PRIM_II  = frozenset("FYWLI")  # aromatic / bulky hydrophobic
_NDEG_SECOND   = frozenset("DE")     # acidic; arginylated
_NDEG_TERTIARY = frozenset("NQ")     # deamidated
_MAP_CLEAVAGE  = frozenset("ACGPSTV")  # penultimate residues triggering MAP


def _classify_ndeg_residue(aa: str) -> str:
    if aa in _NDEG_STABLE:    return "stabilising"
    if aa in _NDEG_PRIM_I:    return "type_I_primary_destabilising"
    if aa in _NDEG_PRIM_II:   return "type_II_primary_destabilising"
    if aa in _NDEG_SECOND:    return "secondary_destabilising"
    if aa in _NDEG_TERTIARY:  return "tertiary_destabilising"
    return "unclassified"


def classify_n_degron(sequence: str) -> Dict:
    """
    Classify the N-terminal residue according to the N-degron (N-end rule) pathway.

    Also predicts Met aminopeptidase (MAP) cleavage: if N-terminal Met is followed
    by a small residue (A/C/G/P/S/T/V) the initiator Met is likely removed,
    exposing the penultimate residue as the effective N-degron.

    Returns a single standard motif dict (start=1, end=1).
    """
    seq = sequence.upper().replace(" ", "").replace("\n", "")
    if not seq:
        return {}

    n_term = seq[0]
    n_class = _classify_ndeg_residue(n_term)

    # MAP cleavage prediction
    map_note = ""
    effective_class = n_class
    if n_term == "M" and len(seq) >= 2 and seq[1] in _MAP_CLEAVAGE:
        exposed = seq[1]
        effective_class = _classify_ndeg_residue(exposed)
        map_note = (
            f" | MAP likely exposes {exposed} → effective class: {effective_class}"
        )

    return dict(
        start  = 1,
        end    = 1,
        motif  = n_term,
        type   = f"n_degron_{effective_class}",
        notes  = (
            f"N-terminal {n_term}: {n_class}{map_note} | "
            "type I (R/H/K)=basic→UBR; type II (F/Y/W/L/I)=hydrophobic→UBR; "
            "secondary (D/E)→ATE1 arginylation; tertiary (N/Q)→deamidation"
        ),
    )


# ═══════════════════════════════════════════════════════════════════════════════
#  C-DEGRON ENGINE
#
#  Source
#  ------
#  Timms RT & Koren I (2020) Biochem Soc Trans DOI: 10.1042/BST20191094
#    APPBP2/VHL-type: Rx(2–4)Gx(0–3) at C-terminus (optimal: RxxGx, RxxGxx)
#    DCAF12-type:     twin C-terminal glutamic acid (−EE)
# ═══════════════════════════════════════════════════════════════════════════════

_CDEG_APPBP2_RE = re.compile(r"R.{2,4}G.{0,3}$")   # Rx(2-4)Gx(0-3) anchored at end
_CDEG_DCAF12_RE = re.compile(r"EE$")                # -EE at extreme C-terminus


def find_cdegron_motifs(sequence: str) -> List[Dict]:
    """
    Scan the C-terminal region for C-degron motifs.

    APPBP2 type: Rx(2–4)Gx(0–3) anchored to C-terminus.
    DCAF12 type: -EE at the very C-terminus.

    Note: C-degrons are end-specific; only terminal residues are checked.
    """
    seq  = sequence.upper().replace(" ", "").replace("\n", "")
    hits = []

    m = _CDEG_APPBP2_RE.search(seq)
    if m:
        hits.append(dict(
            start = m.start() + 1,
            end   = len(seq),
            motif = m.group(0),
            type  = "cdegron_appbp2",
            notes = "Rx(2-4)Gx(0-3) at C-terminus | APPBP2/VHL C-degron (Timms & Koren 2020)"
        ))

    m = _CDEG_DCAF12_RE.search(seq)
    if m:
        hits.append(dict(
            start = m.start() + 1,
            end   = len(seq),
            motif = m.group(0),
            type  = "cdegron_dcaf12",
            notes = "-EE at C-terminus | DCAF12 C-degron (Timms & Koren 2020)"
        ))

    return hits


# ═══════════════════════════════════════════════════════════════════════════════
#  PEST SEQUENCE ENGINE
#
#  Sources
#  -------
#  Rogers S et al. (1986) Science 234:364–368   [original PEST definition]
#  Rechsteiner M & Rogers SW (1996) Trends Biochem Sci 21:267–271
#
#  Structural criteria (Rogers 1986):
#    1. Flanked by K, R, or H at both ends.
#    2. Contains ≥1 P and ≥1 D or E within the region.
#    3. No K, R, or H inside the region.
#    4. Region length ≥ min_length (default 12) residues (exclusive of flanks).
#
#  Score (simplified fraction approach as used in EMBOSS epestfind):
#    PEST fraction = (count of P + E + S + T) / region_length
#    Regions with fraction ≥ 0.30 → "pest_positive"
#    Regions meeting structural criteria but below threshold → "pest_putative"
#
#  Note: the original 1986 formula uses a hydrophilicity-weighted score; this
#  implementation uses the simpler fraction threshold, which is the basis for
#  most modern PEST scanners.
# ═══════════════════════════════════════════════════════════════════════════════

_PEST_FLANK = frozenset("KRH")
_PEST_CHARS = frozenset("PEST")


def find_pest_sequences(
    sequence: str,
    min_length: int = 12,
    max_length: int = 100,
    min_pest_fraction: float = 0.30,
) -> List[Dict]:
    """
    Scan for PEST sequences (signals associated with rapid protein turnover).

    Parameters
    ----------
    min_length       : minimum PEST region length (default 12, exclusive of flanks)
    max_length       : maximum region length to consider (default 100; prevents
                       spanning the whole protein between distant K/R/H residues)
    min_pest_fraction: minimum (P+E+S+T)/length fraction for 'pest_positive'

    Returns list of standard motif dicts. Overlapping/nested regions may appear;
    filter by fraction or length if needed.
    """
    seq   = sequence.upper().replace(" ", "").replace("\n", "")
    flank = [i for i, aa in enumerate(seq) if aa in _PEST_FLANK]
    hits  = []
    seen  = set()

    for li in flank:
        for ri in flank:
            if ri <= li:
                continue
            region = seq[li + 1: ri]   # exclusive of flanking K/R/H
            rlen   = len(region)

            if rlen < min_length or rlen > max_length:
                continue
            if any(aa in _PEST_FLANK for aa in region):
                continue  # K/R/H inside: invalid
            if "P" not in region:
                continue
            if not any(aa in "DE" for aa in region):
                continue

            frac  = sum(1 for aa in region if aa in _PEST_CHARS) / rlen
            ptype = "pest_positive" if frac >= min_pest_fraction else "pest_putative"

            key = (li + 2, ri)
            if key not in seen:
                seen.add(key)
                hits.append(dict(
                    start = li + 2,           # 1-based, first residue after left flank
                    end   = ri,               # 1-based, last residue before right flank
                    motif = region,
                    type  = ptype,
                    notes = (
                        f"PEST fraction={frac:.2f} ({ptype}) | "
                        f"len={rlen} | "
                        f"flanked by {seq[li]}({li+1})…{seq[ri]}({ri+1})"
                    ),
                ))

    return sorted(hits, key=lambda x: x["start"])


# ═══════════════════════════════════════════════════════════════════════════════
#  UNIFIED DEGRADATION MOTIF API
# ═══════════════════════════════════════════════════════════════════════════════

def find_all_degradation_motifs(
    sequence: str,
    allow_kferq_phospho:   bool  = True,
    allow_kferq_acetyl:    bool  = True,
    allow_kferq_n_bearing: bool  = False,
    dbox_extended:         bool  = True,
    pest_min_length:       int   = 12,
    pest_max_length:       int   = 100,
    pest_min_fraction:     float = 0.30,
) -> Dict:
    """
    Run all degradation motif engines and return a unified structured result.

    Returns
    -------
    {
        'kferq'    : List[Dict],   KFERQ-like CMA motifs
        'lir'      : List[Dict],   LIR/AIM autophagy motifs
        'dbox'     : List[Dict],   D-box degrons
        'ken'      : List[Dict],   KEN box degrons
        'n_degron' : Dict,         N-terminal residue classification
        'c_degron' : List[Dict],   C-terminal degrons
        'pest'     : List[Dict],   PEST rapid-turnover sequences
        'summary'  : Dict,         per-engine counts and boolean flags
    }

    All individual motif dicts use the standard format:
        {start, end, motif, type, notes}   (start and end 1-based inclusive)
    """
    kferq  = find_kferq_motifs(sequence,
                                allow_phospho=allow_kferq_phospho,
                                allow_acetyl=allow_kferq_acetyl,
                                allow_n_bearing=allow_kferq_n_bearing)
    lir    = find_lir_motifs(sequence)
    dbox   = find_dbox_motifs(sequence, extended=dbox_extended)
    ken    = find_ken_motifs(sequence)
    n_deg  = classify_n_degron(sequence)
    c_deg  = find_cdegron_motifs(sequence)
    pest   = find_pest_sequences(sequence,
                                  min_length=pest_min_length,
                                  max_length=pest_max_length,
                                  min_pest_fraction=pest_min_fraction)

    cma_cat   = classify_protein_cma(kferq)
    n_destab  = "destabilising" in n_deg.get("type", "") if n_deg else False
    c_destab  = len(c_deg) > 0

    return dict(
        kferq    = kferq,
        lir      = lir,
        dbox     = dbox,
        ken      = ken,
        n_degron = n_deg,
        c_degron = c_deg,
        pest     = pest,
        summary  = dict(
            cma_category             = cma_cat,
            n_kferq                  = len(kferq),
            n_lir                    = len(lir),
            n_dbox                   = len(dbox),
            n_ken                    = len(ken),
            n_pest_positive          = sum(1 for h in pest if h["type"] == "pest_positive"),
            n_pest_putative          = sum(1 for h in pest if h["type"] == "pest_putative"),
            has_destabilising_n_term = n_destab,
            has_destabilising_c_term = c_destab,
        ),
    )


# ═══════════════════════════════════════════════════════════════════════════════
#  FASTA PARSER
# ═══════════════════════════════════════════════════════════════════════════════

def parse_fasta(text: str) -> List[Tuple[str, str]]:
    """Return [(header, sequence), ...] from FASTA-formatted text."""
    records: List[Tuple[str, str]] = []
    header, parts = None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header, parts = line[1:], []
        elif line:
            parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


# ═══════════════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════════════

def _print_all_results(header: str, results: Dict) -> None:
    bar = "=" * 72
    s   = results["summary"]
    print(f"\n{bar}")
    print(f"  {header}")
    print(f"{bar}")
    print(f"  CMA category : {s['cma_category']}")
    print(f"  Motif counts : KFERQ={s['n_kferq']}  LIR={s['n_lir']}  "
          f"D-box={s['n_dbox']}  KEN={s['n_ken']}  "
          f"PEST+={s['n_pest_positive']} PEST?={s['n_pest_putative']}")

    def _show(label, hits, show_notes=True):
        if not hits:
            return
        print(f"\n  ── {label} ──")
        for h in (hits if isinstance(hits, list) else [hits]):
            print(f"    [{h['start']:>5}–{h['end']:<5}]  {h['motif']}  ({h['type']})")
            if show_notes:
                print(f"             {h['notes']}")

    _show("KFERQ-like (CMA)", results["kferq"])
    _show("LIR / AIM",        results["lir"])
    _show("D-box",            results["dbox"])
    _show("KEN box",          results["ken"])
    if results["n_degron"]:
        _show("N-degron",     results["n_degron"])
    _show("C-degron",         results["c_degron"])
    _show("PEST sequences",   results["pest"])

    if s["has_destabilising_n_term"]:
        print("\n  ⚑  Destabilising N-terminal residue detected.")
    if s["has_destabilising_c_term"]:
        print("  ⚑  C-terminal degron detected.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find protein degradation motifs (KFERQ, LIR, D-box, KEN, N-/C-degron, PEST).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("fasta_file", nargs="?",      help="Path to FASTA file")
    src.add_argument("--sequence", "-s", metavar="SEQ", help="Raw amino-acid sequence")

    parser.add_argument("--phospho",     dest="phospho",    action="store_true",  default=True)
    parser.add_argument("--no-phospho",  dest="phospho",    action="store_false")
    parser.add_argument("--acetyl",      dest="acetyl",     action="store_true",  default=True)
    parser.add_argument("--no-acetyl",   dest="acetyl",     action="store_false")
    parser.add_argument("--n-bearing",   dest="n_bearing",  action="store_true",  default=False,
                        help="Include N-bearing (Q→N) KFERQ variants (default OFF)")
    parser.add_argument("--no-extended-dbox", dest="ext_dbox", action="store_false", default=True,
                        help="Report only minimal D-box (RxxL), not extended (RxxLxxxxN)")

    args = parser.parse_args()

    if args.sequence:
        records = [("command-line sequence", args.sequence)]
    else:
        with open(args.fasta_file) as fh:
            records = parse_fasta(fh.read())
        if not records:
            sys.exit(f"ERROR: no FASTA records found in {args.fasta_file}")

    for header, seq in records:
        results = find_all_degradation_motifs(
            seq,
            allow_kferq_phospho=args.phospho,
            allow_kferq_acetyl=args.acetyl,
            allow_kferq_n_bearing=args.n_bearing,
            dbox_extended=args.ext_dbox,
        )
        _print_all_results(header, results)


# ═══════════════════════════════════════════════════════════════════════════════
#  SELF-TEST
#  Run: python3 motif_finder.py --self-test   OR   python3 -c "import motif_finder; motif_finder._self_test()"
# ═══════════════════════════════════════════════════════════════════════════════

def _self_test() -> None:
    """
    Validate all motif engines against known sequences. Raises AssertionError
    on failure; prints PASS / status for each test.
    """
    passed = failed = 0

    def _check(label, condition, detail=""):
        nonlocal passed, failed
        if condition:
            print(f"  PASS  {label}")
            passed += 1
        else:
            print(f"  FAIL  {label}  {detail}")
            failed += 1

    print("─" * 60)
    print("  KFERQ-LIKE MOTIF TESTS")
    print("─" * 60)

    # Canonical: KFERQ — Q-flanked, F(hydro), E(acid), R(basic), K(basic)→core=FERK → 2 basic, 1 hydro? Wait:
    # KFERQ: K at pos 0 (flank candidate) or Q at pos 4 (flank)
    # Window KFERQ: Q at pos 4 → core = KFER → K(basic), F(hydro), E(acid), R(basic) → 2 basic, 1 hydro, 1 acid → VALID canonical
    hits = find_kferq_motifs("KFERQ")
    _check("KFERQ → canonical",     any(h["type"] == "canonical" for h in hits))
    _check("KFERQ Q-flanked",       any(h["motif"] == "KFERQ" for h in hits))

    # ILKEQ: Q at pos 4, core = ILKE → I(hydro), L(hydro), K(basic), E(acid) → 1 basic, 2 hydro, 1 acid → VALID
    hits = find_kferq_motifs("ILKEQ")
    _check("ILKEQ → canonical",     any(h["type"] == "canonical" for h in hits))

    # VKKDQ: Q at pos 4, core = VKKD → V(hydro), K(basic), K(basic), D(acid) → 2 basic, 1 hydro, 1 acid → VALID
    hits = find_kferq_motifs("VKKDQ")
    _check("VKKDQ → canonical",     any(h["type"] == "canonical" for h in hits))

    # LDRLQ: Q at pos 4, core = LDRL → L(hydro), D(acid), R(basic), L(hydro) → 1 basic, 2 hydro, 1 acid → VALID
    hits = find_kferq_motifs("LDRLQ")
    _check("LDRLQ → canonical",     any(h["type"] == "canonical" for h in hits))

    # QKFSR: Q at pos 0, core = KFSR → K(basic), F(hydro), S(phospho), R(basic) → with phospho: 2 basic, 1 hydro, 1 phospho = canonical+phospho → valid phospho; without allow_phospho=True: invalid
    hits_p  = find_kferq_motifs("QKFSR", allow_phospho=True)
    hits_np = find_kferq_motifs("QKFSR", allow_phospho=False)
    _check("QKFSR → phospho_activated (allow_phospho=True)",
           any(h["type"] == "phospho_activated" for h in hits_p))
    _check("QKFSR → no hit (allow_phospho=False)",
           not any(h["type"] == "canonical" for h in hits_np))

    # KFERR: K at pos 0 (acetyl), core = FERR → F(hydro), E(acid), R(basic), R(basic) → 2 basic, 1 hydro, 1 acid → VALID acetyl
    hits = find_kferq_motifs("KFERR", allow_acetyl=True)
    _check("KFERR → acetyl_activated",  any(h["type"] == "acetyl_activated" for h in hits))

    # QAFERQ: Q at pos 0, core = AFER → A(unclassified) → INVALID
    #         Q at pos 5 in window AFERQ → core = QAFE → Q(unclassified) → INVALID
    hits = find_kferq_motifs("QAFERQ")
    _check("QAFERQ → no valid hit",     len(hits) == 0)

    print()
    print("─" * 60)
    print("  LIR / AIM MOTIF TESTS")
    print("─" * 60)

    # p62-like: DDD upstream → xlir_acidic
    hits = find_lir_motifs("DDDWTHLSS")
    _check("DDDWTHLSS → lir_xlir_acidic (DDD upstream of W)",
           any(h["type"] == "lir_xlir_acidic" for h in hits))

    # No upstream context → canonical
    hits = find_lir_motifs("AAWTILAA")
    _check("AAWTILAA → lir_canonical (no acidic/phospho upstream)",
           any(h["type"] == "lir_canonical" for h in hits))

    # S/T upstream only → phospho
    hits = find_lir_motifs("SSWTILAA")
    _check("SSWTILAA → lir_xlir_phospho (SS upstream)",
           any(h["type"] == "lir_xlir_phospho" for h in hits))

    # Both D/E and S/T upstream → optimal
    hits = find_lir_motifs("DSWTILAA")
    _check("DSWTILAA → lir_xlir_optimal (D+S upstream)",
           any(h["type"] == "lir_xlir_optimal" for h in hits))

    # Non-LIR: K at position +3 (not I/L/V)
    hits = find_lir_motifs("WAAAKAA")
    _check("WAAAKAA → no LIR (K at +3 not in ILV)",
           len(hits) == 0)

    print()
    print("─" * 60)
    print("  D-BOX TESTS")
    print("─" * 60)

    # Minimal RxxL
    hits = find_dbox_motifs("AAARTALGAA")
    _check("RTALGAA → dbox_minimal (RxxL)",
           any(h["type"] == "dbox_minimal" and "RTAL" in h["motif"] for h in hits))

    # Extended RxxLxxxxN
    hits = find_dbox_motifs("AAARTALGSSSN")
    _check("RTALGSSSN → dbox_extended (RxxLxxxxN)",
           any(h["type"] == "dbox_extended" for h in hits))

    # K instead of R → no D-box
    hits = find_dbox_motifs("KTALGAA")
    _check("KTALGAA → no D-box (K≠R)",
           not any(h["type"].startswith("dbox") for h in hits))

    print()
    print("─" * 60)
    print("  KEN BOX TESTS")
    print("─" * 60)

    hits = find_ken_motifs("AAKENAAB")
    _check("AAKENAAB → ken_box",         any(h["type"] == "ken_box" for h in hits))
    hits = find_ken_motifs("AAKFNAAB")
    _check("AAKFNAAB → no KEN (K-F-N)",  len(hits) == 0)

    print()
    print("─" * 60)
    print("  N-DEGRON TESTS")
    print("─" * 60)

    r = classify_n_degron("RAAAGGG")
    _check("RAAAGGG → type_I_primary_destabilising",
           "type_I_primary_destabilising" in r.get("type", ""))

    r = classify_n_degron("MAAAGGG")
    _check("MAAAGGG → MAP cleavage → stabilising (M→A)",
           "stabilising" in r.get("type", "") and "MAP" in r.get("notes", ""))

    r = classify_n_degron("DAAAGGG")
    _check("DAAAGGG → secondary_destabilising",
           "secondary_destabilising" in r.get("type", ""))

    r = classify_n_degron("NAAAGGG")
    _check("NAAAGGG → tertiary_destabilising",
           "tertiary_destabilising" in r.get("type", ""))

    print()
    print("─" * 60)
    print("  C-DEGRON TESTS")
    print("─" * 60)

    hits = find_cdegron_motifs("AAAAEE")
    _check("AAAAEE → cdegron_dcaf12 (-EE)",
           any(h["type"] == "cdegron_dcaf12" for h in hits))

    hits = find_cdegron_motifs("AAAARAAGG")
    _check("AAAARAAGG → cdegron_appbp2 (RxxGx)",
           any(h["type"] == "cdegron_appbp2" for h in hits))

    hits = find_cdegron_motifs("AAAAMAAGG")
    _check("AAAAMAAGG → no C-degron (M≠R, last 2 not EE)",
           len(hits) == 0)

    print()
    print("─" * 60)
    print("  PEST SEQUENCE TESTS")
    print("─" * 60)

    # Flanked K…K, region SDEPSTSSTTPPE (13 aa): P=3, E=2, S=4, T=2 → PEST frac=11/13=0.85
    hits = find_pest_sequences("RSDEPSTSSTTPPEK")
    _check("RSDEPSTSSTTPPEK → pest_positive",
           any(h["type"] == "pest_positive" for h in hits))

    # No P → no PEST
    hits = find_pest_sequences("RSDESSTTEEESSK")
    _check("RSDESSTTEEESSK → no PEST (no P inside)",
           len(hits) == 0)

    print()
    print("─" * 60)
    print("  UNIFIED API TEST")
    print("─" * 60)

    seq = "MSQIRLKEFERQSLPKENDDDWTILRAALGSSSNEEMAAEE"
    r = find_all_degradation_motifs(seq)
    _check("Unified API returns all 8 keys",
           set(r.keys()) == {"kferq","lir","dbox","ken","n_degron","c_degron","pest","summary"})
    _check("Unified API: summary is a dict",
           isinstance(r["summary"], dict))
    _check("Unified API: KFERQ count is int",
           isinstance(r["summary"]["n_kferq"], int))

    print()
    bar = "═" * 60
    total = passed + failed
    print(bar)
    if failed == 0:
        print(f"  ALL {total} TESTS PASSED — motif_finder.py v2.0 READY")
    else:
        print(f"  {passed}/{total} passed — {failed} FAILED")
    print(bar)
    if failed:
        raise AssertionError(f"{failed} self-test(s) failed.")


# Allow --self-test flag
if __name__ == "__main__":
    if "--self-test" in sys.argv:
        _self_test()
    else:
        main()