"""
sea_module.py
=============
Super-epitope Architecture (SEA) Scanner

Scores homologous viral/host sequence pairs for their potential to trigger
autoimmune disease via the following coordinated mechanism:

    HINGE  (SS/ST/TT)   →  host-kinase phosphorylation unfolds protein,
                            exposing buried sequences including CMA motifs
    DEGRADATION MOTIF   →  Hsc70 recognizes KFERQ-like motif on unfolded
    (KFERQ-like)           protein → routes to lysosome via LAMP-2A (CMA)
    JAMMER (GA/GR/GK)   →  resists cathepsin proteolysis, protecting the
                            peptide fragment for intact MHC-II loading
    HOMOLOGOUS SEQUENCE →  intact viral peptide on MHC mimics self-protein
                            → cross-reactive autoimmune T/B cell activation

Priority (architecture class, lowest → highest):
    NONE  <  DEGRADATION_ONLY  <  HINGE_ONLY  <  JAMMER_ONLY
          <  HINGE_AND_DEGRADATION / JAMMER_AND_DEGRADATION
          <  SANDWICHED  <  COMPLETE_SEA

Jammer scoring (action-potential amplitude model):
    Single jammer dipeptide  → small baseline contribution
    Regional jammer cluster  → amplified density score
    EBNA-style long repeat   → extreme end, capped at JAMMER_DENSITY_CAP
    This captures both sparse and saturated jammer environments on one
    continuous scale, analogous to growing action-potential amplitude.

Reference framework:
    Cuervo et al. — CMA pathway, KFERQ motif, LAMP-2A translocation
    Levitskaya et al. (1995 Nature) — EBNA-1 GA-repeat proteasome jamming
    Filimonenko et al. — phosphorylation-dependent CMA substrate exposure
"""

import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from enum import Enum


# ══════════════════════════════════════════════════════════════════════════════
#  CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

class SEAConfig:
    """All tunable parameters for the SEA scoring engine."""

    # Residue window on each side of the homologous sequence to scan
    PROXIMITY_WINDOW: int = 15

    # ── Hinge scoring ────────────────────────────────────────────────────────
    HINGE_TIER1_SCORE: float = 1.0   # pure Ser clusters: SS, SSS …
    HINGE_TIER2_SCORE: float = 0.7   # functionally equivalent: ST, TS, TT

    # ── Jammer scoring (per-match, kept for diagnostics) ─────────────────────
    JAMMER_BASE_SCORE:   float = 0.8  # per jammer hit (diagnostic reference)
    JAMMER_REPEAT_BONUS: float = 0.3  # added per repeat unit beyond the first
    JAMMER_MIN_REPEATS:  int   = 1    # minimum repeat units to report

    # ── Jammer density model (action-potential amplitude) ────────────────────
    # Each jammer dipeptide in the broader window contributes a small baseline
    # signal; nearby clusters amplify the score; EBNA-style long repeats
    # saturate it.  Per-residue distance decay prevents distant jammers from
    # dominating.
    JAMMER_DENSITY_WINDOW: int   = 30    # residues each side (2× PROXIMITY_WINDOW)
    JAMMER_DENSITY_WEIGHT: float = 0.12  # contribution per jammer dipeptide
    JAMMER_DENSITY_DECAY:  float = 0.015 # per-residue distance decay factor
    JAMMER_DENSITY_CAP:    float = 3.0   # maximum total jammer density score

    # ── Degradation motif scoring ────────────────────────────────────────────
    DEGRADATION_BASE:  float = 0.6   # per motif found within window
    DEGRADATION_DECAY: float = 0.03  # proximity decay: score × (1 - decay × dist)

    # ── Architecture bonus multipliers ───────────────────────────────────────
    # Applied multiplicatively to the sum of all component scores
    BONUS: Dict[str, float] = {
        "NONE":                   1.0,
        "DEGRADATION_ONLY":       1.5,
        "HINGE_ONLY":             1.5,
        "JAMMER_ONLY":            1.5,
        "HINGE_AND_DEGRADATION":  2.5,
        "JAMMER_AND_DEGRADATION": 2.5,
        "SANDWICHED":             3.0,
        "COMPLETE_SEA":           5.0,
    }

    # ── Cell-type / APC-tropism weights ─────────────────────────────────────
    # Viruses infecting professional APCs (DCs, macrophages, B cells) produce
    # higher-risk exposures because CMA and MHC-II loading are constitutively
    # active in those cells.  Inflammatory stress (kinase activation) is also
    # elevated in APCs during active infection, increasing hinge phosphorylation.
    APC_TROPISM: Dict[str, float] = {
        "EBV":        2.0,   # directly infects B cells (professional APCs)
        "CMV":        1.8,   # monocytes, macrophages, dendritic cells
        "HIV":        1.8,   # CD4+ T cells, macrophages
        "HHV-6":      1.6,   # T cells, macrophages
        "SARS-CoV-2": 1.5,   # macrophages/DCs implicated in long-COVID autoimmunity
        "influenza":  1.3,   # DCs, alveolar macrophages
        "HCV":        1.3,   # hepatic macrophages (Kupffer cells)
        "HSV-1":      1.2,   # DCs, some macrophage tropism
        "default":    1.0,
    }

    # ── KFERQ-like motif biochemical rules ───────────────────────────────────
    # A 5-mer qualifies if it contains Q plus all three of:
    #   positive residue (K/R), negative residue (D/E), hydrophobic (F/I/L/V)
    KFERQ_POSITIVE:    str = "KR"
    KFERQ_NEGATIVE:    str = "DE"
    KFERQ_HYDROPHOBIC: str = "FILV"


# ══════════════════════════════════════════════════════════════════════════════
#  ENUMS
# ══════════════════════════════════════════════════════════════════════════════

class HingeTier(Enum):
    TIER1 = 1   # pure Ser cluster  (SS, SSS …)
    TIER2 = 2   # mixed Ser/Thr     (ST, TS, TT)


class ArchitectureClass(Enum):
    """Ordered by predicted autoimmune risk (higher value = higher risk)."""
    NONE                    = 0
    DEGRADATION_ONLY        = 1
    HINGE_ONLY              = 2
    JAMMER_ONLY             = 3
    HINGE_AND_DEGRADATION   = 4
    JAMMER_AND_DEGRADATION  = 5
    SANDWICHED              = 6
    COMPLETE_SEA            = 7


# ══════════════════════════════════════════════════════════════════════════════
#  DATA CLASSES
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class HingeMatch:
    sequence: str
    position: int        # absolute position in full viral protein
    tier:     HingeTier
    score:    float
    side:     str        # 'N-terminal' | 'C-terminal'
    distance: int        # residues from nearest edge of homologous sequence


@dataclass
class JammerMatch:
    sequence:     str
    position:     int
    repeat_count: int
    score:        float
    side:         str
    distance:     int


@dataclass
class DegradationMotifMatch:
    motif:    str
    position: int
    protein:  str   # 'protein1' (viral) | 'protein2' (host)
    distance: int   # residues from nearest edge of homologous sequence
    source:   str   # 'inline_scan' | 'provided'


@dataclass
class SEAResult:
    # ── Input identifiers ────────────────────────────────────────────────────
    rank:      int
    seq1:      str   # viral (protein1) homologous fragment
    seq2:      str   # host  (protein2) homologous fragment
    position1: int   # start position in viral protein
    position2: int   # start position in host protein

    # ── Detected features ────────────────────────────────────────────────────
    hinges:             List[HingeMatch]            = field(default_factory=list)
    jammers:            List[JammerMatch]           = field(default_factory=list)
    degradation_motifs: List[DegradationMotifMatch] = field(default_factory=list)

    # ── Score components ─────────────────────────────────────────────────────
    base_score:          float = 0.0
    hinge_score:         float = 0.0
    jammer_score:        float = 0.0   # alias for jammer_density_score (backward compat)
    jammer_density_score: float = 0.0  # regional density model (action-potential amplitude)
    degradation_score:   float = 0.0
    architecture_bonus:  float = 1.0
    cell_type_weight:    float = 1.0
    final_sea_score:     float = 0.0

    # ── Classification ───────────────────────────────────────────────────────
    architecture_class: ArchitectureClass = ArchitectureClass.NONE
    is_sandwiched:      bool  = False
    is_complete_sea:    bool  = False

    # ── Phospho-exposure distance ─────────────────────────────────────────────
    # Distance (residues) from the viral fragment centre to the nearest T1
    # hinge found within the proximity window.  None when no T1 hinge is
    # detected.  Used by the Orchestrator to compute phospho_exposure_score.
    phospho_t1_hinge_dist: Optional[int] = None

    notes: List[str] = field(default_factory=list)


# ══════════════════════════════════════════════════════════════════════════════
#  SCANNERS
# ══════════════════════════════════════════════════════════════════════════════

class HingeScanner:
    """
    Detects phosphorylatable Ser/Thr clusters flanking a target sequence.

    Tier 1 — pure Ser runs (SS, SSS …): highest phosphorylation efficiency,
             most likely to function as molecular switches.
    Tier 2 — mixed Ser/Thr (ST, TS, TT): functionally equivalent but lower
             efficiency; treated as a weaker hinge signal.
    """

    _T1 = re.compile(r'S{2,}')                  # Tier 1: SS or longer
    _T2 = re.compile(r'(?=[ST]{2})[ST]+')        # Tier 2: any 2+ S/T stretch

    def scan(
        self,
        protein_seq:  str,
        target_start: int,
        target_end:   int,
        window:       int,
        config:       SEAConfig,
    ) -> List[HingeMatch]:

        matches: List[HingeMatch] = []
        n_region_start = max(0, target_start - window)
        c_region_end   = min(len(protein_seq), target_end + window)

        regions = [
            ('N-terminal', n_region_start,  target_start),
            ('C-terminal', target_end,       c_region_end),
        ]

        for side, r_start, r_end in regions:
            if r_start >= r_end:
                continue
            region = protein_seq[r_start:r_end]

            # Collect Tier-1 positions to avoid double-counting
            t1_positions: set = set()
            for m in self._T1.finditer(region):
                abs_pos  = r_start + m.start()
                dist     = (target_start - (abs_pos + len(m.group()))) \
                           if side == 'N-terminal' \
                           else (abs_pos - target_end)
                t1_positions.update(range(m.start(), m.end()))
                matches.append(HingeMatch(
                    sequence=m.group(),
                    position=abs_pos,
                    tier=HingeTier.TIER1,
                    score=config.HINGE_TIER1_SCORE,
                    side=side,
                    distance=max(0, dist),
                ))

            for m in self._T2.finditer(region):
                # Skip if fully covered by a Tier-1 match
                if all(i in t1_positions for i in range(m.start(), m.end())):
                    continue
                # Skip pure-Ser runs (already Tier 1)
                if set(m.group()) == {'S'}:
                    continue
                abs_pos = r_start + m.start()
                dist    = (target_start - (abs_pos + len(m.group()))) \
                          if side == 'N-terminal' \
                          else (abs_pos - target_end)
                matches.append(HingeMatch(
                    sequence=m.group(),
                    position=abs_pos,
                    tier=HingeTier.TIER2,
                    score=config.HINGE_TIER2_SCORE,
                    side=side,
                    distance=max(0, dist),
                ))

        return matches


class JammerScanner:
    """
    Detects GA/GR/GK aggregating repeats that resist lysosomal cathepsin
    proteolysis, preserving peptide fragments for MHC-II loading.

    Supported core units: GR, RG, GA, AG, GK, KG
    (sourced from EBV EBNA-1 GA-repeat evidence; open to expansion)

    Two complementary methods:
      scan()         — finds discrete repeat runs within PROXIMITY_WINDOW;
                       results stored in SEAResult.jammers for diagnostics.
      density_scan() — regional density model (action-potential amplitude):
                       scans a broader JAMMER_DENSITY_WINDOW, scores each
                       dipeptide with distance decay, sums per side.
                       Single dipeptide → baseline; cluster → amplified;
                       EBNA-style long repeat → saturated/capped.
    """

    _UNITS    = ['GR', 'RG', 'GA', 'AG', 'GK', 'KG']
    _UNIT_SET = {'GR', 'RG', 'GA', 'AG', 'GK', 'KG'}

    def scan(
        self,
        protein_seq:  str,
        target_start: int,
        target_end:   int,
        window:       int,
        config:       SEAConfig,
    ) -> List[JammerMatch]:
        """Discrete repeat-run scan within PROXIMITY_WINDOW (for diagnostics)."""

        matches: List[JammerMatch] = []
        n_start = max(0, target_start - window)
        c_end   = min(len(protein_seq), target_end + window)

        regions = [
            ('N-terminal', n_start,   target_start),
            ('C-terminal', target_end, c_end),
        ]

        for side, r_start, r_end in regions:
            if r_start >= r_end:
                continue
            region = protein_seq[r_start:r_end]
            seen_spans: set = set()

            for unit in self._UNITS:
                pattern = re.compile(f'(?:{re.escape(unit)})+')
                for m in pattern.finditer(region):
                    span = (m.start(), m.end())
                    if span in seen_spans:
                        continue
                    seen_spans.add(span)

                    repeat_count = len(m.group()) // len(unit)
                    if repeat_count < config.JAMMER_MIN_REPEATS:
                        continue

                    score   = (config.JAMMER_BASE_SCORE +
                               config.JAMMER_REPEAT_BONUS * (repeat_count - 1))
                    abs_pos = r_start + m.start()
                    dist    = (target_start - (abs_pos + len(m.group()))) \
                              if side == 'N-terminal' \
                              else (abs_pos - target_end)
                    matches.append(JammerMatch(
                        sequence=m.group(),
                        position=abs_pos,
                        repeat_count=repeat_count,
                        score=score,
                        side=side,
                        distance=max(0, dist),
                    ))

        return matches

    def density_scan(
        self,
        protein_seq:  str,
        target_start: int,
        target_end:   int,
        config:       SEAConfig,
    ) -> Tuple[float, float]:
        """
        Regional density scan over JAMMER_DENSITY_WINDOW residues on each side.

        Slides a 2-residue window through each flanking region; every position
        that matches a jammer dipeptide (GA/AG/GR/RG/GK/KG) contributes
        JAMMER_DENSITY_WEIGHT × (1 − JAMMER_DENSITY_DECAY × dist).

        Returns (n_density, c_density) — uncapped per-side scores.
        Caller applies JAMMER_DENSITY_CAP to their sum.

        Amplitude analogy:
          1 hit  → small baseline (single dipeptide micro-influence)
          ~5 hits clustered → moderate amplitude
          EBNA-style 30-residue GAGAGA… → near-saturated (approaches cap)
        """
        unit_set = self._UNIT_SET
        window   = config.JAMMER_DENSITY_WINDOW
        weight   = config.JAMMER_DENSITY_WEIGHT
        decay    = config.JAMMER_DENSITY_DECAY

        n_density = 0.0
        c_density = 0.0

        # N-terminal side
        n_start = max(0, target_start - window)
        if n_start < target_start:
            region = protein_seq[n_start:target_start]
            for i in range(len(region) - 1):
                if region[i:i+2] in unit_set:
                    abs_pos = n_start + i
                    dist    = max(0, target_start - (abs_pos + 2))
                    n_density += weight * max(0.0, 1.0 - decay * dist)

        # C-terminal side
        c_end = min(len(protein_seq), target_end + window)
        if target_end < c_end:
            region = protein_seq[target_end:c_end]
            for i in range(len(region) - 1):
                if region[i:i+2] in unit_set:
                    dist      = i  # abs_pos - target_end = (target_end+i) - target_end
                    c_density += weight * max(0.0, 1.0 - decay * dist)

        return n_density, c_density


class DegradationMotifScanner:
    """
    Detects KFERQ-like CMA targeting motifs within the flanking window.

    Inline scan: slides a 5-mer window through protein1 flanking region and
    applies the Cuervo biochemical rule (Q + positive + negative + hydrophobic).
    Overlapping hits from the same physical region are deduplicated (greedy,
    lowest distance wins within any 5-residue window) to avoid count inflation.

    Provided list: accepts pre-computed motifs from upstream tools.  Protein1
    motifs within the window are merged (no double-counting with inline hits).
    Protein2 motifs are always included — they indicate host-side CMA activity
    which raises the probability that molecular mimicry leads to autoimmunity.
    """

    def _is_kferq_like(self, pentamer: str, config: SEAConfig) -> bool:
        has_q   = 'Q' in pentamer
        has_pos = any(r in pentamer for r in config.KFERQ_POSITIVE)
        has_neg = any(r in pentamer for r in config.KFERQ_NEGATIVE)
        has_hyd = any(r in pentamer for r in config.KFERQ_HYDROPHOBIC)
        return has_q and has_pos and has_neg and has_hyd

    def scan(
        self,
        protein_seq:        str,
        target_start:       int,
        target_end:         int,
        window:             int,
        config:             SEAConfig,
        provided_motifs:    List[Dict],
        pair_rank:          int,
        host_pair_position: int,
    ) -> List[DegradationMotifMatch]:

        # ── Step 1: Inline scan of protein1 flanking region ──────────────────
        inline_hits: List[DegradationMotifMatch] = []
        scan_start = max(0, target_start - window)
        scan_end   = min(len(protein_seq), target_end + window)
        region     = protein_seq[scan_start:scan_end]

        for i in range(len(region) - 4):
            pentamer = region[i:i+5]
            if self._is_kferq_like(pentamer, config):
                abs_pos = scan_start + i
                dist    = (max(0, target_start - (abs_pos + 5))
                           if abs_pos < target_start
                           else max(0, abs_pos - target_end))
                inline_hits.append(DegradationMotifMatch(
                    motif=pentamer,
                    position=abs_pos,
                    protein='protein1',
                    distance=dist,
                    source='inline_scan',
                ))

        # ── Step 2: Deduplicate inline hits (greedy, best per 5-residue window)
        # Prevents score inflation from overlapping pentamers in the same
        # physical region (e.g., NESRIQREA yields 3 overlapping KFERQ hits).
        inline_hits.sort(key=lambda x: x.position)
        deduped: List[DegradationMotifMatch] = []
        for m in inline_hits:
            if not deduped or m.position - deduped[-1].position >= 5:
                deduped.append(m)
            elif m.distance < deduped[-1].distance:
                # Replace with the closer match (higher proximity score)
                deduped[-1] = m

        # ── Step 3: Merge with provided motifs ───────────────────────────────
        seen_p1_positions = {m.position for m in deduped}
        matches: List[DegradationMotifMatch] = list(deduped)

        for dm in provided_motifs:
            dm_protein = dm.get('protein', 'protein1')
            dm_pos     = dm.get('position', 0)
            dm_motif   = dm.get('motif', '')

            if dm_protein == 'protein1':
                if dm_pos in seen_p1_positions:
                    continue   # already covered by inline scan at same position
                dist = (max(0, target_start - (dm_pos + len(dm_motif)))
                        if dm_pos < target_start
                        else max(0, dm_pos - target_end))
                if dist <= window:
                    seen_p1_positions.add(dm_pos)
                    matches.append(DegradationMotifMatch(
                        motif=dm_motif,
                        position=dm_pos,
                        protein='protein1',
                        distance=dist,
                        source='provided',
                    ))
            else:
                # protein2: always include — indicates host CMA relevance
                dist = abs(dm_pos - host_pair_position)
                matches.append(DegradationMotifMatch(
                    motif=dm_motif,
                    position=dm_pos,
                    protein='protein2',
                    distance=dist,
                    source='provided',
                ))

        return matches


# ══════════════════════════════════════════════════════════════════════════════
#  SCORING ENGINE
# ══════════════════════════════════════════════════════════════════════════════

class SEAScorer:
    """
    Orchestrates all scanners, computes component scores, classifies
    architecture, and produces a final SEA score per homologous pair.
    """

    def __init__(self, config: Optional[SEAConfig] = None):
        self.config   = config or SEAConfig()
        self._hinges  = HingeScanner()
        self._jammers = JammerScanner()
        self._deg     = DegradationMotifScanner()

    def score(
        self,
        pair:               Dict,
        degradation_motifs: List[Dict],
        virus_name:         str = "default",
    ) -> SEAResult:

        cfg          = self.config
        seq1         = pair['seq1']
        pos1         = pair['position1']
        pos2         = pair.get('position2', 0)
        protein1     = pair['protein1']
        base_score   = pair.get('similarity_score', 0.5)
        target_start = pos1
        target_end   = pos1 + len(seq1)

        # ── Run scanners ──────────────────────────────────────────────────
        hinges = self._hinges.scan(
            protein1, target_start, target_end, cfg.PROXIMITY_WINDOW, cfg)

        jammers = self._jammers.scan(
            protein1, target_start, target_end, cfg.PROXIMITY_WINDOW, cfg)

        n_jammer_density, c_jammer_density = self._jammers.density_scan(
            protein1, target_start, target_end, cfg)

        deg_matches = self._deg.scan(
            protein_seq=protein1,
            target_start=target_start,
            target_end=target_end,
            window=cfg.PROXIMITY_WINDOW,
            config=cfg,
            provided_motifs=degradation_motifs,
            pair_rank=pair['rank'],
            host_pair_position=pos2,
        )

        # ── Component scores ──────────────────────────────────────────────
        # Hinge: take best score per side to avoid double-counting
        n_hinge_score = max((h.score for h in hinges if h.side == 'N-terminal'), default=0.0)
        c_hinge_score = max((h.score for h in hinges if h.side == 'C-terminal'), default=0.0)
        hinge_score   = n_hinge_score + c_hinge_score

        # Jammer: regional density model — sum both sides, apply cap
        jammer_density_score = min(
            n_jammer_density + c_jammer_density,
            cfg.JAMMER_DENSITY_CAP,
        )

        # Degradation: additive with proximity decay
        degradation_score = sum(
            cfg.DEGRADATION_BASE * max(0.0, 1.0 - cfg.DEGRADATION_DECAY * dm.distance)
            for dm in deg_matches
        )

        # ── Sandwiched detection ──────────────────────────────────────────
        # A feature on EACH side (hinge or jammer density) is required.
        has_n = (n_hinge_score > 0) or (n_jammer_density > 0)
        has_c = (c_hinge_score > 0) or (c_jammer_density > 0)
        is_sandwiched = has_n and has_c

        # ── Architecture classification ───────────────────────────────────
        has_hinge  = hinge_score > 0
        has_jammer = jammer_density_score > 0
        has_deg    = degradation_score > 0

        is_complete = is_sandwiched and has_deg

        if is_complete:
            arch  = ArchitectureClass.COMPLETE_SEA
            bonus = cfg.BONUS["COMPLETE_SEA"]
        elif is_sandwiched:
            arch  = ArchitectureClass.SANDWICHED
            bonus = cfg.BONUS["SANDWICHED"]
        elif has_hinge and has_deg:
            arch  = ArchitectureClass.HINGE_AND_DEGRADATION
            bonus = cfg.BONUS["HINGE_AND_DEGRADATION"]
        elif has_jammer and has_deg:
            arch  = ArchitectureClass.JAMMER_AND_DEGRADATION
            bonus = cfg.BONUS["JAMMER_AND_DEGRADATION"]
        elif has_hinge:
            arch  = ArchitectureClass.HINGE_ONLY
            bonus = cfg.BONUS["HINGE_ONLY"]
        elif has_jammer:
            arch  = ArchitectureClass.JAMMER_ONLY
            bonus = cfg.BONUS["JAMMER_ONLY"]
        elif has_deg:
            arch  = ArchitectureClass.DEGRADATION_ONLY
            bonus = cfg.BONUS["DEGRADATION_ONLY"]
        else:
            arch  = ArchitectureClass.NONE
            bonus = cfg.BONUS["NONE"]

        # ── Cell-type weight ──────────────────────────────────────────────
        cell_weight = cfg.APC_TROPISM.get(virus_name, cfg.APC_TROPISM['default'])

        # ── Final score ───────────────────────────────────────────────────
        component_sum   = base_score + hinge_score + jammer_density_score + degradation_score
        final_sea_score = round(component_sum * bonus * cell_weight, 4)

        # ── Build notes ───────────────────────────────────────────────────
        notes: List[str] = []
        if is_complete:
            notes.append("COMPLETE-SEA: sequence is sandwiched by structural "
                         "features AND flanked by CMA degradation motif(s)")
        elif is_sandwiched:
            notes.append("Sandwiched: hinge/jammer features on BOTH sides of homologous sequence")
        for h in sorted(hinges, key=lambda x: x.distance):
            notes.append(f"Hinge T{h.tier.value} '{h.sequence}' @ pos {h.position} "
                         f"[{h.side}, dist={h.distance}]")
        notes.append(
            f"Jammer density: {jammer_density_score:.4f} "
            f"(N={n_jammer_density:.4f}, C={c_jammer_density:.4f}, "
            f"cap={cfg.JAMMER_DENSITY_CAP})"
        )
        for j in sorted(jammers, key=lambda x: x.distance):
            notes.append(f"  Jammer run '{j.sequence}' ×{j.repeat_count} @ pos {j.position} "
                         f"[{j.side}, dist={j.distance}]")
        for dm in sorted(deg_matches, key=lambda x: x.distance):
            notes.append(f"Degradation motif '{dm.motif}' in {dm.protein} @ pos {dm.position} "
                         f"[dist={dm.distance}, src={dm.source}]")
        notes.append(f"Cell-type weight: {virus_name} → ×{cell_weight:.1f} (APC tropism)")

        # ── Phospho-exposure distance ─────────────────────────────────────
        # Distance from fragment centre to nearest T1 hinge (within window).
        t1_in_window = [h for h in hinges if h.tier == HingeTier.TIER1]
        phospho_t1d: Optional[int] = (
            min(h.distance for h in t1_in_window) if t1_in_window else None
        )

        return SEAResult(
            rank=pair['rank'],
            seq1=seq1,
            seq2=pair['seq2'],
            position1=pos1,
            position2=pos2,
            hinges=hinges,
            jammers=jammers,
            degradation_motifs=deg_matches,
            base_score=base_score,
            hinge_score=round(hinge_score, 4),
            jammer_score=round(jammer_density_score, 4),          # alias
            jammer_density_score=round(jammer_density_score, 4),  # primary
            degradation_score=round(degradation_score, 4),
            architecture_bonus=bonus,
            cell_type_weight=cell_weight,
            final_sea_score=final_sea_score,
            architecture_class=arch,
            is_sandwiched=is_sandwiched,
            is_complete_sea=is_complete,
            phospho_t1_hinge_dist=phospho_t1d,
            notes=notes,
        )


# ══════════════════════════════════════════════════════════════════════════════
#  TOP-LEVEL MODULE INTERFACE
# ══════════════════════════════════════════════════════════════════════════════

class SEAModule:
    """
    Main entry point for Super-epitope Architecture scanning.

    Parameters
    ----------
    virus_name : str
        Name key for APC-tropism weighting (see SEAConfig.APC_TROPISM).
    config : SEAConfig, optional
        Override default configuration.

    Inputs to run()
    ---------------
    homologous_pairs : list of dict
        Output from upstream homology tool.  Each dict must contain:
          rank             int    — input rank from upstream tool
          seq1             str    — viral protein fragment
          seq2             str    — host protein fragment
          position1        int    — start position of seq1 in protein1
          position2        int    — start position of seq2 in protein2
          protein1         str    — full viral protein sequence
          protein2         str    — full host protein sequence
          similarity_score float  — 0–1 homology score (used as base score)

    autoimmune_diseases : list of str
        Candidate diseases (currently stored in results for downstream use;
        disease-to-protein mapping is a planned extension).

    degradation_motifs : list of dict
        Pre-computed motifs from upstream tool.  Each dict:
          motif    str  — motif sequence (e.g. 'KFERQ')
          position int  — position in the referenced protein
          protein  str  — 'protein1' | 'protein2'
        If absent or empty, the module falls back to inline KFERQ-like scanning.
    """

    def __init__(
        self,
        virus_name: str = "default",
        config: Optional[SEAConfig] = None,
    ):
        self.virus_name = virus_name
        self.scorer     = SEAScorer(config or SEAConfig())

    def run(
        self,
        homologous_pairs:    List[Dict],
        autoimmune_diseases: List[str],
        degradation_motifs:  List[Dict],
    ) -> List[SEAResult]:
        """Score all pairs and return sorted results (highest SEA score first)."""
        results = [
            self.scorer.score(pair, degradation_motifs, self.virus_name)
            for pair in homologous_pairs
        ]
        results.sort(key=lambda r: r.final_sea_score, reverse=True)
        return results

    def run_from_sequences(
        self,
        viral_seq:           str,
        host_seq:            str,
        provided_motifs:     Optional[List[Dict]] = None,
        autoimmune_diseases: Optional[List[str]]  = None,
        min_identity:        float = 0.33,
        window:              int   = 12,
        step:                int   = 4,
        max_pairs:           int   = 40,
    ) -> List["SEAResult"]:
        """
        Convenience entry point: align sequences, then score.

        Calls find_homologous_pairs(viral_seq, host_seq, ...) internally,
        then forwards to run().  Requires BioPython.

        Parameters
        ----------
        viral_seq           : full viral protein sequence
        host_seq            : full host protein sequence
        provided_motifs     : pre-computed degradation motifs from motif_finder
                              (List of dicts with keys: motif, position, protein)
        autoimmune_diseases : optional list of candidate disease names
        min_identity        : passed to find_homologous_pairs (default 0.33)
        window              : sliding window length             (default 12)
        step                : stride between windows            (default 4)
        max_pairs           : cap on accepted pairs             (default 40)

        Returns
        -------
        List[SEAResult] sorted by final_sea_score descending.
        Empty list if no homologous pairs meet the identity threshold.
        """
        pairs = find_homologous_pairs(
            viral_seq    = viral_seq,
            host_seq     = host_seq,
            min_identity = min_identity,
            window       = window,
            step         = step,
            max_pairs    = max_pairs,
        )
        if not pairs:
            return []
        return self.run(
            homologous_pairs    = pairs,
            autoimmune_diseases = autoimmune_diseases or [],
            degradation_motifs  = provided_motifs or [],
        )

    def report(self, results: List[SEAResult], top_n: Optional[int] = None) -> None:
        """Print a formatted SEA scan report to stdout."""
        display = results[:top_n] if top_n else results
        W = 72

        print("\n" + "═" * W)
        print("  COMPLETE-SEA ARCHITECTURE (SEA) SCAN REPORT")
        print(f"  Virus: {self.virus_name}")
        print("═" * W)

        for i, r in enumerate(display, 1):
            if r.is_complete_sea:
                tag = "  ★★ COMPLETE-SEA ★★"
            elif r.is_sandwiched:
                tag = "  ◆ SANDWICHED"
            else:
                tag = ""

            print(f"\n{'─'*W}")
            print(f"  [{i}]  SEA Score: {r.final_sea_score:.4f}{tag}")
            print(f"  Architecture : {r.architecture_class.name}")
            print(f"  Viral seq    : {r.seq1}  (pos {r.position1})")
            print(f"  Host seq     : {r.seq2}  (pos {r.position2})")
            print(f"  Base score   : {r.base_score:.3f}  |  "
                  f"Hinge: {r.hinge_score:.3f}  |  "
                  f"Jammer density: {r.jammer_density_score:.3f}  |  "
                  f"Degradation: {r.degradation_score:.3f}")
            print(f"  Arch bonus ×{r.architecture_bonus:.1f}  |  "
                  f"Cell-type ×{r.cell_type_weight:.1f}")
            for note in r.notes:
                print(f"    → {note}")

        print(f"\n{'═'*W}")
        total  = len(results)
        supers = sum(1 for r in results if r.is_complete_sea)
        sand   = sum(1 for r in results if r.is_sandwiched and not r.is_complete_sea)
        print(f"  Pairs scanned    : {total}")
        print(f"  Complete-SEA hits: {supers}")
        print(f"  Sandwiched       : {sand}")
        print(f"  With degradation : {sum(1 for r in results if r.degradation_score > 0)}")
        print("═" * W + "\n")


# ══════════════════════════════════════════════════════════════════════════════
#  HOMOLOGOUS PAIR FINDER
# ══════════════════════════════════════════════════════════════════════════════

try:
    from Bio.Align import PairwiseAligner as _PairwiseAligner
    _BIOPYTHON_AVAILABLE = True
except ImportError:  # pragma: no cover
    _BIOPYTHON_AVAILABLE = False


def find_homologous_pairs(
    viral_seq:    str,
    host_seq:     str,
    min_identity: float = 0.33,
    window:       int   = 12,
    step:         int   = 4,
    max_pairs:    int   = 40,
) -> List[Dict]:
    """
    Find short homologous sequence pairs between a viral and host protein.

    Slides a window of *window* residues over *host_seq* (step *step*).
    Each window is locally aligned against the full *viral_seq* using
    BioPython PairwiseAligner (match=2, mismatch=-1, gap open=-5,
    gap extend=-0.5).  A pair is kept when:

      - the aligned block is >= 6 residues
      - ungapped identity >= *min_identity*
      - no previously accepted viral hit overlaps by > 50%

    Pairs are sorted by similarity descending; ranks are assigned 1-based.

    Parameters
    ----------
    viral_seq    : full viral protein sequence (query)
    host_seq     : full host protein sequence  (reference window source)
    min_identity : minimum fraction of identical residues (default 0.33)
    window       : sliding window length in residues        (default 12)
    step         : stride between windows                   (default 4)
    max_pairs    : cap on accepted pairs before sorting     (default 40)

    Returns
    -------
    List of pair dicts, each containing:
        seq1             str    viral protein fragment
        seq2             str    host protein fragment
        position1        int    0-based start in viral_seq
        position2        int    0-based start in host_seq
        protein1         str    full viral_seq (reference for SEAScorer)
        similarity_score float  ungapped identity 0–1
        rank             int    1-based rank by similarity_score descending

    Raises
    ------
    ImportError if BioPython is not installed.
    """
    if not _BIOPYTHON_AVAILABLE:
        raise ImportError(
            "BioPython is required for find_homologous_pairs(). "
            "Install it with: pip install biopython"
        )

    aligner = _PairwiseAligner()
    aligner.mode             = 'local'
    aligner.match_score      =  2
    aligner.mismatch_score   = -1
    aligner.open_gap_score   = -5
    aligner.extend_gap_score = -0.5

    pairs       = []
    seen_ranges = []   # (viral_start, viral_end) blocks already accepted

    n_windows = (len(host_seq) - window) // step + 1

    for step_i in range(n_windows):
        host_start = step_i * step
        host_win   = host_seq[host_start : host_start + window]
        if len(host_win) < window:
            break

        try:
            alignment = next(iter(aligner.align(viral_seq, host_win)))
        except StopIteration:
            continue

        blocks = alignment.aligned
        if blocks is None or len(blocks) < 2 or len(blocks[0]) == 0:
            continue

        v_s, v_e = blocks[0][0]
        h_s, h_e = blocks[1][0]

        block_len = min(v_e - v_s, h_e - h_s)
        if block_len < 6:
            continue

        viral_frag = viral_seq[v_s : v_s + block_len]
        host_frag  = host_win[h_s : h_s + block_len]

        matches  = sum(a == b for a, b in zip(viral_frag, host_frag))
        identity = matches / block_len

        if identity < min_identity:
            continue

        overlapping = any(
            max(0, min(v_s + block_len, se) - max(v_s, ss)) > block_len * 0.5
            for ss, se in seen_ranges
        )
        if overlapping:
            continue

        seen_ranges.append((v_s, v_s + block_len))
        pairs.append({
            'seq1':             viral_frag,
            'seq2':             host_frag,
            'position1':        v_s,
            'position2':        host_start + h_s,
            'protein1':         viral_seq,
            'similarity_score': round(identity, 4),
            'rank':             0,
        })

        if len(pairs) >= max_pairs:
            break

    pairs.sort(key=lambda p: p['similarity_score'], reverse=True)
    for i, p in enumerate(pairs):
        p['rank'] = i + 1

    return pairs
