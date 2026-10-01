"""
mclachlan_aligner.py
====================
Gap-tolerant Smith-Waterman local aligner using the McLachlan (1972)
amino-acid similarity matrix (1-6 scale).

Produces hit dicts compatible with orchestrator.mclachlan_to_pairs():
  {
    "motif":             str,   # viral residues in aligned region (no gap chars)
    "p1_pos1":           str,   # "start-end" 1-based in viral_seq
    "p2_window":         str,   # host residues in aligned region (no gap chars)
    "p2_pos":            str,   # "start-end" 1-based in host_seq
    "composite_primary": float, # sum of raw (1-6) McLachlan scores for matched cols
    "layer2_cross_mean": float, # composite_primary / n_matched_cols
  }

Design notes
------------
* Pure Python + stdlib only — BioPython is not required.
* Uses AFFINE gap penalties (gap_open + gap_extend per extension).
* The substitution matrix is CENTERED (background subtracted) for DP scoring
  so that mismatched residues accumulate a small negative contribution, making
  SW correctly prefer compact high-similarity regions over long mediocre ones.
* composite_primary / layer2_cross_mean are computed from RAW (uncentered)
  McLachlan scores so that the existing composite_primary >= 15.0 threshold
  in mclachlan_to_pairs() remains physically meaningful (≈ 2.14 pts/position
  average, same as the gapless v3 calibration).
"""
from __future__ import annotations

import time
from typing import Dict, List, Tuple

# ── McLachlan 1972 matrix (1-6 scale) ────────────────────────────────────────
# Identical values to generate_mclachlan_json.py so this module is self-contained.
_AAS: str = "ARNDCQEGHILKMFPSTWYV"
_IDX: Dict[str, int] = {aa: i for i, aa in enumerate(_AAS)}

_M_RAW: List[List[int]] = [
#    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
    [6, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1],  # A
    [1, 6, 1, 1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1],  # R
    [1, 1, 6, 2, 1, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1],  # N
    [1, 1, 2, 6, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],  # D
    [1, 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],  # C
    [1, 2, 2, 2, 1, 6, 3, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1],  # Q
    [1, 1, 2, 3, 1, 3, 6, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1],  # E
    [2, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1],  # G
    [1, 2, 2, 1, 1, 2, 1, 1, 6, 1, 1, 1, 1, 2, 1, 1, 1, 1, 3, 1],  # H
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 3, 1, 3, 2, 1, 1, 2, 1, 1, 4],  # I
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 6, 1, 3, 2, 1, 1, 1, 1, 1, 2],  # L
    [1, 3, 1, 1, 1, 2, 2, 1, 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1],  # K
    [1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 3, 1, 6, 2, 1, 1, 1, 1, 1, 2],  # M
    [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 6, 1, 1, 1, 3, 4, 1],  # F
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1, 1, 1],  # P
    [2, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 6, 2, 1, 1, 1],  # S
    [2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 6, 1, 1, 2],  # T
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 6, 3, 1],  # W
    [1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 4, 1, 1, 1, 3, 6, 1],  # Y
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 2, 1, 2, 1, 1, 1, 2, 1, 1, 6],  # V
]
_N20 = 20

# ── Background (average off-diagonal score) ───────────────────────────────────
_BG: float = sum(
    _M_RAW[i][j] for i in range(_N20) for j in range(_N20) if i != j
) / (_N20 * _N20 - _N20)

# ── Centered matrix for SW scoring (mismatches become slightly negative) ───────
_M_C: List[List[float]] = [
    [_M_RAW[i][j] - _BG for j in range(_N20)] for i in range(_N20)
]
_BG_C: float = 1.0 - _BG   # centered score for unknown residue

# ── Raw McLachlan lookup (used only for composite_primary computation) ─────────
def _raw(aa1: str, aa2: str) -> float:
    i = _IDX.get(aa1.upper(), -1)
    j = _IDX.get(aa2.upper(), -1)
    return float(_M_RAW[i][j]) if (i >= 0 and j >= 0) else 1.0


# ═════════════════════════════════════════════════════════════════════════════
#  SMITH-WATERMAN DP FILL  (affine gap penalties)
# ═════════════════════════════════════════════════════════════════════════════

def sw_fill(
    seq1: str,
    seq2: str,
    gap_open:   float = -4.0,
    gap_extend: float = -0.5,
):
    """
    Fill Smith-Waterman DP matrices using the centered McLachlan matrix.

    Convention:
      seq1 → rows  (viral,  length n)
      seq2 → cols  (host,   length m)
      H[i][j] = best local-alignment score ending at seq1[i-1] / seq2[j-1]
      E[i][j] = best score ending with a run of gap(s) IN seq1 at column j
                (i.e. seq2[j-1] paired with '-' from seq1)
      F[i][j] = best score ending with a run of gap(s) IN seq2 at row i
                (i.e. seq1[i-1] paired with '-' from seq2)

    Traceback codes
    ---------------
      TB_H[i][j]: 0=zero-stop  1=diagonal  2=came-from-E  3=came-from-F
      TB_E[i][j]: 0=opened-from-H[i][j-1]  1=extended-from-E[i][j-1]
      TB_F[i][j]: 0=opened-from-H[i-1][j]  1=extended-from-F[i-1][j]

    Returns (H, E, F, TB_H, TB_E, TB_F, max_H)
    """
    n = len(seq1)
    m = len(seq2)
    NEG = -1e9

    # Python list-of-lists for fast element access in tight loops
    H    = [[0.0] * (m + 1) for _ in range(n + 1)]
    E    = [[NEG]  * (m + 1) for _ in range(n + 1)]
    F    = [[NEG]  * (m + 1) for _ in range(n + 1)]
    TB_H = [[0]    * (m + 1) for _ in range(n + 1)]
    TB_E = [[0]    * (m + 1) for _ in range(n + 1)]
    TB_F = [[0]    * (m + 1) for _ in range(n + 1)]

    # Pre-index sequences
    idx1 = [_IDX.get(aa.upper(), -1) for aa in seq1]
    idx2 = [_IDX.get(aa.upper(), -1) for aa in seq2]

    # Pre-compute substitution score rows (one per viral position)
    sub: List[List[float]] = []
    for i in range(n):
        ii = idx1[i]
        if ii >= 0:
            row_c = _M_C[ii]
            sub.append([row_c[idx2[j]] if idx2[j] >= 0 else _BG_C for j in range(m)])
        else:
            sub.append([_BG_C] * m)

    max_H = 0.0

    for i in range(1, n + 1):
        Hi     = H[i];   Hi_1   = H[i - 1]
        Ei     = E[i]
        Fi     = F[i];   Fi_1   = F[i - 1]
        TBHi   = TB_H[i]
        TBEi   = TB_E[i]
        TBFi   = TB_F[i]
        sub_i  = sub[i - 1]

        for j in range(1, m + 1):
            s = sub_i[j - 1]

            # ── E: gap in seq1 (advance j, '-' in seq1) ──────────────────────
            eo = Hi[j - 1] + gap_open
            ee = Ei[j - 1] + gap_extend
            if eo >= ee:
                Ei[j] = eo;  TBEi[j] = 0
            else:
                Ei[j] = ee;  TBEi[j] = 1

            # ── F: gap in seq2 (advance i, '-' in seq2) ──────────────────────
            fo = Hi_1[j] + gap_open
            fe = Fi_1[j] + gap_extend
            if fo >= fe:
                Fi[j] = fo;  TBFi[j] = 0
            else:
                Fi[j] = fe;  TBFi[j] = 1

            # ── H ────────────────────────────────────────────────────────────
            diag  = Hi_1[j - 1] + s
            e_val = Ei[j]
            f_val = Fi[j]

            best = 0.0;  src = 0
            if diag  > best:  best = diag;   src = 1
            if e_val > best:  best = e_val;  src = 2
            if f_val > best:  best = f_val;  src = 3

            Hi[j]   = best
            TBHi[j] = src
            if best > max_H:
                max_H = best

    return H, E, F, TB_H, TB_E, TB_F, max_H


# ═════════════════════════════════════════════════════════════════════════════
#  TRACEBACK
# ═════════════════════════════════════════════════════════════════════════════

def traceback(
    seq1: str,
    seq2: str,
    TB_H,
    TB_E,
    TB_F,
    i_end: int,
    j_end: int,
) -> Tuple[str, str, int, int]:
    """
    Traceback from matrix position (i_end, j_end).

    Uses a 3-state machine (H / E / F) so that affine gap runs are
    correctly reconstructed.

    Returns (aln1, aln2, i_start, j_start) where:
      aln1, aln2  aligned strings with '-' for gaps
      i_start     matrix row  where alignment begins (0 = start of seq1)
      j_start     matrix col  where alignment begins (0 = start of seq2)

    Sequence coordinates (1-based, inclusive):
      viral : i_start+1 .. i_end
      host  : j_start+1 .. j_end
    """
    aln1: List[str] = []
    aln2: List[str] = []
    i, j = i_end, j_end
    state = "H"

    while i > 0 or j > 0:
        if state == "H":
            src = TB_H[i][j]
            if src == 0:
                break
            elif src == 1:          # diagonal
                aln1.append(seq1[i - 1])
                aln2.append(seq2[j - 1])
                i -= 1;  j -= 1
                # stay in H
            elif src == 2:          # came from E (gap in seq1)
                state = "E"
                # do NOT advance: E state will consume seq2[j-1]
            else:                   # src == 3, came from F (gap in seq2)
                state = "F"
                # do NOT advance: F state will consume seq1[i-1]

        elif state == "E":
            # gap in seq1: '-' in seq1, seq2[j-1] in seq2
            if j <= 0:
                break
            aln1.append("-")
            aln2.append(seq2[j - 1])
            src = TB_E[i][j]
            j -= 1
            state = "H" if src == 0 else "E"

        else:  # state == "F"
            # gap in seq2: seq1[i-1] in seq1, '-' in seq2
            if i <= 0:
                break
            aln1.append(seq1[i - 1])
            aln2.append("-")
            src = TB_F[i][j]
            i -= 1
            state = "H" if src == 0 else "F"

    i_start = i
    j_start = j
    aln1.reverse()
    aln2.reverse()
    return "".join(aln1), "".join(aln2), i_start, j_start


# ═════════════════════════════════════════════════════════════════════════════
#  HIT SCORING
# ═════════════════════════════════════════════════════════════════════════════

def _hit_scores(aln1: str, aln2: str) -> Tuple[float, float, int]:
    """
    Compute composite_primary and layer2_cross_mean using raw (uncentered)
    McLachlan scores.  Only matched (non-gap) columns are counted.

    Returns (composite_primary, layer2_cross_mean, n_matched).
    """
    raw_sum = 0.0
    n = 0
    for a, b in zip(aln1, aln2):
        if a == "-" or b == "-":
            continue
        raw_sum += _raw(a, b)
        n += 1
    if n == 0:
        return 0.0, 0.0, 0
    l2cm = raw_sum / n
    return raw_sum, l2cm, n


# ═════════════════════════════════════════════════════════════════════════════
#  MAIN ENTRY POINT
# ═════════════════════════════════════════════════════════════════════════════

def hits_from_sw(
    viral_seq:      str,
    host_seq:       str,
    gap_open:       float = -4.0,
    gap_extend:     float = -0.5,
    min_sw_score:   float = 5.0,
    min_composite:  float = 14.5,
    min_matched:    int   = 5,
    max_hits:             int   = 500,
    max_viral_span:       int   = 15,
    max_host_span:        int   = 15,
    max_candidate_trace:  int   = 50_000,
) -> List[Dict]:
    """
    Run Smith-Waterman alignment and return all qualifying local alignment hits.

    Parameters
    ----------
    viral_seq      : full viral protein sequence  (rows in DP)
    host_seq       : full host protein sequence   (cols in DP)
    gap_open       : affine gap-open penalty (≤ 0)
    gap_extend     : affine gap-extend penalty (≤ 0)
    min_sw_score   : minimum centered SW score for a cell to be a traceback endpoint
    min_composite  : minimum raw composite_primary (McLachlan sum) to keep a hit
    min_matched    : minimum number of non-gap aligned columns
    max_hits           : cap on returned hits
    max_viral_span     : maximum allowed span on the viral side (ve - vs + 1).
                         Hits spanning more residues are discarded.  Default 15
                         mirrors the 7-mer gapless window (≤ 7+gaps ≈ 10-15 aa).
    max_host_span      : maximum allowed span on the host side (he - hs + 1).
                         Default 15 for the same reason.
    max_candidate_trace: only traceback the top-N cells by SW score.
                         Prevents O(n·m) traceback work on full proteins.
                         Default 50 000 — captures all good short alignments.

    Returns
    -------
    List of hit dicts (sorted by composite_primary descending) ready for
    orchestrator.mclachlan_to_pairs().  Diagnostic fields prefixed with '_'.
    """
    n = len(viral_seq)
    m = len(host_seq)

    t0 = time.time()
    print(f"  SW-aligner: filling {n}×{m} DP matrix …", flush=True)
    H, E, F, TB_H, TB_E, TB_F, max_H = sw_fill(viral_seq, host_seq, gap_open, gap_extend)
    dt = time.time() - t0
    print(f"  SW-aligner: DP fill done in {dt:.1f}s, max_H = {max_H:.2f}", flush=True)

    if max_H < min_sw_score:
        print("  SW-aligner: no cells above threshold — 0 hits returned")
        return []

    # ── Collect candidate endpoints (cells ≥ min_sw_score) ───────────────────
    candidates: List[Tuple[float, int, int]] = []
    for i in range(1, n + 1):
        Hi = H[i]
        for j in range(1, m + 1):
            v = Hi[j]
            if v >= min_sw_score:
                candidates.append((v, i, j))

    # Sort descending by SW score; cap to avoid O(n·m) traceback work
    candidates.sort(key=lambda x: -x[0])
    candidates = candidates[:max_candidate_trace]
    print(f"  SW-aligner: tracing top {len(candidates):,} candidates (cap={max_candidate_trace:,})", flush=True)

    # ── Greedy non-overlapping traceback (viral-side dedup) ───────────────────
    used_viral: set = set()   # 1-based viral positions already claimed

    hits: List[Dict] = []

    for sw_score, i_end, j_end in candidates:
        if len(hits) >= max_hits:
            break

        # Quick pre-check: is i_end already claimed?
        if i_end in used_viral:
            continue

        aln1, aln2, i_start, j_start = traceback(
            viral_seq, host_seq, TB_H, TB_E, TB_F, i_end, j_end
        )
        if not aln1 or len(aln1) == 0:
            continue

        comp, l2cm, n_matched = _hit_scores(aln1, aln2)
        if comp < min_composite or n_matched < min_matched:
            continue

        # 1-based inclusive positions
        vs = i_start + 1   # viral start
        ve = i_end          # viral end
        hs = j_start + 1   # host  start
        he = j_end          # host  end

        if ve < vs or he < hs:
            continue   # malformed — skip

        # Reject mega-alignments that span more than the allowed window
        if (ve - vs + 1) > max_viral_span or (he - hs + 1) > max_host_span:
            continue

        # Check overlap with already-used viral range
        viral_range = range(vs, ve + 1)
        if any(p in used_viral for p in viral_range):
            continue

        # Mark viral positions as used
        used_viral.update(viral_range)

        viral_res = aln1.replace("-", "")
        host_res  = aln2.replace("-", "")

        if not viral_res or not host_res:
            continue

        hits.append({
            "motif":              viral_res,
            "p1_pos1":            f"{vs}-{ve}",
            "p2_window":          host_res,
            "p2_pos":             f"{hs}-{he}",
            "composite_primary":  round(comp, 4),
            "layer2_cross_mean":  round(l2cm, 4),
            # ── diagnostic (prefixed with '_', ignored by mclachlan_to_pairs) ──
            "_aln1":              aln1,
            "_aln2":              aln2,
            "_n_matched":         n_matched,
            "_sw_score":          round(sw_score, 4),
            "_has_gaps":          ("-" in aln1 or "-" in aln2),
        })

    hits.sort(key=lambda h: -h["composite_primary"])
    n_gapped = sum(1 for h in hits if h.get("_has_gaps", False))
    print(
        f"  SW-aligner: {len(hits)} hits returned "
        f"({n_gapped} gapped, {len(hits) - n_gapped} gapless)",
        flush=True,
    )
    return hits
