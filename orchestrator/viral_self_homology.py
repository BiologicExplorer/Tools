"""
orchestrator/viral_self_homology.py — Viral Self-Homology Fingerprint
=======================================================================
Computes a self-homology fingerprint for a protein sequence: how often
does any k-mer (or 7-mer window scored by the McLachlan 1972 matrix)
re-appear in a non-trivial position elsewhere in the same sequence?

Viral proteins accumulate self-repetition from two evolutionary forces:
  (a) motif reuse — conserved functional k-mers (e.g. transmembrane anchors,
      glycine-rich linkers) recycled throughout the polyprotein.
  (b) autoimmune tuning — sub-sequences that mimic the host's target-cell
      proteome tend to be amplified by positive selection.

This module measures both signals.  Comparing a viral fingerprint against
a host-protein fingerprint makes the differential self-repetition visible.

Metrics produced per sequence
------------------------------
kmer_exact_repeat_fraction  dict  {k: fraction of k-mers that appear ≥2×
                                   at least diagonal_exclusion positions away}
                                  k ∈ {2, 3, 4, 5, 6, 7}
self_hit_count              int   # 7-mer window pairs with McLachlan
                                   composite_primary ≥ min_composite
                                   (excluding the diagonal band)
self_hit_density            float self_hit_count / n_windows
                                   (hits per window position, 0-1 scale)
duplicate_string_index      float fraction of window start positions that
                                   are covered by ≥1 self-hit at
                                   high_composite threshold (default 28.0).
                                   Avoids saturation seen at min_composite.
max_repeat_score            float composite_primary of the single best
                                   non-trivial self-matching pair
repeat_cluster_count        int   # independent repeat families (connected
                                   components of the self-hit graph at
                                   high_composite threshold)
top_hits                    list  top-N hit dicts with positions & scores

Usage
-----
>>> from orchestrator.viral_self_homology import compute_self_homology_fingerprint
>>> from orchestrator.viral_self_homology import compare_fingerprints
>>> from orchestrator.viral_self_homology import write_fingerprint_report

>>> fp_viral = compute_self_homology_fingerprint(viral_seq, "Q9WMX2 HCV polyprotein")
>>> fp_host  = compute_self_homology_fingerprint(host_seq,  "P05181 CYP2E1")
>>> report   = compare_fingerprints([fp_viral, fp_host])
>>> write_fingerprint_report([fp_viral, fp_host], "/home/sandbox/reports/self_homology.md")
"""
from __future__ import annotations

import time
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np

# ── McLachlan 1972 matrix (1-6 scale) — copied from mclachlan_aligner.py ─────
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
_M_NP = np.array(_M_RAW, dtype=np.int32)  # (20, 20)

# Background score = mean off-diagonal entry (≈ 1.28 for McLachlan 1972)
_BG: float = float(
    sum(_M_RAW[i][j] for i in range(20) for j in range(20) if i != j)
) / (20 * 20 - 20)

# Random-protein expected composite_primary for a window of length w:
# E[score] = w * BG  →  used as a normalisation anchor
def _expected_random_score(w: int) -> float:
    return w * _BG


# ════════════════════════════════════════════════════════════════════════════
#  EXACT K-MER REPEAT PROFILE
# ════════════════════════════════════════════════════════════════════════════

def kmer_repeat_profile(
    seq:                str,
    kmin:               int = 2,
    kmax:               int = 7,
    diagonal_exclusion: int = 10,
) -> Dict[int, float]:
    """
    For each k in [kmin, kmax], compute the fraction of distinct k-mers in
    *seq* that appear at least once more at a position ≥ diagonal_exclusion
    away.

    Returns
    -------
    dict {k: repeat_fraction}  where repeat_fraction ∈ [0, 1].

    Example: if 30 of 100 distinct 7-mers appear twice or more (at non-trivial
    positions), repeat_fraction for k=7 is 0.30.
    """
    result: Dict[int, float] = {}
    L = len(seq)
    for k in range(kmin, kmax + 1):
        if L < k + diagonal_exclusion:
            result[k] = 0.0
            continue
        # Map each k-mer to its list of start positions
        pos_map: Dict[str, List[int]] = defaultdict(list)
        for i in range(L - k + 1):
            kmer = seq[i:i + k].upper()
            pos_map[kmer].append(i)

        repeated = 0
        total    = len(pos_map)
        for positions in pos_map.values():
            if len(positions) < 2:
                continue
            # At least one pair with gap >= diagonal_exclusion?
            for a in range(len(positions)):
                found = False
                for b in range(a + 1, len(positions)):
                    if positions[b] - positions[a] >= diagonal_exclusion:
                        found = True
                        break
                if found:
                    repeated += 1
                    break

        result[k] = round(repeated / total, 4) if total else 0.0
    return result


# ════════════════════════════════════════════════════════════════════════════
#  McLACHLAN SELF-SCORE MATRIX
# ════════════════════════════════════════════════════════════════════════════

def _build_index_array(seq: str) -> np.ndarray:
    """Convert sequence to integer index array (unknown residue → -1)."""
    return np.array([_IDX.get(aa.upper(), -1) for aa in seq], dtype=np.int16)


def sliding_self_scores(
    seq:                str,
    window:             int = 7,
    diagonal_exclusion: int = 10,
) -> Tuple[np.ndarray, int]:
    """
    Compute McLachlan composite_primary for ALL pairs of 7-mer windows in
    *seq* that are at least *diagonal_exclusion* positions apart.

    Uses a (n_windows × n_windows) numpy accumulation — O(w * n_windows²)
    time, O(n_windows²) space.  For Q9WMX2 (3010 aa, w=7): ~72 MB, ~0.5 s.

    Returns
    -------
    scores_upper : np.ndarray shape (n_windows, n_windows)
        Upper-triangular matrix; scores[i][j] = composite_primary of windows
        starting at i and j (1-based: positions i+1 and j+1).
        Entries within the diagonal band (|i-j| < diagonal_exclusion) are set
        to 0 so they never trigger a hit.
    n_windows    : int  length of the first dimension
    """
    idx  = _build_index_array(seq)
    L    = len(seq)
    nw   = L - window + 1           # number of window start positions
    if nw <= 0:
        return np.zeros((0, 0), dtype=np.float32), 0

    # Clamp unknown (-1) to index 0 (Alanine) — contributes background score
    idx_clamped = np.where(idx >= 0, idx, 0).astype(np.int32)

    scores = np.zeros((nw, nw), dtype=np.float32)

    for k in range(window):
        # For position k inside the window:
        #   rows[p1] = aa index at seq[p1 + k]
        #   cols[p2] = aa index at seq[p2 + k]
        rows = idx_clamped[k : k + nw]          # shape (nw,)
        cols = idx_clamped[k : k + nw]          # same array (self-comparison)
        # Outer product of McLachlan scores
        contrib = _M_NP[rows[:, np.newaxis], cols[np.newaxis, :]]  # (nw, nw)
        scores += contrib.astype(np.float32)

    # Zero out diagonal band (positions within exclusion zone are trivial)
    diag_offsets = np.abs(
        np.arange(nw)[:, np.newaxis] - np.arange(nw)[np.newaxis, :]
    )
    scores[diag_offsets < diagonal_exclusion] = 0.0

    # Zero out lower triangle so each pair is counted once
    i_lower, j_lower = np.tril_indices(nw, k=-1)
    scores[i_lower, j_lower] = 0.0

    return scores, nw


# ════════════════════════════════════════════════════════════════════════════
#  CLUSTER / COMPONENT COUNTING
# ════════════════════════════════════════════════════════════════════════════

def _count_clusters(positions: List[int], cluster_gap: int = 20) -> int:
    """
    Group sorted positions into clusters where consecutive positions are
    within cluster_gap of each other.  Returns the number of clusters.
    """
    if not positions:
        return 0
    positions = sorted(positions)
    clusters = 1
    for i in range(1, len(positions)):
        if positions[i] - positions[i - 1] > cluster_gap:
            clusters += 1
    return clusters


def _count_components(i_arr: np.ndarray, j_arr: np.ndarray) -> int:
    """
    Count connected components of the self-hit graph using Union-Find.

    Each node is a window start position; each hit pair (i, j) is an edge.
    Returns the number of distinct connected components — i.e. independent
    repeat families.  A higher count means more *independent* repeated motifs
    (not subsumed into one large network).

    This is used for repeat_cluster_count at high_composite threshold:
    low-threshold graphs are near-complete (few components), while
    high-threshold graphs are sparse (many small families).
    """
    if len(i_arr) == 0:
        return 0
    parent: Dict[int, int] = {}

    def find(x: int) -> int:
        root = x
        while parent.get(root, root) != root:
            root = parent.get(root, root)
        # path compression
        while parent.get(x, x) != root:
            nxt = parent.get(x, x)
            parent[x] = root
            x = nxt
        return root

    def union(a: int, b: int) -> None:
        a, b = find(a), find(b)
        if a != b:
            parent[a] = b

    nodes = set(i_arr.tolist()) | set(j_arr.tolist())
    for n in nodes:
        parent[n] = n
    for a, b in zip(i_arr.tolist(), j_arr.tolist()):
        union(a, b)
    return len(set(find(n) for n in nodes))


# ════════════════════════════════════════════════════════════════════════════
#  MAIN FINGERPRINT FUNCTION
# ════════════════════════════════════════════════════════════════════════════

def compute_self_homology_fingerprint(
    seq:                str,
    protein_name:       str,
    window:             int   = 7,
    min_composite:      float = 14.5,
    high_composite:     float = 28.0,
    diagonal_exclusion: int   = 10,
    cluster_gap:        int   = 20,
    top_n:              int   = 20,
    verbose:            bool  = True,
) -> Dict:
    """
    Compute the full self-homology fingerprint for *seq*.

    Parameters
    ----------
    seq               : protein sequence (single-letter amino acids)
    protein_name      : label for this protein (used in reports)
    window            : k-mer window length for McLachlan scoring (default 7)
    min_composite     : McLachlan composite_primary threshold for overall
                        self_hit_count and self_hit_density (default 14.5).
                        Same threshold used in the cross-protein pipeline.
    high_composite    : higher threshold for duplicate_string_index and
                        repeat_cluster_count (default 28.0 = ~4.0/residue).
                        Prevents saturation of coverage metrics at min_composite.
    diagonal_exclusion: minimum position gap to be counted as a non-trivial
                        self-match (default 10)
    cluster_gap       : unused directly; kept for API compatibility.
    top_n             : number of top hits to record verbatim
    verbose           : print progress lines

    Returns
    -------
    dict with keys described in the module docstring.
    """
    L  = len(seq)
    nw = L - window + 1
    t0 = time.time()

    if verbose:
        print(f"\n[SHF] {protein_name}  ({L} aa, window={window})")

    # ── 1. Exact k-mer repeat profile ────────────────────────────────────────
    if verbose:
        print(f"  [SHF] Computing k-mer repeat profile …")
    kmer_profile = kmer_repeat_profile(seq, kmin=2, kmax=7,
                                       diagonal_exclusion=diagonal_exclusion)
    if verbose:
        for k, f in kmer_profile.items():
            print(f"    k={k}  repeat_fraction={f:.4f}")

    # ── 2. McLachlan sliding-window self-scores ───────────────────────────────
    if verbose:
        print(f"  [SHF] Building {nw}×{nw} McLachlan self-score matrix …")
    t1 = time.time()
    scores, _ = sliding_self_scores(seq, window=window,
                                    diagonal_exclusion=diagonal_exclusion)
    if verbose:
        print(f"  [SHF] Score matrix done in {time.time()-t1:.2f}s")

    # ── 3. Collect hit pairs above threshold ─────────────────────────────────
    hit_mask = scores >= min_composite            # boolean (nw, nw)
    i_hits, j_hits = np.where(hit_mask)           # arrays of row, col indices

    n_hits = int(i_hits.shape[0])
    if verbose:
        print(f"  [SHF] {n_hits} self-hits at composite_primary ≥ {min_composite}")

    # ── 4. Aggregate metrics ──────────────────────────────────────────────────
    self_hit_density = round(n_hits / nw, 5) if nw > 0 else 0.0

    # Max score (across all hits at min_composite)
    max_repeat_score = float(np.max(scores)) if n_hits > 0 else 0.0

    # ── 4b. High-threshold metrics (avoid saturation at min_composite) ────────
    # duplicate_string_index: fraction of window positions covered by ≥1
    # high-quality self-hit (scores ≥ high_composite).  Using min_composite
    # saturates to 1.0 for long viral polyproteins because almost every
    # window participates in at least one above-background match.
    i_high, j_high = np.where(scores >= high_composite)
    n_high = int(i_high.shape[0])

    if n_high > 0:
        covered = np.zeros(nw, dtype=bool)
        covered[i_high] = True
        covered[j_high] = True
        duplicate_string_index = round(float(np.sum(covered)) / nw, 4)
    else:
        duplicate_string_index = 0.0

    # repeat_cluster_count: connected components of the high-threshold hit
    # graph.  Each node is a window position; each hit pair is an edge.
    # This counts independent repeat families, not just positional clusters.
    repeat_cluster_count = _count_components(i_high, j_high)

    # Top-N hits
    hit_scores = scores[i_hits, j_hits]
    if n_hits > 0:
        order = np.argsort(-hit_scores)[:top_n]
        top_hits = []
        for idx in order:
            p1 = int(i_hits[idx])   # 0-based window start
            p2 = int(j_hits[idx])
            sc = float(scores[p1, p2])
            top_hits.append({
                "p1_pos":            p1 + 1,        # 1-based
                "p2_pos":            p2 + 1,
                "p1_seq":            seq[p1 : p1 + window],
                "p2_seq":            seq[p2 : p2 + window],
                "composite_primary": round(sc, 2),
                "distance":          p2 - p1,       # always positive (upper tri)
                "is_exact_match":    seq[p1 : p1 + window] == seq[p2 : p2 + window],
            })
    else:
        top_hits = []

    # Score distribution buckets (for histogram-style reporting)
    score_vals = hit_scores if n_hits > 0 else np.array([], dtype=np.float32)
    bucket_edges = [min_composite, 18, 24, 30, 36, 42, 1e9]
    score_distribution: Dict[str, int] = {}
    for lo, hi in zip(bucket_edges[:-1], bucket_edges[1:]):
        label = f"{lo:.0f}–{hi:.0f}" if hi < 1e9 else f"≥{lo:.0f}"
        count = int(np.sum((score_vals >= lo) & (score_vals < hi)))
        score_distribution[label] = count

    dt = time.time() - t0
    if verbose:
        print(f"  [SHF] Fingerprint complete in {dt:.2f}s total")

    return {
        "protein_name":           protein_name,
        "sequence_length":        L,
        "window":                 window,
        "min_composite":          min_composite,
        "high_composite":         high_composite,
        "diagonal_exclusion":     diagonal_exclusion,
        "kmer_repeat_profile":    kmer_profile,
        "self_hit_count":         n_hits,
        "self_hit_density":       self_hit_density,
        "duplicate_string_index": duplicate_string_index,
        "max_repeat_score":       round(max_repeat_score, 2),
        "repeat_cluster_count":   repeat_cluster_count,
        "score_distribution":     score_distribution,
        "top_hits":               top_hits,
        "compute_time_s":         round(dt, 2),
    }


# ════════════════════════════════════════════════════════════════════════════
#  FINGERPRINT COMPARISON
# ════════════════════════════════════════════════════════════════════════════

def compare_fingerprints(fingerprints: List[Dict]) -> Dict:
    """
    Compare a list of fingerprints (≥ 2).  Returns a dict with:
      - per-metric table (list of rows)
      - delta rows (viral fingerprint - host fingerprint for each pair)
      - ratio rows (viral / host)
    """
    if len(fingerprints) < 2:
        raise ValueError("Need at least 2 fingerprints to compare")

    names = [fp["protein_name"] for fp in fingerprints]

    scalar_metrics = [
        "sequence_length",
        "self_hit_count",
        "self_hit_density",
        "duplicate_string_index",
        "max_repeat_score",
        "repeat_cluster_count",
        "high_composite",
    ]

    table = []
    for metric in scalar_metrics:
        row = {"metric": metric}
        for fp in fingerprints:
            row[fp["protein_name"]] = fp.get(metric, "—")
        # ratio: first / second
        v0 = fingerprints[0].get(metric, 0) or 0
        v1 = fingerprints[1].get(metric, 0) or 0
        if isinstance(v0, (int, float)) and isinstance(v1, (int, float)):
            row["ratio_[0]/[1]"] = round(v0 / v1, 3) if v1 else "∞"
        table.append(row)

    # k-mer profile comparison
    kmer_rows = []
    for k in range(2, 8):
        row = {"k": k}
        for fp in fingerprints:
            row[fp["protein_name"]] = fp["kmer_repeat_profile"].get(k, 0.0)
        v0 = fingerprints[0]["kmer_repeat_profile"].get(k, 0.0)
        v1 = fingerprints[1]["kmer_repeat_profile"].get(k, 0.0)
        row["ratio_[0]/[1]"] = round(v0 / v1, 3) if v1 else "∞"
        kmer_rows.append(row)

    return {
        "names":          names,
        "scalar_table":   table,
        "kmer_table":     kmer_rows,
    }


# ════════════════════════════════════════════════════════════════════════════
#  MARKDOWN REPORT
# ════════════════════════════════════════════════════════════════════════════

def write_fingerprint_report(
    fingerprints: List[Dict],
    output_path:  str,
) -> None:
    """
    Write a markdown self-homology fingerprint report comparing all supplied
    fingerprints.  Suitable for an Obsidian vault or GitHub render.
    """
    from datetime import datetime, timezone
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines = []
    lines.append("# Viral Self-Homology Fingerprint Report\n")
    lines.append(f"*Generated {now}*\n")
    lines.append("---\n")

    # ── Executive summary ────────────────────────────────────────────────────
    lines.append("## Executive Summary\n")
    for fp in fingerprints:
        name  = fp["protein_name"]
        L     = fp["sequence_length"]
        n_hit = fp["self_hit_count"]
        dens  = fp["self_hit_density"]
        dsi   = fp["duplicate_string_index"]
        max_s = fp["max_repeat_score"]
        k7    = fp["kmer_repeat_profile"].get(7, 0)
        lines.append(
            f"**{name}** ({L} aa): {n_hit} self-hits, density {dens:.5f}, "
            f"duplicate_string_index {dsi:.4f}, max_score {max_s}, "
            f"7-mer exact repeat fraction {k7:.4f}\n"
        )
    lines.append("\n")

    # ── Scalar metric table ───────────────────────────────────────────────────
    if len(fingerprints) >= 2:
        comp = compare_fingerprints(fingerprints)
        lines.append("## Scalar Metric Comparison\n")

        # header
        header_names = [fp["protein_name"] for fp in fingerprints]
        header = "| Metric | " + " | ".join(header_names) + " | Ratio [0]/[1] |"
        sep    = "| --- |" + " --- |" * (len(header_names) + 1)
        lines.append(header)
        lines.append(sep)
        for row in comp["scalar_table"]:
            vals = [str(row.get(n, "—")) for n in header_names]
            ratio = str(row.get("ratio_[0]/[1]", "—"))
            lines.append("| " + row["metric"] + " | " + " | ".join(vals) + " | " + ratio + " |")
        lines.append("\n")

        # k-mer table
        lines.append("## Exact K-mer Repeat Fraction\n")
        lines.append(
            "Fraction of distinct k-mers that appear ≥2× at non-trivial "
            f"positions (gap ≥ {fingerprints[0]['diagonal_exclusion']})\n"
        )
        header2 = "| k | " + " | ".join(header_names) + " | Ratio [0]/[1] |"
        sep2    = "| --- |" + " --- |" * (len(header_names) + 1)
        lines.append(header2)
        lines.append(sep2)
        for row in comp["kmer_table"]:
            vals = [f"{row.get(n, 0):.4f}" for n in header_names]
            ratio = str(row.get("ratio_[0]/[1]", "—"))
            lines.append("| " + str(row["k"]) + " | " + " | ".join(vals) + " | " + ratio + " |")
        lines.append("\n")

    # ── Score distribution ───────────────────────────────────────────────────
    lines.append("## Score Distribution (composite_primary buckets)\n")
    for fp in fingerprints:
        lines.append(f"### {fp['protein_name']}\n")
        lines.append("| Score range | # hits |")
        lines.append("| --- | --- |")
        for bucket, count in fp["score_distribution"].items():
            lines.append(f"| {bucket} | {count} |")
        lines.append("\n")

    # ── Top hits per protein ─────────────────────────────────────────────────
    for fp in fingerprints:
        lines.append(f"## Top Self-Matching Pairs — {fp['protein_name']}\n")
        if not fp["top_hits"]:
            lines.append("*No self-hits above threshold.*\n")
            continue
        lines.append("| Rank | Pos1 | Seq1 | Pos2 | Seq2 | Score | Dist | Exact? |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for i, h in enumerate(fp["top_hits"], 1):
            exact = "✓" if h["is_exact_match"] else ""
            lines.append(
                f"| {i} | {h['p1_pos']} | `{h['p1_seq']}` | {h['p2_pos']} | "
                f"`{h['p2_seq']}` | {h['composite_primary']} | {h['distance']} | {exact} |"
            )
        lines.append("\n")

    # ── Interpretation notes ──────────────────────────────────────────────────
    lines.append("## Interpretation Notes\n")
    lines.append(
        "- **self_hit_density**: hits per window position at min_composite threshold. "
        "Viral proteins with recycled functional motifs score significantly higher than host proteins.\n"
        "- **duplicate_string_index**: fraction of window positions covered by ≥1 self-hit "
        "at **high_composite** threshold (default 28.0, ~4.0/residue). "
        "This avoids saturation — at the lower min_composite almost every position is covered. "
        "Values > 0.10 indicate that >10% of the sequence participates in a strong repeated element.\n"
        "- **repeat_cluster_count**: independent repeat families (connected components of the "
        "high-composite self-hit graph). More components = more structurally diverse, "
        "independent repeat motifs. A viral protein with protein-wide dense repeats shows "
        "fewer components (positions are inter-connected) while a host protein with sparse, "
        "isolated patches shows more components.\n"
        "- **k-mer exact repeat fraction**: tracks raw string repetition. "
        "For k=7, fractions above 0.05 are unusual in non-viral eukaryotic proteins.\n"
        "- **Ratio [0]/[1]** > 2 for most metrics strongly supports the hypothesis that "
        "the viral protein carries significantly more self-repetition than the host comparator.\n"
        "- Top self-matching pairs reveal which sequence motifs are most duplicated; "
        "these are candidate molecular mimicry seeds if they also score highly in the "
        "cross-protein (viral vs host) comparison.\n"
    )

    with open(output_path, "w") as fh:
        fh.write("\n".join(lines))

    print(f"\n[SHF] Report written → {output_path}")
