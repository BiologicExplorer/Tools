"""
orchestrator/tests/orchestrator_test.py
=======================================
Regression tests for the SEA Orchestrator.

Verified baseline (HCV polyprotein Q9WMX2 vs CYP2E1 P05181):
  - top final_sea_score  ≥ 17.0
  - top architecture     = SUPER_EPITOPE
  - at least 1 super-epitope hit

Run from repo root:
    python orchestrator/tests/orchestrator_test.py

Or with pytest:
    pytest orchestrator/tests/orchestrator_test.py -v
"""

from __future__ import annotations

import json
import os
import sys
import unittest

# ── Path setup: repo root is two levels above this file ──────────────────────
_REPO_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from orchestrator.orchestrator import (
    OrchestratorResult,
    mclachlan_to_pairs,
    _dedup_pairs,
    orchestrate,
)

# ════════════════════════════════════════════════════════════════════════════
#  Sequence & data loading (fetched once via UniProt REST API)
# ════════════════════════════════════════════════════════════════════════════

_DATA_DIR = os.path.join(_REPO_ROOT, "protein_degradation")
_V3_JSON  = os.path.join(_DATA_DIR, "Q9WMX2_vs_P05181_scored_v3.json")

# Module-level cache so API is hit only once per test session
_VIRAL_SEQ: str | None = None
_HOST_SEQ:  str | None = None
_V3_HITS:   list | None = None


def _get_test_data():
    """Return (viral_seq, host_seq, v3_hits).  Sequences fetched from UniProt once."""
    global _VIRAL_SEQ, _HOST_SEQ, _V3_HITS

    if _VIRAL_SEQ is None:
        from protein_degradation.repeat_finder.sequence_loader import load_from_uniprot
        _, _VIRAL_SEQ = load_from_uniprot("Q9WMX2")
        _, _HOST_SEQ  = load_from_uniprot("P05181")

    if _V3_HITS is None:
        with open(_V3_JSON, "r") as fh:
            data = json.load(fh)
        _V3_HITS = data["hits"]

    return _VIRAL_SEQ, _HOST_SEQ, _V3_HITS


# ════════════════════════════════════════════════════════════════════════════
#  Unit: mclachlan_to_pairs()
# ════════════════════════════════════════════════════════════════════════════

class TestMcLachlanToPairs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.viral_seq, cls.host_seq, cls.hits = _get_test_data()

    def test_filter_by_composite(self):
        """All returned pairs must have composite_primary >= threshold."""
        pairs = mclachlan_to_pairs(self.hits, self.viral_seq, min_composite=15.0)
        self.assertGreater(len(pairs), 0, "Expected at least one pair above threshold")
        for p in pairs:
            self.assertEqual(p["source"], "mclachlan")

    def test_positions_are_zero_based(self):
        """position1 and position2 must be 0-based (in valid range for their sequences)."""
        pairs = mclachlan_to_pairs(self.hits, self.viral_seq, min_composite=15.0)
        for p in pairs:
            self.assertGreaterEqual(p["position1"], 0)
            self.assertLess(p["position1"], len(self.viral_seq),
                f"position1={p['position1']} >= viral_seq length {len(self.viral_seq)}")
            self.assertGreaterEqual(p["position2"], 0)
            self.assertLess(p["position2"], len(self.host_seq),
                f"position2={p['position2']} >= host_seq length {len(self.host_seq)}")

    def test_similarity_score_range(self):
        """similarity_score must be in [0, 1]."""
        pairs = mclachlan_to_pairs(self.hits, self.viral_seq, min_composite=15.0)
        for p in pairs:
            self.assertGreaterEqual(p["similarity_score"], 0.0)
            self.assertLessEqual(p["similarity_score"], 1.0)

    def test_max_pairs_cap(self):
        """Result count must not exceed max_pairs."""
        pairs = mclachlan_to_pairs(self.hits, self.viral_seq,
                                   min_composite=0.0, max_pairs=10)
        self.assertLessEqual(len(pairs), 10)

    def test_count_at_threshold_15(self):
        """Verified count: 134 hits at composite_primary >= 15.0."""
        pairs = mclachlan_to_pairs(self.hits, self.viral_seq,
                                   min_composite=15.0, max_pairs=9999)
        self.assertEqual(len(pairs), 134,
            f"Expected 134 pairs at threshold 15.0, got {len(pairs)}")

    def test_protein1_is_full_viral_seq(self):
        """protein1 field must equal the full viral sequence."""
        pairs = mclachlan_to_pairs(self.hits, self.viral_seq, min_composite=15.0)
        for p in pairs:
            self.assertEqual(p["protein1"], self.viral_seq)


# ════════════════════════════════════════════════════════════════════════════
#  Unit: _dedup_pairs()
# ════════════════════════════════════════════════════════════════════════════

class TestDedupPairs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Need a long viral_seq placeholder (3010 chars) to pass protein1 validation
        cls._VIRAL_PLACEHOLDER = "A" * 3010

    def _make_pair(self, pos1, pos2, score=0.5, source="test"):
        return {
            "seq1": "ABCDEF", "seq2": "ABCDEF",
            "position1": pos1, "position2": pos2,
            "protein1": self._VIRAL_PLACEHOLDER,
            "similarity_score": score,
            "rank": 0,
            "source": source,
        }

    def test_identical_positions_deduped(self):
        pairs = [
            self._make_pair(100, 50, score=0.8),
            self._make_pair(100, 50, score=0.6),  # duplicate
        ]
        result = _dedup_pairs(pairs, tolerance=4)
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result[0]["similarity_score"], 0.8)

    def test_within_tolerance_deduped(self):
        """Pairs 3 residues apart in both positions should be deduplicated."""
        pairs = [
            self._make_pair(100, 50, score=0.9),
            self._make_pair(102, 52, score=0.7),  # within tolerance=4
        ]
        result = _dedup_pairs(pairs, tolerance=4)
        self.assertEqual(len(result), 1)

    def test_outside_tolerance_kept(self):
        """Pairs ≥ tolerance residues apart should both be kept."""
        pairs = [
            self._make_pair(100, 50, score=0.9),
            self._make_pair(105, 55, score=0.7),  # pos1 diff = 5 >= tolerance
        ]
        result = _dedup_pairs(pairs, tolerance=4)
        self.assertEqual(len(result), 2)

    def test_mixed_sources_deduped(self):
        """Cross-source near-duplicate should keep higher-score pair."""
        pairs = [
            self._make_pair(200, 80, score=0.85, source="seq_aligner"),
            self._make_pair(201, 81, score=0.70, source="mclachlan"),
        ]
        result = _dedup_pairs(pairs, tolerance=4)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["source"], "seq_aligner")  # higher score wins

    def test_empty_input(self):
        self.assertEqual(_dedup_pairs([]), [])


# ════════════════════════════════════════════════════════════════════════════
#  Integration: orchestrate() — HCV/CYP2E1 regression
# ════════════════════════════════════════════════════════════════════════════

class TestOrchestratorRegression(unittest.TestCase):
    """
    Primary regression suite.  Runs orchestrate() once with both pair sources
    against the verified HCV/CYP2E1 test case.
    """

    @classmethod
    def setUpClass(cls):
        viral_seq, host_seq, hits = _get_test_data()
        cls.viral_seq = viral_seq
        cls.host_seq  = host_seq
        cls.hits      = hits

        cls.result = orchestrate(
            viral_seq               = viral_seq,
            host_seq                = host_seq,
            virus_name              = "HCV",
            viral_protein_name      = "polyprotein Q9WMX2",
            host_protein_name       = "CYP2E1 P05181",
            mclachlan_hits          = hits,
            mclachlan_min_composite = 15.0,
            mclachlan_max_pairs     = 200,
        )

    # ── Result type ──────────────────────────────────────────────────────

    def test_returns_orchestrator_result(self):
        self.assertIsInstance(self.result, OrchestratorResult)

    # ── Key regression targets (from ORCHESTRATOR_CONTEXT.md) ────────────

    def test_top_score_at_least_17(self):
        """Top final_sea_score must be ≥ 17.0 (verified baseline)."""
        top_score = self.result.sea_results[0].final_sea_score
        self.assertGreaterEqual(
            top_score, 17.0,
            f"Expected top score ≥ 17.0, got {top_score:.4f}"
        )

    def test_top_architecture_is_super_epitope(self):
        """Top hit architecture must be SUPER_EPITOPE."""
        top = self.result.sea_results[0]
        arc_name = (
            top.architecture_class.name
            if hasattr(top.architecture_class, "name")
            else str(top.architecture_class)
        )
        self.assertEqual(
            arc_name, "SUPER_EPITOPE",
            f"Expected SUPER_EPITOPE, got {arc_name}"
        )

    def test_at_least_one_super_epitope(self):
        """At least one super-epitope hit must be found."""
        n_super = sum(1 for r in self.result.sea_results if r.is_super_epitope)
        self.assertGreater(n_super, 0, "Expected at least one SUPER_EPITOPE hit")

    # ── Pair source counts ────────────────────────────────────────────────

    def test_pair_source_counts_present(self):
        psc = self.result.pair_source_counts
        for key in ("seq_aligner", "mclachlan", "merged", "after_dedup"):
            self.assertIn(key, psc, f"Missing pair_source_counts key: {key}")

    def test_mclachlan_pairs_contributed(self):
        """McLachlan source should contribute ≥ 100 pairs at threshold 15.0."""
        self.assertGreaterEqual(
            self.result.pair_source_counts["mclachlan"], 100,
            "McLachlan pairs below expected count"
        )

    def test_dedup_reduces_or_equals_merged(self):
        """after_dedup count must be ≤ merged count."""
        psc = self.result.pair_source_counts
        self.assertLessEqual(psc["after_dedup"], psc["merged"])

    # ── Motif profiles ────────────────────────────────────────────────────

    def test_viral_motif_profile_has_summary(self):
        self.assertIn("summary", self.result.viral_motif_profile)

    def test_host_motif_profile_has_kferq(self):
        """CYP2E1 has known KFERQ motifs (at least 5)."""
        n = self.result.host_motif_profile.get("summary", {}).get("n_kferq", 0)
        self.assertGreaterEqual(n, 5, f"Expected ≥5 host KFERQ motifs, got {n}")

    # ── Risk summary ──────────────────────────────────────────────────────

    def test_risk_summary_overall_risk_not_low(self):
        """Given a SUPER_EPITOPE hit, overall_risk must not be LOW."""
        risk = self.result.risk_summary.get("overall_risk")
        self.assertNotEqual(risk, "LOW", f"Unexpected risk=LOW for HCV/CYP2E1 case")

    def test_risk_summary_keys_present(self):
        expected = {
            "top_architecture", "top_final_score", "n_super_epitope",
            "n_sandwiched", "n_pairs_scored", "host_is_cma_member",
            "host_cma_category", "viral_n_kferq", "host_n_kferq", "overall_risk",
        }
        for key in expected:
            self.assertIn(key, self.result.risk_summary, f"Missing risk_summary key: {key}")

    # ── Report generation ────────────────────────────────────────────────

    def test_report_writes_file(self):
        """report() must write a non-empty .md file."""
        out = "/home/sandbox/orchestrator/tests/_regression_report.md"
        written = self.result.report(out)
        self.assertTrue(os.path.exists(written))
        size = os.path.getsize(written)
        self.assertGreater(size, 500, f"Report file too small: {size} bytes")

    def test_report_contains_virus_name(self):
        """Report must mention the virus and host protein names."""
        out = "/home/sandbox/orchestrator/tests/_regression_report.md"
        if not os.path.exists(out):
            self.result.report(out)
        with open(out) as fh:
            content = fh.read()
        self.assertIn("HCV", content)
        self.assertIn("CYP2E1", content)


# ════════════════════════════════════════════════════════════════════════════
#  Integration: McLachlan-only (seq_aligner blocked)
# ════════════════════════════════════════════════════════════════════════════

class TestOrchestratorMcLachlanOnly(unittest.TestCase):
    """Verify that Orchestrator works when mclachlan_hits is supplied but
    the sequence aligner is blocked by a very high min_identity threshold."""

    @classmethod
    def setUpClass(cls):
        viral_seq, host_seq, hits = _get_test_data()
        cls.result = orchestrate(
            viral_seq               = viral_seq,
            host_seq                = host_seq,
            virus_name              = "HCV",
            viral_protein_name      = "polyprotein Q9WMX2",
            host_protein_name       = "CYP2E1 P05181",
            mclachlan_hits          = hits,
            mclachlan_min_composite = 15.0,
            # Require near-perfect identity so seq_aligner finds nothing
            seq_aligner_min_identity = 0.999,
            seq_aligner_max_pairs    = 40,
        )

    def test_produces_results(self):
        """Should still produce SEA results from McLachlan pairs alone."""
        self.assertGreater(len(self.result.sea_results), 0)

    def test_seq_aligner_count_is_zero(self):
        self.assertEqual(self.result.pair_source_counts["seq_aligner"], 0)

    def test_top_score_acceptable(self):
        """Even with seq_aligner blocked, top score should still be ≥ 10.0."""
        top = self.result.sea_results[0].final_sea_score
        self.assertGreaterEqual(top, 10.0, f"McLachlan-only top score too low: {top:.4f}")


# ════════════════════════════════════════════════════════════════════════════
#  Integration: seq_aligner-only (no McLachlan)
# ════════════════════════════════════════════════════════════════════════════

class TestOrchestratorSeqAlignerOnly(unittest.TestCase):
    """Verify that Orchestrator works when mclachlan_hits=None."""

    @classmethod
    def setUpClass(cls):
        viral_seq, host_seq, _ = _get_test_data()
        cls.result = orchestrate(
            viral_seq          = viral_seq,
            host_seq           = host_seq,
            virus_name         = "HCV",
            viral_protein_name = "polyprotein Q9WMX2",
            host_protein_name  = "CYP2E1 P05181",
            mclachlan_hits     = None,   # no McLachlan source
        )

    def test_no_mclachlan_still_runs(self):
        """Pipeline must complete without McLachlan hits."""
        self.assertIsInstance(self.result, OrchestratorResult)

    def test_mclachlan_count_is_zero(self):
        self.assertEqual(self.result.pair_source_counts["mclachlan"], 0)

    def test_has_some_results(self):
        """seq_aligner alone should find at least 1 pair."""
        self.assertGreater(len(self.result.sea_results), 0)


# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    loader = unittest.TestLoader()
    suite  = loader.loadTestsFromModule(sys.modules[__name__])
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
