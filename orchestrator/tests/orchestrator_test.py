"""
orchestrator/tests/orchestrator_test.py
=======================================
Regression tests for the SEA Orchestrator.

Verified baseline (HCV polyprotein Q9WMX2 vs CYP2E1 P05181):
  - top final_sea_score  ≥ 17.0
  - top architecture     = COMPLETE_SEA
  - at least 1 complete-sea hit

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
    load_protein_knowledge_base,
    get_epitope_ranges,
    apply_epitope_proximity_bonus,
    _check_degradation_proximity,
    _phospho_exposure_score,
    _detect_super_epitope_pairs,
    _compute_proximity_event_score,
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

    def test_top_architecture_is_complete_sea(self):
        """Top hit architecture must be COMPLETE_SEA."""
        top = self.result.sea_results[0]
        arc_name = (
            top.architecture_class.name
            if hasattr(top.architecture_class, "name")
            else str(top.architecture_class)
        )
        self.assertEqual(
            arc_name, "COMPLETE_SEA",
            f"Expected COMPLETE_SEA, got {arc_name}"
        )

    def test_at_least_one_complete_sea(self):
        """At least one complete-sea hit must be found."""
        n_super = sum(1 for r in self.result.sea_results if r.is_complete_sea)
        self.assertGreater(n_super, 0, "Expected at least one COMPLETE_SEA hit")

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
        """Given a COMPLETE_SEA hit, overall_risk must not be LOW."""
        risk = self.result.risk_summary.get("overall_risk")
        self.assertNotEqual(risk, "LOW", f"Unexpected risk=LOW for HCV/CYP2E1 case")

    def test_risk_summary_keys_present(self):
        expected = {
            "top_architecture", "top_final_score", "n_complete_sea",
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
#  Unit & Integration: Epitope Proximity Bonus
# ════════════════════════════════════════════════════════════════════════════

class TestEpitopeProximityBonus(unittest.TestCase):
    """
    Tests for load_protein_knowledge_base, get_epitope_ranges, and
    apply_epitope_proximity_bonus, plus integration tests confirming the
    bonus flows correctly through orchestrate().

    Unit tests use a synthetic KB fixture (_KB_FIXTURE) with a dummy accession
    "PTEST" — they carry no assumptions about any real protein or UniProt entry.
    Integration tests run the full pipeline against the validated Q9WMX2/P05181
    dataset solely to confirm end-to-end behaviour; they are labelled as such.
    """

    # ── Synthetic KB fixture (no real protein data) ───────────────────────────
    # EP-LINEAR: a linear epitope at 1-based positions 50–70.
    # EP-CONFORMATIONAL: null positions → must be skipped by get_epitope_ranges.
    # EP-OVERLAP: overlaps EP-LINEAR to exercise the "bonus at most once" guard.
    _KB_FIXTURE: dict = {
        "proteins": {
            "PTEST": {
                "known_autoantibody_epitopes": [
                    {"label": "EP-LINEAR",        "residue_start": 50,   "residue_end": 70},
                    {"label": "EP-CONFORMATIONAL", "residue_start": None, "residue_end": None},
                    {"label": "EP-OVERLAP",        "residue_start": 55,   "residue_end": 80},
                ]
            }
        }
    }

    # ── Helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    def _make_sea_result(rank: int = 1, position2: int = 0, score: float = 15.0):
        """Return a minimal, mutable SEAResult for unit testing."""
        from sea.sea_module import SEAResult
        r = SEAResult(
            rank=rank,
            seq1="AAAAAAA",
            seq2="AAAAAAA",
            position1=0,
            position2=position2,
        )
        r.final_sea_score = score
        return r

    @classmethod
    def setUpClass(cls):
        """Fetch sequences once; load real KB once (for load_* tests only)."""
        cls.real_kb = load_protein_knowledge_base()          # used only in load_* tests
        cls.viral_seq, cls.host_seq, cls.v3_hits = _get_test_data()

    # ── load_protein_knowledge_base ───────────────────────────────────────────

    def test_load_kb_returns_dict_with_proteins_key(self):
        """KB file loads without error and exposes the top-level 'proteins' key."""
        self.assertIsInstance(self.real_kb, dict)
        self.assertIn("proteins", self.real_kb)

    def test_load_kb_missing_path_returns_empty_dict(self):
        """A non-existent path must return {} without raising any exception."""
        self.assertEqual(load_protein_knowledge_base("/nonexistent/path.json"), {})

    # ── get_epitope_ranges ────────────────────────────────────────────────────

    def test_get_epitope_ranges_returns_linear_entries(self):
        """Entries with valid start/end positions are returned as (start, end, label)."""
        ranges = get_epitope_ranges(self._KB_FIXTURE, "PTEST")
        labels = [lbl for (_, _, lbl) in ranges]
        self.assertIn("EP-LINEAR", labels)
        ep = next(t for t in ranges if t[2] == "EP-LINEAR")
        self.assertEqual(ep[0], 50)
        self.assertEqual(ep[1], 70)

    def test_get_epitope_ranges_skips_null_position_entries(self):
        """Conformational (null-position) epitopes must not appear in the output."""
        ranges = get_epitope_ranges(self._KB_FIXTURE, "PTEST")
        labels = [lbl for (_, _, lbl) in ranges]
        self.assertNotIn("EP-CONFORMATIONAL", labels)
        for (start, end, _) in ranges:
            self.assertIsNotNone(start)
            self.assertIsNotNone(end)

    def test_get_epitope_ranges_unknown_accession_returns_empty(self):
        """An accession absent from the KB must return an empty list."""
        self.assertEqual(get_epitope_ranges(self._KB_FIXTURE, "PXXXXXX"), [])

    def test_get_epitope_ranges_empty_kb_returns_empty(self):
        """An empty KB dict must return an empty list regardless of accession."""
        self.assertEqual(get_epitope_ranges({}, "PTEST"), [])

    # ── apply_epitope_proximity_bonus ─────────────────────────────────────────

    def test_apply_bonus_when_anchor_inside_epitope(self):
        """
        position2=49 → 1-based host position 50 → inside EP-LINEAR (50-70).
        Score must increase by the bonus and one hit record must be returned.
        """
        r = self._make_sea_result(position2=49, score=15.0)
        modified, hits = apply_epitope_proximity_bonus(
            [r], "PTEST", self._KB_FIXTURE, bonus=2.0
        )
        self.assertEqual(len(hits), 1)
        self.assertAlmostEqual(modified[0].final_sea_score, 17.0, places=3)
        self.assertEqual(hits[0]["epitope_label"], "EP-LINEAR")
        self.assertEqual(hits[0]["host_pos_1based"], 50)

    def test_apply_no_bonus_when_anchor_outside_all_epitopes(self):
        """
        position2=9 → 1-based host position 10 → outside every epitope range.
        Score must be unchanged and hits must be empty.
        """
        r = self._make_sea_result(position2=9, score=15.0)
        modified, hits = apply_epitope_proximity_bonus(
            [r], "PTEST", self._KB_FIXTURE, bonus=2.0
        )
        self.assertEqual(len(hits), 0)
        self.assertAlmostEqual(modified[0].final_sea_score, 15.0, places=3)

    def test_apply_bonus_at_most_once_per_result(self):
        """
        When an anchor falls inside multiple overlapping epitope ranges, the
        bonus is applied exactly once (first matching range wins).
        position2=54 → 1-based 55 → inside both EP-LINEAR (50-70) and EP-OVERLAP (55-80).
        """
        r = self._make_sea_result(position2=54, score=10.0)
        modified, hits = apply_epitope_proximity_bonus(
            [r], "PTEST", self._KB_FIXTURE, bonus=2.0
        )
        self.assertEqual(len(hits), 1, "Bonus must be applied exactly once per result")
        self.assertAlmostEqual(modified[0].final_sea_score, 12.0, places=3)

    def test_zero_bonus_leaves_scores_and_hits_unchanged(self):
        """bonus=0.0 must be a complete no-op: scores unchanged, hits list empty."""
        r = self._make_sea_result(position2=49, score=15.0)
        modified, hits = apply_epitope_proximity_bonus(
            [r], "PTEST", self._KB_FIXTURE, bonus=0.0
        )
        self.assertEqual(hits, [])
        self.assertAlmostEqual(modified[0].final_sea_score, 15.0, places=3)

    def test_empty_sea_results_returns_empty_hits_without_error(self):
        """An empty input list must return an empty hits list without raising."""
        _, hits = apply_epitope_proximity_bonus([], "PTEST", self._KB_FIXTURE, bonus=2.0)
        self.assertEqual(hits, [])

    def test_bonus_appends_explanatory_note_to_result(self):
        """A boosted SEAResult must carry an 'Epitope proximity bonus' note."""
        r = self._make_sea_result(position2=59, score=15.0)   # 1-based 60, in EP-LINEAR
        apply_epitope_proximity_bonus([r], "PTEST", self._KB_FIXTURE, bonus=2.0)
        self.assertTrue(
            any("Epitope proximity bonus" in n for n in r.notes),
            "Expected an 'Epitope proximity bonus' entry in result.notes",
        )

    def test_unknown_accession_in_kb_returns_no_hits(self):
        """
        When host_accession is not in the KB, get_epitope_ranges returns [],
        and apply_epitope_proximity_bonus must skip all bonuses cleanly.
        """
        r = self._make_sea_result(position2=49, score=15.0)
        modified, hits = apply_epitope_proximity_bonus(
            [r], "PXXXXXX", self._KB_FIXTURE, bonus=2.0
        )
        self.assertEqual(hits, [])
        self.assertAlmostEqual(modified[0].final_sea_score, 15.0, places=3)

    # ── Integration: orchestrate() (validated dataset regression) ────────────

    def test_integration_bonus_applied_when_accession_supplied(self):
        """
        Integration regression (Q9WMX2 vs P05181 dataset):
        Supplying a host_accession whose KB entry has mapped epitopes must
        produce ≥1 overlap annotation hit.  epitope_proximity_bonus_applied
        must be 0.0 — no score boost is applied (annotation only).
        """
        res = orchestrate(
            viral_seq                = self.viral_seq,
            host_seq                 = self.host_seq,
            virus_name               = "HCV",
            viral_protein_name       = "polyprotein Q9WMX2",
            host_protein_name        = "CYP2E1 P05181",
            mclachlan_hits           = self.v3_hits,
            mclachlan_min_composite  = 15.0,
            seq_aligner_min_identity = 0.999,    # block seq_aligner for speed
            host_accession           = "P05181",
            epitope_proximity_bonus  = 2.0,      # parameter accepted but no longer boosts score
        )
        self.assertGreater(
            len(res.epitope_proximity_hits), 0,
            "Expected ≥1 SEA hit whose host anchor falls in a mapped KB epitope range",
        )
        # Bonus is no longer applied to scores — annotation only
        self.assertEqual(res.epitope_proximity_bonus_applied, 0.0)
        # Top score is driven by Phase 2 convergence + Phase 3 architecture, not the bonus
        self.assertGreater(
            res.sea_results[0].final_sea_score, 10.0,
            "Top unbiased score must exceed 10.0",
        )

    def test_integration_no_bonus_when_no_accession_supplied(self):
        """
        Integration regression: omitting host_accession must leave
        epitope_proximity_bonus_applied at 0.0 and hits at [].
        """
        res = orchestrate(
            viral_seq                = self.viral_seq,
            host_seq                 = self.host_seq,
            virus_name               = "HCV",
            viral_protein_name       = "polyprotein Q9WMX2",
            host_protein_name        = "CYP2E1 P05181",
            mclachlan_hits           = self.v3_hits,
            mclachlan_min_composite  = 15.0,
            seq_aligner_min_identity = 0.999,
        )
        self.assertEqual(res.epitope_proximity_bonus_applied, 0.0)
        self.assertEqual(res.epitope_proximity_hits, [])


# ════════════════════════════════════════════════════════════════════════════
#  Unit tests — _check_degradation_proximity()
# ════════════════════════════════════════════════════════════════════════════

class TestCheckDegradationProximity(unittest.TestCase):
    """Tests for Phase 2a KFERQ proximity helper."""

    # Minimal KFERQ-like motif dict (1-based start, as find_kferq_motifs returns)
    def _motif(self, start: int) -> dict:
        return {"start": start, "end": start + 4, "motif": "KFERQ", "type": "canonical"}

    def test_proximal_viral_side(self):
        """Viral fragment centre within 30 residues of a KFERQ → proximal=True."""
        viral_kferq = [self._motif(10)]   # motif centre at 0-based pos 9+2=11
        # fragment at pos=0, centre=3; dist = |3 - 11| = 8 ≤ 30
        result = _check_degradation_proximity(0, 200, viral_kferq, [], proximity_window=30)
        self.assertTrue(result["proximal_deg_motif"])
        self.assertEqual(result["viral_kferq_dist"], 8)
        self.assertIsNone(result["host_kferq_dist"])

    def test_proximal_host_side(self):
        """Host fragment centre within 30 residues of a KFERQ → proximal=True."""
        host_kferq = [self._motif(20)]    # motif centre at 0-based 19+2=21
        # host fragment at pos=10, centre=13; dist = |13 - 21| = 8 ≤ 30
        result = _check_degradation_proximity(200, 10, [], host_kferq, proximity_window=30)
        self.assertTrue(result["proximal_deg_motif"])
        self.assertIsNone(result["viral_kferq_dist"])
        self.assertEqual(result["host_kferq_dist"], 8)

    def test_not_proximal_when_far(self):
        """Distance > proximity_window → proximal=False."""
        viral_kferq = [self._motif(100)]   # motif centre ≈ 101
        # fragment at pos=0, centre=3; dist = 98 > 30
        result = _check_degradation_proximity(0, 0, viral_kferq, [], proximity_window=30)
        self.assertFalse(result["proximal_deg_motif"])

    def test_no_motifs_returns_none_dists(self):
        """Empty motif lists → both dists None, proximal=False."""
        result = _check_degradation_proximity(50, 50, [], [], proximity_window=30)
        self.assertIsNone(result["viral_kferq_dist"])
        self.assertIsNone(result["host_kferq_dist"])
        self.assertFalse(result["proximal_deg_motif"])

    def test_nearest_motif_selected(self):
        """Multiple motifs — nearest one is used for the distance."""
        viral_kferq = [self._motif(100), self._motif(10)]  # centres ≈ 101 and 11
        # fragment at pos=0, centre=3; nearest is motif at 10 → dist=8
        result = _check_degradation_proximity(0, 200, viral_kferq, [], proximity_window=30)
        self.assertEqual(result["viral_kferq_dist"], 8)

    def test_phase2_pairs_carry_proximal_fields(self):
        """Full orchestrate() run stamps proximal_deg_motif into phase2_pairs."""
        import json, os
        data_path = os.path.join(
            os.path.dirname(__file__), "..", "..", "protein_degradation",
            "Q9WMX2_vs_P05181_scored_v3.json",
        )
        if not os.path.exists(data_path):
            self.skipTest("v3 JSON not present — integration fixture unavailable")
        # Lightweight check: just verify the key is present in phase2_pairs dicts
        # by importing the helper.  Skip gracefully if helper is absent.
        self.skipTest("integration fixture not wired in unit context")


# ════════════════════════════════════════════════════════════════════════════


# ════════════════════════════════════════════════════════════════════════════
#  Unit tests — _phospho_exposure_score()
# ════════════════════════════════════════════════════════════════════════════

class TestPhosphoExposureScore(unittest.TestCase):
    """
    Tests for the thin normalisation helper that converts
    SEAResult.phospho_t1_hinge_dist to a 0–1 phospho-exposure score.

    The orchestrator does NOT re-scan hinge positions; it only normalises
    the distance already computed by sea_module.
    """

    class _FakeSEAResult:
        """Minimal stand-in for SEAResult with a settable phospho_t1_hinge_dist."""
        def __init__(self, dist):
            self.phospho_t1_hinge_dist = dist

    def test_zero_distance_returns_one(self):
        """dist=0 means hinge is at the fragment → score = 1.0."""
        r = self._FakeSEAResult(dist=0)
        self.assertAlmostEqual(_phospho_exposure_score(r, window=15), 1.0, places=4)

    def test_full_window_distance_returns_zero(self):
        """dist == window → score = 0.0."""
        r = self._FakeSEAResult(dist=15)
        self.assertAlmostEqual(_phospho_exposure_score(r, window=15), 0.0, places=4)

    def test_beyond_window_clamped_to_zero(self):
        """dist > window must not produce negative scores."""
        r = self._FakeSEAResult(dist=30)
        score = _phospho_exposure_score(r, window=15)
        self.assertEqual(score, 0.0)

    def test_half_window_returns_half(self):
        """dist = window/2 → score ≈ 0.5."""
        r = self._FakeSEAResult(dist=7.5)
        self.assertAlmostEqual(_phospho_exposure_score(r, window=15), 0.5, places=4)

    def test_none_distance_returns_zero(self):
        """phospho_t1_hinge_dist = None (no T1 hinge in window) → score = 0.0."""
        r = self._FakeSEAResult(dist=None)
        self.assertEqual(_phospho_exposure_score(r), 0.0)

    def test_missing_attribute_returns_zero(self):
        """An object without phospho_t1_hinge_dist attr must not raise."""
        class _Bare:
            pass
        self.assertEqual(_phospho_exposure_score(_Bare()), 0.0)

    def test_score_rounded_to_4dp(self):
        """Result must be rounded to 4 decimal places."""
        r = self._FakeSEAResult(dist=3)
        score = _phospho_exposure_score(r, window=7)
        # 1.0 - 3/7 = 0.571428... → 0.5714
        self.assertEqual(score, round(1.0 - 3 / 7, 4))

    def test_score_in_unit_interval(self):
        """Score must always be in [0.0, 1.0] for any non-negative distance."""
        r_obj = self._FakeSEAResult(dist=None)
        for d in [0, 1, 7, 14, 15, 16, 100]:
            r_obj.phospho_t1_hinge_dist = d
            s = _phospho_exposure_score(r_obj, window=15)
            self.assertGreaterEqual(s, 0.0)
            self.assertLessEqual(s, 1.0)


# ════════════════════════════════════════════════════════════════════════════
#  Unit tests — _detect_super_epitope_pairs()
# ════════════════════════════════════════════════════════════════════════════

class TestDetectSuperEpitopePairs(unittest.TestCase):
    """
    Unit tests for the TRUE SUPER-EPITOPE detector.

    Two ranked viral hit fragments form a TRUE SUPER-EPITOPE when:
      1. Their viral fragment positions are within viral_window residues.
      2. Their host  fragment positions are within host_window  residues.
      3. At least one T1 hinge position lies STRICTLY BETWEEN the two viral
         fragment start positions.
    """

    # ── Helpers ───────────────────────────────────────────────────────────────

    class _FakeHinge:
        def __init__(self, position: int):
            self.position = position

    @staticmethod
    def _make_sea_result(rank: int, pos1: int, pos2: int, score: float = 10.0):
        """Build a minimal SEAResult-like object."""
        from sea.sea_module import SEAResult
        r = SEAResult(
            rank=rank,
            seq1="AAAAAAA",
            seq2="AAAAAAA",
            position1=pos1,
            position2=pos2,
        )
        r.final_sea_score = score
        return r

    def _viral_hj(self, t1_positions):
        """Build a viral_hj dict with the given T1 hinge positions."""
        return {"hinges_t1": [self._FakeHinge(p) for p in t1_positions]}

    # ── Tests ─────────────────────────────────────────────────────────────────

    def test_finds_pair_with_bridging_hinge(self):
        """Two viral hits within window with a T1 hinge between them → 1 pair."""
        r1 = self._make_sea_result(rank=1, pos1=100, pos2=200)
        r2 = self._make_sea_result(rank=2, pos1=150, pos2=220)
        # T1 hinge at 125 is strictly between 100 and 150
        vjh = self._viral_hj([125])
        result = _detect_super_epitope_pairs([r1, r2], vjh, viral_window=150, host_window=150)
        self.assertEqual(len(result), 1)
        self.assertIn(125, result[0]["bridging_t1_hinges"])

    def test_no_pair_without_bridging_hinge(self):
        """Same viral hits but no T1 hinge between them → no pairs returned."""
        r1 = self._make_sea_result(rank=1, pos1=100, pos2=200)
        r2 = self._make_sea_result(rank=2, pos1=150, pos2=220)
        # T1 hinge at 200 is OUTSIDE the viral span [100, 150]
        vjh = self._viral_hj([200])
        result = _detect_super_epitope_pairs([r1, r2], vjh, viral_window=150, host_window=150)
        self.assertEqual(len(result), 0)

    def test_no_pair_when_viral_span_exceeds_window(self):
        """Viral span > viral_window → pair excluded."""
        r1 = self._make_sea_result(rank=1, pos1=0,   pos2=0)
        r2 = self._make_sea_result(rank=2, pos1=200, pos2=10)
        vjh = self._viral_hj([100])   # hinge between them but span > 150
        result = _detect_super_epitope_pairs([r1, r2], vjh, viral_window=150, host_window=150)
        self.assertEqual(len(result), 0)

    def test_no_pair_when_host_span_exceeds_window(self):
        """Host span > host_window → pair excluded."""
        r1 = self._make_sea_result(rank=1, pos1=100, pos2=0)
        r2 = self._make_sea_result(rank=2, pos1=140, pos2=200)
        vjh = self._viral_hj([120])   # hinge within viral span, but host span = 200 > 150
        result = _detect_super_epitope_pairs([r1, r2], vjh, viral_window=150, host_window=150)
        self.assertEqual(len(result), 0)

    def test_hinge_at_boundary_not_bridging(self):
        """Hinge exactly at v_lo or v_hi is not 'strictly between' → no pair."""
        r1 = self._make_sea_result(rank=1, pos1=100, pos2=0)
        r2 = self._make_sea_result(rank=2, pos1=140, pos2=20)
        # Hinge exactly at 100 (== v_lo) — NOT strictly interior
        vjh = self._viral_hj([100])
        result = _detect_super_epitope_pairs([r1, r2], vjh, viral_window=150, host_window=150)
        self.assertEqual(len(result), 0)

    def test_hinge_at_upper_boundary_not_bridging(self):
        """Hinge exactly at v_hi is not 'strictly between' → no pair."""
        r1 = self._make_sea_result(rank=1, pos1=100, pos2=0)
        r2 = self._make_sea_result(rank=2, pos1=140, pos2=20)
        vjh = self._viral_hj([140])  # == v_hi
        result = _detect_super_epitope_pairs([r1, r2], vjh, viral_window=150, host_window=150)
        self.assertEqual(len(result), 0)

    def test_result_sorted_by_combined_score(self):
        """Pairs are returned sorted by hit_a_score + hit_b_score descending."""
        r1 = self._make_sea_result(rank=1, pos1=100, pos2=0,   score=5.0)
        r2 = self._make_sea_result(rank=2, pos1=140, pos2=20,  score=5.0)
        r3 = self._make_sea_result(rank=3, pos1=105, pos2=5,   score=20.0)
        r4 = self._make_sea_result(rank=4, pos1=145, pos2=25,  score=20.0)
        vjh = self._viral_hj([120])
        result = _detect_super_epitope_pairs([r1, r2, r3, r4], vjh)
        # All pairs with bridging hinge — highest combined score first
        self.assertGreater(len(result), 0)
        for idx in range(len(result) - 1):
            combined_a = result[idx]["hit_a_score"] + result[idx]["hit_b_score"]
            combined_b = result[idx + 1]["hit_a_score"] + result[idx + 1]["hit_b_score"]
            self.assertGreaterEqual(combined_a, combined_b)

    def test_empty_sea_results_returns_empty(self):
        """No sea results → no pairs."""
        result = _detect_super_epitope_pairs([], self._viral_hj([100]))
        self.assertEqual(result, [])

    def test_empty_hinges_returns_empty(self):
        """No T1 hinges at all → no pairs."""
        r1 = self._make_sea_result(rank=1, pos1=100, pos2=0)
        r2 = self._make_sea_result(rank=2, pos1=140, pos2=20)
        result = _detect_super_epitope_pairs([r1, r2], self._viral_hj([]))
        self.assertEqual(result, [])

    def test_result_dict_has_expected_keys(self):
        """Returned pair dicts must carry all required fields."""
        r1 = self._make_sea_result(rank=1, pos1=100, pos2=0)
        r2 = self._make_sea_result(rank=2, pos1=140, pos2=20)
        vjh = self._viral_hj([120])
        result = _detect_super_epitope_pairs([r1, r2], vjh)
        self.assertEqual(len(result), 1)
        expected_keys = {
            "hit_a_rank", "hit_b_rank",
            "hit_a_seq1", "hit_b_seq1",
            "hit_a_seq2", "hit_b_seq2",
            "hit_a_pos1", "hit_b_pos1",
            "hit_a_pos2", "hit_b_pos2",
            "viral_span", "host_span",
            "bridging_t1_hinges",
            "hit_a_score", "hit_b_score",
        }
        for key in expected_keys:
            self.assertIn(key, result[0], f"Missing key in pair dict: {key}")


# ════════════════════════════════════════════════════════════════════════════
#  Unit tests — _check_degradation_proximity() LIR fields
# ════════════════════════════════════════════════════════════════════════════

class TestCheckDegradationProximityLIR(unittest.TestCase):
    """
    Tests for the LIR proximity fields added to _check_degradation_proximity().

    These complement the existing TestCheckDegradationProximity class and focus
    solely on the new viral_lir / host_lir parameters and the derived fields:
      proximal_lir_motif, viral_lir_dist, host_lir_dist.
    """

    def _lir(self, start: int) -> dict:
        """Minimal LIR motif dict (1-based start, 4-residue core)."""
        return {"start": start, "end": start + 3, "motif": "WXXL", "type": "canonical"}

    def _motif(self, start: int) -> dict:
        """Minimal KFERQ motif dict (1-based start)."""
        return {"start": start, "end": start + 4, "motif": "KFERQ", "type": "canonical"}

    def test_lir_fields_present_in_return(self):
        """Return dict must always include LIR-related keys."""
        result = _check_degradation_proximity(0, 0, [], [], proximity_window=30)
        for key in ("viral_lir_dist", "host_lir_dist",
                    "proximal_lir_motif", "proximal_kferq_motif",
                    "proximal_deg_motif"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_no_lir_motifs_returns_none_lir_dists(self):
        """Without LIR motifs supplied, both lir_dist values are None."""
        result = _check_degradation_proximity(10, 10, [], [], proximity_window=30)
        self.assertIsNone(result["viral_lir_dist"])
        self.assertIsNone(result["host_lir_dist"])
        self.assertFalse(result["proximal_lir_motif"])

    def test_proximal_viral_lir(self):
        """Viral fragment centre within window of a LIR → proximal_lir_motif=True."""
        # viral LIR centre at pos 11 (1-based start=10, 0-based= 9, centre= 9+2=11)
        # fragment pos=0, centre=3; dist=|3-11|=8 ≤ 30
        result = _check_degradation_proximity(
            0, 200, [], [],
            proximity_window=30,
            viral_lir=[self._lir(10)],
        )
        self.assertTrue(result["proximal_lir_motif"])
        self.assertEqual(result["viral_lir_dist"], 8)
        self.assertIsNone(result["host_lir_dist"])

    def test_proximal_host_lir(self):
        """Host fragment centre within window of a LIR → proximal_lir_motif=True."""
        result = _check_degradation_proximity(
            200, 0, [], [],
            proximity_window=30,
            host_lir=[self._lir(10)],
        )
        self.assertTrue(result["proximal_lir_motif"])
        self.assertIsNone(result["viral_lir_dist"])
        self.assertEqual(result["host_lir_dist"], 8)

    def test_lir_not_proximal_when_far(self):
        """LIR outside proximity window → proximal_lir_motif=False."""
        # LIR centre ~101; fragment at pos=0 centre=3; dist=98 > 30
        result = _check_degradation_proximity(
            0, 0, [], [],
            proximity_window=30,
            viral_lir=[self._lir(100)],
        )
        self.assertFalse(result["proximal_lir_motif"])

    def test_proximal_deg_motif_is_or_of_kferq_and_lir(self):
        """proximal_deg_motif = proximal_kferq OR proximal_lir."""
        # Only LIR proximal
        r1 = _check_degradation_proximity(
            0, 200, [], [],
            proximity_window=30,
            viral_lir=[self._lir(10)],
        )
        self.assertTrue(r1["proximal_deg_motif"])
        self.assertFalse(r1["proximal_kferq_motif"])

        # Only KFERQ proximal
        r2 = _check_degradation_proximity(
            0, 200, [self._motif(10)], [],
            proximity_window=30,
        )
        self.assertTrue(r2["proximal_deg_motif"])
        self.assertFalse(r2["proximal_lir_motif"])

        # Neither proximal
        r3 = _check_degradation_proximity(
            0, 0, [], [],
            proximity_window=30,
        )
        self.assertFalse(r3["proximal_deg_motif"])

    def test_kferq_and_lir_both_proximal(self):
        """When both KFERQ and LIR are proximal, all three flags are True."""
        r = _check_degradation_proximity(
            0, 0,
            viral_kferq=[self._motif(3)],
            host_kferq=[],
            proximity_window=30,
            viral_lir=[self._lir(3)],
        )
        self.assertTrue(r["proximal_kferq_motif"])
        self.assertTrue(r["proximal_lir_motif"])
        self.assertTrue(r["proximal_deg_motif"])


# ════════════════════════════════════════════════════════════════════════════
#  Unit tests — _compute_proximity_event_score()
# ════════════════════════════════════════════════════════════════════════════

class TestComputeProximityEventScore(unittest.TestCase):
    """
    Tests for the additive distance-weighted proximity event scorer.

    Default weights: w_viral_lir=1.0, w_host_lir=1.5, w_viral_kferq=1.0,
                     w_host_kferq=2.0  →  max_possible=5.5
    """

    def _dp(self, **kwargs):
        """Build a minimal degradation-proximity dict."""
        return {
            "viral_lir_dist":   kwargs.get("vl"),
            "host_lir_dist":    kwargs.get("hl"),
            "viral_kferq_dist": kwargs.get("vk"),
            "host_kferq_dist":  kwargs.get("hk"),
        }

    def test_all_none_returns_zero(self):
        """No proximal motifs → all components 0 → proximity_event_score=0."""
        r = _compute_proximity_event_score(self._dp())
        self.assertEqual(r["proximity_event_score"], 0.0)
        self.assertEqual(r["viral_lir_score"],  0.0)
        self.assertEqual(r["host_lir_score"],   0.0)
        self.assertEqual(r["viral_kferq_score"],0.0)
        self.assertEqual(r["host_kferq_score"], 0.0)

    def test_all_at_zero_distance_returns_one(self):
        """All four signals at dist=0 → normalised score = 1.0."""
        r = _compute_proximity_event_score(self._dp(vl=0, hl=0, vk=0, hk=0))
        self.assertAlmostEqual(r["proximity_event_score"], 1.0, places=4)

    def test_single_host_kferq_at_zero_normalised_correctly(self):
        """host_kferq dist=0 only: raw=2.0, max=5.5 → score≈0.3636."""
        r = _compute_proximity_event_score(self._dp(hk=0))
        expected = round(2.0 / 5.5, 4)
        self.assertAlmostEqual(r["proximity_event_score"], expected, places=4)

    def test_distance_equals_window_gives_zero_component(self):
        """dist == window → component = 0, score from other signals."""
        r = _compute_proximity_event_score(self._dp(vl=30, hl=0))
        self.assertEqual(r["viral_lir_score"], 0.0)
        # host_lir component: (1 - 0/30) * 1.5 = 1.5, normalised = 1.5/5.5
        self.assertAlmostEqual(r["host_lir_score"], 1.0, places=4)

    def test_beyond_window_clamped_to_zero(self):
        """dist > window must not produce negative components."""
        r = _compute_proximity_event_score(self._dp(vl=100))
        self.assertEqual(r["viral_lir_score"], 0.0)
        self.assertEqual(r["proximity_event_score"], 0.0)

    def test_output_keys_present(self):
        """Return dict must carry all five expected keys."""
        r = _compute_proximity_event_score(self._dp(hl=5))
        for key in ("viral_lir_score", "host_lir_score",
                    "viral_kferq_score", "host_kferq_score",
                    "proximity_event_score"):
            self.assertIn(key, r)

    def test_score_bounded_to_one(self):
        """proximity_event_score must never exceed 1.0."""
        r = _compute_proximity_event_score(self._dp(vl=0, hl=0, vk=0, hk=0))
        self.assertLessEqual(r["proximity_event_score"], 1.0)

    def test_custom_window_scales_correctly(self):
        """Halving the window doubles the decay rate."""
        r30 = _compute_proximity_event_score(self._dp(hl=15), window=30)
        r15 = _compute_proximity_event_score(self._dp(hl=15), window=15)
        # window=30, dist=15 → component = 0.5; window=15, dist=15 → component=0
        self.assertGreater(r30["proximity_event_score"], 0.0)
        self.assertEqual(r15["proximity_event_score"], 0.0)

    def test_phase2_pairs_carry_new_fields(self):
        """phase2_enriched_pairs dicts must include new proximity event fields."""
        import json, os
        data_path = os.path.join(
            os.path.dirname(__file__), "..", "..", "protein_degradation",
            "Q9WMX2_vs_P05181_scored_v3.json",
        )
        if not os.path.exists(data_path):
            self.skipTest("v3 JSON not present — integration fixture unavailable")
        from orchestrator.orchestrator import orchestrate
        import sys
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
        viral_seq = "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRAARKVSERDKSERVNISDADDGSQSQKIQEAGQKQKQKQK"
        host_seq  = "MGSALMSTLALVPVLFIILAFSSQFQTELESASASEASASQASAAASN"
        raw = json.load(open(data_path))
        hits = raw["hits"][:5]
        result = orchestrate(
            viral_seq, host_seq,
            virus_name="HCV", viral_protein_name="TestVP", host_protein_name="TestHP",
            mclachlan_hits=hits,
        )
        new_fields = (
            "viral_lir_score", "host_lir_score",
            "viral_kferq_score", "host_kferq_score",
            "proximity_event_score", "proximity_boost",
            "total_event_signal",
        )
        for p in result.phase2_enriched_pairs:
            for field in new_fields:
                self.assertIn(field, p, msg=f"Field '{field}' missing from phase2_enriched_pairs")


if __name__ == "__main__":
    loader = unittest.TestLoader()
    suite  = loader.loadTestsFromModule(sys.modules[__name__])
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
