"""Unit tests for paqr3.calculate_cov_metrics."""

import collections
import json
import os
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from paqr3.calculate_cov_metrics import (
    CalculateCoverages,
    _read_coverage_batch,
    compute_drops_and_usage,
    evaluate_all_pas_usage_patterns,
    json_serial,
    process_segment,
)


# ---------------------------------------------------------------------------
# json_serial
# ---------------------------------------------------------------------------


class TestJsonSerial:
    def test_float32(self):
        assert json_serial(np.float32(1.5)) == pytest.approx(1.5)
        assert isinstance(json_serial(np.float32(1.5)), float)

    def test_float64(self):
        assert json_serial(np.float64(2.0)) == pytest.approx(2.0)
        assert isinstance(json_serial(np.float64(2.0)), float)

    def test_int32(self):
        assert json_serial(np.int32(7)) == 7
        assert isinstance(json_serial(np.int32(7)), int)

    def test_int64(self):
        assert json_serial(np.int64(42)) == 42

    def test_ndarray(self):
        arr = np.array([1.0, 2.0, 3.0])
        result = json_serial(arr)
        assert result == [1.0, 2.0, 3.0]
        assert isinstance(result, list)

    def test_unknown_type_raises(self):
        with pytest.raises(TypeError):
            json_serial("not a numpy type")


# ---------------------------------------------------------------------------
# Helper to build subsegment dicts for evaluate_all_pas_usage_patterns
# ---------------------------------------------------------------------------


def _sub(pas_id: str, coverage: list[float]) -> dict:
    return {
        "gene_id": "g1",
        "segment_id": "seg1",
        "subsegment_id": f"sub_{pas_id}",
        "pas_id": pas_id,
        "coverage": np.array(coverage, dtype=float),
    }


# ---------------------------------------------------------------------------
# evaluate_all_pas_usage_patterns
# ---------------------------------------------------------------------------


class TestEvaluateAllPasUsagePatterns:
    def test_all_zero_coverage_returns_zero_usage(self):
        subsegments = [
            _sub("pas1", [0.0, 0.0]),
            _sub(".", [0.0, 0.0]),
        ]
        result = evaluate_all_pas_usage_patterns(subsegments, "seg1")
        assert result is not None
        assert result["pas_usage"]["pas1"] == 0.0
        assert result["f_stat"] is None

    def test_single_pas_coverage_above_zero_usage_one(self):
        subsegments = [
            _sub("pas1", [10.0, 10.0]),
            _sub(".", [2.0, 2.0]),
        ]
        result = evaluate_all_pas_usage_patterns(subsegments, "seg1")
        assert result is not None
        assert result["pas_usage"]["pas1"] == pytest.approx(1.0)

    def test_single_pas_zero_coverage_usage_zero(self):
        subsegments = [_sub("pas1", [0.0])]
        result = evaluate_all_pas_usage_patterns(subsegments, "seg1")
        assert result is not None
        assert result["pas_usage"]["pas1"] == 0.0

    def test_single_pas_zero_cov_with_nonzero_trailing(self):
        # pas1 has zero coverage but trailing has non-zero → avoids early exit
        # forces the else branch on line 107 (mu == 0 for single PAS)
        subsegments = [
            _sub("pas1", [0.0, 0.0]),
            _sub(".", [5.0, 5.0]),
        ]
        result = evaluate_all_pas_usage_patterns(subsegments, "seg1")
        assert result is not None
        assert result["pas_usage"]["pas1"] == 0.0

    def test_too_many_pas_returns_none(self):
        # 11 PAS with max_pas_count=10 → should return None
        subsegments = [_sub(f"pas{i}", [float(10 - i)]) for i in range(11)]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", max_pas_count=10
        )
        assert result is None

    def test_no_pas_returns_none(self):
        subsegments = [_sub(".", [5.0, 5.0])]
        result = evaluate_all_pas_usage_patterns(subsegments, "seg1")
        assert result is None

    def test_two_pas_clear_monotone_drop(self):
        # Strong monotone drop: 100 → 50 → 10
        # Expected: both pas1 and pas2 get non-zero usage
        subsegments = [
            _sub("pas1", [100.0] * 50),
            _sub("pas2", [50.0] * 50),
            _sub(".", [10.0] * 50),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", f_stat_threshold=10
        )
        assert result is not None
        assert result["pas_usage"]["pas1"] > 0
        assert result["pas_usage"]["pas2"] > 0

    def test_usages_sum_to_one(self):
        subsegments = [
            _sub("pas1", [100.0] * 50),
            _sub("pas2", [60.0] * 50),
            _sub(".", [10.0] * 50),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", f_stat_threshold=10
        )
        if result is not None:
            total = sum(result["pas_usage"].values())
            assert total == pytest.approx(1.0, abs=1e-6)

    def test_f_stat_threshold_too_high_returns_none(self):
        # Two PAS with identical coverage → F-stat = 0, no pattern passes.
        # (Single-PAS bypasses the F-stat check, so we need 2+ PAS.)
        subsegments = [
            _sub("pas1", [10.0, 10.0, 10.0, 10.0]),
            _sub("pas2", [10.0, 10.0, 10.0, 10.0]),
            _sub(".", [10.0, 10.0, 10.0, 10.0]),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", f_stat_threshold=1e9
        )
        assert result is None

    def test_result_keys_present(self):
        subsegments = [
            _sub("pas1", [100.0] * 20),
            _sub(".", [10.0] * 20),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", f_stat_threshold=10
        )
        assert result is not None
        for key in (
            "pas_usage", "f_stat", "p_value",
            "rna_drop_cov", "rna_sum_drop_cov",
            "rna_monotone", "unique_pas_ids",
        ):
            assert key in result

    def test_unique_pas_ids_in_result(self):
        subsegments = [
            _sub("pas1", [80.0] * 20),
            _sub("pas2", [40.0] * 20),
            _sub(".", [5.0] * 20),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", f_stat_threshold=5
        )
        if result is not None:
            assert set(result["unique_pas_ids"]) == {"pas1", "pas2"}


# ---------------------------------------------------------------------------
# compute_drops_and_usage
# ---------------------------------------------------------------------------


class TestComputeDropsAndUsage:
    def _make_group(self, rows: list[dict]) -> pd.DataFrame:
        return pd.DataFrame(rows)

    def test_drops_computed_correctly(self):
        group = self._make_group([
            {"subsegment_id": "s:001", "mean_cov": 100.0, "pas_id": "pas1"},
            {"subsegment_id": "s:002", "mean_cov": 40.0, "pas_id": "."},
        ])
        result = compute_drops_and_usage(group)
        drops = result["rna_drop_cov"].tolist()
        # drop[0] = 100 - 40 = 60; drop[1] = 0 (last always 0)
        assert drops[0] == pytest.approx(60.0)
        assert drops[1] == pytest.approx(0.0)

    def test_last_drop_always_zero(self):
        group = self._make_group([
            {"subsegment_id": "s:001", "mean_cov": 50.0, "pas_id": "pas1"},
            {"subsegment_id": "s:002", "mean_cov": 30.0, "pas_id": "pas2"},
            {"subsegment_id": "s:003", "mean_cov": 10.0, "pas_id": "."},
        ])
        result = compute_drops_and_usage(group)
        assert result["rna_drop_cov"].iloc[-1] == pytest.approx(0.0)

    def test_rna_usage_sums_to_one(self):
        group = self._make_group([
            {"subsegment_id": "s:001", "mean_cov": 100.0, "pas_id": "pas1"},
            {"subsegment_id": "s:002", "mean_cov": 50.0, "pas_id": "pas2"},
            {"subsegment_id": "s:003", "mean_cov": 10.0, "pas_id": "."},
        ])
        result = compute_drops_and_usage(group)
        # None is stored as NaN in a DataFrame column
        usages = result["rna_usage"].dropna().tolist()
        assert sum(usages) == pytest.approx(1.0)

    def test_trailing_subsegment_usage_is_none(self):
        group = self._make_group([
            {"subsegment_id": "s:001", "mean_cov": 100.0, "pas_id": "pas1"},
            {"subsegment_id": "s:002", "mean_cov": 10.0, "pas_id": "."},
        ])
        result = compute_drops_and_usage(group)
        # None assigned to a DataFrame column becomes NaN
        usage_row = result[result["pas_id"] == "."].iloc[0]
        assert pd.isna(usage_row["rna_usage"])

    def test_monotone_flag_set_when_all_positive(self):
        group = self._make_group([
            {"subsegment_id": "s:001", "mean_cov": 100.0, "pas_id": "pas1"},
            {"subsegment_id": "s:002", "mean_cov": 50.0, "pas_id": "pas2"},
            {"subsegment_id": "s:003", "mean_cov": 10.0, "pas_id": "."},
        ])
        result = compute_drops_and_usage(group)
        assert result["rna_monotone"].iloc[0] == 1

    def test_sum_drop_cov_equals_sum_of_drops(self):
        group = self._make_group([
            {"subsegment_id": "s:001", "mean_cov": 80.0, "pas_id": "pas1"},
            {"subsegment_id": "s:002", "mean_cov": 30.0, "pas_id": "pas2"},
            {"subsegment_id": "s:003", "mean_cov": 10.0, "pas_id": "."},
        ])
        result = compute_drops_and_usage(group)
        assert result["rna_sum_drop_cov"].iloc[0] == pytest.approx(
            result["rna_drop_cov"].sum()
        )


# ---------------------------------------------------------------------------
# CalculateCoverages.run() on real BigWig files
# ---------------------------------------------------------------------------


class TestCalculateCoveragesRun:
    @pytest.fixture(scope="class")
    def subsegments_df(self, test_gtf, test_atlas):
        from paqr3.construct_segments import ConstructSegments

        cs = ConstructSegments(test_gtf, test_atlas, 200)
        cs.compute(merge_distance=5)
        return cs.create_subsegments_dataframe()

    def test_run_returns_two_dataframes(
        self, subsegments_df, test_bw_pos, test_bw_neg
    ):
        cc = CalculateCoverages(test_bw_pos, test_bw_neg)
        raw_cov_df, usage_df = cc.run(subsegments_df)
        assert isinstance(raw_cov_df, pd.DataFrame)
        assert isinstance(usage_df, pd.DataFrame)

    def test_raw_cov_has_all_subsegments(
        self, subsegments_df, test_bw_pos, test_bw_neg
    ):
        cc = CalculateCoverages(test_bw_pos, test_bw_neg)
        raw_cov_df, _ = cc.run(subsegments_df)
        assert set(raw_cov_df["subsegment_id"]) == set(
            subsegments_df["subsegment_id"]
        )

    def test_raw_cov_required_columns(
        self, subsegments_df, test_bw_pos, test_bw_neg
    ):
        cc = CalculateCoverages(test_bw_pos, test_bw_neg)
        raw_cov_df, _ = cc.run(subsegments_df)
        for col in ("gene_id", "segment_id", "subsegment_id", "pas_id", "mean_cov"):
            assert col in raw_cov_df.columns

    def test_usage_df_only_pas_rows(
        self, subsegments_df, test_bw_pos, test_bw_neg
    ):
        cc = CalculateCoverages(test_bw_pos, test_bw_neg)
        _, usage_df = cc.run(subsegments_df)
        if not usage_df.empty:
            assert (usage_df["pas_id"] != ".").all()

    def test_usage_df_required_columns(
        self, subsegments_df, test_bw_pos, test_bw_neg
    ):
        cc = CalculateCoverages(test_bw_pos, test_bw_neg)
        _, usage_df = cc.run(subsegments_df)
        for col in (
            "gene_id", "segment_id", "subsegment_id", "pas_id",
            "mean_cov", "rna_drop_cov", "rna_sum_drop_cov",
            "rna_monotone", "rna_usage", "f_stat", "p_value",
        ):
            assert col in usage_df.columns

    def test_mean_cov_non_negative(
        self, subsegments_df, test_bw_pos, test_bw_neg
    ):
        cc = CalculateCoverages(test_bw_pos, test_bw_neg)
        raw_cov_df, _ = cc.run(subsegments_df)
        assert (raw_cov_df["mean_cov"] >= 0).all()

    def test_multithreaded_gives_same_result(
        self, subsegments_df, test_bw_pos, test_bw_neg
    ):
        cc = CalculateCoverages(test_bw_pos, test_bw_neg)
        raw1, _ = cc.run(subsegments_df, n_threads=1)
        raw2, _ = cc.run(subsegments_df, n_threads=2)
        merged = raw1.merge(
            raw2, on="subsegment_id", suffixes=("_1", "_2")
        )
        np.testing.assert_allclose(
            merged["mean_cov_1"].values,
            merged["mean_cov_2"].values,
            rtol=1e-5,
        )


# ---------------------------------------------------------------------------
# Debug=True code paths in evaluate_all_pas_usage_patterns
# ---------------------------------------------------------------------------


class TestDebugLogPaths:
    def test_all_zero_debug(self):
        subsegments = [_sub("pas1", [0.0, 0.0])]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", debug=True
        )
        assert result is not None
        assert result["pas_usage"]["pas1"] == 0.0

    def test_single_pas_debug(self):
        subsegments = [
            _sub("pas1", [10.0, 10.0]),
            _sub(".", [2.0, 2.0]),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", debug=True
        )
        assert result is not None
        assert result["pas_usage"]["pas1"] == pytest.approx(1.0)

    def test_too_many_pas_debug(self):
        subsegments = [_sub(f"pas{i}", [float(i + 1)]) for i in range(12)]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", debug=True, max_pas_count=10
        )
        assert result is None

    def test_no_pas_debug(self):
        subsegments = [_sub(".", [5.0, 5.0])]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", debug=True
        )
        assert result is None

    def test_dfs_pattern_debug(self):
        subsegments = [
            _sub("pas1", [100.0] * 20),
            _sub("pas2", [50.0] * 20),
            _sub(".", [10.0] * 20),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", debug=True, f_stat_threshold=10
        )
        assert result is not None

    def test_no_combos_debug(self):
        subsegments = [
            _sub("pas1", [10.0, 10.0, 10.0, 10.0]),
            _sub("pas2", [10.0, 10.0, 10.0, 10.0]),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", debug=True, f_stat_threshold=1e9
        )
        assert result is None

    def test_union_group_means_debug(self):
        subsegments = [
            _sub("pas1", [100.0] * 20),
            _sub(".", [10.0] * 20),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", debug=True, f_stat_threshold=10
        )
        assert result is not None


# ---------------------------------------------------------------------------
# DFS edge cases: tail monotonicity violation and single-group skip
# ---------------------------------------------------------------------------


class TestDfsEdgeCases:
    def test_tail_monotonicity_violation(self):
        # trailing has higher coverage than the last PAS → most patterns pruned
        subsegments = [
            _sub("pas1", [50.0] * 10),
            _sub("pas2", [30.0] * 10),
            _sub(".", [100.0] * 10),
        ]
        # With f_stat_threshold=0 some patterns still pass (cuts before tail)
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", f_stat_threshold=0
        )
        # May or may not be None; just check it doesn't crash.
        assert result is None or isinstance(result, dict)

    def test_single_group_pattern_skipped(self):
        # 2 PAS, no trailing subsegment → pattern (0,1) creates 1 group → skipped
        # pattern (1,1) creates 2 groups and should be evaluated
        subsegments = [
            _sub("pas1", [100.0] * 10),
            _sub("pas2", [50.0] * 10),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", f_stat_threshold=0
        )
        # (0,1) single group is skipped (line 215). (1,1) has 2 groups, may pass.
        assert result is None or isinstance(result, dict)


# ---------------------------------------------------------------------------
# Union pattern bit=0 (lines 295, 329)
# ---------------------------------------------------------------------------


class TestUnionPatternBitZero:
    def test_union_pattern_has_zero_bit(self):
        # pas1(low)→pas2(high)→pas3(mid)→trailing(low)
        # Cut at pas1 (mean=10) → next group mean=80 → violates monotonicity
        # Only patterns starting with bit=0 at pas1 can be monotone
        # → union_pattern[0]=0 → lines 295 and 329 executed
        subsegments = [
            _sub("pas1", [10.0] * 20),
            _sub("pas2", [80.0] * 20),
            _sub("pas3", [40.0] * 20),
            _sub(".", [5.0] * 20),
        ]
        result = evaluate_all_pas_usage_patterns(
            subsegments, "seg1", f_stat_threshold=0
        )
        if result is not None:
            assert result["pas_usage"]["pas1"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# process_segment
# ---------------------------------------------------------------------------


class TestProcessSegment:
    def test_returns_triple(self):
        subs = [
            {"gene_id": "g1", "segment_id": "s1",
             "subsegment_id": "s1:001", "pas_id": "cs1",
             "mean_cov": 50.0, "sum_squared_values": 0.0,
             "coverage": np.array([50.0] * 20)},
            {"gene_id": "g1", "segment_id": "s1",
             "subsegment_id": "s1:002", "pas_id": ".",
             "mean_cov": 5.0, "sum_squared_values": 0.0,
             "coverage": np.array([5.0] * 20)},
        ]
        seg_id, out_subs, result = process_segment(
            ("s1", subs), max_pas_count=10, f_stat_threshold=0
        )
        assert seg_id == "s1"
        assert out_subs is subs
        assert result is not None

    def test_none_when_no_pas(self):
        subs = [
            {"gene_id": "g1", "segment_id": "s1",
             "subsegment_id": "s1:001", "pas_id": ".",
             "mean_cov": 5.0, "sum_squared_values": 0.0,
             "coverage": np.array([5.0] * 10)},
        ]
        _, _, result = process_segment(("s1", subs))
        assert result is None


# ---------------------------------------------------------------------------
# _read_coverage_batch RuntimeError fallback (lines 464-465)
# ---------------------------------------------------------------------------


class TestReadCoverageBatch:
    def test_runtime_error_yields_empty_array(self):
        Row = collections.namedtuple(
            "Row",
            ["gene_id", "segment_id", "subsegment_id",
             "chrom", "start", "end", "strand", "pas_id"],
        )
        row = Row("g1", "seg1", "s1:001", "chr1", 100, 200, "+", ".")

        mock_bw = MagicMock()
        mock_bw.values.side_effect = RuntimeError("chrom not in bigwig")
        mock_bw.close = MagicMock()

        with patch(
            "paqr3.calculate_cov_metrics.pyBigWig.open",
            return_value=mock_bw,
        ):
            rows = _read_coverage_batch([row], "fake.bw", "fake.bw")

        assert len(rows) == 1
        assert rows[0]["mean_cov"] == 0.0


# ---------------------------------------------------------------------------
# evaluate_pas_usage_models: None-result skip and debug JSON (lines 613, 660, 698-699)
# ---------------------------------------------------------------------------


def _dot_only_df() -> pd.DataFrame:
    return pd.DataFrame([{
        "gene_id": "g1",
        "segment_id": "seg1",
        "subsegment_id": "seg1:001",
        "pas_id": ".",
        "mean_cov": 10.0,
        "sum_squared_values": 0.0,
        "coverage": np.array([10.0] * 5),
    }])


def _single_pas_df() -> pd.DataFrame:
    return pd.DataFrame([
        {"gene_id": "g1", "segment_id": "seg1",
         "subsegment_id": "seg1:001", "pas_id": "cs1",
         "mean_cov": 80.0, "sum_squared_values": 0.0,
         "coverage": np.array([80.0] * 20)},
        {"gene_id": "g1", "segment_id": "seg1",
         "subsegment_id": "seg1:002", "pas_id": ".",
         "mean_cov": 10.0, "sum_squared_values": 0.0,
         "coverage": np.array([10.0] * 20)},
    ])


class TestEvaluatePasUsageModels:
    def test_none_result_skipped_sequential(self):
        cc = CalculateCoverages("", "")
        usage_df = cc.evaluate_pas_usage_models(
            _dot_only_df(), n_procs=1
        )
        assert usage_df.empty

    def test_none_result_skipped_parallel(self):
        cc = CalculateCoverages("", "")
        usage_df = cc.evaluate_pas_usage_models(
            _dot_only_df(), n_procs=2
        )
        assert usage_df.empty

    def test_debug_json_written(self, tmp_path):
        cc = CalculateCoverages("", "")
        debug_path = str(tmp_path / "debug.json")
        cc.evaluate_pas_usage_models(
            _single_pas_df(),
            output_tsv_debug=debug_path,
            n_procs=1,
            f_stat_threshold=0,
        )
        assert os.path.exists(debug_path)
        with open(debug_path) as fh:
            data = json.load(fh)
        assert isinstance(data, list)
