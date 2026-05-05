"""Unit tests for paqr3.calculate_posterior_usage."""

import pandas as pd
import pytest

from paqr3.calculate_posterior_usage import CalculatePosteriorUsage


def _df(**kwargs) -> pd.DataFrame:
    """Build a minimal usage DataFrame from keyword columns."""
    return pd.DataFrame(kwargs)


class TestCalculatePosteriorUsage:
    def test_observed_rpm_normalised_to_1e6(self):
        df = _df(
            segment_id=["s1"],
            rna_drop_cov=[500.0],
            atlas_rpm=[0.0],
        )
        result = CalculatePosteriorUsage(weight=1.0).compute(df)
        # Single row: observed_rpm must equal 1e6
        assert result["observed_rpm"].iloc[0] == pytest.approx(1e6)

    def test_observed_rpm_proportional(self):
        df = _df(
            segment_id=["s1", "s1"],
            rna_drop_cov=[100.0, 300.0],
            atlas_rpm=[0.0, 0.0],
        )
        result = CalculatePosteriorUsage(weight=1.0).compute(df)
        rpm = result["observed_rpm"].tolist()
        assert rpm[1] == pytest.approx(3 * rpm[0])

    def test_zero_total_drop_gives_zero_observed_rpm(self):
        df = _df(
            segment_id=["s1", "s1"],
            rna_drop_cov=[0.0, 0.0],
            atlas_rpm=[50.0, 50.0],
        )
        result = CalculatePosteriorUsage(weight=0.5).compute(df)
        assert (result["observed_rpm"] == 0.0).all()

    def test_weight_zero_is_pure_atlas(self):
        df = _df(
            segment_id=["s1", "s1"],
            rna_drop_cov=[200.0, 100.0],
            atlas_rpm=[300.0, 100.0],
        )
        result = CalculatePosteriorUsage(weight=0.0).compute(df)
        pd.testing.assert_series_equal(
            result["posterior_rpm"].reset_index(drop=True),
            result["atlas_rpm"].reset_index(drop=True),
            check_names=False,
        )

    def test_weight_one_is_pure_observed(self):
        df = _df(
            segment_id=["s1", "s1"],
            rna_drop_cov=[200.0, 100.0],
            atlas_rpm=[0.0, 0.0],
        )
        result = CalculatePosteriorUsage(weight=1.0).compute(df)
        pd.testing.assert_series_equal(
            result["posterior_rpm"].reset_index(drop=True),
            result["observed_rpm"].reset_index(drop=True),
            check_names=False,
        )

    def test_posterior_rpm_blends_correctly(self):
        df = _df(
            segment_id=["s1"],
            rna_drop_cov=[1000.0],
            atlas_rpm=[200.0],
        )
        w = 0.4
        result = CalculatePosteriorUsage(weight=w).compute(df)
        expected = w * 1e6 + (1 - w) * 200.0
        assert result["posterior_rpm"].iloc[0] == pytest.approx(expected)

    def test_posterior_rel_usage_sums_to_one_per_segment(self):
        df = _df(
            segment_id=["s1", "s1", "s1"],
            rna_drop_cov=[100.0, 200.0, 300.0],
            atlas_rpm=[10.0, 20.0, 30.0],
        )
        result = CalculatePosteriorUsage(weight=0.1).compute(df)
        for _, grp in result.groupby("segment_id"):
            assert grp["posterior_rel_usage"].sum() == pytest.approx(1.0)

    def test_atlas_rel_usage_sums_to_one_per_segment(self):
        df = _df(
            segment_id=["s1", "s1"],
            rna_drop_cov=[100.0, 100.0],
            atlas_rpm=[200.0, 400.0],
        )
        result = CalculatePosteriorUsage(weight=0.1).compute(df)
        assert result["atlas_rel_usage"].sum() == pytest.approx(1.0)
        assert result["atlas_rel_usage"].iloc[1] == pytest.approx(
            result["atlas_rel_usage"].iloc[0] * 2
        )

    def test_multi_segment_normalisation_is_independent(self):
        df = _df(
            segment_id=["s1", "s1", "s2", "s2"],
            rna_drop_cov=[100.0, 300.0, 50.0, 50.0],
            atlas_rpm=[10.0, 30.0, 5.0, 5.0],
        )
        result = CalculatePosteriorUsage(weight=0.1).compute(df)
        for seg, grp in result.groupby("segment_id"):
            assert grp["posterior_rel_usage"].sum() == pytest.approx(
                1.0
            ), f"segment {seg} doesn't sum to 1"

    def test_output_has_required_columns(self):
        df = _df(
            segment_id=["s1"],
            rna_drop_cov=[100.0],
            atlas_rpm=[50.0],
        )
        result = CalculatePosteriorUsage().compute(df)
        for col in (
            "observed_rpm",
            "posterior_rpm",
            "posterior_rel_usage",
            "atlas_rel_usage",
        ):
            assert col in result.columns

    def test_atlas_nan_treated_as_zero(self):
        df = _df(
            segment_id=["s1"],
            rna_drop_cov=[100.0],
            atlas_rpm=[float("nan")],
        )
        result = CalculatePosteriorUsage(weight=0.5).compute(df)
        assert not result["posterior_rpm"].isna().any()

    def test_input_dataframe_not_mutated(self):
        df = _df(
            segment_id=["s1"],
            rna_drop_cov=[100.0],
            atlas_rpm=[50.0],
        )
        original_cols = set(df.columns)
        CalculatePosteriorUsage().compute(df)
        assert set(df.columns) == original_cols
