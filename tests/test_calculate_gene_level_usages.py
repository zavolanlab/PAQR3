"""Unit tests for paqr3.calculate_gene_level_usages."""

import pandas as pd
import pytest

from paqr3.calculate_gene_level_usages import CalculateGeneLevelUsage


def _df(**kwargs) -> pd.DataFrame:
    return pd.DataFrame(kwargs)


class TestCalculateGeneLevelUsage:
    def test_single_gene_single_pas_gives_one(self):
        df = _df(
            gene_id=["g1"],
            subsegment_id=["g1:T001:001"],
            posterior_rpm=[100.0],
        )
        result = CalculateGeneLevelUsage(df).compute()
        assert result["gene_level_usage"].iloc[0] == pytest.approx(1.0)

    def test_single_gene_two_pas_proportional(self):
        df = _df(
            gene_id=["g1", "g1"],
            subsegment_id=["g1:T001:001", "g1:T001:002"],
            posterior_rpm=[100.0, 300.0],
        )
        result = CalculateGeneLevelUsage(df).compute()
        idx = result.set_index("subsegment_id")
        assert idx.loc["g1:T001:001", "gene_level_usage"] == pytest.approx(0.25)
        assert idx.loc["g1:T001:002", "gene_level_usage"] == pytest.approx(0.75)

    def test_gene_usages_sum_to_one(self):
        df = _df(
            gene_id=["g1", "g1", "g1"],
            subsegment_id=["g1:T001:001", "g1:T001:002", "g1:T002:001"],
            posterior_rpm=[200.0, 300.0, 500.0],
        )
        result = CalculateGeneLevelUsage(df).compute()
        assert result["gene_level_usage"].sum() == pytest.approx(1.0)

    def test_multi_gene_normalised_independently(self):
        df = _df(
            gene_id=["g1", "g1", "g2", "g2"],
            subsegment_id=[
                "g1:T001:001",
                "g1:T001:002",
                "g2:T001:001",
                "g2:T001:002",
            ],
            posterior_rpm=[200.0, 200.0, 100.0, 400.0],
        )
        result = CalculateGeneLevelUsage(df).compute()
        idx = result.set_index("subsegment_id")
        # g1: both equal
        assert idx.loc["g1:T001:001", "gene_level_usage"] == pytest.approx(0.5)
        assert idx.loc["g1:T001:002", "gene_level_usage"] == pytest.approx(0.5)
        # g2: 1:4 split
        assert idx.loc["g2:T001:001", "gene_level_usage"] == pytest.approx(0.2)
        assert idx.loc["g2:T001:002", "gene_level_usage"] == pytest.approx(0.8)

    def test_zero_total_rpm_fills_zero(self):
        df = _df(
            gene_id=["g1", "g1"],
            subsegment_id=["g1:T001:001", "g1:T001:002"],
            posterior_rpm=[0.0, 0.0],
        )
        result = CalculateGeneLevelUsage(df).compute()
        assert (result["gene_level_usage"] == 0.0).all()

    def test_output_columns(self):
        df = _df(
            gene_id=["g1"],
            subsegment_id=["g1:T001:001"],
            posterior_rpm=[1.0],
        )
        result = CalculateGeneLevelUsage(df).compute()
        assert list(result.columns) == ["subsegment_id", "gene_level_usage"]

    def test_output_row_count_matches_input(self):
        df = _df(
            gene_id=["g1", "g1", "g2"],
            subsegment_id=["g1:T001:001", "g1:T001:002", "g2:T001:001"],
            posterior_rpm=[100.0, 200.0, 300.0],
        )
        result = CalculateGeneLevelUsage(df).compute()
        assert len(result) == len(df)

    def test_input_dataframe_not_mutated(self):
        df = _df(
            gene_id=["g1"],
            subsegment_id=["g1:T001:001"],
            posterior_rpm=[100.0],
        )
        original_cols = set(df.columns)
        CalculateGeneLevelUsage(df).compute()
        assert set(df.columns) == original_cols
