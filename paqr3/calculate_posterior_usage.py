"""Posterior PAS usage blending for PAQR3."""

import logging

import pandas as pd  # type: ignore

logger = logging.getLogger(__name__)


class CalculatePosteriorUsage:
    """Blend atlas RPM with observed drop RPM to produce posterior usage.

    posterior_rpm = weight * observed_rpm + (1 - weight) * atlas_rpm,
    then normalised within each segment.

    Args:
        weight: Mixing weight for observed RPM. Must be in [0, 1].
    """

    def __init__(self, weight: float = 0.1) -> None:
        self.weight = weight

    def compute(self, usage_df: pd.DataFrame) -> pd.DataFrame:
        """Compute posterior relative usage.

        Converts rna_drop_cov to observed_rpm (normalised to 1e6),
        blends with atlas_rpm using self.weight, then normalises
        posterior_rpm and atlas_rpm within each segment_id.

        Args:
            usage_df: Per-PAS DataFrame with rna_drop_cov, atlas_rpm,
                and segment_id columns. atlas_rpm must already be
                present (carried through from the segments TSV).

        Returns:
            Copy of usage_df with added columns observed_rpm,
            posterior_rpm, posterior_rel_usage, atlas_rel_usage.
        """
        df = usage_df.copy()

        total_drop = df["rna_drop_cov"].sum()
        df["observed_rpm"] = (
            df["rna_drop_cov"] / total_drop * 1e6 if total_drop > 0 else 0.0
        )

        w = self.weight
        atlas = df["atlas_rpm"].fillna(0.0)
        df["posterior_rpm"] = w * df["observed_rpm"] + (1 - w) * atlas

        df["posterior_rel_usage"] = df.groupby("segment_id")[
            "posterior_rpm"
        ].transform(lambda x: x / x.sum() if x.sum() > 0 else 0.0)

        df["atlas_rel_usage"] = df.groupby("segment_id")[
            "atlas_rpm"
        ].transform(lambda x: x / x.sum() if x.sum() > 0 else 0.0)

        return df
