import logging
import pandas as pd  # type: ignore


def log_message(msg):
    logging.info(msg)


class CalculatePosteriorUsage:
    """
    Blend atlas‐derived RPMs (from pas_df) with observed drop‐derived RPMs
    (from usage_df) to get posterior relative usage per PAS.
    """

    def __init__(self, atlas_df: pd.DataFrame):
        # atlas_df must have columns: ['pas_id','atlas_rpm']
        self.atlas_df = atlas_df.copy()
        # Ensure atlas pas_id is integer
        self.atlas_df["pas_id"] = self.atlas_df["pas_id"].astype(int)

    def compute(self, usage_df: pd.DataFrame) -> pd.DataFrame:
        df = usage_df.copy()

        # Coerce usage_df['pas_id'] to int so it matches atlas_df
        df["pas_id"] = df["pas_id"].astype(int)

        # 1) observed RPM from drop
        total_drop = df["rna_drop_cov"].sum()
        log_message(f"Total rna_drop_cov sum = {total_drop:.3f}")
        df["observed_rpm"] = (
            df["rna_drop_cov"] / total_drop * 1e6 if total_drop > 0 else 0.0
        )

        # 2) merge atlas RPM & compute posterior RPM
        merged = df.merge(
            self.atlas_df[["pas_id", "atlas_rpm"]], on="pas_id", how="left"
        )
        merged["atlas_rpm"] = merged["atlas_rpm"].fillna(0.0)
        merged["posterior_rpm"] = merged["atlas_rpm"] + merged["observed_rpm"]

        # 3) normalize within each segment
        merged["posterior_rel_usage"] = merged.groupby("segment_id")[
            "posterior_rpm"
        ].transform(lambda x: x / x.sum())

        return merged
