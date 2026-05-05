"""Gene-level PAS usage computation for PAQR3."""

import pandas as pd  # type: ignore


class CalculateGeneLevelUsage:
    """Compute per-PAS gene-level usage fraction.

    gene_level_usage = posterior_rpm / sum(posterior_rpm for all PAS
    in the same gene). For single-segment genes this equals
    posterior_rel_usage; for multi-segment genes it correctly weights
    expression across segments.

    Args:
        posterior_df: Per-PAS DataFrame with gene_id, subsegment_id,
            and posterior_rpm columns.
    """

    def __init__(self, posterior_df: pd.DataFrame):
        self.posterior_df = posterior_df.copy()

    def compute(self) -> pd.DataFrame:
        """Compute gene_level_usage for each PAS.

        Returns:
            DataFrame with columns subsegment_id and gene_level_usage.
        """
        df = self.posterior_df[
            ["gene_id", "subsegment_id", "posterior_rpm"]
        ].copy()

        gene_total = df.groupby("gene_id")["posterior_rpm"].transform("sum")
        df["gene_level_usage"] = (
            df["posterior_rpm"].div(gene_total).fillna(0.0)
        )

        return df[["subsegment_id", "gene_level_usage"]]
