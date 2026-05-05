import pandas as pd  # type: ignore


class CalculateGeneLevelUsage:
    """Compute each PAS's share of total gene-level posterior expression.

    For each PAS: ``gene_level_usage = posterior_rpm / sum(posterior_rpm
    across all PAS in the same gene)``.  For a single-segment gene this is
    identical to ``posterior_rel_usage``; for multi-segment genes it correctly
    weights each segment's expression level when comparing PAS across segments.

    Returns a DataFrame with columns ``subsegment_id`` and
    ``gene_level_usage``.
    """

    def __init__(self, posterior_df: pd.DataFrame):
        self.posterior_df = posterior_df.copy()

    def compute(self) -> pd.DataFrame:
        df = self.posterior_df[
            ["gene_id", "subsegment_id", "posterior_rpm"]
        ].copy()

        gene_total = df.groupby("gene_id")["posterior_rpm"].transform("sum")
        df["gene_level_usage"] = (
            df["posterior_rpm"].div(gene_total).fillna(0.0)
        )

        return df[["subsegment_id", "gene_level_usage"]]
