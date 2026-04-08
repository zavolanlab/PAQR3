import pandas as pd  # type: ignore


class CalculateGeneLevelUsage:
    """
    Aggregate PAS-level posterior relative usage to gene-level usage.

    For each gene:
      - Multiply PAS posterior_rel_usage by its RNA-seq drop coverage (rna_drop_cov)
      - Sum this product over the gene
      - Normalize by total rna_drop_cov per gene

    Output: gene_id, gene_weighted_usage (float between 0-1)
    """

    def __init__(self, posterior_df: pd.DataFrame):
        self.posterior_df = posterior_df.copy()

    def compute(self) -> pd.DataFrame:
        df = self.posterior_df.copy()

        # Ensure numeric and non-null values
        df = df[
            df[["gene_id", "rna_drop_cov", "posterior_rel_usage"]]
            .notnull()
            .all(axis=1)
        ].copy()

        # Multiply posterior usage by drop coverage
        df["weighted_usage"] = df["posterior_rel_usage"] * df["rna_drop_cov"]

        # Group by gene_id and compute totals
        gene_df = (
            df.groupby("gene_id")
            .agg(
                total_weighted_usage=("weighted_usage", "sum"),
                total_drop_cov=("rna_drop_cov", "sum"),
            )
            .reset_index()
        )

        # Compute gene-level weighted usage
        gene_df["gene_weighted_usage"] = gene_df.apply(
            lambda row: (
                row["total_weighted_usage"] / row["total_drop_cov"]
                if row["total_drop_cov"] > 0
                else 0.0
            ),
            axis=1,
        )

        # Keep only final output columns
        result = gene_df[["gene_id", "gene_weighted_usage"]].copy()
        return result
