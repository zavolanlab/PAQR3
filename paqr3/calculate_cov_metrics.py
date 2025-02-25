import logging
import pandas as pd  # type: ignore
import pyBigWig  # type: ignore
import numpy as np


def log_message(message):
    logging.info(message)


class CalculateCoverages:
    def __init__(self, coverage_bw_pos, coverage_bw_neg):
        self.coverage_bw_pos = coverage_bw_pos
        self.coverage_bw_neg = coverage_bw_neg

    def calculate_coverage_metrics(self, subsegments_df):
        """Calculates mean coverage and sum of squared coverage values."""

        bw_pos = pyBigWig.open(self.coverage_bw_pos)
        bw_neg = pyBigWig.open(self.coverage_bw_neg)

        results = []

        for _, row in subsegments_df.iterrows():
            subsegment_id = f"{row['gene_unique_id']}.{row['segment_number']}.{row['subsegment_number']}"
            chrom = row["chrom"]
            start = row["start"]
            end = row["end"]
            strand = row["strand"]

            if strand == "+":
                bw = bw_pos
            elif strand == "-":
                bw = bw_neg
            else:
                raise ValueError(f"Invalid strand: {strand}")

            if start < end:
                try:
                    coverage = bw.values(chrom, start, end, numpy=True)
                except ValueError:
                    coverage = np.array([])
            else:
                coverage = np.array([])

            mean_coverage = np.nanmean(coverage) if len(coverage) > 0 else 0
            sum_squared_values = (
                np.nansum(coverage**2) if len(coverage) > 0 else 0
            )

            results.append([subsegment_id, mean_coverage, sum_squared_values])

        bw_pos.close()
        bw_neg.close()

        results_df = pd.DataFrame(
            results,
            columns=["subsegment_id", "mean_cov", "sum_squared_values"],
        )
        return results_df

    def write_coverage_results(self, results_df, output_tsv):
        """Writes the coverage results DataFrame to a TSV file."""

        results_df.to_csv(output_tsv, sep="\t", index=False)
        log_message(f"Coverage results written to {output_tsv}")

    def run(self, subsegments_df, output_coverage_tsv):
        log_message(
            "Calculating mean coverage and sum of squared values for subsegments..."
        )
        coverage_results = self.calculate_coverage_metrics(subsegments_df)

        log_message("Writing coverage results to TSV...")
        self.write_coverage_results(coverage_results, output_coverage_tsv)
