from datetime import datetime
import pyBigWig  # type: ignore
import pandas as pd  # type: ignore
import numpy as np


def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


def calculate_mean_coverage(subsegments_df, coverage_bw_pos, coverage_bw_neg):
    """Calculates mean coverage and returns a DataFrame."""

    log_message("Calculating mean coverage...")

    bw_pos = pyBigWig.open(coverage_bw_pos)
    bw_neg = pyBigWig.open(coverage_bw_neg)

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

        results.append([subsegment_id, mean_coverage])

    bw_pos.close()
    bw_neg.close()

    results_df = pd.DataFrame(results, columns=["subsegment_id", "mean_cov"])
    return results_df


def write_coverage_results(results_df, output_tsv):
    """Writes the coverage results DataFrame to a TSV file."""

    results_df.to_csv(output_tsv, sep="\t", index=False)
    log_message(f"Coverage results written to {output_tsv}")
