from datetime import datetime
import pyBigWig
import pandas as pd
import numpy as np

def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def calculate_mean_coverage(genes, coverage_bw):
    """Calculates mean squared coverage for PAS sites using bigWig (corrected bounds)."""

    log_message("Calculating mean squared coverage for PAS sites (using bigWig)...")

    bw = pyBigWig.open(coverage_bw)

    results = []

    for gene in genes.values():
        for segment in gene.segments:
            if "overlapping_pas" in segment.attributes:
                pas_ids = segment.attributes["overlapping_pas"]

                for pas_id in pas_ids:
                    pas_info = pas_id.split(":")
                    pas_chrom = pas_info[0]
                    pas_pos = int(pas_info[1])
                    pas_strand = pas_info[2]

                    segment_start = segment.start
                    segment_end = segment.end

                    if pas_strand == "+":
                        upstream_start = segment_start
                        upstream_end = pas_pos
                        downstream_start = pas_pos + 1  # Corrected downstream start
                        downstream_end = segment_end
                    elif pas_strand == "-":
                        upstream_start = pas_pos
                        upstream_end = segment_end
                        downstream_start = segment_start
                        downstream_end = pas_pos - 1 # Corrected downstream end

                    # Ensure valid intervals (start < end)
                    if upstream_start < upstream_end:
                        try:
                            upstream_coverage = bw.values(pas_chrom, upstream_start, upstream_end, numpy=True)
                        except ValueError:  # Chromosome not in bigWig or invalid interval
                            upstream_coverage = np.array([])
                    else:
                        upstream_coverage = np.array([])

                    if downstream_start < downstream_end: # Check for correct interval bounds
                        try:
                            downstream_coverage = bw.values(pas_chrom, downstream_start, downstream_end, numpy=True)
                        except ValueError:  # Chromosome not in bigWig or invalid interval
                            downstream_coverage = np.array([])
                    else:
                        downstream_coverage = np.array([])

                    upstream_mean = np.nanmean(upstream_coverage) if len(upstream_coverage) > 0 else 0  # Use nanmean to ignore NAs
                    downstream_mean = np.nanmean(downstream_coverage) if len(downstream_coverage) > 0 else 0

                    upstream_mean_sq = upstream_mean**2
                    downstream_mean_sq = downstream_mean**2

                    results.append(
                        [
                            gene.gene_id,
                            pas_id,
                            segment.attributes["segment_number"],
                            f"{segment.chrom}:{segment.start}-{segment.end}:{segment.strand}",
                            upstream_mean_sq,
                            downstream_mean_sq,
                        ]
                    )

    bw.close()

    results_df = pd.DataFrame(
        results,
        columns=[
            "gene_id",
            "pas_id",
            "segment_number",
            "segment_coordinates",
            "upstream_mean_sq_coverage",
            "downstream_mean_sq_coverage",
        ],
    )

    return results_df


def write_coverage_results(results_df, output_file):
    """Writes the coverage results to a TSV file."""
    results_df.to_csv(output_file, sep="\t", index=False)
    log_message(f"Coverage results written to {output_file}")
