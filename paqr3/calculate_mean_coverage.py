from datetime import datetime
import pyBigWig
import pandas as pd
import numpy as np

def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def calculate_mean_coverage(genes, coverage_bw_pos, coverage_bw_neg): #Added genes parameter
    """Calculates mean coverage for subsegments using bigWig."""

    log_message("Calculating mean coverage for subsegments (using bigWig)...")

    bw_pos = pyBigWig.open(coverage_bw_pos)
    bw_neg = pyBigWig.open(coverage_bw_neg)

    results = []

    for gene in genes.values():
        for segment in gene.segments:
            if hasattr(segment, "subsegments"):
                for subsegment in segment.subsegments:
                    subsegment_start = subsegment.start
                    subsegment_end = subsegment.end
                    subsegment_strand = subsegment.strand

                    if subsegment_strand == "+":
                        bw = bw_pos
                    elif subsegment_strand == "-":
                        bw = bw_neg
                    else:
                        raise ValueError(f"Invalid strand: {subsegment_strand}")

                    if subsegment_start < subsegment_end:
                        try:
                            coverage = bw.values(subsegment.chrom, subsegment_start, subsegment_end, numpy=True)
                        except ValueError:
                            coverage = np.array([])
                    else:
                        coverage = np.array([])

                    mean_coverage = np.nanmean(coverage) if len(coverage) > 0 else 0

                    results.append(
                        [
                            gene.gene_id,
                            segment.attributes["segment_number"],
                            subsegment.attributes["subsegment_number"],
                            f"{subsegment.chrom}:{subsegment.start}-{subsegment.end}:{subsegment.strand}",
                            mean_coverage,
                        ]
                    )

    bw_pos.close()
    bw_neg.close()

    results_df = pd.DataFrame(
        results,
        columns=[
            "gene_id",
            "segment_number",
            "subsegment_number",
            "subsegment_coordinates",
            "mean_coverage",
        ],
    )

    return results_df


def write_coverage_results(genes, results_df, output_genes_tsv, output_segments_tsv, output_subsegments_tsv): #Added genes parameter
    """Writes the coverage results to TSV files."""

    if results_df.empty:  # Check if DataFrame is empty
        genes_df = pd.DataFrame(columns=["unique_gene_id", "gene_id", "gene_name", "chrom", "start", "end"])
        segments_df = pd.DataFrame(columns=["gene_unique_id", "segment_number", "segment_coordinates"])
        subsegments_df = pd.DataFrame(columns=["gene_unique_id", "segment_number", "subsegment_number", "subsegment_coordinates", "mean_coverage"])
    else:
        # Create unique gene IDs and DataFrames
        gene_mapping = {}
        unique_gene_id_counter = 1
        genes_data = []
        segments_data = []
        subsegments_data = []

        for _, row in results_df.iterrows():
            gene_id = row["gene_id"]
            if gene_id not in gene_mapping:
                gene_mapping[gene_id] = unique_gene_id_counter
                unique_gene_id_counter += 1

            unique_gene_id = gene_mapping[gene_id]

            gene = next((g for g in genes.values() if g.gene_id == gene_id), None)
            if gene:
                # Get chromosome from the first transcript's first region (assuming all transcripts are on the same chromosome)
                first_transcript = next(iter(gene.transcripts.values()), None)
                chrom = first_transcript.regions[0].chrom if first_transcript and first_transcript.regions else ""

                genes_data.append([unique_gene_id, gene_id, gene.attributes.get("gene_name", ""), chrom, min(list(r.start for t in gene.transcripts.values() for r in t.regions)), max(list(r.end for t in gene.transcripts.values() for r in t.regions))])  # Corrected: Added chrom
            else:
                genes_data.append([unique_gene_id, gene_id, "", "", "", ""])

            segments_data.append([unique_gene_id, row["segment_number"], row["subsegment_coordinates"].rsplit(":", 1)[0]])
            subsegments_data.append([unique_gene_id, row["segment_number"], row["subsegment_number"], row["subsegment_coordinates"], row["mean_coverage"]])

        genes_df = pd.DataFrame(genes_data, columns=["unique_gene_id", "gene_id", "gene_name", "chrom", "start", "end"])
        segments_df = pd.DataFrame(segments_data, columns=["gene_unique_id", "segment_number", "segment_coordinates"])
        subsegments_df = pd.DataFrame(subsegments_data, columns=["gene_unique_id", "segment_number", "subsegment_number", "subsegment_coordinates", "mean_coverage"])

    genes_df.to_csv(output_genes_tsv, sep="\t", index=False)
    segments_df.to_csv(output_segments_tsv, sep="\t", index=False)
    subsegments_df.to_csv(output_subsegments_tsv, sep="\t", index=False)

    log_message(f"Coverage results written to {output_genes_tsv}, {output_segments_tsv}, and {output_subsegments_tsv}")
