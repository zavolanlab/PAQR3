from datetime import datetime
import pyBigWig  # type: ignore
import pandas as pd  # type: ignore
import numpy as np


def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


def calculate_mean_coverage(genes, coverage_bw_pos, coverage_bw_neg):
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
                            coverage = bw.values(
                                subsegment.chrom,
                                subsegment_start,
                                subsegment_end,
                                numpy=True,
                            )
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


def write_coverage_results(
    genes, results_df, out_genes_tsv, out_segments_tsv, out_subsegments_tsv
):
    """Writes the coverage results to TSV files."""

    if results_df.empty:  # Check if DataFrame is empty
        genes_df = pd.DataFrame(
            columns=[
                "unique_gene_id",
                "gene_id",
                "gene_name",
                "chrom",
                "start",
                "end",
                "strand",
            ]
        )
        segments_df = pd.DataFrame(
            columns=[
                "gene_unique_id",
                "segment_number",
                "chrom",
                "start",
                "end",
                "strand",
                "overlapping_pas",
            ]
        )
        subsegments_df = pd.DataFrame(
            columns=[
                "gene_unique_id",
                "segment_number",
                "subsegment_number",
                "chrom",
                "start",
                "end",
                "strand",
                "mean_coverage",
            ]
        )
    else:
        gene_mapping = {}
        unique_gene_id_counter = 1
        segments_data = []
        subsegments_data = []

        unique_genes = {}
        processed_segments = set()

        for _, row in results_df.iterrows():
            gene_id = row["gene_id"]
            if gene_id not in gene_mapping:
                gene_mapping[gene_id] = unique_gene_id_counter
                unique_gene_id_counter += 1

            unique_gene_id = gene_mapping[gene_id]

            gene = next((g for g in genes.values() if g.gene_id == gene_id), None)
            if gene:
                first_transcript = next(iter(gene.transcripts.values()), None)
                chrom = (
                    first_transcript.regions[0].chrom
                    if first_transcript and first_transcript.regions
                    else ""
                )
                strand = first_transcript.strand if first_transcript else ""

                unique_genes[unique_gene_id] = [
                    unique_gene_id,
                    gene_id,
                    gene.attributes.get("gene_name", ""),
                    chrom,
                    min(
                        list(
                            r.start
                            for t in gene.transcripts.values()
                            for r in t.regions
                        )
                    ),
                    max(
                        list(
                            r.end for t in gene.transcripts.values() for r in t.regions
                        )
                    ),
                    strand,
                ]
            else:
                unique_genes[unique_gene_id] = [
                    unique_gene_id,
                    gene_id,
                    "",
                    "",
                    "",
                    "",
                    "",
                ]

            segment_key = (unique_gene_id, row["segment_number"])
            segment = next(
                (
                    s
                    for s in gene.segments
                    if s.attributes["segment_number"] == row["segment_number"]
                ),
                None,
            )
            if segment and segment_key not in processed_segments:
                overlapping_pas_str = ",".join(
                    segment.attributes.get("overlapping_pas", [])
                )
                segments_data.append(
                    [
                        unique_gene_id,
                        row["segment_number"],
                        segment.chrom,
                        segment.start,
                        segment.end,
                        segment.strand,
                        overlapping_pas_str,
                    ]
                )
                processed_segments.add(segment_key)

            if (
                segment
                and "overlapping_pas" in segment.attributes
                and segment.attributes["overlapping_pas"]
            ):
                subsegment_coords = row["subsegment_coordinates"].split(":")
                subsegment_chrom = (
                    subsegment_coords[0] if len(subsegment_coords) > 0 else ""
                )
                subsegment_start_end = (
                    subsegment_coords[1] if len(subsegment_coords) > 1 else ""
                )
                subsegment_strand = (
                    subsegment_coords[2] if len(subsegment_coords) > 2 else ""
                )

                subsegment_start = (
                    int(subsegment_start_end.split("-")[0])
                    if subsegment_start_end and "-" in subsegment_start_end
                    else ""
                )
                subsegment_end = (
                    int(subsegment_start_end.split("-")[1])
                    if subsegment_start_end and "-" in subsegment_start_end
                    else ""
                )

                subsegments_data.append(
                    [
                        unique_gene_id,
                        row["segment_number"],
                        row["subsegment_number"],
                        subsegment_chrom,
                        subsegment_start,
                        subsegment_end,
                        subsegment_strand,
                        row["mean_coverage"],
                    ]
                )

        genes_df = pd.DataFrame(
            list(unique_genes.values()),
            columns=[
                "unique_gene_id",
                "gene_id",
                "gene_name",
                "chrom",
                "start",
                "end",
                "strand",
            ],
        )
        segments_df = pd.DataFrame(
            segments_data,
            columns=[
                "gene_unique_id",
                "segment_number",
                "chrom",
                "start",
                "end",
                "strand",
                "overlapping_pas",
            ],
        )
        subsegments_df = pd.DataFrame(
            subsegments_data,
            columns=[
                "gene_unique_id",
                "segment_number",
                "subsegment_number",
                "chrom",
                "start",
                "end",
                "strand",
                "mean_coverage",
            ],
        )

    genes_df.to_csv(out_genes_tsv, sep="\t", index=False)
    segments_df.to_csv(out_segments_tsv, sep="\t", index=False)
    subsegments_df.to_csv(out_subsegments_tsv, sep="\t", index=False)

    log_message(
        f"Coverage results written to {out_genes_tsv}, {out_segments_tsv}, and {out_subsegments_tsv}"
    )
