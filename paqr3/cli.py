import argparse
import os
from datetime import datetime
from gtfparse import read_gtf  # type: ignore
from paqr3.version import __version__
from paqr3.construct_segments import (
    extend_gene_coordinates,
    extend_exon_downstream,
    define_exons_introns,
    construct_segments,
    identify_pas_in_segments,
    write_segments_pas_to_tsv,
    write_segments_pas_to_gtf,
    create_subsegments_dataframe,
)
from paqr3.calculate_coverages import (
    calculate_mean_coverage,
    write_coverage_results,
)
from paqr3.models import Gene, Transcript, Region
import pandas as pd  # type: ignore


# Utility function for timestamped messages
def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


def parse_gtf_to_genes(gtf_data):
    """Parses GTF data (pandas DataFrame) into a dictionary of Gene objects."""
    genes = {}
    for _, row in gtf_data.iterrows():
        if row["feature"] == "gene":
            gene_id = row["gene_id"]
            gene_attributes = {
                k: v
                for k, v in row.items()
                if k
                not in [
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                ]
            }
            gene = Gene(gene_id, attributes=gene_attributes)
            genes[gene_id] = gene
        elif row["feature"] == "transcript":
            gene_id = row["gene_id"]
            transcript_id = row["transcript_id"]
            transcript_attributes = {
                k: v
                for k, v in row.items()
                if k
                not in [
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "gene_id",
                ]
            }

            if gene_id in genes:
                transcript = Transcript(
                    transcript_id,
                    row["strand"],
                    attributes=transcript_attributes,
                )
                genes[gene_id].add_transcript(transcript)
        elif row["feature"] == "exon":
            gene_id = row["gene_id"]
            transcript_id = row["transcript_id"]
            exon_attributes = {
                k: v
                for k, v in row.items()
                if k
                not in [
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "gene_id",
                    "transcript_id",
                ]
            }
            if (
                gene_id in genes
                and transcript_id in genes[gene_id].transcripts
            ):
                exon = Region(
                    region_type="exon",
                    chrom=row["seqname"],
                    start=row["start"],
                    end=row["end"],
                    strand=row["strand"],
                    attributes=exon_attributes,
                )
                genes[gene_id].transcripts[transcript_id].add_region(exon)

    return genes


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter and construct genomic segments from RNA-Seq annotations."
        )
    )
    parser.add_argument(
        "-de",
        "--downstream_exon_extension",
        type=int,
        default=200,
        help="Number of bases to extend terminal exons. Default: 200.",
    )
    parser.add_argument(
        "-a",
        "--annotation",
        type=str,
        required=True,
        help="Path to the annotation GTF file.",
    )
    parser.add_argument(
        "-pa",
        "--pas_atlas",
        type=str,
        required=True,
        help="Path to the PAS atlas BED file.",
    )
    parser.add_argument(
        "-c",
        "--coverage",
        type=str,
        required=True,
        nargs=2,
        help=(
            "Paths to the positive and negative strand coverage bigWig files"
            "(e.g., pos.bw neg.bw)."
        ),
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to the output directory.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=(
            f"PAQR3 v{__version__},"
            "(c) 2025 by Zavolab (zavolab-biozentrum@unibas.ch)"
        ),
        help="show version information and exit",
    )

    args = parser.parse_args()

    log_message("Starting genomic segment construction pipeline...")

    # Extract sample name from coverage files
    sample_name_pos = os.path.basename(args.coverage[0]).split("_")[0]
    sample_name_neg = os.path.basename(args.coverage[1]).split("_")[0]

    if sample_name_pos != sample_name_neg:
        raise ValueError("Coverage files must have the same sample name.")

    sample_name = sample_name_pos

    output_dir = os.path.join(args.output_dir, f"{sample_name}_results")
    os.makedirs(output_dir, exist_ok=True)

    output_gtf = os.path.join(output_dir, f"{sample_name}.gtf")
    output_genes_tsv = os.path.join(output_dir, f"{sample_name}_genes.tsv")
    output_segments_tsv = os.path.join(
        output_dir, f"{sample_name}_segments.tsv"
    )
    output_subsegments_tsv = os.path.join(
        output_dir, f"{sample_name}_subsegments.tsv"
    )
    output_coverage_tsv = os.path.join(
        output_dir, f"{sample_name}_coverage.tsv"
    )

    # Step 1: Read annotation
    log_message("Reading annotation GTF...")
    gtf_data = read_gtf(args.annotation, result_type="pandas")

    # Step 2: Parse GTF to Gene/Transcript/Region objects
    log_message("Parsing GTF to Gene/Transcript/Region objects...")
    genes = parse_gtf_to_genes(gtf_data)

    # Step 2.5: Extend gene coordinates
    log_message("Extending gene coordinates...")
    genes = extend_gene_coordinates(
        genes, args.downstream_exon_extension
    )  # Extend gene coordinates first

    # Step 3: Extend exons (now only extends exons,
    # using extended gene coordinates)
    log_message("Extending terminal exons...")
    genes = extend_exon_downstream(genes, args.downstream_exon_extension)

    # Step 4: Define exons and introns
    log_message("Defining exons and introns...")
    genes = define_exons_introns(genes)

    # Step 5: Construct segments
    log_message("Constructing segments...")
    genes = construct_segments(genes)

    # Step 6: Identify PAS in segments and construct subsegments
    log_message(
        "Identifying PAS sites in segments and constructing subsegments..."
    )
    genes = identify_pas_in_segments(genes, args.pas_atlas)

    # Step 7: Write to TSV (including segments and subsegments)
    log_message("Writing genes, segments and subsegments to TSV...")
    gene_mapping = write_segments_pas_to_tsv(
        genes,
        output_genes_tsv,
        output_segments_tsv,
        output_subsegments_tsv,
    )
    # Step 8: Write to GTF (including segments and subsegments)
    log_message("Writing genes, segments and subsegments to GTF...")
    write_segments_pas_to_gtf(genes, output_gtf)

    # Step 9: Calculate mean coverage for subsegments
    log_message("Calculating mean coverage for subsegments...")

    subsegments_df = create_subsegments_dataframe(genes, gene_mapping)

    coverage_results = calculate_mean_coverage(
        subsegments_df, args.coverage[0], args.coverage[1]
    )  # Pass DataFrame

    # Step 10: Write coverage results to TSV
    log_message("Writing coverage results to TSV...")
    write_coverage_results(
        coverage_results, output_coverage_tsv
    )  # Write to the new TSV

    log_message(
        "Genomic segment construction, PAS identification,"
        "subsegment construction and coverage calculation pipeline completed."
    )


if __name__ == "__main__":
    main()
