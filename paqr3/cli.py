import argparse
from datetime import datetime
from paqr3.version import __version__
from paqr3.filter_annotation import filter_annotation_by_gene_type, resolve_gene_overlaps
from paqr3.construct_segments import extend_gene_coordinates, extend_exon_downstream, define_exons_introns, construct_segments
from paqr3.detect_pas import identify_pas_in_segments, write_segments_pas_to_gtf
from paqr3.calculate_mean_coverage import calculate_mean_coverage, write_coverage_results
from paqr3.models import Gene, Transcript, Region

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
            gene_attributes = {k: v for k, v in row.items() if k not in ["seqname", "source", "feature", "start", "end", "score", "strand"]}
            gene = Gene(gene_id, attributes=gene_attributes)
            genes[gene_id] = gene
        elif row["feature"] == "transcript":
            gene_id = row["gene_id"]
            transcript_id = row["transcript_id"]
            transcript_attributes = {k: v for k, v in row.items() if k not in ["seqname", "source", "feature", "start", "end", "score", "strand", "gene_id"]}

            if gene_id in genes:
                transcript = Transcript(transcript_id, row["strand"], attributes=transcript_attributes)
                genes[gene_id].add_transcript(transcript)
        elif row["feature"] == "exon":
            gene_id = row["gene_id"]
            transcript_id = row["transcript_id"]
            exon_attributes = {k: v for k, v in row.items() if k not in ["seqname", "source", "feature", "start", "end", "score", "strand", "gene_id", "transcript_id"]}
            if gene_id in genes and transcript_id in genes[gene_id].transcripts:
                exon = Region(
                    region_type="exon",
                    chrom=row["seqname"],
                    start=row["start"],
                    end=row["end"],
                    strand=row["strand"],
                    attributes=exon_attributes
                )
                genes[gene_id].transcripts[transcript_id].add_region(exon)

    return genes

def main():
    parser = argparse.ArgumentParser(description="Filter and construct genomic segments from RNA-Seq annotations.")
    parser.add_argument("--annotation", "-a", type=str, required=True, help="Path to the annotation GTF file.")
    parser.add_argument("--strandedness", type=str, default="true", choices=["true", "false"], help="Whether the data is stranded (true) or unstranded (false).")
    parser.add_argument("--downstream_exon_extension", type=int, default=200, help="Number of bases to extend terminal exons.")
    parser.add_argument("--pas_atlas", "-pa", type=str, required=True, help="Path to the PAS atlas BED file.")
    parser.add_argument("--output_regions", "-or", type=str, required=True, help="Path to the output GTF file for regions and segments.")
    parser.add_argument("--coverage", "-c", type=str, required=True, help="Path to the coverage BED file.")
    parser.add_argument("--output_coverage", "-oc", type=str, required=True, help="Path to output TSV file for coverage results.")
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=f"PAQR3 v{__version__}, (c) 2025 by Zavolab (zavolab-biozentrum@unibas.ch)",
        help="show version information and exit",
    )


    args = parser.parse_args()

    log_message("Starting genomic segment construction pipeline...")

    # Step 1: Filter annotation
    filtered_gtf = filter_annotation_by_gene_type(args.annotation)
    stranded = args.strandedness.lower() == "true"
    log_message(f"Resolving gene overlaps (strandedness: {stranded})...")
    filtered_gtf = resolve_gene_overlaps(filtered_gtf, strandedness=stranded)

    # Step 2: Parse GTF to Gene/Transcript/Region objects
    log_message("Parsing GTF to Gene/Transcript/Region objects...")
    genes = parse_gtf_to_genes(filtered_gtf)

    # Step 2.5: Extend gene coordinates
    log_message("Extending gene coordinates...")
    genes = extend_gene_coordinates(genes, args.downstream_exon_extension) # Extend gene coordinates first

    # Step 3: Extend exons (now only extends exons, using extended gene coordinates)
    log_message("Extending terminal exons...")
    genes = extend_exon_downstream(genes, args.downstream_exon_extension)

    # Step 4: Define exons and introns
    log_message("Defining exons and introns...")
    genes = define_exons_introns(genes)

    # Step 5: Construct segments
    log_message("Constructing segments...")
    genes = construct_segments(genes)

    # Step 6: Identify PAS in segments
    log_message("Identifying PAS sites in segments...")
    genes = identify_pas_in_segments(genes, args.pas_atlas)

    # Step 7: Write to GTF (including PAS info)
    log_message("Writing genes, segments and overlapping PAS to GTF...")
    write_segments_pas_to_gtf(genes, args.output_regions)

    # Step 8: Calculate mean coverage
    log_message("Calculating mean coverage...")
    coverage_results = calculate_mean_coverage(genes, args.coverage)

    # Step 9: Write coverage results
    log_message("Writing coverage results...")
    write_coverage_results(coverage_results, args.output_coverage)

    log_message("Genomic segment construction, PAS identification and coverage calculation pipeline completed.")


if __name__ == "__main__":
    main()
