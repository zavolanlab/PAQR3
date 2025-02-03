
import argparse
from datetime import datetime
from filter_gtf import filter_gtf_by_gene_type, resolve_strand_overlaps, write_filtered_gtf

# Utility function for timestamped messages
def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def main():
    parser = argparse.ArgumentParser(description="Filter and quantify PAS usage from RNA-Seq alignments.")
    parser.add_argument(
        "--annotation",
        type=str,
        required=True,
        help="Path to the annotation GTF file."
    )
    parser.add_argument(
        "--pas_atlas",
        type=str,
        required=True,
        help="Path to the PAS atlas BED file."
    )
    parser.add_argument(
        "--coverage",
        type=str,
        required=False,
        help="Path to the RNA-Seq coverage file (e.g., from samtools coverage)."
    )
    parser.add_argument(
        "--output_regions",
        type=str,
        required=False,
        help="Path to the output filtered GTF file (optional)."
    )
    parser.add_argument(
        "--output_pas_usage",
        type=str,
        required=False,
        help="Path to the output file for PAS usage (optional)."
    )
    parser.add_argument(
        "--filter_strand_overlaps",
        action="store_true",
        help="Enable filtering of overlapping genes based on strand."
    )

    args = parser.parse_args()

    log_message("Starting GTF filtering pipeline...")

    # Step 1: Filter the GTF file by gene type
    filtered_gtf = filter_gtf_by_gene_type(args.annotation)

    # Step 2: Resolve strand overlaps (optional)
    if args.filter_strand_overlaps:
        log_message("Resolving strand overlaps...")
        filtered_gtf = resolve_strand_overlaps(filtered_gtf)

    # Step 3: Write the filtered GTF to a file (optional)
    if args.output_regions:
        write_filtered_gtf(filtered_gtf, args.output_regions)

    log_message("GTF filtering pipeline completed.")

if __name__ == "__main__":
    main()