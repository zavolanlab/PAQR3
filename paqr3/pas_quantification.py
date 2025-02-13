import argparse
import pybedtools
from gtfparse import read_gtf

# Define a class to hold information about each region in the transcript
# Add a string representation for debugging
class Region:
    def __init__(self, region_type, chrom, start, end, strand, exon_number=None, intron_number=None, attributes=None):
        self.region_type = region_type  # exon or intron
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.exon_number = exon_number  # For exons
        self.intron_number = intron_number
        self.attributes = attributes or {}  # GTF attributes

    def to_gtf_format(self):
        """Convert the Region object to a GTF format string."""
        attr_str = "; ".join([f'{key} "{value}"' for key, value in self.attributes.items()]) + ";"
        return f"{self.chrom}\t{self.attributes.get('source', 'unknown')}\t{self.region_type}\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{attr_str}"

    def __str__(self):
        """String representation for debugging."""
        return f"Region({self.region_type}, {self.chrom}:{self.start}-{self.end}, {self.strand}, attributes={self.attributes})"

    def __repr__(self):
        return self.__str__()
        
# Define a class to hold transcript information
class Transcript:
    def __init__(self, transcript_id, strand, attributes=None):
        self.transcript_id = transcript_id
        self.strand = strand
        self.regions = []
        self.attributes = attributes or {}  # GTF attributes

    def add_region(self, region):
        self.regions.append(region)

    def to_gtf_format(self):
        """Convert the Transcript object to a GTF format string."""
        attr_str = "; ".join([f'{key} "{value}"' for key, value in self.attributes.items()]) + ";"
        start = min(region.start for region in self.regions)
        end = max(region.end for region in self.regions)
        return f"{self.regions[0].chrom}\t{self.attributes.get('source', 'unknown')}\ttranscript\t{start}\t{end}\t.\t{self.strand}\t.\t{attr_str}"

    def __repr__(self):
        return f"Transcript: {self.transcript_id}, Regions: {self.regions}"


# Define a class to hold gene information and associated transcripts
class Gene:
    def __init__(self, gene_id, attributes=None):
        self.gene_id = gene_id
        self.transcripts = {}
        self.attributes = attributes or {}  # GTF attributes

    def add_transcript(self, transcript):
        self.transcripts[transcript.transcript_id] = transcript

    def to_gtf_format(self):
        """Convert the Gene object to a GTF format string."""
        attr_str = "; ".join([f'{key} "{value}"' for key, value in self.attributes.items()]) + ";"
        start = min(region.start for transcript in self.transcripts.values() for region in transcript.regions)
        end = max(region.end for transcript in self.transcripts.values() for region in transcript.regions)
        strand = next(iter(self.transcripts.values())).strand  # Use strand of the first transcript
        chrom = next(iter(self.transcripts.values())).regions[0].chrom  # Use chrom of the first region
        return f"{chrom}\t{self.attributes.get('source', 'unknown')}\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}"

    def __repr__(self):
        return f"Gene: {self.gene_id}, Transcripts: {self.transcripts}"
    
# Function to read the GTF file and create gene and transcript objects with exon regions
def read_gtf_file(gtf_file):
    # Read GTF file into a pandas dataframe
    gtf_data = read_gtf(gtf_file, result_type="pandas")

    # Filter for exons only
    exons = gtf_data[gtf_data["feature"] == "exon"]

    # Create a dictionary to hold genes
    genes = {}

    # Group exons by gene and transcript
    for _, row in exons.iterrows():
        gene_id = row["gene_id"]
        transcript_id = row["transcript_id"]
        exon_number = int(row.get("exon_number", -1))  # Default to -1 if missing
        chrom = row["seqname"]  # Extract chromosome name
        start = row["start"]
        end = row["end"]
        strand = row["strand"]

        # Extract attributes for genes, transcripts, and exons
        gene_attributes = {key: row[key] for key in row.index if "gene" in key}
        transcript_attributes = {key: row[key] for key in row.index if "transcript" in key}
        exon_attributes = {key: row[key] for key in row.index if "exon" in key}
        exon_attributes.update(gene_attributes)
        exon_attributes.update(transcript_attributes)

        # Get or create the gene
        if gene_id not in genes:
            genes[gene_id] = Gene(gene_id, gene_attributes)
        gene = genes[gene_id]

        # Get or create the transcript
        if transcript_id not in gene.transcripts:
            transcript = Transcript(transcript_id, strand, transcript_attributes)
            gene.add_transcript(transcript)
        else:
            transcript = gene.transcripts[transcript_id]

        # Create the exon region and add to the transcript
        exon_region = Region("exon", chrom, start, end, strand, exon_number, attributes=exon_attributes)
        transcript.add_region(exon_region)

    return genes

def add_intronic_regions(genes):
    """
    Add intronic regions to the existing annotation.

    Args:
        genes (dict): A dictionary of Gene objects.

    Returns:
        dict: Updated dictionary of Gene objects with intronic regions added.
    """
    for gene in genes.values():
        for transcript in gene.transcripts.values():
            # Sort exons by their start positions
            transcript.regions.sort(key=lambda r: r.start)

            # Add intronic regions between exons
            reconstructed_regions = []
            intron_count = 0
            for i, region in enumerate(transcript.regions):
                reconstructed_regions.append(region)

                # Add intron if not the last exon
                if i < len(transcript.regions) - 1:
                    next_region = transcript.regions[i + 1]
                    intron_start = region.end + 1
                    intron_end = next_region.start - 1
                    if intron_start < intron_end:
                        intron_count += 1
                        intron_region = Region(
                            region_type="intron",
                            chrom=region.chrom,
                            start=intron_start,
                            end=intron_end,
                            strand=region.strand,
                            intron_number=intron_count,
                            attributes={"gene_id": gene.gene_id, "transcript_id": transcript.transcript_id}
                        )
                        reconstructed_regions.append(intron_region)

            # Update transcript regions
            transcript.regions = reconstructed_regions

    return genes

def construct_contiguous_regions(genes):
    """
    Construct contiguous, non-overlapping regions within each gene based on exon and intron boundaries.

    Args:
        genes (dict): A dictionary of Gene objects with exon and intron annotations.

    Returns:
        dict: Updated dictionary of Gene objects with constructed regions.
    """
    for gene in genes.values():
        # Collect all unique boundaries (start and end) across all transcripts
        boundaries = set()
        for transcript in gene.transcripts.values():
            for region in transcript.regions:
                boundaries.add(region.start)
                boundaries.add(region.end + 1)  # Ensure the end is exclusive for the next region

        # Sort boundaries
        sorted_boundaries = sorted(boundaries)

        # Construct contiguous, non-overlapping regions
        regions = []
        for i in range(len(sorted_boundaries) - 1):
            start = sorted_boundaries[i]
            end = sorted_boundaries[i + 1] - 1  # Adjust to make the region inclusive
            if start <= end:  # Only include valid regions
                # Ensure there are no gaps between regions
                if regions and start != regions[-1][1] + 1:
                    start = regions[-1][1] + 1  # Adjust to make regions contiguous
                regions.append((start, end))

        # Replace transcript regions with the newly constructed regions
        gene_regions = []
        for i, (start, end) in enumerate(regions, start=1):
            # Create a Region object for each new region
            new_region = Region(
                region_type="gene_fragment",
                chrom=gene.transcripts[next(iter(gene.transcripts))].regions[0].chrom,  # Use chrom of the first transcript
                start=start,
                end=end,
                strand=gene.transcripts[next(iter(gene.transcripts))].strand,  # Use strand of the first transcript
                attributes={"gene_id": gene.gene_id, "region_id": i},
            )
            gene_regions.append(new_region)

        # Assign regions at the gene level (regions are shared by all transcripts of the gene)
        for transcript in gene.transcripts.values():
            transcript.regions = gene_regions

    return genes

# Add a function to write the regions to a GTF file
def write_regions_to_gtf(genes, output_file):
    """
    Write the gene and region information to a GTF file.

    Args:
        genes (dict): A dictionary of Gene objects with constructed regions.
        output_file (str): Path to the output GTF file.
    """
    with open(output_file, 'w') as f:
        for gene in genes.values():
            # Write gene-level GTF entry
            f.write(gene.to_gtf_format() + '\n')

            # Write region-level GTF entries
            unique_regions = set()  # Ensure no duplicate regions
            for transcript in gene.transcripts.values():
                for region in transcript.regions:
                    region_key = (region.chrom, region.start, region.end, region.strand)
                    if region_key not in unique_regions:
                        unique_regions.add(region_key)
                        f.write(region.to_gtf_format() + '\n')

    print(f"Regions written to {output_file}")


def filter_regions_with_pas_and_track_overlaps(genes, pas_bed_file):
    """
    Filter regions based on overlap with PAS (polyadenylation sites) from a BED file and track overlaps.

    Args:
        genes (dict): A dictionary of Gene objects with constructed regions.
        pas_bed_file (str): Path to the BED file containing PAS sites.

    Returns:
        dict: A dictionary where keys are Region objects with at least one overlapping PAS,
              and values are lists of overlapping PAS IDs.
    """
    # Load PAS sites from the BED file
    pas_bed = pybedtools.BedTool(pas_bed_file)

    region_pas_map = {}

    for gene in genes.values():
        for transcript in gene.transcripts.values():
            regions = transcript.regions

            for i, region in enumerate(regions):
                # Skip first region (positive strand) or last region (negative strand)
                if (transcript.strand == "+" and i == 0) or (transcript.strand == "-" and i == len(regions) - 1):
                    continue

                # Convert region to a BedTool-compatible format
                region_bed = pybedtools.create_interval_from_list([
                    region.chrom,
                    str(region.start - 1),  # BED is 0-based
                    str(region.end),       # BED is 1-based
                    region.attributes["region_id"],
                    "0",                   # Dummy score
                    region.strand
                ])

                # Check for overlaps with PAS
                overlaps = pas_bed.intersect(pybedtools.BedTool([region_bed]), wa=True)

                valid_pas_ids = set()  # Use a set to avoid duplicate PAS IDs
                for pas_interval in overlaps:
                    pas_start = pas_interval.start + 1  # BED interval is 0-based, GTF is 1-based
                    pas_end = pas_interval.end
                    pas_id = pas_interval[3]  # 4th column of the BED file

                    # Exclude PAS overlapping boundaries or the entire region
                    if pas_start <= region.start or pas_end >= region.end:
                        continue  # Skip this PAS
                    else:
                        valid_pas_ids.add(pas_id)

                if valid_pas_ids:
                    # Add region and its PAS IDs to the dictionary
                    region_pas_map[region] = list(valid_pas_ids)  # Convert set to list

    return region_pas_map

# Add an argument parser for the script
def main():
    parser = argparse.ArgumentParser(description="Quantify PAS usage from RNA-Seq alignments.")
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
        required=True,
        help="Path to the output GTF file for regions."
    )
    parser.add_argument(
        "--output_pas_usage",
        type=str,
        required=False,
        help="Path to the output file for PAS usage (optional)."
    )

    args = parser.parse_args()

    # Step 1: Load the genes and regions
    print("Reading the annotation GTF file...")
    genes = read_gtf_file(args.annotation)
    genes_with_introns = add_intronic_regions(genes)
    genes_with_regions = construct_contiguous_regions(genes_with_introns)

    # Step 2: Filter regions based on PAS overlaps
    print("Filtering regions with PAS overlaps...")
    region_pas_map = filter_regions_with_pas_and_track_overlaps(genes_with_regions, args.pas_atlas)

    # Step 3: Write regions to GTF for visualization
    print("Writing regions to output GTF...")
    write_regions_to_gtf(genes_with_regions, args.output_regions)

    # Step 4 (optional): Quantify PAS usage (if coverage file is provided)
    if args.coverage:
        print("Calculating PAS usage...")
        usage_results = calculate_usage(region_pas_map, args.coverage)
        if args.output_pas_usage:
            print(f"Writing PAS usage results to {args.output_pas_usage}...")
            with open(args.output_pas_usage, 'w') as f:
                for region, metrics in usage_results.items():
                    f.write(f"{region}: {metrics}\n")

        # Display a sample of the results
        print("Sample PAS usage results:")
        for region, metrics in list(usage_results.items())[:5]:
            print(f"{region}: {metrics}")

if __name__ == "__main__":
    main()