from datetime import datetime
from pybedtools import BedTool
from intervaltree import Interval, IntervalTree


def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


def identify_pas_in_segments(genes, pas_atlas_bed, strandedness=True):
    """Identifies PAS sites overlapping with segments using an interval tree."""

    log_message("Reading PAS atlas...")
    pas_bed = BedTool(pas_atlas_bed)

    log_message("Building PAS interval tree...")
    pas_trees = {}  # Dictionary to store interval trees per chromosome and strand

    for pas in pas_bed:
        chrom = pas.chrom
        strand = pas.strand if strandedness else "+" # Use dummy strand if unstranded
        if (chrom, strand) not in pas_trees:
            pas_trees[(chrom, strand)] = IntervalTree()
        pas_trees[(chrom, strand)].add(Interval(int(pas.start), int(pas.end), pas.fields))

    log_message("Identifying PAS overlaps with segments...")

    for gene in genes.values():
        for segment in gene.segments:
            chrom = segment.chrom
            strand = segment.strand
            segment_start = segment.start
            segment_end = segment.end

            overlapping_pas = []
            if (chrom, strand) in pas_trees:
                for interval in pas_trees[(chrom, strand)][segment_start:segment_end]:
                    pas_start, pas_end, pas_fields = interval.begin, interval.end, interval.data
                    pas_id = pas_fields[3] # Assuming PAS ID is in the 4th field

                    # Check for boundary overlaps (same as before)
                    if (int(segment_start) < int(pas_start) < int(segment_end)) or (int(segment_start) < int(pas_end) < int(segment_end)):
                        overlapping_pas.append(pas_id)

            segment.attributes["overlapping_pas"] = overlapping_pas

    return genes


def write_segments_pas_to_gtf(genes, output_file):
    """Writes genes, segments and overlapping PAS to a GTF file."""

    with open(output_file, "w") as f:
        for gene in genes.values():
            f.write(gene.to_gtf_format() + "\n")  # Write gene information
            for segment in gene.segments:
                segment_line = segment.to_gtf_format()
                if "overlapping_pas" in segment.attributes:
                    pas_list = ",".join(segment.attributes["overlapping_pas"])
                    segment_line = segment_line.replace(";", f'; overlapping_pas "{pas_list}";')
                f.write(segment_line + "\n")  # Write segment and PAS information

    log_message(f"Genes, segments and overlapping PAS written to {output_file}")