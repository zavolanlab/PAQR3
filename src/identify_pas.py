import argparse
from datetime import datetime
from pybedtools import BedTool
from intervaltree import Interval, IntervalTree

def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def identify_pas_in_segments(genes, pas_atlas_bed):
    """Identifies PAS sites overlapping with segments."""

    log_message("Reading PAS atlas...")
    pas_bed = BedTool(pas_atlas_bed)

    log_message("Identifying PAS overlaps with segments...")

    for gene in genes.values():
        for segment in gene.segments:
            segment_bed = BedTool(
                [(segment.chrom, segment.start, segment.end, gene.gene_id, segment.attributes.get("segment_number"), segment.strand)]
            )
            overlaps = segment_bed.intersect(pas_bed, wa=True, wb=True)
            overlapping_pas = []

            for overlap in overlaps:
                seg_chrom, seg_start, seg_end, gene_id, seg_num, seg_strand = overlap[:6]
                pas_chrom, pas_start, pas_end, pas_id, pas_score, pas_strand = overlap[6:]

                # Check for boundary overlaps
                if (int(seg_start) < int(pas_start) < int(seg_end)) or (int(seg_start) < int(pas_end) < int(seg_end)):
                    overlapping_pas.append(pas_id)

            segment.attributes["overlapping_pas"] = overlapping_pas  # Store PAS IDs in segment attributes

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