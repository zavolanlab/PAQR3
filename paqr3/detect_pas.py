from datetime import datetime
from pybedtools import BedTool  # type: ignore
from intervaltree import Interval, IntervalTree  # type: ignore
from paqr3.models import Region


def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


def identify_pas_in_segments(genes, pas_atlas_bed):
    """
    Identifies PAS sites overlapping with segments
    and constructs subsegments.
    """

    log_message("Reading PAS atlas...")
    pas_bed = BedTool(pas_atlas_bed)

    log_message("Building PAS interval tree...")
    pas_trees = {}

    for pas in pas_bed:
        chrom = pas.chrom
        strand = pas.strand
        if (chrom, strand) not in pas_trees:
            pas_trees[(chrom, strand)] = IntervalTree()
        pas_trees[(chrom, strand)].add(
            Interval(int(pas.start), int(pas.end), pas.fields)
        )

    log_message(
        "Identifying PAS overlaps with segments"
        "and constructing subsegments..."
    )

    for gene in genes.values():
        for segment in gene.segments:
            chrom = segment.chrom
            strand = segment.strand
            segment_start = segment.start
            segment_end = segment.end

            overlapping_pas = []
            subsegments = []
            if (chrom, strand) in pas_trees:
                overlapping_intervals = pas_trees[(chrom, strand)][
                    segment_start:segment_end
                ]
                sorted_intervals = sorted(
                    overlapping_intervals, key=lambda x: x.begin
                )  # Sort by start position

                current_subsegment_start = segment_start

                for interval in sorted_intervals:
                    pas_start, pas_end, pas_fields = (
                        interval.begin,
                        interval.end,
                        interval.data,
                    )
                    pas_id = pas_fields[3]

                    # Create subsegment up to PAS start
                    if pas_start > current_subsegment_start:
                        subsegments.append(
                            Region(
                                region_type="subsegment",
                                chrom=chrom,
                                start=current_subsegment_start,
                                end=pas_start - 1,
                                strand=strand,
                                attributes={
                                    "gene_id": gene.gene_id,
                                    "segment_id": segment.attributes[
                                        "segment_number"
                                    ],
                                    "subsegment_number": len(subsegments) + 1,
                                    "strand": strand,
                                },  # Added strand
                            )
                        )

                    # Create subsegment for PAS region
                    subsegments.append(
                        Region(
                            region_type="subsegment",
                            chrom=chrom,
                            start=pas_start,
                            end=pas_end,
                            strand=strand,
                            attributes={
                                "gene_id": gene.gene_id,
                                "segment_id": segment.attributes[
                                    "segment_number"
                                ],
                                "subsegment_number": len(subsegments) + 1,
                                "strand": strand,
                            },  # Added strand
                        )
                    )

                    current_subsegment_start = (
                        pas_end + 1
                    )  # Start next subsegment after PAS end
                    overlapping_pas.append(pas_id)

                # Create any remaining subsegment after the last PAS
                if current_subsegment_start <= segment_end:
                    subsegments.append(
                        Region(
                            region_type="subsegment",
                            chrom=chrom,
                            start=current_subsegment_start,
                            end=segment_end,
                            strand=strand,
                            attributes={
                                "gene_id": gene.gene_id,
                                "segment_id": segment.attributes[
                                    "segment_number"
                                ],
                                "subsegment_number": len(subsegments) + 1,
                                "strand": strand,
                            },  # Added strand
                        )
                    )

            segment.attributes["overlapping_pas"] = overlapping_pas
            segment.subsegments = subsegments  # Add subsegments to the segment

    return genes


def write_segments_pas_to_gtf(genes, output_file):
    """Writes genes and segments to a GTF file."""

    with open(output_file, "w") as f:
        for gene in genes.values():
            for segment in gene.segments:
                segment_line = segment.to_gtf_format()
                if "overlapping_pas" in segment.attributes:
                    pas_list = ",".join(segment.attributes["overlapping_pas"])
                    segment_line = segment_line.replace(
                        ";", f'; overlapping_pas "{pas_list}";'
                    )
                f.write(segment_line + "\n")

    log_message(f"Genes and segments written to {output_file}")
