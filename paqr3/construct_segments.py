import logging
import pandas as pd  # type: ignore
from pybedtools import BedTool  # type: ignore
from intervaltree import Interval, IntervalTree  # type: ignore
from paqr3.models import Region, Gene, Transcript
from gtfparse import read_gtf  # type: ignore


def log_message(message):
    logging.info(message)


class ConstructSegments:
    def __init__(
        self, annotation_file, pas_atlas_file, downstream_exon_extension
    ):
        self.annotation_file = annotation_file
        self.pas_atlas_file = pas_atlas_file
        self.downstream_exon_extension = downstream_exon_extension
        self.genes = {}
        self.gene_mapping = {}

    def parse_gtf_to_genes(self, gtf_data):
        """Parses GTF data (pandas DataFrame) into a dictionary of Gene objects."""
        genes = {}
        for _, row in gtf_data.iterrows():
            if row["feature"] == "gene":
                gene_id = row["gene_id"]
                if pd.isna(row["seqname"]):  # Check if chromosome name is NaN
                    continue  # Skip genes with missing chromosome names

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

    def extend_gene_coordinates(self):
        """Extends gene coordinates (correctly handles negative strand)."""

        sorted_genes = sorted(
            self.genes.values(),
            key=lambda g: (
                next(
                    (
                        r.chrom
                        for t in g.transcripts.values()
                        for r in t.regions
                    ),
                    None,
                ),
                (
                    min(
                        list(
                            r.start
                            for t in g.transcripts.values()
                            for r in t.regions
                        )
                    )
                    if g.transcripts
                    else float("inf")
                ),
            ),
        )

        for i, gene in enumerate(sorted_genes):
            gene_strand = next(
                (t.strand for t in gene.transcripts.values()), None
            )
            if gene_strand is None:
                continue

            gene_start = (
                min(
                    list(
                        r.start
                        for t in gene.transcripts.values()
                        for r in t.regions
                    )
                )
                if gene.transcripts
                else float("inf")
            )
            gene_end = (
                max(
                    list(
                        r.end
                        for t in gene.transcripts.values()
                        for r in t.regions
                    )
                )
                if gene.transcripts
                else float("-inf")
            )

            if gene_strand == "+":
                potential_end = gene_end + self.downstream_exon_extension
                next_gene_start = potential_end

                for other_gene in sorted_genes[i + 1 :]:
                    other_gene_strand = next(
                        (t.strand for t in other_gene.transcripts.values()),
                        None,
                    )
                    if (
                        next(
                            (
                                r.chrom
                                for t in other_gene.transcripts.values()
                                for r in t.regions
                            ),
                            None,
                        )
                        == next(
                            (
                                r.chrom
                                for t in gene.transcripts.values()
                                for r in t.regions
                            ),
                            None,
                        )
                        and other_gene_strand == gene_strand
                        and min(
                            list(
                                r.start
                                for t in other_gene.transcripts.values()
                                for r in t.regions
                            )
                        )
                        > gene_end
                    ):
                        next_gene_start = min(
                            potential_end,
                            min(
                                list(
                                    r.start
                                    for t in other_gene.transcripts.values()
                                    for r in t.regions
                                )
                            )
                            - 1,
                        )
                        break

                gene_end = min(potential_end, next_gene_start)

            elif gene_strand == "-":
                potential_start = gene_start - self.downstream_exon_extension
                previous_gene_end = potential_start

                for other_gene in reversed(sorted_genes[:i]):
                    other_gene_strand = next(
                        (t.strand for t in other_gene.transcripts.values()),
                        None,
                    )
                    if (
                        next(
                            (
                                r.chrom
                                for t in other_gene.transcripts.values()
                                for r in t.regions
                            ),
                            None,
                        )
                        == next(
                            (
                                r.chrom
                                for t in gene.transcripts.values()
                                for r in t.regions
                            ),
                            None,
                        )
                        and other_gene_strand == gene_strand
                        and max(
                            list(
                                r.end
                                for t in other_gene.transcripts.values()
                                for r in t.regions
                            )
                        )
                        < gene_start
                    ):
                        previous_gene_end = max(
                            potential_start,
                            max(
                                list(
                                    r.end
                                    for t in other_gene.transcripts.values()
                                    for r in t.regions
                                )
                            )
                            + 1,
                        )
                        break

                gene_start = max(potential_start, previous_gene_end)

            gene.attributes["start"] = gene_start
            gene.attributes["end"] = gene_end

    def extend_exon_downstream(self):
        """Extends terminal exons (CORRECTLY handles negative strand)."""

        for gene in self.genes.values():
            gene_strand = next(
                (t.strand for t in gene.transcripts.values()), None
            )
            if gene_strand is None:
                continue

            for transcript in gene.transcripts.values():
                if transcript.regions:

                    if (
                        gene_strand == "+"
                    ):  # Positive strand - same logic as before
                        terminal_exon = max(
                            transcript.regions,
                            key=lambda r: (
                                r.end
                                if r.region_type == "exon" and r.strand == "+"
                                else (
                                    r.start
                                    if r.region_type == "exon"
                                    and r.strand == "-"
                                    else -1
                                )
                            ),
                        )
                        if terminal_exon.strand == "+":
                            terminal_exon.end = min(
                                terminal_exon.end
                                + self.downstream_exon_extension,
                                gene.attributes["end"],
                            )
                        elif terminal_exon.strand == "-":
                            terminal_exon.start = max(
                                terminal_exon.start
                                - self.downstream_exon_extension,
                                gene.attributes["start"],
                            )

                    elif gene_strand == "-":  # Negative strand
                        # Find the exon with the lowest exon number
                        terminal_exon = min(
                            transcript.regions,
                            key=lambda r: (
                                r.start  # Use start for negative strand sorting
                            ),
                        )

                        if terminal_exon.strand == "-":
                            terminal_exon.start = max(
                                terminal_exon.start
                                - self.downstream_exon_extension,
                                gene.attributes["start"],
                            )  # Extend the START on the negative strand
                        elif terminal_exon.strand == "+":
                            terminal_exon.end = min(
                                terminal_exon.end
                                + self.downstream_exon_extension,
                                gene.attributes["end"],
                            )

    def define_exons_introns(self):
        """Constructs exons and introns for all transcripts."""

        for gene in self.genes.values():
            for transcript in gene.transcripts.values():
                transcript.regions.sort(
                    key=lambda r: r.start
                )  # Ensure regions are sorted

                exons = [
                    r for r in transcript.regions if r.region_type == "exon"
                ]
                introns = []

                for i in range(len(exons) - 1):
                    exon1 = exons[i]
                    exon2 = exons[i + 1]
                    intron_start = exon1.end + 1
                    intron_end = exon2.start - 1

                    if (
                        intron_start <= intron_end
                    ):  # Only create if it's a valid intron
                        intron = Region(
                            region_type="intron",
                            chrom=exon1.chrom,
                            start=intron_start,
                            end=intron_end,
                            strand=exon1.strand,
                            intron_number=i
                            + 1,  # Intron numbering starts from 1
                            attributes={
                                "gene_id": gene.gene_id,
                                "transcript_id": transcript.transcript_id,
                            },
                        )
                        introns.append(intron)

                transcript.regions = (
                    exons + introns
                )  # Combine exons and introns
                transcript.regions.sort(
                    key=lambda r: r.start
                )  # Sort all regions

    def construct_segments(self):
        """Constructs segments."""

        for gene in self.genes.values():
            strand = next((t.strand for t in gene.transcripts.values()), None)
            if strand is None:
                continue

            all_regions = []

            for transcript in gene.transcripts.values():
                for region in transcript.regions:
                    all_regions.append((region.start, region.end, region))

            all_regions.sort(key=lambda x: (x[0], x[1]))

            segment_boundaries = set()

            for start, end, _ in all_regions:
                segment_boundaries.add(start)
                segment_boundaries.add(end + 1)

            gene_start = (
                min(
                    list(
                        r.start
                        for t in gene.transcripts.values()
                        for r in t.regions
                    )
                )
                if gene.transcripts
                else float("inf")
            )
            gene_end = (
                max(
                    list(
                        r.end
                        for t in gene.transcripts.values()
                        for r in t.regions
                    )
                )
                if gene.transcripts
                else float("-inf")
            )

            terminal_exons_ends = set()
            terminal_exons_starts = set()

            for transcript in gene.transcripts.values():
                for region in transcript.regions:
                    if region.region_type == "exon":
                        is_terminal = (
                            region.strand == "+"
                            and region.end
                            == max(
                                r.end
                                for r in transcript.regions
                                if r.region_type == "exon"
                            )
                        ) or (
                            region.strand == "-"
                            and region.start
                            == min(
                                r.start
                                for r in transcript.regions
                                if r.region_type == "exon"
                            )
                        )
                        if is_terminal:
                            if region.strand == "+":
                                terminal_exons_ends.add(region.end)
                            else:
                                terminal_exons_starts.add(region.start)

            # Remove invalid terminal exon boundaries
            valid_segment_boundaries = set()
            for boundary in segment_boundaries:
                if strand == "+":
                    if (
                        boundary - 1 in terminal_exons_ends
                        and boundary - 1 != gene_end
                    ):
                        continue  # Skip this boundary
                    else:
                        valid_segment_boundaries.add(boundary)
                elif strand == "-":
                    if (
                        boundary in terminal_exons_starts
                        and boundary != gene_start
                    ):
                        continue  # Skip this boundary
                    else:
                        valid_segment_boundaries.add(boundary)

            valid_segment_boundaries = sorted(list(valid_segment_boundaries))
            segments = []
            for i in range(len(valid_segment_boundaries) - 1):
                start = valid_segment_boundaries[i]
                end = valid_segment_boundaries[i + 1] - 1

                # Prevent zero-length segments
                if start < end:  # Only create segments if start < end
                    segments.append(
                        Region(
                            region_type="segment",
                            chrom=all_regions[0][2].chrom,
                            start=start,
                            end=end,
                            strand=strand,
                            attributes={
                                "gene_id": gene.gene_id,
                                "segment_strand": strand,
                            },
                        )
                    )

            if strand == "-":
                segments.reverse()

            for i, segment in enumerate(segments):
                segment.attributes["segment_number"] = i + 1

            gene.segments = segments

    def identify_pas_in_segments(self):
        """
        Identifies PAS sites overlapping with segments
        and constructs subsegments.
        """

        log_message("Reading PAS atlas...")
        pas_bed = BedTool(self.pas_atlas_file)

        log_message("Building PAS interval tree...")
        pas_trees = {}

        # Create PAS mapping and data in identify_pas_in_segments
        self.pas_mapping = {}
        self.pas_data = []  # Corrected line
        unique_pas_id_counter = 1

        for pas in pas_bed:
            pas_id = f"{pas.chrom}:{pas.start}:{pas.end}:{pas.strand}"
            self.pas_mapping[pas_id] = unique_pas_id_counter
            self.pas_data.append(
                [
                    unique_pas_id_counter,
                    pas.chrom,
                    pas.start,
                    pas.end,
                    pas.strand,
                    pas.name,
                ]
            )
            unique_pas_id_counter += 1

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

        for gene in self.genes.values():
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

                    if (
                        sorted_intervals
                    ):  # Only proceed if there are overlapping intervals.
                        for interval in sorted_intervals:
                            pas_start, pas_end, pas_fields = (
                                interval.begin,
                                interval.end,
                                interval.data,
                            )
                            pas_id = f"{pas_fields[0]}:{pas_fields[1]}:{pas_fields[2]}:{pas_fields[5]}"
                            unique_pas_id = self.pas_mapping.get(pas_id)

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
                                            "subsegment_number": len(
                                                subsegments
                                            )
                                            + 1,
                                            "strand": strand,
                                        },
                                    )
                                )

                            current_subsegment_start = pas_end + 1
                            overlapping_pas.append(unique_pas_id)

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
                                        "subsegment_number": len(subsegments)
                                        + 1,
                                        "strand": strand,
                                    },
                                )
                            )

                segment.attributes["overlapping_pas"] = overlapping_pas
                segment.subsegments = subsegments

    def write_segments_pas_to_tsv(
        self, out_genes_tsv, out_segments_tsv, out_subsegments_tsv, out_pas_tsv
    ):
        """Writes genes, segments, subsegments, and PAS to TSV files."""

        gene_mapping = {}
        unique_gene_id_counter = 1
        genes_data = []
        segments_data = []
        subsegments_data = []

        pas_df = pd.DataFrame(
            self.pas_data,
            columns=[
                "unique_pas_id",
                "chrom",
                "start",
                "end",
                "strand",
                "pas_id",
            ],
        )
        pas_df.to_csv(out_pas_tsv, sep="\t", index=False)

        for gene in self.genes.values():
            first_transcript = next(iter(gene.transcripts.values()), None)
            chrom = (
                first_transcript.regions[0].chrom
                if first_transcript and first_transcript.regions
                else ""
            )
            strand = first_transcript.strand if first_transcript else ""

            unique_gene_id = unique_gene_id_counter
            gene_mapping[gene.gene_id] = unique_gene_id
            unique_gene_id_counter += 1

            genes_data.append(
                [
                    unique_gene_id,
                    gene.gene_id,
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
                            r.end
                            for t in gene.transcripts.values()
                            for r in t.regions
                        )
                    ),
                    strand,
                ]
            )

            for segment in gene.segments:
                overlapping_pas_ids = []
                if "overlapping_pas" in segment.attributes:
                    for pas_id in segment.attributes["overlapping_pas"]:
                        if pas_id:
                            overlapping_pas_ids.append(str(pas_id))
                overlapping_pas_str = ",".join(overlapping_pas_ids)

                segments_data.append(
                    [
                        unique_gene_id,
                        segment.attributes["segment_number"],
                        segment.chrom,
                        segment.start,
                        segment.end,
                        segment.strand,
                        overlapping_pas_str,
                    ]
                )

                if hasattr(segment, "subsegments"):
                    for subsegment in segment.subsegments:
                        subsegments_data.append(
                            [
                                unique_gene_id,
                                segment.attributes["segment_number"],
                                subsegment.attributes["subsegment_number"],
                                subsegment.chrom,
                                subsegment.start,
                                subsegment.end,
                                subsegment.strand,
                            ]
                        )

        genes_df = pd.DataFrame(
            genes_data,
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
            ],
        )

        genes_df.to_csv(out_genes_tsv, sep="\t", index=False)
        segments_df.to_csv(out_segments_tsv, sep="\t", index=False)
        subsegments_df.to_csv(out_subsegments_tsv, sep="\t", index=False)

        log_message(
            "PAS, genes, segments, and subsegments written to "
            f"{out_pas_tsv}, {out_genes_tsv}, {out_segments_tsv}, and {out_subsegments_tsv}"
        )
        self.gene_mapping = gene_mapping

    def create_subsegments_dataframe(self):
        subsegments_data = []
        for gene in self.genes.values():
            for segment in gene.segments:
                if hasattr(segment, "subsegments"):
                    for subsegment in segment.subsegments:
                        # Filter out subsegments that overlap PAS sites
                        if subsegment.end - subsegment.start > 0:
                            subsegments_data.append(
                                [
                                    self.gene_mapping[gene.gene_id],
                                    segment.attributes["segment_number"],
                                    subsegment.attributes["subsegment_number"],
                                    subsegment.chrom,
                                    subsegment.start,
                                    subsegment.end,
                                    subsegment.strand,
                                ]
                            )
        return pd.DataFrame(
            subsegments_data,
            columns=[
                "gene_unique_id",
                "segment_number",
                "subsegment_number",
                "chrom",
                "start",
                "end",
                "strand",
            ],
        )

    def write_segments_pas_to_gtf(self, output_file):
        """Writes segments and PAS to GTF file."""

        with open(output_file, "w") as f:
            for gene in self.genes.values():
                # Write gene
                if "start" in gene.attributes:
                    gene_start = gene.attributes["start"]
                    gene_end = gene.attributes["end"]
                else:
                    gene_start = min(
                        list(
                            r.start
                            for t in gene.transcripts.values()
                            for r in t.regions
                        )
                    )
                    gene_end = max(
                        list(
                            r.end
                            for t in gene.transcripts.values()
                            for r in t.regions
                        )
                    )
                gene_line = (
                    f"{next(iter(gene.transcripts.values())).regions[0].chrom}\t"
                    f"paqr3\tgene\t{gene_start}\t{gene_end}\t.\t"
                    f"{next(iter(gene.transcripts.values())).strand}\t.\t"
                    f'gene_id "{gene.gene_id}";'
                )
                f.write(gene_line + "\n")

                # Write segments
                for segment in gene.segments:
                    segment_line = (
                        f"{segment.chrom}\tpaqr3\tsegment\t{segment.start}\t"
                        f"{segment.end}\t.\t{segment.strand}\t.\t"
                        f'gene_id "{gene.gene_id}"; segment_number '
                        f'"{segment.attributes["segment_number"]}";'
                    )
                    f.write(segment_line + "\n")

                    # Write subsegments
                    if hasattr(segment, "subsegments"):
                        for subsegment in segment.subsegments:
                            subsegment_line = (
                                f"{subsegment.chrom}\tpaqr3\tsubsegment\t"
                                f"{subsegment.start}\t{subsegment.end}\t.\t"
                                f"{subsegment.strand}\t.\tgene_id "
                                f'"{gene.gene_id}"; segment_number '
                                f'"{segment.attributes["segment_number"]}"; '
                                f'subsegment_number "{subsegment.attributes["subsegment_number"]}";'
                            )
                            f.write(subsegment_line + "\n")

    def run(
        self,
        output_gtf,
        output_genes_tsv,
        output_segments_tsv,
        output_subsegments_tsv,
        output_pas_tsv,
    ):
        log_message("Reading annotation GTF...")
        gtf_data = read_gtf(self.annotation_file, result_type="pandas")

        log_message("Parsing GTF to Gene/Transcript/Region objects...")
        self.genes = self.parse_gtf_to_genes(gtf_data)

        log_message("Extending gene coordinates...")
        self.extend_gene_coordinates()

        log_message("Extending terminal exons...")
        self.extend_exon_downstream()

        log_message("Defining exons and introns...")
        self.define_exons_introns()

        log_message("Constructing segments...")
        self.construct_segments()

        log_message(
            "Identifying PAS sites in segments and constructing subsegments..."
        )
        self.identify_pas_in_segments()

        log_message("Writing genes, segments, subsegments, and PAS to TSV...")
        self.write_segments_pas_to_tsv(
            output_genes_tsv,
            output_segments_tsv,
            output_subsegments_tsv,
            output_pas_tsv,
        )

        log_message("Writing genes, segments and subsegments to GTF...")
        self.write_segments_pas_to_gtf(output_gtf)
