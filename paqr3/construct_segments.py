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
        # After process_pas_atlas, this will hold per‐PAS info including merged RPM
        self.pas_df = None

    def parse_gtf_to_genes(self, gtf_data):
        """Parses GTF data (pandas DataFrame) into a dictionary of Gene objects."""
        genes = {}
        for _, row in gtf_data.iterrows():
            if row["feature"] == "gene":
                gene_id = row["gene_id"]
                if pd.isna(row["seqname"]):
                    continue

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
                        start=row["start"] - 1,
                        end=row["end"],
                        strand=row["strand"],
                        attributes=exon_attributes,
                    )
                    genes[gene_id].transcripts[transcript_id].add_region(exon)

        return genes

    def extend_gene_coordinates(self):
        """Extends gene coordinates (handles strand) for downstream exons."""
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
        """Extends terminal exons (handles strand)."""
        for gene in self.genes.values():
            gene_strand = next(
                (t.strand for t in gene.transcripts.values()), None
            )
            if gene_strand is None:
                continue

            for transcript in gene.transcripts.values():
                if transcript.regions:
                    if gene_strand == "+":
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

                    elif gene_strand == "-":
                        terminal_exon = min(
                            transcript.regions,
                            key=lambda r: r.start,
                        )
                        if terminal_exon.strand == "-":
                            terminal_exon.start = max(
                                terminal_exon.start
                                - self.downstream_exon_extension,
                                gene.attributes["start"],
                            )
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
                transcript.regions.sort(key=lambda r: r.start)
                exons = [
                    r for r in transcript.regions if r.region_type == "exon"
                ]
                introns = []

                for i in range(len(exons) - 1):
                    exon1 = exons[i]
                    exon2 = exons[i + 1]
                    intron_start = exon1.end
                    intron_end = exon2.start
                    if intron_start <= intron_end:
                        intron = Region(
                            region_type="intron",
                            chrom=exon1.chrom,
                            start=intron_start,
                            end=intron_end,
                            strand=exon1.strand,
                            intron_number=i + 1,
                            attributes={
                                "gene_id": gene.gene_id,
                                "transcript_id": transcript.transcript_id,
                            },
                        )
                        introns.append(intron)

                transcript.regions = exons + introns
                transcript.regions.sort(key=lambda r: r.start)

    def construct_segments(self):
        """Constructs segments for each gene."""
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
                segment_boundaries.add(end)

            gene_start = (
                min(
                    r.start
                    for t in gene.transcripts.values()
                    for r in t.regions
                )
                if gene.transcripts
                else float("inf")
            )
            gene_end = (
                max(
                    r.end for t in gene.transcripts.values() for r in t.regions
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

            valid_segment_boundaries = set()
            for boundary in segment_boundaries:
                if strand == "+":
                    if (
                        boundary in terminal_exons_ends
                        and boundary != gene_end
                    ):
                        continue
                    else:
                        valid_segment_boundaries.add(boundary)
                elif strand == "-":
                    if (
                        boundary in terminal_exons_starts
                        and boundary != gene_start
                    ):
                        continue
                    else:
                        valid_segment_boundaries.add(boundary)

            valid_segment_boundaries = sorted(valid_segment_boundaries)
            segments = []
            for i in range(len(valid_segment_boundaries) - 1):
                start = valid_segment_boundaries[i]
                end = valid_segment_boundaries[i + 1]
                if start < end:
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

    def process_pas_atlas(self, merge_distance=5):
        """
        Read the PAS‐atlas BED file, merge sites within merge_distance bp,
        compute mean RPM for merged sites, assign new IDs, and build an interval tree.
        Stores resulting PAS info in self.pas_df (columns: pas_id, chrom, start, end, strand, atlas_rpm).
        """
        log_message("Reading PAS atlas...")
        pas_bed = BedTool(self.pas_atlas_file)

        pas_sorted = sorted(
            pas_bed, key=lambda p: (p.chrom, p.strand, int(p.start))
        )

        log_message(f"Merging PAS sites within {merge_distance} bp...")
        merged_sites = []
        current_chrom = None
        current_strand = None
        current_start = None
        current_end = None
        rpm_values = []
        for pas in pas_sorted:
            chrom = pas.chrom
            start = int(pas.start)
            end = int(pas.end)
            strand = pas.strand
            # 5th column is RPM
            try:
                rpm = float(pas.fields[4])
            except (ValueError, IndexError):
                rpm = 0.0

            if current_start is None:
                # initialize first cluster
                current_chrom = chrom
                current_strand = strand
                current_start = start
                current_end = end
                rpm_values = [rpm]
            elif (
                chrom == current_chrom
                and strand == current_strand
                and start - current_end <= merge_distance
            ):
                # extend cluster
                current_end = max(current_end, end)
                rpm_values.append(rpm)
            else:
                # finalize previous cluster
                mean_rpm = sum(rpm_values) / len(rpm_values)
                merged_sites.append(
                    (
                        current_chrom,
                        current_start,
                        current_end,
                        current_strand,
                        mean_rpm,
                    )
                )
                # start new cluster
                current_chrom = chrom
                current_strand = strand
                current_start = start
                current_end = end
                rpm_values = [rpm]

        # finalize last
        if current_start is not None:
            mean_rpm = sum(rpm_values) / len(rpm_values)
            merged_sites.append(
                (
                    current_chrom,
                    current_start,
                    current_end,
                    current_strand,
                    mean_rpm,
                )
            )

        # Build interval trees and pas_data list
        pas_trees = {}
        pas_records = (
            []
        )  # list of [pas_id, chrom, start, end, strand, atlas_rpm]
        unique_pas_id = 1
        for chrom, start, end, strand, atlas_rpm in merged_sites:
            pas_id = unique_pas_id
            pas_records.append([pas_id, chrom, start, end, strand, atlas_rpm])

            if (chrom, strand) not in pas_trees:
                pas_trees[(chrom, strand)] = IntervalTree()
            pas_trees[(chrom, strand)].add(Interval(start, end, data=pas_id))

            unique_pas_id += 1

        # Save DataFrame for downstream use
        self.pas_df = pd.DataFrame(
            pas_records,
            columns=["pas_id", "chrom", "start", "end", "strand", "atlas_rpm"],
        )

        return pas_trees

    def identify_pas_in_segments(self, pas_trees):
        """
        Identifies merged PAS sites overlapping segments and constructs subsegments.
        Assigns the merged pas_id to each subsegment.
        """
        log_message(
            "Identifying PAS overlaps with segments and constructing subsegments..."
        )

        for gene in self.genes.values():
            for segment in gene.segments:
                chrom = segment.chrom
                strand = segment.strand
                seg_start = segment.start
                seg_end = segment.end

                overlapping_pas_ids = []
                subsegments = []

                if (chrom, strand) in pas_trees:
                    intervals = pas_trees[(chrom, strand)][seg_start:seg_end]
                    sorted_intervals = sorted(
                        [
                            iv
                            for iv in intervals
                            if iv.begin >= seg_start and iv.end <= seg_end
                        ],
                        key=lambda x: x.begin,
                    )

                    current_start = seg_start
                    for iv in sorted_intervals:
                        iv_start = iv.begin
                        iv_end = iv.end
                        pas_id = iv.data

                        if iv_start > current_start:
                            subsegments.append(
                                Region(
                                    region_type="subsegment",
                                    chrom=chrom,
                                    start=current_start,
                                    end=iv_start,
                                    strand=strand,
                                    attributes={
                                        "gene_id": gene.gene_id,
                                        "segment_id": segment.attributes[
                                            "segment_number"
                                        ],
                                        "strand": strand,
                                        "pas_id": None,
                                    },
                                )
                            )
                        current_start = iv_end
                        overlapping_pas_ids.append(pas_id)

                    if current_start < seg_end:
                        subsegments.append(
                            Region(
                                region_type="subsegment",
                                chrom=chrom,
                                start=current_start,
                                end=seg_end,
                                strand=strand,
                                attributes={
                                    "gene_id": gene.gene_id,
                                    "segment_id": segment.attributes[
                                        "segment_number"
                                    ],
                                    "strand": strand,
                                    "pas_id": None,
                                },
                            )
                        )

                if strand == "-":
                    subsegments.reverse()

                for i, sub in enumerate(subsegments):
                    sub.attributes["subsegment_number"] = i + 1
                    # assign pas_id in order
                    pas_index = (
                        i
                        if strand == "+"
                        else len(overlapping_pas_ids) - 1 - i
                    )
                    if 0 <= pas_index < len(overlapping_pas_ids):
                        sub.attributes["pas_id"] = overlapping_pas_ids[
                            pas_index
                        ]
                    else:
                        sub.attributes["pas_id"] = "."

                segment.attributes["overlapping_pas"] = overlapping_pas_ids
                segment.subsegments = subsegments

    def write_segments_pas_to_bed(
        self, out_genes_bed, out_segments_bed, out_subsegments_bed, out_pas_bed
    ):
        """Writes genes, segments, subsegments, and merged PAS to BED files."""
        gene_mapping = {}
        unique_gene_id = 1
        genes_data = []
        segments_data = []
        subsegments_data = []

        pas_bed_df = self.pas_df[
            ["chrom", "start", "end", "pas_id", "atlas_rpm", "strand"]
        ]
        # output columns: chrom, start, end, pas_id, atlas_rpm, strand
        pas_bed_df.to_csv(out_pas_bed, sep="\t", index=False, header=False)

        for gene in self.genes.values():
            first_tx = next(iter(gene.transcripts.values()), None)
            chrom = (
                first_tx.regions[0].chrom
                if first_tx and first_tx.regions
                else ""
            )
            strand = first_tx.strand if first_tx else ""

            gid = unique_gene_id
            gene_mapping[gene.gene_id] = gid
            unique_gene_id += 1

            genes_data.append(
                [
                    chrom,
                    min(
                        r.start
                        for t in gene.transcripts.values()
                        for r in t.regions
                    ),
                    max(
                        r.end
                        for t in gene.transcripts.values()
                        for r in t.regions
                    ),
                    gid,
                    gene.gene_id,
                    strand,
                ]
            )

            for segment in gene.segments:
                overlapping = segment.attributes.get("overlapping_pas", [])
                overlaps = (
                    ",".join(str(x) for x in overlapping)
                    if overlapping
                    else "."
                )
                segments_data.append(
                    [
                        segment.chrom,
                        segment.start,
                        segment.end,
                        f"{gid}.{segment.attributes['segment_number']}",
                        overlaps,
                        segment.strand,
                    ]
                )
                if hasattr(segment, "subsegments"):
                    for sub in segment.subsegments:
                        subsegments_data.append(
                            [
                                sub.chrom,
                                sub.start,
                                sub.end,
                                f"{gid}.{segment.attributes['segment_number']}.{sub.attributes['subsegment_number']}",
                                ".",
                                sub.strand,
                            ]
                        )

        genes_df = pd.DataFrame(
            genes_data,
            columns=[
                "chrom",
                "start",
                "end",
                "unique_gene_id",
                "gene_id",
                "strand",
            ],
        )
        segments_df = pd.DataFrame(
            segments_data,
            columns=[
                "chrom",
                "start",
                "end",
                "gene_segment_id",
                "overlapping_pas",
                "strand",
            ],
        )
        subsegments_df = pd.DataFrame(
            subsegments_data,
            columns=[
                "chrom",
                "start",
                "end",
                "gene_segment_subsegment_id",
                "score",
                "strand",
            ],
        )

        genes_df.to_csv(out_genes_bed, sep="\t", index=False, header=False)
        segments_df.to_csv(
            out_segments_bed, sep="\t", index=False, header=False
        )
        subsegments_df.to_csv(
            out_subsegments_bed, sep="\t", index=False, header=False
        )

        log_message(
            f"Written PAS to {out_pas_bed}, genes to {out_genes_bed}, segments to {out_segments_bed}, and subsegments to {out_subsegments_bed}"
        )
        self.gene_mapping = gene_mapping

    def create_subsegments_dataframe(self):
        """Turn in‐memory subsegments into a DataFrame for coverage calculations."""
        rows = []
        for gene in self.genes.values():
            for segment in gene.segments:
                if hasattr(segment, "subsegments"):
                    for sub in segment.subsegments:
                        if sub.end - sub.start > 0:
                            rows.append(
                                [
                                    self.gene_mapping[gene.gene_id],
                                    segment.attributes["segment_number"],
                                    sub.attributes["subsegment_number"],
                                    sub.chrom,
                                    sub.start,
                                    sub.end,
                                    sub.strand,
                                    sub.attributes.get("pas_id", "."),
                                ]
                            )
        return pd.DataFrame(
            rows,
            columns=[
                "gene_unique_id",
                "segment_number",
                "subsegment_number",
                "chrom",
                "start",
                "end",
                "strand",
                "pas_id",
            ],
        )

    def run(
        self,
        out_genes_bed,
        out_segments_bed,
        out_subsegments_bed,
        out_pas_bed,
        merge_distance=5,
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
            "Processing PAS atlas and building merged PAS interval tree..."
        )
        pas_trees = self.process_pas_atlas(merge_distance)

        log_message(
            "Identifying PAS in segments and constructing subsegments..."
        )
        self.identify_pas_in_segments(pas_trees)

        log_message(
            "Writing genes, segments, subsegments, and merged PAS to BED/TSV..."
        )
        self.write_segments_pas_to_bed(
            out_genes_bed,
            out_segments_bed,
            out_subsegments_bed,
            out_pas_bed,
        )
