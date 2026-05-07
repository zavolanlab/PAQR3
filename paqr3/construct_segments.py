"""Construct genomic segments and sub-segments for PAQR3 analysis."""

import bisect
import logging
from collections import defaultdict

import pandas as pd  # type: ignore
from gtfparse import read_gtf  # type: ignore
from intervaltree import Interval, IntervalTree  # type: ignore
from pybedtools import BedTool  # type: ignore

from paqr3.models import Gene, Region, Transcript

logger = logging.getLogger(__name__)

# 0-based index of the RPM column in the PAS BED file (column 5 in
# 1-based BED notation). Configurable here if the atlas format changes.
_PAS_RPM_COL: int = 4

# Maps segment_origin values to single-letter type prefixes used in IDs.
# "EX" → regular exon, "IN" → intron, "TE" → terminal exon.
_ORIGIN_TO_PREFIX: dict[str, str] = {"EX": "E", "IN": "I", "TE": "T"}

# Columns excluded when building GTF attribute dicts.
_GENE_SKIP_COLS: frozenset[str] = frozenset(
    {"seqname", "source", "feature", "start", "end", "score", "strand"}
)
_TX_SKIP_COLS: frozenset[str] = _GENE_SKIP_COLS | {"gene_id"}
_EXON_SKIP_COLS: frozenset[str] = _TX_SKIP_COLS | {"transcript_id"}


class ConstructSegments:
    """Build genomic segments and sub-segments from a GTF annotation.

    Parses a GTF into Gene objects, extends terminal exons, splits genes
    into non-overlapping segments, and intersects with a PAS atlas to
    define sub-segments.

    Attributes:
        annotation_file: Path to the GTF annotation file.
        pas_atlas_file: Path to the PAS atlas BED file.
        downstream_exon_extension: Bases to extend terminal exons and
            gene boundaries downstream.
        genes: Mapping from gene ID to Gene after parsing.
        gene_mapping: Original gene ID to integer unique gene ID.
        pas_df: DataFrame of merged PAS sites from process_pas_atlas().
    """

    def __init__(
        self,
        annotation_file: str,
        pas_atlas_file: str,
        downstream_exon_extension: int,
    ) -> None:
        """Args:
            annotation_file: Path to the GTF annotation file.
            pas_atlas_file: Path to the PAS atlas BED file.
            downstream_exon_extension: Bases to extend terminal exons
                and gene boundaries in the 3′ direction.
        """
        self.annotation_file = annotation_file
        self.pas_atlas_file = pas_atlas_file
        self.downstream_exon_extension = downstream_exon_extension
        self.genes: dict[str, Gene] = {}
        self.gene_mapping: dict[str, int] = {}
        self.pas_df: pd.DataFrame | None = None

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _gene_extent(self, gene: Gene) -> tuple[int, int]:
        """Return (min_start, max_end) across all regions of gene, or (0, 0)."""
        starts = [
            r.start for t in gene.transcripts.values() for r in t.regions
        ]
        ends = [r.end for t in gene.transcripts.values() for r in t.regions]
        if not starts:
            return (0, 0)
        return min(starts), max(ends)

    def _gene_chrom(self, gene: Gene) -> str | None:
        """Return the chromosome of the first region in gene, or None."""
        return next(
            (r.chrom for t in gene.transcripts.values() for r in t.regions),
            None,
        )

    def _gene_strand(self, gene: Gene) -> str | None:
        """Return the strand of the first transcript in gene, or None."""
        return next(
            (t.strand for t in gene.transcripts.values()),
            None,
        )

    # ------------------------------------------------------------------
    # Pipeline steps
    # ------------------------------------------------------------------

    def parse_gtf_to_genes(self, gtf_data: pd.DataFrame) -> dict[str, Gene]:
        """Parse a GTF DataFrame into a gene_id → Gene mapping.

        Processes features in gene → transcript → exon order regardless
        of file ordering.

        Args:
            gtf_data: DataFrame from gtfparse.read_gtf with at minimum
                seqname, source, feature, start, end, score, strand,
                gene_id, and transcript_id columns.

        Returns:
            Dict mapping gene_id to Gene.
        """
        feature_groups = {k: v for k, v in gtf_data.groupby("feature")}
        gene_rows = feature_groups.get("gene", pd.DataFrame())
        tx_rows = feature_groups.get("transcript", pd.DataFrame())
        exon_rows = feature_groups.get("exon", pd.DataFrame())

        genes: dict[str, Gene] = {}

        for row in gene_rows.itertuples(index=False):
            row_dict = row._asdict()
            gene_id = row_dict["gene_id"]
            if pd.isna(row_dict["seqname"]):
                continue
            gene_attributes = {
                k: v for k, v in row_dict.items() if k not in _GENE_SKIP_COLS
            }
            genes[gene_id] = Gene(gene_id, attributes=gene_attributes)

        for row in tx_rows.itertuples(index=False):
            row_dict = row._asdict()
            gene_id = row_dict["gene_id"]
            if gene_id not in genes:
                continue
            transcript_id = row_dict["transcript_id"]
            tx_attributes = {
                k: v for k, v in row_dict.items() if k not in _TX_SKIP_COLS
            }
            transcript = Transcript(
                transcript_id,
                row_dict["strand"],
                attributes=tx_attributes,
            )
            genes[gene_id].add_transcript(transcript)

        for row in exon_rows.itertuples(index=False):
            row_dict = row._asdict()
            gene_id = row_dict["gene_id"]
            transcript_id = row_dict["transcript_id"]
            if (
                gene_id not in genes
                or transcript_id not in genes[gene_id].transcripts
            ):
                continue
            exon_attributes = {
                k: v for k, v in row_dict.items() if k not in _EXON_SKIP_COLS
            }
            exon = Region(
                region_type="exon",
                chrom=row_dict["seqname"],
                start=row_dict["start"] - 1,
                end=row_dict["end"],
                strand=row_dict["strand"],
                attributes=exon_attributes,
            )
            genes[gene_id].transcripts[transcript_id].add_region(exon)

        return genes

    def extend_gene_coordinates(self) -> None:
        """Extend each gene boundary downstream, capped at the nearest gene.

        Uses bisect for O(n log n) neighbour lookup. Updated start/end
        values are written into gene.attributes.
        """
        sorted_genes = sorted(
            self.genes.values(),
            key=lambda g: (
                self._gene_chrom(g) or "",
                (self._gene_extent(g)[0] if g.transcripts else float("inf")),
            ),
        )

        # Build per-(chrom, strand) sorted index once.
        # Each entry is (gene_start, Gene), sorted by gene_start.
        chrom_strand_genes: dict[tuple[str, str], list[tuple[int, Gene]]] = (
            defaultdict(list)
        )
        for gene in sorted_genes:
            strand = self._gene_strand(gene)
            chrom = self._gene_chrom(gene)
            if strand and chrom and gene.transcripts:
                gene_start, _ = self._gene_extent(gene)
                chrom_strand_genes[(chrom, strand)].append((gene_start, gene))

        # Pre-extract plain start lists for bisect (bisect operates on
        # a sequence of comparable values, not tuples).
        chrom_strand_starts: dict[tuple[str, str], list[int]] = {
            key: [s for s, _ in entries]
            for key, entries in chrom_strand_genes.items()
        }

        for gene in sorted_genes:
            gene_strand = self._gene_strand(gene)
            if gene_strand is None:
                continue
            gene_chrom = self._gene_chrom(gene)
            if gene_chrom is None or not gene.transcripts:
                continue

            gene_start, gene_end = self._gene_extent(gene)
            key = (gene_chrom, gene_strand)
            starts = chrom_strand_starts.get(key, [])
            entries = chrom_strand_genes.get(key, [])

            if gene_strand == "+":
                potential_end = gene_end + self.downstream_exon_extension
                next_gene_start = potential_end
                # bisect_right finds the first entry with start >
                # gene_end — i.e. the nearest downstream gene.
                idx = bisect.bisect_right(starts, gene_end)
                if idx < len(entries):
                    next_gene_start = min(potential_end, entries[idx][0] - 1)
                gene_end = min(potential_end, next_gene_start)

            elif gene_strand == "-":
                potential_start = gene_start - self.downstream_exon_extension
                previous_gene_end = potential_start
                # bisect_left finds the insertion point for gene_start;
                # going one back gives the nearest upstream gene.
                idx = bisect.bisect_left(starts, gene_start) - 1
                while idx >= 0:
                    _, prev_gene = entries[idx]
                    _, prev_end = self._gene_extent(prev_gene)
                    if prev_end < gene_start:
                        previous_gene_end = max(potential_start, prev_end + 1)
                        break
                    idx -= 1
                gene_start = max(potential_start, previous_gene_end)

            gene.attributes["start"] = gene_start
            gene.attributes["end"] = gene_end

    def extend_exon_downstream(self) -> None:
        """Extend the terminal exon of each transcript downstream.

        Capped at the gene boundary set by extend_gene_coordinates().
        """
        for gene in self.genes.values():
            gene_strand = self._gene_strand(gene)
            if gene_strand is None:
                continue

            for transcript in gene.transcripts.values():
                if not transcript.regions:
                    continue

                if gene_strand == "+":
                    terminal_exon = max(
                        transcript.regions,
                        key=lambda r: (
                            r.end
                            if r.region_type == "exon" and r.strand == "+"
                            else (
                                r.start
                                if r.region_type == "exon" and r.strand == "-"
                                else -1
                            )
                        ),
                    )
                    gene_end: int = gene.attributes["end"]
                    gene_start_attr: int = gene.attributes["start"]
                    if terminal_exon.strand == "+":
                        terminal_exon.end = min(
                            terminal_exon.end + self.downstream_exon_extension,
                            gene_end,
                        )
                    elif terminal_exon.strand == "-":
                        terminal_exon.start = max(
                            terminal_exon.start
                            - self.downstream_exon_extension,
                            gene_start_attr,
                        )

                elif gene_strand == "-":
                    terminal_exon = min(
                        transcript.regions,
                        key=lambda r: r.start,
                    )
                    gene_end = gene.attributes["end"]
                    gene_start_attr = gene.attributes["start"]
                    if terminal_exon.strand == "-":
                        terminal_exon.start = max(
                            terminal_exon.start
                            - self.downstream_exon_extension,
                            gene_start_attr,
                        )
                    elif terminal_exon.strand == "+":
                        terminal_exon.end = min(
                            terminal_exon.end + self.downstream_exon_extension,
                            gene_end,
                        )

    def define_exons_introns(self) -> None:
        """Infer intron regions and mark terminal exons for all transcripts.

        Introns are gaps between consecutive exons. The 3′-most exon per
        transcript gets is_terminal_exon=True in its attributes.
        """
        for gene in self.genes.values():
            for transcript in gene.transcripts.values():
                transcript.regions.sort(key=lambda r: r.start)
                exons = [
                    r for r in transcript.regions if r.region_type == "exon"
                ]

                # Mark the terminal exon per transcript
                # (post-extension, strand-aware).
                if exons:
                    if transcript.strand == "+":
                        max_end = max(e.end for e in exons)
                        for e in exons:
                            e.attributes["is_terminal_exon"] = e.end == max_end
                    else:
                        min_start = min(e.start for e in exons)
                        for e in exons:
                            e.attributes["is_terminal_exon"] = (
                                e.start == min_start
                            )

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
                                "transcript_id": (transcript.transcript_id),
                            },
                        )
                        introns.append(intron)

                transcript.regions = exons + introns
                transcript.regions.sort(key=lambda r: r.start)

    def construct_segments(self) -> None:
        """Split each gene into non-overlapping genomic segments.

        Boundaries inside terminal exons (except at the gene boundary)
        are suppressed to keep the terminal-exon region intact.
        Each segment gets a segment_number (1-based) and segment_origin
        ("IN", "TE", or "EX") in its attributes. Results stored in
        gene.segments.
        """
        for gene in self.genes.values():
            strand = self._gene_strand(gene)
            if strand is None or not gene.transcripts:
                continue

            gene_start, gene_end = self._gene_extent(gene)

            all_regions: list[tuple[int, int, Region]] = []
            for transcript in gene.transcripts.values():
                for region in transcript.regions:
                    all_regions.append((region.start, region.end, region))

            all_regions.sort(key=lambda x: (x[0], x[1]))
            segment_boundaries: set[int] = set()
            for start, end, _ in all_regions:
                segment_boundaries.add(start)
                segment_boundaries.add(end)

            # Collect terminal-exon boundaries to suppress mid-TE
            # splits while keeping the gene-boundary split.
            terminal_exons_ends: set[int] = set()
            terminal_exons_starts: set[int] = set()
            for transcript in gene.transcripts.values():
                for region in transcript.regions:
                    if region.region_type == "exon":
                        is_terminal = bool(
                            region.attributes.get("is_terminal_exon", False)
                        )
                        if is_terminal:
                            if region.strand == "+":
                                terminal_exons_ends.add(region.end)
                            else:
                                terminal_exons_starts.add(region.start)

            valid_segment_boundaries: set[int] = set()
            for boundary in segment_boundaries:
                if strand == "+":
                    if (
                        boundary in terminal_exons_ends
                        and boundary != gene_end
                    ):
                        continue
                    valid_segment_boundaries.add(boundary)
                elif strand == "-":
                    if (
                        boundary in terminal_exons_starts
                        and boundary != gene_start
                    ):
                        continue
                    valid_segment_boundaries.add(boundary)

            sorted_boundaries = sorted(valid_segment_boundaries)
            segments: list[Region] = []
            for i in range(len(sorted_boundaries) - 1):
                start = sorted_boundaries[i]
                end = sorted_boundaries[i + 1]
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

            # Pre-collect region lists for origin labelling.
            gene_introns: list[Region] = []
            gene_exons: list[Region] = []
            gene_terminal_exons: list[Region] = []
            for transcript in gene.transcripts.values():
                for r in transcript.regions:
                    if r.region_type == "intron":
                        gene_introns.append(r)
                    elif r.region_type == "exon":
                        gene_exons.append(r)
                        if r.attributes.get("is_terminal_exon", False):
                            gene_terminal_exons.append(r)

            def _overlaps(
                a_start: int,
                a_end: int,
                b_start: int,
                b_end: int,
            ) -> bool:
                return (a_start < b_end) and (b_start < a_end)

            for i, segment in enumerate(segments):
                segment.attributes["segment_number"] = i + 1
                seg_s, seg_e = segment.start, segment.end

                is_intron = any(
                    _overlaps(seg_s, seg_e, r.start, r.end)
                    for r in gene_introns
                )
                if is_intron:
                    segment.attributes["segment_origin"] = "IN"
                else:
                    is_te = any(
                        _overlaps(seg_s, seg_e, r.start, r.end)
                        for r in gene_terminal_exons
                    )
                    if is_te:
                        segment.attributes["segment_origin"] = "TE"
                    else:
                        segment.attributes["segment_origin"] = "EX"

            gene.segments = segments

    def process_pas_atlas(
        self, merge_distance: int = 5
    ) -> dict[tuple[str, str], IntervalTree]:
        """Read and merge the PAS atlas, then build interval trees.

        Sites within merge_distance on the same chrom/strand are merged
        into a cluster; the cluster stores the mean RPM and the
        highest-RPM site as rep_cs. Results are stored in self.pas_df.

        Args:
            merge_distance: Max gap in bp between sites to merge.

        Returns:
            Dict mapping (chrom, strand) to an IntervalTree of merged
            PAS intervals; each interval's data is the integer pas_id.
        """
        logger.info("Reading PAS atlas...")
        pas_bed = BedTool(self.pas_atlas_file)

        pas_sorted = sorted(
            pas_bed, key=lambda p: (p.chrom, p.strand, int(p.start))
        )

        logger.info("Merging PAS sites within %d bp...", merge_distance)
        merged_sites: list[tuple[str, int, int, str, float, str]] = []
        current_chrom: str | None = None
        current_strand: str | None = None
        current_start: int | None = None
        current_end: int | None = None
        rpm_values: list[float] = []
        rep_start: int | None = None
        rep_end: int | None = None
        rep_rpm: float | None = None
        rep_name: str | None = None

        for pas in pas_sorted:
            chrom = pas.chrom
            start = int(pas.start)
            end = int(pas.end)
            strand = pas.strand
            try:
                rpm = float(pas.fields[_PAS_RPM_COL])
            except ValueError:
                rpm = 0.0
            except IndexError:
                logger.warning(
                    "PAS entry %s:%s-%s has no column %d; "
                    "defaulting RPM to 0.0",
                    pas.chrom,
                    pas.start,
                    pas.end,
                    _PAS_RPM_COL + 1,
                )
                rpm = 0.0
            try:
                name: str | None = pas.fields[3] or None
            except IndexError:
                name = None

            if current_start is None:
                current_chrom = chrom
                current_strand = strand
                current_start = start
                current_end = end
                rpm_values = [rpm]
                rep_start, rep_end, rep_rpm, rep_name = start, end, rpm, name
            elif (
                current_end is not None
                and chrom == current_chrom
                and strand == current_strand
                and start - current_end <= merge_distance
            ):
                current_end = max(current_end, end)
                rpm_values.append(rpm)
                if rep_rpm is None or rpm > rep_rpm:
                    rep_start, rep_end, rep_rpm, rep_name = (
                        start, end, rpm, name
                    )
            else:
                assert current_chrom is not None
                assert current_strand is not None
                assert current_start is not None
                assert current_end is not None
                assert rep_start is not None
                assert rep_end is not None
                mean_rpm = sum(rpm_values) / len(rpm_values)
                if rep_name and rep_name != ".":
                    rep_cs = rep_name
                else:
                    rep_pos = (
                        rep_start
                        if rep_end == rep_start
                        else (rep_start + rep_end) // 2
                    )
                    rep_cs = f"{current_chrom}:{rep_pos}:{current_strand}"
                merged_sites.append(
                    (
                        current_chrom,
                        current_start,
                        current_end,
                        current_strand,
                        mean_rpm,
                        rep_cs,
                    )
                )
                current_chrom = chrom
                current_strand = strand
                current_start = start
                current_end = end
                rpm_values = [rpm]
                rep_start, rep_end, rep_rpm, rep_name = start, end, rpm, name

        if current_start is not None:
            assert current_chrom is not None
            assert current_strand is not None
            assert current_end is not None
            assert rep_start is not None
            assert rep_end is not None
            mean_rpm = sum(rpm_values) / len(rpm_values)
            if rep_name and rep_name != ".":
                rep_cs = rep_name
            else:
                rep_pos = (
                    rep_start
                    if rep_end == rep_start
                    else (rep_start + rep_end) // 2
                )
                rep_cs = f"{current_chrom}:{rep_pos}:{current_strand}"
            merged_sites.append(
                (
                    current_chrom,
                    current_start,
                    current_end,
                    current_strand,
                    mean_rpm,
                    rep_cs,
                )
            )

        pas_trees: dict[tuple[str, str], IntervalTree] = {}
        pas_records: list[list] = []
        unique_pas_id = 1
        for chrom, start, end, strand, atlas_rpm, rep_cs in merged_sites:
            pas_id = unique_pas_id
            pas_records.append(
                [pas_id, chrom, start, end, strand, atlas_rpm, rep_cs]
            )
            if (chrom, strand) not in pas_trees:
                pas_trees[(chrom, strand)] = IntervalTree()
            pas_trees[(chrom, strand)].add(Interval(start, end, data=pas_id))
            unique_pas_id += 1

        self.pas_df = pd.DataFrame(
            pas_records,
            columns=[
                "pas_id",
                "chrom",
                "start",
                "end",
                "strand",
                "atlas_rpm",
                "rep_cs",
            ],
        )
        return pas_trees

    def identify_pas_in_segments(
        self,
        pas_trees: dict[tuple[str, str], IntervalTree],
    ) -> None:
        """Overlap segments with PAS sites and build sub-segments.

        Splits each segment at PAS boundaries; each sub-segment gets
        the ID of the PAS it precedes, with "." for the trailing piece.
        Sub-segment order follows transcript direction. Results stored
        in segment.subsegments.

        Args:
            pas_trees: Interval trees from process_pas_atlas, keyed by
                (chrom, strand).
        """
        logger.info(
            "Identifying PAS overlaps with segments and "
            "constructing subsegments..."
        )

        for gene in self.genes.values():
            for segment in gene.segments:
                chrom = segment.chrom
                strand = segment.strand
                seg_start = segment.start
                seg_end = segment.end

                overlapping_pas_ids: list[int] = []
                subsegments: list[Region] = []

                # Cache the lookup key to avoid constructing the tuple
                # twice (once for `in`, once for indexing).
                key = (chrom, strand)
                if key in pas_trees:
                    intervals = pas_trees[key][seg_start:seg_end]
                    contained = [
                        iv
                        for iv in intervals
                        if iv.begin >= seg_start and iv.end <= seg_end
                    ]
                    if contained:
                        sorted_intervals = sorted(
                            contained, key=lambda x: x.begin
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
                                            "segment_id": (
                                                segment.attributes[
                                                    "segment_number"
                                                ]
                                            ),
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
                                        "segment_id": (
                                            segment.attributes[
                                                "segment_number"
                                            ]
                                        ),
                                        "strand": strand,
                                        "pas_id": None,
                                    },
                                )
                            )

                if not subsegments:
                    segment.attributes["overlapping_pas"] = []
                    continue

                if strand == "-":
                    subsegments.reverse()
                    overlapping_pas_ids.reverse()

                for i, sub in enumerate(subsegments):
                    sub.attributes["subsegment_number"] = i + 1
                    sub.attributes["pas_id"] = (
                        overlapping_pas_ids[i]
                        if i < len(overlapping_pas_ids)
                        else "."
                    )

                segment.attributes["overlapping_pas"] = overlapping_pas_ids
                segment.subsegments = subsegments

    def _assign_segment_ids(self) -> dict[tuple[str, int], str]:
        """Map (gene_id, segment_number) to a segment ID string.

        IDs take the form {gene_id}:{prefix}{index:03d} where prefix is
        E (exon), I (intron), or T (terminal exon), and index is a
        1-based per-type counter within the gene.

        Returns:
            Dict mapping (gene_id, segment_number) to a segment ID
            like "ENSG00000186092.1:T001".
        """
        seg_id_map: dict[tuple[str, int], str] = {}
        for gene in self.genes.values():
            segs = sorted(
                gene.segments,
                key=lambda s: s.attributes["segment_number"],
            )
            counters: dict[str, int] = {}
            for seg in segs:
                origin = seg.attributes.get("segment_origin", "EX")
                prefix = _ORIGIN_TO_PREFIX.get(origin, "E")
                counters[prefix] = counters.get(prefix, 0) + 1
                seg_num = seg.attributes["segment_number"]
                seg_id_map[(gene.gene_id, seg_num)] = (
                    f"{gene.gene_id}:{prefix}{counters[prefix]:03d}"
                )
        return seg_id_map

    def create_segments_tsv(
        self,
        out_tsv: str | None = None,
        out_debug_prefix: str | None = None,
    ) -> None:
        """Write the segments TSV and optionally debug BED files.

        One row per sub-segment in segments with at least one PAS;
        segments without PAS overlap are omitted. If out_debug_prefix
        is given, four headerless BED files are written for genes,
        segments, subsegments, and PAS sites.

        Args:
            out_tsv: Output path for the TSV; pass None to skip.
            out_debug_prefix: Path prefix for debug BEDs (no extension;
                category suffix appended automatically).
        """
        assert self.pas_df is not None, (
            "process_pas_atlas() must be called before "
            "create_segments_tsv()"
        )

        # Build lookups: integer pas_id → rep_cs, atlas_rpm, cluster coords.
        id_to_rep: dict[int, str] = {}
        id_to_rpm: dict[int, float] = {}
        id_to_pas_start: dict[int, int] = {}
        id_to_pas_end: dict[int, int] = {}
        # Also: rep_cs → (chrom, start, end, strand) for debug BED.
        rep_to_coords: dict[str, tuple[str, int, int, str]] = {}
        for row in self.pas_df.itertuples(index=False):
            pid = int(row.pas_id)
            id_to_rep[pid] = row.rep_cs
            id_to_rpm[pid] = float(row.atlas_rpm)
            id_to_pas_start[pid] = int(row.start)
            id_to_pas_end[pid] = int(row.end)
            rep_to_coords[row.rep_cs] = (
                row.chrom, int(row.start), int(row.end), row.strand
            )

        seg_id_map = self._assign_segment_ids()

        tsv_rows: list[dict] = []
        debug_rows: list[dict] = []
        seen_pas: set[str] = set()

        for gene in self.genes.values():
            first_tx = next(iter(gene.transcripts.values()), None)
            g_chrom = (
                first_tx.regions[0].chrom
                if first_tx and first_tx.regions
                else ""
            )
            g_strand = first_tx.strand if first_tx else ""
            g_start, g_end = (
                self._gene_extent(gene)
                if gene.transcripts
                else (0, 0)
            )

            if out_debug_prefix:
                debug_rows.append({
                    "chrom": g_chrom, "start": g_start,
                    "end": g_end, "id": gene.gene_id,
                    "type": "gene", "strand": g_strand,
                })

            for segment in gene.segments:
                seg_num = segment.attributes["segment_number"]
                seg_id = seg_id_map.get(
                    (gene.gene_id, seg_num),
                    f"{gene.gene_id}:E001",
                )

                if out_debug_prefix:
                    debug_rows.append({
                        "chrom": segment.chrom,
                        "start": segment.start,
                        "end": segment.end,
                        "id": seg_id,
                        "type": "segment",
                        "strand": segment.strand,
                    })

                if not segment.subsegments:
                    continue  # no PAS overlap → skip from TSV

                for sub in segment.subsegments:
                    sub_num = sub.attributes["subsegment_number"]
                    sub_id = f"{seg_id}:{sub_num:03d}"
                    int_pid = sub.attributes.get("pas_id")
                    if int_pid is not None and int_pid != ".":
                        rep_cs = id_to_rep.get(int(int_pid), ".")
                        rpm = id_to_rpm.get(int(int_pid), 0.0)
                    else:
                        rep_cs = "."
                        rpm = 0.0

                    if int_pid is not None and int_pid != ".":
                        pas_start: int | None = id_to_pas_start.get(
                            int(int_pid)
                        )
                        pas_end: int | None = id_to_pas_end.get(int(int_pid))
                    else:
                        pas_start = None
                        pas_end = None
                    tsv_rows.append({
                        "chrom": sub.chrom,
                        "start": sub.start,
                        "end": sub.end,
                        "subsegment_id": sub_id,
                        "strand": sub.strand,
                        "overlapping_pas_id": rep_cs,
                        "atlas_rpm": rpm,
                        "pas_start": pas_start,
                        "pas_end": pas_end,
                    })

                    if out_debug_prefix:
                        debug_rows.append({
                            "chrom": sub.chrom,
                            "start": sub.start,
                            "end": sub.end,
                            "id": sub_id,
                            "type": "subsegment",
                            "strand": sub.strand,
                        })
                        if rep_cs != "." and rep_cs not in seen_pas:
                            seen_pas.add(rep_cs)
                            coords = rep_to_coords.get(rep_cs)
                            if coords:
                                pc, ps, pe, pstr = coords
                                debug_rows.append({
                                    "chrom": pc,
                                    "start": ps,
                                    "end": pe,
                                    "id": rep_cs,
                                    "type": "pas",
                                    "strand": pstr,
                                })

        tsv_df = pd.DataFrame(
            tsv_rows,
            columns=[
                "chrom", "start", "end", "subsegment_id",
                "strand", "overlapping_pas_id", "atlas_rpm",
                "pas_start", "pas_end",
            ],
        )
        if out_tsv is not None:
            tsv_df.to_csv(out_tsv, sep="\t", index=False)
            logger.info(
                "Segments TSV written to %s (%d subsegments)",
                out_tsv, len(tsv_df),
            )

        if out_debug_prefix and debug_rows:
            debug_df = pd.DataFrame(
                debug_rows,
                columns=["chrom", "start", "end", "id", "type", "strand"],
            )
            for cat, suffix in (
                ("gene", "genes"),
                ("segment", "segments"),
                ("subsegment", "subsegments"),
                ("pas", "pas"),
            ):
                sub_df = (
                    debug_df[debug_df["type"] == cat]
                    .drop(columns=["type"])
                )
                path = f"{out_debug_prefix}_{suffix}.bed"
                sub_df.to_csv(path, sep="\t", index=False, header=False)
                logger.info("Debug BED written to %s", path)

    def write_segments_pas_to_bed(
        self,
        out_genes_bed: str,
        out_segments_bed: str,
        out_subsegments_bed: str,
        out_pas_bed: str,
    ) -> None:
        """Write genes, segments, sub-segments, and PAS to BED files.

        Args:
            out_genes_bed: Output path for the genes BED file.
            out_segments_bed: Output path for the segments BED file.
            out_subsegments_bed: Output path for the sub-segments BED.
            out_pas_bed: Output path for the merged PAS BED file.
        """
        gene_mapping: dict[str, int] = {}
        unique_gene_id = 1
        genes_data: list[dict] = []
        segments_data: list[dict] = []
        subsegments_data: list[dict] = []

        assert (
            self.pas_df is not None
        ), "process_pas_atlas() must be called before write_segments_pas_to_bed()"
        pas_bed_df = self.pas_df[
            ["chrom", "start", "end", "pas_id", "rep_cs", "strand"]
        ]
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

            if gene.transcripts:
                gene_start, gene_end = self._gene_extent(gene)
            else:
                gene_start, gene_end = 0, 0

            genes_data.append(
                {
                    "chrom": chrom,
                    "start": gene_start,
                    "end": gene_end,
                    "unique_gene_id": gid,
                    "gene_id": gene.gene_id,
                    "strand": strand,
                }
            )

            for segment in gene.segments:
                overlapping = segment.attributes.get("overlapping_pas", [])
                overlaps = (
                    ",".join(str(x) for x in overlapping)
                    if overlapping
                    else "."
                )
                origin = segment.attributes.get("segment_origin", "EX")
                seg_num = segment.attributes["segment_number"]
                segments_data.append(
                    {
                        "chrom": segment.chrom,
                        "start": segment.start,
                        "end": segment.end,
                        "gene_segment_id": f"{gid}.{seg_num}",
                        "overlapping_pas": overlaps,
                        "strand": segment.strand,
                        "segment_origin": origin,
                    }
                )
                if segment.subsegments:
                    for sub in segment.subsegments:
                        sub_num = sub.attributes["subsegment_number"]
                        sub_pas_id = sub.attributes.get("pas_id", ".")
                        subsegments_data.append(
                            {
                                "chrom": sub.chrom,
                                "start": sub.start,
                                "end": sub.end,
                                "gene_segment_subsegment_id": (
                                    f"{gid}.{seg_num}.{sub_num}"
                                ),
                                "pas_id": sub_pas_id,
                                "strand": sub.strand,
                                "segment_origin": origin,
                            }
                        )

        genes_df = pd.DataFrame(genes_data)
        segments_df = pd.DataFrame(segments_data)
        subsegments_df = pd.DataFrame(subsegments_data)

        genes_df.to_csv(out_genes_bed, sep="\t", index=False, header=False)
        segments_df.to_csv(
            out_segments_bed, sep="\t", index=False, header=False
        )
        subsegments_df.to_csv(
            out_subsegments_bed, sep="\t", index=False, header=False
        )

        logger.info(
            "Written PAS to %s, genes to %s, segments to %s, "
            "and subsegments to %s",
            out_pas_bed,
            out_genes_bed,
            out_segments_bed,
            out_subsegments_bed,
        )
        self.gene_mapping = gene_mapping

    def create_subsegments_dataframe(self) -> pd.DataFrame:
        """Build a flat DataFrame of all sub-segments.

        Collects non-zero-length sub-segments with human-readable IDs.
        pas_id is the rep_cs string for PAS-terminating sub-segments,
        or "." for trailing ones.

        Returns:
            DataFrame with columns gene_id, segment_id, subsegment_id,
            chrom, start, end, strand, pas_id.
        """
        assert self.pas_df is not None, (
            "process_pas_atlas() must be called before "
            "create_subsegments_dataframe()"
        )
        id_to_rep: dict[int, str] = {
            int(row.pas_id): row.rep_cs
            for row in self.pas_df.itertuples(index=False)
        }
        seg_id_map = self._assign_segment_ids()

        rows: list[dict] = []
        for gene in self.genes.values():
            for segment in gene.segments:
                if not segment.subsegments:
                    continue
                seg_num = segment.attributes["segment_number"]
                seg_id = seg_id_map.get(
                    (gene.gene_id, seg_num),
                    f"{gene.gene_id}:E001",
                )
                for sub in segment.subsegments:
                    if sub.end - sub.start <= 0:
                        continue
                    sub_num = sub.attributes["subsegment_number"]
                    raw_pid = sub.attributes.get("pas_id")
                    if raw_pid is not None and raw_pid != ".":
                        pas_id = id_to_rep.get(int(raw_pid), ".")
                    else:
                        pas_id = "."
                    rows.append(
                        {
                            "gene_id": gene.gene_id,
                            "segment_id": seg_id,
                            "subsegment_id": f"{seg_id}:{sub_num:03d}",
                            "chrom": sub.chrom,
                            "start": sub.start,
                            "end": sub.end,
                            "strand": sub.strand,
                            "pas_id": pas_id,
                        }
                    )
        return pd.DataFrame(rows)

    def compute(self, merge_distance: int = 5) -> None:
        """Run all in-memory segmentation steps without writing files.

        Call this before create_segments_tsv or
        create_subsegments_dataframe.

        Args:
            merge_distance: Max gap in bp for merging adjacent PAS
                sites in the atlas.
        """
        logger.info("Reading annotation GTF...")
        gtf_data = read_gtf(self.annotation_file, result_type="pandas")

        logger.info("Parsing GTF to Gene/Transcript/Region objects...")
        self.genes = self.parse_gtf_to_genes(gtf_data)

        logger.info("Extending gene coordinates...")
        self.extend_gene_coordinates()

        logger.info("Extending terminal exons...")
        self.extend_exon_downstream()

        logger.info("Defining exons and introns...")
        self.define_exons_introns()

        logger.info("Constructing segments...")
        self.construct_segments()

        logger.info(
            "Processing PAS atlas and building merged PAS interval tree..."
        )
        pas_trees = self.process_pas_atlas(merge_distance)

        logger.info(
            "Identifying PAS in segments and constructing subsegments..."
        )
        self.identify_pas_in_segments(pas_trees)

    def run(
        self,
        out_tsv: str | None = None,
        merge_distance: int = 5,
        out_debug_prefix: str | None = None,
    ) -> None:
        """Run the full segmentation pipeline and write output files.

        Calls compute(), then writes the segments TSV if out_tsv is
        given and debug BEDs if out_debug_prefix is given.

        Args:
            out_tsv: Output path for the segments TSV; None to skip.
            merge_distance: Max gap in bp for merging adjacent PAS.
            out_debug_prefix: Path prefix for debug BED files.
        """
        self.compute(merge_distance)
        self.create_segments_tsv(out_tsv, out_debug_prefix=out_debug_prefix)
