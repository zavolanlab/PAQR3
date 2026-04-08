"""Construct genomic segments and sub-segments for PAQR3 analysis.

This module defines :class:`ConstructSegments`, which parses a GTF
annotation file into gene/transcript/region objects, extends gene and
terminal-exon coordinates downstream, splits genes into non-overlapping
segments at exon/intron boundaries, overlaps those segments with a
poly-A site (PAS) atlas, and writes the resulting BED-format output
files.
"""

import bisect
import logging
from collections import defaultdict

import pandas as pd  # type: ignore
from gtfparse import read_gtf  # type: ignore
from intervaltree import Interval, IntervalTree  # type: ignore
from pybedtools import BedTool  # type: ignore

from paqr3.models import Gene, Region, Transcript

logger = logging.getLogger(__name__)

# Columns excluded when building GTF attribute dicts.
_GENE_SKIP_COLS: frozenset[str] = frozenset(
    {"seqname", "source", "feature", "start", "end", "score", "strand"}
)
_TX_SKIP_COLS: frozenset[str] = _GENE_SKIP_COLS | {"gene_id"}
_EXON_SKIP_COLS: frozenset[str] = _TX_SKIP_COLS | {"transcript_id"}


class ConstructSegments:
    """Pipeline for constructing genomic segments from a GTF annotation.

    The pipeline reads a GTF annotation and a poly-A site (PAS) atlas,
    then:

    1. Parses genes, transcripts, and exons into
       :class:`~paqr3.models.Gene` objects.
    2. Extends gene boundaries and terminal exons downstream by a fixed
       number of bases.
    3. Splits each gene into non-overlapping segments at exon/intron
       boundaries.
    4. Overlaps segments with merged PAS sites to define sub-segments.
    5. Writes BED-format output files for downstream analysis.

    Attributes:
        annotation_file: Path to the input GTF annotation file.
        pas_atlas_file: Path to the PAS atlas BED file.
        downstream_exon_extension: Number of bases to extend terminal
            exons and gene boundaries downstream.
        genes: Mapping from gene ID to :class:`~paqr3.models.Gene`
            after parsing.
        gene_mapping: Mapping from original gene ID to integer unique
            gene ID used in output files.
        pas_df: DataFrame of merged PAS sites produced by
            :meth:`process_pas_atlas`.
    """

    def __init__(
        self,
        annotation_file: str,
        pas_atlas_file: str,
        downstream_exon_extension: int,
    ) -> None:
        """Initialise the pipeline with input files and extension distance.

        Args:
            annotation_file: Path to the GTF annotation file.
            pas_atlas_file: Path to the PAS atlas BED file.
            downstream_exon_extension: Number of bases to extend
                terminal exons and gene boundaries in the 3′ direction.
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
        """Return the (start, end) span of all regions in *gene*.

        Computes the minimum start and maximum end across every region
        of every transcript. The gene must have at least one transcript
        with at least one region.

        Args:
            gene: A Gene with at least one region across its
                transcripts.

        Returns:
            A tuple ``(min_start, max_end)`` over all regions.
        """
        starts = [
            r.start for t in gene.transcripts.values() for r in t.regions
        ]
        ends = [r.end for t in gene.transcripts.values() for r in t.regions]
        return min(starts), max(ends)

    def _gene_chrom(self, gene: Gene) -> str | None:
        """Return the chromosome of the first region in *gene*.

        Args:
            gene: A Gene object.

        Returns:
            Chromosome name, or ``None`` if the gene has no regions.
        """
        return next(
            (r.chrom for t in gene.transcripts.values() for r in t.regions),
            None,
        )

    def _gene_strand(self, gene: Gene) -> str | None:
        """Return the strand of the first transcript in *gene*.

        Args:
            gene: A Gene object.

        Returns:
            ``"+"`` or ``"-"``, or ``None`` if the gene has no
            transcripts.
        """
        return next(
            (t.strand for t in gene.transcripts.values()),
            None,
        )

    # ------------------------------------------------------------------
    # Pipeline steps
    # ------------------------------------------------------------------

    def parse_gtf_to_genes(self, gtf_data: pd.DataFrame) -> dict[str, Gene]:
        """Parse a GTF DataFrame into a mapping of gene ID to Gene.

        Rows are grouped by feature type (``"gene"``,
        ``"transcript"``, ``"exon"``) up front, then each group is
        iterated with :meth:`~pandas.DataFrame.itertuples` rather than
        :meth:`~pandas.DataFrame.iterrows`.  This avoids constructing a
        full ``Series`` per row and gives a 3–5× speed improvement.

        Processing order is always gene → transcript → exon, regardless
        of the order features appear in the file, making the parser
        more robust to non-standard GTF orderings.

        Args:
            gtf_data: A DataFrame produced by ``gtfparse.read_gtf``,
                with at minimum the columns ``seqname``, ``source``,
                ``feature``, ``start``, ``end``, ``score``,
                ``strand``, ``gene_id``, and ``transcript_id``.

        Returns:
            A dict mapping ``gene_id`` to the corresponding
            :class:`~paqr3.models.Gene` object.
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
        """Extend each gene boundary downstream by the configured offset.

        For ``+``-strand genes the right (3′) boundary is extended; for
        ``-``-strand genes the left (3′) boundary is extended.
        Extension is capped so that it never overlaps the nearest
        downstream gene on the same chromosome and strand.

        **Complexity improvement** — the original implementation used a
        nested loop that scanned all remaining genes for each gene,
        giving O(n²) comparisons.  This method instead builds a
        per-``(chrom, strand)`` sorted list of ``(start, Gene)`` pairs
        once, then uses :func:`bisect.bisect_right` and
        :func:`bisect.bisect_left` to locate the neighbouring gene in
        O(log n) per query, reducing overall complexity to O(n log n).

        Updated ``"start"`` and ``"end"`` values are written into
        ``gene.attributes`` for use by downstream steps.
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

        For ``+``-strand transcripts the rightmost exon's end is
        extended; for ``-``-strand transcripts the leftmost exon's
        start is extended.  Extension is capped at the gene boundary
        set by :meth:`extend_gene_coordinates`.
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
        """Build intron regions and mark terminal exons for all transcripts.

        For each transcript, introns are inferred as the gaps between
        consecutive exons.  The terminal exon (3′-most relative to
        strand) is flagged with ``"is_terminal_exon": True`` in its
        attributes.  The transcript's region list is replaced with the
        combined and position-sorted list of exons and introns.
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

        Segment boundaries are placed at every exon/intron start and
        end coordinate within the gene.  Boundaries that fall inside a
        terminal exon (but not at the gene boundary) are suppressed to
        avoid splitting the terminal-exon region.

        Each resulting segment is annotated with:

        - ``"segment_number"``: 1-based index along the transcript
          direction.
        - ``"segment_origin"``: ``"IN"`` (intronic), ``"TE"``
          (terminal-exon), or ``"EX"`` (internal exon).

        Segments are stored in ``gene.segments``.
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

        PAS sites on the same chromosome and strand that are within
        *merge_distance* bases of each other are merged into a single
        cluster.  For each cluster, the mean RPM across constituent
        sites is stored and the site with the highest RPM is chosen as
        the representative cleavage site (``rep_cs``).

        Results are stored in ``self.pas_df`` and also returned as a
        nested interval-tree structure for fast overlap queries.

        Args:
            merge_distance: Maximum gap (in bp) between two PAS sites
                that are still merged into the same cluster.

        Returns:
            A dict mapping ``(chrom, strand)`` to an
            :class:`~intervaltree.IntervalTree` of merged PAS
            intervals.  Each interval's ``data`` attribute holds the
            integer ``pas_id``.
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

        for pas in pas_sorted:
            chrom = pas.chrom
            start = int(pas.start)
            end = int(pas.end)
            strand = pas.strand
            try:
                rpm = float(pas.fields[4])
            except (ValueError, IndexError):
                rpm = 0.0

            if current_start is None:
                current_chrom = chrom
                current_strand = strand
                current_start = start
                current_end = end
                rpm_values = [rpm]
                rep_start, rep_end, rep_rpm = start, end, rpm
            elif (
                current_end is not None
                and chrom == current_chrom
                and strand == current_strand
                and start - current_end <= merge_distance
            ):
                current_end = max(current_end, end)
                rpm_values.append(rpm)
                if rep_rpm is None or rpm > rep_rpm:
                    rep_start, rep_end, rep_rpm = start, end, rpm
            else:
                assert current_chrom is not None
                assert current_strand is not None
                assert current_start is not None
                assert current_end is not None
                assert rep_start is not None
                assert rep_end is not None
                mean_rpm = sum(rpm_values) / len(rpm_values)
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
                rep_start, rep_end, rep_rpm = start, end, rpm

        if current_start is not None:
            assert current_chrom is not None
            assert current_strand is not None
            assert current_end is not None
            assert rep_start is not None
            assert rep_end is not None
            mean_rpm = sum(rpm_values) / len(rpm_values)
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

        For each segment, PAS intervals fully contained within the
        segment are identified.  The segment is then split into
        sub-segments at the PAS boundaries: each sub-segment preceding
        a PAS receives that PAS's ID, and the final tail sub-segment
        receives ``"."``.  On the ``"-"`` strand, sub-segment order is
        reversed to follow transcript direction.

        Sub-segments are stored in ``segment.subsegments``.

        Args:
            pas_trees: Interval trees returned by
                :meth:`process_pas_atlas`, keyed by
                ``(chrom, strand)``.
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

                if (chrom, strand) in pas_trees:
                    intervals = pas_trees[(chrom, strand)][seg_start:seg_end]
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

    def write_segments_pas_to_bed(
        self,
        out_genes_bed: str,
        out_segments_bed: str,
        out_subsegments_bed: str,
        out_pas_bed: str,
    ) -> None:
        """Write genes, segments, sub-segments, and PAS to BED files.

        Output format (all tab-separated, no header):

        - *out_genes_bed*: chrom, start, end, unique_gene_id,
          gene_id, strand.
        - *out_segments_bed*: chrom, start, end, gene_segment_id,
          overlapping_pas, strand, segment_origin.
        - *out_subsegments_bed*: chrom, start, end,
          gene_segment_subsegment_id, pas_id, strand,
          segment_origin.
        - *out_pas_bed*: chrom, start, end, pas_id, rep_cs, strand.

        Args:
            out_genes_bed: Output path for the genes BED file.
            out_segments_bed: Output path for the segments BED file.
            out_subsegments_bed: Output path for the sub-segments BED
                file.
            out_pas_bed: Output path for the merged PAS BED file.
        """
        gene_mapping: dict[str, int] = {}
        unique_gene_id = 1
        genes_data: list[list] = []
        segments_data: list[list] = []
        subsegments_data: list[list] = []

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
                [chrom, gene_start, gene_end, gid, gene.gene_id, strand]
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
                    [
                        segment.chrom,
                        segment.start,
                        segment.end,
                        f"{gid}.{seg_num}",
                        overlaps,
                        segment.strand,
                        origin,
                    ]
                )
                if segment.subsegments:
                    for sub in segment.subsegments:
                        sub_num = sub.attributes["subsegment_number"]
                        sub_pas_id = sub.attributes.get("pas_id", ".")
                        subsegments_data.append(
                            [
                                sub.chrom,
                                sub.start,
                                sub.end,
                                f"{gid}.{seg_num}.{sub_num}",
                                sub_pas_id,
                                sub.strand,
                                origin,
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
                "segment_origin",
            ],
        )
        subsegments_df = pd.DataFrame(
            subsegments_data,
            columns=[
                "chrom",
                "start",
                "end",
                "gene_segment_subsegment_id",
                "pas_id",
                "strand",
                "segment_origin",
            ],
        )

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
        """Build a DataFrame of all sub-segments for coverage calculations.

        Iterates over ``gene.segments`` and collects every sub-segment
        with a non-zero length into a flat table.

        Returns:
            A DataFrame with columns ``gene_unique_id``,
            ``segment_number``, ``subsegment_number``, ``chrom``,
            ``start``, ``end``, ``strand``, ``pas_id``.
        """
        rows: list[list] = []
        for gene in self.genes.values():
            for segment in gene.segments:
                if not segment.subsegments:
                    continue
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
        out_genes_bed: str,
        out_segments_bed: str,
        out_subsegments_bed: str,
        out_pas_bed: str,
        merge_distance: int = 5,
    ) -> None:
        """Execute the full segment-construction pipeline.

        Steps, in order:

        1. Read and parse the GTF annotation.
        2. Extend gene coordinate boundaries downstream.
        3. Extend terminal exons downstream.
        4. Infer intron regions and mark terminal exons.
        5. Construct non-overlapping segments per gene.
        6. Read and merge the PAS atlas.
        7. Identify PAS within segments and build sub-segments.
        8. Write all output BED files.

        Args:
            out_genes_bed: Output path for the genes BED file.
            out_segments_bed: Output path for the segments BED file.
            out_subsegments_bed: Output path for the sub-segments BED
                file.
            out_pas_bed: Output path for the merged PAS BED file.
            merge_distance: Maximum gap (in bp) for merging adjacent
                PAS sites in the atlas.
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

        logger.info(
            "Writing genes, segments, subsegments, and merged PAS "
            "to BED/TSV..."
        )
        self.write_segments_pas_to_bed(
            out_genes_bed,
            out_segments_bed,
            out_subsegments_bed,
            out_pas_bed,
        )
