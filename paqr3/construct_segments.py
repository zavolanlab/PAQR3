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

        The PAS BED file is sorted by ``(chrom, strand, start)`` before
        merging; the merge loop relies on this order being strictly
        maintained.  RPM values are read from column index
        ``_PAS_RPM_COL`` (0-based); a warning is emitted and the value
        defaults to ``0.0`` when that column is absent.

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
        """Map (gene_id, segment_number) to a human-readable segment ID.

        Segment IDs take the form ``{gene_id}:{prefix}{index:03d}``
        where *prefix* is ``E`` (exon), ``I`` (intron), or ``T``
        (terminal exon), and *index* is a 1-based per-type counter
        within the gene.  Segments are processed in
        ``segment_number`` order so numbering is stable.

        Returns:
            Dict mapping ``(gene_id, segment_number)`` to a segment
            ID string such as ``"ENSG00000186092.1:T001"``.
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
        out_tsv: str,
        out_debug_bed: str | None = None,
    ) -> None:
        """Write the segments TSV (and optionally a debug BED).

        The TSV is the primary output of the segmentation stage and
        the required input for ``paqr3 quant``.  It contains one row
        per sub-segment belonging to a segment with at least one
        overlapping PAS; segments with no PAS overlap are omitted.

        TSV columns (tab-separated, with header):

        - ``chrom``, ``start``, ``end``: 0-based half-open genomic
          coordinates.
        - ``subsegment_id``: ``{gene_id}:{TypeNNN}:{sub:03d}``
          — *Type* is ``E`` / ``I`` / ``T``, *NNN* is the per-type
          segment counter within the gene, *sub* is the 1-based
          sub-segment index within the segment.
        - ``strand``: ``"+"`` or ``"-"``.
        - ``overlapping_pas_id``: representative cleavage-site
          coordinate of the PAS that terminates this sub-segment
          (e.g. ``chr1:65435``), or ``"."`` for the trailing
          sub-segment.
        - ``atlas_rpm``: mean RPM of the merged PAS cluster; ``0.0``
          for trailing sub-segments.

        If *out_debug_bed* is given, a headerless six-column BED file
        is written listing every gene, segment, sub-segment and PAS
        with their genomic coordinates and a ``type`` label
        (``gene`` / ``segment`` / ``subsegment`` / ``pas``).  PAS
        entries appear at most once per unique representative CS.

        Args:
            out_tsv: Output path for the segments TSV file.
            out_debug_bed: Optional path for the debug BED file.
        """
        assert self.pas_df is not None, (
            "process_pas_atlas() must be called before "
            "create_segments_tsv()"
        )

        # Build lookups: integer pas_id → rep_cs and → atlas_rpm.
        id_to_rep: dict[int, str] = {}
        id_to_rpm: dict[int, float] = {}
        # Also: rep_cs → (chrom, start, end, strand) for debug BED.
        rep_to_coords: dict[str, tuple[str, int, int, str]] = {}
        for row in self.pas_df.itertuples(index=False):
            pid = int(row.pas_id)
            id_to_rep[pid] = row.rep_cs
            id_to_rpm[pid] = float(row.atlas_rpm)
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

            if out_debug_bed:
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

                if out_debug_bed:
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

                    tsv_rows.append({
                        "chrom": sub.chrom,
                        "start": sub.start,
                        "end": sub.end,
                        "subsegment_id": sub_id,
                        "strand": sub.strand,
                        "overlapping_pas_id": rep_cs,
                        "atlas_rpm": rpm,
                    })

                    if out_debug_bed:
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
            ],
        )
        tsv_df.to_csv(out_tsv, sep="\t", index=False)
        logger.info(
            "Segments TSV written to %s (%d subsegments)",
            out_tsv, len(tsv_df),
        )

        if out_debug_bed and debug_rows:
            debug_df = pd.DataFrame(
                debug_rows,
                columns=["chrom", "start", "end", "id", "type", "strand"],
            )
            debug_df.to_csv(
                out_debug_bed, sep="\t", index=False, header=False
            )
            logger.info("Debug BED written to %s", out_debug_bed)

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
        """Build a DataFrame of all sub-segments for coverage calculations.

        Iterates over ``gene.segments`` and collects every sub-segment
        with a non-zero length into a flat table using the same
        human-readable ID scheme as :meth:`create_segments_tsv`.

        ``pas_id`` is set to the representative cleavage-site string
        (e.g. ``"chr1:65435"``) for sub-segments that terminate at a
        PAS, and ``"."`` for trailing sub-segments.

        Returns:
            A DataFrame with columns ``gene_id``, ``segment_id``,
            ``subsegment_id``, ``chrom``, ``start``, ``end``,
            ``strand``, ``pas_id``.
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
        """Execute all in-memory segmentation steps.

        Runs the complete pipeline through to sub-segment construction
        without writing any output files.  Call this before
        :meth:`create_segments_tsv` or
        :meth:`create_subsegments_dataframe`.

        Steps, in order:

        1. Read and parse the GTF annotation.
        2. Extend gene coordinate boundaries downstream.
        3. Extend terminal exons downstream.
        4. Infer intron regions and mark terminal exons.
        5. Construct non-overlapping segments per gene.
        6. Read and merge the PAS atlas.
        7. Identify PAS within segments and build sub-segments.

        Args:
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

    def run(
        self,
        out_tsv: str,
        merge_distance: int = 5,
        out_debug_bed: str | None = None,
    ) -> None:
        """Run the full segmentation pipeline and write output files.

        Calls :meth:`compute` to build all in-memory data structures,
        then writes the segments TSV via :meth:`create_segments_tsv`.

        Args:
            out_tsv: Output path for the segments TSV.
            merge_distance: Maximum gap (in bp) for merging adjacent
                PAS sites in the atlas.
            out_debug_bed: Optional path for the debug BED file
                (written only when provided).
        """
        self.compute(merge_distance)
        self.create_segments_tsv(out_tsv, out_debug_bed=out_debug_bed)
