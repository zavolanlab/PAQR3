"""Unit tests for paqr3.construct_segments."""

import os

import pandas as pd
import pytest

from paqr3.construct_segments import ConstructSegments
from paqr3.models import Gene, Region, Transcript


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_gene(
    gene_id: str,
    exon_intervals: list[tuple[int, int]],
    strand: str = "+",
    chrom: str = "chr1",
) -> Gene:
    """Build a Gene with one transcript and the given exon intervals."""
    gene = Gene(gene_id, attributes={"gene_id": gene_id})
    tx = Transcript(f"{gene_id}_tx1", strand)
    for i, (start, end) in enumerate(exon_intervals):
        tx.add_region(
            Region(
                region_type="exon",
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                exon_number=i + 1,
                attributes={
                    "gene_id": gene_id,
                    "transcript_id": f"{gene_id}_tx1",
                },
            )
        )
    gene.add_transcript(tx)
    return gene


def _make_cs(**kwargs) -> ConstructSegments:
    """Return a ConstructSegments instance with dummy file paths."""
    return ConstructSegments(
        annotation_file=kwargs.get("annotation_file", ""),
        pas_atlas_file=kwargs.get("pas_atlas_file", ""),
        downstream_exon_extension=kwargs.get("downstream_exon_extension", 200),
    )


def _write_bed(path: str, rows: list) -> None:
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")


# ---------------------------------------------------------------------------
# Private helpers: _gene_extent, _gene_chrom, _gene_strand
# ---------------------------------------------------------------------------


class TestGeneHelpers:
    def test_gene_extent_single_exon(self):
        cs = _make_cs()
        gene = _make_gene("g1", [(100, 200)])
        assert cs._gene_extent(gene) == (100, 200)

    def test_gene_extent_multiple_exons(self):
        cs = _make_cs()
        gene = _make_gene("g1", [(100, 200), (400, 600)])
        assert cs._gene_extent(gene) == (100, 600)

    def test_gene_chrom(self):
        cs = _make_cs()
        gene = _make_gene("g1", [(100, 200)], chrom="chrX")
        assert cs._gene_chrom(gene) == "chrX"

    def test_gene_strand_plus(self):
        cs = _make_cs()
        gene = _make_gene("g1", [(100, 200)], strand="+")
        assert cs._gene_strand(gene) == "+"

    def test_gene_strand_minus(self):
        cs = _make_cs()
        gene = _make_gene("g1", [(100, 200)], strand="-")
        assert cs._gene_strand(gene) == "-"

    def test_gene_chrom_no_transcripts_returns_none(self):
        cs = _make_cs()
        assert cs._gene_chrom(Gene("g_empty")) is None

    def test_gene_strand_no_transcripts_returns_none(self):
        cs = _make_cs()
        assert cs._gene_strand(Gene("g_empty")) is None


# ---------------------------------------------------------------------------
# define_exons_introns
# ---------------------------------------------------------------------------


class TestDefineExonsIntrons:
    def _run(self, gene):
        cs = _make_cs()
        cs.genes = {gene.gene_id: gene}
        cs.define_exons_introns()
        return cs

    def test_intron_created_between_two_exons(self):
        gene = _make_gene("g1", [(1000, 1500), (1700, 2000)])
        self._run(gene)
        tx = gene.transcripts["g1_tx1"]
        introns = [r for r in tx.regions if r.region_type == "intron"]
        assert len(introns) == 1
        assert introns[0].start == 1500
        assert introns[0].end == 1700

    def test_three_exons_produce_two_introns(self):
        gene = _make_gene("g1", [(100, 200), (300, 400), (500, 600)])
        self._run(gene)
        tx = gene.transcripts["g1_tx1"]
        introns = [r for r in tx.regions if r.region_type == "intron"]
        assert len(introns) == 2

    def test_single_exon_no_intron(self):
        gene = _make_gene("g1", [(1000, 2000)])
        self._run(gene)
        tx = gene.transcripts["g1_tx1"]
        introns = [r for r in tx.regions if r.region_type == "intron"]
        assert len(introns) == 0

    def test_terminal_exon_on_plus_strand(self):
        gene = _make_gene("g1", [(1000, 1500), (1700, 2000)], strand="+")
        self._run(gene)
        tx = gene.transcripts["g1_tx1"]
        exons = [r for r in tx.regions if r.region_type == "exon"]
        terminals = [e for e in exons if e.attributes.get("is_terminal_exon")]
        assert len(terminals) == 1
        assert terminals[0].end == 2000

    def test_terminal_exon_on_minus_strand(self):
        gene = _make_gene("g1", [(1000, 1500), (1700, 2000)], strand="-")
        self._run(gene)
        tx = gene.transcripts["g1_tx1"]
        exons = [r for r in tx.regions if r.region_type == "exon"]
        terminals = [e for e in exons if e.attributes.get("is_terminal_exon")]
        assert len(terminals) == 1
        assert terminals[0].start == 1000

    def test_non_terminal_exons_marked_false(self):
        gene = _make_gene("g1", [(1000, 1500), (1700, 2000)], strand="+")
        self._run(gene)
        tx = gene.transcripts["g1_tx1"]
        exons = sorted(
            [r for r in tx.regions if r.region_type == "exon"],
            key=lambda r: r.start,
        )
        assert exons[0].attributes["is_terminal_exon"] is False
        assert exons[1].attributes["is_terminal_exon"] is True


# ---------------------------------------------------------------------------
# construct_segments
# ---------------------------------------------------------------------------


class TestConstructSegments:
    def _setup_and_run(
        self,
        gene_id: str,
        exon_intervals: list[tuple[int, int]],
        strand: str = "+",
        ext: int = 200,
    ):
        gene = _make_gene(gene_id, exon_intervals, strand=strand)
        # Set gene boundary attributes as extend_gene_coordinates would.
        lo = min(s for s, e in exon_intervals)
        hi = max(e for s, e in exon_intervals)
        gene.attributes["start"] = lo - ext
        gene.attributes["end"] = hi + ext
        cs = _make_cs(downstream_exon_extension=ext)
        cs.genes = {gene_id: gene}
        cs.define_exons_introns()
        cs.construct_segments()
        return cs, gene

    def test_segments_non_overlapping(self):
        _, gene = self._setup_and_run("g1", [(1000, 1500), (1700, 2000)])
        segs = sorted(gene.segments, key=lambda s: s.start)
        for i in range(len(segs) - 1):
            assert segs[i].end <= segs[i + 1].start

    def test_segment_origins_labelled(self):
        _, gene = self._setup_and_run("g1", [(1000, 1500), (1700, 2000)])
        for seg in gene.segments:
            assert "segment_origin" in seg.attributes
            assert seg.attributes["segment_origin"] in ("EX", "IN", "TE")

    def test_segment_numbers_sequential(self):
        _, gene = self._setup_and_run("g1", [(1000, 1500), (1700, 2000)])
        nums = [s.attributes["segment_number"] for s in gene.segments]
        assert nums == list(range(1, len(nums) + 1))

    def test_terminal_exon_is_one_segment(self):
        _, gene = self._setup_and_run("g1", [(1000, 1500), (1700, 2000)])
        te_segs = [
            s
            for s in gene.segments
            if s.attributes.get("segment_origin") == "TE"
        ]
        assert len(te_segs) == 1

    def test_minus_strand_segments_descending(self):
        _, gene = self._setup_and_run(
            "g1", [(1000, 1500), (1700, 2000)], strand="-"
        )
        starts = [s.start for s in gene.segments]
        assert starts == sorted(starts, reverse=True)

    def test_single_exon_gene_gets_one_te_segment(self):
        _, gene = self._setup_and_run("g1", [(1000, 2000)])
        origins = [s.attributes.get("segment_origin") for s in gene.segments]
        assert all(o == "TE" for o in origins)


# ---------------------------------------------------------------------------
# process_pas_atlas
# ---------------------------------------------------------------------------


class TestProcessPasAtlas:
    def test_basic_loading(self, tmp_path):
        bed = str(tmp_path / "atlas.bed")
        _write_bed(bed, [
            ["chr1", 1000, 1010, "cs1", 50.0, "+"],
            ["chr1", 2000, 2010, "cs2", 30.0, "+"],
        ])
        cs = ConstructSegments("", bed, 200)
        trees = cs.process_pas_atlas(merge_distance=5)

        assert cs.pas_df is not None
        assert len(cs.pas_df) == 2
        assert ("chr1", "+") in trees

    def test_nearby_sites_merged(self, tmp_path):
        bed = str(tmp_path / "atlas.bed")
        # Gap of 3 bp, within merge_distance=5
        _write_bed(bed, [
            ["chr1", 1000, 1005, "cs1", 40.0, "+"],
            ["chr1", 1008, 1012, "cs2", 60.0, "+"],
        ])
        cs = ConstructSegments("", bed, 200)
        cs.process_pas_atlas(merge_distance=5)
        assert len(cs.pas_df) == 1

    def test_distant_sites_not_merged(self, tmp_path):
        bed = str(tmp_path / "atlas.bed")
        # Gap of 10 bp, beyond merge_distance=5
        _write_bed(bed, [
            ["chr1", 1000, 1005, "cs1", 40.0, "+"],
            ["chr1", 1015, 1020, "cs2", 60.0, "+"],
        ])
        cs = ConstructSegments("", bed, 200)
        cs.process_pas_atlas(merge_distance=5)
        assert len(cs.pas_df) == 2

    def test_different_strands_not_merged(self, tmp_path):
        bed = str(tmp_path / "atlas.bed")
        _write_bed(bed, [
            ["chr1", 1000, 1005, "cs1", 40.0, "+"],
            ["chr1", 1000, 1005, "cs2", 60.0, "-"],
        ])
        cs = ConstructSegments("", bed, 200)
        trees = cs.process_pas_atlas(merge_distance=5)
        assert len(cs.pas_df) == 2
        assert ("chr1", "+") in trees
        assert ("chr1", "-") in trees

    def test_rep_cs_uses_name_field(self, tmp_path):
        bed = str(tmp_path / "atlas.bed")
        _write_bed(bed, [["chr1", 1000, 1005, "my_site", 40.0, "+"]])
        cs = ConstructSegments("", bed, 200)
        cs.process_pas_atlas()
        assert cs.pas_df.iloc[0]["rep_cs"] == "my_site"

    def test_rep_cs_fallback_for_dot_name(self, tmp_path):
        bed = str(tmp_path / "atlas.bed")
        _write_bed(bed, [["chr1", 1000, 1010, ".", 40.0, "+"]])
        cs = ConstructSegments("", bed, 200)
        cs.process_pas_atlas()
        rep_cs = cs.pas_df.iloc[0]["rep_cs"]
        assert rep_cs.startswith("chr1:")

    def test_mean_rpm_on_merge(self, tmp_path):
        bed = str(tmp_path / "atlas.bed")
        _write_bed(bed, [
            ["chr1", 1000, 1005, ".", 40.0, "+"],
            ["chr1", 1007, 1012, ".", 60.0, "+"],
        ])
        cs = ConstructSegments("", bed, 200)
        cs.process_pas_atlas(merge_distance=5)
        assert cs.pas_df.iloc[0]["atlas_rpm"] == pytest.approx(50.0)

    def test_highest_rpm_site_chosen_as_rep(self, tmp_path):
        bed = str(tmp_path / "atlas.bed")
        _write_bed(bed, [
            ["chr1", 1000, 1005, "cs_low", 20.0, "+"],
            ["chr1", 1007, 1012, "cs_high", 80.0, "+"],
        ])
        cs = ConstructSegments("", bed, 200)
        cs.process_pas_atlas(merge_distance=5)
        assert cs.pas_df.iloc[0]["rep_cs"] == "cs_high"

    def test_pas_df_columns(self, tmp_path):
        bed = str(tmp_path / "atlas.bed")
        _write_bed(bed, [["chr1", 1000, 1010, "cs1", 50.0, "+"]])
        cs = ConstructSegments("", bed, 200)
        cs.process_pas_atlas()
        assert set(cs.pas_df.columns) == {
            "pas_id", "chrom", "start", "end", "strand", "atlas_rpm", "rep_cs"
        }


# ---------------------------------------------------------------------------
# _assign_segment_ids
# ---------------------------------------------------------------------------


class TestAssignSegmentIds:
    def _gene_with_segments(self, origins: list[str]) -> Gene:
        gene = Gene("ENSG1")
        for i, origin in enumerate(origins):
            seg = Region(
                region_type="segment",
                chrom="chr1",
                start=i * 100,
                end=(i + 1) * 100,
                strand="+",
                attributes={
                    "segment_number": i + 1,
                    "segment_origin": origin,
                },
            )
            gene.segments.append(seg)
        return gene

    def test_exon_prefix(self):
        cs = _make_cs()
        cs.genes = {"ENSG1": self._gene_with_segments(["EX"])}
        seg_id = cs._assign_segment_ids()[("ENSG1", 1)]
        assert ":E" in seg_id

    def test_intron_prefix(self):
        cs = _make_cs()
        cs.genes = {"ENSG1": self._gene_with_segments(["IN"])}
        seg_id = cs._assign_segment_ids()[("ENSG1", 1)]
        assert ":I" in seg_id

    def test_terminal_exon_prefix(self):
        cs = _make_cs()
        cs.genes = {"ENSG1": self._gene_with_segments(["TE"])}
        seg_id = cs._assign_segment_ids()[("ENSG1", 1)]
        assert ":T" in seg_id

    def test_counter_increments_per_type(self):
        cs = _make_cs()
        cs.genes = {"ENSG1": self._gene_with_segments(["EX", "IN", "EX", "TE"])}
        seg_id_map = cs._assign_segment_ids()
        ids = list(seg_id_map.values())
        e_ids = [i for i in ids if ":E" in i]
        assert len(e_ids) == 2
        assert any("E001" in i for i in e_ids)
        assert any("E002" in i for i in e_ids)

    def test_three_digit_zero_padded(self):
        cs = _make_cs()
        cs.genes = {"ENSG1": self._gene_with_segments(["EX"])}
        seg_id = cs._assign_segment_ids()[("ENSG1", 1)]
        assert "E001" in seg_id


# ---------------------------------------------------------------------------
# Full compute() on real test data
# ---------------------------------------------------------------------------


class TestComputeOnRealData:
    @pytest.fixture(scope="class")
    def computed_cs(self, test_gtf, test_atlas):
        cs = ConstructSegments(
            test_gtf, test_atlas, downstream_exon_extension=200
        )
        cs.compute(merge_distance=5)
        return cs

    def test_genes_parsed(self, computed_cs):
        assert len(computed_cs.genes) > 0

    def test_pas_df_populated(self, computed_cs):
        assert computed_cs.pas_df is not None
        assert len(computed_cs.pas_df) > 0

    def test_pas_df_required_columns(self, computed_cs):
        for col in (
            "pas_id", "chrom", "start", "end",
            "strand", "atlas_rpm", "rep_cs",
        ):
            assert col in computed_cs.pas_df.columns

    def test_genes_have_segments(self, computed_cs):
        genes_with_segments = [
            g for g in computed_cs.genes.values() if g.segments
        ]
        assert len(genes_with_segments) > 0

    def test_subsegments_dataframe_not_empty(self, computed_cs):
        df = computed_cs.create_subsegments_dataframe()
        assert not df.empty

    def test_subsegments_dataframe_columns(self, computed_cs):
        df = computed_cs.create_subsegments_dataframe()
        for col in (
            "gene_id", "segment_id", "subsegment_id",
            "chrom", "start", "end", "strand", "pas_id",
        ):
            assert col in df.columns

    def test_subsegment_ids_unique(self, computed_cs):
        df = computed_cs.create_subsegments_dataframe()
        assert df["subsegment_id"].nunique() == len(df)

    def test_segment_id_is_prefix_of_subsegment_id(self, computed_cs):
        df = computed_cs.create_subsegments_dataframe()
        bad = df[
            ~df.apply(
                lambda row: row["subsegment_id"].startswith(row["segment_id"]),
                axis=1,
            )
        ]
        assert len(bad) == 0

    def test_run_writes_tsv(self, tmp_path, test_gtf, test_atlas):
        cs = ConstructSegments(test_gtf, test_atlas, downstream_exon_extension=200)
        out_tsv = str(tmp_path / "segments.tsv")
        cs.run(out_tsv=out_tsv, merge_distance=5)

        assert os.path.exists(out_tsv)
        df = pd.read_csv(out_tsv, sep="\t")
        assert not df.empty
        for col in (
            "chrom", "start", "end", "subsegment_id",
            "strand", "overlapping_pas_id", "atlas_rpm",
        ):
            assert col in df.columns

    def test_run_without_tsv_still_populates(self, test_gtf, test_atlas):
        cs = ConstructSegments(test_gtf, test_atlas, downstream_exon_extension=200)
        cs.run(out_tsv=None, merge_distance=5)
        assert cs.pas_df is not None
        assert len(cs.genes) > 0
