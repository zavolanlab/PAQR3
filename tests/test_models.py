"""Unit tests for paqr3.models."""

import pytest

from paqr3.models import Gene, Region, Transcript


# ---------------------------------------------------------------------------
# Region
# ---------------------------------------------------------------------------


class TestRegion:
    def test_defaults(self):
        r = Region(
            region_type="exon", chrom="chr1", start=100, end=200, strand="+"
        )
        assert r.region_type == "exon"
        assert r.chrom == "chr1"
        assert r.start == 100
        assert r.end == 200
        assert r.strand == "+"
        assert r.exon_number is None
        assert r.intron_number is None
        assert r.attributes == {}
        assert r.subsegments == []

    def test_optional_fields(self):
        r = Region(
            region_type="exon",
            chrom="chr1",
            start=0,
            end=10,
            strand="-",
            exon_number=3,
            intron_number=2,
            attributes={"gene_id": "g1"},
        )
        assert r.exon_number == 3
        assert r.intron_number == 2
        assert r.attributes["gene_id"] == "g1"

    def test_to_gtf_format_fields(self):
        r = Region(
            region_type="exon",
            chrom="chr1",
            start=100,
            end=200,
            strand="+",
            attributes={"source": "HAVANA", "gene_id": "g1"},
        )
        fields = r.to_gtf_format().split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "HAVANA"
        assert fields[2] == "exon"
        assert fields[3] == "100"
        assert fields[4] == "200"
        assert fields[5] == "."
        assert fields[6] == "+"
        assert fields[7] == "."
        assert 'gene_id "g1"' in fields[8]
        assert fields[8].endswith(";")

    def test_to_gtf_format_unknown_source(self):
        r = Region(region_type="exon", chrom="chr1", start=1, end=10, strand="+")
        fields = r.to_gtf_format().split("\t")
        assert fields[1] == "unknown"

    def test_to_gtf_format_multiple_attributes(self):
        r = Region(
            region_type="exon",
            chrom="chr1",
            start=1,
            end=100,
            strand="+",
            attributes={
                "source": "TEST",
                "gene_id": "g1",
                "transcript_id": "tx1",
            },
        )
        attr_str = r.to_gtf_format().split("\t")[8]
        assert 'gene_id "g1"' in attr_str
        assert 'transcript_id "tx1"' in attr_str

    def test_str_contains_key_info(self):
        r = Region(
            region_type="intron", chrom="chr2", start=500, end=600, strand="-"
        )
        s = str(r)
        assert "intron" in s
        assert "chr2" in s
        assert "500" in s
        assert "600" in s
        assert "-" in s

    def test_repr_equals_str(self):
        r = Region(region_type="exon", chrom="chr1", start=0, end=10, strand="+")
        assert repr(r) == str(r)

    def test_subsegments_list_is_independent(self):
        r1 = Region(region_type="exon", chrom="chr1", start=0, end=10, strand="+")
        r2 = Region(region_type="exon", chrom="chr1", start=0, end=10, strand="+")
        r1.subsegments.append(r2)
        # r2's subsegment list must not be shared with r1's
        assert r2.subsegments == []


# ---------------------------------------------------------------------------
# Transcript
# ---------------------------------------------------------------------------


def _exon(start: int, end: int, strand: str = "+") -> Region:
    return Region(
        region_type="exon", chrom="chr1", start=start, end=end, strand=strand
    )


class TestTranscript:
    def test_defaults(self):
        t = Transcript("tx1", "+")
        assert t.transcript_id == "tx1"
        assert t.strand == "+"
        assert t.regions == []
        assert t.attributes == {}

    def test_add_region(self):
        t = Transcript("tx1", "+")
        r = _exon(100, 200)
        t.add_region(r)
        assert len(t.regions) == 1
        assert t.regions[0] is r

    def test_add_multiple_regions(self):
        t = Transcript("tx1", "+")
        t.add_region(_exon(100, 200))
        t.add_region(_exon(400, 600))
        assert len(t.regions) == 2

    def test_to_gtf_format_feature_type(self):
        t = Transcript(
            "tx1",
            "+",
            attributes={"source": "ENSEMBL", "gene_id": "g1", "transcript_id": "tx1"},
        )
        t.add_region(_exon(100, 200))
        fields = t.to_gtf_format().split("\t")
        assert fields[2] == "transcript"

    def test_to_gtf_format_coords_span_all_regions(self):
        t = Transcript(
            "tx1",
            "+",
            attributes={"source": "ENSEMBL"},
        )
        t.add_region(_exon(100, 200))
        t.add_region(_exon(400, 600))
        fields = t.to_gtf_format().split("\t")
        assert fields[3] == "100"
        assert fields[4] == "600"

    def test_to_gtf_format_strand(self):
        t = Transcript("tx1", "-", attributes={"source": "TEST"})
        t.add_region(_exon(10, 50, strand="-"))
        fields = t.to_gtf_format().split("\t")
        assert fields[6] == "-"

    def test_repr_contains_id(self):
        t = Transcript("tx_abc", "+")
        assert "tx_abc" in repr(t)


# ---------------------------------------------------------------------------
# Gene
# ---------------------------------------------------------------------------


def _transcript_with_exon(tx_id: str, start: int, end: int) -> Transcript:
    t = Transcript(tx_id, "+")
    t.add_region(_exon(start, end))
    return t


class TestGene:
    def test_defaults(self):
        g = Gene("g1")
        assert g.gene_id == "g1"
        assert g.transcripts == {}
        assert g.attributes == {}
        assert g.segments == []

    def test_add_transcript(self):
        g = Gene("g1")
        t = _transcript_with_exon("tx1", 100, 200)
        g.add_transcript(t)
        assert "tx1" in g.transcripts
        assert g.transcripts["tx1"] is t

    def test_add_multiple_transcripts(self):
        g = Gene("g1")
        g.add_transcript(_transcript_with_exon("tx1", 100, 200))
        g.add_transcript(_transcript_with_exon("tx2", 500, 800))
        assert len(g.transcripts) == 2

    def test_to_gtf_format_feature_type(self):
        g = Gene("g1", attributes={"source": "ENSEMBL", "gene_id": "g1"})
        g.add_transcript(_transcript_with_exon("tx1", 100, 200))
        fields = g.to_gtf_format().split("\t")
        assert fields[2] == "gene"

    def test_to_gtf_format_spans_all_transcripts(self):
        g = Gene("g1", attributes={"source": "ENSEMBL", "gene_id": "g1"})
        g.add_transcript(_transcript_with_exon("tx1", 100, 200))
        g.add_transcript(_transcript_with_exon("tx2", 500, 800))
        fields = g.to_gtf_format().split("\t")
        assert fields[3] == "100"
        assert fields[4] == "800"

    def test_to_gtf_format_uses_first_transcript_strand(self):
        g = Gene("g1", attributes={"source": "TEST", "gene_id": "g1"})
        t = Transcript("tx1", "-")
        t.add_region(_exon(100, 200, strand="-"))
        g.add_transcript(t)
        fields = g.to_gtf_format().split("\t")
        assert fields[6] == "-"

    def test_repr_contains_gene_id(self):
        g = Gene("ENSG00000001")
        assert "ENSG00000001" in repr(g)

    def test_segments_list_is_independent(self):
        g1 = Gene("g1")
        g2 = Gene("g2")
        r = Region(region_type="segment", chrom="chr1", start=0, end=10, strand="+")
        g1.segments.append(r)
        assert g2.segments == []
