"""Unit tests for paqr3.paqr3 (PAQR3 class and helpers)."""

import gzip
import os
from unittest.mock import patch

import pandas as pd
import pyBigWig
import pytest

from paqr3.paqr3 import PAQR3, _load_chrom_sizes, _resolve_emit, _write_bigwig, _write_tsv


# ---------------------------------------------------------------------------
# _resolve_emit
# ---------------------------------------------------------------------------


class TestResolveEmit:
    def test_none_returns_empty_set(self):
        assert _resolve_emit(None) == set()

    def test_empty_list_returns_empty_set(self):
        assert _resolve_emit([]) == set()

    def test_single_mode_passed_through(self):
        assert _resolve_emit(["mean_cov"]) == {"mean_cov"}

    def test_multiple_modes(self):
        assert _resolve_emit(["mean_cov", "observed"]) == {"mean_cov", "observed"}

    def test_all_expands_to_bw_modes(self):
        result = _resolve_emit(["all"])
        assert result == {"mean_cov", "observed", "posterior", "atlas"}

    def test_debug_expands_to_all_plus_debug_tokens(self):
        result = _resolve_emit(["debug"])
        assert "mean_cov" in result
        assert "observed" in result
        assert "posterior" in result
        assert "atlas" in result
        assert "debug_json" in result
        assert "debug" in result

    def test_all_and_extra_token(self):
        result = _resolve_emit(["all", "debug_json"])
        assert "debug_json" in result
        assert "mean_cov" in result

    def test_debug_is_superset_of_all(self):
        all_result = _resolve_emit(["all"])
        debug_result = _resolve_emit(["debug"])
        assert all_result.issubset(debug_result)


# ---------------------------------------------------------------------------
# _load_chrom_sizes
# ---------------------------------------------------------------------------


class TestLoadChromSizes:
    def test_basic_parsing(self, tmp_path):
        path = str(tmp_path / "chrom.sizes")
        with open(path, "w") as fh:
            fh.write("chr1\t248956422\n")
            fh.write("chr2\t242193529\n")
        sizes = _load_chrom_sizes(path)
        assert sizes == {"chr1": 248956422, "chr2": 242193529}

    def test_blank_lines_skipped(self, tmp_path):
        path = str(tmp_path / "chrom.sizes")
        with open(path, "w") as fh:
            fh.write("chr1\t100\n\nchr2\t200\n")
        sizes = _load_chrom_sizes(path)
        assert len(sizes) == 2

    def test_returns_int_values(self, tmp_path):
        path = str(tmp_path / "chrom.sizes")
        with open(path, "w") as fh:
            fh.write("chr1\t1000\n")
        sizes = _load_chrom_sizes(path)
        assert isinstance(sizes["chr1"], int)


# ---------------------------------------------------------------------------
# _write_tsv
# ---------------------------------------------------------------------------


class TestWriteTsv:
    def _sample_df(self):
        return pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})

    def test_plain_tsv_written(self, tmp_path):
        path = str(tmp_path / "out.tsv")
        _write_tsv(self._sample_df(), path, compress=False)
        assert os.path.exists(path)
        df = pd.read_csv(path, sep="\t")
        assert list(df.columns) == ["a", "b"]
        assert len(df) == 2

    def test_gzip_tsv_written(self, tmp_path):
        path = str(tmp_path / "out.tsv.gz")
        _write_tsv(self._sample_df(), path, compress=True)
        assert os.path.exists(path)
        with gzip.open(path, "rt") as fh:
            header = fh.readline().strip()
        assert header == "a\tb"

    def test_plain_tsv_correct_content(self, tmp_path):
        df = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})
        path = str(tmp_path / "out.tsv")
        _write_tsv(df, path, compress=False)
        result = pd.read_csv(path, sep="\t")
        pd.testing.assert_frame_equal(result, df)


# ---------------------------------------------------------------------------
# PAQR3 helper methods
# ---------------------------------------------------------------------------


class TestPAQR3Helpers:
    def _make_paqr3(self, bw_pos="sample.pos.bw", **kwargs):
        return PAQR3(
            annotation_file=kwargs.get("annotation_file", ""),
            pas_atlas_file=kwargs.get("pas_atlas_file", ""),
            coverage_bw_pos=bw_pos,
            coverage_bw_neg=kwargs.get("bw_neg", "sample.neg.bw"),
            output_dir=kwargs.get("output_dir", "/tmp"),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
        )

    def test_sample_name_from_bw_pos(self):
        paqr3 = self._make_paqr3("my_sample.positive.bw")
        assert paqr3._sample_name() == "my_sample"

    def test_sample_name_stops_at_first_dot(self):
        paqr3 = self._make_paqr3("SRR12345.chr1.bw")
        assert paqr3._sample_name() == "SRR12345"

    def test_make_results_dir_creates_directory(self, tmp_path):
        paqr3 = self._make_paqr3(output_dir=str(tmp_path))
        results_dir = paqr3._make_results_dir("my_sample")
        assert os.path.isdir(results_dir)
        assert results_dir.endswith("my_sample_results")


# ---------------------------------------------------------------------------
# PAQR3.run_segment on real test data
# ---------------------------------------------------------------------------


class TestRunSegment:
    def test_writes_segments_tsv(self, tmp_path, test_gtf, test_atlas):
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos="",
            coverage_bw_neg="",
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
        )
        out = paqr3.run_segment()
        assert os.path.exists(out)
        assert out.endswith("output_segments.tsv")

    def test_segments_tsv_has_required_columns(self, tmp_path, test_gtf, test_atlas):
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos="",
            coverage_bw_neg="",
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
        )
        out = paqr3.run_segment()
        df = pd.read_csv(out, sep="\t")
        for col in (
            "chrom", "start", "end", "subsegment_id",
            "strand", "overlapping_pas_id", "atlas_rpm",
        ):
            assert col in df.columns

    def test_returns_path_string(self, tmp_path, test_gtf, test_atlas):
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos="",
            coverage_bw_neg="",
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
        )
        result = paqr3.run_segment()
        assert isinstance(result, str)


# ---------------------------------------------------------------------------
# PAQR3.run_quant on real test data
# ---------------------------------------------------------------------------


class TestRunQuant:
    @pytest.fixture(scope="class")
    def segments_tsv(self, tmp_path_factory, test_gtf, test_atlas):
        out_dir = tmp_path_factory.mktemp("seg")
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos="",
            coverage_bw_neg="",
            output_dir=str(out_dir),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
        )
        return paqr3.run_segment()

    def test_writes_three_output_tsvs(
        self, tmp_path, segments_tsv, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file="",
            pas_atlas_file="",
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=False,
        )
        paqr3.run_quant(segments_tsv, sample_name="test_sample")
        results_dir = tmp_path / "test_sample_results"
        assert (results_dir / "test_sample_segment_results.tsv").exists()
        assert (results_dir / "test_sample_subsegment_results.tsv").exists()
        assert (results_dir / "test_sample_pas_results.tsv").exists()

    def test_segment_results_columns(
        self, tmp_path, segments_tsv, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file="",
            pas_atlas_file="",
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=False,
        )
        paqr3.run_quant(segments_tsv, sample_name="test_sample")
        df = pd.read_csv(
            tmp_path / "test_sample_results" / "test_sample_segment_results.tsv",
            sep="\t",
        )
        for col in ("chr", "start", "end", "strand", "segment_id",
                    "rna_sum_drop_cov", "f_stat", "p_value"):
            assert col in df.columns

    def test_subsegment_results_columns(
        self, tmp_path, segments_tsv, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file="",
            pas_atlas_file="",
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=False,
        )
        paqr3.run_quant(segments_tsv, sample_name="test_sample")
        df = pd.read_csv(
            tmp_path / "test_sample_results" / "test_sample_subsegment_results.tsv",
            sep="\t",
        )
        for col in ("chr", "start", "end", "strand", "subsegment_id", "pas_id",
                    "mean_cov", "atlas_rpm", "observed_rpm", "posterior_rpm"):
            assert col in df.columns

    def test_pas_results_columns(
        self, tmp_path, segments_tsv, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file="",
            pas_atlas_file="",
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=False,
        )
        paqr3.run_quant(segments_tsv, sample_name="test_sample")
        df = pd.read_csv(
            tmp_path / "test_sample_results" / "test_sample_pas_results.tsv",
            sep="\t",
        )
        for col in ("chr", "start", "end", "strand", "pas_id", "subsegment_id",
                    "atlas_usage", "observed_usage", "posterior_usage",
                    "gene_level_usage"):
            assert col in df.columns

    def test_mismatched_bw_stems_raises(
        self, tmp_path, segments_tsv, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file="",
            pas_atlas_file="",
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg="different_name.bw",
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
        )
        with pytest.raises(ValueError, match="sample name"):
            paqr3.run_quant(segments_tsv)


# ---------------------------------------------------------------------------
# PAQR3.run_full on real test data
# ---------------------------------------------------------------------------


class TestRunFull:
    def test_writes_three_output_tsvs(
        self, tmp_path, test_gtf, test_atlas, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=False,
        )
        paqr3.run_full(sample_name="test_sample")
        results_dir = tmp_path / "test_sample_results"
        assert (results_dir / "test_sample_segment_results.tsv").exists()
        assert (results_dir / "test_sample_subsegment_results.tsv").exists()
        assert (results_dir / "test_sample_pas_results.tsv").exists()

    def test_no_intermediate_segments_tsv(
        self, tmp_path, test_gtf, test_atlas, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=False,
        )
        paqr3.run_full(sample_name="test_sample")
        # run_full must NOT write output_segments.tsv to output_dir
        assert not (tmp_path / "output_segments.tsv").exists()

    def test_gzip_output(
        self, tmp_path, test_gtf, test_atlas, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=True,
        )
        paqr3.run_full(sample_name="test_sample")
        results_dir = tmp_path / "test_sample_results"
        assert (results_dir / "test_sample_pas_results.tsv.gz").exists()

    def test_inferred_sample_name_from_stems(
        self, tmp_path, test_gtf, test_atlas, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=False,
        )
        paqr3.run_full()  # no sample_name → inferred from BigWig stem
        results_dir = tmp_path / "test_data_results"
        assert results_dir.exists()

    def test_mismatched_stems_raises(
        self, tmp_path, test_gtf, test_atlas, test_bw_pos
    ):
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg="different_name.bw",
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
        )
        with pytest.raises(ValueError, match="sample name"):
            paqr3.run_full()


# ---------------------------------------------------------------------------
# _write_tsv gzip fallback (lines 98-99)
# ---------------------------------------------------------------------------


class TestWriteTsvGzipFallback:
    def test_fallback_to_stdlib_gzip(self, tmp_path):
        df = pd.DataFrame({"x": [1], "y": ["a"]})
        path = str(tmp_path / "out.tsv.gz")
        with patch("paqr3.paqr3.shutil.which", return_value=None):
            _write_tsv(df, path, compress=True)
        with gzip.open(path, "rt") as fh:
            header = fh.readline().strip()
        assert header == "x\ty"


# ---------------------------------------------------------------------------
# _write_bigwig (lines 119-139)
# ---------------------------------------------------------------------------


class TestWriteBigwig:
    def test_empty_df_not_written(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "start": [100],
            "end": [200],
            "val": [float("nan")],
        })
        path = str(tmp_path / "out.bw")
        _write_bigwig(df, "val", {"chr1": 1000}, path)
        assert not os.path.exists(path)

    def test_nonempty_df_written(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "start": [100],
            "end": [200],
            "val": [5.0],
        })
        path = str(tmp_path / "out.bw")
        _write_bigwig(df, "val", {"chr1": 1000}, path)
        assert os.path.exists(path)
        bw = pyBigWig.open(path)
        assert bw.stats("chr1", 100, 200)[0] == pytest.approx(5.0)
        bw.close()


# ---------------------------------------------------------------------------
# _emit_bigwigs + run_quant inferred sample name + pas_coords=None
# (lines 211-232, 288, 465, 507)
# ---------------------------------------------------------------------------


class TestRunQuantExtra:
    @pytest.fixture(scope="class")
    def segments_tsv(self, tmp_path_factory, test_gtf, test_atlas):
        out_dir = tmp_path_factory.mktemp("seg2")
        paqr3 = PAQR3(
            annotation_file=test_gtf,
            pas_atlas_file=test_atlas,
            coverage_bw_pos="",
            coverage_bw_neg="",
            output_dir=str(out_dir),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
        )
        return paqr3.run_segment()

    def test_inferred_sample_name(
        self, tmp_path, segments_tsv, test_bw_pos, test_bw_neg
    ):
        paqr3 = PAQR3(
            annotation_file="",
            pas_atlas_file="",
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=False,
        )
        paqr3.run_quant(segments_tsv)  # no sample_name
        assert (tmp_path / "test_data_results").exists()

    def test_run_quant_without_pas_coords(
        self, tmp_path, segments_tsv, test_bw_pos, test_bw_neg
    ):
        seg_df = pd.read_csv(segments_tsv, sep="\t")
        no_coords = str(tmp_path / "no_coords.tsv")
        seg_df.drop(
            columns=["pas_start", "pas_end"], errors="ignore"
        ).to_csv(no_coords, sep="\t", index=False)

        paqr3 = PAQR3(
            annotation_file="",
            pas_atlas_file="",
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            use_gzip=False,
        )
        paqr3.run_quant(no_coords, sample_name="noseg_sample")
        results_dir = tmp_path / "noseg_sample_results"
        assert (results_dir / "noseg_sample_pas_results.tsv").exists()

    def _chr_sizes_path(self, tmp_path, bw_path):
        bw = pyBigWig.open(bw_path)
        chroms = bw.chroms()
        bw.close()
        path = str(tmp_path / "chrom.sizes")
        with open(path, "w") as fh:
            for chrom, size in chroms.items():
                fh.write(f"{chrom}\t{size}\n")
        return path

    def test_emit_bigwigs_mean_cov(
        self, tmp_path, segments_tsv, test_bw_pos, test_bw_neg
    ):
        chr_sizes_path = self._chr_sizes_path(tmp_path, test_bw_pos)
        paqr3 = PAQR3(
            annotation_file="",
            pas_atlas_file="",
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            emit=["mean_cov"],
            chr_sizes_file=chr_sizes_path,
            use_gzip=False,
        )
        paqr3.run_quant(segments_tsv, sample_name="bw_sample")
        results_dir = tmp_path / "bw_sample_results"
        assert (
            (results_dir / "bw_sample_mean_cov.pos.bw").exists()
            or (results_dir / "bw_sample_mean_cov.neg.bw").exists()
        )

    def test_emit_bigwigs_observed_uses_pas_coords(
        self, tmp_path, segments_tsv, test_bw_pos, test_bw_neg
    ):
        # emit="observed" uses coord_type="pas" → line 223 in _emit_bigwigs
        chr_sizes_path = self._chr_sizes_path(tmp_path, test_bw_pos)
        paqr3 = PAQR3(
            annotation_file="",
            pas_atlas_file="",
            coverage_bw_pos=test_bw_pos,
            coverage_bw_neg=test_bw_neg,
            output_dir=str(tmp_path),
            downstream_exon_extension=200,
            merge_distance=5,
            max_pas_count=10,
            emit=["observed"],
            chr_sizes_file=chr_sizes_path,
            use_gzip=False,
        )
        paqr3.run_quant(segments_tsv, sample_name="obs_sample")
        # BigWig may or may not exist depending on whether any PAS pass,
        # but the code path (line 223) must be executed without error.
