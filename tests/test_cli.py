"""Unit tests for paqr3.cli."""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from paqr3.cli import main


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _run(args: list[str]) -> None:
    """Run the CLI with the given argument list."""
    with patch.object(sys, "argv", ["paqr3"] + args):
        main()


def _run_exits(args: list[str]) -> pytest.ExceptionInfo:
    """Assert the CLI exits (SystemExit) and return the ExceptionInfo."""
    with pytest.raises(SystemExit) as exc_info:
        _run(args)
    return exc_info


# ---------------------------------------------------------------------------
# Argument parsing — segment subcommand
# ---------------------------------------------------------------------------


class TestSegmentArgParsing:
    def test_missing_annotation_exits(self, test_atlas, tmp_path):
        exc = _run_exits([
            "segment",
            "--pas-atlas", test_atlas,
            "--output-dir", str(tmp_path),
        ])
        assert exc.value.code != 0

    def test_missing_atlas_exits(self, test_gtf, tmp_path):
        exc = _run_exits([
            "segment",
            "--annotation", test_gtf,
            "--output-dir", str(tmp_path),
        ])
        assert exc.value.code != 0

    def test_missing_output_dir_exits(self, test_gtf, test_atlas):
        exc = _run_exits([
            "segment",
            "--annotation", test_gtf,
            "--pas-atlas", test_atlas,
        ])
        assert exc.value.code != 0

    def test_nonexistent_annotation_exits(self, test_atlas, tmp_path):
        exc = _run_exits([
            "segment",
            "--annotation", "/nonexistent/file.gtf",
            "--pas-atlas", test_atlas,
            "--output-dir", str(tmp_path),
        ])
        assert exc.value.code != 0

    def test_nonexistent_atlas_exits(self, test_gtf, tmp_path):
        exc = _run_exits([
            "segment",
            "--annotation", test_gtf,
            "--pas-atlas", "/nonexistent/atlas.bed",
            "--output-dir", str(tmp_path),
        ])
        assert exc.value.code != 0

    def test_valid_args_calls_run_segment(self, test_gtf, test_atlas, tmp_path):
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "segment",
                "--annotation", test_gtf,
                "--pas-atlas", test_atlas,
                "--output-dir", str(tmp_path),
            ])
            mock_instance.run_segment.assert_called_once()

    def test_downstream_extension_passed(self, test_gtf, test_atlas, tmp_path):
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "segment",
                "--annotation", test_gtf,
                "--pas-atlas", test_atlas,
                "--output-dir", str(tmp_path),
                "--downstream-exon-extension", "500",
            ])
            _, kwargs = mock_cls.call_args
            assert kwargs["downstream_exon_extension"] == 500

    def test_merge_distance_passed(self, test_gtf, test_atlas, tmp_path):
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "segment",
                "--annotation", test_gtf,
                "--pas-atlas", test_atlas,
                "--output-dir", str(tmp_path),
                "--merge-distance", "10",
            ])
            _, kwargs = mock_cls.call_args
            assert kwargs["merge_distance"] == 10


# ---------------------------------------------------------------------------
# Argument parsing — quant subcommand
# ---------------------------------------------------------------------------


class TestQuantArgParsing:
    def _seg_tsv(self, tmp_path: Path) -> str:
        path = tmp_path / "output_segments.tsv"
        path.write_text("chrom\tstart\tend\n")
        return str(path)

    def test_missing_segments_tsv_exits(
        self, test_bw_pos, test_bw_neg, tmp_path
    ):
        exc = _run_exits([
            "quant",
            "--coverage-pos", test_bw_pos,
            "--coverage-neg", test_bw_neg,
            "--output-dir", str(tmp_path),
        ])
        assert exc.value.code != 0

    def test_missing_coverage_pos_exits(
        self, test_bw_neg, tmp_path
    ):
        seg = self._seg_tsv(tmp_path)
        exc = _run_exits([
            "quant",
            "--segments-tsv", seg,
            "--coverage-neg", test_bw_neg,
            "--output-dir", str(tmp_path),
        ])
        assert exc.value.code != 0

    def test_bw_emit_without_chr_sizes_exits(
        self, test_bw_pos, test_bw_neg, tmp_path
    ):
        seg = self._seg_tsv(tmp_path)
        exc = _run_exits([
            "quant",
            "--segments-tsv", seg,
            "--coverage-pos", test_bw_pos,
            "--coverage-neg", test_bw_neg,
            "--output-dir", str(tmp_path),
            "--emit", "mean_cov",
        ])
        assert exc.value.code != 0

    def test_valid_args_calls_run_quant(
        self, test_bw_pos, test_bw_neg, tmp_path
    ):
        seg = self._seg_tsv(tmp_path)
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "quant",
                "--segments-tsv", seg,
                "--coverage-pos", test_bw_pos,
                "--coverage-neg", test_bw_neg,
                "--output-dir", str(tmp_path),
            ])
            mock_instance.run_quant.assert_called_once()

    def test_sample_id_forwarded(self, test_bw_pos, test_bw_neg, tmp_path):
        seg = self._seg_tsv(tmp_path)
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "quant",
                "--segments-tsv", seg,
                "--coverage-pos", test_bw_pos,
                "--coverage-neg", test_bw_neg,
                "--output-dir", str(tmp_path),
                "--sample-id", "my_run",
            ])
            args, _ = mock_instance.run_quant.call_args
            assert "my_run" in args or any("my_run" == v for v in _.values())

    def test_f_stat_threshold_forwarded(
        self, test_bw_pos, test_bw_neg, tmp_path
    ):
        seg = self._seg_tsv(tmp_path)
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "quant",
                "--segments-tsv", seg,
                "--coverage-pos", test_bw_pos,
                "--coverage-neg", test_bw_neg,
                "--output-dir", str(tmp_path),
                "--f-stat-threshold", "50",
            ])
            _, kwargs = mock_cls.call_args
            assert kwargs["f_stat_threshold"] == 50

    def test_no_gzip_flag(self, test_bw_pos, test_bw_neg, tmp_path):
        seg = self._seg_tsv(tmp_path)
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "quant",
                "--segments-tsv", seg,
                "--coverage-pos", test_bw_pos,
                "--coverage-neg", test_bw_neg,
                "--output-dir", str(tmp_path),
                "--no-gzip",
            ])
            _, kwargs = mock_cls.call_args
            assert kwargs["use_gzip"] is False


# ---------------------------------------------------------------------------
# Argument parsing — full subcommand
# ---------------------------------------------------------------------------


class TestFullArgParsing:
    def test_valid_args_calls_run_full(
        self, test_gtf, test_atlas, test_bw_pos, test_bw_neg, tmp_path
    ):
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "full",
                "--annotation", test_gtf,
                "--pas-atlas", test_atlas,
                "--coverage-pos", test_bw_pos,
                "--coverage-neg", test_bw_neg,
                "--output-dir", str(tmp_path),
            ])
            mock_instance.run_full.assert_called_once()

    def test_posterior_weight_forwarded(
        self, test_gtf, test_atlas, test_bw_pos, test_bw_neg, tmp_path
    ):
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "full",
                "--annotation", test_gtf,
                "--pas-atlas", test_atlas,
                "--coverage-pos", test_bw_pos,
                "--coverage-neg", test_bw_neg,
                "--output-dir", str(tmp_path),
                "--posterior-usage-weight", "0.3",
            ])
            _, kwargs = mock_cls.call_args
            assert kwargs["posterior_usage_weight"] == pytest.approx(0.3)

    def test_threads_forwarded(
        self, test_gtf, test_atlas, test_bw_pos, test_bw_neg, tmp_path
    ):
        with patch("paqr3.cli.PAQR3") as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance
            _run([
                "full",
                "--annotation", test_gtf,
                "--pas-atlas", test_atlas,
                "--coverage-pos", test_bw_pos,
                "--coverage-neg", test_bw_neg,
                "--output-dir", str(tmp_path),
                "--threads", "4",
            ])
            _, kwargs = mock_cls.call_args
            assert kwargs["n_threads"] == 4

    def test_no_subcommand_exits(self):
        exc = _run_exits([])
        assert exc.value.code != 0

    def test_version_flag_exits_zero(self):
        exc = _run_exits(["--version"])
        assert exc.value.code == 0
