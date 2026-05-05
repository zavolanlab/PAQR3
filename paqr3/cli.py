"""Command-line interface for PAQR3."""

import argparse
import logging
import sys
from pathlib import Path

from paqr3.version import __version__
from paqr3.paqr3 import PAQR3

logging.basicConfig(
    format="[{asctime}] {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)
logging.getLogger().handlers.clear()
_handler = logging.StreamHandler()
_handler.setFormatter(
    logging.Formatter(
        "[{asctime}] {message}", style="{", datefmt="%Y-%m-%d %H:%M:%S"
    )
)
logging.getLogger().addHandler(_handler)
logging.getLogger().setLevel(logging.INFO)

logger = logging.getLogger(__name__)


def _require_file(path: str, flag: str) -> None:
    if not Path(path).is_file():
        print(f"Error: {flag} file not found: {path}", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Shared argument groups
# ---------------------------------------------------------------------------


_EMIT_CHOICES = [
    "mean_cov",
    "observed",
    "posterior",
    "atlas",
    "all",
    "debug",
]


def _add_common_args(p: argparse.ArgumentParser) -> None:
    """Arguments required by every sub-command."""
    p.add_argument(
        "--output-dir",
        "-o",
        required=True,
        help="Path to the output directory.",
    )
    p.add_argument(
        "--verbosity",
        "-v",
        choices=["INFO", "DEBUG"],
        default="INFO",
        help="Log verbosity level (default: INFO).",
    )
    p.add_argument(
        "--emit",
        "-e",
        nargs="+",
        choices=_EMIT_CHOICES,
        default=None,
        metavar="MODE",
        help=(
            "Additional outputs to emit.  Choices: "
            + ", ".join(_EMIT_CHOICES)
            + ".  'observed' → observed_usage.bw + observed_rpm.bw;"
            " 'posterior' → posterior_rpm.bw + posterior_usage.bw;"
            " 'atlas' → atlas_rpm.bw + atlas_usage.bw;"
            " 'mean_cov' → mean_cov.bw (subsegment-level);"
            " 'all' → all four BigWig groups;"
            " 'debug' → all BigWigs + JSON + debug BEDs."
            " Default: none (only the three results TSVs are written)."
        ),
    )
    p.add_argument(
        "--gzip",
        "-g",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Gzip-compress output TSV files (default: true). "
            "Use --no-gzip to write plain TSVs."
        ),
    )


def _add_segment_args(p: argparse.ArgumentParser) -> None:
    """Add arguments specific to the segmentation stage."""
    p.add_argument(
        "--annotation",
        "-a",
        required=True,
        help="Path to the annotation GTF file.",
    )
    p.add_argument(
        "--pas-atlas",
        "-pa",
        required=True,
        help="Path to the PAS atlas BED file.",
    )
    p.add_argument(
        "--downstream-exon-extension",
        "-de",
        type=int,
        default=200,
        help="Bases to extend terminal exons downstream (default: 200).",
    )
    p.add_argument(
        "--merge-distance",
        "-md",
        type=int,
        default=5,
        help="Merge PAS sites within this distance in bp (default: 5; 0 to disable).",
    )


def _add_quant_args(p: argparse.ArgumentParser) -> None:
    """Add arguments specific to the quantification stage."""
    p.add_argument(
        "--coverage-pos",
        "-c-pos",
        required=True,
        help="Path to the positive-strand coverage BigWig file.",
    )
    p.add_argument(
        "--coverage-neg",
        "-c-neg",
        required=True,
        help="Path to the negative-strand coverage BigWig file.",
    )
    p.add_argument(
        "--sample-id",
        "-sid",
        default=None,
        help=(
            "Sample identifier used for output filenames and directory. "
            "Defaults to the BigWig filename stem (before the first '.')."
        ),
    )
    p.add_argument(
        "--max-pas-count",
        "-mpc",
        type=int,
        default=10,
        help="Skip segments with more PAS than this (default: 10).",
    )
    p.add_argument(
        "--bam",
        "-b",
        default=None,
        help="Path to the aligned BAM file for expression rank statistics (optional).",
    )
    p.add_argument(
        "--f-stat-threshold",
        "-fst",
        type=int,
        default=100,
        help="F-statistic threshold for PAS usage (default: 100).",
    )
    p.add_argument(
        "--posterior-usage-weight",
        "-puw",
        type=float,
        default=0.1,
        help="Weight for observed RPM in posterior blending (default: 0.1).",
    )
    p.add_argument(
        "--threads",
        "-t",
        type=int,
        default=1,
        help="Threads for BigWig reading and F-stat worker processes (default: 1).",
    )
    p.add_argument(
        "--chr-sizes",
        "-cs",
        default=None,
        help=(
            "Path to a two-column chromosome-sizes file (chrom TAB size). "
            "Required when any BigWig emit mode is requested "
            "(final, mean_cov, rna_u, obs_rpm, post_rpm, all, debug)."
        ),
    )


# ---------------------------------------------------------------------------
# Sub-command handlers
# ---------------------------------------------------------------------------


_BW_EMIT_TOKENS = frozenset(
    {"mean_cov", "observed", "posterior", "atlas", "all", "debug"}
)


def _set_verbosity(args: argparse.Namespace) -> None:
    """Apply the --verbosity choice to the root logger."""
    level = logging.DEBUG if args.verbosity == "DEBUG" else logging.INFO
    logging.getLogger().setLevel(level)


def _require_chr_sizes(args: argparse.Namespace) -> None:
    """Exit with an error if BigWig modes are requested without --chr-sizes."""
    emit = set(args.emit or [])
    if emit & _BW_EMIT_TOKENS:
        if not args.chr_sizes:
            print(
                "Error: --chr-sizes is required when BigWig emit modes "
                f"are requested ({', '.join(sorted(emit & _BW_EMIT_TOKENS))}).",
                file=sys.stderr,
            )
            sys.exit(1)
        _require_file(args.chr_sizes, "--chr-sizes")


def _run_segment(args: argparse.Namespace) -> None:
    """Handle the segment sub-command."""
    _set_verbosity(args)
    _require_file(args.annotation, "--annotation")
    _require_file(args.pas_atlas, "--pas-atlas")
    logger.info("Starting PAQR3 segmentation...")
    paqr3 = PAQR3(
        annotation_file=args.annotation,
        pas_atlas_file=args.pas_atlas,
        coverage_bw_pos="",  # not used in segment mode
        coverage_bw_neg="",
        output_dir=args.output_dir,
        downstream_exon_extension=args.downstream_exon_extension,
        merge_distance=args.merge_distance,
        max_pas_count=0,
        emit=args.emit,
        use_gzip=args.gzip,
    )
    paqr3.run_segment()


def _run_quant(args: argparse.Namespace) -> None:
    """Handle the quant sub-command."""
    _set_verbosity(args)
    _require_chr_sizes(args)
    _require_file(args.segments_tsv, "--segments-tsv")
    _require_file(args.coverage_pos, "--coverage-pos")
    _require_file(args.coverage_neg, "--coverage-neg")
    logger.info("Starting PAQR3 quantification...")
    paqr3 = PAQR3(
        annotation_file="",  # not used in quant mode
        pas_atlas_file="",
        coverage_bw_pos=args.coverage_pos,
        coverage_bw_neg=args.coverage_neg,
        output_dir=args.output_dir,
        downstream_exon_extension=0,
        merge_distance=0,
        max_pas_count=args.max_pas_count,
        bam_file=args.bam,
        f_stat_threshold=args.f_stat_threshold,
        posterior_usage_weight=args.posterior_usage_weight,
        n_threads=args.threads,
        emit=args.emit,
        chr_sizes_file=args.chr_sizes,
        use_gzip=args.gzip,
    )
    paqr3.run_quant(args.segments_tsv, sample_name=args.sample_id)


def _run_full(args: argparse.Namespace) -> None:
    """Handle the full sub-command."""
    _set_verbosity(args)
    _require_chr_sizes(args)
    _require_file(args.annotation, "--annotation")
    _require_file(args.pas_atlas, "--pas-atlas")
    _require_file(args.coverage_pos, "--coverage-pos")
    _require_file(args.coverage_neg, "--coverage-neg")
    logger.info("Starting PAQR3 full pipeline...")
    paqr3 = PAQR3(
        annotation_file=args.annotation,
        pas_atlas_file=args.pas_atlas,
        coverage_bw_pos=args.coverage_pos,
        coverage_bw_neg=args.coverage_neg,
        output_dir=args.output_dir,
        downstream_exon_extension=args.downstream_exon_extension,
        merge_distance=args.merge_distance,
        max_pas_count=args.max_pas_count,
        bam_file=args.bam,
        f_stat_threshold=args.f_stat_threshold,
        posterior_usage_weight=args.posterior_usage_weight,
        n_threads=args.threads,
        emit=args.emit,
        chr_sizes_file=args.chr_sizes,
        use_gzip=args.gzip,
    )
    paqr3.run_full(sample_name=args.sample_id)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Entry point for the ``paqr3`` CLI."""
    parser = argparse.ArgumentParser(
        prog="paqr3",
        description=(
            "PAQR3: Poly(A) site quantification on standard RNA-Seq data."
        ),
    )
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=(
            f"PAQR3 v{__version__}, "
            "(c) 2025 by Zavolab (zavolab-biozentrum@unibas.ch)"
        ),
    )

    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    # ---- segment ----
    p_seg = sub.add_parser(
        "segment",
        help="Run annotation-based segmentation only.",
        description=(
            "Parse a GTF annotation, extend terminal exons, merge PAS "
            "sites and build sub-segments.  Produces a segments TSV "
            "(plus a debug BED when --emit debug is set)."
        ),
    )
    _add_common_args(p_seg)
    _add_segment_args(p_seg)
    p_seg.set_defaults(func=_run_segment)

    # ---- quant ----
    p_quant = sub.add_parser(
        "quant",
        help="Run quantification from a segments TSV.",
        description=(
            "Read a segments TSV produced by 'paqr3 segment', compute "
            "coverage metrics, evaluate PAS usage models and output "
            "posterior usage and gene-level summaries."
        ),
    )
    _add_common_args(p_quant)
    p_quant.add_argument(
        "--segments-tsv",
        "-s",
        required=True,
        help="Path to the segments TSV produced by 'paqr3 segment'.",
    )
    _add_quant_args(p_quant)
    p_quant.set_defaults(func=_run_quant)

    # ---- full ----
    p_full = sub.add_parser(
        "full",
        help="Run segmentation and quantification end-to-end.",
        description=(
            "Run the complete PAQR3 pipeline: segmentation followed by "
            "quantification.  Equivalent to 'paqr3 segment' then "
            "'paqr3 quant'."
        ),
    )
    _add_common_args(p_full)
    _add_segment_args(p_full)
    _add_quant_args(p_full)
    p_full.set_defaults(func=_run_full)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
