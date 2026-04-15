"""Command-line interface for PAQR3.

Three sub-commands are exposed:

- ``paqr3 segment`` — run only the segmentation stage.
- ``paqr3 quant``   — run only the quantification stage.
- ``paqr3 full``    — run both stages end-to-end.
"""

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


def _add_common_args(p: argparse.ArgumentParser) -> None:
    """Arguments required by every sub-command."""
    p.add_argument(
        "--output_dir",
        "-o",
        required=True,
        help="Path to the output directory.",
    )
    p.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug output files (default: off).",
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
        "--pas_atlas",
        "-pa",
        required=True,
        help="Path to the PAS atlas BED file.",
    )
    p.add_argument(
        "--downstream_exon_extension",
        "-de",
        type=int,
        default=200,
        help="Bases to extend terminal exons downstream (default: 200).",
    )
    p.add_argument(
        "--merge_distance",
        "-md",
        type=int,
        default=5,
        help="Merge PAS sites within this distance in bp (default: 5; 0 to disable).",
    )


def _add_quant_args(p: argparse.ArgumentParser) -> None:
    """Add arguments specific to the quantification stage."""
    p.add_argument(
        "--coverage",
        "-c",
        required=True,
        nargs=2,
        metavar=("POS_BW", "NEG_BW"),
        help="Paths to the positive- and negative-strand coverage BigWig files.",
    )
    p.add_argument(
        "--max_pas_count",
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
        "--f_stat_threshold",
        "-fst",
        type=int,
        default=100,
        help="F-statistic threshold for PAS usage (default: 100).",
    )
    p.add_argument(
        "--posterior_usage_weight",
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


# ---------------------------------------------------------------------------
# Sub-command handlers
# ---------------------------------------------------------------------------


def _run_segment(args: argparse.Namespace) -> None:
    _require_file(args.annotation, "--annotation")
    _require_file(args.pas_atlas, "--pas_atlas")
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
        debug=args.debug,
    )
    paqr3.run_segment()


def _run_quant(args: argparse.Namespace) -> None:
    _require_file(args.segments_tsv, "--segments_tsv")
    _require_file(args.coverage[0], "--coverage")
    _require_file(args.coverage[1], "--coverage")
    logger.info("Starting PAQR3 quantification...")
    paqr3 = PAQR3(
        annotation_file="",  # not used in quant mode
        pas_atlas_file="",
        coverage_bw_pos=args.coverage[0],
        coverage_bw_neg=args.coverage[1],
        output_dir=args.output_dir,
        downstream_exon_extension=0,
        merge_distance=0,
        max_pas_count=args.max_pas_count,
        bam_file=args.bam,
        f_stat_threshold=args.f_stat_threshold,
        posterior_usage_weight=args.posterior_usage_weight,
        n_threads=args.threads,
        debug=args.debug,
    )
    paqr3.run_quant(args.segments_tsv)


def _run_full(args: argparse.Namespace) -> None:
    _require_file(args.annotation, "--annotation")
    _require_file(args.pas_atlas, "--pas_atlas")
    _require_file(args.coverage[0], "--coverage")
    _require_file(args.coverage[1], "--coverage")
    logger.info("Starting PAQR3 full pipeline...")
    paqr3 = PAQR3(
        annotation_file=args.annotation,
        pas_atlas_file=args.pas_atlas,
        coverage_bw_pos=args.coverage[0],
        coverage_bw_neg=args.coverage[1],
        output_dir=args.output_dir,
        downstream_exon_extension=args.downstream_exon_extension,
        merge_distance=args.merge_distance,
        max_pas_count=args.max_pas_count,
        bam_file=args.bam,
        f_stat_threshold=args.f_stat_threshold,
        posterior_usage_weight=args.posterior_usage_weight,
        n_threads=args.threads,
        debug=args.debug,
    )
    paqr3.run_full()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> None:
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
            "(plus a debug BED when --debug is set)."
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
        "--segments_tsv",
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
