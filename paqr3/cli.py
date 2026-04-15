import argparse
import sys
import logging
from pathlib import Path
from paqr3.version import __version__
from paqr3.paqr3 import PAQR3

# ——————————— Logging setup ———————————
logging.basicConfig(
    format="[{asctime}] {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)
logging.getLogger().handlers.clear()
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "[{asctime}] {message}", style="{", datefmt="%Y-%m-%d %H:%M:%S"
)
handler.setFormatter(formatter)
logging.getLogger().addHandler(handler)
logging.getLogger().setLevel(logging.INFO)

logger = logging.getLogger(__name__)


def _require_file(path: str, flag: str) -> None:
    if not Path(path).is_file():
        print(f"Error: {flag} file not found: {path}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=("Construct genomic segments from RNA-Seq annotations.")
    )
    parser.add_argument(
        "--annotation",
        "-a",
        type=str,
        required=True,
        help="Path to the annotation GTF file.",
    )
    parser.add_argument(
        "--downstream-exon-extension",
        "-de",
        type=int,
        default=200,
        help="Number of bases to extend terminal exons.",
    )
    parser.add_argument(
        "--pas_atlas",
        "-pa",
        type=str,
        required=True,
        help="Path to the PAS atlas BED file.",
    )
    parser.add_argument(
        "--merge-distance",
        "-md",
        type=int,
        default=5,
        help="Merge PAS sites that are within this distance (in base pairs). Set to 0 to disable merging.",
    )
    parser.add_argument(
        "--coverage",
        "-c",
        type=str,
        required=True,
        nargs=2,
        help=(
            "Paths to the positive and negative strand coverage bigWig files"
            " (e.g., pos.bw neg.bw)."
        ),
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=str,
        required=True,
        help="Path to the output directory.",
    )
    parser.add_argument(
        "--max_pas_count",
        "-mpc",
        type=int,
        default=10,
        help="Only evaluate segments with this number or fewer PAS.",
    )
    parser.add_argument(
        "--bam",
        "-b",
        type=str,
        default=None,
        help="Path to the aligned RNA-Seq BAM file for segment read counts (optional).",
    )
    parser.add_argument(
        "--f_stat_threshold",
        "-fst",
        type=int,
        default=100,
        help="F-statistic threshold for PAS usage evaluation.",
    )
    parser.add_argument(
        "--posterior_usage_weight",
        "-puw",
        type=float,
        default=0.1,
        help="Weight to combine RNA-seq and atlas RPMs for posterior PAS usage (default: 0.1)",
    )
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=1,
        help=(
            "Number of parallel threads for BigWig reading and worker"
            " processes for F-statistic evaluation (default: 1)."
        ),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug output (default: off)",
    )
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=(
            f"PAQR3 v{__version__},"
            "(c) 2025 by Zavolab (zavolab-biozentrum@unibas.ch)"
        ),
        help="Show version information and exit",
    )

    args = parser.parse_args()

    _require_file(args.annotation, "--annotation")
    _require_file(args.pas_atlas, "--pas-atlas")
    _require_file(args.coverage[0], "--coverage")
    _require_file(args.coverage[1], "--coverage")

    logger.info("Starting PAQR3 pipeline... - new version")

    paqr3_instance = PAQR3(
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
    paqr3_instance.run()


if __name__ == "__main__":
    main()
