import argparse
from paqr3.version import __version__
from paqr3.paqr3 import PAQR3
from paqr3.construct_segments import log_message


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
        "--downstream_exon_extension",
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
        "--output_dir",
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

    log_message("Starting PAQR3 pipeline...")

    paqr3_instance = PAQR3(
        args.annotation,
        args.pas_atlas,
        args.coverage[0],
        args.coverage[1],
        args.output_dir,
        args.downstream_exon_extension,
        args.merge_distance,
        args.max_pas_count,
        args.bam,
        args.f_stat_threshold,
        args.posterior_usage_weight,
        args.debug,
    )
    paqr3_instance.run()


if __name__ == "__main__":
    main()
