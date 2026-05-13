"""PAQR3 pipeline orchestration."""

import gzip
import logging
import os
import shutil
import subprocess

import pandas as pd  # type: ignore
import pyBigWig  # type: ignore

from paqr3.construct_segments import ConstructSegments
from paqr3.calculate_cov_metrics import CalculateCoverages
from paqr3.calculate_posterior_usage import CalculatePosteriorUsage
from paqr3.calculate_gene_level_usages import CalculateGeneLevelUsage

logger = logging.getLogger(__name__)

# Emit tokens that produce BigWig output.
_EMIT_BW_MODES: frozenset[str] = frozenset(
    {"mean_cov", "observed", "posterior", "atlas"}
)

# Maps emit token → list of (source key, column, coord type, file suffix).
# Source keys:  "raw"  = raw_cov_df,  "post" = posterior_df.
# Coord types:  "subseg" = subsegment intervals; "pas" = PAS cluster coords.
_EMIT_BW_SPEC: dict[str, list[tuple[str, str, str, str]]] = {
    "mean_cov": [
        ("raw", "mean_cov", "subseg", "mean_cov"),
    ],
    "observed": [
        ("post", "rna_usage", "pas", "observed_usage"),
        ("post", "observed_rpm", "pas", "observed_rpm"),
    ],
    "posterior": [
        ("post", "posterior_rpm", "pas", "posterior_rpm"),
        ("post", "posterior_rel_usage", "pas", "posterior_usage"),
    ],
    "atlas": [
        ("post", "atlas_rpm", "pas", "atlas_rpm"),
        ("post", "atlas_rel_usage", "pas", "atlas_usage"),
    ],
}


def _resolve_emit(emit: list[str] | None) -> set[str]:
    """Expand shortcut emit tokens into the full set of requested outputs.

    "all" expands to all BigWig modes; "debug" expands to all BigWig
    modes plus debug_json and debug (which also triggers debug BEDs).
    """
    if not emit:
        return set()
    resolved: set[str] = set()
    for token in emit:
        if token == "all":
            resolved |= _EMIT_BW_MODES
        elif token == "debug":
            resolved |= _EMIT_BW_MODES
            resolved.add("debug_json")
            resolved.add("debug")
        else:
            resolved.add(token)
    return resolved


def _load_chrom_sizes(path: str) -> dict[str, int]:
    """Parse a two-column chrom-sizes file into a {chrom: size} dict."""
    sizes: dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, size = line.split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def _write_tsv(
    df: pd.DataFrame,
    path: str,
    n_threads: int = 1,
    compress: bool = True,
) -> None:
    """Write df to a TSV, optionally gzip-compressed via pigz or gzip."""
    if compress:
        csv_bytes = df.to_csv(sep="\t", index=False).encode()
        if shutil.which("pigz"):
            with open(path, "wb") as fh:
                subprocess.run(
                    ["pigz", "-p", str(n_threads), "-c"],
                    input=csv_bytes,
                    stdout=fh,
                    check=True,
                )
        else:
            with gzip.open(path, "wb") as fh:
                fh.write(csv_bytes)
    else:
        df.to_csv(path, sep="\t", index=False)
    logger.info("Written: %s", path)


def _write_bigwig(
    df: pd.DataFrame,
    value_col: str,
    chrom_sizes: dict[str, int],
    path: str,
) -> None:
    """Write a BigWig with one constant-value entry per interval.

    Args:
        df: DataFrame with columns chrom, start, end, and value_col.
        value_col: Column to use as the signal value.
        chrom_sizes: {chrom: size} dict for the BigWig header.
        path: Output file path.
    """
    mask = df[value_col].notna() & df["chrom"].isin(chrom_sizes)
    df = df[mask].sort_values(["chrom", "start"]).reset_index(drop=True)

    if df.empty:
        logger.warning("No data for BigWig '%s' — file not written.", path)
        return

    bw = pyBigWig.open(path, "w")
    bw.addHeader(sorted(chrom_sizes.items()))

    for chrom, grp in df.groupby("chrom", sort=True):
        grp = grp.sort_values("start")
        bw.addEntries(
            chroms=[str(chrom)] * len(grp),
            starts=grp["start"].astype(int).tolist(),
            ends=grp["end"].astype(int).tolist(),
            values=grp[value_col].astype(float).tolist(),
        )

    bw.close()
    logger.info("Written: %s", path)


class PAQR3:
    """Orchestrate the PAQR3 segmentation and quantification pipeline."""

    def __init__(
        self,
        annotation_file: str,
        pas_atlas_file: str,
        coverage_bw_pos: str,
        coverage_bw_neg: str,
        output_dir: str,
        downstream_exon_extension: int,
        merge_distance: int,
        max_pas_count: int,
        f_stat_threshold: int = 100,
        posterior_usage_weight: float = 0.1,
        n_threads: int = 1,
        emit: list[str] | None = None,
        chr_sizes_file: str | None = None,
        use_gzip: bool = True,
    ) -> None:
        self.annotation_file = annotation_file
        self.pas_atlas_file = pas_atlas_file
        self.coverage_bw_pos = coverage_bw_pos
        self.coverage_bw_neg = coverage_bw_neg
        self.output_dir = output_dir
        self.downstream_exon_extension = downstream_exon_extension
        self.merge_distance = merge_distance
        self.max_pas_count = max_pas_count
        self.f_stat_threshold = f_stat_threshold
        self.posterior_usage_weight = posterior_usage_weight
        self.n_threads = n_threads
        self.emit = emit
        self.chr_sizes_file = chr_sizes_file
        self.gzip = use_gzip

    def _sample_name(self) -> str:
        return os.path.basename(self.coverage_bw_pos).split(".")[0]

    def _make_results_dir(self, sample_name: str) -> str:
        results_dir = os.path.join(self.output_dir, f"{sample_name}_results")
        os.makedirs(results_dir, exist_ok=True)
        return results_dir

    def _emit_bigwigs(
        self,
        effective_emit: set[str],
        results_dir: str,
        sname: str,
        raw_cov_df: pd.DataFrame,
        posterior_df: pd.DataFrame,
        subsegments_df: pd.DataFrame,
        pas_coords: pd.DataFrame | None = None,
    ) -> None:
        """Write BigWig files for all requested emit modes.

        Args:
            effective_emit: Resolved set of emit tokens.
            results_dir: Output directory for this sample.
            sname: Sample name prefix for file names.
            raw_cov_df: Raw per-subsegment coverage metrics.
            posterior_df: Per-PAS posterior usage DataFrame.
            subsegments_df: Sub-segment coordinate table.
            pas_coords: PAS cluster coordinates (pas_id, chrom, start,
                end) for PAS coord-type modes.
        """
        bw_modes = effective_emit & _EMIT_BW_MODES
        if not bw_modes:
            return

        assert self.chr_sizes_file is not None
        chrom_sizes = _load_chrom_sizes(self.chr_sizes_file)

        coords = subsegments_df[
            ["subsegment_id", "chrom", "strand", "start", "end"]
        ].drop_duplicates("subsegment_id")

        for mode in bw_modes:
            for source_key, col, coord_type, suffix in _EMIT_BW_SPEC[mode]:
                source_df = raw_cov_df if source_key == "raw" else posterior_df
                if coord_type == "pas" and pas_coords is not None:
                    df_with_coords = source_df.merge(
                        pas_coords[["pas_id", "chrom", "strand", "start", "end"]],
                        on="pas_id",
                        how="left",
                    )
                else:
                    df_with_coords = source_df.merge(
                        coords, on="subsegment_id", how="left"
                    )
                for strand_char, strand_tag in [("+", "pos"), ("-", "neg")]:
                    out_path = os.path.join(
                        results_dir, f"{sname}_{suffix}.{strand_tag}.bw"
                    )
                    _write_bigwig(
                        df_with_coords[
                            df_with_coords["strand"] == strand_char
                        ],
                        col,
                        chrom_sizes,
                        out_path,
                    )

    def run_segment(self, output_dir: str | None = None) -> str:
        """Run the segmentation stage and write the segments TSV.

        Args:
            output_dir: Override the instance output directory.

        Returns:
            Path to the written segments TSV file.
        """
        out_dir = output_dir or self.output_dir
        os.makedirs(out_dir, exist_ok=True)
        segments_tsv = os.path.join(out_dir, "output_segments.tsv")
        effective_emit = _resolve_emit(self.emit)
        debug_prefix = (
            os.path.join(out_dir, "output_segments_debug")
            if "debug" in effective_emit
            else None
        )
        cs = ConstructSegments(
            self.annotation_file,
            self.pas_atlas_file,
            self.downstream_exon_extension,
        )
        cs.run(
            segments_tsv,
            merge_distance=self.merge_distance,
            out_debug_prefix=debug_prefix,
        )
        logger.info("Segmentation complete. Output: %s", segments_tsv)
        return segments_tsv

    def run_quant(
        self,
        segments_tsv: str,
        output_dir: str | None = None,
        sample_name: str | None = None,
    ) -> None:
        """Run the quantification stage from a segments TSV.

        Args:
            segments_tsv: Path to the segments TSV from segmentation.
            output_dir: Override the instance output directory.
            sample_name: Override the sample name from the BigWig
                filename.
        """
        if sample_name is not None:
            sname = sample_name
        else:
            pos_stem = self._sample_name()
            neg_stem = os.path.basename(self.coverage_bw_neg).split(".")[0]
            if pos_stem != neg_stem:
                raise ValueError(
                    "Coverage files must share the same sample name."
                )
            sname = pos_stem
        results_dir = output_dir or self._make_results_dir(sname)
        os.makedirs(results_dir, exist_ok=True)

        effective_emit = _resolve_emit(self.emit)
        ext = ".tsv.gz" if self.gzip else ".tsv"

        debug_json = (
            os.path.join(results_dir, f"{sname}_debug.json")
            if "debug_json" in effective_emit
            else None
        )

        # Load segments TSV; derive segment_id and gene_id from
        # subsegment_id ({gene_id}:{TypeNNN}:{sub:03d}).
        seg_df = pd.read_csv(segments_tsv, sep="\t")
        seg_df["segment_id"] = (
            seg_df["subsegment_id"].str.rsplit(":", n=1).str[0]
        )
        seg_df["gene_id"] = seg_df["subsegment_id"].str.split(":").str[0]
        subsegments_df = seg_df.rename(
            columns={"overlapping_pas_id": "pas_id"}
        )

        cc = CalculateCoverages(
            self.coverage_bw_pos,
            self.coverage_bw_neg,
            f_stat_threshold=self.f_stat_threshold,
        )
        raw_cov_df, usage_df = cc.run(
            subsegments_df,
            output_debug_json=debug_json,
            n_threads=self.n_threads,
            max_pas_count=self.max_pas_count,
            f_stat_threshold=self.f_stat_threshold,
            debug=logger.isEnabledFor(logging.DEBUG),
        )

        # Map atlas_rpm via subsegment_id.
        atlas_map = (
            seg_df[["subsegment_id", "atlas_rpm"]]
            .drop_duplicates("subsegment_id")
            .set_index("subsegment_id")["atlas_rpm"]
        )
        usage_df["atlas_rpm"] = (
            usage_df["subsegment_id"].map(atlas_map).fillna(0.0)
        )

        cpu = CalculatePosteriorUsage(weight=self.posterior_usage_weight)
        posterior_df = cpu.compute(usage_df)

        # Segment coordinates (min start / max end over subsegments).
        seg_coords = (
            subsegments_df.groupby("segment_id", sort=False)
            .agg(
                chrom=("chrom", "first"),
                strand=("strand", "first"),
                start=("start", "min"),
                end=("end", "max"),
            )
            .reset_index()
        )

        gene_usage_df = CalculateGeneLevelUsage(posterior_df).compute()

        # _segment_results.tsv
        seg_stats = (
            usage_df[
                [
                    "segment_id",
                    "rna_sum_drop_cov",
                    "f_stat",
                    "p_value",
                ]
            ]
            .drop_duplicates("segment_id")
            .reset_index(drop=True)
        )
        segment_results = seg_coords.merge(
            seg_stats, on="segment_id", how="inner"
        )[
            [
                "chrom",
                "start",
                "end",
                "strand",
                "segment_id",
                "rna_sum_drop_cov",
                "f_stat",
                "p_value",
            ]
        ].rename(
            columns={"chrom": "chr"}
        )
        _write_tsv(
            segment_results,
            os.path.join(results_dir, f"{sname}_segment_results{ext}"),
            n_threads=self.n_threads,
            compress=self.gzip,
        )

        # _subsegment_results.tsv
        subseg_results = (
            subsegments_df[
                [
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "subsegment_id",
                    "pas_id",
                    "atlas_rpm",
                ]
            ]
            .merge(
                raw_cov_df[["subsegment_id", "mean_cov"]].drop_duplicates(
                    "subsegment_id"
                ),
                on="subsegment_id",
                how="left",
            )
            .merge(
                posterior_df[
                    [
                        "subsegment_id",
                        "observed_rpm",
                        "posterior_rpm",
                    ]
                ].drop_duplicates("subsegment_id"),
                on="subsegment_id",
                how="left",
            )[
                [
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "subsegment_id",
                    "pas_id",
                    "mean_cov",
                    "atlas_rpm",
                    "observed_rpm",
                    "posterior_rpm",
                ]
            ]
            .rename(columns={"chrom": "chr"})
        )
        _write_tsv(
            subseg_results,
            os.path.join(results_dir, f"{sname}_subsegment_results{ext}"),
            n_threads=self.n_threads,
            compress=self.gzip,
        )

        # PAS cluster coordinates (with strand) from segments TSV.
        if "pas_start" in seg_df.columns:
            pas_coords: pd.DataFrame | None = (
                seg_df[
                    [
                        "chrom",
                        "strand",
                        "overlapping_pas_id",
                        "pas_start",
                        "pas_end",
                    ]
                ]
                .query("overlapping_pas_id != '.'")
                .drop_duplicates("overlapping_pas_id")
                .rename(
                    columns={
                        "overlapping_pas_id": "pas_id",
                        "pas_start": "start",
                        "pas_end": "end",
                    }
                )
            )
        else:
            pas_coords = None

        # _pas_results.tsv
        pas_usage = posterior_df[
            [
                "pas_id",
                "subsegment_id",
                "rna_usage",
                "atlas_rel_usage",
                "posterior_rel_usage",
            ]
        ].rename(
            columns={
                "rna_usage": "observed_usage",
                "atlas_rel_usage": "atlas_usage",
                "posterior_rel_usage": "posterior_usage",
            }
        )
        if pas_coords is not None:
            pas_results = (
                pas_usage.merge(
                    pas_coords[["pas_id", "chrom", "start", "end", "strand"]],
                    on="pas_id",
                    how="left",
                )
                .merge(gene_usage_df, on="subsegment_id", how="left")[
                    [
                        "chrom",
                        "start",
                        "end",
                        "strand",
                        "pas_id",
                        "subsegment_id",
                        "atlas_usage",
                        "observed_usage",
                        "posterior_usage",
                        "gene_level_usage",
                    ]
                ]
                .rename(columns={"chrom": "chr"})
            )
        else:
            pas_results = pas_usage.merge(
                gene_usage_df, on="subsegment_id", how="left"
            )[
                [
                    "pas_id",
                    "subsegment_id",
                    "atlas_usage",
                    "observed_usage",
                    "posterior_usage",
                    "gene_level_usage",
                ]
            ]
        _write_tsv(
            pas_results,
            os.path.join(results_dir, f"{sname}_pas_results{ext}"),
            n_threads=self.n_threads,
            compress=self.gzip,
        )

        self._emit_bigwigs(
            effective_emit,
            results_dir,
            sname,
            raw_cov_df,
            posterior_df,
            subsegments_df,
            pas_coords=pas_coords,
        )

        logger.info("Quantification complete.")

    def run_full(self, sample_name: str | None = None) -> None:
        """Run segmentation and quantification end-to-end.

        Args:
            sample_name: Override the sample name from the BigWig
                filename.
        """
        if sample_name is not None:
            sname = sample_name
        else:
            sname = self._sample_name()
            neg_stem = os.path.basename(self.coverage_bw_neg).split(".")[0]
            if sname != neg_stem:
                raise ValueError(
                    "Coverage files must share the same sample name."
                )
        results_dir = self._make_results_dir(sname)

        effective_emit = _resolve_emit(self.emit)
        ext = ".tsv.gz" if self.gzip else ".tsv"

        debug_prefix = (
            os.path.join(results_dir, f"{sname}_debug")
            if "debug" in effective_emit
            else None
        )

        # Stage 1: segmentation (TSV skipped; in-memory DataFrame used).
        cs = ConstructSegments(
            self.annotation_file,
            self.pas_atlas_file,
            self.downstream_exon_extension,
        )
        cs.run(
            merge_distance=self.merge_distance,
            out_debug_prefix=debug_prefix,
        )
        assert cs.pas_df is not None, "ConstructSegments.run() must set pas_df"

        # Stage 2: in-memory DataFrame (avoids re-reading the TSV).
        subsegments_df = cs.create_subsegments_dataframe()

        debug_json = (
            os.path.join(results_dir, f"{sname}_debug.json")
            if "debug_json" in effective_emit
            else None
        )

        cc = CalculateCoverages(
            self.coverage_bw_pos,
            self.coverage_bw_neg,
            f_stat_threshold=self.f_stat_threshold,
        )
        raw_cov_df, usage_df = cc.run(
            subsegments_df,
            output_debug_json=debug_json,
            n_threads=self.n_threads,
            max_pas_count=self.max_pas_count,
            f_stat_threshold=self.f_stat_threshold,
            debug=logger.isEnabledFor(logging.DEBUG),
        )

        # Map atlas_rpm from in-memory pas_df (rep_cs → atlas_rpm).
        atlas_map = (
            cs.pas_df[["rep_cs", "atlas_rpm"]]
            .drop_duplicates("rep_cs")
            .set_index("rep_cs")["atlas_rpm"]
        )
        usage_df["atlas_rpm"] = usage_df["pas_id"].map(atlas_map).fillna(0.0)
        subsegments_df["atlas_rpm"] = (
            subsegments_df["pas_id"].map(atlas_map).fillna(0.0)
        )

        cpu = CalculatePosteriorUsage(weight=self.posterior_usage_weight)
        posterior_df = cpu.compute(usage_df)

        gene_usage_df = CalculateGeneLevelUsage(posterior_df).compute()

        # Segment coordinates (min start / max end over subsegments).
        seg_coords = (
            subsegments_df.groupby("segment_id", sort=False)
            .agg(
                chrom=("chrom", "first"),
                strand=("strand", "first"),
                start=("start", "min"),
                end=("end", "max"),
            )
            .reset_index()
        )

        # _segment_results.tsv
        seg_stats = (
            usage_df[
                [
                    "segment_id",
                    "rna_sum_drop_cov",
                    "f_stat",
                    "p_value",
                ]
            ]
            .drop_duplicates("segment_id")
            .reset_index(drop=True)
        )
        segment_results = seg_coords.merge(
            seg_stats, on="segment_id", how="inner"
        )[
            [
                "chrom",
                "start",
                "end",
                "strand",
                "segment_id",
                "rna_sum_drop_cov",
                "f_stat",
                "p_value",
            ]
        ].rename(
            columns={"chrom": "chr"}
        )
        _write_tsv(
            segment_results,
            os.path.join(results_dir, f"{sname}_segment_results{ext}"),
            n_threads=self.n_threads,
            compress=self.gzip,
        )

        # _subsegment_results.tsv
        subseg_results = (
            subsegments_df[
                [
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "subsegment_id",
                    "pas_id",
                    "atlas_rpm",
                ]
            ]
            .merge(
                raw_cov_df[["subsegment_id", "mean_cov"]].drop_duplicates(
                    "subsegment_id"
                ),
                on="subsegment_id",
                how="left",
            )
            .merge(
                posterior_df[
                    [
                        "subsegment_id",
                        "observed_rpm",
                        "posterior_rpm",
                    ]
                ].drop_duplicates("subsegment_id"),
                on="subsegment_id",
                how="left",
            )[
                [
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "subsegment_id",
                    "pas_id",
                    "mean_cov",
                    "atlas_rpm",
                    "observed_rpm",
                    "posterior_rpm",
                ]
            ]
            .rename(columns={"chrom": "chr"})
        )
        _write_tsv(
            subseg_results,
            os.path.join(results_dir, f"{sname}_subsegment_results{ext}"),
            n_threads=self.n_threads,
            compress=self.gzip,
        )

        # PAS cluster coordinates (with strand) from in-memory pas_df.
        pas_coords = (
            cs.pas_df[["rep_cs", "chrom", "strand", "start", "end"]]
            .rename(columns={"rep_cs": "pas_id"})
            .drop_duplicates("pas_id")
        )

        # _pas_results.tsv
        pas_usage = posterior_df[
            [
                "pas_id",
                "subsegment_id",
                "rna_usage",
                "atlas_rel_usage",
                "posterior_rel_usage",
            ]
        ].rename(
            columns={
                "rna_usage": "observed_usage",
                "atlas_rel_usage": "atlas_usage",
                "posterior_rel_usage": "posterior_usage",
            }
        )
        pas_results = (
            pas_usage.merge(
                pas_coords[["pas_id", "chrom", "start", "end", "strand"]],
                on="pas_id",
                how="left",
            )
            .merge(gene_usage_df, on="subsegment_id", how="left")[
                [
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "pas_id",
                    "subsegment_id",
                    "atlas_usage",
                    "observed_usage",
                    "posterior_usage",
                    "gene_level_usage",
                ]
            ]
            .rename(columns={"chrom": "chr"})
        )
        _write_tsv(
            pas_results,
            os.path.join(results_dir, f"{sname}_pas_results{ext}"),
            n_threads=self.n_threads,
            compress=self.gzip,
        )

        self._emit_bigwigs(
            effective_emit,
            results_dir,
            sname,
            raw_cov_df,
            posterior_df,
            subsegments_df,
            pas_coords=pas_coords,
        )

        logger.info("PAQR3 pipeline completed.")
