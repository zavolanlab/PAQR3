"""PAQR3 pipeline orchestration.

Exposes three run modes:

- :meth:`PAQR3.run_segment`: annotation-based segmentation only;
  produces the segments TSV (and optionally a debug BED).
- :meth:`PAQR3.run_quant`: quantification only; reads the segments TSV
  and produces PAS usage tables.
- :meth:`PAQR3.run_full`: runs both stages end-to-end.
"""

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
    {"final", "mean_cov", "rna_u", "obs_rpm", "post_rpm"}
)

# Maps emit token → (source DataFrame key, column name).
# Source keys: "raw" = raw_cov_df, "post" = posterior_df.
_EMIT_BW_SPEC: dict[str, tuple[str, str]] = {
    "mean_cov": ("raw", "mean_cov"),
    "rna_u": ("post", "rna_usage"),
    "obs_rpm": ("post", "observed_rpm"),
    "post_rpm": ("post", "posterior_rpm"),
    "final": ("post", "posterior_rel_usage"),
}


def _resolve_emit(emit: list[str] | None) -> set[str]:
    """Expand shortcut emit tokens into the full set of requested outputs.

    ``all``   expands to all BigWig modes.
    ``debug`` expands to all BigWig modes plus ``debug_json`` and
              ``debug`` (the last flag is also used to trigger debug
              BEDs in segment mode).
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


def _write_tsv_gz(df: pd.DataFrame, path: str, n_threads: int = 1) -> None:
    """Write *df* to a gzip-compressed TSV, using pigz when available."""
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
    logger.info("Written: %s", path)


def _write_bigwig(
    df: pd.DataFrame,
    value_col: str,
    chrom_sizes: dict[str, int],
    path: str,
) -> None:
    """Write a BigWig file with one constant-value entry per subsegment.

    Args:
        df: DataFrame with columns ``chrom``, ``start``, ``end``, and
            *value_col*.
        value_col: Column in *df* to use as the BigWig signal value.
        chrom_sizes: ``{chrom: size}`` dict used for the BigWig header.
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
        bam_file: str | None = None,
        f_stat_threshold: int = 100,
        posterior_usage_weight: float = 0.1,
        n_threads: int = 1,
        emit: list[str] | None = None,
        chr_sizes_file: str | None = None,
    ) -> None:
        self.annotation_file = annotation_file
        self.pas_atlas_file = pas_atlas_file
        self.coverage_bw_pos = coverage_bw_pos
        self.coverage_bw_neg = coverage_bw_neg
        self.output_dir = output_dir
        self.downstream_exon_extension = downstream_exon_extension
        self.merge_distance = merge_distance
        self.max_pas_count = max_pas_count
        self.bam_file = bam_file
        self.f_stat_threshold = f_stat_threshold
        self.posterior_usage_weight = posterior_usage_weight
        self.n_threads = n_threads
        self.emit = emit
        self.chr_sizes_file = chr_sizes_file

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
    ) -> None:
        """Write BigWig files for all requested emit modes.

        Args:
            effective_emit: Resolved set of emit tokens.
            results_dir: Output directory for this sample.
            sname: Sample name prefix for file names.
            raw_cov_df: Raw per-subsegment coverage metrics (all
                subsegments, including non-PAS; no coordinate columns).
            posterior_df: Per-PAS posterior usage DataFrame.
            subsegments_df: Sub-segment coordinate table used to attach
                ``chrom``/``start``/``end`` to PAS-only DataFrames.
        """
        bw_modes = effective_emit & _EMIT_BW_MODES
        if not bw_modes:
            return

        assert self.chr_sizes_file is not None
        chrom_sizes = _load_chrom_sizes(self.chr_sizes_file)

        coords = subsegments_df[
            ["subsegment_id", "chrom", "start", "end"]
        ].drop_duplicates("subsegment_id")

        for mode in bw_modes:
            source_key, col = _EMIT_BW_SPEC[mode]
            out_path = os.path.join(results_dir, f"{sname}_{mode}.bw")
            source_df = raw_cov_df if source_key == "raw" else posterior_df
            merged = source_df.merge(coords, on="subsegment_id", how="left")
            _write_bigwig(merged, col, chrom_sizes, out_path)

    def run_segment(self, output_dir: str | None = None) -> str:
        """Run the segmentation stage and write the segments TSV.

        Produces:

        - ``output_segments.tsv`` — primary output; input for
          :meth:`run_quant`.
        - ``output_segments_debug.bed`` — debug BED (only when
          ``debug`` is in ``emit``).

        Args:
            output_dir: Override the instance output directory.

        Returns:
            Path to the segments TSV file.
        """
        out_dir = output_dir or self.output_dir
        os.makedirs(out_dir, exist_ok=True)
        segments_tsv = os.path.join(out_dir, "output_segments.tsv")
        effective_emit = _resolve_emit(self.emit)
        debug_bed = (
            os.path.join(out_dir, "output_segments_debug.bed")
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
            out_debug_bed=debug_bed,
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

        Always writes ``{sname}_posterior_usage.tsv.gz``.  Additional
        outputs are controlled by ``self.emit``:

        - ``segment``  — ``{sname}_segment.tsv.gz``.
        - ``debug``    — ``{sname}_debug.json``.
        - ``mean_cov`` — ``{sname}_mean_cov.bw`` (mean BigWig coverage
          per subsegment, all subsegments including non-PAS).
        - ``rna_u``    — ``{sname}_rna_u.bw`` (F-stat usage fraction).
        - ``obs_rpm``  — ``{sname}_obs_rpm.bw``.
        - ``post_rpm`` — ``{sname}_post_rpm.bw``.
        - ``final``    — ``{sname}_final.bw`` (posterior relative usage).
        - ``all``      — all BigWig modes above.

        Args:
            segments_tsv: Path to the segments TSV from the
                segmentation stage.
            output_dir: Override the instance output directory.
            sample_name: Override the sample name derived from the
                BigWig filename.
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

        posterior_tsv = os.path.join(
            results_dir, f"{sname}_posterior_usage.tsv.gz"
        )
        debug_json = (
            os.path.join(results_dir, f"{sname}_debug.json")
            if "debug_json" in effective_emit
            else None
        )
        segment_tsv = (
            os.path.join(results_dir, f"{sname}_segment.tsv.gz")
            if "segment" in effective_emit
            else None
        )

        # Load the segments TSV; derive segment_id and gene_id from
        # the subsegment_id field ({gene_id}:{TypeNNN}:{sub:03d}).
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
            bam_file=self.bam_file,
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

        if segment_tsv is not None:
            segment_df = (
                usage_df[
                    [
                        "gene_id",
                        "segment_id",
                        "rna_sum_drop_cov",
                        "f_stat",
                        "p_value",
                    ]
                ]
                .drop_duplicates("segment_id")
                .reset_index(drop=True)
            )
            _write_tsv_gz(segment_df, segment_tsv, n_threads=self.n_threads)

        cpu = CalculatePosteriorUsage(weight=self.posterior_usage_weight)
        posterior_df = cpu.compute(usage_df)

        gene_usage_df = CalculateGeneLevelUsage(posterior_df).compute()

        # Merge gene_weighted_usage into posterior output; drop redundant
        # gene_id and segment_id columns (subsegment_id encodes both).
        posterior_out = posterior_df.merge(
            gene_usage_df, on="gene_id", how="left"
        ).drop(columns=["gene_id", "segment_id"])
        _write_tsv_gz(posterior_out, posterior_tsv, n_threads=self.n_threads)

        self._emit_bigwigs(
            effective_emit,
            results_dir,
            sname,
            raw_cov_df,
            posterior_df,
            subsegments_df,
        )

        logger.info("Quantification complete.")

    def run_full(self, sample_name: str | None = None) -> None:
        """Run segmentation and quantification end-to-end.

        Args:
            sample_name: Override the sample name derived from the
                BigWig filename.
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

        segments_tsv = os.path.join(results_dir, f"{sname}_segments.tsv")
        debug_bed = (
            os.path.join(results_dir, f"{sname}_debug.bed")
            if "debug" in effective_emit
            else None
        )

        # Stage 1: segmentation.
        cs = ConstructSegments(
            self.annotation_file,
            self.pas_atlas_file,
            self.downstream_exon_extension,
        )
        cs.run(
            segments_tsv,
            merge_distance=self.merge_distance,
            out_debug_bed=debug_bed,
        )
        assert cs.pas_df is not None, "ConstructSegments.run() must set pas_df"

        # Stage 2: use in-memory DataFrame (avoids re-reading the TSV).
        subsegments_df = cs.create_subsegments_dataframe()

        posterior_tsv = os.path.join(
            results_dir, f"{sname}_posterior_usage.tsv.gz"
        )
        debug_json = (
            os.path.join(results_dir, f"{sname}_debug.json")
            if "debug_json" in effective_emit
            else None
        )
        segment_tsv = (
            os.path.join(results_dir, f"{sname}_segment.tsv.gz")
            if "segment" in effective_emit
            else None
        )

        cc = CalculateCoverages(
            self.coverage_bw_pos,
            self.coverage_bw_neg,
            bam_file=self.bam_file,
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

        # Map atlas_rpm from the in-memory pas_df (rep_cs → atlas_rpm).
        atlas_map = (
            cs.pas_df[["rep_cs", "atlas_rpm"]]
            .drop_duplicates("rep_cs")
            .set_index("rep_cs")["atlas_rpm"]
        )
        usage_df["atlas_rpm"] = usage_df["pas_id"].map(atlas_map).fillna(0.0)

        if segment_tsv is not None:
            segment_df = (
                usage_df[
                    [
                        "gene_id",
                        "segment_id",
                        "rna_sum_drop_cov",
                        "f_stat",
                        "p_value",
                    ]
                ]
                .drop_duplicates("segment_id")
                .reset_index(drop=True)
            )
            _write_tsv_gz(segment_df, segment_tsv, n_threads=self.n_threads)

        cpu = CalculatePosteriorUsage(weight=self.posterior_usage_weight)
        posterior_df = cpu.compute(usage_df)

        gene_usage_df = CalculateGeneLevelUsage(posterior_df).compute()

        # Merge gene_weighted_usage into posterior output; drop redundant
        # gene_id and segment_id columns (subsegment_id encodes both).
        posterior_out = posterior_df.merge(
            gene_usage_df, on="gene_id", how="left"
        ).drop(columns=["gene_id", "segment_id"])
        _write_tsv_gz(posterior_out, posterior_tsv, n_threads=self.n_threads)

        self._emit_bigwigs(
            effective_emit,
            results_dir,
            sname,
            raw_cov_df,
            posterior_df,
            subsegments_df,
        )

        logger.info("PAQR3 pipeline completed.")
