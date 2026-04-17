"""PAQR3 pipeline orchestration.

Exposes three run modes:

- :meth:`PAQR3.run_segment`: annotation-based segmentation only;
  produces the segments TSV (and optionally a debug BED).
- :meth:`PAQR3.run_quant`: quantification only; reads the segments TSV
  and produces PAS usage tables.
- :meth:`PAQR3.run_full`: runs both stages end-to-end.
"""

import logging
import os

import pandas as pd  # type: ignore

from paqr3.construct_segments import ConstructSegments
from paqr3.calculate_cov_metrics import CalculateCoverages
from paqr3.calculate_posterior_usage import CalculatePosteriorUsage
from paqr3.calculate_gene_level_usages import CalculateGeneLevelUsage

logger = logging.getLogger(__name__)


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
        debug: bool = False,
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
        self.debug = debug

    def _sample_name(self) -> str:
        return os.path.basename(self.coverage_bw_pos).split(".")[0]

    def _make_results_dir(self, sample_name: str) -> str:
        results_dir = os.path.join(
            self.output_dir, f"{sample_name}_results"
        )
        os.makedirs(results_dir, exist_ok=True)
        return results_dir

    def run_segment(self, output_dir: str | None = None) -> str:
        """Run the segmentation stage and write the segments TSV.

        Produces:

        - ``{sample}_segments.tsv`` — primary output; input for
          :meth:`run_quant`.
        - ``{sample}_debug.bed`` — debug BED (only when debug=True).

        Args:
            output_dir: Override the instance output directory.

        Returns:
            Path to the segments TSV file.
        """
        out_dir = output_dir or self.output_dir
        os.makedirs(out_dir, exist_ok=True)
        segments_tsv = os.path.join(out_dir, "output_segments.tsv")
        debug_bed = (
            os.path.join(out_dir, "output_segments_debug.bed")
            if self.debug
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

        Args:
            segments_tsv: Path to the segments TSV from the
                segmentation stage.
            output_dir: Override the instance output directory.
            sample_name: Override the sample name derived from the
                BigWig filename.
        """
        pos_stem = self._sample_name()
        neg_stem = os.path.basename(self.coverage_bw_neg).split(".")[0]
        if pos_stem != neg_stem:
            raise ValueError(
                "Coverage files must share the same sample name."
            )
        sname = sample_name or pos_stem
        results_dir = (
            output_dir or self._make_results_dir(sname)
        )
        os.makedirs(results_dir, exist_ok=True)

        usage_tsv = os.path.join(results_dir, f"{sname}_usage.tsv")
        posterior_tsv = os.path.join(
            results_dir, f"{sname}_posterior_usage.tsv"
        )
        coverage_tsv = os.path.join(results_dir, f"{sname}_coverage.tsv")
        debug_json = os.path.join(results_dir, f"{sname}_debug.json")

        # Load the segments TSV; derive segment_id and gene_id from
        # the subsegment_id field ({gene_id}:{TypeNNN}:{sub:03d}).
        seg_df = pd.read_csv(segments_tsv, sep="\t")
        seg_df["segment_id"] = seg_df["subsegment_id"].str.rsplit(
            ":", n=1
        ).str[0]
        seg_df["gene_id"] = seg_df["subsegment_id"].str.split(
            ":"
        ).str[0]
        subsegments_df = seg_df.rename(
            columns={"overlapping_pas_id": "pas_id"}
        )

        cc = CalculateCoverages(
            self.coverage_bw_pos,
            self.coverage_bw_neg,
            bam_file=self.bam_file,
            f_stat_threshold=self.f_stat_threshold,
        )
        usage_df = cc.run(
            subsegments_df,
            output_raw_tsv=coverage_tsv if self.debug else None,
            output_final_tsv=usage_tsv if self.debug else None,
            output_debug_json=debug_json if self.debug else None,
            n_threads=self.n_threads,
            max_pas_count=self.max_pas_count,
            f_stat_threshold=self.f_stat_threshold,
            debug=self.debug,
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

        gene_usage_df = CalculateGeneLevelUsage(posterior_df).compute()

        # Merge gene_weighted_usage into posterior output; drop redundant
        # gene_id and segment_id columns (subsegment_id encodes both).
        posterior_out = posterior_df.merge(
            gene_usage_df, on="gene_id", how="left"
        ).drop(columns=["gene_id", "segment_id"])
        posterior_out.to_csv(posterior_tsv, sep="\t", index=False)
        logger.info("Posterior usage written to %s", posterior_tsv)
        logger.info("Quantification complete.")

    def run_full(self) -> None:
        """Run segmentation and quantification end-to-end."""
        sname = self._sample_name()
        neg_stem = os.path.basename(self.coverage_bw_neg).split(".")[0]
        if sname != neg_stem:
            raise ValueError(
                "Coverage files must share the same sample name."
            )
        results_dir = self._make_results_dir(sname)

        segments_tsv = os.path.join(results_dir, f"{sname}_segments.tsv")
        debug_bed = (
            os.path.join(results_dir, f"{sname}_debug.bed")
            if self.debug
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

        # Stage 2: use in-memory DataFrame (avoids re-reading the TSV).
        subsegments_df = cs.create_subsegments_dataframe()

        usage_tsv = os.path.join(results_dir, f"{sname}_usage.tsv")
        posterior_tsv = os.path.join(
            results_dir, f"{sname}_posterior_usage.tsv"
        )
        coverage_tsv = os.path.join(results_dir, f"{sname}_coverage.tsv")
        debug_json = os.path.join(results_dir, f"{sname}_debug.json")

        cc = CalculateCoverages(
            self.coverage_bw_pos,
            self.coverage_bw_neg,
            bam_file=self.bam_file,
            f_stat_threshold=self.f_stat_threshold,
        )
        usage_df = cc.run(
            subsegments_df,
            output_raw_tsv=coverage_tsv if self.debug else None,
            output_final_tsv=usage_tsv if self.debug else None,
            output_debug_json=debug_json if self.debug else None,
            n_threads=self.n_threads,
            max_pas_count=self.max_pas_count,
            f_stat_threshold=self.f_stat_threshold,
            debug=self.debug,
        )

        # Map atlas_rpm from the in-memory pas_df (rep_cs → atlas_rpm).
        atlas_map = (
            cs.pas_df[["rep_cs", "atlas_rpm"]]
            .drop_duplicates("rep_cs")
            .set_index("rep_cs")["atlas_rpm"]
        )
        usage_df["atlas_rpm"] = (
            usage_df["pas_id"].map(atlas_map).fillna(0.0)
        )

        cpu = CalculatePosteriorUsage(weight=self.posterior_usage_weight)
        posterior_df = cpu.compute(usage_df)

        gene_usage_df = CalculateGeneLevelUsage(posterior_df).compute()

        # Merge gene_weighted_usage into posterior output; drop redundant
        # gene_id and segment_id columns (subsegment_id encodes both).
        posterior_out = posterior_df.merge(
            gene_usage_df, on="gene_id", how="left"
        ).drop(columns=["gene_id", "segment_id"])
        posterior_out.to_csv(posterior_tsv, sep="\t", index=False)
        logger.info("Posterior usage written to %s", posterior_tsv)
        logger.info("PAQR3 pipeline completed.")
