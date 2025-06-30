import logging
import os

from paqr3.construct_segments import ConstructSegments, log_message
from paqr3.calculate_cov_metrics import CalculateCoverages
from paqr3.calculate_posterior_usage import CalculatePosteriorUsage

# ———————————————— Logging setup ————————————————
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


class PAQR3:
    def __init__(
        self,
        annotation_file,
        pas_atlas_file,
        coverage_bw_pos,
        coverage_bw_neg,
        output_dir,
        downstream_exon_extension,
        merge_distance,
        max_pas_count,
        bam_file=None,
        f_stat_threshold=100,
        posterior_usage_weight=0.1,
        debug=False,
    ):
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
        self.debug = debug

    def run(self):
        # 1) Sample name / output folder
        sample_name = os.path.basename(self.coverage_bw_pos).split(".")[0]
        if sample_name != os.path.basename(self.coverage_bw_neg).split(".")[0]:
            raise ValueError("Coverage files must have the same sample name.")

        results_dir = os.path.join(self.output_dir, f"{sample_name}_results")
        os.makedirs(results_dir, exist_ok=True)

        # 2) Paths for intermediate outputs
        genes_bed = os.path.join(results_dir, f"{sample_name}_genes.bed")
        segments_bed = os.path.join(results_dir, f"{sample_name}_segments.bed")
        subsegments_bed = os.path.join(
            results_dir, f"{sample_name}_subsegments.bed"
        )
        pas_bed = os.path.join(results_dir, f"{sample_name}_PAS.bed")
        coverage_tsv = os.path.join(results_dir, f"{sample_name}_coverage.tsv")
        usage_tsv = os.path.join(results_dir, f"{sample_name}_usage.tsv")
        debug_json = os.path.join(results_dir, f"{sample_name}_debug.json")

        # 3) Build segments + subsegments
        cs = ConstructSegments(
            self.annotation_file,
            self.pas_atlas_file,
            self.downstream_exon_extension,
        )
        cs.run(
            genes_bed,
            segments_bed,
            subsegments_bed,
            pas_bed,
            self.merge_distance,
        )
        subsegments_df = cs.create_subsegments_dataframe()

        # 4) Compute coverage + F-stat usage
        cc = CalculateCoverages(
            self.coverage_bw_pos,
            self.coverage_bw_neg,
            bam_file=self.bam_file,
            f_stat_threshold=self.f_stat_threshold,
        )
        refined_usage_df = cc.run(
            subsegments_df,
            output_raw_tsv=coverage_tsv,
            output_final_tsv=usage_tsv,
            output_debug_json=debug_json if self.debug else None,
            max_pas_count=self.max_pas_count,
            f_stat_threshold=self.f_stat_threshold,
            debug=self.debug,
        )

        # 5) Grab the merged PAS DataFrame (with atlas RPM) from ConstructSegments
        atlas_df = cs.pas_df[["pas_id", "atlas_rpm"]].copy()

        # 6) Compute Bayesian‐style posterior usage
        cpu = CalculatePosteriorUsage(
            atlas_df, weight=self.posterior_usage_weight
        )
        posterior_df = cpu.compute(refined_usage_df)

        # 7) Write out posterior usage for comparison
        posterior_tsv = os.path.join(
            results_dir, f"{sample_name}_posterior_usage.tsv"
        )
        posterior_df.to_csv(posterior_tsv, sep="\t", index=False)
        log_message(f"Posterior usage written to {posterior_tsv}")

        log_message("PAQR3 pipeline completed.")
