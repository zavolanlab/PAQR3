import logging
import os
from paqr3.construct_segments import ConstructSegments
from paqr3.calculate_cov_metrics import CalculateCoverages

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


def log_message(message):
    logging.info(message)


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
    ):
        self.annotation_file = annotation_file
        self.pas_atlas_file = pas_atlas_file
        self.coverage_bw_pos = coverage_bw_pos
        self.coverage_bw_neg = coverage_bw_neg
        self.output_dir = output_dir
        self.downstream_exon_extension = downstream_exon_extension
        self.merge_distance = merge_distance

    def run(self):
        sample_name_pos = os.path.basename(self.coverage_bw_pos).split(".")[0]
        sample_name_neg = os.path.basename(self.coverage_bw_neg).split(".")[0]

        if sample_name_pos != sample_name_neg:
            raise ValueError("Coverage files must have the same sample name.")

        sample_name = sample_name_pos
        output_dir = os.path.join(self.output_dir, f"{sample_name}_results")
        os.makedirs(output_dir, exist_ok=True)

        # Output paths
        output_genes_bed = os.path.join(output_dir, f"{sample_name}_genes.bed")
        output_segments_bed = os.path.join(
            output_dir, f"{sample_name}_segments.bed"
        )
        output_subsegments_bed = os.path.join(
            output_dir, f"{sample_name}_subsegments.bed"
        )
        output_pas_bed = os.path.join(output_dir, f"{sample_name}_PAS.bed")
        output_coverage_tsv = os.path.join(
            output_dir, f"{sample_name}_coverage.tsv"
        )
        output_usage_tsv = os.path.join(output_dir, f"{sample_name}_usage.tsv")
        output_debug_json = os.path.join(
            output_dir, f"{sample_name}_debug.json"
        )

        # Construct segments
        construct_segments = ConstructSegments(
            self.annotation_file,
            self.pas_atlas_file,
            self.downstream_exon_extension,
        )
        construct_segments.run(
            output_genes_bed,
            output_segments_bed,
            output_subsegments_bed,
            output_pas_bed,
            self.merge_distance,
        )

        # Calculate coverages and usage
        calculate_coverages = CalculateCoverages(
            self.coverage_bw_pos, self.coverage_bw_neg
        )
        subsegments_df = construct_segments.create_subsegments_dataframe()

        calculate_coverages.run(
            subsegments_df,
            output_raw_tsv=output_coverage_tsv,
            output_final_tsv=output_usage_tsv,
            output_debug_json=output_debug_json,
        )

        log_message("PAQR3 pipeline completed.")
