import os
from datetime import datetime
from paqr3.construct_segments import ConstructSegments
from paqr3.calculate_coverages import CalculateCoverages


def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


class PAQR3:
    def __init__(
        self,
        annotation_file,
        pas_atlas_file,
        coverage_bw_pos,
        coverage_bw_neg,
        output_dir,
        downstream_exon_extension,
    ):
        self.annotation_file = annotation_file
        self.pas_atlas_file = pas_atlas_file
        self.coverage_bw_pos = coverage_bw_pos
        self.coverage_bw_neg = coverage_bw_neg
        self.output_dir = output_dir
        self.downstream_exon_extension = downstream_exon_extension

    def run(self):
        # Extract sample name from coverage files
        sample_name_pos = os.path.basename(self.coverage_bw_pos).split("_")[0]
        sample_name_neg = os.path.basename(self.coverage_bw_neg).split("_")[0]

        if sample_name_pos != sample_name_neg:
            raise ValueError("Coverage files must have the same sample name.")

        sample_name = sample_name_pos

        output_dir = os.path.join(self.output_dir, f"{sample_name}_results")
        os.makedirs(output_dir, exist_ok=True)

        output_gtf = os.path.join(output_dir, f"{sample_name}.gtf")
        output_genes_tsv = os.path.join(output_dir, f"{sample_name}_genes.tsv")
        output_segments_tsv = os.path.join(
            output_dir, f"{sample_name}_segments.tsv"
        )
        output_subsegments_tsv = os.path.join(
            output_dir, f"{sample_name}_subsegments.tsv"
        )
        output_coverage_tsv = os.path.join(
            output_dir, f"{sample_name}_coverage.tsv"
        )

        # Construct segments
        construct_segments = ConstructSegments(
            self.annotation_file,
            self.pas_atlas_file,
            self.downstream_exon_extension,
        )
        construct_segments.run(
            output_gtf,
            output_genes_tsv,
            output_segments_tsv,
            output_subsegments_tsv,
        )

        # Calculate coverages
        calculate_coverages = CalculateCoverages(
            self.coverage_bw_pos, self.coverage_bw_neg
        )
        subsegments_df = construct_segments.create_subsegments_dataframe()

        calculate_coverages.run(subsegments_df, output_coverage_tsv)

        log_message(
            "Genomic segment construction, PAS identification, "
            "subsegment construction and coverage calculation pipeline completed."
        )
