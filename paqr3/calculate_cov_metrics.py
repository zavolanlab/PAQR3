import logging
import pandas as pd  # type: ignore
import pyBigWig  # type: ignore
import numpy as np
import itertools
from collections import defaultdict


def log_message(message):
    logging.info(message)


def evaluate_all_pas_usage_patterns(subsegments, segment_id, debug=False):
    log = logging.info if debug else lambda *args, **kwargs: None

    mean_covs = [s["mean_cov"] for s in subsegments]
    pas_ids = [str(s["pas_id"]) for s in subsegments]

    if all(mu == 0 for mu in mean_covs):
        log(
            f"[{segment_id}] All subsegment mean coverages are 0 → Assign usage = 0 to all PAS"
        )
        return {
            "pas_usage": {pid: 0.0 for pid in pas_ids if pid != "."},
            "rna_monotone": 1,
            "f_stat": 0.0,
            "used_combination": "none",
            "rna_drop_cov": [0.0] * len(subsegments),
            "rna_sum_drop_cov": 0.0,
            "debug_info": "All coverages zero, no PAS used.",
        }

    unique_pas_ids = [pid for pid in pas_ids if pid != "."]
    if not unique_pas_ids:
        log(f"[{segment_id}] No PAS found in this segment, skipping.")
        return None

    n = len(unique_pas_ids)
    best_f_stat = -np.inf
    best_usage = None
    best_monotonic = False
    best_combo = None
    best_drop_cov = []
    best_drop_sum = 0.0
    debug_summary = []

    for binary_pattern in itertools.product([0, 1], repeat=n):
        if sum(binary_pattern) == 0:
            continue

        used_pas = [
            pid for pid, used in zip(unique_pas_ids, binary_pattern) if used
        ]
        log(f"[{segment_id}] Testing PAS combination: {used_pas}")

        groups = []
        current_group = []
        current_pas_idx = 0

        for i, pid in enumerate(pas_ids):
            if pid == ".":
                current_group.append(i)
            else:
                if unique_pas_ids[current_pas_idx] in used_pas:
                    if current_group:
                        groups.append(current_group)
                    current_group = [i]
                else:
                    current_group.append(i)
                current_pas_idx += 1
        if current_group:
            groups.append(current_group)

        group_means = [
            np.mean([mean_covs[i] for i in g]) if g else 0.0 for g in groups
        ]
        is_monotonic = all(
            group_means[i] >= group_means[i + 1]
            for i in range(len(group_means) - 1)
        )

        all_covs = []
        labels = []
        for idx, group in enumerate(groups):
            for i in group:
                all_covs.append(mean_covs[i])
                labels.append(idx)

        group_means_array = [
            np.mean([c for j, c in enumerate(all_covs) if labels[j] == g])
            for g in set(labels)
        ]
        overall_mean = np.mean(all_covs)
        ss_between = sum(
            len([1 for j in labels if j == g]) * (gm - overall_mean) ** 2
            for g, gm in enumerate(group_means_array)
        )
        ss_within = sum(
            (all_covs[i] - group_means_array[labels[i]]) ** 2
            for i in range(len(all_covs))
        )

        df_between = len(set(labels)) - 1
        df_within = len(all_covs) - len(set(labels))
        f_stat = (
            (ss_between / df_between) / (ss_within / df_within)
            if df_within != 0
            else 0.0
        )

        drop_cov = [
            group_means[i] - group_means[i + 1]
            for i in range(len(group_means) - 1)
        ]
        drop_cov.append(0.0)
        drop_sum = sum(drop_cov)

        # 💡 FIX: Usage for single PAS
        if sum(binary_pattern) == 1:
            usage = [1.0]
        else:
            usage = [
                (d / drop_sum if drop_sum > 0 else 0) for d in drop_cov[:-1]
            ]

        usage_dict = {
            pid: u
            for pid, u, used in zip(unique_pas_ids, usage, binary_pattern)
            if used
        }

        debug_summary.append(
            {
                "combo": binary_pattern,
                "f_stat": f_stat,
                "group_means": group_means,
                "monotonic": is_monotonic,
                "groups": groups,
            }
        )

        if f_stat > best_f_stat or (
            f_stat == best_f_stat and is_monotonic and not best_monotonic
        ):
            best_f_stat = f_stat
            best_usage = usage_dict
            best_monotonic = is_monotonic
            best_combo = binary_pattern
            best_drop_cov = drop_cov
            best_drop_sum = drop_sum

    if best_usage is None:
        best_usage = {pid: 0.0 for pid in unique_pas_ids}
    return {
        "pas_usage": best_usage,
        "rna_monotone": int(best_monotonic),
        "f_stat": best_f_stat,
        "used_combination": best_combo,
        "rna_drop_cov": best_drop_cov,
        "rna_sum_drop_cov": best_drop_sum,
        "debug_info": debug_summary,
    }


class CalculateCoverages:
    def __init__(self, coverage_bw_pos, coverage_bw_neg):
        self.coverage_bw_pos = coverage_bw_pos
        self.coverage_bw_neg = coverage_bw_neg

    def calculate_coverage_metrics(self, subsegments_df):
        bw_pos = pyBigWig.open(self.coverage_bw_pos)
        bw_neg = pyBigWig.open(self.coverage_bw_neg)

        all_rows = []

        for _, row in subsegments_df.iterrows():
            gene_id = row["gene_unique_id"]
            segment_number = row["segment_number"]
            subsegment_number = row["subsegment_number"]
            chrom = row["chrom"]
            start = row["start"]
            end = row["end"]
            strand = row["strand"]
            pas_id = row.get("pas_id", ".")

            subsegment_id = f"{gene_id}.{segment_number}.{subsegment_number}"
            segment_id = f"{gene_id}.{segment_number}"

            if strand == "+":
                bw = bw_pos
            elif strand == "-":
                bw = bw_neg
            else:
                raise ValueError(f"Invalid strand: {strand}")

            try:
                coverage = bw.values(chrom, start, end, numpy=True)
                coverage = np.nan_to_num(coverage, nan=0.0)
            except RuntimeError:
                coverage = np.array([])

            mean_cov = np.nanmean(coverage) if len(coverage) > 0 else 0.0
            sum_squared = np.nansum(coverage**2) if len(coverage) > 0 else 0.0

            all_rows.append(
                {
                    "gene_id": gene_id,
                    "segment_id": segment_id,
                    "subsegment_id": subsegment_id,
                    "pas_id": pas_id,
                    "mean_cov": mean_cov,
                    "sum_squared_values": sum_squared,
                }
            )

        bw_pos.close()
        bw_neg.close()

        df = pd.DataFrame(all_rows)

        def check_monotonicity_and_drops(group):
            means = group["mean_cov"].values
            drops = [means[i] - means[i + 1] for i in range(len(means) - 1)]
            drops.append(0.0)

            group = group.copy()
            group["rna_drop_cov"] = drops
            rna_sum_drop_cov = sum(drops)
            group["rna_sum_drop_cov"] = rna_sum_drop_cov
            group["rna_monotone"] = (
                int(all(d > 0 for d in drops[:-1]))
                if rna_sum_drop_cov > 0
                else 0
            )

            return group

        df = df.groupby(["gene_id", "segment_id"], group_keys=False).apply(
            check_monotonicity_and_drops
        )
        return df

    def evaluate_pas_usage_models(self, raw_cov_df, output_tsv_debug=None):
        usage_rows = []
        debug_rows = []

        grouped = raw_cov_df.groupby("segment_id")
        for segment_id, group in grouped:
            subsegments = group.to_dict("records")
            result = evaluate_all_pas_usage_patterns(
                subsegments, segment_id, debug=True
            )

            if result is None:
                continue
            assert (
                result is not None
            ), f"Result became None unexpectedly at segment {segment_id}"
            print(
                f"Evaluating segment {segment_id} → result keys: {result.keys() if result else 'None'}"
            )
            for idx, sub in enumerate(subsegments):
                pid = str(sub["pas_id"])
                if pid == ".":
                    continue
                if pid in result["pas_usage"]:
                    usage_rows.append(
                        {
                            "gene_id": sub["gene_id"],
                            "segment_id": sub["segment_id"],
                            "subsegment_id": sub["subsegment_id"],
                            "pas_id": pid,
                            "mean_cov": sub["mean_cov"],
                            "sum_squared_values": sub["sum_squared_values"],
                            "rna_usage": result["pas_usage"][pid],
                            "rna_monotone": result["rna_monotone"],
                            "rna_drop_cov": (
                                result["rna_drop_cov"][idx]
                                if idx < len(result["rna_drop_cov"])
                                else 0.0
                            ),
                            "rna_sum_drop_cov": result["rna_sum_drop_cov"],
                            "f_stat": result["f_stat"],
                        }
                    )

            if output_tsv_debug:
                debug_rows.append(
                    {
                        "segment_id": segment_id,
                        "used_combination": result["used_combination"],
                        "debug_info": result["debug_info"],
                    }
                )

        usage_df = pd.DataFrame(usage_rows)

        if output_tsv_debug and debug_rows:
            pd.DataFrame(debug_rows).to_json(
                output_tsv_debug, orient="records", lines=True
            )

        return usage_df

    def write_coverage_results(self, results_df, output_tsv):
        results_df.to_csv(output_tsv, sep="\t", index=False)
        log_message(f"Results written to {output_tsv}")

    def run(
        self,
        subsegments_df,
        output_raw_tsv,
        output_final_tsv,
        output_debug_json=None,
    ):
        log_message("Step 1: Calculating initial coverage metrics...")
        raw_cov_df = self.calculate_coverage_metrics(subsegments_df)
        self.write_coverage_results(raw_cov_df, output_raw_tsv)

        log_message("Step 2: Evaluating PAS usage models via F-statistics...")
        refined_usage_df = self.evaluate_pas_usage_models(
            raw_cov_df, output_debug_json
        )
        self.write_coverage_results(refined_usage_df, output_final_tsv)

        log_message("PAS usage evaluation complete.")
