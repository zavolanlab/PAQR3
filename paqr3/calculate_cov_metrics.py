import logging
import pandas as pd  # type: ignore
import pyBigWig  # type: ignore
import numpy as np
import itertools
import json
from scipy.stats import f  # type: ignore


def log_message(message):
    logging.info(message)


def evaluate_all_pas_usage_patterns(subsegments, segment_id, debug=True):
    # Use logging.info as a default logger; also use logging.info for detailed messages.
    log = (
        logging.info if debug is False else logging.info
    )  # We always log info messages
    # (But we'll add debug messages with logging.info below.)

    mean_covs = [s["mean_cov"] for s in subsegments]
    pas_ids = [str(s["pas_id"]) for s in subsegments]

    # Early exit if all coverages are zero.
    if all(mu == 0 for mu in mean_covs):
        msg = f"[{segment_id}] All subsegment mean coverages are 0 → Assign usage = 0 to all PAS"
        log(msg)
        if debug:
            logging.info("DEBUG: " + msg)
        unique_pas_ids = [pid for pid in pas_ids if pid != "."]
        return {
            "pas_usage": {pid: 0.0 for pid in unique_pas_ids},
            "rna_monotone": 0,
            "f_stat": 0.0,
            "p_value": None,
            "used_combination": "none",
            "rna_drop_cov": [0.0] * len(subsegments),
            "rna_sum_drop_cov": 0.0,
            "debug_info": "All coverages zero, no PAS used.",
        }

    unique_pas_ids = [pid for pid in pas_ids if pid != "."]
    if debug:
        logging.info(
            f"DEBUG: Processing segment {segment_id} with PAS: {unique_pas_ids}"
        )
    if not unique_pas_ids:
        msg = f"[{segment_id}] No PAS found in this segment, skipping."
        log(msg)
        if debug:
            logging.info("DEBUG: " + msg)
        return None

    best_f_stat = -np.inf
    best_usage = None
    best_monotonic = False
    best_combo = None
    best_drop_cov = []
    best_drop_sum = 0.0
    best_p_value = None
    debug_summary = []

    m = len(unique_pas_ids)
    # Iterate over all possible combinations of active/inactive for each PAS.
    for binary_pattern in itertools.product([0, 1], repeat=m):
        if sum(binary_pattern) == 0:
            continue

        used_pas = [
            pid
            for pid, flag in zip(unique_pas_ids, binary_pattern)
            if flag == 1
        ]

        # Modified grouping: when an active PAS is encountered, include the current subsegment
        # in the current group before splitting.
        groups = []
        current_group = []
        current_pas_idx = 0
        for i, pid in enumerate(pas_ids):
            if pid == ".":
                current_group.append(i)
            else:
                if unique_pas_ids[current_pas_idx] in used_pas:
                    current_group.append(i)
                    groups.append(current_group)
                    current_group = []
                else:
                    current_group.append(i)
                current_pas_idx += 1
        if current_group:
            groups.append(current_group)

        # Compute group means.
        group_means = [
            np.mean([mean_covs[i] for i in group]) if group else 0.0
            for group in groups
        ]

        # One-way ANOVA: compute F-statistic across groups.
        all_covs = []
        labels = []
        for idx, group in enumerate(groups):
            for i in group:
                all_covs.append(mean_covs[i])
                labels.append(idx)
        if len(set(labels)) > 1:
            group_means_array = [
                np.mean([c for j, c in enumerate(all_covs) if labels[j] == g])
                for g in sorted(set(labels))
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
            if df_within == 0 or ss_within == 0:
                f_stat = 0.0
                p_value = None
            else:
                f_stat = (ss_between / df_between) / (ss_within / df_within)
                p_value = 1 - f.cdf(f_stat, df_between, df_within)
        else:
            f_stat = 0.0
            p_value = None

        # Check monotonicity.
        is_monotonic = all(
            group_means[i] >= group_means[i + 1]
            for i in range(len(group_means) - 1)
        )

        # Always compute drop_list from the groups.
        if len(groups) >= 2:
            drop_list = [
                group_means[i] - group_means[i + 1]
                for i in range(len(groups) - 1)
            ]
        else:
            drop_list = [group_means[0]]
        drop_list.append(0.0)
        drop_sum = sum(drop_list)

        usage_values = []
        drop_idx = 0
        for flag in binary_pattern:
            if flag == 1:
                usage_values.append(
                    (drop_list[drop_idx] / drop_sum) if drop_sum > 0 else 0.0
                )
                drop_idx += 1
            else:
                usage_values.append(0.0)

        usage_dict = {
            pid: u
            for pid, u, flag in zip(
                unique_pas_ids, usage_values, binary_pattern
            )
            if flag == 1
        }

        # Log debug info for this combination.
        if debug:
            logging.info(
                f"DEBUG: Segment {segment_id} combo {binary_pattern}: f_stat={f_stat}, p_value={p_value}, group_means={group_means}, drop_list={drop_list}, usage={usage_dict}"
            )

        debug_summary.append(
            {
                "combo": binary_pattern,
                "f_stat": f_stat,
                "p_value": p_value,
                "group_means": group_means,
                "monotonic": is_monotonic,
                "groups": groups,
                "usage": usage_dict,
                "drop_list": drop_list,
                "drop_sum": drop_sum,
            }
        )

        # Choose the best combination.
        if f_stat > best_f_stat or (
            f_stat == best_f_stat and is_monotonic and not best_monotonic
        ):
            best_f_stat = f_stat
            best_usage = usage_dict
            best_monotonic = is_monotonic
            best_combo = binary_pattern
            best_drop_cov = drop_list
            best_drop_sum = drop_sum
            best_p_value = p_value

    if best_usage is None:
        best_usage = {pid: 0.0 for pid in unique_pas_ids}

    if debug:
        logging.info(
            f"DEBUG: Segment {segment_id} best combo: {best_combo}, best usage: {best_usage}, f_stat: {best_f_stat}, p_value: {best_p_value}"
        )

    return {
        "pas_usage": best_usage,
        "rna_monotone": int(best_monotonic),
        "f_stat": best_f_stat,
        "p_value": best_p_value,
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
                    "segment_number": segment_number,
                }
            )

        bw_pos.close()
        bw_neg.close()

        df = pd.DataFrame(all_rows)

        def compute_drops_and_usage(group):
            group = group.sort_values("subsegment_id")
            mean_covs = group["mean_cov"].values
            pas_ids = group["pas_id"].astype(str).values

            drops = [
                mean_covs[i] - mean_covs[i + 1]
                for i in range(len(mean_covs) - 1)
            ]
            drops.append(0.0)
            sum_drops = sum(drops)

            rna_usage = [
                (
                    (drop / sum_drops)
                    if pid != "." and sum_drops > 0
                    else (None if pid == "." else 0.0)
                )
                for drop, pid in zip(drops, pas_ids)
            ]

            group["rna_drop_cov"] = drops
            group["rna_sum_drop_cov"] = sum_drops
            group["rna_monotone"] = (
                int(all(d > 0 for d in drops[:-1])) if sum_drops > 0 else 1
            )
            group["rna_usage"] = rna_usage

            return group

        df = df.groupby(["gene_id", "segment_id"], group_keys=False).apply(
            compute_drops_and_usage
        )
        df = df.drop(columns=["segment_number"])

        return df

    def evaluate_pas_usage_models(self, raw_cov_df, output_tsv_debug=None):
        usage_rows = []
        debug_rows = []
        # Use sort=False to preserve the original order.
        grouped = raw_cov_df.groupby("segment_id", sort=False)
        for segment_id, group in grouped:
            if output_tsv_debug:
                logging.info(f"DEBUG: Evaluating segment {segment_id}")
            subsegments = group.to_dict("records")
            result = evaluate_all_pas_usage_patterns(
                subsegments, segment_id, debug=True
            )
            if result is None:
                continue
            unique_pas_ids = []
            for sub in subsegments:
                pid = str(sub["pas_id"])
                if pid != "." and pid not in unique_pas_ids:
                    unique_pas_ids.append(pid)
            for pid in unique_pas_ids:
                usage_val = result["pas_usage"].get(pid, 0.0)
                for sub in subsegments:
                    if str(sub["pas_id"]) == pid:
                        usage_rows.append(
                            {
                                "gene_id": sub["gene_id"],
                                "segment_id": sub["segment_id"],
                                "subsegment_id": sub["subsegment_id"],
                                "pas_id": pid,
                                "mean_cov": sub["mean_cov"],
                                "sum_squared_values": sub[
                                    "sum_squared_values"
                                ],
                                "rna_drop_cov": (
                                    result["rna_drop_cov"][
                                        unique_pas_ids.index(pid)
                                    ]
                                    if unique_pas_ids.index(pid)
                                    < len(result["rna_drop_cov"])
                                    else 0.0
                                ),
                                "rna_sum_drop_cov": result["rna_sum_drop_cov"],
                                "rna_monotone": result["rna_monotone"],
                                "rna_usage": usage_val,
                                "f_stat": result["f_stat"],
                                "p_value": result["p_value"],
                            }
                        )
                        break
            # Append debug information for the current segment.
            debug_rows.append(
                {
                    "segment_id": segment_id,
                    "used_combination": result.get("used_combination", []),
                    "debug_info": result.get("debug_info", {}),
                }
            )
        if output_tsv_debug:
            with open(output_tsv_debug, "w") as f:
                json.dump(debug_rows, f, indent=4)
        usage_df = pd.DataFrame(usage_rows)
        # Reorder columns to match _coverage.tsv (with f_stat and p_value appended at the end).
        ordered_columns = [
            "gene_id",
            "segment_id",
            "subsegment_id",
            "pas_id",
            "mean_cov",
            "sum_squared_values",
            "rna_drop_cov",
            "rna_sum_drop_cov",
            "rna_monotone",
            "rna_usage",
            "f_stat",
            "p_value",
        ]
        usage_df = usage_df[ordered_columns]
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
