import logging
import pandas as pd  # type: ignore
import pyBigWig  # type: ignore
import numpy as np
import itertools
import json
from scipy.stats import f  # type: ignore
import time
import concurrent.futures
from functools import partial


def log_message(message):
    logging.info(message)


def json_serial(obj):
    """JSON serializer for objects not serializable by default JSON code"""
    if isinstance(obj, (np.float32, np.float64)):
        return float(obj)
    if isinstance(obj, (np.int32, np.int64)):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Type {type(obj)} not serializable")


def evaluate_all_pas_usage_patterns(
    subsegments, segment_id, debug=True, max_pas_count=10
):
    log = logging.info

    # Retrieve the raw coverage arrays and PAS IDs.
    coverage_arrays = [s["coverage"] for s in subsegments]
    pas_ids = [str(s["pas_id"]) for s in subsegments]
    computed_means = [
        np.mean(cov) if len(cov) > 0 else 0.0 for cov in coverage_arrays
    ]

    # Early exit if all coverage arrays are zero.
    if all(mu == 0 for mu in computed_means):
        msg = f"[{segment_id}] All subsegment coverages are 0 → Assign usage = 0 to all PAS"
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

    # New check: only evaluate segments with max_pas_count or fewer PAS
    if len(unique_pas_ids) > max_pas_count:
        msg = (
            f"[{segment_id}] Too many PAS ({len(unique_pas_ids)}) in this segment; "
            f"only segments with {max_pas_count} or fewer PAS are evaluated. Skipping."
        )
        log(msg)
        if debug:
            logging.info("DEBUG: " + msg)
        return None

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

    best_mono = None  # best monotonic candidate info
    best_mono_f_stat = -np.inf

    m = len(unique_pas_ids)
    # Phase 1: Evaluate only monotonic combinations
    for binary_pattern in itertools.product([0, 1], repeat=m):
        if sum(binary_pattern) == 0:
            continue

        used_pas = [
            pid
            for pid, flag in zip(unique_pas_ids, binary_pattern)
            if flag == 1
        ]

        # Create groups: when an active PAS is encountered, include the current subsegment then split.
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

        # Recalculate group means from the raw coverage values.
        group_means = [
            (
                np.mean(np.concatenate([coverage_arrays[i] for i in group]))
                if group
                else 0.0
            )
            for group in groups
        ]
        # Compute drop_list and drop_sum
        if len(group_means) >= 2:
            drop_list = [
                group_means[i] - group_means[i + 1]
                for i in range(len(group_means) - 1)
            ]
        else:
            drop_list = [group_means[0]]
        drop_list.append(0.0)
        drop_sum = sum(drop_list)
        # Check monotonicity (if drop_sum > 0, then check if non-increasing)
        is_monotonic = False
        if drop_sum > 0:
            is_monotonic = all(
                group_means[i] >= group_means[i + 1]
                for i in range(len(group_means) - 1)
            )
        if not is_monotonic:
            continue  # Skip combinations that are not monotonic

        # For monotonic combinations, compute F-statistics
        all_covs = np.concatenate(
            [
                np.concatenate([coverage_arrays[i] for i in group])
                for group in groups
            ]
        )
        overall_mean = np.mean(all_covs)
        ss_between = sum(
            len(np.concatenate([coverage_arrays[i] for i in group]))
            * (gm - overall_mean) ** 2
            for group, gm in zip(groups, group_means)
        )
        ss_within = sum(
            np.sum(
                (np.concatenate([coverage_arrays[i] for i in group]) - gm) ** 2
            )
            for group, gm in zip(groups, group_means)
        )
        df_between = len(groups) - 1
        df_within = len(all_covs) - len(groups)
        if df_within == 0 or ss_within == 0:
            f_stat = 0.0
            p_value = None
        else:
            f_stat = (ss_between / df_between) / (ss_within / df_within)
            p_value = 1 - f.cdf(f_stat, df_between, df_within)

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

        debug_info = {
            "combo": binary_pattern,
            "f_stat": f_stat,
            "p_value": p_value,
            "group_means": group_means,
            "monotonic": True,
            "groups": groups,
            "usage": usage_dict,
            "drop_list": drop_list,
            "drop_sum": drop_sum,
        }

        if debug:
            logging.info(
                f"DEBUG: Segment {segment_id} monotonic combo {binary_pattern}: f_stat={f_stat}, p_value={p_value}, "
                f"group_means={group_means}, drop_list={drop_list}, usage={usage_dict}"
            )

        if f_stat > best_mono_f_stat:
            best_mono_f_stat = f_stat
            best_mono = {
                "pas_usage": usage_dict,
                "f_stat": f_stat,
                "p_value": p_value,
                "used_combination": binary_pattern,
                "rna_drop_cov": drop_list,
                "rna_sum_drop_cov": drop_sum,
                "debug_info": debug_info,
                "rna_monotone": 1,
            }

    # If no monotonic combination was found, skip this segment.
    if best_mono is None:
        msg = f"[{segment_id}] No monotonic PAS combinations found; skipping segment."
        log(msg)
        if debug:
            logging.info("DEBUG: " + msg)
        return None

    return best_mono


def compute_drops_and_usage(group):
    """Compute rna_drop_cov, rna_sum_drop_cov, rna_monotone, and rna_usage based on precomputed mean_cov values."""
    group = group.sort_values("subsegment_id")
    mean_covs = group["mean_cov"].values
    pas_ids = group["pas_id"].astype(str).values

    drops = [
        mean_covs[i] - mean_covs[i + 1] for i in range(len(mean_covs) - 1)
    ]
    drops.append(0.0)
    sum_drops = sum(drops)
    mono = int(all(d > 0 for d in drops[:-1])) if sum_drops > 0 else 0

    rna_usage = [
        (drop / sum_drops if sum_drops > 0 else 0.0) if pid != "." else None
        for drop, pid in zip(drops, pas_ids)
    ]

    group["rna_drop_cov"] = drops
    group["rna_sum_drop_cov"] = sum_drops
    group["rna_monotone"] = mono
    group["rna_usage"] = rna_usage

    return group


def process_segment(segment_tuple, debug=True, max_pas_count=10):
    segment_id, group = segment_tuple
    subsegments = group.to_dict("records")
    result = evaluate_all_pas_usage_patterns(
        subsegments, segment_id, debug=debug, max_pas_count=max_pas_count
    )
    return segment_id, subsegments, result


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
            seg_num = row["segment_number"]
            subseg_num = row["subsegment_number"]
            chrom = row["chrom"]
            start = row["start"]
            end = row["end"]
            strand = row["strand"]
            pas_id = row.get("pas_id", ".")
            segment_id = f"{gene_id}.{seg_num}"
            subsegment_id = f"{gene_id}.{seg_num}.{subseg_num}"

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
                    "coverage": coverage,  # Keep raw coverage for downstream calculations.
                }
            )

        bw_pos.close()
        bw_neg.close()

        df = pd.DataFrame(all_rows)
        return df

    def evaluate_pas_usage_models(
        self, raw_cov_df, output_tsv_debug=None, n_procs=8, max_pas_count=10
    ):
        usage_rows = []
        debug_rows = []
        grouped = list(raw_cov_df.groupby("segment_id", sort=False))
        if n_procs > 1:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=n_procs
            ) as executor:
                process_func = partial(
                    process_segment, debug=True, max_pas_count=max_pas_count
                )
                results = executor.map(process_func, grouped)
                for segment_id, subsegments, result in results:
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
                                        "rna_sum_drop_cov": result[
                                            "rna_sum_drop_cov"
                                        ],
                                        "rna_monotone": result["rna_monotone"],
                                        "rna_usage": usage_val,
                                        "f_stat": result["f_stat"],
                                        "p_value": result["p_value"],
                                    }
                                )
                                break
                    debug_rows.append(
                        {
                            "segment_id": segment_id,
                            "used_combination": result.get(
                                "used_combination", []
                            ),
                            "debug_info": result.get("debug_info", {}),
                        }
                    )
        else:
            for segment_id, group in grouped:
                logging.info(f"DEBUG: Evaluating segment {segment_id}")
                subsegments = group.to_dict("records")
                result = evaluate_all_pas_usage_patterns(
                    subsegments,
                    segment_id,
                    debug=True,
                    max_pas_count=max_pas_count,
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
                                    "rna_sum_drop_cov": result[
                                        "rna_sum_drop_cov"
                                    ],
                                    "rna_monotone": result["rna_monotone"],
                                    "rna_usage": usage_val,
                                    "f_stat": result["f_stat"],
                                    "p_value": result["p_value"],
                                }
                            )
                            break
                debug_rows.append(
                    {
                        "segment_id": segment_id,
                        "used_combination": result.get("used_combination", []),
                        "debug_info": result.get("debug_info", {}),
                    }
                )
        if output_tsv_debug:
            with open(output_tsv_debug, "w") as f:
                json.dump(debug_rows, f, indent=4, default=json_serial)
        usage_df = pd.DataFrame(usage_rows)
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
        # Drop "coverage" column before writing output.
        if "coverage" in results_df.columns:
            results_df = results_df.drop(columns=["coverage"])
        results_df.to_csv(output_tsv, sep="\t", index=False)
        log_message(f"Results written to {output_tsv}")

    def run(
        self,
        subsegments_df,
        output_raw_tsv,
        output_final_tsv,
        output_debug_json=None,
        n_procs=8,
        max_pas_count=10,
    ):
        start_time = time.time()
        log_message("Step 1: Calculating initial coverage metrics...")
        raw_cov_df = self.calculate_coverage_metrics(subsegments_df)
        # For _coverage.tsv output, compute extra columns using compute_drops_and_usage.
        cov_df = raw_cov_df.groupby(
            ["gene_id", "segment_id"], group_keys=False
        ).apply(compute_drops_and_usage)
        self.write_coverage_results(cov_df.copy(), output_raw_tsv)

        log_message("Step 2: Evaluating PAS usage models via F-statistics...")
        refined_usage_df = self.evaluate_pas_usage_models(
            raw_cov_df,
            output_debug_json,
            n_procs=n_procs,
            max_pas_count=max_pas_count,
        )
        self.write_coverage_results(refined_usage_df, output_final_tsv)

        end_time = time.time()
        log_message("PAS usage evaluation complete.")
        log_message(f"Total runtime: {end_time - start_time:.2f} seconds")
