import logging
import pandas as pd  # type: ignore
import pyBigWig  # type: ignore
import numpy as np
import itertools
import json
import time
import concurrent.futures
import pysam  # type: ignore
from functools import partial
from scipy.stats import f  # type: ignore
from pybedtools import BedTool  # type: ignore


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
    subsegments,
    segment_id,
    debug=True,
    max_pas_count=10,
    f_stat_threshold=100,
):
    """
    Select all monotonic PAS combinations with F-stat >= threshold,
    then build a single 'union' pattern of all PAS used in any of them,
    and compute drops, usage, and monotonicity from that union pattern.
    For segments with exactly one PAS, assigns usage=1.0 if coverage>0.
    """
    log = logging.info

    coverage_arrays = [s["coverage"] for s in subsegments]
    pas_ids = [str(s["pas_id"]) for s in subsegments]
    computed_means = [
        np.mean(cov) if len(cov) > 0 else 0.0 for cov in coverage_arrays
    ]
    unique_pas_ids = [pid for pid in pas_ids if pid != "."]

    # 1) Early exit if all coverage zero
    if all(mu == 0 for mu in computed_means):
        msg = f"[{segment_id}] All subsegment coverages are 0 → Assign usage = 0 to all PAS"
        log(msg)
        if debug:
            logging.info("DEBUG: " + msg)
        zero_drops = [0.0] * len(unique_pas_ids)
        return {
            "pas_usage": {pid: 0.0 for pid in unique_pas_ids},
            "f_stat": None,
            "p_value": None,
            "rna_drop_cov": zero_drops,
            "rna_sum_drop_cov": 0.0,
            "rna_monotone": 0,
            "used_combos": [],
            "debug_info": "All coverages zero, no PAS used.",
            "unique_pas_ids": unique_pas_ids,
        }

    # 2) Single‐PAS segment: force usage = 1 if coverage > 0
    if len(unique_pas_ids) == 1:
        pid = unique_pas_ids[0]
        # find first coverage for that PAS
        idx = next(i for i, x in enumerate(pas_ids) if x == pid)
        mu = computed_means[idx]
        if mu > 0:
            usage = 1.0
            drop = mu
            sum_drop = mu
            mono = 1
        else:
            usage = 0.0
            drop = 0.0
            sum_drop = 0.0
            mono = 0
        if debug:
            logging.info(
                f"DEBUG: [{segment_id}] Single PAS '{pid}', mean_cov={mu:.3f} → usage={usage}"
            )
        return {
            "pas_usage": {pid: usage},
            "f_stat": None,
            "p_value": None,
            "rna_drop_cov": [drop],
            "rna_sum_drop_cov": sum_drop,
            "rna_monotone": mono,
            "used_combos": [],
            "debug_info": "Single-PAS segment; forced usage.",
            "unique_pas_ids": unique_pas_ids,
        }

    # 3) Too many PAS? skip
    if len(unique_pas_ids) > max_pas_count:
        msg = (
            f"[{segment_id}] Too many PAS ({len(unique_pas_ids)}) in this segment; "
            f"only segments with {max_pas_count} or fewer PAS are evaluated. Skipping."
        )
        log(msg)
        if debug:
            logging.info("DEBUG: " + msg)
        return None

    # 4) No PAS at all? skip
    if not unique_pas_ids:
        msg = f"[{segment_id}] No PAS found in this segment, skipping."
        log(msg)
        if debug:
            logging.info("DEBUG: " + msg)
        return None

    # 5) Evaluate every monotonic combo, collect those above threshold
    combos_info = []
    m = len(unique_pas_ids)
    for pattern in itertools.product([0, 1], repeat=m):
        if sum(pattern) == 0:
            continue

        groups, current, idx = [], [], 0
        for cov_idx, pid in enumerate(pas_ids):
            if pid == ".":
                current.append(cov_idx)
            else:
                if pattern[idx] == 1:
                    current.append(cov_idx)
                    groups.append(current)
                    current = []
                else:
                    current.append(cov_idx)
                idx += 1
        if current:
            groups.append(current)

        group_means = [
            (
                np.mean(np.concatenate([coverage_arrays[i] for i in g]))
                if g
                else 0.0
            )
            for g in groups
        ]

        if len(group_means) >= 2:
            drops = [
                group_means[i] - group_means[i + 1]
                for i in range(len(group_means) - 1)
            ]
        else:
            drops = [group_means[0]]
        drops.append(0.0)
        total_drops = sum(drops)
        if not (
            total_drops > 0
            and all(
                group_means[i] >= group_means[i + 1]
                for i in range(len(group_means) - 1)
            )
        ):
            continue

        all_covs = np.concatenate(
            [np.concatenate([coverage_arrays[i] for i in g]) for g in groups]
        )
        overall_mean = np.mean(all_covs)
        ss_between = sum(
            len(np.concatenate([coverage_arrays[i] for i in g]))
            * (gm - overall_mean) ** 2
            for g, gm in zip(groups, group_means)
        )
        ss_within = sum(
            np.sum((np.concatenate([coverage_arrays[i] for i in g]) - gm) ** 2)
            for g, gm in zip(groups, group_means)
        )
        df_b, df_w = len(groups) - 1, len(all_covs) - len(groups)
        if df_w == 0 or ss_within == 0:
            f_stat, p_val = 0.0, None
        else:
            f_stat = (ss_between / df_b) / (ss_within / df_w)
            p_val = 1 - f.cdf(f_stat, df_b, df_w)

        if debug:
            logging.info(
                f"DEBUG: Segment {segment_id} pattern {pattern}: f_stat={f_stat:.3f}"
            )

        if f_stat >= f_stat_threshold:
            combos_info.append((pattern, f_stat, p_val))

    # 6) If none passed, skip
    if not combos_info:
        msg = f"[{segment_id}] No PAS combinations met the F-stat threshold; skipping segment."
        log(msg)
        if debug:
            logging.info("DEBUG: " + msg)
        return None

    # 7) Build the union pattern and compute final drops & usage
    union_pattern = [
        int(any(pat[i] for pat, *_ in combos_info)) for i in range(m)
    ]
    groups, current, idx = [], [], 0
    for cov_idx, pid in enumerate(pas_ids):
        if pid == ".":
            current.append(cov_idx)
        else:
            if union_pattern[idx] == 1:
                current.append(cov_idx)
                groups.append(current)
                current = []
            else:
                current.append(cov_idx)
            idx += 1
    if current:
        groups.append(current)

    group_means = [
        np.mean(np.concatenate([coverage_arrays[i] for i in g])) if g else 0.0
        for g in groups
    ]
    if debug:
        logging.info(
            f"DEBUG: Segment {segment_id} union group means: {group_means}"
        )

    if len(group_means) >= 2:
        union_drops = [
            group_means[i] - group_means[i + 1]
            for i in range(len(group_means) - 1)
        ]
    else:
        union_drops = [group_means[0]]
    union_drops.append(0.0)
    sum_union_drops = sum(union_drops)

    union_mono = int(
        sum_union_drops > 0 and all(d >= 0 for d in union_drops[:-1])
    )

    drop_per_pas, drop_idx = [], 0
    for flag in union_pattern:
        if flag:
            drop_per_pas.append(union_drops[drop_idx])
            drop_idx += 1
        else:
            drop_per_pas.append(0.0)

    usage_per_pas = [
        (d / sum_union_drops if sum_union_drops > 0 else 0.0)
        for d in drop_per_pas
    ]
    best_pattern, best_f, best_p = max(combos_info, key=lambda x: x[1])

    if debug:
        logging.info(
            f"DEBUG: Segment {segment_id} union_pattern: {union_pattern}"
        )
        logging.info(f"DEBUG: Segment {segment_id} union_drops: {union_drops}")
        logging.info(
            f"DEBUG: Segment {segment_id} union_monotonic: {union_mono}"
        )
        logging.info(
            f"DEBUG: Segment {segment_id} usage_per_pas: {usage_per_pas}"
        )

    return {
        "pas_usage": {
            pid: usage for pid, usage in zip(unique_pas_ids, usage_per_pas)
        },
        "f_stat": best_f,
        "p_value": best_p,
        "rna_drop_cov": drop_per_pas,
        "rna_sum_drop_cov": sum_union_drops,
        "rna_monotone": union_mono,
        "used_combos": [pat for pat, *_ in combos_info],
        "debug_info": {
            "threshold": f_stat_threshold,
            "num_kept_combos": len(combos_info),
        },
        "unique_pas_ids": unique_pas_ids,
    }


def compute_drops_and_usage(group):
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


def process_segment(
    segment_tuple,
    debug=True,
    max_pas_count=10,
    f_stat_threshold=100,
):
    segment_id, group = segment_tuple
    subsegments = group.to_dict("records")
    result = evaluate_all_pas_usage_patterns(
        subsegments,
        segment_id,
        debug=debug,
        max_pas_count=max_pas_count,
        f_stat_threshold=f_stat_threshold,
    )
    return segment_id, subsegments, result


class CalculateCoverages:
    def __init__(
        self,
        coverage_bw_pos,
        coverage_bw_neg,
        bam_file=None,
        f_stat_threshold=100,
    ):
        self.coverage_bw_pos = coverage_bw_pos
        self.coverage_bw_neg = coverage_bw_neg
        self.bam_file = bam_file
        self.f_stat_threshold = f_stat_threshold

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

            bw = bw_pos if strand == "+" else bw_neg

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
                    "coverage": coverage,
                }
            )

        bw_pos.close()
        bw_neg.close()
        return pd.DataFrame(all_rows)

    def evaluate_pas_usage_models(
        self,
        raw_cov_df,
        output_tsv_debug=None,
        n_procs=8,
        max_pas_count=10,
        f_stat_threshold=100,
    ):
        usage_rows = []
        debug_rows = []
        grouped = list(raw_cov_df.groupby("segment_id", sort=False))

        if n_procs > 1:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=n_procs
            ) as executor:
                func = partial(
                    process_segment,
                    debug=True,
                    max_pas_count=max_pas_count,
                    f_stat_threshold=f_stat_threshold,
                )
                for segment_id, subsegments, result in executor.map(
                    func, grouped
                ):
                    if not result:
                        continue
                    uids = result["unique_pas_ids"]
                    for pid, usage in result["pas_usage"].items():
                        drop_idx = uids.index(pid)
                        for s in subsegments:
                            if str(s["pas_id"]) == pid:
                                usage_rows.append(
                                    {
                                        "gene_id": s["gene_id"],
                                        "segment_id": s["segment_id"],
                                        "subsegment_id": s["subsegment_id"],
                                        "pas_id": pid,
                                        "mean_cov": s["mean_cov"],
                                        "sum_squared_values": s[
                                            "sum_squared_values"
                                        ],
                                        "rna_drop_cov": result["rna_drop_cov"][
                                            drop_idx
                                        ],
                                        "rna_sum_drop_cov": result[
                                            "rna_sum_drop_cov"
                                        ],
                                        "rna_monotone": result["rna_monotone"],
                                        "rna_usage": usage,
                                        "f_stat": result["f_stat"],
                                        "p_value": result["p_value"],
                                    }
                                )
                                break
                    debug_rows.append(
                        {
                            "segment_id": segment_id,
                            "used_combos": result["used_combos"],
                            "debug_info": result["debug_info"],
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
                    f_stat_threshold=f_stat_threshold,
                )
                if not result:
                    continue
                uids = result["unique_pas_ids"]
                for pid, usage in result["pas_usage"].items():
                    drop_idx = uids.index(pid)
                    for s in subsegments:
                        if str(s["pas_id"]) == pid:
                            usage_rows.append(
                                {
                                    "gene_id": s["gene_id"],
                                    "segment_id": s["segment_id"],
                                    "subsegment_id": s["subsegment_id"],
                                    "pas_id": pid,
                                    "mean_cov": s["mean_cov"],
                                    "sum_squared_values": s[
                                        "sum_squared_values"
                                    ],
                                    "rna_drop_cov": result["rna_drop_cov"][
                                        drop_idx
                                    ],
                                    "rna_sum_drop_cov": result[
                                        "rna_sum_drop_cov"
                                    ],
                                    "rna_monotone": result["rna_monotone"],
                                    "rna_usage": usage,
                                    "f_stat": result["f_stat"],
                                    "p_value": result["p_value"],
                                }
                            )
                            break
                debug_rows.append(
                    {
                        "segment_id": segment_id,
                        "used_combos": result["used_combos"],
                        "debug_info": result["debug_info"],
                    }
                )

        if output_tsv_debug:
            with open(output_tsv_debug, "w") as f:
                json.dump(debug_rows, f, indent=4, default=json_serial)

        usage_df = pd.DataFrame(usage_rows)
        cols = [
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
        return usage_df[cols]

    def write_coverage_results(self, results_df, output_tsv):
        if "coverage" in results_df.columns:
            results_df = results_df.drop(columns=["coverage"])
        results_df.to_csv(output_tsv, sep="\t", index=False)
        log_message(f"Results written to {output_tsv}")

    def calculate_segment_expression_stats(self, usage_df, subsegments_df):
        """
        Calculate per-segment:
          - rpm = read_count / segment_length
          - rpm_rank (dense rank, highest=1)
          - median_cov across the segment (from BigWig)
          - median_rank (dense rank, highest=1)
        """

        # 1) Build segment metadata
        seg_meta = subsegments_df.groupby(
            ["gene_unique_id", "segment_number"], as_index=False
        ).agg(
            chrom=("chrom", "first"),
            strand=("strand", "first"),
            start=("start", "min"),
            end=("end", "max"),
        )
        seg_meta["segment_id"] = (
            seg_meta["gene_unique_id"].astype(str)
            + "."
            + seg_meta["segment_number"].astype(str)
        )

        # 2) Count reads per segment via pysam (streaming, low memory)
        if not self.bam_file:
            raise ValueError("BAM file is required to compute RPM.")
        bam = pysam.AlignmentFile(self.bam_file, "rb")
        read_counts = []
        for _, row in seg_meta.iterrows():
            count = bam.count(
                contig=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
            )
            read_counts.append(count)
        bam.close()

        seg_meta["read_count"] = read_counts

        # 3) Compute rpm and median coverage from BigWig
        bw_pos = pyBigWig.open(self.coverage_bw_pos)
        bw_neg = pyBigWig.open(self.coverage_bw_neg)

        def _get_median_cov(r):
            bw = bw_pos if r["strand"] == "+" else bw_neg
            vals = bw.values(
                r["chrom"], int(r["start"]), int(r["end"]), numpy=True
            )
            return float(np.nanmedian(np.nan_to_num(vals, nan=0.0)))

        seg_meta["segment_length"] = seg_meta["end"] - seg_meta["start"]
        seg_meta["rpm"] = seg_meta["read_count"] / seg_meta["segment_length"]
        seg_meta["median_cov"] = seg_meta.apply(_get_median_cov, axis=1)

        bw_pos.close()
        bw_neg.close()

        # 4) Dense ranking (highest value → rank 1)
        seg_meta["rpm_rank"] = (
            seg_meta["rpm"].rank(method="dense", ascending=False).astype(int)
        )
        seg_meta["median_rank"] = (
            seg_meta["median_cov"]
            .rank(method="dense", ascending=False)
            .astype(int)
        )

        # 5) Merge back onto the PAS‐level usage_df
        merged = usage_df.merge(
            seg_meta[
                ["segment_id", "rpm", "rpm_rank", "median_cov", "median_rank"]
            ],
            on="segment_id",
            how="left",
        )

        # 6) Enforce column order
        cols = [
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
            "rpm",
            "rpm_rank",
            "median_cov",
            "median_rank",
        ]
        return merged[cols]

    def run(
        self,
        subsegments_df,
        output_raw_tsv,
        output_final_tsv,
        output_debug_json=None,
        n_procs=8,
        max_pas_count=10,
        f_stat_threshold=None,
    ):
        start_time = time.time()

        # Step 1: initial coverage metrics
        log_message("Step 1: Calculating initial coverage metrics...")
        raw_cov_df = self.calculate_coverage_metrics(subsegments_df)
        cov_df = raw_cov_df.groupby(
            ["gene_id", "segment_id"], group_keys=False
        ).apply(compute_drops_and_usage)
        self.write_coverage_results(cov_df.copy(), output_raw_tsv)

        # Step 2: PAS usage via F-statistics
        log_message("Step 2: Evaluating PAS usage models via F-statistics...")
        # prefer the passed-in threshold, else use the one from __init__
        thr = (
            f_stat_threshold
            if f_stat_threshold is not None
            else self.f_stat_threshold
        )
        refined_usage_df = self.evaluate_pas_usage_models(
            raw_cov_df,
            output_debug_json,
            n_procs=n_procs,
            max_pas_count=max_pas_count,
            f_stat_threshold=thr,
        )

        # free memory from raw coverage arrays
        del raw_cov_df
        import gc

        gc.collect()

        # Step 3: per-segment expression stats (if BAM provided)
        if self.bam_file:
            refined_usage_df = self.calculate_segment_expression_stats(
                refined_usage_df, subsegments_df
            )

        # Final write
        self.write_coverage_results(refined_usage_df, output_final_tsv)

        end_time = time.time()
        log_message("PAS usage evaluation complete.")
        log_message(f"Total runtime: {end_time - start_time:.2f} seconds")

        return refined_usage_df
