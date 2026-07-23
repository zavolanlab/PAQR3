"""Coverage metrics and PAS usage evaluation for PAQR3."""

import concurrent.futures
import json
import logging
import time
import warnings
from functools import partial

import numpy as np
import pandas as pd  # type: ignore
import pyBigWig  # type: ignore
from pybedtools import BedTool  # type: ignore
from scipy.stats import f  # type: ignore

warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")


logger = logging.getLogger(__name__)


def json_serial(obj):
    """JSON serialiser for NumPy types; pass as json.dump default.

    Args:
        obj: Object to serialise.

    Returns:
        float, int, or list.

    Raises:
        TypeError: If obj is not a recognised NumPy type.
    """
    if isinstance(obj, (np.float32, np.float64)):
        return float(obj)
    if isinstance(obj, (np.int32, np.int64)):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Type {type(obj)} not serializable")


def _assign_cluster_ids(
    pas_centers: list,
    cluster_distance: float,
) -> list[int]:
    """Assign a cluster ID to each PAS based on centre-to-centre distance.

    Consecutive PAS whose centres are within cluster_distance bp of
    each other receive the same cluster ID.  IDs are 0-based integers
    that increase whenever a gap exceeds the threshold.

    Args:
        pas_centers: Centre coordinates (float or None) for each PAS
            in gene order.  None entries force a cluster break.
        cluster_distance: Maximum centre-to-centre distance to merge.

    Returns:
        List of integer cluster IDs, same length as pas_centers.
    """
    if not pas_centers:
        return []
    ids = [0]
    current = 0
    for k in range(1, len(pas_centers)):
        p1, p2 = pas_centers[k - 1], pas_centers[k]
        if p1 is None or p2 is None or abs(p2 - p1) > cluster_distance:
            current += 1
        ids.append(current)
    return ids


def evaluate_all_pas_usage_patterns(
    subsegments,
    segment_id,
    debug=False,
    max_pas_count=10,
    f_stat_threshold=100,
    cluster_distance=0,
):
    """Find PAS usage patterns satisfying monotonicity and F-stat criteria.

    Uses a pruned DFS over binary assignment patterns on PAS clusters,
    keeping those where coverage is non-increasing across groups and
    whose F-statistic meets the threshold.  Nearby PAS within
    cluster_distance are merged into a single unit for the ANOVA; their
    individual drops are then distributed proportionally by atlas_rpm.

    Args:
        subsegments: List of sub-segment dicts with coverage (numpy
            array), pas_id ("." for non-PAS), and optionally
            pas_start, pas_end, atlas_rpm keys.
        segment_id: Used in debug log messages.
        debug: Emit verbose per-pattern log messages.
        max_pas_count: Skip segments with more clusters than this.
        f_stat_threshold: Minimum F-statistic to retain a pattern.
        cluster_distance: Maximum centre-to-centre distance (bp) to
            merge consecutive PAS into one cluster.  0 = no clustering.

    Returns:
        Dict with keys pas_usage, f_stat, p_value, rna_drop_cov,
        rna_sum_drop_cov, rna_monotone, used_combos, debug_info,
        unique_pas_ids. Returns None if no patterns pass or the
        segment is skipped.
    """
    coverage_arrays = [s["coverage"] for s in subsegments]
    pas_ids = [str(s["pas_id"]) for s in subsegments]
    computed_means = [
        np.mean(cov) if len(cov) > 0 else 0.0 for cov in coverage_arrays
    ]
    unique_pas_ids = [pid for pid in pas_ids if pid != "."]

    # 1) Early exit if all coverage zero
    if all(mu == 0 for mu in computed_means):
        if debug:
            logger.debug(
                "DEBUG: [%s] All subsegment coverages are 0 "
                "→ Assign usage = 0 to all PAS",
                segment_id,
            )
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

    # 2) Single-PAS segment: force usage = 1 if coverage > 0
    if len(unique_pas_ids) == 1:
        pid = unique_pas_ids[0]
        idx = next(i for i, x in enumerate(pas_ids) if x == pid)
        mu = computed_means[idx]
        if mu > 0:
            usage, drop, sum_drop, mono = 1.0, mu, mu, 1
        else:
            usage, drop, sum_drop, mono = 0.0, 0.0, 0.0, 0
        if debug:
            logger.debug(
                "DEBUG: [%s] Single PAS '%s', mean_cov=%.3f → usage=%.3f",
                segment_id, pid, mu, usage,
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

    # 3) No PAS at all? skip
    if not unique_pas_ids:
        if debug:
            logger.debug(
                "DEBUG: [%s] No PAS found in this segment, skipping.",
                segment_id,
            )
        return None

    # Extract atlas_rpm and centre positions for each PAS (in PAS order).
    atlas_rpms: list[float] = []
    pas_centers: list = []
    for s in subsegments:
        if str(s["pas_id"]) == ".":
            continue
        ar = s.get("atlas_rpm", 0.0)
        try:
            ar = float(ar)
            if np.isnan(ar):
                ar = 0.0
        except (TypeError, ValueError):
            ar = 0.0
        atlas_rpms.append(ar)

        ps = s.get("pas_start")
        pe = s.get("pas_end")
        try:
            center = (float(ps) + float(pe)) / 2.0
            if np.isnan(center):
                center = None
        except (TypeError, ValueError):
            center = None
        pas_centers.append(center)

    # Build per-PAS slices: slices[k] = coverage indices from after
    # PAS k-1 up to and including PAS k.  tail = indices after last PAS.
    m_orig = len(unique_pas_ids)
    slices: list[list[int]] = []
    current_slice: list[int] = []
    for cov_idx, pid in enumerate(pas_ids):
        current_slice.append(cov_idx)
        if pid != ".":
            slices.append(current_slice)
            current_slice = []
    tail = current_slice

    # Cluster consecutive PAS whose centres are ≤ cluster_distance apart,
    # then merge their slices so the DFS operates on clusters.
    if cluster_distance > 0 and m_orig > 1:
        cluster_id_per_pas = _assign_cluster_ids(pas_centers, cluster_distance)
        n_clusters = cluster_id_per_pas[-1] + 1
        virtual_slices: list[list[int]] = [[] for _ in range(n_clusters)]
        for i, cid in enumerate(cluster_id_per_pas):
            virtual_slices[cid].extend(slices[i])
        slices = virtual_slices
        if debug and n_clusters < m_orig:
            logger.debug(
                "DEBUG: [%s] Clustered %d PAS → %d clusters "
                "(cluster_distance=%d bp)",
                segment_id, m_orig, n_clusters, cluster_distance,
            )
    else:
        cluster_id_per_pas = list(range(m_orig))
        n_clusters = m_orig

    m = n_clusters  # DFS dimensionality = number of clusters

    # 4) Too many clusters? skip
    if m > max_pas_count:
        if debug:
            logger.debug(
                "DEBUG: [%s] Too many clusters (%d); "
                "only segments with %d or fewer are evaluated. Skipping.",
                segment_id, m, max_pas_count,
            )
        return None

    # 5) Evaluate every monotone cluster pattern via pruned DFS.
    #
    # slices[k] now holds concatenated coverage indices for all PAS in
    # cluster k.  The DFS explores cut/no-cut at each cluster boundary.
    _concat_cache: dict = {}

    def _group_array(indices):
        if len(indices) == 1:
            return coverage_arrays[indices[0]]
        key = frozenset(indices)
        if key not in _concat_cache:
            _concat_cache[key] = np.concatenate(
                [coverage_arrays[i] for i in indices]
            )
        return _concat_cache[key]

    def _group_mean(indices):
        arr = _group_array(indices)
        return float(np.mean(arr)) if len(arr) > 0 else 0.0

    combos_info = []

    def _dfs(k, g_start, prev_mean, groups, bits):
        """DFS over binary cut/no-cut decisions at each cluster position.

        Prunes any branch where the newly closed group's mean exceeds
        the previous group's mean, cutting off non-monotone sub-trees
        before their F-statistics are computed.

        Args:
            k: Index of the cluster slice currently being considered.
            g_start: Slice index at which the current open group began.
            prev_mean: Mean of the last closed group; float("inf") at
                the root so any first group is accepted.
            groups: Closed groups accumulated so far (mutated in place;
                restored on backtrack).
            bits: Pattern bits accumulated so far (mutated in place;
                restored on backtrack).
        """
        if k == len(slices):
            if not groups:
                return  # all-zeros pattern — skip
            last_group = [i for s in slices[g_start:] for i in s] + tail
            if last_group:
                if _group_mean(last_group) > prev_mean:
                    return  # tail violates monotonicity
            complete = groups + ([last_group] if last_group else [])
            if len(complete) < 2:
                return  # single group — df_between = 0
            all_means = [_group_mean(g) for g in complete]
            group_arrays = [_group_array(g) for g in complete]
            all_covs = np.concatenate(group_arrays)
            overall_mean = float(np.mean(all_covs))
            ss_between = sum(
                len(arr) * (gm - overall_mean) ** 2
                for arr, gm in zip(group_arrays, all_means)
            )
            ss_within = sum(
                float(np.sum((arr - gm) ** 2))
                for arr, gm in zip(group_arrays, all_means)
            )
            df_b = len(complete) - 1
            df_w = len(all_covs) - len(complete)
            if (
                df_b == 0
                or df_w == 0
                or ss_within == 0
                or np.isnan(ss_within)
                or np.isnan(ss_between)
            ):
                f_stat_val = 0.0
                p_val = None
            else:
                f_stat_val = (ss_between / df_b) / (ss_within / df_w)
                p_val = float(1 - f.cdf(f_stat_val, df_b, df_w))
            if debug:
                logger.debug(
                    "DEBUG: Segment %s pattern %s: f_stat=%.3f",
                    segment_id, tuple(bits), f_stat_val,
                )
            if f_stat_val >= f_stat_threshold:
                combos_info.append((tuple(bits), f_stat_val, p_val))
            return

        # Option A: cut at cluster k (bit = 1).
        cut_indices = [i for s in slices[g_start : k + 1] for i in s]
        cut_mean = _group_mean(cut_indices)
        if cut_mean <= prev_mean:
            groups.append(cut_indices)
            bits.append(1)
            _dfs(k + 1, k + 1, cut_mean, groups, bits)
            groups.pop()
            bits.pop()

        # Option B: no cut at cluster k (bit = 0).
        bits.append(0)
        _dfs(k + 1, g_start, prev_mean, groups, bits)
        bits.pop()

    _dfs(0, 0, float("inf"), [], [])

    # 6) If none passed, skip
    if not combos_info:
        if debug:
            logger.debug(
                "DEBUG: [%s] No PAS combinations met the F-stat "
                "threshold; skipping segment.",
                segment_id,
            )
        return None

    # 7) Build the union pattern (cluster-level) and compute drops.
    union_pattern = [
        int(any(pat[i] for pat, *_ in combos_info)) for i in range(m)
    ]

    # Reconstruct groups by scanning subsegments in order.  A group is
    # closed when the last PAS of a cut-point cluster is encountered.
    groups, current, pas_idx = [], [], 0
    for cov_idx, pid in enumerate(pas_ids):
        if pid == ".":
            current.append(cov_idx)
        else:
            cid = cluster_id_per_pas[pas_idx]
            is_last_in_cluster = (
                pas_idx == m_orig - 1
                or cluster_id_per_pas[pas_idx + 1] != cid
            )
            current.append(cov_idx)
            if is_last_in_cluster and union_pattern[cid] == 1:
                groups.append(current)
                current = []
            pas_idx += 1
    if current:
        groups.append(current)

    group_means = [_group_mean(g) if g else 0.0 for g in groups]
    if debug:
        logger.debug(
            "DEBUG: Segment %s union group means: %s",
            segment_id, group_means,
        )

    if len(group_means) >= 2:
        cluster_drops = [
            group_means[i] - group_means[i + 1]
            for i in range(len(group_means) - 1)
        ]
    else:
        cluster_drops = [group_means[0]]  # pragma: no cover
    cluster_drops.append(0.0)
    sum_cluster_drops = sum(cluster_drops)

    union_mono = int(
        sum_cluster_drops > 0 and all(d >= 0 for d in cluster_drops[:-1])
    )

    # Map each cut-point cluster to its drop value.
    drop_per_cluster: dict[int, float] = {}
    drop_idx = 0
    for cid, flag in enumerate(union_pattern):
        if flag:
            drop_per_cluster[cid] = cluster_drops[drop_idx]
            drop_idx += 1
        else:
            drop_per_cluster[cid] = 0.0

    # Distribute each cluster's drop to its individual PAS proportionally
    # by atlas_rpm.  Fall back to equal shares when all atlas_rpm = 0.
    drop_per_pas: list[float] = []
    for i, cid in enumerate(cluster_id_per_pas):
        members = [j for j, c in enumerate(cluster_id_per_pas) if c == cid]
        total_atlas = sum(atlas_rpms[j] for j in members)
        weight = (
            atlas_rpms[i] / total_atlas
            if total_atlas > 0
            else 1.0 / len(members)
        )
        drop_per_pas.append(drop_per_cluster[cid] * weight)

    sum_union_drops = sum(drop_per_pas)
    usage_per_pas = [
        (d / sum_union_drops if sum_union_drops > 0 else 0.0)
        for d in drop_per_pas
    ]
    _, best_f, best_p = max(combos_info, key=lambda x: x[1])

    logger.debug(
        "DEBUG: Segment %s union_pattern (clusters): %s",
        segment_id, union_pattern,
    )
    logger.debug(
        "DEBUG: Segment %s cluster_drops: %s", segment_id, cluster_drops,
    )
    logger.debug(
        "DEBUG: Segment %s drop_per_pas: %s", segment_id, drop_per_pas,
    )
    logger.debug(
        "DEBUG: Segment %s union_monotonic: %s", segment_id, union_mono,
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
        "used_combos": [
            {"pattern": list(pat), "f_stat": f_val, "p_value": p_val}
            for pat, f_val, p_val in combos_info
        ],
        "debug_info": {
            "threshold": f_stat_threshold,
            "num_kept_combos": len(combos_info),
            "n_original_pas": m_orig,
            "n_clusters": n_clusters,
            "cluster_assignments": cluster_id_per_pas,
        },
        "unique_pas_ids": unique_pas_ids,
    }


def compute_drops_and_usage(group):
    """Compute coverage drops and RNA usage fractions for one segment.

    Drop at position i is mean_cov[i] - mean_cov[i+1]; the last is 0.
    Usage is each drop's share of the total drop.

    Args:
        group: DataFrame slice for one segment with subsegment_id,
            mean_cov, and pas_id columns.

    Returns:
        Input DataFrame with added columns rna_drop_cov,
        rna_sum_drop_cov, rna_monotone, rna_usage (None for ".").
    """
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
    debug=False,
    max_pas_count=10,
    f_stat_threshold=100,
    cluster_distance=0,
):
    """Evaluate PAS usage patterns for one segment.

    Intended for ProcessPoolExecutor workers; accepts plain dicts to
    avoid DataFrame pickling overhead.

    Args:
        segment_tuple: (segment_id, subsegments) pair.
        debug: Forwarded to evaluate_all_pas_usage_patterns.
        max_pas_count: Forwarded to evaluate_all_pas_usage_patterns.
        f_stat_threshold: Forwarded to evaluate_all_pas_usage_patterns.
        cluster_distance: Forwarded to evaluate_all_pas_usage_patterns.

    Returns:
        (segment_id, subsegments, result) triple.
    """
    segment_id, subsegments = segment_tuple
    result = evaluate_all_pas_usage_patterns(
        subsegments,
        segment_id,
        debug=debug,
        max_pas_count=max_pas_count,
        f_stat_threshold=f_stat_threshold,
        cluster_distance=cluster_distance,
    )
    return segment_id, subsegments, result


def _read_coverage_batch(batch, bw_pos_path, bw_neg_path):
    """Read BigWig coverage for a batch of sub-segment rows.

    Opens strand-specific handles once, reads coverage for each row,
    then closes. Safe for ThreadPoolExecutor workers (each thread owns
    its handles; BigWig reads release the GIL).

    Args:
        batch: Iterable of namedtuples with gene_id, segment_id,
            subsegment_id, chrom, start, end, strand, pas_id fields.
        bw_pos_path: Path to the positive-strand BigWig file.
        bw_neg_path: Path to the negative-strand BigWig file.

    Returns:
        List of dicts with gene_id, segment_id, subsegment_id, pas_id,
        mean_cov, sum_squared_values, coverage.
    """
    bw_pos = pyBigWig.open(bw_pos_path)
    bw_neg = pyBigWig.open(bw_neg_path)
    rows_out = []
    for row in batch:
        chrom = row.chrom
        start = row.start
        end = row.end
        strand = row.strand
        bw = bw_pos if strand == "+" else bw_neg
        try:
            coverage = bw.values(chrom, start, end, numpy=True)
            coverage = np.nan_to_num(coverage, nan=0.0)
        except RuntimeError:
            coverage = np.array([])

        mean_cov = np.nanmean(coverage) if len(coverage) > 0 else 0.0
        sum_squared = np.nansum(coverage**2) if len(coverage) > 0 else 0.0
        rows_out.append(
            {
                "gene_id": row.gene_id,
                "segment_id": row.segment_id,
                "subsegment_id": row.subsegment_id,
                "pas_id": row.pas_id,
                "mean_cov": mean_cov,
                "sum_squared_values": sum_squared,
                "coverage": coverage,
            }
        )
    bw_pos.close()
    bw_neg.close()
    return rows_out


class CalculateCoverages:
    """Compute per-sub-segment coverage metrics and PAS usage fractions.

    Reads strand-specific BigWig tracks, evaluates PAS usage models via
    monotone ANOVA F-statistics, and optionally enriches results with
    per-segment expression ranks from a BAM file.

    Attributes:
        coverage_bw_pos: Path to the positive-strand BigWig file.
        coverage_bw_neg: Path to the negative-strand BigWig file.
        f_stat_threshold: Default F-statistic threshold for run().
    """

    def __init__(
        self,
        coverage_bw_pos,
        coverage_bw_neg,
        f_stat_threshold=100,
    ):
        """Args:
            coverage_bw_pos: Path to the positive-strand BigWig file.
            coverage_bw_neg: Path to the negative-strand BigWig file.
            f_stat_threshold: Default minimum F-statistic for run().
        """
        self.coverage_bw_pos = coverage_bw_pos
        self.coverage_bw_neg = coverage_bw_neg
        self.f_stat_threshold = f_stat_threshold

    def calculate_coverage_metrics(self, subsegments_df, n_threads=1):
        """Read BigWig coverage for every sub-segment.

        When n_threads > 1, splits work across threads (BigWig reads
        release the GIL, so threads provide real I/O parallelism).

        Args:
            subsegments_df: DataFrame with gene_id, segment_id,
                subsegment_id, chrom, start, end, strand, pas_id.
            n_threads: Number of threads (1 = sequential).

        Returns:
            DataFrame with one row per sub-segment and columns gene_id,
            segment_id, subsegment_id, pas_id, mean_cov,
            sum_squared_values, coverage.
        """
        rows = list(subsegments_df.itertuples(index=False))

        if n_threads <= 1:
            return pd.DataFrame(
                _read_coverage_batch(
                    rows, self.coverage_bw_pos, self.coverage_bw_neg
                )
            )

        # Split rows into n_threads chunks; submit each to a thread.
        chunk_size = max(1, (len(rows) + n_threads - 1) // n_threads)
        chunks = [
            rows[i : i + chunk_size] for i in range(0, len(rows), chunk_size)
        ]
        all_rows = []
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=n_threads
        ) as executor:
            futures = [
                executor.submit(
                    _read_coverage_batch,
                    chunk,
                    self.coverage_bw_pos,
                    self.coverage_bw_neg,
                )
                for chunk in chunks
            ]
            for fut in futures:  # iterate in submission order → stable output
                all_rows.extend(fut.result())
        return pd.DataFrame(all_rows)

    def evaluate_pas_usage_models(
        self,
        raw_cov_df,
        output_tsv_debug=None,
        n_procs=16,
        max_pas_count=10,
        f_stat_threshold=100,
        debug=False,
        cluster_distance=0,
    ):
        """Evaluate PAS usage patterns across all segments.

        Groups sub-segments by segment_id and runs
        evaluate_all_pas_usage_patterns on each. When n_procs > 1,
        dispatches to a ProcessPoolExecutor; groups are converted to
        dicts first to reduce pickle overhead.

        Args:
            raw_cov_df: DataFrame from calculate_coverage_metrics.
            output_tsv_debug: Optional path to write a JSON debug dump.
            n_procs: Number of worker processes.
            max_pas_count: Skip segments with more clusters than this.
            f_stat_threshold: Minimum F-statistic to retain a pattern.
            debug: Emit verbose per-pattern debug log messages.
            cluster_distance: Forwarded to evaluate_all_pas_usage_patterns.

        Returns:
            DataFrame with columns gene_id, segment_id, subsegment_id,
            pas_id, mean_cov, sum_squared_values, rna_drop_cov,
            rna_sum_drop_cov, rna_monotone, rna_usage, f_stat, p_value.
        """
        usage_rows = []
        debug_rows = []

        # Convert each group to a list of dicts up front to reduce
        # per-process pickle overhead (DataFrame → list[dict]).
        grouped = [
            (seg_id, grp.to_dict("records"))
            for seg_id, grp in raw_cov_df.groupby("segment_id", sort=False)
        ]

        if n_procs > 1:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=n_procs
            ) as executor:
                func = partial(
                    process_segment,
                    debug=debug,
                    max_pas_count=max_pas_count,
                    f_stat_threshold=f_stat_threshold,
                    cluster_distance=cluster_distance,
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
            for segment_id, subsegments in grouped:
                logger.debug("DEBUG: Evaluating segment %s", segment_id)
                result = evaluate_all_pas_usage_patterns(
                    subsegments,
                    segment_id,
                    debug=debug,
                    max_pas_count=max_pas_count,
                    f_stat_threshold=f_stat_threshold,
                    cluster_distance=cluster_distance,
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
            with open(output_tsv_debug, "w") as fh:
                json.dump(debug_rows, fh, indent=4, default=json_serial)

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
        usage_df = pd.DataFrame(usage_rows)
        return usage_df.reindex(columns=cols)

    def run(
        self,
        subsegments_df,
        output_debug_json=None,
        n_threads=1,
        max_pas_count=10,
        f_stat_threshold=None,
        debug=False,
        cluster_distance=0,
    ):
        """Run the full coverage and PAS usage pipeline.

        Computes coverage metrics and evaluates PAS usage models. Returns
        a slimmed coverage table (coverage arrays dropped) alongside the
        usage DataFrame.

        Args:
            subsegments_df: Sub-segment coordinate table.
            output_debug_json: Optional path for a JSON debug dump.
            n_threads: Threads for BigWig reading and F-stat workers.
            max_pas_count: Skip segments with more clusters than this.
            f_stat_threshold: F-statistic threshold; falls back to
                self.f_stat_threshold when None.
            debug: Emit verbose per-pattern debug log messages.
            cluster_distance: Maximum centre-to-centre distance (bp) to
                merge consecutive PAS into one cluster for ANOVA.
                0 = no clustering.

        Returns:
            (raw_cov_slim, usage_df) tuple.
        """
        start_time = time.time()

        # Step 1: initial coverage metrics
        logger.info("Step 1: Calculating initial coverage metrics...")
        raw_cov_df = self.calculate_coverage_metrics(
            subsegments_df, n_threads=n_threads
        )

        # Merge PAS coordinates and atlas_rpm so the DFS can cluster PAS
        # and distribute within-cluster drops proportionally.
        coord_cols = [
            c for c in ("pas_start", "pas_end", "atlas_rpm")
            if c in subsegments_df.columns
        ]
        if coord_cols:
            raw_cov_df = raw_cov_df.merge(
                subsegments_df[["subsegment_id"] + coord_cols]
                .drop_duplicates("subsegment_id"),
                on="subsegment_id",
                how="left",
            )

        # Step 2: PAS usage via F-statistics
        if cluster_distance > 0:
            logger.info(
                "Step 2: Evaluating PAS usage models via F-statistics "
                "(clustering PAS within %d bp)...",
                cluster_distance,
            )
        else:
            logger.info(
                "Step 2: Evaluating PAS usage models via F-statistics..."
            )
        thr = (
            f_stat_threshold
            if f_stat_threshold is not None
            else self.f_stat_threshold
        )
        refined_usage_df = self.evaluate_pas_usage_models(
            raw_cov_df,
            output_debug_json,
            n_procs=n_threads,
            max_pas_count=max_pas_count,
            f_stat_threshold=thr,
            debug=debug,
            cluster_distance=cluster_distance,
        )

        # Return a slimmed raw-coverage table (no large coverage arrays)
        # for callers that need per-subsegment mean coverage (e.g. BigWig
        # writing).  Drop the array column before freeing the full frame.
        raw_cov_slim = raw_cov_df.drop(columns=["coverage"], errors="ignore")
        del raw_cov_df

        end_time = time.time()
        logger.info("PAS usage evaluation complete.")
        logger.info("Total runtime: %.2f seconds", end_time - start_time)

        return raw_cov_slim, refined_usage_df
