"""Coverage metric calculation for PAQR3 poly-A site usage analysis.

This module defines :class:`CalculateCoverages`, which reads RNA-seq
coverage BigWig files, computes mean and sum-squared coverage per
sub-segment, evaluates poly-A site (PAS) usage models using an F-
statistic criterion, and writes per-PAS usage tables.

Key functions and methods:

- :func:`evaluate_all_pas_usage_patterns`: searches all monotone PAS
  assignment patterns using a pruned DFS, selecting those that meet
  an F-statistic threshold.
- :meth:`~CalculateCoverages.calculate_coverage_metrics`: reads
  BigWig coverage per sub-segment.
- :meth:`~CalculateCoverages.evaluate_pas_usage_models`: parallelises
  pattern evaluation across segments.
- :meth:`~CalculateCoverages.run`: full pipeline entry point.
"""

import concurrent.futures
import json
import logging
import time
import warnings
from functools import partial

import numpy as np
import pandas as pd  # type: ignore
import pyBigWig  # type: ignore
import pysam  # type: ignore
from pybedtools import BedTool  # type: ignore
from scipy.stats import f  # type: ignore

warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")


logger = logging.getLogger(__name__)


def json_serial(obj):
    """JSON serialiser for NumPy scalar and array types.

    Used as the ``default`` argument to :func:`json.dump` when writing
    debug output that may contain NumPy values.

    Args:
        obj: Object to serialise.

    Returns:
        A JSON-compatible Python object (``float``, ``int``, or
        ``list``).

    Raises:
        TypeError: If *obj* is not a recognised NumPy type.
    """
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
    debug=False,
    max_pas_count=10,
    f_stat_threshold=100,
):
    """Select PAS usage patterns satisfying monotonicity and F-stat criteria.

    Searches all binary patterns that assign each PAS to a group
    boundary, keeping only those where coverage is non-increasing
    across groups (monotone) and whose one-way ANOVA F-statistic is at
    least *f_stat_threshold*.  A "union" pattern is built from the OR
    of all passing patterns, and per-PAS usage fractions are derived
    from the resulting coverage drops.

    **Performance improvements**

    - *Pruned DFS* (point 14): rather than enumerating all 2^m binary
      strings with :func:`itertools.product`, a depth-first search
      abandons any branch the moment the current group's mean exceeds
      the previous group's mean.  In practice most patterns are non-
      monotone, so the majority of the search tree is pruned before
      any F-stat is computed.
    - *Memoised concatenation* (point 15): group concatenated-array
      and mean values are cached by ``frozenset`` of coverage indices;
      groups shared across different patterns are computed only once.

    Args:
        subsegments: List of sub-segment records, each a dict with at
            least keys ``"coverage"`` (1-D NumPy array) and
            ``"pas_id"`` (string; ``"."`` for non-PAS positions).
        segment_id: Identifier used in debug log messages.
        debug: If ``True``, emit detailed per-pattern log messages.
        max_pas_count: Maximum number of unique PAS allowed; segments
            with more are skipped (returns ``None``).
        f_stat_threshold: Minimum F-statistic for a pattern to be
            retained.

    Returns:
        A result dict with keys:

        - ``"pas_usage"``: ``{pas_id: usage_fraction}``
        - ``"f_stat"``: best F-statistic across retained patterns
        - ``"p_value"``: corresponding p-value
        - ``"rna_drop_cov"``: coverage drop per PAS position
        - ``"rna_sum_drop_cov"``: total coverage drop
        - ``"rna_monotone"``: 1 if the union pattern is monotone
        - ``"used_combos"``: list of retained patterns (tuples)
        - ``"debug_info"``: summary string or dict
        - ``"unique_pas_ids"``: ordered list of PAS IDs

        Returns ``None`` if no pattern passes the threshold or the
        segment should be skipped.
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
                "DEBUG: [%s] Single PAS '%s', mean_cov=%.3f "
                "→ usage=%.3f",
                segment_id,
                pid,
                mu,
                usage,
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
        if debug:
            logger.debug(
                "DEBUG: [%s] Too many PAS (%d) in this segment; "
                "only segments with %d or fewer PAS are evaluated. "
                "Skipping.",
                segment_id,
                len(unique_pas_ids),
                max_pas_count,
            )
        return None

    # 4) No PAS at all? skip
    if not unique_pas_ids:
        if debug:
            logger.debug(
                "DEBUG: [%s] No PAS found in this segment, skipping.",
                segment_id,
            )
        return None

    # 5) Evaluate every monotone pattern via pruned DFS.
    #
    # Build per-PAS "slices": slices[k] contains coverage indices from
    # immediately after PAS k-1 up to and including PAS k.  tail holds
    # any indices that follow the last PAS.  A group is a contiguous
    # run of slices, so its concatenated array and mean can be memoised
    # by frozenset of coverage indices — groups shared across different
    # patterns are computed only once.
    m = len(unique_pas_ids)
    slices = []
    current_slice = []
    for cov_idx, pid in enumerate(pas_ids):
        current_slice.append(cov_idx)
        if pid != ".":
            slices.append(current_slice)
            current_slice = []
    tail = current_slice

    # Cache: frozenset(cov_indices) → concatenated coverage array.
    # Single-element groups bypass the cache to avoid dict overhead.
    _concat_cache = {}

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
        """DFS over binary cut/no-cut decisions at each PAS position.

        Prunes any branch where the newly closed group's mean exceeds
        the previous group's mean, cutting off non-monotone sub-trees
        before their F-statistics are computed.

        Args:
            k: Index of the PAS slice currently being considered.
            g_start: Slice index at which the current open group
                began.
            prev_mean: Mean of the last closed group; ``float("inf")``
                at the root so any first group is accepted.
            groups: Closed groups accumulated so far (mutated in
                place; restored on backtrack).
            bits: Pattern bits accumulated so far (mutated in place;
                restored on backtrack).
        """
        if k == len(slices):
            if not groups:
                return  # all-zeros pattern — skip
            last_group = [
                i for s in slices[g_start:] for i in s
            ] + tail
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
                    segment_id,
                    tuple(bits),
                    f_stat_val,
                )
            if f_stat_val >= f_stat_threshold:
                combos_info.append((tuple(bits), f_stat_val, p_val))
            return

        # Option A: cut at PAS k (bit = 1).
        # Closes the current group as the union of slices[g_start..k].
        cut_indices = [i for s in slices[g_start : k + 1] for i in s]
        cut_mean = _group_mean(cut_indices)
        if cut_mean <= prev_mean:  # monotone prefix — explore
            groups.append(cut_indices)
            bits.append(1)
            _dfs(k + 1, k + 1, cut_mean, groups, bits)
            groups.pop()
            bits.pop()

        # Option B: no cut at PAS k (bit = 0).
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

    # Reuse the memoised _group_mean for the union groups.
    group_means = [_group_mean(g) if g else 0.0 for g in groups]
    if debug:
        logger.debug(
            "DEBUG: Segment %s union group means: %s",
            segment_id,
            group_means,
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
    _, best_f, best_p = max(combos_info, key=lambda x: x[1])

    logger.debug(
        "DEBUG: Segment %s union_pattern: %s", segment_id, union_pattern
    )
    logger.debug(
        "DEBUG: Segment %s union_drops: %s", segment_id, union_drops
    )
    logger.debug(
        "DEBUG: Segment %s union_monotonic: %s", segment_id, union_mono
    )
    logger.debug(
        "DEBUG: Segment %s usage_per_pas: %s", segment_id, usage_per_pas
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
    """Compute per-sub-segment coverage drops and RNA usage fractions.

    Rows are sorted by ``subsegment_id`` before processing.  The
    *drop* at position *i* is ``mean_cov[i] - mean_cov[i+1]``; the
    last position always has a drop of 0.  Usage for each sub-segment
    is that drop's share of the total drop.

    Args:
        group: DataFrame slice for a single segment, containing at
            minimum columns ``subsegment_id``, ``mean_cov``,
            ``pas_id``.

    Returns:
        The input DataFrame augmented with columns ``rna_drop_cov``,
        ``rna_sum_drop_cov``, ``rna_monotone``, and ``rna_usage``.
        Sub-segments whose ``pas_id`` is ``"."`` receive
        ``rna_usage = None``.
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
):
    """Evaluate PAS usage patterns for a single segment.

    Designed to be called from a
    :class:`~concurrent.futures.ProcessPoolExecutor` worker.  Accepts
    pre-converted sub-segment records (plain dicts with NumPy arrays)
    rather than a :class:`~pandas.DataFrame` to avoid the significant
    overhead of pickling a DataFrame across process boundaries.

    Args:
        segment_tuple: A ``(segment_id, subsegments)`` pair where
            *subsegments* is a list of dicts, each containing at least
            the keys ``"coverage"``, ``"pas_id"``, ``"gene_id"``,
            ``"segment_id"``, ``"subsegment_id"``, ``"mean_cov"``, and
            ``"sum_squared_values"``.
        debug: Forwarded to
            :func:`evaluate_all_pas_usage_patterns`.
        max_pas_count: Forwarded to
            :func:`evaluate_all_pas_usage_patterns`.
        f_stat_threshold: Forwarded to
            :func:`evaluate_all_pas_usage_patterns`.

    Returns:
        A ``(segment_id, subsegments, result)`` triple where *result*
        is the dict from :func:`evaluate_all_pas_usage_patterns` or
        ``None``.
    """
    segment_id, subsegments = segment_tuple
    result = evaluate_all_pas_usage_patterns(
        subsegments,
        segment_id,
        debug=debug,
        max_pas_count=max_pas_count,
        f_stat_threshold=f_stat_threshold,
    )
    return segment_id, subsegments, result


class CalculateCoverages:
    """Compute RNA-seq coverage metrics and PAS usage fractions.

    Reads strand-specific BigWig coverage tracks, computes per-sub-
    segment mean and sum-squared coverage, evaluates poly-A site (PAS)
    usage models via monotone ANOVA F-statistics, and optionally
    enriches results with per-segment expression rank information from
    a BAM file.

    Attributes:
        coverage_bw_pos: Path to the positive-strand BigWig file.
        coverage_bw_neg: Path to the negative-strand BigWig file.
        bam_file: Optional path to a BAM file for read-count-based
            expression statistics.
        f_stat_threshold: Default F-statistic threshold used when
            :meth:`run` is called without an explicit threshold.
    """

    def __init__(
        self,
        coverage_bw_pos,
        coverage_bw_neg,
        bam_file=None,
        f_stat_threshold=100,
    ):
        """Initialise with BigWig paths and optional BAM file.

        Args:
            coverage_bw_pos: Path to the positive-strand BigWig file.
            coverage_bw_neg: Path to the negative-strand BigWig file.
            bam_file: Optional BAM file path; required only when
                :meth:`calculate_segment_expression_stats` is called.
            f_stat_threshold: Default minimum F-statistic for
                :meth:`run`.
        """
        self.coverage_bw_pos = coverage_bw_pos
        self.coverage_bw_neg = coverage_bw_neg
        self.bam_file = bam_file
        self.f_stat_threshold = f_stat_threshold

    def calculate_coverage_metrics(self, subsegments_df):
        """Read BigWig coverage for every sub-segment.

        Opens the strand-specific BigWig files once, then reads
        coverage values for each sub-segment row using
        :meth:`pyBigWig.pyBigWig.values`.  Rows are iterated with
        :meth:`~pandas.DataFrame.itertuples` (3–5× faster than
        :meth:`~pandas.DataFrame.iterrows`, which boxes each row into
        a :class:`~pandas.Series` object).

        Args:
            subsegments_df: DataFrame with columns
                ``gene_unique_id``, ``segment_number``,
                ``subsegment_number``, ``chrom``, ``start``, ``end``,
                ``strand``, ``pas_id``.

        Returns:
            A new DataFrame with one row per sub-segment and columns
            ``gene_id``, ``segment_id``, ``subsegment_id``,
            ``pas_id``, ``mean_cov``, ``sum_squared_values``,
            ``coverage``.
        """
        bw_pos = pyBigWig.open(self.coverage_bw_pos)
        bw_neg = pyBigWig.open(self.coverage_bw_neg)

        all_rows = []
        for row in subsegments_df.itertuples(index=False):
            gene_id = row.gene_unique_id
            seg_num = row.segment_number
            subseg_num = row.subsegment_number
            chrom = row.chrom
            start = row.start
            end = row.end
            strand = row.strand
            pas_id = row.pas_id
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
        debug=False,
    ):
        """Evaluate PAS usage patterns across all segments.

        Groups sub-segments by ``segment_id``, then evaluates
        :func:`evaluate_all_pas_usage_patterns` for each group.  When
        ``n_procs > 1``, work is dispatched to a
        :class:`~concurrent.futures.ProcessPoolExecutor`.

        Groups are converted to plain lists of dicts before
        dispatching to reduce pickle overhead.  Serialising a
        :class:`~pandas.DataFrame` across process boundaries is
        significantly more expensive than serialising a list of dicts
        with NumPy arrays.

        .. note::
            For further optimisation, consider splitting I/O-bound
            BigWig reads (:meth:`calculate_coverage_metrics`) onto a
            :class:`~concurrent.futures.ThreadPoolExecutor` (BigWig
            reads release the GIL) and reserving
            :class:`~concurrent.futures.ProcessPoolExecutor` only for
            the CPU-bound F-statistic computation here.

        Args:
            raw_cov_df: DataFrame as returned by
                :meth:`calculate_coverage_metrics`, with a
                ``"segment_id"`` column.
            output_tsv_debug: Optional path to write a JSON debug
                dump of all evaluated patterns.
            n_procs: Number of worker processes.
            max_pas_count: Forwarded to
                :func:`evaluate_all_pas_usage_patterns`.
            f_stat_threshold: Forwarded to
                :func:`evaluate_all_pas_usage_patterns`.
            debug: Emit verbose per-pattern debug log messages.

        Returns:
            A DataFrame with columns ``gene_id``, ``segment_id``,
            ``subsegment_id``, ``pas_id``, ``mean_cov``,
            ``sum_squared_values``, ``rna_drop_cov``,
            ``rna_sum_drop_cov``, ``rna_monotone``, ``rna_usage``,
            ``f_stat``, ``p_value``.
        """
        usage_rows = []
        debug_rows = []

        # Convert each group to a list of dicts up front to reduce
        # per-process pickle overhead (DataFrame → list[dict]).
        grouped = [
            (seg_id, grp.to_dict("records"))
            for seg_id, grp in raw_cov_df.groupby(
                "segment_id", sort=False
            )
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
                                        "subsegment_id": s[
                                            "subsegment_id"
                                        ],
                                        "pas_id": pid,
                                        "mean_cov": s["mean_cov"],
                                        "sum_squared_values": s[
                                            "sum_squared_values"
                                        ],
                                        "rna_drop_cov": result[
                                            "rna_drop_cov"
                                        ][drop_idx],
                                        "rna_sum_drop_cov": result[
                                            "rna_sum_drop_cov"
                                        ],
                                        "rna_monotone": result[
                                            "rna_monotone"
                                        ],
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
                logger.debug(
                    "DEBUG: Evaluating segment %s", segment_id
                )
                result = evaluate_all_pas_usage_patterns(
                    subsegments,
                    segment_id,
                    debug=debug,
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
                                    "rna_drop_cov": result[
                                        "rna_drop_cov"
                                    ][drop_idx],
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
        """Write a coverage results DataFrame to a tab-separated file.

        The ``"coverage"`` column is dropped before writing to keep
        output file sizes manageable.

        Args:
            results_df: DataFrame to write; may contain a
                ``"coverage"`` column.
            output_tsv: Destination file path.
        """
        if "coverage" in results_df.columns:
            results_df = results_df.drop(columns=["coverage"])
        results_df.to_csv(output_tsv, sep="\t", index=False)
        logger.info("Results written to %s", output_tsv)

    def calculate_segment_expression_stats(self, usage_df, subsegments_df):
        """Compute per-segment expression ranks and merge onto usage table.

        Derives four expression metrics per segment:

        - **rpm**: read count divided by segment length (from BAM).
        - **rpm_rank**: dense rank of ``rpm`` (1 = highest).
        - **median_cov**: median BigWig coverage across the segment.
        - **median_rank**: dense rank of ``median_cov`` (1 = highest).

        Args:
            usage_df: Per-PAS usage DataFrame as returned by
                :meth:`evaluate_pas_usage_models`.
            subsegments_df: Sub-segment coordinate table used to
                derive segment boundaries.

        Returns:
            *usage_df* merged with four new columns: ``rpm``,
            ``rpm_rank``, ``median_cov``, ``median_rank``.

        Raises:
            ValueError: If ``self.bam_file`` is not set.
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
            seg_meta["rpm"]
            .rank(method="dense", ascending=False)
            .astype(int)
        )
        seg_meta["median_rank"] = (
            seg_meta["median_cov"]
            .rank(method="dense", ascending=False)
            .astype(int)
        )

        # 5) Merge back onto the PAS-level usage_df
        merged = usage_df.merge(
            seg_meta[
                [
                    "segment_id",
                    "rpm",
                    "rpm_rank",
                    "median_cov",
                    "median_rank",
                ]
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
        debug=False,
    ):
        """Execute the full coverage and PAS usage pipeline.

        Steps, in order:

        1. Compute per-sub-segment coverage metrics from BigWig files.
        2. Apply simple drop-based PAS usage and write the raw TSV.
        3. Evaluate PAS usage models via F-statistics (parallelised
           across segments when ``n_procs > 1``).
        4. Optionally compute per-segment expression ranks from BAM.
        5. Write the final per-PAS usage TSV.

        Args:
            subsegments_df: Sub-segment coordinate table from
                :meth:`~paqr3.construct_segments.ConstructSegments\
.create_subsegments_dataframe`.
            output_raw_tsv: Path for the intermediate coverage TSV.
            output_final_tsv: Path for the final per-PAS usage TSV.
            output_debug_json: Optional path for a JSON debug dump of
                all evaluated patterns.
            n_procs: Number of worker processes for parallel
                F-statistic evaluation.
            max_pas_count: Passed to
                :func:`evaluate_all_pas_usage_patterns`; segments with
                more PAS are skipped.
            f_stat_threshold: F-statistic threshold; falls back to
                ``self.f_stat_threshold`` when ``None``.
            debug: Emit verbose per-pattern debug log messages.

        Returns:
            The final per-PAS usage DataFrame.
        """
        start_time = time.time()

        # Step 1: initial coverage metrics
        logger.info("Step 1: Calculating initial coverage metrics...")
        raw_cov_df = self.calculate_coverage_metrics(subsegments_df)

        grouped = raw_cov_df.groupby(
            ["gene_id", "segment_id"], group_keys=False
        )

        results = [
            compute_drops_and_usage(group)
            for _, group in grouped
            if not group["mean_cov"].isna().all()
        ]

        # Filter out empty or all-NA DataFrames before concatenation
        cov_df = (
            pd.concat(
                [
                    r
                    for r in results
                    if not r.empty and not r.isna().all().all()
                ],
                ignore_index=True,
            )
            if results
            else pd.DataFrame(
                columns=raw_cov_df.columns.tolist()
                + [
                    "rna_drop_cov",
                    "rna_sum_drop_cov",
                    "rna_monotone",
                    "rna_usage",
                ]
            )
        )

        self.write_coverage_results(cov_df.copy(), output_raw_tsv)

        # Step 2: PAS usage via F-statistics
        logger.info(
            "Step 2: Evaluating PAS usage models via F-statistics..."
        )
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
            debug=debug,
        )

        # free memory from raw coverage arrays
        del raw_cov_df
        # Step 3: per-segment expression stats (if BAM provided)
        if self.bam_file:
            refined_usage_df = self.calculate_segment_expression_stats(
                refined_usage_df, subsegments_df
            )

        # Final write
        self.write_coverage_results(refined_usage_df, output_final_tsv)

        end_time = time.time()
        logger.info("PAS usage evaluation complete.")
        logger.info("Total runtime: %.2f seconds", end_time - start_time)

        return refined_usage_df
