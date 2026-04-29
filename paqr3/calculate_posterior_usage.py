"""Posterior PAS usage computation for PAQR3.

Blends RNA-seq-derived drop coverage with atlas RPM values to produce
a smoothed posterior relative usage estimate per PAS per segment.
"""

import logging

import pandas as pd  # type: ignore

logger = logging.getLogger(__name__)


class CalculatePosteriorUsage:
    """Blend atlas RPM with observed drop RPM for posterior PAS usage.

    ``atlas_rpm`` values are expected to already be present as a column
    in the usage DataFrame (carried through from the segments TSV).
    The ``pas_id`` column contains representative cleavage-site strings
    (e.g. ``"chr1:65435"``), not integers.

    Attributes:
        weight: Mixing weight for the observed RPM.  The posterior RPM
            is ``weight * observed_rpm + (1 - weight) * atlas_rpm``.
    """

    def __init__(self, weight: float = 0.1) -> None:
        """Initialise with a mixing weight.

        Args:
            weight: Fraction of the posterior assigned to the
                observed (RNA-seq drop) RPM.  Must be in [0, 1].
        """
        self.weight = weight

    def compute(self, usage_df: pd.DataFrame) -> pd.DataFrame:
        """Compute posterior relative usage from a usage DataFrame.

        Expects *usage_df* to contain columns ``rna_drop_cov``,
        ``atlas_rpm``, and ``segment_id``.  ``pas_id`` is a string
        (rep_cs or ``"."``) and is preserved as-is.

        Steps:

        1. Convert ``rna_drop_cov`` to observed RPM (per-total-drop
           normalisation scaled to 1 × 10⁶).
        2. Blend with ``atlas_rpm`` using the configured weight.
        3. Normalise within each segment so usage fractions sum to 1.

        Args:
            usage_df: Per-PAS usage DataFrame as returned by
                :meth:`~paqr3.calculate_cov_metrics\
.CalculateCoverages.evaluate_pas_usage_models`, augmented with an
                ``atlas_rpm`` column.

        Returns:
            A copy of *usage_df* with three additional columns:
            ``observed_rpm``, ``posterior_rpm``,
            ``posterior_rel_usage``.
        """
        df = usage_df.copy()

        # 1) Observed RPM from drop coverage.
        total_drop = df["rna_drop_cov"].sum()
        df["observed_rpm"] = (
            df["rna_drop_cov"] / total_drop * 1e6
            if total_drop > 0
            else 0.0
        )

        # 2) Weighted posterior RPM.
        w = self.weight
        atlas = df["atlas_rpm"].fillna(0.0)
        df["posterior_rpm"] = w * df["observed_rpm"] + (1 - w) * atlas

        # 3) Normalise within each segment.
        df["posterior_rel_usage"] = df.groupby("segment_id")[
            "posterior_rpm"
        ].transform(lambda x: x / x.sum() if x.sum() > 0 else 0.0)

        # 4) Atlas relative usage per segment.
        df["atlas_rel_usage"] = df.groupby("segment_id")[
            "atlas_rpm"
        ].transform(lambda x: x / x.sum() if x.sum() > 0 else 0.0)

        return df
