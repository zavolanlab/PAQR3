# PAQR3: Poly(A) site Quantification on RNA-Seq data

[![license][badge-license]][badge-url-license]

PAQR3 is a command-line tool for quantifying poly(A) sites from standard RNA-Seq data. 
It processes annotation files (GTF), PAS atlases (BED), and coverage data (BED) to identify and quantify poly(A) sites associated with genomic segments.

## Installation

### Using conda (recommended)

```bash
conda env create -f install/environment.yml
conda activate paqr3
```

This will install PAQR3 and its dependencies.


## Usage

PAQR3 exposes three sub-commands with higly customizable parameters.

### `paqr3 segment` — segmentation only

Parses the GTF annotation, extends terminal exons, merges nearby PAS sites, and writes a segments TSV used as input for `paqr3 quant`.

```
paqr3 segment \
  --annotation/-a                   GTF annotation file (required) \
  --pas-atlas/-pa                   PAS atlas BED file (required) \
  --output-dir/-o                   Output directory (required) \
  --downstream-exon-extension/-de   Bases to extend terminal exons (default: 200) \
  --merge-distance/-md              Merge PAS within this distance in bp (default: 5) \
  --emit/-e                         Additional outputs (see below) \
  --gzip/-g / --no-gzip             Gzip-compress TSV outputs (default: true) \
  --verbosity/-v                    Log level: INFO or DEBUG (default: INFO)
```

### `paqr3 quant` — quantification from an existing segments TSV

```
paqr3 quant \
  --segments-tsv/-s                 Segments TSV from paqr3 segment (required) \
  --coverage-pos/-c-pos             Positive-strand coverage BigWig (required) \
  --coverage-neg/-c-neg             Negative-strand coverage BigWig (required) \
  --output-dir/-o                   Output directory (required) \
  --sample-id/-sid                  Sample name prefix (default: BigWig filename stem) \
  --max-pas-count/-mpc              Skip segments with more PAS than this (default: 10) \
  --f-stat-threshold/-fst           F-statistic threshold for PAS usage (default: 100) \
  --posterior-usage-weight/-puw     Weight for observed RPM in posterior blend (default: 0.1) \
  --bam/-b                          BAM file for expression rank statistics (optional) \
  --threads/-t                      Threads for BigWig reading and F-stat workers (default: 1) \
  --chr-sizes/-cs                   Chromosome-sizes file (required when any BigWig emit mode is set) \
  --emit/-e                         Additional outputs (see below) \
  --gzip/-g / --no-gzip             Gzip-compress TSV outputs (default: true) \
  --verbosity/-v                    Log level: INFO or DEBUG (default: INFO)
```

### `paqr3 full` — segmentation and quantification end-to-end

Accepts all arguments from both `segment` and `quant` (except `--segments-tsv`).

```
paqr3 full \
  --annotation/-a  --pas-atlas/-pa \
  --coverage-pos/-c-pos  --coverage-neg/-c-neg \
  --output-dir/-o  [all other quant/segment options]
```

### `--emit` modes

The `--emit` flag controls which additional outputs are written. Multiple modes can be combined.

| Token | Files written |
|---|---|
| `mean_cov` | `{sample}_mean_cov.bw` — mean coverage per subsegment |
| `observed` | `{sample}_observed_usage.bw` + `{sample}_observed_rpm.bw` |
| `posterior` | `{sample}_posterior_rpm.bw` + `{sample}_posterior_usage.bw` |
| `atlas` | `{sample}_atlas_rpm.bw` + `{sample}_atlas_usage.bw` |
| `all` | All four BigWig groups above |
| `debug` | All BigWigs + debug JSON + four debug BED files |

BigWig emission requires `--chr-sizes`.

## Output files

### `paqr3 segment`

Outputs are written directly to `--output-dir`:

| File | Description |
|---|---|
| `output_segments.tsv` | One row per sub-segment; input for `paqr3 quant`. Columns: `chrom`, `start`, `end`, `subsegment_id`, `strand`, `overlapping_pas_id` (rep CS string or `.`), `atlas_rpm`, `pas_start`, `pas_end` (merged PAS cluster boundaries). |
| `output_segments_debug_{genes,segments,subsegments,pas}.bed` | Four headerless BED files (chrom, start, end, id, strand) — written only when `--emit debug` is set. |

### `paqr3 quant` / `paqr3 full`

Outputs are written to `{output-dir}/{sample}_results/`:

| File | Description |
|---|---|
| `{sample}_segment_results.tsv[.gz]` | One row per evaluated segment. Columns: `chrom`, `start`, `end`, `strand`, `segment_id`, `rna_sum_drop_cov`, `f_stat`, `p_value`. |
| `{sample}_subsegment_results.tsv[.gz]` | One row per subsegment (including trailing). Columns: `chrom`, `start`, `end`, `subsegment_id`, `strand`, `pas_id`, `mean_cov`, `atlas_rpm`, `observed_rpm`, `posterior_rpm`. |
| `{sample}_pas_results.tsv[.gz]` | One row per PAS. Columns: `chrom`, `start`, `end`, `strand`, `pas_id`, `subsegment_id`, `atlas_usage`, `observed_usage`, `posterior_usage`. |
| `{sample}_mean_cov.bw` | Mean RNA-seq coverage per subsegment. Written when `--emit mean_cov`, `all`, or `debug`. |
| `{sample}_observed_usage.bw` | F-stat derived PAS usage fraction, at PAS cluster coordinates. |
| `{sample}_observed_rpm.bw` | Observed coverage-drop RPM per PAS. |
| `{sample}_posterior_rpm.bw` | Posterior blended RPM per PAS. |
| `{sample}_posterior_usage.bw` | Posterior relative usage per PAS (normalised within segment). |
| `{sample}_atlas_rpm.bw` | Atlas RPM per PAS. |
| `{sample}_atlas_usage.bw` | Atlas relative usage per PAS (normalised within segment). |
| `{sample}_debug.json` | Per-segment pattern evaluation dump. Written when `--emit debug`. |
| `{sample}_debug_{genes,segments,subsegments,pas}.bed` | Debug BED files. Written when `--emit debug`. |

`{sample}` defaults to the stem of the positive-strand BigWig filename (everything before the first `.`) and can be overridden with `--sample-id`.  The `.gz` suffix appears when `--gzip` is active (default).

## Contributing

Contributions are welcome! Please see the contributing guidelines for more information.

## License

GNU General Public License v2.0

[badge-license]: <https://img.shields.io/badge/License-GPL_v2-blue.svg>
[badge-url-license]: <https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html>