# PAQR3: Poly(A) site Quantification on RNA-Seq data

[![license][badge-license]][badge-url-license]

PAQR3 is a command-line tool for quantifying poly(A) sites from standard RNA-Seq data. 
It processes annotation files (GTF), PAS atlases (BED), and coverage data (BED) to identify and quantify poly(A) sites associated with genomic segments.

## Installation

### Using conda (recommended)

```bash
conda create -f install/environment.yml
conda activate paqr3
```

This will install PAQR3 and its dependencies.


## Usage

The main entry point for PAQR3 is the `paqr3` command.  Here's how to use it:

```bash
usage: paqr3 [-h] --annotation ANNOTATION [--downstream_exon_extension DOWNSTREAM_EXON_EXTENSION] --pas_atlas PAS_ATLAS --coverage COVERAGE COVERAGE --output_dir OUTPUT_DIR [--version]

Filter and construct genomic segments from RNA-Seq annotations.

options:
  -h, --help            show this help message and exit
  -v, --version,        show version information and exit
  -de, --downstream_exon_extension INT
                        Number of bases to extend terminal exons. (default:
                        200)
  -a, --annotation ANNOTATION
                        Path to the annotation GTF file.
  -pa, --pas_atlas PAS_ATLAS   Path to the PAS atlas BED file.
  -c, --coverage COVERAGE   Path to the coverage BigWig files. For every sample, a positive and negative stranded BW is required.
  -o, --output_dir OUTPUT_DIRECTORY
                        Path to the output directory.
  -md, --merge-distance INT
                        Merge PAS that are within this distance in base pairs (default: 5). Set to 0 to disable merging.
```

Output Files

    output_regions.gtf: A GTF file containing the processed genomic regions, constructed segments, and identified PAS sites.
    output_genes.tsv: A TSV file containing the unique gene ID for every gene in the input annotation.
    output_PAS.tsv: A TSV file containing the unique gene ID for every PAS in the input PAS atlas.
    output_segments.tsv: A TSV file contining the segments constructed based on the transcripts isoforms for each gene.
    output_subsegments.tsv: A TSV file containing the subsegments constructed based on segments with overlapping PAS.
    output_coverage.tsv: A TSV file containing the calculated mean coverage and the sum of squared coverage values for each segment.

## Contributing

Contributions are welcome! Please see the contributing guidelines for more information.

## License

GNU General Public License v2.0

[badge-license]: <https://img.shields.io/badge/License-GPL_v2-blue.svg>
[badge-url-license]: <https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html>