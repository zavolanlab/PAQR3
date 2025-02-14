# PAQR3: Poly(A) site Quantification on RNA-Seq data

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
usage: paqr3 [-h] [--strandedness {true,false}] [--downstream_exon_extension INT]
             [--version]
             --annotation ANNOTATION --pas_atlas PAS_ATLAS --output_regions OUTPUT_REGIONS
             --coverage COVERAGE --output_coverage OUTPUT_COVERAGE

Filter and construct genomic segments from RNA-Seq annotations.

options:
  -h, --help            show this help message and exit
  -v, --version,        show version information and exit
  --strandedness {true,false}
                        Whether the data is stranded (true) or unstranded
                        (false). (default: true)
  --downstream_exon_extension INT
                        Number of bases to extend terminal exons. (default:
                        200)
  --annotation ANNOTATION
                        Path to the annotation GTF file.
  --pas_atlas PAS_ATLAS   Path to the PAS atlas BED file.
  --output_regions OUTPUT_REGIONS
                        Path to the output GTF file for regions and segments.
  --coverage COVERAGE   Path to the coverage BED file.
  --output_coverage OUTPUT_COVERAGE
                        Path to output TSV file for coverage results.
```

Output Files

    output_regions.gtf: A GTF file containing the processed genomic regions, constructed segments, and identified PAS sites.
    output_coverage.tsv: A TSV file containing the calculated mean coverage for each segment.

## Contributing

Contributions are welcome! Please see the contributing guidelines for more information.

## License

GNU General Public License v2.0