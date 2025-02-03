from models import Gene, Transcript, Region
from pybedtools import BedTool
from datetime import datetime

# Utility function for timestamped messages
def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def add_intronic_regions(genes):
    """
    Add intronic regions to the existing annotation.

    Args:
        genes (dict): A dictionary of Gene objects.

    Returns:
        dict: Updated dictionary of Gene objects with intronic regions added.
    """
    for gene in genes.values():
        for transcript in gene.transcripts.values():
            transcript.regions.sort(key=lambda r: r.start)  # Ensure exons are sorted
            reconstructed_regions = []
            intron_count = 0

            for i, region in enumerate(transcript.regions):
                reconstructed_regions.append(region)
                if i < len(transcript.regions) - 1:
                    next_region = transcript.regions[i + 1]
                    intron_start = region.end + 1
                    intron_end = next_region.start - 1
                    if intron_start < intron_end:
                        intron_count += 1
                        intron_region = Region(
                            region_type="intron",
                            chrom=region.chrom,
                            start=intron_start,
                            end=intron_end,
                            strand=region.strand,
                            intron_number=intron_count,
                            attributes={"gene_id": gene.gene_id, "transcript_id": transcript.transcript_id},
                        )
                        reconstructed_regions.append(intron_region)

            transcript.regions = reconstructed_regions

    return genes

def construct_regions(genes, extension=200):
    """
    Construct contiguous regions for each gene based on exon and intron boundaries.

    Args:
        genes (dict): Dictionary of Gene objects.
        extension (int): Number of bases to extend terminal exons, default is 200.

    Returns:
        dict: Updated dictionary of Gene objects with constructed regions.
    """
    log_message("Constructing regions with terminal exon extensions...")

    for gene in genes.values():
        all_regions = []
        for transcript in gene.transcripts.values():
            for region in transcript.regions:
                all_regions.append((region.start, region.end, region))

        all_regions = sorted(all_regions, key=lambda x: (x[0], x[1]))  # Sort by start, then end
        constructed_regions = []
        for start, end, region in all_regions:
            constructed_regions.append(region)

        # Extend terminal exon
        terminal_exon = max(constructed_regions, key=lambda r: r.end if r.region_type == "exon" else -1)
        if terminal_exon:
            extended_end = terminal_exon.end + extension

            # Ensure no overlap with the start of the next gene
            next_gene_start = min(
                [
                    r.start for g in genes.values() for r in g.transcripts.values() for r in r.regions
                    if r.start > terminal_exon.end and r.chrom == terminal_exon.chrom
                ],
                default=float("inf"),
            )
            extended_end = min(extended_end, next_gene_start - 1)
            terminal_exon.end = extended_end

        gene.transcripts[next(iter(gene.transcripts))].regions = constructed_regions

    return genes