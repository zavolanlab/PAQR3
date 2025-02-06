from models import Region
from datetime import datetime

def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def extend_exon_downstream(genes, downstream_exon_extension=200):
    """Constructs regions, extending terminal exons per transcript."""

    log_message(f"Constructing regions with terminal exon extensions ({downstream_exon_extension} bp)...")

    for gene in genes.values():
        # Find the longest transcript for boundary checks
        longest_transcript = max(gene.transcripts.values(), key=lambda t: max((r.end - r.start) for r in t.regions) if t.regions else 0) if gene.transcripts else None
        longest_transcript_start = min(r.start for r in longest_transcript.regions) if longest_transcript else float('inf')
        longest_transcript_end = max(r.end for r in longest_transcript.regions) if longest_transcript else float('-inf')

        for transcript in gene.transcripts.values():
            # Extend terminal exon for the current transcript
            if transcript.regions:  # Check if there are any regions (transcripts can be empty)
                terminal_exon = max(transcript.regions, key=lambda r: r.end if r.region_type == "exon" else -1)

                if terminal_exon.strand == "+":
                    extended_end = terminal_exon.end + downstream_exon_extension
                    extended_end = min(extended_end, longest_transcript_end)  # Boundary check
                    terminal_exon.end = extended_end
                elif terminal_exon.strand == "-":
                    extended_start = terminal_exon.start - downstream_exon_extension
                    extended_start = max(extended_start, longest_transcript_start)  # Boundary check
                    terminal_exon.start = extended_start
        # Update gene coordinates based on the *extended* transcripts
        gene_start = min(r.start for t in gene.transcripts.values() for r in t.regions) if gene.transcripts else float('inf')
        gene_end = max(r.end for t in gene.transcripts.values() for r in t.regions) if gene.transcripts else float('-inf')
        for t in gene.transcripts.values():
            for r in t.regions:
                r.attributes["gene_id"] = gene.gene_id
                r.attributes["transcript_id"] = t.transcript_id
        gene.attributes["start"] = gene_start
        gene.attributes["end"] = gene_end

    return genes


def define_exons_introns(genes):
    """Constructs exons and introns for all transcripts."""

    log_message("Defining exons and introns...")

    for gene in genes.values():
        for transcript in gene.transcripts.values():
            transcript.regions.sort(key=lambda r: r.start)  # Ensure regions are sorted

            exons = [r for r in transcript.regions if r.region_type == "exon"]
            introns = []

            for i in range(len(exons) - 1):
                exon1 = exons[i]
                exon2 = exons[i + 1]
                intron_start = exon1.end + 1
                intron_end = exon2.start - 1

                if intron_start <= intron_end:  # Only create if it's a valid intron
                    intron = Region(
                        region_type="intron",
                        chrom=exon1.chrom,
                        start=intron_start,
                        end=intron_end,
                        strand=exon1.strand,
                        intron_number=i + 1,  # Intron numbering starts from 1
                        attributes={"gene_id": gene.gene_id, "transcript_id": transcript.transcript_id},
                    )
                    introns.append(intron)

            transcript.regions = exons + introns  # Combine exons and introns
            transcript.regions.sort(key=lambda r: r.start) # Sort all regions

    return genes


def construct_segments(genes):
    """Constructs segments based on all transcripts of a gene."""

    log_message("Constructing segments...")

    for gene in genes.values():
        # Get strand from any of the transcripts (assuming all transcripts of a gene are on the same strand)
        strand = next((t.strand for t in gene.transcripts.values()), None) # Handle the case where there are no transcripts
        if strand is None: # If there are no transcripts, skip this gene
            continue

        all_regions = []
        for transcript in gene.transcripts.values():
            for region in transcript.regions:  # Use the regions AFTER extension!
                all_regions.append((region.start, region.end, region))

        all_regions.sort(key=lambda x: (x[0], x[1])) 

        segment_boundaries = sorted(list(set([r[0] for r in all_regions] + [r[1] + 1 for r in all_regions])))

        segments = []
        for i in range(len(segment_boundaries) - 1):
            start = segment_boundaries[i]
            end = segment_boundaries[i + 1] - 1
            if start <= end:  # Only create if it's a valid segment
                segments.append(
                    Region(
                        region_type="segment",
                        chrom=all_regions[0][2].chrom,  # Use chromosome from any region
                        start=start,
                        end=end,
                        strand=strand,  # Use strand from transcript
                        attributes={"gene_id": gene.gene_id},
                    )
                )

        if strand == "-":  # Use the strand variable
            segments.reverse()  # Reverse for negative strand

        for i, segment in enumerate(segments):
            segment.attributes["segment_number"] = i + 1

        gene.segments = segments  # Add segments to gene object

    return genes


def write_segments_to_gtf(genes, output_file):
    """Writes genes and segments to a GTF file."""

    with open(output_file, "w") as f:
        for gene in genes.values():
            f.write(gene.to_gtf_format() + "\n")  # Write gene information
            for segment in gene.segments:
                f.write(segment.to_gtf_format() + "\n")  # Write segment information

    log_message(f"Genes and segments written to {output_file}")