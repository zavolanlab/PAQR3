from paqr3.models import Region
from datetime import datetime

def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def extend_gene_coordinates(genes, gene_extension_length=200):
    """Extends gene coordinates (simplified and robust - FIXED MIN/MAX)."""

    log_message(f"Extending gene coordinates ({gene_extension_length} bp)...")

    sorted_genes = sorted(genes.values(), key=lambda g: (next((r.chrom for t in g.transcripts.values() for r in t.regions), None), min(list(r.start for t in g.transcripts.values() for r in t.regions)) if g.transcripts else float('inf')))

    for i, gene in enumerate(sorted_genes):
        gene_strand = next((t.strand for t in gene.transcripts.values()), None)
        if gene_strand is None:
            continue

        gene_start = min(list(r.start for t in gene.transcripts.values() for r in t.regions)) if gene.transcripts else float('inf')
        gene_end = max(list(r.end for t in gene.transcripts.values() for r in t.regions)) if gene.transcripts else float('-inf')

        if gene_strand == "+":
            potential_end = gene_end + gene_extension_length
            next_gene_start = potential_end

            for other_gene in sorted_genes[i+1:]:
                if next((r.chrom for t in other_gene.transcripts.values() for r in t.regions), None) == next((r.chrom for t in gene.transcripts.values() for r in t.regions), None) and min(list(r.start for t in other_gene.transcripts.values() for r in t.regions)) > gene_end:
                    next_gene_start = min(potential_end, min(list(r.start for t in other_gene.transcripts.values() for r in t.regions)) - 1)
                    break

            gene_end = min(potential_end, next_gene_start)

        elif gene_strand == "-":
            potential_start = gene_start - gene_extension_length
            previous_gene_end = potential_start

            for other_gene in reversed(sorted_genes[:i]):
                if next((r.chrom for t in other_gene.transcripts.values() for r in t.regions), None) == next((r.chrom for t in gene.transcripts.values() for r in t.regions), None) and max(list(r.end for t in other_gene.transcripts.values() for r in t.regions)) < gene_start:
                    previous_gene_end = max(potential_start, max(list(r.end for t in other_gene.transcripts.values() for r in t.regions)) + 1)
                    break

            gene_start = max(potential_start, previous_gene_end)

        gene.attributes["start"] = gene_start
        gene.attributes["end"] = gene_end

    return genes


def extend_exon_downstream(genes, downstream_exon_extension=200):
    """Extends terminal exons, using pre-extended gene coordinates."""

    log_message(f"Extending terminal exons ({downstream_exon_extension} bp)...")

    for gene in genes.values():
        gene_strand = next((t.strand for t in gene.transcripts.values()), None)
        if gene_strand is None:
            continue

        for transcript in gene.transcripts.values():
            if transcript.regions:
                terminal_exon = max(transcript.regions, key=lambda r: r.end if r.region_type == "exon" and r.strand == "+" else r.start if r.region_type == "exon" and r.strand == "-" else -1)

                if terminal_exon.strand == "+":
                    terminal_exon.end = min(terminal_exon.end + downstream_exon_extension, gene.attributes["end"])  # Cap at extended gene end
                elif terminal_exon.strand == "-":
                    terminal_exon.start = max(terminal_exon.start - downstream_exon_extension, gene.attributes["start"])  # Cap at extended gene start

    return genes


def define_exons_introns(genes):
    """Constructs exons and introns for all transcripts."""

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
    """Constructs segments (correctly handles multiple transcripts and terminal exons)."""

    for gene in genes.values():
        strand = next((t.strand for t in gene.transcripts.values()), None)
        if strand is None:
            continue

        all_regions = []
        terminal_exon_starts = set()
        terminal_exon_ends = set()

        for transcript in gene.transcripts.values():
            for region in transcript.regions:
                all_regions.append((region.start, region.end, region))
                if region.region_type == "exon":
                    is_terminal = (region.strand == "+" and region.end == max(r.end for r in transcript.regions if r.region_type == "exon")) or \
                                  (region.strand == "-" and region.start == min(r.start for r in transcript.regions if r.region_type == "exon"))
                    if is_terminal:
                        if region.strand == "+":
                            terminal_exon_ends.add(region.end)
                        else:
                            terminal_exon_starts.add(region.start)

        all_regions.sort(key=lambda x: (x[0], x[1]))

        segment_boundaries = set()  # Use a set for efficient boundary management

        for start, end, _ in all_regions:
            segment_boundaries.add(start)
            segment_boundaries.add(end + 1)

        # Remove boundaries that are ONLY terminal exon starts/ends
        segment_boundaries -= terminal_exon_ends
        segment_boundaries -= terminal_exon_starts


        segment_boundaries = sorted(list(segment_boundaries))

        segments = []
        for i in range(len(segment_boundaries) - 1):
            start = segment_boundaries[i]
            end = segment_boundaries[i + 1] - 1

            if start <= end:
                segments.append(
                    Region(
                        region_type="segment",
                        chrom=all_regions[0][2].chrom,
                        start=start,
                        end=end,
                        strand=strand,
                        attributes={"gene_id": gene.gene_id, "segment_strand": strand},
                    )
                )

        if strand == "-":
            segments.reverse()

        for i, segment in enumerate(segments):
            segment.attributes["segment_number"] = i + 1

        gene.segments = segments

    return genes
