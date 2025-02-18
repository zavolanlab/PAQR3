class Region:
    def __init__(
        self,
        region_type,
        chrom,
        start,
        end,
        strand,
        exon_number=None,
        intron_number=None,
        attributes=None,
    ):
        self.region_type = region_type  # exon or intron
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.exon_number = exon_number  # For exons
        self.intron_number = intron_number
        self.attributes = attributes or {}  # GTF attributes

    def to_gtf_format(self):
        """Convert the Region object to a GTF format string."""
        attr_str = (
            "; ".join([f'{key} "{value}"' for key, value in self.attributes.items()])
            + ";"
        )
        return f"{self.chrom}\t{self.attributes.get('source', 'unknown')}\t{self.region_type}\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{attr_str}"

    def __str__(self):
        """String representation for debugging."""
        return f"Region({self.region_type}, {self.chrom}:{self.start}-{self.end}, {self.strand}, attributes={self.attributes})"

    def __repr__(self):
        return self.__str__()


class Transcript:
    def __init__(self, transcript_id, strand, attributes=None):
        self.transcript_id = transcript_id
        self.strand = strand
        self.regions = []
        self.attributes = attributes or {}  # GTF attributes

    def add_region(self, region):
        self.regions.append(region)

    def to_gtf_format(self):
        """Convert the Transcript object to a GTF format string."""
        attr_str = (
            "; ".join([f'{key} "{value}"' for key, value in self.attributes.items()])
            + ";"
        )
        start = min(region.start for region in self.regions)
        end = max(region.end for region in self.regions)
        return f"{self.regions[0].chrom}\t{self.attributes.get('source', 'unknown')}\ttranscript\t{start}\t{end}\t.\t{self.strand}\t.\t{attr_str}"

    def __repr__(self):
        return f"Transcript: {self.transcript_id}, Regions: {self.regions}"


class Gene:
    def __init__(self, gene_id, attributes=None):
        self.gene_id = gene_id
        self.transcripts = {}
        self.attributes = attributes or {}  # GTF attributes

    def add_transcript(self, transcript):
        self.transcripts[transcript.transcript_id] = transcript

    def to_gtf_format(self):
        """Convert the Gene object to a GTF format string."""
        attr_str = (
            "; ".join([f'{key} "{value}"' for key, value in self.attributes.items()])
            + ";"
        )
        start = min(
            region.start
            for transcript in self.transcripts.values()
            for region in transcript.regions
        )
        end = max(
            region.end
            for transcript in self.transcripts.values()
            for region in transcript.regions
        )
        strand = next(
            iter(self.transcripts.values())
        ).strand  # Use strand of the first transcript
        chrom = (
            next(iter(self.transcripts.values())).regions[0].chrom
        )  # Use chrom of the first region
        return f"{chrom}\t{self.attributes.get('source', 'unknown')}\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}"

    def __repr__(self):
        return f"Gene: {self.gene_id}, Transcripts: {self.transcripts}"
