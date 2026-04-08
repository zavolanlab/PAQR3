"""Data models for genomic annotation objects.

This module defines the core data structures used to represent genomic
features parsed from GTF files: Region (exon or intron), Transcript,
and Gene. Each class can serialise itself back to a GTF format string
via ``to_gtf_format()``.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass(slots=True)
class Region:
    """A single genomic region (exon or intron).

    Attributes:
        region_type: Feature type, either ``"exon"`` or ``"intron"``.
        chrom: Chromosome or contig name.
        start: 1-based start coordinate (GTF convention).
        end: 1-based end coordinate, inclusive (GTF convention).
        strand: Strand, either ``"+"`` or ``"-"``.
        exon_number: Exon index within the transcript, if applicable.
        intron_number: Intron index within the transcript, if
            applicable.
        attributes: Key-value pairs from the GTF attributes field.
    """

    region_type: str  # "exon" or "intron"
    chrom: str
    start: int
    end: int
    strand: str
    exon_number: Optional[int] = None
    intron_number: Optional[int] = None
    attributes: dict[str, str] = field(default_factory=dict)

    def to_gtf_format(self) -> str:
        """Convert the Region to a GTF format string.

        Returns:
            A tab-delimited GTF line representing this region.
        """
        attr_str = (
            "; ".join(
                [f'{key} "{value}"' for key, value in self.attributes.items()]
            )
            + ";"
        )
        source = self.attributes.get("source", "unknown")
        return (
            f"{self.chrom}\t{source}\t{self.region_type}\t"
            f"{self.start}\t{self.end}\t.\t{self.strand}\t.\t{attr_str}"
        )

    def __str__(self) -> str:
        """Return a human-readable string representation.

        Returns:
            A string describing the region type, coordinates, strand,
            and attributes.
        """
        return (
            f"Region({self.region_type}, "
            f"{self.chrom}:{self.start}-{self.end}, "
            f"{self.strand}, attributes={self.attributes})"
        )

    def __repr__(self) -> str:
        """Return the string representation of the Region.

        Returns:
            Same as __str__.
        """
        return self.__str__()


@dataclass(slots=True)
class Transcript:
    """A transcript composed of an ordered list of Region objects.

    Attributes:
        transcript_id: Unique identifier for this transcript.
        strand: Strand, either ``"+"`` or ``"-"``.
        regions: Ordered list of exon/intron regions belonging to this
            transcript.
        attributes: Key-value pairs from the GTF attributes field.
    """

    transcript_id: str
    strand: str
    regions: list[Region] = field(default_factory=list)
    attributes: dict[str, str] = field(default_factory=dict)

    def add_region(self, region: Region) -> None:
        """Append a Region to this transcript's region list.

        Args:
            region: The Region object to add.
        """
        self.regions.append(region)

    def to_gtf_format(self) -> str:
        """Convert the Transcript to a GTF format string.

        The start and end coordinates are derived from the minimum and
        maximum coordinates of all contained regions.

        Returns:
            A tab-delimited GTF line representing this transcript.
        """
        attr_str = (
            "; ".join(
                [f'{key} "{value}"' for key, value in self.attributes.items()]
            )
            + ";"
        )
        start = min(region.start for region in self.regions)
        end = max(region.end for region in self.regions)
        chrom = self.regions[0].chrom
        source = self.attributes.get("source", "unknown")
        return (
            f"{chrom}\t{source}\ttranscript\t"
            f"{start}\t{end}\t.\t{self.strand}\t.\t{attr_str}"
        )

    def __repr__(self) -> str:
        """Return a string representation of the Transcript.

        Returns:
            A string with the transcript ID and its regions.
        """
        return f"Transcript: {self.transcript_id}, Regions: {self.regions}"


@dataclass(slots=True)
class Gene:
    """A gene composed of one or more Transcript objects.

    Attributes:
        gene_id: Unique identifier for this gene.
        transcripts: Mapping from transcript ID to Transcript object.
        attributes: Key-value pairs from the GTF attributes field.
    """

    gene_id: str
    transcripts: dict[str, Transcript] = field(default_factory=dict)
    attributes: dict[str, str] = field(default_factory=dict)

    def add_transcript(self, transcript: Transcript) -> None:
        """Add a Transcript to this gene, keyed by transcript ID.

        Args:
            transcript: The Transcript object to add.
        """
        self.transcripts[transcript.transcript_id] = transcript

    def to_gtf_format(self) -> str:
        """Convert the Gene to a GTF format string.

        The start and end coordinates are derived from the minimum and
        maximum coordinates across all regions of all transcripts.
        Strand and chromosome are taken from the first transcript.

        Returns:
            A tab-delimited GTF line representing this gene.
        """
        attr_str = (
            "; ".join(
                [f'{key} "{value}"' for key, value in self.attributes.items()]
            )
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
        first_transcript = next(iter(self.transcripts.values()))
        strand = first_transcript.strand
        chrom = first_transcript.regions[0].chrom
        source = self.attributes.get("source", "unknown")
        return (
            f"{chrom}\t{source}\tgene\t"
            f"{start}\t{end}\t.\t{strand}\t.\t{attr_str}"
        )

    def __repr__(self) -> str:
        """Return a string representation of the Gene.

        Returns:
            A string with the gene ID and its transcripts.
        """
        return f"Gene: {self.gene_id}, Transcripts: {self.transcripts}"
