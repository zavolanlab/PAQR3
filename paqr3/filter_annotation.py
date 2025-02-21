from gtfparse import read_gtf  # type: ignore
from pybedtools import BedTool  # type: ignore
from datetime import datetime


def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


def filter_annotation_by_gene_type(gtf_file, gene_types=["protein_coding"]):
    """Filters GTF by gene type, handling case variations for lncRNA."""

    log_message("Reading the GTF file...")
    gtf_data = read_gtf(gtf_file, result_type="pandas")

    log_message("Filtering genes by type...")

    # Standardize lncRNA casing for filtering
    gene_types_lower = [gt.lower() for gt in gene_types]
    gtf_data["gene_type_lower"] = gtf_data["gene_type"].str.lower()

    valid_genes = gtf_data[
        (gtf_data["feature"] == "gene")
        & (gtf_data["gene_type_lower"].isin(gene_types_lower))
    ]["gene_id"]

    log_message(f"Found {len(valid_genes)} genes matching the criteria.")

    filtered_gtf = gtf_data[gtf_data["gene_id"].isin(valid_genes)]
    filtered_gtf = filtered_gtf.drop(
        columns="gene_type_lower"
    )  # Drop temporary column

    return filtered_gtf


def resolve_gene_overlaps(filtered_gtf, strandedness=True):
    """
    Resolves gene overlaps based on strandedness,
    updating exon coordinates.
    """

    log_message("Resolving overlaps...")
    gene_lengths = {}
    gene_coordinates = {}  # Store original gene coordinates
    for _, row in filtered_gtf[filtered_gtf["feature"] == "gene"].iterrows():
        gene_id = row["gene_id"]
        gene_lengths[gene_id] = row["end"] - row["start"]
        gene_coordinates[gene_id] = (row["start"], row["end"])

    # Create BedTool from GENE coordinates ONLY
    gene_bed = BedTool.from_dataframe(
        filtered_gtf[filtered_gtf["feature"] == "gene"][
            ["seqname", "start", "end", "gene_id", "score", "strand"]
        ]
    )

    log_message("Checking for overlaps...")
    if strandedness:
        overlaps = gene_bed.intersect(gene_bed, wo=True, s=True)
    else:
        overlaps = gene_bed.intersect(gene_bed, wo=True)

    log_message("Resolving overlaps...")
    genes_to_remove = set()
    regions_to_update = {}

    overlap_groups = {}
    for line in overlaps:
        fields = line.fields
        gene1 = fields[3]
        gene2 = fields[9]
        start1, end1 = int(fields[1]), int(fields[2])
        start2, end2 = int(fields[7]), int(fields[8])
        chrom = fields[0]
        strand = fields[5] if strandedness else "+"

        group_key = (chrom, strand) if strandedness else (chrom,)
        if group_key not in overlap_groups:
            overlap_groups[group_key] = []
        overlap_groups[group_key].append(
            (gene1, gene2, start1, end1, start2, end2)
        )

    for group_key, overlaps_in_group in overlap_groups.items():
        genes_in_group = set()
        for gene1, gene2, _, _, _, _ in overlaps_in_group:
            genes_in_group.add(gene1)
            genes_in_group.add(gene2)

        genes_to_remove_in_group = set()
        regions_to_update_in_group = {}

        overlaps_in_group.sort(
            key=lambda x: x[5] - x[4] if x[5] > x[4] else x[3] - x[2],
            reverse=True,
        )

        updated_genes = set()

        for gene1, gene2, start1, end1, start2, end2 in overlaps_in_group:
            if gene1 == gene2:
                continue

            if (
                gene1 in genes_to_remove_in_group
                or gene2 in genes_to_remove_in_group
            ):
                continue

            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            overlap_len = overlap_end - overlap_start

            if (
                start1 <= start2 and end1 >= end2
            ):  # Gene1 completely overlaps Gene2
                genes_to_remove_in_group.add(
                    gene2
                    if gene_lengths[gene1] >= gene_lengths[gene2]
                    else gene1
                )

            elif (
                start2 <= start1 and end2 >= end1
            ):  # Gene2 completely overlaps Gene1
                genes_to_remove_in_group.add(
                    gene1
                    if gene_lengths[gene2] >= gene_lengths[gene1]
                    else gene2
                )

            elif overlap_len > 0:  # Partial overlap
                current_start1 = (
                    regions_to_update_in_group.get(gene1, (start1, end1))[0]
                    if gene1 in updated_genes
                    else start1
                )
                current_end1 = (
                    regions_to_update_in_group.get(gene1, (start1, end1))[1]
                    if gene1 in updated_genes
                    else end1
                )
                current_start2 = (
                    regions_to_update_in_group.get(gene2, (start2, end2))[0]
                    if gene2 in updated_genes
                    else start2
                )
                current_end2 = (
                    regions_to_update_in_group.get(gene2, (start2, end2))[1]
                    if gene2 in updated_genes
                    else end2
                )

                new_start1 = (
                    current_start1
                    if current_start1 < overlap_start
                    else overlap_end + 1
                )
                new_end1 = (
                    current_end1
                    if current_end1 > overlap_end
                    else overlap_start - 1
                )
                new_start2 = (
                    current_start2
                    if current_start2 < overlap_start
                    else overlap_end + 1
                )
                new_end2 = (
                    current_end2
                    if current_end2 > overlap_end
                    else overlap_start - 1
                )

                if new_start1 <= new_end1:
                    regions_to_update_in_group[gene1] = (new_start1, new_end1)
                    updated_genes.add(gene1)
                if new_start2 <= new_end2:
                    regions_to_update_in_group[gene2] = (new_start2, new_end2)
                    updated_genes.add(gene2)

        genes_to_remove.update(genes_to_remove_in_group)
        regions_to_update.update(regions_to_update_in_group)

    log_message(
        f"Removing {len(genes_to_remove)} completely overlapping genes."
    )
    filtered_gtf = filtered_gtf[~filtered_gtf["gene_id"].isin(genes_to_remove)]

    # Update exon coordinates based on resolved gene overlaps
    log_message("Updating exon coordinates...")
    updated_gtf = filtered_gtf.copy()

    for gene_id, (new_start, new_end) in regions_to_update.items():
        updated_gtf.loc[
            (updated_gtf["gene_id"] == gene_id)
            & (updated_gtf["feature"] != "gene"),
            "start",
        ] = updated_gtf.loc[
            (updated_gtf["gene_id"] == gene_id)
            & (updated_gtf["feature"] != "gene"),
            "start",
        ].apply(
            lambda x: max(x, new_start)
        )
        updated_gtf.loc[
            (updated_gtf["gene_id"] == gene_id)
            & (updated_gtf["feature"] != "gene"),
            "end",
        ] = updated_gtf.loc[
            (updated_gtf["gene_id"] == gene_id)
            & (updated_gtf["feature"] != "gene"),
            "end",
        ].apply(
            lambda x: min(x, new_end)
        )

    return updated_gtf
