import pandas as pd
from gtfparse import read_gtf
from pybedtools import BedTool
from datetime import datetime

# Utility function for timestamped messages
def log_message(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def filter_gtf_by_gene_type(gtf_file, gene_types=["protein_coding", "lncrna"]):
    """
    Filter the GTF file for genes and associated entries of specified types.

    Args:
        gtf_file (str): Path to the input GTF file.
        gene_types (list): List of gene types to retain (default: ["protein_coding", "lncrna"]).

    Returns:
        pd.DataFrame: Filtered GTF data.
    """
    log_message("Reading the GTF file...")
    gtf_data = read_gtf(gtf_file, result_type="pandas")

    log_message("Filtering genes by type...")
    valid_genes = gtf_data[
        (gtf_data["feature"] == "gene") & (gtf_data["gene_biotype"].isin(gene_types))
    ]["gene_id"]

    log_message(f"Found {len(valid_genes)} genes matching the criteria.")

    filtered_gtf = gtf_data[gtf_data["gene_id"].isin(valid_genes)]
    return filtered_gtf

def resolve_strand_overlaps(filtered_gtf):
    """
    Resolve strand overlaps by modifying or removing regions based on overlap conditions.

    Args:
        filtered_gtf (pd.DataFrame): Filtered GTF data.

    Returns:
        pd.DataFrame: Updated GTF data with resolved overlaps.
    """
    log_message("Converting GTF data to BedTool format...")
    gtf_bed = BedTool.from_dataframe(
        filtered_gtf[["seqname", "start", "end", "gene_id", "score", "strand"]]
    )

    log_message("Checking for overlaps...")
    overlaps = gtf_bed.intersect(gtf_bed, wo=True, s=True)

    log_message("Resolving overlaps...")
    to_remove = set()
    to_update = {}

    for line in overlaps:
        fields = line.fields
        gene1 = fields[3]  # First gene ID
        gene2 = fields[9]  # Second gene ID
        start1, end1 = int(fields[1]), int(fields[2])
        start2, end2 = int(fields[7]), int(fields[8])
        overlap_length = int(fields[-1])

        if gene1 == gene2:
            continue

        if start1 <= start2 and end1 >= end2:  # Gene1 completely overlaps Gene2
            to_remove.add(gene2 if end1 - start1 > end2 - start2 else gene1)

        elif start2 <= start1 and end2 >= end1:  # Gene2 completely overlaps Gene1
            to_remove.add(gene1 if end2 - start2 > end1 - start1 else gene2)

        else:  # Partial overlap
            new_start1, new_end1 = max(start1, end2 + 1), min(end1, start2 - 1)
            new_start2, new_end2 = max(start2, end1 + 1), min(end2, start1 - 1)

            if new_start1 < new_end1:
                to_update[gene1] = (new_start1, new_end1)
            if new_start2 < new_end2:
                to_update[gene2] = (new_start2, new_end2)

    log_message(f"Removing {len(to_remove)} completely overlapping genes.")
    filtered_gtf = filtered_gtf[~filtered_gtf["gene_id"].isin(to_remove)]

    log_message(f"Updating {len(to_update)} genes with modified coordinates.")
    for gene_id, (new_start, new_end) in to_update.items():
        filtered_gtf.loc[filtered_gtf["gene_id"] == gene_id, "start"] = new_start
        filtered_gtf.loc[filtered_gtf["gene_id"] == gene_id, "end"] = new_end

    return filtered_gtf

def write_filtered_gtf(filtered_gtf, output_file):
    """
    Write the filtered GTF data to a new file.

    Args:
        filtered_gtf (pd.DataFrame): Filtered GTF data.
        output_file (str): Path to the output GTF file.
    """
    filtered_gtf.to_csv(output_file, sep="\t", header=False, index=False, quoting=None)
    log_message(f"Filtered GTF written to {output_file}")
