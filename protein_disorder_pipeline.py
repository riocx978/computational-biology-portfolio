#!/usr/bin/env python3
"""
Protein Disorder Annotation & Dataset Splitting Pipeline
Author: Rhea Charles | University of South Florida

Description:
    Processes UniProt and DisProt databases to extract intrinsically disordered
    protein (IDP) region annotations, then splits the resulting dataset into
    train/test sets using CD-HIT sequence clusters to prevent data leakage.

    This pipeline was developed as part of a deep learning project to predict
    intrinsically disordered protein linkers.

Workflow:
    1. Parse UniProt JSON — extract linker/disordered region annotations
    2. Parse DisProt JSON — encode disorder state per residue (D=1, S=2, unknown=0)
    3. Parse CD-HIT clusters — separate empty vs non-empty clusters
    4. Split non-empty clusters 70/30 (train/test) to avoid sequence similarity leakage

Input:
    uniprot.json                                    — UniProt API export
    DisProt release_2023_06 with_ambiguous_evidences.json — DisProt database release
    Disprot_output.fasta.clstr                      — CD-HIT clustering output

Output:
    output_with_sequences.json                      — UniProt annotations with sequences
    DisProt_disorder_annotation_ambiguous.txt       — per-residue disorder encoding
    Disprot_output.fasta.emptyclstrs                — singleton clusters
    Disprot_output.fasta.nonemptyclstrs             — multi-member clusters
    Disprotclusters_70_percent.clstr                — training split
    Disprotclusters_30_percent.clstr                — test split
"""

import os
import json
import random
import numpy as np


# -----------------------------------------------------------------------------
# Configuration — update these paths for your environment
# -----------------------------------------------------------------------------

UNIPROT_DIR   = "data/uniprot"
DISPROT_DIR   = "data/disprot"
SPLIT_DIR     = "data/disprot/split"

UNIPROT_FILE  = "uniprot.json"
DISPROT_FILE  = "DisProt release_2023_06 with_ambiguous_evidences.json"
CDHIT_FILE    = "Disprot_output.fasta.clstr"
SEQUENCES_FILE = "DisProt_sequences.fasta"

TRAIN_SPLIT   = 0.7   # Fraction of clusters for training set
RANDOM_SEED   = 42    # For reproducibility

random.seed(RANDOM_SEED)


# =============================================================================
# Step 1. Parse UniProt — extract linker and disordered region annotations
# =============================================================================

def parse_uniprot(filepath):
    """
    Extracts protein features annotated as linker, binding, or disordered
    regions from a UniProt JSON export.

    Returns a list of dicts with protein metadata, region coordinates,
    and full amino acid sequence.
    """
    with open(filepath, "r") as f:
        data = json.load(f)

    records = []
    target_terms = {"linker", "binding", "Disordered"}

    for entry in data["results"]:
        protein_regions = {}

        for feature in entry.get("features", []):
            description = feature.get("description", "")

            if not any(term in description for term in target_terms):
                continue

            feature_type = feature.get("type")
            start = feature["location"]["start"]["value"]
            end   = feature["location"]["end"]["value"]

            if feature_type not in protein_regions:
                protein_regions[feature_type] = []
            protein_regions[feature_type].append({"start": start, "end": end})

        protein_desc = entry.get("proteinDescription", {})

        record = {
            "uniProtkbId":        entry.get("uniProtkbId", ""),
            "uniprotID":          entry.get("primaryAccession", ""),
            "scientificName":     entry.get("organism", {}).get("scientificName", ""),
            "taxonId":            entry.get("organism", {}).get("taxonId", ""),
            "regions":            protein_regions,
            "recommendedName":    protein_desc.get("recommendedName", {})
                                              .get("fullName", {})
                                              .get("value", ""),
            "alternativeNames":   [n.get("fullName", {}).get("value", "")
                                   for n in protein_desc.get("alternativeNames", [])],
            "geneNames":          [g.get("geneName", {}).get("value", "")
                                   for g in entry.get("genes", [])],
            "sequence":           entry.get("sequence", {}).get("value", "")
        }
        records.append(record)

    return records


os.chdir(UNIPROT_DIR)
uniprot_records = parse_uniprot(UNIPROT_FILE)

with open("output_with_sequences.json", "w") as f:
    json.dump(uniprot_records, f, indent=4)

print(f"UniProt: {len(uniprot_records)} records saved to output_with_sequences.json")


# =============================================================================
# Step 2. Parse DisProt — encode disorder state per residue
# =============================================================================
# Each residue is encoded as:
#   0 = unknown / unannotated
#   1 = disordered (D)
#   2 = structured (S)
#
# Output format (FASTA-like):
#   > uniprotID_disprotID
#   <amino acid sequence>
#   <disorder encoding string>

def parse_disprot(filepath, encoding="utf-8"):
    """
    Reads a DisProt JSON release and produces per-residue disorder annotations.
    Returns a list of formatted annotation strings.
    """
    with open(filepath, "rt", encoding=encoding) as f:
        data = json.load(f)

    annotations = []

    for entry in data["data"]:
        uniprot_acc = entry["acc"]
        disprot_id  = entry["disprot_id"]
        sequence    = entry["sequence"]

        # Initialize all residues as unknown (0)
        disorder_encoding = np.zeros(len(sequence), dtype=np.uint8)

        structural_states = entry.get("disprot_consensus", {}) \
                                 .get("Structural state", [])

        for region in structural_states:
            start = region["start"] - 1  # Convert to 0-indexed
            end   = region["end"]

            if region["type"] == "D":
                disorder_encoding[start:end] = 1
            elif region["type"] == "S":
                disorder_encoding[start:end] = 2

        annotation = (
            f">{uniprot_acc}_{disprot_id}\n"
            f"{sequence.strip()}\n"
            f"{''.join(map(str, disorder_encoding))}"
        )
        annotations.append(annotation)

    return annotations


os.chdir(DISPROT_DIR)
disprot_annotations = parse_disprot(DISPROT_FILE)

with open("DisProt_disorder_annotation_ambiguous.txt", "wt") as f:
    f.write("\n".join(disprot_annotations))

print(f"DisProt: {len(disprot_annotations)} sequences annotated")


# =============================================================================
# Step 3. Parse CD-HIT clusters — separate singletons from multi-member clusters
# =============================================================================
# Singletons (only one sequence in a cluster) are separated from multi-member
# clusters. Only multi-member clusters are used for train/test splitting to
# ensure meaningful sequence diversity in both sets.

def parse_cdhit_clusters(clstr_filepath, out_empty, out_nonempty):
    """
    Splits a CD-HIT .clstr file into singleton and multi-member cluster files.

    Args:
        clstr_filepath: path to .clstr file from CD-HIT
        out_empty:      output path for singleton clusters
        out_nonempty:   output path for multi-member clusters
    """
    with open(out_empty, "w") as f_single, \
         open(out_nonempty, "w") as f_multi, \
         open(clstr_filepath, "r") as f_in:

        current_name  = None
        current_lines = []
        member_count  = 0

        for line in f_in:
            if line.startswith(">"):
                # Flush previous cluster
                if current_name:
                    target = f_single if member_count == 1 else f_multi
                    target.writelines(current_lines)

                current_name  = line.strip()
                current_lines = [line]
                member_count  = 0
            else:
                current_lines.append(line)
                member_count += 1

        # Flush final cluster
        if current_name:
            target = f_single if member_count == 1 else f_multi
            target.writelines(current_lines)


parse_cdhit_clusters(
    CDHIT_FILE,
    out_empty    = "Disprot_output.fasta.emptyclstrs",
    out_nonempty = "Disprot_output.fasta.nonemptyclstrs"
)
print("CD-HIT clusters parsed into singleton and multi-member files")


# =============================================================================
# Step 4. Split clusters 70/30 — train/test without sequence similarity leakage
# =============================================================================
# Splitting at the cluster level (not sequence level) ensures that sequences
# with high similarity do not appear in both train and test sets — a critical
# step for unbiased model evaluation in protein ML tasks.

def split_clusters(clstr_filepath, train_path, test_path, train_fraction=0.7):
    """
    Randomly splits CD-HIT clusters into train and test sets.

    Splitting by cluster (not by sequence) prevents data leakage from
    highly similar sequences appearing in both sets.
    """
    with open(clstr_filepath, "r") as f:
        lines = f.read().splitlines()

    # Group lines into clusters
    clusters = []
    current  = []

    for line in lines:
        if line.startswith(">Cluster"):
            if current:
                clusters.append(current)
            current = [line]
        else:
            current.append(line)

    if current:
        clusters.append(current)

    random.shuffle(clusters)

    n_train = int(len(clusters) * train_fraction)
    train_clusters = clusters[:n_train]
    test_clusters  = clusters[n_train:]

    with open(train_path, "w") as f:
        for cluster in train_clusters:
            f.write("\n".join(cluster) + "\n")

    with open(test_path, "w") as f:
        for cluster in test_clusters:
            f.write("\n".join(cluster) + "\n")

    return len(clusters), n_train


os.chdir(SPLIT_DIR)

total, n_train = split_clusters(
    f"../{'/'.join(['Disprot_output.fasta.nonemptyclstrs'])}",
    train_path     = "Disprotclusters_70_percent.clstr",
    test_path      = "Disprotclusters_30_percent.clstr",
    train_fraction = TRAIN_SPLIT
)

print(f"\nDataset split complete:")
print(f"  Total clusters : {total}")
print(f"  Training (70%) : {n_train}")
print(f"  Test (30%)     : {total - n_train}")
