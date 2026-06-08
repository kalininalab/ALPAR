#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np

def calculate_median_length(fasta_files):
    lengths = []
    for fasta in fasta_files:
        length = 0
        with open(fasta, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    length += len(line.strip())
        lengths.append(length)
    return np.median(lengths), lengths

def check_genome_quality(fasta_path, median_length=None, length_threshold=0.1, max_contigs=500):
    """
    Explicit QC check for a genome.
    - length_threshold: fraction of median length allowed (default 0.1 for 90%-110%)
    - max_contigs: maximum allowed number of contigs
    """
    length = 0
    contigs = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contigs += 1
            else:
                length += len(line.strip())
    
    reasons = []
    if median_length is not None:
        lower_bound = median_length * (1 - length_threshold)
        upper_bound = median_length * (1 + length_threshold)
        if length < lower_bound or length > upper_bound:
            reasons.append(f"Length {length} is outside bounds [{lower_bound:.0f}, {upper_bound:.0f}]")
    
    if contigs > max_contigs:
        reasons.append(f"Contig count {contigs} exceeds maximum {max_contigs}")
        
    if reasons:
        return False, "; ".join(reasons)
    return True, "Pass"

def run_qc_pipeline(input_strains, output_folder, length_threshold=0.1, max_contigs=500):
    print("Running Quality Control on input genomes...")
    median_length, lengths = calculate_median_length(input_strains)
    
    qc_results = []
    passed_strains = []
    
    for strain in input_strains:
        passed, reason = check_genome_quality(strain, median_length, length_threshold, max_contigs)
        qc_results.append({
            'strain': os.path.basename(strain),
            'path': strain,
            'passed': passed,
            'reason': reason
        })
        if passed:
            passed_strains.append(strain)
            
    qc_df = pd.DataFrame(qc_results)
    qc_df.to_csv(os.path.join(output_folder, "qc_report.tsv"), sep="\t", index=False)
    
    print(f"QC finished: {len(passed_strains)}/{len(input_strains)} genomes passed.")
    return passed_strains
