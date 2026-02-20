#!/usr/bin/env python3

import csv
import os
import sys

import warnings

warnings.filterwarnings("ignore")

csv.field_size_limit(sys.maxsize)

def binary_table_threshold_with_percentage(binary_table, output_folder, threshold_percentage):
    # 1. First Pass: Count '1's for each column
    col_counts = []
    num_strains = 0
    
    with open(binary_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        mutation_names = headers[1:]
        # Initialize counts with zeros for each mutation
        col_counts = [0] * len(mutation_names)
        
        for row in reader:
            num_strains += 1
            # Using enumerate is faster than dictionary lookups here
            for i, val in enumerate(row[1:]):
                if val == '1' or val == '1.0':
                    col_counts[i] += 1

    # 2. Calculate Threshold
    threshold_value = max(1, int((threshold_percentage / 100) * num_strains))
    
    # Identify indices of columns to KEEP
    keep_indices = [i for i, count in enumerate(col_counts) if count > threshold_value]
    
    print(f"Total strains: {num_strains}")
    print(f"Threshold: {threshold_value}")
    print(f"Mutations dropped: {len(mutation_names) - len(keep_indices)}")

    # 3. Second Pass: Filter and Write
    file_suffix = "with_gene_presence_absence" if "with_gene_presence_absence" in binary_table else ""
    output_filename = f"binary_mutation_table_{file_suffix}_threshold_{threshold_percentage}_percent.tsv"
    output_file = os.path.join(output_folder, output_filename.replace("__", "_"))

    with open(binary_table, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        # Write filtered headers
        original_headers = next(reader)
        # Strain column + only the mutations that passed the threshold
        new_headers = [original_headers[0]] + [mutation_names[i] for i in keep_indices]
        writer.writerow(new_headers)
        
        # Write filtered rows
        for row in reader:
            # Construct new row by grabbing only the kept indices (+1 to offset Strain)
            new_row = [row[0]] + [row[i+1] for i in keep_indices]
            writer.writerow(new_row)

    return output_file