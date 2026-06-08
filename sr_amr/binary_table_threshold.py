#!/usr/bin/env python3

import csv
import os
import sys

import warnings

warnings.filterwarnings("ignore")

csv.field_size_limit(sys.maxsize)

def binary_table_threshold_with_percentage(binary_table, output_folder, threshold_percentage):
    col_counts = []
    num_strains = 0
    
    with open(binary_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        mutation_names = headers[1:]
        col_counts = [0] * len(mutation_names)
        
        for row in reader:
            num_strains += 1
            for i, val in enumerate(row[1:]):
                if val == '1' or val == '1.0':
                    col_counts[i] += 1

    threshold_value = max(1, int((threshold_percentage / 100) * num_strains))
    
    keep_indices = [i for i, count in enumerate(col_counts) if count > threshold_value]
    
    print(f"Total strains: {num_strains}")
    print(f"Threshold: {threshold_value}")
    print(f"Mutations dropped: {len(mutation_names) - len(keep_indices)}")

    file_suffix = "with_gene_presence_absence" if "with_gene_presence_absence" in binary_table else ""
    output_filename = f"binary_mutation_table_{file_suffix}_threshold_{threshold_percentage}_percent.tsv"
    output_file = os.path.join(output_folder, output_filename.replace("__", "_"))

    with open(binary_table, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        original_headers = next(reader)
        new_headers = [original_headers[0]] + [mutation_names[i] for i in keep_indices]
        writer.writerow(new_headers)

        for row in reader:
            new_row = [row[0]] + [row[i+1] for i in keep_indices]
            writer.writerow(new_row)

    return output_file