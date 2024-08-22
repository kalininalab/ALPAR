#!/usr/bin/env python3

import csv
import os
import sys

import warnings

warnings.filterwarnings("ignore")

csv.field_size_limit(sys.maxsize)

def binary_table_threshold_with_percentage(binary_table, output_folder, threshold_percentage):

    with open(binary_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        header_indices = {name: index for index, name in enumerate(headers[1:], start=1)}
        binary_table_dict = {}
        for row in reader:
            strain = row[0]
            mutations = row[1:]
            binary_table_dict[strain] = {mutation_name: mutations[header_indices[mutation_name]-1] for mutation_name in headers[1:]}
            
    print(f"Number of mutations in the table: {len(headers[1:])}")

    cols_to_be_dropped = []

    amount_of_strains = len(binary_table_dict)

    threshold_value = (threshold_percentage / 100) * amount_of_strains

    threshold_value = int(threshold_value)

    if threshold_value < 1:
        threshold_value = 1

    print(f"Threshold value: {threshold_value}")

    cols_to_be_dropped = []

    for col in headers[1:]:
        count_ones = sum(binary_table_dict[strain][col] == '1' for strain in binary_table_dict)
        if count_ones <= threshold_value:
            cols_to_be_dropped.append(col)

    print(f"Number of mutations to be dropped: {len(cols_to_be_dropped)}")

    cols_to_be_dropped_set = set(cols_to_be_dropped)

    for col in cols_to_be_dropped:
        for strain in binary_table_dict.keys():
            del binary_table_dict[strain][col]

    headers = [header for header in headers if header not in cols_to_be_dropped_set]

    print(f"Number of mutations in the table after dropping: {len(headers) - 1}")

    if "with_gene_presence_absence" in binary_table:
        output_file = os.path.join(output_folder, f"binary_mutation_table_with_gene_presence_absence_threshold_{threshold_percentage}_percent.tsv")
    else:
        output_file = os.path.join(output_folder, f"binary_mutation_table_threshold_{threshold_percentage}_percent.tsv")

    with open(f"{output_file}", 'w') as file:
        headers = ['Strain'] + list(next(iter(binary_table_dict.values())).keys())
        file.write('\t'.join(headers) + '\n')
        
        for strain, mutations in binary_table_dict.items():
            row = [strain] + [mutations[mutation] for mutation in headers[1:]]
            file.write('\t'.join(row) + '\n')

    return output_file