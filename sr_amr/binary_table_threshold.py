#!/usr/bin/env python3

import csv
import os
import sys

csv.field_size_limit(sys.maxsize)

def binary_table_threshold_with_percentage(binary_table, output_folder, threshold_percentage):

    with open(binary_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        binary_table_dict = {rows[0]: rows[1:] for rows in reader}

    print(f"Number of mutations in the table: {len(headers[1:])}")

    cols_to_be_dropped = []

    amount_of_strains = len(binary_table_dict)

    threshold_value = (threshold_percentage / 100) * amount_of_strains

    col_indices = {col: idx for idx, col in enumerate(headers)}

    cols_to_be_dropped = []

    for col in headers[1:]:
        col_index = col_indices[col] - 1
        count_ones = sum(binary_table_dict[strain][col_index] == '1' for strain in binary_table_dict)
        if count_ones <= threshold_value:
            cols_to_be_dropped.append(col)

    print(f"Number of mutations to be dropped: {len(cols_to_be_dropped)}")

    cols_to_be_dropped_set = set(cols_to_be_dropped)

    indices_to_keep = [index for index, header in enumerate(headers) if header not in cols_to_be_dropped_set]
    
    adjusted_indices_to_keep = [0] + [index - 1 for index in indices_to_keep[1:]]

    for strain in binary_table_dict:
        binary_table_dict[strain] = [binary_table_dict[strain][index] for index in adjusted_indices_to_keep]

    print(f"Number of mutations in the table after dropping: {len(headers[1:]) - len(cols_to_be_dropped)}")

    if "with_gene_presence_absence" in binary_table:
        output_file = os.path.join(output_folder, f"binary_mutation_table_with_gene_presence_absence_threshold_{threshold_percentage}_percent.tsv")
    else:
        output_file = os.path.join(output_folder, f"binary_mutation_table_threshold_{threshold_percentage}_percent_2.tsv")

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        for key, value in binary_table_dict.items():
            writer.writerow([key] + value)

    return output_file