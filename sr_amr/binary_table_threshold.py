#!/usr/bin/env python3

import csv
import os

def binary_table_threshold_with_percentage(binary_table, output_folder, threshold_percentage):

    with open(binary_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        binary_table_dict = {rows[0]: rows[1:] for rows in reader}

    print(f"Number of mutations in the table: {len(headers[1:])}")

    cols_to_be_dropped = []

    amount_of_strains = len(binary_table_dict)

    threshold_value = (threshold_percentage / 100) * amount_of_strains

    for col in headers[1:]:
        if sum(int(binary_table_dict[strain][headers.index(col)-1] == '1') for strain in binary_table_dict) <= threshold_value:
            cols_to_be_dropped.append(col)

    print(f"Number of mutations to be dropped: {len(cols_to_be_dropped)}")

    for strain in binary_table_dict:
        for col in cols_to_be_dropped:
            del binary_table_dict[strain][headers.index(col)-1]

    print(f"Number of mutations in the table after dropping: {len(headers[1:]) - len(cols_to_be_dropped)}")

    if "with_gene_presence_absence" in binary_table:
        output_file = os.path.join(output_folder, f"binary_mutation_table_with_gene_presence_absence_threshold_{threshold_percentage}_percent.tsv")
    else:
        output_file = os.path.join(output_folder, f"binary_mutation_table_threshold_{threshold_percentage}_percent.tsv")

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        for key, value in binary_table_dict.items():
            writer.writerow([key] + value)

    return output_file