#!/usr/bin/env python3

import os
import pandas as pd


def binary_table_threshold_with_percentage(binary_table, output_folder, threshold_percentage):

    df = pd.read_csv(binary_table, sep="\t")

    print(f"Number of mutations in the table: {len(df.columns[1:])}")

    cols_to_be_dropped = []

    amount_of_strains = int(df.shape[0])

    threshold_value = (threshold_percentage / 100) * amount_of_strains

    for col in df.columns[1:]:
        if int(df[col].sum()) <= int(threshold_value):
            cols_to_be_dropped.append(col)

    print(f"Number of mutations to be dropped: {len(cols_to_be_dropped)}")

    dropped_df = df.drop(cols_to_be_dropped, axis=1)

    print(
        f"Number of mutations in the table after dropping: {len(dropped_df.columns[1:])}")
    
    if f"with_gene_presence_absence" in binary_table:
        dropped_df.to_csv(os.path.join(
            output_folder, f"binary_mutation_table_with_gene_presence_absence_threshold_{threshold_percentage}_percent.tsv"), sep="\t", index=False)

    else:
        dropped_df.to_csv(os.path.join(
            output_folder, f"binary_mutation_table_threshold_{threshold_percentage}_percent.tsv"), sep="\t", index=False)
