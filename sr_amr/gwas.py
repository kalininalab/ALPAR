#!/usr/bin/env python3

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
import pathlib
import csv
from sr_amr.ml import decision_tree
import sys

csv.field_size_limit(sys.maxsize)

# Get the path of the script
PATH_OF_SCRIPT = pathlib.Path(__file__).parent.resolve()

# This function creates a genotype matrix from a binary mutation table.
# It reads the table into a pandas DataFrame, fills any missing values with "0",
# transposes the DataFrame, and then writes it to a file.

def pyseer_genotype_matrix_creator(binary_mutation_table, output_file):

    # Read the binary mutation table into a list of lists
    with open(binary_mutation_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        data = list(reader)

    # Fill any missing values with "0"
    for row in data:
        row += ["0"] * (len(data[0]) - len(row))

    # Transpose the list of lists
    data_transposed = list(map(list, zip(*data)))

    # Set the index name to "Gene"
    data_transposed[0][0] = "Gene"

    # Write the list of lists to a file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(data_transposed)

# This function creates a phenotype file for each antibiotic in the phenotype file.
# It reads the phenotype file into a pandas DataFrame, then for each column (antibiotic),
# it creates a new DataFrame with just that column, sets the index name to "samples",
# and writes it to a file with a name based on the antibiotic.
def pyseer_phenotype_file_creator(phenotype_file, output_file_directory):

    # Read the phenotype file into a DataFrame
    phenotype_df = pd.read_csv(phenotype_file, sep="\t", index_col=0)

    # For each antibiotic in the DataFrame
    for antibiotic in phenotype_df.columns:
        # Create a new DataFrame with just that column
        df_column = phenotype_df[[antibiotic]]

        # Set the index name to "samples"
        df_column.index.name = "samples"

        # Create a file name based on the column name
        output_file = f"{output_file_directory}/{antibiotic}.pheno"

        # Write the DataFrame to a file
        df_column.to_csv(output_file, sep="\t")

    for antibiotic in phenotype_df.columns:
        with open(f"{output_file_directory}/{antibiotic}.pheno", "r") as infile:
            lines = infile.readlines()
            with open(f"{output_file_directory}/{antibiotic}.pheno", "w") as ofile:
                ofile.write(lines[0])
                for line in lines[1:]:
                    if line.split("\t")[1].strip() == "1" or line.split("\t")[1].strip() == "0":
                        ofile.write(line)

def pyseer_individual_genotype_creator(pyseer_genotype_matrix, pyseer_phenotype_matrix, output_folder):

    # Read the genotype matrix into a dictionary
    with open(pyseer_genotype_matrix, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        genotype_header = next(reader)
        genotype_dict = {rows[0]: rows[1:] for rows in reader}

    # Read the phenotype matrix into a dictionary
    with open(pyseer_phenotype_matrix, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        phenotype_header = next(reader)
        phenotype_dict = {rows[0]: rows[1:] for rows in reader}

    strains_to_be_dropped = []

    with open(pyseer_phenotype_matrix) as infile:
        lines = infile.readlines()
        for line in lines[1:]:
            splitted = line.split("\t")
            if int(splitted[1].strip()) == 2:
                strains_to_be_dropped.append(splitted[0].strip())

    # Drop the strains from the genotype and phenotype data
    for strain in strains_to_be_dropped:
        del genotype_dict[strain]
        del phenotype_dict[strain]

    os.mkdir(os.path.join(output_folder, "genotype_files_individual"))
    os.mkdir(os.path.join(output_folder, "phenotype_files_individual"))

    # Write the genotype and phenotype data back to files
    with open(f"{output_folder}/genotype_files_individual/{os.path.basename(pyseer_phenotype_matrix)[:-6]}.tsv", 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(genotype_header)
        for key, value in genotype_dict.items():
            writer.writerow([key] + value)

    with open(f"{output_folder}/phenotype_files_individual/{os.path.basename(pyseer_phenotype_matrix)[:-6]}.pheno", 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(phenotype_header)
        for key, value in phenotype_dict.items():
            writer.writerow([key] + value)

# This function creates a similarity matrix from a phylogenetic tree.
# It runs a script with the phylogenetic tree as an argument and writes the output to a file.
def pyseer_similarity_matrix_creator(phylogenetic_tree, output_file):

    # Define the command to run the script
    script_command = f"python {PATH_OF_SCRIPT}/phylogeny_distance.py --lmm {phylogenetic_tree} > {output_file}"

    # Run the command
    os.system(script_command)

# This function runs the pyseer tool for each phenotype in the phenotype file path.
# It constructs a command to run pyseer with the appropriate arguments for each phenotype,
# checks if the output directory exists and creates it if not, then runs the command.
def pyseer_runner(genotype_file_path, phenotype_file_path, similarity_matrix, output_file_directory, threads):

    # Get a list of phenotypes
    phenotypes = os.listdir(f"{phenotype_file_path}")

    # For each phenotype
    for phenotype in phenotypes:

        # Construct the command to run pyseer
        script_command = f"pyseer --lmm --phenotypes {phenotype_file_path}/{phenotype} --pres {genotype_file_path} --similarity {similarity_matrix} --cpu {threads} > {output_file_directory}/{phenotype}.tsv"

        # If the output directory doesn't exist, create it
        if not os.path.exists(f"{output_file_directory}"):
            os.mkdir(f"{output_file_directory}")

        # Run the command
        os.system(script_command)

# This function processes the output of pyseer.
# It gets a list of result files, then for each file, it sorts the file based on the fourth column,
# writes the sorted file to a new file, then reads the sorted file and writes a cleaned version
# to another new file, removing any lines where the last column is "bad-chisq".
def pyseer_post_processor(pyseer_output_folder, output_folder):

    # Get a list of result files
    gwas_results_files = os.listdir(pyseer_output_folder)

    # For each result file
    for gwas_result_file in gwas_results_files:
        # Construct the command to sort the file based on the fourth column
        script_command = f"sort -g -k4,4 {os.path.join(pyseer_output_folder, gwas_result_file)}" + " > " + os.path.join(output_folder,'sorted', f"{gwas_result_file[:-4]}_sorted.tsv")

        # Run the command
        os.system(script_command)

        # Open a new file to write the cleaned version
        with open(os.path.join(output_folder, 'sorted_cleaned', f"{gwas_result_file[:-4]}_sorted.tsv"), "w") as ofile:
            # Open the sorted file to read
            with open(os.path.join(output_folder,'sorted', f"{gwas_result_file[:-4]}_sorted.tsv"), "r") as infile:
                # Read all lines
                lines = infile.readlines()
                # For each line
                for line in lines:
                    # Split the line into columns
                    splitted = line.split("\t")
                    # If the last column is not "bad-chisq", write the line to the cleaned file
                    if not splitted[-1].strip() == "bad-chisq":
                        ofile.write(f"{line.strip()}\n")


def pyseer_plot_file_creator(input_file, output_file):

    with open(input_file) as infile:
        lines = infile.readlines()

    with open(output_file, "w") as ofile:
        ofile.write("#CHR\tSNP\tBP\tminLOG10(P)\tlog10(p)\tr^2\n")
        for line in lines[1:]:
            splitted = line.split("\t")

            try:
                mut_position = int(
                    splitted[0].strip().split(",")[0].strip("'"))

            except:
                mut_position = int(0)

            lrt_p_val = float(splitted[3].strip())

            log_val = -math.log(lrt_p_val)

            ofile.write(
                f"26\t{splitted[0].strip()}\t{mut_position}\t{log_val}\t{log_val}\t0\n")


def pyseer_gwas_graph_creator(pyseer_output_folder, output_folder):

    gwas_sorted_files = os.listdir(
        os.path.join(pyseer_output_folder, "sorted_cleaned"))

    sorted_cleaned_files_path = os.path.join(pyseer_output_folder, "sorted_cleaned")

    raw_gwas_output_path = os.path.join(pyseer_output_folder, "gwas_results")

    for gwas_sorted_file in gwas_sorted_files:

        pyseer_plot_file_creator(os.path.join(sorted_cleaned_files_path, gwas_sorted_file), os.path.join(
            output_folder, f"{gwas_sorted_file[:-4]}.plot"))

        with open(os.path.join(raw_gwas_output_path, f"{gwas_sorted_file.split('_')[0]}.tsv")) as raw_gwas_file:
            lines = raw_gwas_file.readlines()
            threshold_denominator = len(lines)-1

        df = pd.read_csv(
            f"{os.path.join(output_folder, gwas_sorted_file[:-4])}.plot", sep='\t')

        bonferini_adjusted_threshold = 0.05 / threshold_denominator

        threshold = -(math.log(bonferini_adjusted_threshold))

        grid = sns.relplot(data=df, x='BP', y='log10(p)',
                           hue='log10(p)', palette='RdYlGn_r', aspect=1)

        grid.ax.set_xlabel("Position")
        grid.ax.set_ylabel("-log10(p-value)")
        grid.ax.axhline(threshold, linestyle='--', linewidth='1')
        plt.savefig(
            f"{os.path.join(output_folder, gwas_sorted_file[:-4])}.jpg", dpi=1200)

        os.remove(os.path.join(output_folder, f"{gwas_sorted_file[:-4]}.plot"))

def decision_tree_input_creator(binary_table, phenotype_file_path, pyseer_output_folder, output_folder):
    
    gwas_sorted_files = os.listdir(
        os.path.join(pyseer_output_folder, "sorted_cleaned"))

    sorted_cleaned_files_path = os.path.join(pyseer_output_folder, "sorted_cleaned")

    raw_gwas_output_path = os.path.join(pyseer_output_folder, "gwas_results")

    for gwas_sorted_file in gwas_sorted_files:

        pyseer_plot_file_creator(os.path.join(sorted_cleaned_files_path, gwas_sorted_file), os.path.join(
            output_folder, f"{gwas_sorted_file[:-4]}.plot"))

        with open(os.path.join(raw_gwas_output_path, f"{gwas_sorted_file.split('_')[0]}.tsv")) as raw_gwas_file:
            lines = raw_gwas_file.readlines()
            threshold_denominator = len(lines)-1

        bonferini_adjusted_threshold = 0.05 / threshold_denominator

        threshold = -(math.log(bonferini_adjusted_threshold))

        top_gwas_results = []
        
        with open(f"{os.path.join(output_folder, gwas_sorted_file[:-4])}.plot", "r") as plot_file:
            lines = plot_file.readlines()
            for line in lines[1:]:
                splitted = line.split("\t")
                if float(splitted[3].strip()) > threshold:
                    top_gwas_results.append(splitted[1].strip())

        # We need at least 2 features for decision trees
        if len(top_gwas_results) < 2:
            print(f"No significant SNPs found for {gwas_sorted_file.split('.')[0]}")
            print("Decision tree will not be generated for this antibiotic")
            os.remove(os.path.join(output_folder, f"{gwas_sorted_file[:-4]}.plot"))
           
        else:
            # Read the binary table into a dictionary
            with open(binary_table, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                headers = next(reader)
                header_indices = {name: index for index, name in enumerate(headers[1:], start=1)}
                binary_table_dict = {}
                for row in reader:
                    strain = row[0]
                    mutations = row[1:]
                    binary_table_dict[strain] = {mutation_name: mutations[header_indices[mutation_name]-1] for mutation_name in headers[1:]}
            
            cols_to_be_dropped = []

            for col in headers[1:]:
                if col not in top_gwas_results:
                    cols_to_be_dropped.append(col)
            
            for col in cols_to_be_dropped:
                for strain in binary_table_dict.keys():
                    del binary_table_dict[strain][col]

            headers = [header for header in headers if header not in cols_to_be_dropped]

            with open(os.path.join(output_folder, 'gwas_top_results.tsv'), 'w') as ofile:
                headers = ['Strain'] + list(next(iter(binary_table_dict.values())).keys())
                ofile.write('\t'.join(headers) + '\n')
                
                for strain, mutations in binary_table_dict.items():
                    row = [strain] + [mutations[mutation] for mutation in headers[1:]]
                    ofile.write('\t'.join(row) + '\n')

            os.makedirs(os.path.join(output_folder, 'decision_tree'), exist_ok=True)
            decision_tree(os.path.join(output_folder, 'gwas_top_results.tsv'), phenotype_file_path, gwas_sorted_file.split('.')[0], 42, 0.2, os.path.join(output_folder, 'decision_tree'))