#!/usr/bin/env python3

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
import pathlib

# Get the path of the script
PATH_OF_SCRIPT = pathlib.Path(__file__).parent.resolve()

def pyseer_genotype_matrix_creator(binary_mutation_table, output_file):

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t", index_col=0)

    genotype_df_transposed = genotype_df.transpose()

    genotype_df_transposed.index.name = "Gene"
    
    genotype_df_transposed.to_csv(output_file, sep="\t")


def pyseer_phenotype_file_creator(phenotype_file, output_file_directory):

    phenotype_df = pd.read_csv(phenotype_file, sep="\t", index_col=0)

    for antibiotic in phenotype_df.columns:
        df_column = phenotype_df[[antibiotic]]

        df_column.index.name = "samples"

        # Create a file name based on the column name
        output_file = f"{output_file_directory}/{antibiotic}.pheno"

        # Write the DataFrame to a file
        df_column.to_csv(output_file, sep="\t")


def pyseer_individual_genotype_creator(pyseer_genotype_matrix, pyseer_phenotype_matrix, output_folder):

    genotype_df = pd.read_csv(pyseer_genotype_matrix, sep="\t", index_col=0)

    phenotype_df = pd.read_csv(pyseer_phenotype_matrix, sep="\t", index_col=0)

    strains_to_be_dropped = []

    with open(pyseer_phenotype_matrix) as infile:
        lines = infile.readlines()
        for line in lines[1:]:
            splitted = line.split("\t")
            if int(splitted[1].strip()) == 2:
                strains_to_be_dropped.append(splitted[0].strip())
    
    genotype_df_dropped = genotype_df.drop(strains_to_be_dropped, axis=1)

    phenotype_df_dropped = phenotype_df.drop(strains_to_be_dropped, axis=0)

    os.mkdir(os.path.join(output_folder, "genotype_files_individual"))

    os.mkdir(os.path.join(output_folder, "phenotype_files_individual"))

    genotype_df_dropped.to_csv(f"{output_folder}/genotype_files_individual/{pyseer_phenotype_matrix.split('/')[-1][:-6]}.tsv", sep="\t")

    phenotype_df_dropped.to_csv(f"{output_folder}/phenotype_files_individual/{pyseer_phenotype_matrix.split('/')[-1][:-6]}.pheno", sep="\t")


def pyseer_similarity_matrix_creator(phylogenetic_tree, output_file):
    
    script_command = f"python {PATH_OF_SCRIPT}/phylogeny_distance.py --lmm {phylogenetic_tree} > {output_file}"

    os.system(script_command)


def pyseer_runner(genotype_file_path, phenotype_file_path, similarity_matrix, output_file_directory):

    phenotypes = os.listdir(f"{phenotype_file_path}")

    genotypes = os.listdir(f"{genotype_file_path}")

    for phenotype in phenotypes:
        script_command = f"pyseer --lmm --phenotypes {phenotype_file_path}/{phenotype} --pres {genotypes}/{phenotype[:-6]}.tsv --similarity {similarity_matrix} > {output_file_directory}/{phenotype}.tsv"
                
        if not os.path.exists(f"{output_file_directory}"):
            os.mkdir(f"{output_file_directory}")

        os.system(script_command)


def pyseer_post_processor(pyseer_output_folder, output_folder):

    gwas_results_files = os.listdir(pyseer_output_folder)

    for gwas_result_file in gwas_results_files:
        script_command = f"sort -g -k4,4 {pyseer_output_folder}/{gwas_result_file} > {output_folder}/sorted/{gwas_result_file[:-4]}_sorted.tsv"
    
        os.system(script_command)

        with open (f"{output_folder}/sorted_cleaned/{gwas_result_file[:-4]}_sorted.tsv", "w") as ofile:
            with open(f"{output_folder}/sorted/{gwas_result_file[:-4]}_sorted.tsv", "r") as infile:
                lines = infile.readlines()
                for line in lines:
                    splitted = line.split("\t")
                    if not splitted[-1].strip() == "bad-chisq":
                        ofile.write(f"{line.strip()}\n")


def pyseer_plot_file_creator(input_file, output_file):

    with open(input_file) as infile:
        lines = infile.readlines()
    
    with open(output_file, "w") as ofile:
        ofile.write("#CHR\tSNP\tBP\tminLOG10(P)\tlog10(p)\tr^2\n")
        for line in lines[1:]:
            splitted = line.split("\t")

            mut_position = int(splitted[0].strip().split(",")[0].strip("'"))

            lrt_p_val = float(splitted[3].strip())

            log_val = -math.log(lrt_p_val)

            ofile.write(f"26\t{splitted[0].strip()}\t{mut_position}\t{log_val}\t{log_val}\t0\n")


def pyseer_gwas_graph_creator(pyseer_output_folder, output_folder):

    gwas_sorted_files = os.listdir(os.path.join(pyseer_output_folder,"sorted"))

    sorted_files_path = os.path.join(pyseer_output_folder,"sorted")

    raw_gwas_output_path = os.path.join(pyseer_output_folder, "gwas_results")

    for gwas_sorted_file in gwas_sorted_files:

        pyseer_plot_file_creator(os.path.join(sorted_files_path, gwas_sorted_file), os.path.join(output_folder, f"{gwas_sorted_file[:-4]}.plot"))

        with open(os.path.join(raw_gwas_output_path, f"{gwas_sorted_file.split('_')[0]}.tsv")) as raw_gwas_file:
            lines = raw_gwas_file.readlines()
            threshold_denominator= len(lines)-1

        df = pd.read_csv(f"{os.path.join(output_folder, gwas_sorted_file[:-4])}.plot", sep='\t')

        bonferini_adjusted_threshold = 0.05 / threshold_denominator

        threshold = -(math.log(bonferini_adjusted_threshold))

        grid = sns.relplot(data=df, x = 'BP', y = 'log10(p)', hue='log10(p)', palette = 'RdYlGn_r', aspect=1)

        grid.ax.set_xlabel("Position")
        grid.ax.set_ylabel("-log10(p-value)")
        grid.ax.axhline(threshold, linestyle='--', linewidth='1')
        plt.savefig(f"{os.path.join(output_folder, gwas_sorted_file[:-4])}.jpg", dpi = 1200)
