import os
import sys
import argparse
import pathlib
import subprocess
import random
import string
import shutil
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from datetime import datetime
from joblib import Parallel, delayed
import time
import copy
from utils import snippy_runner, prokka_runner, random_name_giver, panaroo_input_creator, panaroo_runner, binary_table_creator, binary_mutation_table_gpa_information_adder, phenotype_dataframe_creator, panacota_pipeline, pyseer_runner, pyseer_similarity_matrix_creator, pyseer_phenotype_file_creator, pyseer_genotype_matrix_creator

# SNIPPY VCF EMPTY ISSUE SOLUTION = conda install snippy vt=0.57721

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Single reference AMR is a tool to get mutation and gene presence absence information from genome sequences.")

    parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')

    subparsers = parser.add_subparsers(help='Suggested pipeline = create_binary_tables -> panacota -> gwas')

    parser_main_pipeline = subparsers.add_parser('create_binary_tables', help='from genomic files, create binary mutation table and phenotype table')
    parser_main_pipeline.add_argument('-i', '--input', type=str, nargs=1, help='txt file that contains path of each strain per line or input folder path (check folder structure)', required=True)
    parser_main_pipeline.add_argument('-o', '--output', type=str, nargs=1, help='path of the output folder', required=True)
    parser_main_pipeline.add_argument('--reference', type=str, nargs=1, help='path of the reference file', required=True)
    parser_main_pipeline.add_argument('--temp', type=str, nargs=1, help='path of the temporary directory', required=True)
    parser_main_pipeline.add_argument('--override', action='store_true', help='override the output and temp folder if exists')
    parser_main_pipeline.add_argument('--keep-temp-files', action='store_true', help='keep the temporary files')
    parser_main_pipeline.add_argument('--cpus', type=int, nargs=1, help='number of cpus to use', default=1)
    parser_main_pipeline.add_argument('--ram', type=int, nargs=1, help='amount of ram to use in GB', default=4)
    parser_main_pipeline.add_argument('--create_phenotype_from_folder', action='store_true', help='create phenotype file from the folders that contains genomic files, folder path should be given with --phenotype_folder option')
    parser_main_pipeline.add_argument('--phenotype_folder', type=str, nargs=1, help='folder path to create phenotype file')
    parser_main_pipeline.set_defaults(func=binary_table_pipeline)

    parser_panacota = subparsers.add_parser('panacota', help='run panacota analysis')
    parser_panacota.add_argument('-i', '--input', type=str, nargs=1, help='txt file that contains path of each strain per line or input folder path, can be found create_binary_tables output path as strains.txt', required=True)
    parser_panacota.add_argument('--reference', type=str, nargs=1, help='reference genome path, either gbk or gbff', required=True)
    parser_panacota.add_argument('-o', '--output', type=str, nargs=1, help='path of the output folder', required=True)
    parser_panacota.add_argument('--override', action='store_true', help='override the output folder if exists')
    parser_panacota.add_argument('--cpus', type=int, nargs=1, help='number of cpus to use', default=1)
    parser_panacota.add_argument('--name', type=str, nargs=1, help='name of the analysis', default="WIBI")
    parser_panacota.add_argument('--min_seq_id', type=float, nargs=1, help='Minimum sequence identity to be considered in the same cluster (float between 0 and 1). Default is 0.8', default=0.8)
    parser_panacota.add_argument('--clustering_mode', type=int, nargs=1, help='Choose the clustering mode: 0 for set cover, 1 for single-linkage, 2 for CD-Hit. Default is single-linkage (1)', default=1)
    parser_panacota.set_defaults(func=panacota_pipeline)
    
    parser_gwas = subparsers.add_parser('gwas', help='run gwas analysis')
    parser_gwas.add_argument('-i', '--input', type=str, nargs=1, help='binary mutation table path', required=True)
    parser_gwas.add_argument('-p', '--phenotype', type=str, nargs=1, help='phenotype table path', required=True)
    parser_gwas.add_argument('-t', '--tree', type=str, nargs=1, help='phylogenetic tree path', required=True)
    parser_gwas.add_argument('-o', '--output', type=str, nargs=1, help='path of the output folder', required=True)
    parser_gwas.add_argument('--override', action='store_true', help='override the output folder if exists')
    parser_gwas.set_defaults(func=gwas_pipeline)

    # Parse the arguments
    args = parser.parse_args()

    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


def binary_table_pipeline(args):
    
    # Sanity checks

    # Check the arguments
    if args.input is None:
        print("Error: Input file is required.")
        sys.exit(1)
    
    # Check if the phenotype file will be created from the folder and folder is given
    if args.create_phenotype_from_folder:
        if args.phenotype_folder is None:
            print("Error: phenotype_folder is required.")
            sys.exit(1)

    # Check if the input is a folder or a file
    if os.path.isdir(args.input[0]):
        input_folder = args.input[0]
        input_file = None
    else:
        input_folder = None
        input_file = args.input[0]
    
    # Check if reference file exists and correct extension  
    accepted_reference_extensions = ['.gbk', '.gbff']
    if not os.path.exists(args.reference[0]):
        print("Error: Reference file does not exist.")
        sys.exit(1)
    else:
        if pathlib.Path(args.reference[0]).suffix not in accepted_reference_extensions:
            print("Error: Reference file extension is not accepted.")
            sys.exit(1)

    # Check if output folder empty
    if os.path.exists(args.output[0]) and os.path.isdir(args.output[0]) and os.listdir(args.output[0]):
        if not args.override:
            print("Error: Output folder is not empty.")
            sys.exit(1)
    
    # Check if temp folder empty
    if os.path.exists(args.temp[0]) and os.path.isdir(args.temp[0]) and os.listdir(args.temp[0]):
        if not args.override:
            print("Error: Temp folder is not empty.")
            sys.exit(1)

    # Check if cpus is positive
    if args.cpus is not None:
        if args.cpus[0] <= 0:
            print("Error: Number of cpus should be positive.")
            sys.exit(1)
    else:
        args.cpus = [1]

    # Check if ram is positive
    if args.ram is not None:
        if args.ram[0] <= 0:
            print("Error: Amount of ram should be positive.")
            sys.exit(1)
    else:
        args.ram = [4]

    # Check if phenotype folder exists
    if args.create_phenotype_from_folder:
        if not os.path.exists(args.phenotype_folder[0]):
            print("Error: Phenotype folder does not exist.")
            sys.exit(1)

    # Create the output folder
    if not os.path.exists(args.output[0]):
        os.mkdir(args.output[0])

    if not os.path.exists(f"{args.output[0]}/snippy"):
        os.mkdir(f"{args.output[0]}/snippy")
    if not os.path.exists(f"{args.output[0]}/prokka"):
        os.mkdir(f"{args.output[0]}/prokka")
    if not os.path.exists(f"{args.output[0]}/panaroo"):
        os.mkdir(f"{args.output[0]}/panaroo")

    snippy_output = f"{args.output[0]}/snippy"
    prokka_output = f"{args.output[0]}/prokka"
    panaroo_output = f"{args.output[0]}/panaroo"

    # Create the temp folder
    if not os.path.exists(args.temp[0]):
        os.mkdir(args.temp[0])

    # Create the temp folder for the panaroo input
    if not os.path.exists(f"{args.temp[0]}/panaroo"):
        os.mkdir(f"{args.temp[0]}/panaroo")

    if input_folder is not None:
        antibiotics = os.listdir(input_folder)
        for antibiotic in antibiotics:
            antibiotic_path = os.path.join(input_folder, antibiotic)
            status = os.listdir(antibiotic_path)
            if not 'Resistant' in status:
                print(f"Error: {antibiotic} folder does not contain resistant folder.")
                sys.exit(1)
            if not 'Susceptible' in status:
                print(f"Error: {antibiotic} folder does not contain susceptible folder.")
                sys.exit(1)
            
            resistant_path = os.path.join(antibiotic_path, 'Resistant')
            susceptible_path = os.path.join(antibiotic_path, 'Susceptible')

            resistant_strains = os.listdir(resistant_path)
            susceptible_strains = os.listdir(susceptible_path)

            with open(f"{args.temp[0]}/strains.txt", "w") as outfile:
                for strain in resistant_strains:
                    # Make sure path is same in both Windows and Linux 
                    strain_path = os.path.join(resistant_path, strain)
                    strain_path = strain_path.replace("\\", "/")
                    outfile.write(f"{strain_path}\n")
                for strain in susceptible_strains:
                    strain_path = os.path.join(susceptible_path, strain)
                    strain_path = strain_path.replace("\\", "/")
                    outfile.write(f"{strain_path}\n")

        input_file = f"{args.temp[0]}/strains.txt"
    
    if input_file is not None:
        random_names = random_name_giver(input_file, f"{args.temp[0]}/random_names.txt")

    strain_list = []

    with open(input_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            strain_list.append(line.strip())
 
    # Run snippy and prokka
    
    print(f"Number of strains to be processed: {len(strain_list)}")
    print("Running snippy and prokka...")
        
    for strain in strain_list:
        # input, output, reference, cpus = 1, memory = 4, parallel_run = False
        snippy_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[0]] ,snippy_output, args.reference[0], f"{args.temp[0]}/snippy_log.txt" ,args.cpus[0], args.ram[0])
        # input, output, reference, cpus = 1, parallel_run = False
        prokka_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[0]], prokka_output, args.reference[0], f"{args.temp[0]}/prokka_log.txt", args.cpus[0])
        #break

    print("Creating panaroo input...")
    # Create the panaroo input
    panaroo_input_creator(f"{args.temp[0]}/random_names.txt", prokka_output, f"{args.temp[0]}/panaroo")

    print("Running panaroo...")
    # Run panaroo
    panaroo_runner(f"{args.temp[0]}/panaroo", panaroo_output, f"{args.temp[0]}/panaroo_log.txt")

    print("Creating binary mutation table...")
    # Create the binary table
    binary_table_creator (snippy_output, f"{args.output[0]}/binary_mutation_table.tsv", args.cpus[0])

    print("Adding gene presence absence information to the binary table...")
    # Add gene presence absence information to the binary table
    binary_mutation_table_gpa_information_adder(f"{args.output[0]}/binary_mutation_table.tsv", f"{panaroo_output}/gene_presence_absence.csv", f"{args.output[0]}/binary_mutation_table_with_gene_presence_absence.tsv")

    if args.create_phenotype_from_folder:
        print("Creating phenotype dataframe...")
        # Create the phenotype dataframe
        phenotype_dataframe_creator(input_folder, f"{args.output[0]}/phenotype_table.tsv", random_names)


def panacota_pipeline(args):

    panacota_pipeline(args.input[0], args.reference[0], args.output[0], args.name[0], args.cpus[0])


def gwas_pipeline(args):

    # Sanity checks

    # Check the output_folder
    if not os.path.exists(args.output[0]):
        os.mkdir(args.output[0])
        os.mkdir(f"{args.output[0]}/gwas_output")
        os.mkdir(f"{args.output[0]}/pyseer_phenotypes")

    # Check if output folder empty
    if os.path.exists(args.output[0]) and os.path.isdir(args.output[0]) and os.listdir(args.output[0]):
        if not args.override:
            print("Error: Output folder is not empty.")
            sys.exit(1)

    # pyseer_genotype_matrix_creator(binary_mutation_table, output_file):
    pyseer_genotype_matrix_creator(args.input[0], f"{args.output[0]}/genotype_matrix.tsv")
    # pyseer_phenotype_file_creator(phenotype_file, output_file_directory):
    pyseer_phenotype_file_creator(args.phenotype[0], f"{args.output[0]}/pyseer_phenotypes/")
    # pyseer_similarity_matrix_creator(phylogenetic_tree, output_file):
    pyseer_similarity_matrix_creator(args.tree[0], f"{args.output[0]}/similarity_matrix.tsv")
    # pyseer_runner(genotype_file_path, phenotype_file_path, similarity_matrix, output_file_directory, cpus):
    pyseer_runner(f"{args.output[0]}/genotype_matrix.tsv", f"{args.output[0]}/pyseer_phenotypes/", f"{args.output[0]}/similarity_matrix.tsv", f"{args.output[0]}/gwas_output")



if __name__ == "__main__":
    main()
