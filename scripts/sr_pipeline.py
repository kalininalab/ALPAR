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
from utils import snippy_runner, prokka_runner, random_name_giver, panaroo_input_creator, panaroo_runner, binary_table_creator, binary_mutation_table_gpa_information_adder, phenotype_dataframe_creator

# SNIPPY VCF EMPTY ISSUE SOLUTION = conda install snippy vt=0.57721

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Single reference AMR is a tool to get mutation and gene presence absence information from genome sequences.")

    # Add the arguments
    parser.add_argument('-i', '--input', type=str, nargs=1, help='txt file that contains path of each strain per line or input folder path (check folder structure)', required=True)

    parser.add_argument('-o', '--output', type=str, nargs=1, help='path of the output folder', required=True)

    parser.add_argument('--reference', type=str, nargs=1, help='path of the reference file', required=True)

    parser.add_argument('--temp', type=str, nargs=1, help='path of the temporary directory', required=True)

    parser.add_argument('--override', action='store_true', help='override the output and temp folder if exists')
    
    parser.add_argument('--keep-temp-files', action='store_true', help='keep the temporary files')

    parser.add_argument('--cpus', type=int, nargs=1, help='number of cpus to use', default=1)

    parser.add_argument('--ram', type=int, nargs=1, help='amount of ram to use in GB', default=4)

    parser.add_argument('--create_phenotype_from_folder', action='store_true', help='create phenotype file from the folders that contains genomic files, folder path should be given with --phenotype_folder option')

    parser.add_argument('--phenotype_folder', type=str, nargs=1, help='folder path to create phenotype file')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('--gwas', action='store_true', help='run gwas analysis')

    # Parse the arguments
    args = parser.parse_args()

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

    # # Run gwas analysis if argument is given
    # if args.gwas:


if __name__ == "__main__":
    main()
