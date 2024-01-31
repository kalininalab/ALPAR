import os
import random
import pandas as pd
#import pyarrow as pa
import shutil
import string
import subprocess
import sys
import time
import argparse
import pathlib
from snippy_runner import snippy_runner
from prokka_runner import prokka_runner
from rename_files_with_random_names import random_name_giver


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

    parser.add_argument('-c', '--cpus', type=int, nargs=1, help='number of cpus to use', default=1)

    parser.add_argument('-r', '--ram', type=int, nargs=1, help='amount of ram to use in GB', default=4)

    parser.add_argument('--create_phenotype_from_folder', action='store_true', help='create phenotype file from the folders that contains genomic files, folder path should be given with --phenotype_folder option')

    parser.add_argument('--phenotype_folder', type=str, nargs=1, help='folder path to create phenotype file')

    # Parse the arguments
    args = parser.parse_args()

    # # Test prints
    # print(args.input)
    # print(args.output)
    # print(args.reference)
    # print(args.temp)
    # print(args.override)

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
    if args.cpus <= 0:
        print("Error: Number of cpus should be positive.")
        sys.exit(1)

    # Check if ram is positive
    if args.ram <= 0:
        print("Error: Amount of ram should be positive.")
        sys.exit(1)

    # Check if phenotype folder exists
    if args.create_phenotype_from_folder:
        if not os.path.exists(args.phenotype_folder[0]):
            print("Error: Phenotype folder does not exist.")
            sys.exit(1)

    # Create the output folder
    if not os.path.exists(args.output[0]):
        os.mkdir(args.output[0])

    # Create the temp folder
    if not os.path.exists(args.temp[0]):
        os.mkdir(args.temp[0])

    # Create the temp folder for the snippy output
    if not os.path.exists(f"{args.temp[0]}/snippy"):
        os.mkdir(f"{args.temp[0]}/snippy")
    
    snippy_temp = f"{args.temp[0]}/snippy"

    # Create the temp folder for the prokka output
    if not os.path.exists(f"{args.temp[0]}/prokka"):
        os.mkdir(f"{args.temp[0]}/prokka")

    prokka_temp = f"{args.temp[0]}/prokka"

    if input_folder is not None:
        antibiotics = os.listdir(input_folder)
        for antibiotic in antibiotics:
            antibiotic_path = os.path.join(input_folder, antibiotic)
            status = os.listdir(antibiotic_path)
            if not 'resistant' in status:
                print(f"Error: {antibiotic} folder does not contain resistant folder.")
                sys.exit(1)
            if not 'susceptible' in status:
                print(f"Error: {antibiotic} folder does not contain susceptible folder.")
                sys.exit(1)
            
            resistant_path = os.path.join(antibiotic_path, 'resistant')
            susceptible_path = os.path.join(antibiotic_path, 'susceptible')

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

        
    # Run snippy




if __name__ == "__main__":
    main()
