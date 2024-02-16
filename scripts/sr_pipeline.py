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
import multiprocessing
from utils import is_tool_installed, snippy_runner, prokka_runner, random_name_giver, panaroo_input_creator, panaroo_runner, binary_table_creator, binary_mutation_table_gpa_information_adder, phenotype_dataframe_creator, panacota_pipeline_runner, pyseer_runner, pyseer_similarity_matrix_creator, pyseer_phenotype_file_creator, pyseer_genotype_matrix_creator, panacota_pre_processor, panacota_post_processor, temp_folder_remover, binary_table_threshold_with_percentage
from prps import PRPS_runner
from ml import rf_auto_ml, svm, rf, svm_cv, prps_ml_preprecessor, gb_auto_ml, gb

# SNIPPY VCF EMPTY ISSUE SOLUTION = conda install snippy vt=0.57721
def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Single reference AMR is a tool to get mutation and gene presence absence information from genome sequences.")

    parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')

    subparsers = parser.add_subparsers(help='For suggested pipeline, check out github page: https://github.com/kalininalab/SR-AMR')

    parser_main_pipeline = subparsers.add_parser('create_binary_tables', help='from genomic files, create binary mutation table and phenotype table')
    parser_main_pipeline.add_argument('-i', '--input', type=str, help='txt file that contains path of each strain per line or input folder path (check folder structure)', required=True)
    parser_main_pipeline.add_argument('-o', '--output', type=str, help='path of the output folder', required=True)
    parser_main_pipeline.add_argument('--reference', type=str, help='path of the reference file', required=True)
    parser_main_pipeline.add_argument('--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_main_pipeline.add_argument('--overwrite', action='store_true', help='overwrite the output and temp folder if exists, default=False')
    parser_main_pipeline.add_argument('--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_main_pipeline.add_argument('--cpus', type=int, help='number of cpus to use, default=1', default=1)
    parser_main_pipeline.add_argument('--ram', type=int, help='amount of ram to use in GB, default=4', default=4)
    parser_main_pipeline.add_argument('--create_phenotype_from_folder', action='store_true', help='create phenotype file from the folders that contains genomic files, folder path should be given with --phenotype_folder option, default=False')
    parser_main_pipeline.add_argument('--phenotype_folder', type=str, help='folder path to create phenotype file, default=None')
    parser_main_pipeline.set_defaults(func=binary_table_pipeline)

    parser_binary_tables_threshold = subparsers.add_parser('binary_tables_threshold', help='apply threshold to binary mutation table, drops columns that has less than threshold percentage, default=0.2')
    parser_binary_tables_threshold.add_argument('-i', '--input', type=str, help='binary mutation table path', required=True)
    parser_binary_tables_threshold.add_argument('-o', '--output', type=str, help='path of the output folder', required=True)
    parser_binary_tables_threshold.add_argument('--threshold_percentage', type=float, help='threshold percentage value to apply, default=0.2', default=0.2)
    parser_binary_tables_threshold.add_argument('--overwrite', action='store_true', help='overwrite the output folder if exists, default=False')
    parser_binary_tables_threshold.add_argument('--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_binary_tables_threshold.add_argument('--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_binary_tables_threshold.set_defaults(func=binary_table_threshold)

    parser_panacota = subparsers.add_parser('panacota', help='run panacota analysis')
    parser_panacota.add_argument('-i', '--input', type=str, help='txt file that contains path of each strain per line or input folder path, can be found create_binary_tables output path as strains.txt', required=True)
    parser_panacota.add_argument('-o', '--output', type=str, help='path of the output folder', required=True)
    parser_panacota.add_argument('--overwrite', action='store_true', help='overwrite the output folder if exists, default=False')
    parser_panacota.add_argument('--cpus', type=int, help='number of cpus to use, default=1', default=1)
    parser_panacota.add_argument('--name', type=str, help='name of the analysis, default=WIBI', default="WIBI")
    parser_panacota.add_argument('--min_seq_id', type=float, help='Minimum sequence identity to be considered in the same cluster (float between 0 and 1). Default is 0.8', default=0.8)
    parser_panacota.add_argument('--clustering_mode', type=int, help='Choose the clustering mode: 0 for set cover, 1 for single-linkage, 2 for CD-Hit. Default is single-linkage (1)', default=1)
    parser_panacota.add_argument('--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_panacota.add_argument('--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_panacota.add_argument('--random_names_dict', type=str, help='random names dictionary path')
    parser_panacota.add_argument('--data_type', type=str, help='data type of the input, either "nucl" or "prot", default=nucl', default="nucl")
    parser_panacota.set_defaults(func=panacota_pipeline)
    
    parser_gwas = subparsers.add_parser('gwas', help='run gwas analysis')
    parser_gwas.add_argument('-i', '--input', type=str, help='binary mutation table path', required=True)
    parser_gwas.add_argument('-p', '--phenotype', type=str, help='phenotype table path', required=True)
    parser_gwas.add_argument('-t', '--tree', type=str, help='phylogenetic tree path', required=True)
    parser_gwas.add_argument('-o', '--output', type=str, help='path of the output folder', required=True)
    parser_gwas.add_argument('--overwrite', action='store_true', help='overwrite the output folder if exists, default=False')
    parser_gwas.set_defaults(func=gwas_pipeline)

    parser_prps = subparsers.add_parser('prps', help='run prps analysis')
    parser_prps.add_argument('-t', '--tree', type=str, help='phylogenetic tree file path', required=True)
    parser_prps.add_argument('-b', '--binary_mutation_file', type=str, help='binary mutation file path', required=True)
    parser_prps.add_argument('-o', '--output', type=str, help='path of the output folder', required=True)
    parser_prps.add_argument('--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_prps.add_argument('--overwrite', action='store_true', help='overwrite the output and temp folder if exists, default=False')
    parser_prps.add_argument('--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_prps.add_argument('--cpus', type=int, help='number of cpus to use, default=1', default=1)
    parser_prps.set_defaults(func=prps_pipeline)

    parser_ml = subparsers.add_parser('ml', help='run machine learning analysis')
    parser_ml.add_argument('-i', '--input', type=str, help='binary mutation table path', required=True)
    parser_ml.add_argument('-p', '--phenotype', type=str, help='phenotype table path', required=True)
    parser_ml.add_argument('-o', '--output', type=str, help='path of the output folder', required=True)
    parser_ml.add_argument('-a', '--antibiotic', type=str, help='antibiotic name', required=True)
    parser_ml.add_argument('--prps', type=str, help='prps score file path')
    parser_ml.add_argument('--prps_percentage', type=int, help='percentage of the top scores of prps to be used, should be used with --prps option, default=30', default=30)
    parser_ml.add_argument('--overwrite', action='store_true', help='overwrite the output folder if exists, default=False')
    parser_ml.add_argument('--cpus', type=int, help='number of cpus to use, default=1', default=1)
    parser_ml.add_argument('--ram', type=int, help='amount of ram to use in GB, default=4', default=4)
    parser_ml.add_argument('--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_ml.add_argument('--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_ml.add_argument('--ml_algorithm', type=str, help='classification algorithm to be used, available selections: [rf, svm, gb], default=rf', default="rf")
    parser_ml.add_argument('--test_train_split', type=float, help='test train split ratio, default=0.20', default=0.20)
    parser_ml.add_argument('--random_state', type=int, help='random state, default=42', default=42)
    parser_ml.add_argument('--n_estimators', type=int, help='number of estimators for random forest, default=100', default=100)
    parser_ml.add_argument('--max_depth', type=int, help='max depth for random forest, default=10', default=10)
    parser_ml.add_argument('--min_samples_split', type=int, help='min samples split for random forest, default=2', default=2)
    parser_ml.add_argument('--min_samples_leaf', type=int, help='min samples leaf for random forest, default=1', default=1)
    parser_ml.add_argument('--max_features', type=str, help='max features for random forest, default=auto', default="auto")
    parser_ml.add_argument('--resampling_strategy', type=str, help='resampling strategy for ml, available selections: [holdout, cv], default=holdout', default="holdout")
    #parser_ml.add_argument('--bootstrap', action='store_true', help='bootstrap for random forest, default=True')
    parser_ml.add_argument('--parameter_optimization', action='store_true', help='runs parameter optimization, default=False')
    parser_ml.add_argument('--n_jobs', type=int, help='number of jobs for random forest, default=1', default=1)
    parser_ml.add_argument('--cv', type=int, help='applies Cross-Validation with given number of splits, default=4', default=4)
    parser_ml.add_argument('--scoring', type=str, help='scoring method for cross-validation, available selections: [MCC,accuracy,f1,roc_auc], default=MCC', default="MCC")
    parser_ml.add_argument('--save_model', action='store_true', help='save the ml model, default=False')
    parser_ml.add_argument('--feature_importance_analysis', action='store_true', help='analyze feature importance, default=False')
    parser_ml.add_argument('--feature_importance_analysis_number_of_repeats', type=int, help='number of repeats for feature importance analysis should be given with --feature_importance_analysis option, default=5', default=5)
    parser_ml.add_argument('--optimization_time_limit', type=int, help='time limit for parameter optimization with AutoML, default=3600', default=3600)
    parser_ml.add_argument('--svm_kernel', type=str, help='kernel for svm, available selections: [linear, poly, rbf, sigmoid], default=linear', default="linear")

    parser_ml.set_defaults(func=ml_pipeline)


    # Parse the arguments
    args = parser.parse_args()

    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


def run_snippy_and_prokka(strain, random_names, snippy_output, prokka_output, args):
    snippy_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[0]], snippy_output, args.reference, f"{args.temp}/snippy_log.txt", 1, args.ram)
    prokka_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[0]], prokka_output, args.reference, f"{args.temp}/prokka_log.txt", 1)


def binary_table_pipeline(args):
    
    # Sanity checks

    # Check if the tools are installed
    tool_list = ["snippy", "prokka", "panaroo"]

    for tool in tool_list:
        if not is_tool_installed(tool):
            print(f"Error: {tool} is not installed.")
            sys.exit(1)

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
    if os.path.isdir(args.input):
        input_folder = args.input
        input_file = None
    else:
        input_folder = None
        input_file = args.input
    
    # Check if reference file exists and correct extension  
    accepted_reference_extensions = ['.gbk', '.gbff']
    if not os.path.exists(args.reference):
        print("Error: Reference file does not exist.")
        sys.exit(1)
    else:
        if pathlib.Path(args.reference).suffix not in accepted_reference_extensions:
            print("Error: Reference file extension is not accepted.")
            sys.exit(1)
    
    # Check if output folder exists and create if not
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    # Check if output folder empty
    if os.path.exists(args.output) and os.path.isdir(args.output):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)
    
    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        os.mkdir(args.temp)
    
    # Check if temp folder empty
    if os.path.exists(args.temp) and os.path.isdir(args.temp):
        if not args.overwrite:
            print("Error: Temp folder is not empty.")
            sys.exit(1)

    # Check if cpus is positive
    if args.cpus is not None:
        if args.cpus <= 0:
            print("Error: Number of cpus should be positive.")
            sys.exit(1)
    else:
        args.cpus = 1

    # Check if ram is positive
    if args.ram is not None:
        if args.ram <= 0:
            print("Error: Amount of ram should be positive.")
            sys.exit(1)
    else:
        args.ram = 4

    # Check if phenotype folder exists
    if args.create_phenotype_from_folder:
        if not os.path.exists(args.phenotype_folder):
            print("Error: Phenotype folder does not exist.")
            sys.exit(1)

    # Create the output folder
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(os.path.join(args.output, "snippy")):
        os.mkdir(os.path.join(args.output, "snippy"))
    if not os.path.exists(os.path.join(args.output, "prokka")):
        os.mkdir(os.path.join(args.output, "prokka"))
    if not os.path.exists(os.path.join(args.output, "panaroo")):
        os.mkdir(os.path.join(args.output, "panaroo"))

    snippy_output = os.path.join(args.output, "snippy")
    prokka_output = os.path.join(args.output, "prokka")
    panaroo_output = os.path.join(args.output, "panaroo")

    # Create the temp folder
    if not os.path.exists(args.temp):
        os.mkdir(args.temp)

    # Create the temp folder for the panaroo input
    if not os.path.exists(os.path.join(args.temp, "panaroo")):
        os.mkdir(os.path.join(args.temp, "panaroo"))

    if os.path.exists(os.path.join(args.output, "strains.txt")):
        if not args.overwrite:
            print("Error: strains.txt already exists.")
            sys.exit(1)
        else:
            print("Warning: strains.txt already exists. Old one will be deleted..")
            os.remove(os.path.join(args.output, "strains.txt"))

    accepted_fasta_file_extensions = [".fna", ".fasta", ".faa"]

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

            # Checking if folders contain fasta files that are accepted

            files_in_resistant_path = os.listdir(resistant_path)
            files_in_susceptible_path = os.listdir(susceptible_path)

            resistant_strains = []
            susceptible_strains = []
            
            for file in files_in_resistant_path:
                if pathlib.Path(file).suffix in accepted_fasta_file_extensions:
                    resistant_strains.append(file)

            for file in files_in_susceptible_path:
                if pathlib.Path(file).suffix in accepted_fasta_file_extensions:
                    susceptible_strains.append(file)

            
            with open(os.path.join(args.output, "strains.txt"), "a") as outfile:
                for strain in resistant_strains:
                    # Make sure path is same in both Windows and Linux 
                    strain_path = os.path.join(resistant_path, strain)
                    strain_path = strain_path.replace("\\", "/")
                    outfile.write(f"{strain_path}\n")
                for strain in susceptible_strains:
                    strain_path = os.path.join(susceptible_path, strain)
                    strain_path = strain_path.replace("\\", "/")
                    outfile.write(f"{strain_path}\n")

        input_file = os.path.join(args.output, "strains.txt")
    
    if input_file is not None:
        random_names = random_name_giver(input_file, os.path.join(args.output, "random_names.txt"))

    strain_list = []

    # TODO check if strain is added, if yes, don't add and process it again!!!
    with open(input_file, "r") as infile:
        added_strains = []
        lines = infile.readlines()
        for line in lines:
            if os.path.splitext(line.split("/")[-1].strip())[0] not in added_strains:
                added_strains.append(os.path.splitext(line.split("/")[-1].strip())[0])
                strain_list.append(line.strip())
 
    # Run snippy and prokka
    
    print(f"Number of strains to be processed: {len(strain_list)}")
    print("Running snippy and prokka...")

    num_parallel_tasks = args.cpus

    params = [(strain, random_names, snippy_output, prokka_output, args) for strain in strain_list]

    with multiprocessing.Pool(num_parallel_tasks) as pool:
        pool.starmap(run_snippy_and_prokka, params)

    strains_to_be_processed = []

    prokka_output_strains = os.listdir(prokka_output)
    snippy_output_strains = os.listdir(snippy_output)

    strains_to_be_skiped = []

    for strain in random_names.keys():
        if random_names[strain] in prokka_output_strains and random_names[strain] in snippy_output_strains:
            strains_to_be_processed.append(random_names[strain])
        else:
            strains_to_be_skiped.append(random_names[strain])

    print(f"Number of strains processed: {len(strains_to_be_processed)}")
    print(f"Number of strains skipped: {len(strains_to_be_skiped)}")

    print("Creating panaroo input...")
    # Create the panaroo input
    panaroo_input_creator(os.path.join(args.output, "random_names.txt"), prokka_output, os.path.join(args.temp, "panaroo"), strains_to_be_processed)

    print("Running panaroo...")
    # Run panaroo
    panaroo_runner(os.path.join(args.temp, "panaroo"), panaroo_output, os.path.join(args.temp, "panaroo_log.txt"))

    print("Creating binary mutation table...")
    # Create the binary table
    binary_table_creator(snippy_output, os.path.join(args.output, "binary_mutation_table.tsv"), args.cpus, strains_to_be_processed)

    print("Adding gene presence absence information to the binary table...")
    # Add gene presence absence information to the binary table
    binary_mutation_table_gpa_information_adder(os.path.join(args.output, "binary_mutation_table.tsv"), os.path.join(panaroo_output, "gene_presence_absence.csv"), os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv"))

    if args.create_phenotype_from_folder:
        print("Creating phenotype dataframe...")
        # Create the phenotype dataframe
        phenotype_dataframe_creator(input_folder, os.path.join(args.output, "phenotype_table.tsv"), random_names)

    if not args.keep_temp_files:
        temp_folder_remover(os.path.join(args.temp, "panaroo"))


def panacota_pipeline(args):

    tool_list = ["PanACoTA"]

    for tool in tool_list:
        if not is_tool_installed(tool):
            print(f"Error: {tool} is not installed.")
            sys.exit(1)

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        if not os.path.exists(args.temp):
            os.mkdir(args.temp)

    panacota_output = os.path.join(args.output,"panacota")
    panacota_temp = os.path.join(args.temp,"panacota")

    if not os.path.exists(panacota_temp):
            os.mkdir(panacota_temp)
    
    if not os.path.exists(panacota_output):
            os.mkdir(panacota_output)
    
    # Check if temp folder empty
    if os.path.exists(panacota_temp) and os.path.isdir(panacota_temp):
        if not args.overwrite:
            print("Error: Temp folder is not empty.")
            sys.exit(1)
        
    # Check if output folder empty
    if os.path.exists(panacota_output) and os.path.isdir(panacota_output):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)

    panacota_log_file = os.path.join(panacota_output, "panacota_log.txt")

    print(f"Running PanACoTA pipeline pre-precessor...")
    
    panacota_pre_processor(args.input, panacota_temp, panacota_output, args.random_names_dict)

    print(f"Running PanACoTA pipeline with {args.cpus} cores...")

    panacota_pipeline_runner(os.path.join(panacota_output, "panacota_input.lst"), panacota_temp, panacota_output, args.name, args.cpus, panacota_log_file, type=args.data_type, min_seq_id=args.min_seq_id, mode=args.clustering_mode)

    print(f"Running PanACoTA pipeline post-precessor...")

    panacota_post_processor(panacota_output, args.name, args.data_type)

    if not args.keep_temp_files:
        temp_folder_remover(panacota_temp)


def gwas_pipeline(args):

    # Sanity checks

    tool_list = ["pyseer"]

    for tool in tool_list:
        if not is_tool_installed(tool):
            print(f"Error: {tool} is not installed.")
            sys.exit(1)

    # Check the output_folder
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    
    gwas_output = os.path.join(args.output, "gwas")

    # Check if output folder empty
    if os.path.exists(gwas_output) and os.path.isdir(gwas_output):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)
    else:
        os.mkdir(gwas_output)
    
    if not os.path.exists(os.path.join(gwas_output, "gwas_results")):
        if not args.overwrite:
            os.mkdir(os.path.join(gwas_output, "gwas_results"))

    if not os.path.exists(os.path.join(gwas_output, "pyseer_phenotypes")):
        if not args.overwrite:
            os.mkdir(os.path.join(gwas_output, "pyseer_phenotypes"))

    # pyseer_genotype_matrix_creator(binary_mutation_table, output_file):
    pyseer_genotype_matrix_creator(args.input, os.path.join(gwas_output, "genotype_matrix.tsv"))
    # pyseer_phenotype_file_creator(phenotype_file, output_file_directory):
    pyseer_phenotype_file_creator(args.phenotype, os.path.join(gwas_output, "pyseer_phenotypes"))
    # pyseer_similarity_matrix_creator(phylogenetic_tree, output_file):
    pyseer_similarity_matrix_creator(args.tree, os.path.join(gwas_output, "similarity_matrix.tsv"))
    # pyseer_runner(genotype_file_path, phenotype_file_path, similarity_matrix, output_file_directory, cpus):
    pyseer_runner(os.path.join(gwas_output, "genotype_matrix.tsv"), os.path.join(gwas_output, "pyseer_phenotypes"), os.path.join(gwas_output, "similarity_matrix.tsv"), os.path.join(gwas_output, "gwas_results"))


def prps_pipeline(args):

    # Sanity checks

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        if not os.path.exists(args.temp):
            os.mkdir(args.temp)

    prps_output = os.path.join(args.output,"prps")
    prps_temp = os.path.join(args.temp,"prps")

    # Check if output folder empty
    if os.path.exists(prps_output) and os.path.isdir(prps_output):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)
    else:
        os.mkdir(prps_output)

    if not os.path.exists(prps_temp):
        os.mkdir(prps_temp)

    PRPS_runner(args.tree, args.binary_mutation_file, prps_output, prps_temp)

    if not args.keep_temp_files:
        temp_folder_remover(prps_temp)


def ml_pipeline(args):
    
    # Sanity checks

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        if not os.path.exists(args.temp):
            os.mkdir(args.temp)

    ml_output = os.path.join(args.output,"ml")
    ml_temp = os.path.join(args.temp,"ml")

    # Check if output folder empty
    if os.path.exists(ml_output) and os.path.isdir(ml_output):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)
    else:
        os.mkdir(ml_output)

    if not os.path.exists(ml_temp):
        os.mkdir(ml_temp)

    accepted_ml_algorithms = ["rf", "svm"]

    if args.ml_algorithm not in accepted_ml_algorithms:
        print("Error: ML algorithm is not accepted.")
        sys.exit(1)

    accepted_scoring_methods = ["MCC", "accuracy", "f1", "roc_auc"]

    if args.scoring not in accepted_scoring_methods:
        print("Error: Scoring method is not accepted.")
        sys.exit(1)

    PRPS_flag = False

    if args.prps is not None:
        if os.path.exists(args.prps):
            PRPS_flag = True
        else:
            print("Error: PRPS output does not exist.")
            sys.exit(1)
    
    PRPS_percentage = 30

    if args.prps_percentage:

        PRPS_percentage = args.prps_percentage

        if args.prps_percentage < 0 or args.prps_percentage > 100:
            print("Error: PRPS percentage should be between 0 and 100.")
            sys.exit(1)
    
    binary_mutation_table_path = args.input

    if PRPS_flag:
        prps_ml_preprecessor(binary_mutation_table_path, args.prps, PRPS_percentage, ml_temp)

        binary_mutation_table_path = os.path.join(ml_temp, "prps_filtered_table.tsv")

    if float(args.test_train_split) < 0 or float(args.test_train_split) > 1:
        print("Error: Test train split should be between 0 and 1.")
        sys.exit(1)

    if int(args.random_state) < 0:
        print("Error: Random state should be positive.")
        sys.exit(1)

    if int(args.cpus) < 1:
        print("Error: Number of cpus should be positive.")
        sys.exit(1)
    
    if int(args.ram) < 1:
        print("Error: Amount of ram should be positive.")
        sys.exit(1)

    if int(args.feature_importance_analysis_number_of_repeats) < 1:
        print("Error: Number of repeats for feature importance analysis should be positive and bigger than 0.")
        sys.exit(1)

    if args.ml_algorithm == "rf":

        if args.parameter_optimization:
            rf_auto_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.cpus, ml_temp, args.ram, args.optimization_time_limit, args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5)
        else:
            rf(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.cpus, args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf, min_samples_split=args.min_samples_split)

    elif args.ml_algorithm == "svm":

        if args.resampling_strategy == "cv":
            svm_cv(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.test_train_split, ml_output, args.cpus, args.feature_importance_analysis, args.save_model, resampling_strategy="cv", fia_repeats=5, optimization=False)
        elif args.resampling_strategy == "holdout":
            svm(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.test_train_split, ml_output, args.cpus, args.feature_importance_analysis, args.save_model, resampling_strategy="holdout", fia_repeats=5, optimization=False)

    elif args.ml_algorithm == "gb":

        if args.parameter_optimization:
            gb_auto_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.cpus, ml_temp, args.ram, args.optimization_time_limit, args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5)
        else:
            gb(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.cpus, args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf, min_samples_split=args.min_samples_split)

    if not args.keep_temp_files:
        temp_folder_remover(ml_temp)


def binary_table_threshold(args):
    
    # Check the arguments
    if args.input is None:
        print("Error: Input file is required.")
        sys.exit(1)
    
    # Check if threshold is between 0 and 1
    if float(args.threshold_percentage) < 0 or float(args.threshold_percentage) > 100:
        print("Error: Threshold percentage should be between 0 and 100.")
        sys.exit(1)

    # Check if output folder exists and create if not
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    binary_table_threshold_temp = os.path.join(args.temp, "binary_table_threshold")
    binary_table_threshold_output = os.path.join(args.output, "binary_table_threshold")

    # Check if output folder empty
    if os.path.exists(binary_table_threshold_output) and os.path.isdir(binary_table_threshold_output):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        os.mkdir(args.temp)

    # Check if temp folder empty
    if os.path.exists(args.temp) and os.path.isdir(args.temp):
        if not args.overwrite:
            print("Error: Temp folder is not empty.")
            sys.exit(1)

    # Create the temp folder
    if not os.path.exists(binary_table_threshold_temp):
        os.mkdir(binary_table_threshold_temp)

    # Create the output folder
    if not os.path.exists(binary_table_threshold_output):
        os.mkdir(binary_table_threshold_output)

    threshold_percentage_float = float(args.threshold_percentage)
    
    binary_table_threshold_with_percentage(args.input, binary_table_threshold_output, threshold_percentage_float)
    

if __name__ == "__main__":
    main()
