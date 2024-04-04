#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import os
import sys
import argparse
import pathlib
import contextlib
import time
import multiprocessing

from utils import is_tool_installed, temp_folder_remover, time_function

try:
    from panacota import panacota_pre_processor, panacota_post_processor, panacota_pipeline_runner
    from gwas import pyseer_runner, pyseer_similarity_matrix_creator, pyseer_phenotype_file_creator, pyseer_genotype_matrix_creator, pyseer_post_processor, pyseer_gwas_graph_creator
    from binary_tables import snippy_runner, prokka_runner, random_name_giver, panaroo_input_creator, panaroo_runner, binary_table_creator, binary_mutation_table_gpa_information_adder, phenotype_dataframe_creator, phenotype_dataframe_creator_post_processor, prokka_create_database, snippy_processed_file_creator
    from binary_table_threshold import binary_table_threshold_with_percentage
    from phylogeny_tree import mash_preprocessor, mash_distance_runner
    from prps import PRPS_runner
    from ds import datasail_runner, datasail_pre_precessor
    from ml import rf_auto_ml, svm, rf, svm_cv, prps_ml_preprecessor, gb_auto_ml, gb
    isLite = False
    # print("Full version is running.")

except ImportError as e:
    print(e)
    from binary_tables import snippy_runner, prokka_runner, random_name_giver, panaroo_input_creator, panaroo_runner, binary_table_creator, binary_mutation_table_gpa_information_adder, phenotype_dataframe_creator, phenotype_dataframe_creator_post_processor, prokka_create_database, snippy_processed_file_creator
    from binary_table_threshold import binary_table_threshold_with_percentage
    from ml_lite import svm, rf, svm_cv, prps_ml_preprecessor, gb
    isLite = True
    # print("Lite version is running.")


def main():
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Single reference AMR is a tool to get mutation and gene presence absence information from genome sequences.")

    if isLite:
        parser.add_argument('--version', action='version',
                            version='%(prog)s 0.1.1 Lite')
    else:
        parser.add_argument('--version', action='version',
                            version='%(prog)s 0.1.1')

    subparsers = parser.add_subparsers(
        help='For suggested pipeline, check out our github page: https://github.com/kalininalab/SR-AMR')

    parser_main_pipeline = subparsers.add_parser(
        'create_binary_tables', help='from genomic files, create binary mutation table and phenotype table')
    parser_main_pipeline.add_argument(
        '-i', '--input', type=str, help='txt file that contains path of each strain per line or input folder path (check folder structure)', required=True)
    parser_main_pipeline.add_argument(
        '-o', '--output', type=str, help='path of the output folder', required=True)
    parser_main_pipeline.add_argument(
        '--reference', type=str, help='path of the reference file', required=True)
    parser_main_pipeline.add_argument('--create_phenotype_from_folder', type=str,
                                      help='create phenotype file from the folders that contains genomic files, folder path should be given with the option, default=None')
    parser_main_pipeline.add_argument('--custom_database', type=str, nargs=2,
                                      help='creates and uses custom database for prokka, require path of the fasta file and genus name, default=None')
    parser_main_pipeline.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_main_pipeline.add_argument('--overwrite', action='store_true',
                                      help='overwrite the output and temp folder if exists, default=False')
    parser_main_pipeline.add_argument(
        '--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_main_pipeline.add_argument(
        '--threads', type=int, help='number of threads to use, default=1', default=1)
    parser_main_pipeline.add_argument(
        '--ram', type=int, help='amount of ram to use in GB, default=4', default=4)
    parser_main_pipeline.add_argument('--no_gene_presence_absence', action='store_true',
                                      help='do not run gene presence absence functions, default=False')
    parser_main_pipeline.add_argument(
        '--no_gene_annotation', action='store_true', help='do not run gene annotation, default=False')
    parser_main_pipeline.set_defaults(func=binary_table_pipeline)

    parser_phenotype_table = subparsers.add_parser(
        'create_phenotype_table', help='from specific folder structure, create phenotype table, does not needed to be run if create_phenotype_from_folder used in create_binary_tables')
    parser_phenotype_table.add_argument(
        '-i', '--input', type=str, help='input folder path (check folder structure)', required=True)
    parser_phenotype_table.add_argument(
        '-o', '--output', type=str, help='path of the output folder', required=True)
    parser_phenotype_table.add_argument(
        '--random_names_dict', type=str, help='random names dictionary path')
    parser_phenotype_table.add_argument(
        '--overwrite', action='store_true', help='overwrite the output and temp folder if exists, default=False')
    parser_phenotype_table.set_defaults(func=phenotype_table_pipeline)

    parser_binary_tables_threshold = subparsers.add_parser(
        'binary_table_threshold', help='apply threshold to binary mutation table, drops columns that has less than threshold percentage, default=0.2')
    parser_binary_tables_threshold.add_argument(
        '-i', '--input', type=str, help='binary mutation table path', required=True)
    parser_binary_tables_threshold.add_argument(
        '-o', '--output', type=str, help='path of the output folder', required=True)
    parser_binary_tables_threshold.add_argument(
        '--threshold_percentage', type=float, help='threshold percentage value to apply, default=0.2', default=0.2)
    parser_binary_tables_threshold.add_argument(
        '--overwrite', action='store_true', help='overwrite the output folder if exists, default=False')
    parser_binary_tables_threshold.add_argument(
        '--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_binary_tables_threshold.set_defaults(func=binary_table_threshold)

    parser_panacota = subparsers.add_parser(
        'panacota', help='run panacota analysis')
    parser_panacota.add_argument(
        '-i', '--input', type=str, help='txt file that contains path of each strain per line or input folder path, can be found create_binary_tables output path as strains.txt', required=True)
    parser_panacota.add_argument(
        '-o', '--output', type=str, help='path of the output folder', required=True)
    parser_panacota.add_argument(
        '--random_names_dict', type=str, help='random names dictionary path')
    parser_panacota.add_argument(
        '--data_type', type=str, help='data type of the input, either "nucl" or "prot", default=nucl', default="nucl")
    parser_panacota.add_argument('--overwrite', action='store_true',
                                 help='overwrite the output folder if exists, default=False')
    parser_panacota.add_argument(
        '--threads', type=int, help='number of threads to use, default=1', default=1)
    parser_panacota.add_argument(
        '--name', type=str, help='name of the analysis, default=WIBI', default="WIBI")
    parser_panacota.add_argument(
        '--min_seq_id', type=float, help='Minimum sequence identity to be considered in the same cluster (float between 0 and 1). Default is 0.8', default=0.8)
    parser_panacota.add_argument('--core_genome_percentage', type=float,
                                 help='Percentage of core genome to be considered as core genome, default=1', default=1)
    parser_panacota.add_argument('--clustering_mode', type=int,
                                 help='Choose the clustering mode: 0 for set cover, 1 for single-linkage, 2 for CD-Hit. Default is single-linkage (1)', default=1)
    parser_panacota.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_panacota.add_argument(
        '--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_panacota.set_defaults(func=panacota_pipeline)

    parser_phylogenetic_tree = subparsers.add_parser(
        'phylogenetic_tree', help='create phylogenetic tree from genomic fasta files via mashtree')
    parser_phylogenetic_tree.add_argument(
        '-i', '--input', type=str, help='txt file that contains path of each strain per line or input folder path, can be found create_binary_tables output path as strains.txt', required=True)
    parser_phylogenetic_tree.add_argument(
        '-o', '--output', type=str, help='path of the output folder', required=True)
    parser_phylogenetic_tree.add_argument(
        '--random_names_dict', type=str, help='random names dictionary path')
    parser_phylogenetic_tree.add_argument(
        '--overwrite', action='store_true', help='overwrite the output folder if exists, default=False')
    parser_phylogenetic_tree.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_phylogenetic_tree.add_argument(
        '--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_phylogenetic_tree.set_defaults(func=phylogenetic_tree_pipeline)

    parser_gwas = subparsers.add_parser('gwas', help='run gwas analysis')
    parser_gwas.add_argument('-i', '--input', type=str,
                             help='binary mutation table path', required=True)
    parser_gwas.add_argument('-o', '--output', type=str,
                             help='path of the output folder', required=True)
    parser_gwas.add_argument(
        '-p', '--phenotype', type=str, help='phenotype table path', required=True)
    parser_gwas.add_argument('-t', '--tree', type=str,
                             help='phylogenetic tree path', required=True)
    parser_gwas.add_argument('--overwrite', action='store_true',
                             help='overwrite the output folder if exists, default=False')
    parser_gwas.set_defaults(func=gwas_pipeline)

    parser_prps = subparsers.add_parser('prps', help='run prps analysis')
    parser_prps.add_argument('-i', '--input', type=str,
                             help='binary mutation file path', required=True)
    parser_prps.add_argument('-o', '--output', type=str,
                             help='path of the output folder', required=True)
    parser_prps.add_argument('-t', '--tree', type=str,
                             help='phylogenetic tree file path', required=True)
    parser_prps.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_prps.add_argument('--overwrite', action='store_true',
                             help='overwrite the output and temp folder if exists, default=False')
    parser_prps.add_argument(
        '--threads', type=int, help='number of threads to use, default=1', default=1)
    parser_prps.add_argument('--keep_temp_files', action='store_true',
                             help='keep the temporary files, default=False')
    parser_prps.set_defaults(func=prps_pipeline)

    parser_ml = subparsers.add_parser(
        'ml', help='run machine learning analysis')
    parser_ml.add_argument('-i', '--input', type=str,
                           help='binary mutation table path', required=True)
    parser_ml.add_argument('-o', '--output', type=str,
                           help='path of the output folder', required=True)
    parser_ml.add_argument('-p', '--phenotype', type=str,
                           help='phenotype table path', required=True)
    parser_ml.add_argument('-a', '--antibiotic', type=str,
                           help='antibiotic name', required=True)
    parser_ml.add_argument('--prps', type=str, help='prps score file path')
    parser_ml.add_argument('--prps_percentage', type=int,
                           help='percentage of the top scores of prps to be used, should be used with --prps option, default=30', default=30)
    parser_ml.add_argument('--sail', type=str, help='split against information leakage, requires txt file that contains path of each strain per line or input folder path, can be found create_binary_tables output path as strains.txt, default=None')
    parser_ml.add_argument('--overwrite', action='store_true',
                           help='overwrite the output folder if exists, default=False')
    parser_ml.add_argument(
        '--threads', type=int, help='number of threads to use, default=1', default=1)
    parser_ml.add_argument(
        '--ram', type=int, help='amount of ram to use in GB, default=4', default=4)
    parser_ml.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_ml.add_argument('--keep_temp_files', action='store_true',
                           help='keep the temporary files, default=False')
    parser_ml.add_argument('--ml_algorithm', type=str,
                           help='classification algorithm to be used, available selections: [rf, svm, gb], default=rf', default="rf")
    parser_ml.add_argument('--test_train_split', type=float,
                           help='test train split ratio, default=0.20', default=0.20)
    parser_ml.add_argument('--random_state', type=int,
                           help='random state, default=42', default=42)
    parser_ml.add_argument('--n_estimators', type=int,
                           help='number of estimators for random forest, default=100', default=100)
    parser_ml.add_argument('--max_depth', type=int,
                           help='max depth for random forest, default=10', default=10)
    parser_ml.add_argument('--min_samples_split', type=int,
                           help='min samples split for random forest, default=2', default=2)
    parser_ml.add_argument('--min_samples_leaf', type=int,
                           help='min samples leaf for random forest, default=1', default=1)
    parser_ml.add_argument('--max_features', type=str,
                           help='max features for random forest, default=auto', default="auto")
    parser_ml.add_argument('--resampling_strategy', type=str,
                           help='resampling strategy for ml, available selections: [holdout, cv], default=holdout', default="holdout")
    parser_ml.add_argument('--parameter_optimization', action='store_true',
                           help='runs parameter optimization, default=False')
    parser_ml.add_argument(
        '--cv', type=int, help='applies Cross-Validation with given number of splits, default=4', default=4)
    parser_ml.add_argument(
        '--scoring', type=str, help='scoring method for cross-validation, available selections: [MCC,accuracy,f1,roc_auc], default=MCC', default="MCC")
    parser_ml.add_argument('--save_model', action='store_true',
                           help='save the ml model, default=False')
    parser_ml.add_argument('--feature_importance_analysis', action='store_true',
                           help='analyze feature importance, default=False')
    parser_ml.add_argument('--feature_importance_analysis_number_of_repeats', type=int,
                           help='number of repeats for feature importance analysis should be given with --feature_importance_analysis option, default=5', default=5)
    parser_ml.add_argument('--optimization_time_limit', type=int,
                           help='time limit for parameter optimization with AutoML, default=3600', default=3600)
    parser_ml.add_argument('--svm_kernel', type=str,
                           help='kernel for svm, available selections: [linear, poly, rbf, sigmoid], default=linear', default="linear")

    parser_ml.set_defaults(func=ml_pipeline)

    # Parse the arguments
    args = parser.parse_args()

    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


def run_snippy_and_prokka(strain, random_names, snippy_output, prokka_output, args, snippy_flag, prokka_flag, custom_db=None):
    if snippy_flag:
        snippy_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[
                      0]], snippy_output, args.reference, f"{args.temp}/snippy_log.txt", 1, args.ram)
    if prokka_flag:
        prokka_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[
                      0]], prokka_output, args.reference, f"{args.temp}/prokka_log.txt", 1, custom_db)


def binary_table_pipeline(args):

    start_time = time.time()

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

    # Check if output folder empty
    if os.path.exists(args.output) and os.path.isdir(args.output):
        if os.listdir(args.output) and not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)

    # Check if output folder exists and create if not
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    temp_folder_created = False

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        if not os.path.exists(args.temp):
            os.mkdir(args.temp)
            temp_folder_created = True

    if not temp_folder_created:
        # Check if temp folder empty
        if os.path.exists(args.temp) and os.path.isdir(args.temp):
            if os.listdir(args.temp) and not args.overwrite:
                print("Error: Temp folder is not empty.")
                sys.exit(1)

    # Check if threads is positive
    if args.threads is not None:
        if args.threads <= 0:
            print("Error: Number of threads should be positive.")
            sys.exit(1)
    else:
        args.threads = 1

    # Check if ram is positive
    if args.ram is not None:
        if args.ram <= 0:
            print("Error: Amount of ram should be positive.")
            sys.exit(1)
    else:
        args.ram = 4

    # Check if phenotype folder exists
    if args.create_phenotype_from_folder:
        if not os.path.exists(args.create_phenotype_from_folder):
            print("Error: Phenotype folder does not exist.")
            sys.exit(1)

    if args.no_gene_annotation:
        if not args.no_gene_presence_absence:
            print(
                "Error: If gene annotation is not run, gene presence absence can not be run.")
            sys.exit(1)

    snippy_flag = True
    prokka_flag = True
    panaroo_flag = True

    if args.no_gene_presence_absence:
        panaroo_flag = False
    if args.no_gene_annotation:
        prokka_flag = False

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
            if antibiotic.startswith("."):
                continue
            antibiotic_path = os.path.join(input_folder, antibiotic)
            status = os.listdir(antibiotic_path)
            if not 'Resistant' in status:
                print(
                    f"Error: {antibiotic} folder does not contain resistant folder.")
                sys.exit(1)
            if not 'Susceptible' in status:
                print(
                    f"Error: {antibiotic} folder does not contain susceptible folder.")
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
        random_names = random_name_giver(
            input_file, os.path.join(args.output, "random_names.txt"))

    strain_list = []

    with open(input_file, "r") as infile:
        added_strains = []
        lines = infile.readlines()
        for line in lines:
            if os.path.splitext(line.split("/")[-1].strip())[0] not in added_strains:
                added_strains.append(os.path.splitext(
                    line.split("/")[-1].strip())[0])
                strain_list.append(line.strip())

    if args.custom_database:
        if len(args.custom_database) != 2:
            print("Error: Custom database option should have two arguments.")
            sys.exit(1)

        if not os.path.exists(args.custom_database[0]):
            print("Error: Custom database fasta file does not exist.")
            sys.exit(1)

        print("Creating custom database...")
        prokka_create_database(
            args.custom_database[0], args.custom_database[1], args.temp, args.threads, args.ram)
        print("Custom database created.")

    # Run snippy and prokka

    print(f"Number of strains to be processed: {len(strain_list)}")
    print("Running snippy and prokka...")

    num_parallel_tasks = args.threads

    params = [(strain, random_names, snippy_output, prokka_output,
               args, snippy_flag, prokka_flag) for strain in strain_list]

    if args.custom_database:
        params = [(strain, random_names, snippy_output, prokka_output, args,
                   snippy_flag, prokka_flag, args.custom_database[1]) for strain in strain_list]

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

    if panaroo_flag:

        print("Creating panaroo input...")
        # Create the panaroo input
        panaroo_input_creator(os.path.join(args.output, "random_names.txt"), prokka_output, os.path.join(
            args.temp, "panaroo"), strains_to_be_processed)

        print("Running panaroo...")
        # Run panaroo
        panaroo_runner(os.path.join(args.temp, "panaroo"), panaroo_output, os.path.join(
            args.temp, "panaroo_log.txt"), args.threads)

        print("Creating binary mutation table...")
        # Create the binary table
        binary_table_creator(snippy_output, os.path.join(
            args.output, "binary_mutation_table.tsv"), args.threads, strains_to_be_processed)

        print("Adding gene presence absence information to the binary table...")
        # Add gene presence absence information to the binary table
        binary_mutation_table_gpa_information_adder(os.path.join(args.output, "binary_mutation_table.tsv"), os.path.join(
            panaroo_output, "gene_presence_absence.csv"), os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv"))

        if args.create_phenotype_from_folder:
            print("Creating phenotype dataframe...")
            # Create the phenotype dataframe
            phenotype_dataframe_creator(args.create_phenotype_from_folder, os.path.join(
                args.output, "phenotype_table.tsv"), random_names)

            phenotype_dataframe_creator_post_processor(os.path.join(
                args.output, "binary_mutation_table.tsv"), os.path.join(args.output, "phenotype_table.tsv"))

    snippy_processed_file_creator(snippy_output, os.path.join(
        args.output, "snippy_processed_strains.txt"))

    if args.keep_temp_files:
        print("Warning, temp files will be kept this might take up space.")

    if not args.keep_temp_files:
        temp_folder_remover(os.path.join(args.temp))
        temp_folder_remover(os.path.join(args.output, "snippy"))
        temp_folder_remover(os.path.join(args.output, "prokka"))
        temp_folder_remover(os.path.join(args.output, "panaroo"))

    print("Done")

    end_time = time.time()

    print(time_function(start_time, end_time))


def panacota_pipeline(args):

    start_time = time.time()

    tool_list = ["PanACoTA"]

    for tool in tool_list:
        if not is_tool_installed(tool):
            print(f"Error: {tool} is not installed.")
            sys.exit(1)

    temp_folder_created = False

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        if not os.path.exists(args.temp):
            os.mkdir(args.temp)
            temp_folder_created = True

    panacota_output = os.path.join(args.output, "panacota")
    panacota_temp = os.path.join(args.temp, "panacota")

    if not os.path.exists(panacota_temp):
        os.mkdir(panacota_temp)
        temp_folder_created = True

    output_folder_created = False

    if not os.path.exists(panacota_output):
        os.mkdir(panacota_output)
        output_folder_created = True

    if not temp_folder_created:
        # Check if temp folder empty
        if os.path.exists(panacota_temp) and os.path.isdir(panacota_temp):
            if not args.overwrite:
                print("Error: Temp folder is not empty.")
                sys.exit(1)

    if not output_folder_created:
        # Check if output folder empty
        if os.path.exists(panacota_output) and os.path.isdir(panacota_output):
            if not args.overwrite:
                print("Error: Output folder is not empty.")
                sys.exit(1)

    panacota_log_file = os.path.join(panacota_output, "panacota_log.txt")

    print(f"Running PanACoTA pipeline pre-precessor...")

    panacota_pre_processor(args.input, panacota_temp,
                           panacota_output, args.random_names_dict)

    print(f"Running PanACoTA pipeline with {args.threads} cores...")

    panacota_pipeline_runner(os.path.join(panacota_output, "panacota_input.lst"), panacota_temp, panacota_output, args.name, args.threads, panacota_log_file,
                             type=args.data_type, min_seq_id=args.min_seq_id, mode=args.clustering_mode, core_genome_percentage=args.core_genome_percentage)

    print(f"Running PanACoTA pipeline post-precessor...")

    panacota_post_processor(panacota_output, args.name,
                            args.output, args.data_type)

    if not args.keep_temp_files:
        print("Removing temp folder...")
        temp_folder_remover(panacota_temp)

    print("Done")

    end_time = time.time()

    print(time_function(start_time, end_time))


def gwas_pipeline(args):

    start_time = time.time()

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
        os.mkdir(os.path.join(gwas_output, "gwas_results"))

    if not os.path.exists(os.path.join(gwas_output, "pyseer_phenotypes")):
        os.mkdir(os.path.join(gwas_output, "pyseer_phenotypes"))

    if not os.path.exists(os.path.join(gwas_output, "graphs")):
        os.mkdir(os.path.join(gwas_output, "graphs"))

    # pyseer_genotype_matrix_creator(binary_mutation_table, output_file):
    pyseer_genotype_matrix_creator(
        args.input, os.path.join(gwas_output, "genotype_matrix.tsv"))
    # pyseer_phenotype_file_creator(phenotype_file, output_file_directory):
    pyseer_phenotype_file_creator(
        args.phenotype, os.path.join(gwas_output, "pyseer_phenotypes"))
    # pyseer_similarity_matrix_creator(phylogenetic_tree, output_file):
    pyseer_similarity_matrix_creator(
        args.tree, os.path.join(gwas_output, "similarity_matrix.tsv"))
    # pyseer_runner(genotype_file_path, phenotype_file_path, similarity_matrix, output_file_directory, threads):
    pyseer_runner(os.path.join(gwas_output, "genotype_matrix.tsv"), os.path.join(gwas_output, "pyseer_phenotypes"),
                  os.path.join(gwas_output, "similarity_matrix.tsv"), os.path.join(gwas_output, "gwas_results"))

    if not os.path.exists(os.path.join(gwas_output, "sorted")):
        os.mkdir(os.path.join(gwas_output, "sorted"))

    if not os.path.exists(os.path.join(gwas_output, "sorted_cleaned")):
        os.mkdir(os.path.join(gwas_output, "sorted_cleaned"))

    pyseer_post_processor(os.path.join(
        gwas_output, "gwas_results"), gwas_output)

    pyseer_gwas_graph_creator(gwas_output, os.path.join(gwas_output, "graphs"))

    print("Done")

    end_time = time.time()

    print(time_function(start_time, end_time))


def prps_pipeline(args):

    start_time = time.time()

    # Sanity checks

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        if not os.path.exists(args.temp):
            os.mkdir(args.temp)

    prps_output = os.path.join(args.output, "prps")
    prps_temp = os.path.join(args.temp, "prps")

    # Check if output folder empty
    if os.path.exists(prps_output) and os.path.isdir(prps_output):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)
    else:
        os.mkdir(prps_output)

    if not os.path.exists(prps_temp):
        os.mkdir(prps_temp)

    print("Running PRPS...")

    PRPS_runner(args.tree, args.input, prps_output, prps_temp)

    if not args.keep_temp_files:
        print("Removing temp folder...")
        temp_folder_remover(prps_temp)

    print("Done")

    end_time = time.time()

    print(time_function(start_time, end_time))


def ml_pipeline(args):

    start_time = time.time()

    # Sanity checks

    # Check if output folder exists and create if not
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        if not os.path.exists(args.temp):
            os.mkdir(args.temp)

    ml_output = os.path.join(args.output, "ml")
    ml_temp = os.path.join(args.temp, "ml")

    # Check if output folder empty
    if os.path.exists(ml_output) and os.path.isdir(ml_output):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)
    else:
        os.mkdir(ml_output)

    if not os.path.exists(ml_temp):
        os.mkdir(ml_temp)

    accepted_ml_algorithms = ["rf", "svm", "gb"]

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

        print("Running PRPS ML pre-precessor...")
        prps_ml_preprecessor(binary_mutation_table_path,
                             args.prps, PRPS_percentage, ml_temp)

        binary_mutation_table_path = os.path.join(
            ml_temp, "prps_filtered_table.tsv")

    if float(args.test_train_split) < 0 or float(args.test_train_split) > 1:
        print("Error: Test train split should be between 0 and 1.")
        sys.exit(1)

    if int(args.random_state) < 0:
        print("Error: Random state should be positive.")
        sys.exit(1)

    if int(args.threads) < 1:
        print("Error: Number of threads should be positive.")
        sys.exit(1)

    if int(args.ram) < 1:
        print("Error: Amount of ram should be positive.")
        sys.exit(1)

    if int(args.feature_importance_analysis_number_of_repeats) < 1:
        print("Error: Number of repeats for feature importance analysis should be positive and bigger than 0.")
        sys.exit(1)

    train_strains = []
    test_strains = []

    if args.sail:

        datasail_temp = os.path.join(args.temp, "datasail")
        datasail_output = os.path.join(args.output, "datasail")

        if os.path.exists(os.path.join(datasail_output, "splits.tsv")):
            print("Warning: Split file already exists, it will be used for calculations. If you want to re-run the datasail, please remove the splits.tsv file from the output folder.")

        else:
            # Check if output folder empty
            if os.path.exists(datasail_output) and os.path.isdir(datasail_output):
                if not args.overwrite:
                    print("Error: Output folder is not empty.")
                    sys.exit(1)

            # Create the temp folder
            if not os.path.exists(datasail_temp):
                os.mkdir(datasail_temp)

            # Create the output folder
            if not os.path.exists(datasail_output):
                os.mkdir(datasail_output)

            if os.path.exists(f"{os.path.dirname(args.sail)}/random_names.txt"):
                random_names_dict = f"{os.path.dirname(args.sail)}/random_names.txt"
            else:
                random_names_dict = None

            if args.test_train_split:
                train_test = [float(1-args.test_train_split),
                              float(args.test_train_split)]

            print("Creating distance matrix...")

            datasail_pre_precessor(
                args.sail, datasail_temp, random_names_dict, datasail_output, args.threads)

            print("Running datasail...")

            distance_matrix = os.path.join(
                datasail_output, "distance_matrix.tsv")

            datasail_runner(distance_matrix, datasail_output,
                            splits=train_test, threads=args.threads)

            if not args.keep_temp_files:
                print(f"Removing temp folder {datasail_temp}...")
                temp_folder_remover(datasail_temp)

            if not os.path.exists(os.path.join(datasail_output, "splits.tsv")):
                print(
                    "Error: Splits file does not exist. Check the datasail output folder.")
                sys.exit(1)

        with open(os.path.join(datasail_output, "splits.tsv")) as splits_file:
            lines = splits_file.readlines()
            for line in lines:
                splitted = line.split("\t")
                if splitted[1].strip() == "train":
                    train_strains.append(splitted[0].strip())
                elif splitted[1].strip() == "test":
                    test_strains.append(splitted[0].strip())

    if args.ml_algorithm == "rf":

        if args.parameter_optimization:

            same_setup_run_count = 1

            while True:

                if same_setup_run_count == 99:
                    print("Error: Same setup run count reached 99.")
                    print("Please change the output folder name.")
                    sys.exit(1)
                if args.sail:
                    output_folder_name = f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_RF_AutoML_{same_setup_run_count}_sail"
                else:
                    output_folder_name = f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_RF_AutoML_{same_setup_run_count}"
                if not os.path.exists(os.path.join(ml_output, output_folder_name)):
                    os.mkdir(os.path.join(ml_output, output_folder_name))
                    ml_output = os.path.join(ml_output, output_folder_name)
                    break
                else:
                    same_setup_run_count += 1

            with open(os.path.join(ml_output, "log_file.txt"), "w") as log_file:
                with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
                    rf_auto_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.optimization_time_limit,
                               args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, train=train_strains, test=test_strains, same_setup_run_count=same_setup_run_count)

        else:
            with open(os.path.join(ml_output, "log_file.txt"), "w") as log_file:
                with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
                    rf(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy,
                       custom_scorer="MCC", fia_repeats=5, n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf, min_samples_split=args.min_samples_split, train=train_strains, test=test_strains)

    elif args.ml_algorithm == "svm":

        ml_log_name = f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_SVM"

        if args.resampling_strategy == "cv":
            with open(os.path.join(ml_output, f"{ml_log_name}_log_file.txt"), "w") as log_file:
                with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
                    svm_cv(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.test_train_split, ml_output, args.threads,
                           args.feature_importance_analysis, args.save_model, resampling_strategy="cv", fia_repeats=5, optimization=False, train=train_strains, test=test_strains)
        elif args.resampling_strategy == "holdout":
            with open(os.path.join(ml_output, f"{ml_log_name}_log_file.txt"), "w") as log_file:
                with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
                    svm(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.test_train_split, ml_output, args.threads,
                        args.feature_importance_analysis, args.save_model, resampling_strategy="holdout", fia_repeats=5, optimization=False, train=train_strains, test=test_strains)

    elif args.ml_algorithm == "gb":

        if args.parameter_optimization:

            same_setup_run_count = 1

            while True:

                if same_setup_run_count == 99:
                    print("Error: Same setup run count reached 99.")
                    print("Please change the output folder name.")
                    sys.exit(1)

                if not os.path.exists(os.path.join(ml_output, f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_GB_AutoML_{same_setup_run_count}")):
                    os.mkdir(os.path.join(
                        ml_output, f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_GB_AutoML_{same_setup_run_count}"))
                    ml_output = os.path.join(
                        ml_output, f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_GB_AutoML_{same_setup_run_count}")
                    break
                else:
                    same_setup_run_count += 1
            with open(os.path.join(ml_output, "log_file.txt"), "w") as log_file:
                with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
                    gb_auto_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.optimization_time_limit,
                               args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, train=train_strains, test=test_strains)
        else:
            with open(os.path.join(ml_output, "log_file.txt"), "w") as log_file:
                with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
                    gb(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy,
                       custom_scorer="MCC", fia_repeats=5, n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf, min_samples_split=args.min_samples_split, train=train_strains, test=test_strains)

    if not args.keep_temp_files:
        print(f"Removing temp folder {ml_temp}...")
        temp_folder_remover(ml_temp)

    print("Done")

    end_time = time.time()

    print(time_function(start_time, end_time))


def binary_table_threshold(args):

    start_time = time.time()

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

    binary_table_threshold_output = os.path.join(
        args.output, "binary_table_threshold")

    # Check if output folder empty
    if os.path.exists(binary_table_threshold_output) and os.path.isdir(binary_table_threshold_output):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)

    # Create the output folder
    if not os.path.exists(binary_table_threshold_output):
        os.mkdir(binary_table_threshold_output)

    threshold_percentage_float = float(args.threshold_percentage)

    print("Creating binary table with threshold...")

    binary_table_threshold_with_percentage(
        args.input, binary_table_threshold_output, threshold_percentage_float)

    print(f"Binary table with threshold {threshold_percentage_float} is created. Can be found in: {binary_table_threshold_output}/binary_mutation_table_threshold_{threshold_percentage_float}_percent.tsv")

    with open(args.input) as infile:
        line = infile.readline()
        splitted = line.split("\t")
        print(
            f"Number of mutations in the original binary table: {len(splitted) - 1}")
        original_table_mutations = len(splitted) - 1

    with open(os.path.join(binary_table_threshold_output, f"binary_mutation_table_threshold_{threshold_percentage_float}_percent.tsv")) as infile:
        line = infile.readline()
        splitted = line.split("\t")
        print(
            f"Number of mutations in the thresholded binary table: {len(splitted) - 1}")
        thresholded_table_mutations = len(splitted) - 1

    print(
    f"Percentage of mutations kept: {thresholded_table_mutations/original_table_mutations * 100:.2f}%")

    print("Done")

    end_time = time.time()

    print(time_function(start_time, end_time))


def phenotype_table_pipeline(args):

    start_time = time.time()

    # Sanity checks

    # Check the arguments
    if args.input is None:
        print("Error: Input folder path is required.")
        sys.exit(1)

    # Check if the input is a folder or a file
    if os.path.isdir(args.input):
        input_folder = args.input
    else:
        print("Error: Input should be a folder.")
        sys.exit(1)

    phenotype_output = os.path.join(args.output, "phenotype_table")

    # Check if output folder empty
    if os.path.exists(phenotype_output) and os.path.isdir(phenotype_output):
        if len(os.listdir(phenotype_output)) == 0 and not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)

    # Check if output folder exists and create if not
    if not os.path.exists(phenotype_output):
        os.mkdir(phenotype_output)

    accepted_fasta_file_extensions = [".fna", ".fasta", ".faa"]

    random_names_will_be_used = False

    if args.random_names_dict is not None:

        random_names_will_be_used = True
        print("Random names will be used...")

        random_names = {}
        with open(args.random_names_dict, "r") as infile:
            lines = infile.readlines()
            for line in lines:
                splitted = line.split("\t")
                random_names[splitted[0].strip()] = splitted[1].strip()

    if input_folder is not None:
        antibiotics = os.listdir(input_folder)
        for antibiotic in antibiotics:
            if antibiotic.startswith("."):
                continue
            antibiotic_path = os.path.join(input_folder, antibiotic)
            status = os.listdir(antibiotic_path)
            if not 'Resistant' in status:
                print(
                    f"Error: {antibiotic} folder does not contain resistant folder.")
                sys.exit(1)
            if not 'Susceptible' in status:
                print(
                    f"Error: {antibiotic} folder does not contain susceptible folder.")
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

            with open(os.path.join(phenotype_output, "strains.txt"), "a") as outfile:
                for strain in resistant_strains:
                    # Make sure path is same in both Windows and Linux
                    strain_path = os.path.join(resistant_path, strain)
                    strain_path = strain_path.replace("\\", "/")
                    outfile.write(f"{strain_path}\n")
                for strain in susceptible_strains:
                    # Make sure path is same in both Windows and Linux
                    strain_path = os.path.join(susceptible_path, strain)
                    strain_path = strain_path.replace("\\", "/")
                    outfile.write(f"{strain_path}\n")

        input_file = os.path.join(phenotype_output, "strains.txt")

    if input_file is not None:
        if not random_names_will_be_used:
            random_names = random_name_giver(
                input_file, os.path.join(phenotype_output, "random_names.txt"))

    print("Creating phenotype dataframe...")
    # Create the phenotype dataframe
    phenotype_dataframe_creator(args.input, os.path.join(
        phenotype_output, "phenotype_table.tsv"), random_names)

    print("Done")

    end_time = time.time()

    print(time_function(start_time, end_time))


def phylogenetic_tree_pipeline(args):

    start_time = time.time()

    # Sanity checks

    # Check the arguments
    if args.input is None:
        print("Error: Input file is required.")
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

    # Check if temp folder exists and create if not
    if not os.path.exists(args.temp):
        os.mkdir(args.temp)

    # Check if temp folder empty
    if os.path.exists(os.path.join(args.temp, "phylogeny")) and os.path.isdir(os.path.join(args.temp, "phylogeny")):
        if not args.overwrite:
            print("Error: Temp folder is not empty.")
            sys.exit(1)
    else:
        os.mkdir(os.path.join(args.temp, "phylogeny"))

    mash_temp = os.path.join(args.temp, "phylogeny")
    mash_output = os.path.join(args.output, "phylogeny")

    if not os.path.exists(mash_output):
        os.mkdir(mash_output)

    print("Running phylogenetic tree pipeline...")
    mash_preprocessor(args.input, mash_output,
                      mash_temp, args.random_names_dict)

    mash_distance_runner(mash_output, mash_temp)

    print("Done")

    if not args.keep_temp_files:
        print("Removing temp folder...")
        temp_folder_remover(mash_temp)

    end_time = time.time()

    print(time_function(start_time, end_time))


if __name__ == "__main__":
    main()
