#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import os
import sys
import argparse
import pathlib
import contextlib
import time
import multiprocessing
import shutil

from sr_amr.utils import is_tool_installed, temp_folder_remover, time_function, copy_and_zip_file
from sr_amr.version import __version__


from sr_amr.panacota import panacota_pre_processor, panacota_post_processor, panacota_pipeline_runner
from sr_amr.gwas import pyseer_runner, pyseer_similarity_matrix_creator, pyseer_phenotype_file_creator, pyseer_genotype_matrix_creator, pyseer_post_processor, pyseer_gwas_graph_creator, decision_tree_input_creator
from sr_amr.binary_tables import snippy_runner, prokka_runner, random_name_giver, panaroo_input_creator, panaroo_runner, binary_table_creator, binary_mutation_table_gpa_information_adder, binary_mutation_table_gpa_information_adder_panaroo, phenotype_dataframe_creator, phenotype_dataframe_creator_post_processor, prokka_create_database, snippy_processed_file_creator, annotation_file_from_snippy, cdhit_preprocessor, cdhit_runner, gene_presence_absence_file_creator
from sr_amr.binary_table_threshold import binary_table_threshold_with_percentage
from sr_amr.phylogeny_tree import mash_preprocessor, mash_distance_runner
from sr_amr.prps import PRPS_runner
from sr_amr.ds import datasail_runner, datasail_pre_precessor
from sr_amr.ml import prps_ml_preprecessor, combined_ml
from sr_amr.full_automatix import automatix_runner
from sr_amr.ml_common_files import fia_file_annotation
from sr_amr.structman import structman_input_creator, annotation_function
from sr_amr.prediction import process_data_for_prediction, predict, equalize_columns

import warnings

warnings.filterwarnings("ignore")


def main():
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Automated Learning Pipeline for Antimicrobial Resistance (ALPAR)")

    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + __version__)

    subparsers = parser.add_subparsers(
        help='For suggested pipeline, check out our github page: https://github.com/kalininalab/ALPAR')
    
    ##############################################################################################################
    # Fully automated pipeline

    parser_automatix = subparsers.add_parser(
        'automatix', help='run automated pipeline')
    parser_automatix.add_argument(
        '-i', '--input', type=str, help='input folder path (check folder structure)', required=True)
    parser_automatix.add_argument(
        '-o', '--output', type=str, help='path of the output folder', required=True)
    parser_automatix.add_argument(
        '--reference', type=str, help='path of the reference file', required=True)
    parser_automatix.add_argument('--custom_database', type=str,
                                  help='creates and uses custom database for prokka, require path of the fasta file, default=None')
    parser_automatix.add_argument('--just_mutations', action='store_true',
                                  help='only creates binary mutation table with mutations, without gene presence absence information, default=False')
    parser_automatix.add_argument('--no_feature_importance_analysis', action='store_true',help='do not run feature importance analysis on ML step, default=False')
    parser_automatix.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_automatix.add_argument(
        '--threads', type=int, help='number of threads to use, default=1', default=1)
    parser_automatix.add_argument(
        '--ram', type=int, help='amount of ram to use in GB, default=8', default=8)
    parser_automatix.add_argument(
        '--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_automatix.add_argument('--overwrite', action='store_true',
                                  help='overwrite the output and temp folder if exists, default=False')
    parser_automatix.add_argument('--ml_algorithm', nargs='+',
                              help='classification algorithm to be used, available selections: [rf, svm, gb], default=[rf, svm, gb]', default=["rf", "svm", "gb"])
    parser_automatix.add_argument('--no_ml', action='store_true', help='do not run machine learning analysis, default=False')
    parser_automatix.add_argument('--fast', action='store_true', help='fast mode, does not run PanACoTA pipeline for phylogenetic tree analysis, default=False')
    parser_automatix.add_argument('--checkpoint', action='store_true',
                                  help='continues run from the checkpoint, default=False')
    parser_automatix.add_argument('--use_panaroo', action='store_true',help='use panaroo for gene presence absence analysis, WARNING: REQUIRES A LOT OF MEMORY, default=False')
    parser_automatix.add_argument('--no_datasail', action='store_true', help='splits data randomly instead of using genomic distances, default=False')
    parser_automatix.add_argument('--verbosity', type=int,
                                  help='verbosity level, default=1', default=1)
    parser_automatix.set_defaults(func=fully_automated_pipeline)

    ##############################################################################################################
    # Create binary tables

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
    parser_main_pipeline.add_argument(
        '--use_panaroo', action='store_true', help='use panaroo for gene presence absence analysis, WARNING: REQUIRES A LOT OF MEMORY, default=False')
    parser_main_pipeline.add_argument('--checkpoint', action='store_true',
                                      help='continues run from the checkpoint, default=False')
    parser_main_pipeline.add_argument('--verbosity', type=int,
                                      help='verbosity level, default=1', default=1)
    parser_main_pipeline.set_defaults(func=binary_table_pipeline)

    ##############################################################################################################
    # Create phenotype table

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
    parser_phenotype_table.add_argument('--verbosity', type=int,
                                        help='verbosity level, default=1', default=1)
    parser_phenotype_table.set_defaults(func=phenotype_table_pipeline)

    ##############################################################################################################
    # Create thresholded binary tables

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
    parser_binary_tables_threshold.add_argument('--verbosity', type=int,
                                                help='verbosity level, default=1', default=1)
    parser_binary_tables_threshold.set_defaults(func=binary_table_threshold)

    ##############################################################################################################
    # Create phylogenetic tree with PanACoTA pipeline

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
    parser_panacota.add_argument('--verbosity', type=int,
                                 help='verbosity level, default=1', default=1)
    parser_panacota.set_defaults(func=panacota_pipeline)

    ##############################################################################################################
    # Create phylogenetic tree with Mashtree

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
    parser_phylogenetic_tree.add_argument('--verbosity', type=int,
                                          help='verbosity level, default=1', default=1)
    parser_phylogenetic_tree.set_defaults(func=phylogenetic_tree_pipeline)

    ##############################################################################################################
    # Run GWAS analysis

    parser_gwas = subparsers.add_parser('gwas', help='run gwas analysis')
    parser_gwas.add_argument('-i', '--input', type=str,
                             help='binary mutation table path', required=True)
    parser_gwas.add_argument('-o', '--output', type=str,
                             help='path of the output folder', required=True)
    parser_gwas.add_argument(
        '-p', '--phenotype', type=str, help='phenotype table path', required=True)
    parser_gwas.add_argument('-t', '--tree', type=str,
                             help='phylogenetic tree path', required=True)
    parser_gwas.add_argument('--threads', type=int,
                             help='number of threads to use, default=1', default=1)
    parser_gwas.add_argument('--overwrite', action='store_true',
                             help='overwrite the output folder if exists, default=False')
    parser_gwas.add_argument('--verbosity', type=int,
                             help='verbosity level, default=1', default=1)
    parser_gwas.set_defaults(func=gwas_pipeline)

    ##############################################################################################################
    # Run PRPS analysis

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
    parser_prps.add_argument('--verbosity', type=int,
                             help='verbosity level, default=1', default=1)
    parser_prps.set_defaults(func=prps_pipeline)

    ##############################################################################################################
    # Run Machine Learning analysis

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
    parser_ml.add_argument('--annotation', type=str, help='annotation file path')
    parser_ml.add_argument('--prps', type=str, help='prps score file path')
    parser_ml.add_argument('--prps_percentage', type=int,
                           help='percentage of the top scores of prps to be used, should be used with --prps option, default=30', default=30)
    parser_ml.add_argument('--sail', type=str, help='split against information leakage, requires txt file that contains path of each strain per line or input folder path, can be found create_binary_tables output path as strains.txt, default=None')
    parser_ml.add_argument('--sail_epsilon', type=str, help='epsilon value for datasail, default=0.1',
                           default=0.1)
    parser_ml.add_argument('--sail_delta', type=str, help='delta value for datasail, default=0.1')
    parser_ml.add_argument('--sail_solver', type=str, help='solver for datasail, available selections: [SCIP, MOSEK, GUROBI, CPLEX], check https://datasail.readthedocs.io/en/latest/workflow/solvers.html for more information default=SCIP')
    parser_ml.add_argument('--train_strains_file', type=str,
                           help='train strains file path', default=None)
    parser_ml.add_argument('--test_strains_file', type=str,
                           help='test strains file path', default=None)
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
                           help='max depth for random forest, default=None', default=None)
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
    parser_ml.add_argument('--important_feature_limit', type=int,
                           help='number of reported maximum number of features in FIA file, default=10', default=10)
    parser_ml.add_argument('--feature_importance_analysis_number_of_repeats', type=int,
                           help='number of repeats for feature importance analysis should be given with --feature_importance_analysis option, default=5', default=5)
    parser_ml.add_argument('--feature_importance_analysis_strategy', type=str,
                           help='strategy for feature importance analysis, available selections: [gini, permutation_importance], default=gini', default="gini")
    parser_ml.add_argument('--optimization_time_limit', type=int,
                           help='time limit for parameter optimization with AutoML, default=3600', default=3600)
    parser_ml.add_argument('--svm_kernel', type=str,
                           help='kernel for svm, available selections: [linear, poly, rbf, sigmoid], default=linear', default="linear")
    parser_ml.add_argument('--no_stratify_split', type=str,
                           help='if given, does not uses stratify in random split', default="False")
    parser_ml.add_argument('--verbosity', type=int,
                           help='verbosity level, default=1', default=1)
    parser_ml.set_defaults(func=ml_pipeline)

    ##############################################################################################################
    # Create StructMAn input

    parser_structman = subparsers.add_parser(
        'structman', help='structman input creator from feature importance analysis results')
    parser_structman.add_argument(
        '-i', '--input', type=str, help='output file path or txt file that contains path of each feature importance analysis result per line', required=True)
    parser_structman.add_argument(
        '-o', '--output', type=str, help='path of the output folder', required=True)
    parser_structman.add_argument(
        '--annotation', type=str, help='path of the annotation file, can be found in binary tables output as "mutations_annotations.tsv"', required=True)
    parser_structman.add_argument(
        '--reference', type=str, help='path of the reference file', required=True)
    parser_structman.add_argument('--overwrite', action='store_true',
                                  help='overwrite the output folder if exists, default=False')
    parser_structman.add_argument(
        '--name', type=str, help='name of the structman analysis, default=WIBI', default="WIBI")
    parser_structman.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_structman.add_argument(
        '--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_structman.add_argument('--verbosity', type=int,
                                  help='verbosity level, default=1', default=1)
    parser_structman.set_defaults(func=structman_pipeline)

    ##############################################################################################################
    # Create prediction table & predict unknown strains

    parser_prediction = subparsers.add_parser(
        'prediction', help='run prediction analysis')
    parser_prediction.add_argument('-i', '--input', type=str,
                                   help='txt file that contains path of each strain per line or input folder path')
    parser_prediction.add_argument('-p', '--prediction_table', type=str, help='prediction table path')
    parser_prediction.add_argument('-b', '--model_binary_table', type=str, help='binary mutation table path that is used to create model', required=True)
    parser_prediction.add_argument('-o', '--output', type=str,
                                   help='path of the output folder', required=True)
    parser_prediction.add_argument('--model', type=str,
                                   help='path of the ml model', required=True)
    parser_prediction.add_argument(
        '--reference', type=str, help='path of the reference file')
    parser_prediction.add_argument('--custom_database', type=str, nargs=2,
                                    help='creates and uses custom database for prokka, require path of the fasta file and genus name, default=None')
    parser_prediction.add_argument('--no_gene_presence_absence', action='store_true',
                                    help='do not run gene presence absence functions, default=False')
    parser_prediction.add_argument('--prps', type=str, nargs=2,
                                    help='use prps scores, prps score table path and percentage should be given, percentage default value = 30, default = None, 30', default=[None, 30])
    parser_prediction.add_argument(
        '--no_gene_annotation', action='store_true', help='do not run gene annotation, default=False')
    parser_prediction.add_argument('--use_panaroo', action='store_true', help='use panaroo for gene presence absence analysis, WARNING: REQUIRES A LOT OF MEMORY, default=False')
    parser_prediction.add_argument('--threads', type=int,
                                   help='number of threads to use, default=1', default=1)
    parser_prediction.add_argument('--ram', type=int, help='amount of ram to use in GB, default=4', default=4)
    parser_prediction.add_argument('--temp', type=str,
                                   help='path of the temporary directory, default=output_folder/temp')
    parser_prediction.add_argument(
        '--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_prediction.add_argument('--overwrite', action='store_true',
                                   help='overwrite the output folder if exists, default=False')
    parser_prediction.add_argument('--verbosity', type=int,
                                   help='verbosity level, default=1', default=1)
    
    parser_prediction.set_defaults(func=prediction_pipeline)

    # Parse the arguments
    args = parser.parse_args()

    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


def run_snippy_and_prokka(strain, random_names, snippy_output, prokka_output, args, snippy_flag, prokka_flag, custom_db=None):
    # if args.ram / args.threads < 8:
    #     print("Warning: Not enough ram for the processes. Minimum 8 GB of ram per thread is recommended.")
    
    # Snippy creates issue with high memory usage, so it is limited to 100 GB
    if args.ram > 100:
        #print(f"Warning: Due to restrictions of snippy, ram is limited to 100 GB for snippy run only.")
        args.ram = 100

    if snippy_flag:
        snippy_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[
                      0]], snippy_output, args.reference, f"{args.temp}/snippy_log.txt", 1, args.ram)
    if prokka_flag:
        prokka_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[
                      0]], prokka_output, args.reference, f"{args.temp}/prokka_log.txt", 1, custom_db)


def binary_table_pipeline(args):

    start_time = time.time()

    # Sanity checks

    status = -1

    if args.verbosity > 3:
        print("Checking the installed tools...")
    # Check if the tools are installed
    tool_list = ["snippy", "prokka", "panaroo"]

    for tool in tool_list:
        if not is_tool_installed(tool):
            print(f"Error: {tool} is not installed.")
            sys.exit(1)

    if args.verbosity > 3:
        print("Checking the input...")
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
            print("Accepted extensions: .gbk, .gbff")
            sys.exit(1)

    if args.custom_database:
        if len(args.custom_database) != 2:
            print("Error: Custom database option should have two arguments.")
            sys.exit(1)

        if not os.path.exists(args.custom_database[0]):
            print("Error: Custom database fasta file does not exist.")
            sys.exit(1)

        if pathlib.Path(args.custom_database[0]).suffix != ".fasta":
            print("Error: Custom database file extension is not accepted.")
            print("Accepted extension: .fasta")
            sys.exit(1)

    # Check if output folder empty
    if os.path.exists(args.output) and os.path.isdir(args.output):
        if len(os.listdir(args.output)) > 0 and args.checkpoint:
            print("Warning: Output folder is not empty.")
            print("Checking if there is previous run files to continue.")
            if args.temp is None:
                print(
                    "Warning: Temp folder is not given. Checking output folder for temp folder.")
                temp_folder = os.path.join(args.output, "temp")
            else:
                temp_folder = args.temp
            if os.path.exists(temp_folder):
                if os.path.exists(os.path.join(temp_folder, "status.txt")):
                    with open(os.path.join(temp_folder, "status.txt"), "r") as infile:
                        line = infile.readline()
                        status = int(line.strip())
                        print(f"Previous run found. Continuing from step {status}")
                        with open(os.path.join(args.output, "random_names.txt")) as random_names_file:
                            random_names = {}
                            for random_names_file_line in random_names_file.readlines():
                                random_names[random_names_file_line.split(
                                    "\t")[0].strip()] = random_names_file_line.split("\t")[1].strip()
                        
                        if args.verbosity > 4:
                            print(f"Length of random names: {len(random_names)}")

        elif len(os.listdir(args.output)) > 0 and not args.overwrite:
            print("Error: Output folder is not empty.")
            print("If you want to overwrite the output folder, use --overwrite option.")
            print("If you want to continue previous run, use --checkpoint option.")
            sys.exit(1)

    # Check if output folder exists and create if not
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    temp_folder_created = False
    continue_flag = False

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
        if not os.path.exists(args.temp):
            os.mkdir(args.temp)
            temp_folder_created = True
        else:
            print("Warning: Temp folder already exists. Will be used for the run.")
            temp_folder_created = True

    else:
        if not os.path.exists(args.temp):
            os.mkdir(args.temp)
            temp_folder_created = True
        else:
            print("Warning: Temp folder already exists. Will be used for the run.")
            temp_folder_created = True

    if not temp_folder_created:
        # Check if temp folder empty
        if os.path.exists(args.temp) and os.path.isdir(args.temp):
            if os.listdir(args.temp) and not args.overwrite:
                print(
                    "Error: Temp folder is not empty. Please remove the temp folder or use --overwrite option.")
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
    panaroo_flag = False
    gpa_flag = True

    if args.no_gene_presence_absence:
        panaroo_flag = False
        gpa_flag = False
    if args.no_gene_annotation:
        prokka_flag = False
    if args.use_panaroo:
        panaroo_flag = True

    if args.verbosity > 3:
        print("Will run snippy: ", snippy_flag)
        print("Will run prokka: ", prokka_flag)
        print("Will run panaroo: ", panaroo_flag)
        print("Will run gene presence absence: ", gpa_flag)

    # Create the output folder
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(os.path.join(args.output, "snippy")):
        os.mkdir(os.path.join(args.output, "snippy"))
    if not os.path.exists(os.path.join(args.output, "prokka")):
        os.mkdir(os.path.join(args.output, "prokka"))
    if not os.path.exists(os.path.join(args.output, "panaroo")):
        os.mkdir(os.path.join(args.output, "panaroo"))
    if not os.path.exists(os.path.join(args.output, "cd-hit")):
        os.mkdir(os.path.join(args.output, "cd-hit"))

    snippy_output = os.path.join(args.output, "snippy")
    prokka_output = os.path.join(args.output, "prokka")
    panaroo_output = os.path.join(args.output, "panaroo")

    if args.verbosity > 3:
        print(f"Output folder created: {args.output}")
        print(f"Snippy output folder created: {snippy_output}")
        print(f"Prokka output folder created: {prokka_output}")
        print(f"Panaroo output folder created: {panaroo_output}")

    # Create the temp folder for the panaroo input
    if not os.path.exists(os.path.join(args.temp, "panaroo")):
        os.mkdir(os.path.join(args.temp, "panaroo"))
    # Create the temp folder for the cdhit input
    if not os.path.exists(os.path.join(args.temp, "cdhit")):
        os.mkdir(os.path.join(args.temp, "cdhit"))

    if status < 0:

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
                # Checking for Mac OS hidden files
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

                if args.verbosity > 3:
                    print(f"Checking {antibiotic} folder...")
                    print(f"Resistant folder: {resistant_path}")
                    print(f"Amount of files in resistant folder: {len(files_in_resistant_path)}")
                    print(f"Susceptible folder: {susceptible_path}")
                    print(f"Amount of files in susceptible folder: {len(files_in_susceptible_path)}")

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
                        strain_path = os.path.abspath(strain_path)
                        strain_path = strain_path.replace("\\", "/")
                        outfile.write(f"{strain_path}\n")
                    for strain in susceptible_strains:
                        strain_path = os.path.join(susceptible_path, strain)
                        strain_path = os.path.abspath(strain_path)
                        strain_path = strain_path.replace("\\", "/")
                        outfile.write(f"{strain_path}\n")

            input_file = os.path.join(args.output, "strains.txt")

            if args.verbosity > 3:
                print(f"Strains.txt file created: {input_file}")

        if input_file is not None:
            random_names = random_name_giver(
                input_file, os.path.join(args.output, "random_names.txt"))
            if args.verbosity > 3:
                print(f"Random names file created: {os.path.join(args.output, 'random_names.txt')}")

        with open(os.path.join(args.temp, "status.txt"), "w") as outfile:
            outfile.write(f"0")
            status = 0

    with open(os.path.join(args.output, "random_names.txt")) as random_names_file:
        random_names = {}
        for random_names_file_line in random_names_file.readlines():
            random_names[random_names_file_line.split(
                "\t")[0].strip()] = random_names_file_line.split("\t")[1].strip()

    input_file = os.path.join(args.output, "strains.txt")

    if status < 1:

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

            with open(os.path.join(args.temp, "status.txt"), "w") as outfile:
                outfile.write(f"1")
                status = 1

    if status < 2:
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

        # We will use a status file to indicate checkpoints

        with open(os.path.join(args.temp, "status.txt"), "w") as outfile:
            outfile.write(f"2")
            status = 2

    if status < 3:

        strains_to_be_processed = []

        prokka_output_strains = os.listdir(prokka_output)
        snippy_output_strains = os.listdir(snippy_output)

        strains_to_be_skiped = []

        for strain in random_names.keys():
            if prokka_flag and snippy_flag:
                if random_names[strain] in prokka_output_strains and random_names[strain] in snippy_output_strains:
                    strains_to_be_processed.append(random_names[strain])
                else:
                    strains_to_be_skiped.append(random_names[strain])
            elif prokka_flag and not snippy_flag:
                if random_names[strain] in prokka_output_strains:
                    strains_to_be_processed.append(random_names[strain])
                else:
                    strains_to_be_skiped.append(random_names[strain])
            elif snippy_flag and not prokka_flag:
                if random_names[strain] in snippy_output_strains:
                    strains_to_be_processed.append(random_names[strain])
                else:
                    strains_to_be_skiped.append(random_names[strain])

        print(f"Number of strains processed: {len(strains_to_be_processed)}")
        print(f"Number of strains skipped: {len(strains_to_be_skiped)}")

        snippy_processed_file_creator(snippy_output, os.path.join(
            args.output, "snippy_processed_strains.txt"))

        print("Creating binary mutation table...")
        # Create the binary table
        binary_table_creator(snippy_output, os.path.join(
            args.output, "binary_mutation_table.tsv"), args.threads, strains_to_be_processed)

        print("Creating annotation table...")

        annotation_file_from_snippy(snippy_output, args.output)

        with open(os.path.join(args.temp, "status.txt"), "w") as outfile:
            outfile.write(f"3")
            status = 3

    if status < 4:

        if panaroo_flag:

            try:

                print("Creating panaroo input...")
                # Create the panaroo input
                panaroo_input_creator(os.path.join(args.output, "random_names.txt"), prokka_output, os.path.join(
                    args.temp, "panaroo"), strains_to_be_processed)

                print("Running panaroo...")
                # Run panaroo
                panaroo_runner(os.path.join(args.temp, "panaroo"), panaroo_output, os.path.join(
                    args.temp, "panaroo_log.txt"), args.threads)

                print("Adding gene presence absence information to the binary table...")
                # Add gene presence absence information to the binary table

                if not os.path.exists(os.path.join(panaroo_output, "gene_presence_absence.csv")):
                    print("Warning: Gene presence absence file does not exist.")
                    print(
                        "Gene presence absence information will not be added to the binary table.")
                    do_not_remove_temp = True

                else:
                    binary_mutation_table_gpa_information_adder_panaroo(os.path.join(args.output, "binary_mutation_table.tsv"), os.path.join(
                        panaroo_output, "gene_presence_absence.csv"), os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv"))
                    do_not_remove_temp = False

            except Exception as e:
                print("Error: Panaroo could not be run.")
                print(e)
                do_not_remove_temp = True

        if gpa_flag:

            try: 

                print("Creating gene presence absence input...")
                # Create the gene presence absence information
                print(f"CD-HIT preprecessor is running...")
                cdhit_preprocessor(os.path.join(args.output, "random_names.txt"), prokka_output, os.path.join(args.temp, "cdhit"), strains_to_be_processed)

                shutil.copy(os.path.join(args.temp, 'cdhit', 'protein_positions.csv'), os.path.join(args.output, 'cd-hit', 'protein_positions.csv'))

                print(f"CD-HIT is running...")
                cdhit_runner(os.path.join(args.temp, "cdhit", "combined_proteins.faa"), os.path.join(args.output, "cd-hit", "cdhit_output.txt"), n_cpu=args.threads)

                print(f"Gene presence-absence matrix is being created...")
                gene_presence_absence_file_creator(os.path.join(args.output, "cd-hit", "cdhit_output.txt.clstr"), strains_to_be_processed, os.path.join(args.temp, "cdhit"))

                print("Adding gene presence absence information to the binary table...")
                # Add gene presence absence information to the binary table

                if not os.path.exists(os.path.join(args.temp, "cdhit", "gene_presence_absence_matrix.csv")):
                    print("Warning: Gene presence absence file does not exist.")
                    print(
                        "Gene presence absence information will not be added to the binary table.")
                    do_not_remove_temp = True
                
                else:
                    binary_mutation_table_gpa_information_adder(os.path.join(args.output, "binary_mutation_table.tsv"), os.path.join(args.temp, "cdhit", "gene_presence_absence_matrix.csv"), os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv"))
                    do_not_remove_temp = False

            except Exception as e:
                print("Error: CD-HIT could not be run.")
                print(e)
                do_not_remove_temp = True

        if args.create_phenotype_from_folder:
            print("Creating phenotype dataframe...")
            # Create the phenotype dataframe
            phenotype_dataframe_creator(args.create_phenotype_from_folder, os.path.join(
                args.output, "phenotype_table.tsv"), random_names)

            phenotype_dataframe_creator_post_processor(os.path.join(
                args.output, "binary_mutation_table.tsv"), os.path.join(args.output, "phenotype_table.tsv"))

        with open(os.path.join(args.temp, "status.txt"), "w") as outfile:
            outfile.write(f"4")
            status = 4

    if args.keep_temp_files:
        print("Warning, temp files will be kept this might take up space.")

    if not args.keep_temp_files:
        if not do_not_remove_temp:
            print("Removing temp folder...")
            temp_folder_remover(os.path.join(args.temp))
            temp_folder_remover(os.path.join(args.output, "snippy"))
            temp_folder_remover(os.path.join(args.output, "prokka"))
            temp_folder_remover(os.path.join(args.output, "panaroo"))

        else:
            print(
                "Warning: Temp folder will not be removed because gene presence absence file does not exist.")
            print(
                "Temp folder can be used for re-run panaroo again for gene presence absence information.")
            print("You can remove the temp folder manually. Temp folder path: ",
                  os.path.join(args.temp))

    print(f"Binary tables are created, can be found in {args.output}")

    end_time = time.time()

    print(time_function(start_time, end_time))


def panacota_pipeline(args):

    start_time = time.time()

    tool_list = ["PanACoTA"]

    for tool in tool_list:
        if not is_tool_installed(tool):
            print(f"Error: {tool} is not installed.")
            print("Please install the tool via `pip install panacota` and try again.")
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
            if len(os.listdir(panacota_temp)) > 0 and not args.overwrite:
                print("Error: Temp folder is not empty.")
                sys.exit(1)

    if not output_folder_created:
        # Check if output folder empty
        if os.path.exists(panacota_output) and os.path.isdir(panacota_output):
            if len(os.listdir(panacota_output)) > 0 and not args.overwrite:
                print("Error: Output folder is not empty.")
                sys.exit(1)

    panacota_log_file = os.path.join(panacota_output, "panacota_log.txt")

    try:

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

        print(
            f"PanACoTA pipeline is finished, results can be found in the {panacota_output}")
        
    except Exception as e:
        print("Error: PanACoTA pipeline could not be run.")
        print(e)

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
        if len(os.listdir(gwas_output)) > 0 and not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)
        else:
            if args.overwrite:
                print("Warning: Output folder is not empty. Old files will be deleted.")
                temp_folder_remover(gwas_output)
                os.makedirs(gwas_output, exist_ok=True)
    else:
        os.mkdir(gwas_output)

    if not os.path.exists(os.path.join(gwas_output, "gwas_results")):
        os.mkdir(os.path.join(gwas_output, "gwas_results"))

    if not os.path.exists(os.path.join(gwas_output, "pyseer_phenotypes")):
        os.mkdir(os.path.join(gwas_output, "pyseer_phenotypes"))

    if not os.path.exists(os.path.join(gwas_output, "graphs")):
        os.mkdir(os.path.join(gwas_output, "graphs"))

    try:

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
                    os.path.join(gwas_output, "similarity_matrix.tsv"), os.path.join(gwas_output, "gwas_results"), args.threads)

        if not os.path.exists(os.path.join(gwas_output, "sorted")):
            os.mkdir(os.path.join(gwas_output, "sorted"))

        if not os.path.exists(os.path.join(gwas_output, "sorted_cleaned")):
            os.mkdir(os.path.join(gwas_output, "sorted_cleaned"))

        pyseer_post_processor(os.path.join(
            gwas_output, "gwas_results"), gwas_output)

        pyseer_gwas_graph_creator(gwas_output, os.path.join(gwas_output, "graphs"))

        #def decision_tree_input_creator(binary_table, phenotype_file_path, pyseer_output_folder, output_folder):
        decision_tree_input_creator(args.input, args.phenotype, gwas_output, gwas_output)
    
    except Exception as e:
        print("Error: GWAS pipeline could not be run.")
        print(e)

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
        if len(os.listdir(prps_output)) > 0 and not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)
        else:
            if args.overwrite:
                print("Warning: Output folder is not empty. Old files will be deleted.")
                temp_folder_remover(prps_output)
                os.makedirs(prps_output, exist_ok=True)
    else:
        os.makedirs(prps_output, exist_ok=True)

    if not os.path.exists(prps_temp):
        os.mkdir(prps_temp)

    else:
        if args.overwrite:
            print("Warning: Temp folder is not empty. Old files will be deleted.")
            temp_folder_remover(prps_temp)
            os.makedirs(prps_temp, exist_ok=True)

    try:
        print("Running PRPS...")

        PRPS_runner(args.tree, args.input, prps_output, prps_temp)

        if not args.keep_temp_files:
            print("Removing temp folder...")
            temp_folder_remover(prps_temp)

        print(
            f"PRPS is finished, results can be found in the {os.path.abspath(prps_output)}")
    
    except Exception as e:
        print("Error: PRPS could not be run.")
        print(e)

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
    svm_output = os.path.join(ml_output, "svm")
    rf_output = os.path.join(ml_output, "rf")
    gb_output = os.path.join(ml_output, "gb")
    hist_gb_output = os.path.join(ml_output, "hist_gb")

    # Check if output folder empty
    if os.path.exists(ml_output) and os.path.isdir(ml_output):
        if len(os.listdir(ml_output)) > 0 and not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)
        else:
            if args.overwrite:
                print("Warning: Output folder is not empty. Old files will be deleted.")
                temp_folder_remover(ml_output)
                os.makedirs(ml_output, exist_ok=True)
    else:
        os.mkdir(ml_output)

    if os.path.exists(ml_temp) and os.path.isdir(ml_temp):
        if len(os.listdir(ml_temp)) > 0 and not args.overwrite:
            print("Error: Temp folder is not empty.")
            sys.exit(1)
        else:
            if args.overwrite:
                print("Warning: Temp folder is not empty. Old files will be deleted.")
                temp_folder_remover(ml_temp)
                os.makedirs(ml_temp, exist_ok=True)
    else:
        os.makedirs(ml_temp, exist_ok=True)

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

    binary_table_path_used_for_model = copy_and_zip_file(binary_mutation_table_path, ml_output, "model_binary_mutation_table")

    print(f"Binary mutation table that will be used for model can be found in the {os.path.join(ml_output,binary_table_path_used_for_model)}")

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

    if args.feature_importance_analysis_strategy not in ["gini", "permutation_importance"]:
        print("Error: Feature importance analysis strategy should be either gini or permutation_importance.")
        sys.exit(1)

    train_strains = []
    test_strains = []

    if args.sail:

        if args.train_strains_file or args.test_strains_file:
            print("Error: If you want to use DataSAIL as data split, train and test strains files should not be provided.")
            sys.exit(1)

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
            
            if not args.sail_epsilon:
                args.sail_epsilon = 0.1
            if not args.sail_delta:
                args.sail_delta = 0.1
            if not args.sail_solver:
                args.sail_solver = "SCIP"

            datasail_runner(distance_matrix, datasail_output,
                            splits=train_test, cpus=args.threads, epsilon=args.sail_epsilon, delta=args.sail_delta, solver=args.sail_solver)

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

    if args.train_strains_file and args.test_strains_file:
        if os.path.exists(args.train_strains_file):
            with open(args.train_strains_file) as train_file:
                train_strains_lines = train_file.readlines()
                for train_strain_line in train_strains_lines:
                    train_strains.append(train_strain_line.strip())
        else:
            print("Error: Train strains file does not exist.")
            sys.exit(1)

        if os.path.exists(args.test_strains_file):
            with open(args.test_strains_file) as test_file:
                test_strains_lines = test_file.readlines()
                for test_strain_line in test_strains_lines:
                    test_strains.append(test_strain_line.strip())
        else:
            print("Error: Test strains file does not exist.")
            sys.exit(1)

        if args.verbosity > 3:
            print(f"Train strains: {train_strains}")
            print(f"Test strains: {test_strains}")

    if args.train_strains_file and not args.test_strains_file:
        print("Error: Test strains file is missing.")
        sys.exit(1)

    if not args.train_strains_file and args.test_strains_file:
        print("Error: Train strains file is missing.")
        sys.exit(1)

    if args.no_stratify_split:
        stratiy_random_split = False
    else:
        stratiy_random_split = True

    if args.ml_algorithm == "rf":

        if os.path.exists(rf_output):
            if args.overwrite:
                print("Warning: Output folder is not empty. Old files will be deleted.")
                temp_folder_remover(rf_output)
                os.makedirs(rf_output, exist_ok=True)
            else:
                print("Error: Output folder is not empty.")
                print(
                    "If you want to overwrite the output folder, use --overwrite option.")
                sys.exit(1)
        else:
            os.makedirs(rf_output, exist_ok=True)

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

                    fia_file = combined_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.optimization_time_limit, "rf_auto_ml", args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, train=train_strains, test=test_strains, same_setup_run_count=same_setup_run_count, stratify=stratiy_random_split, feature_importance_analysis_strategy=args.feature_importance_analysis_strategy, important_feature_limit=args.important_feature_limit) 

        else:
            with open(os.path.join(ml_output, "log_file.txt"), "w") as log_file:
                with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):

                    fia_file = combined_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.optimization_time_limit, "rf", args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf, min_samples_split=args.min_samples_split, train=train_strains, test=test_strains, stratify=stratiy_random_split, feature_importance_analysis_strategy=args.feature_importance_analysis_strategy, important_feature_limit=args.important_feature_limit) 

    elif args.ml_algorithm == "svm":

        if os.path.exists(svm_output):
            if args.overwrite:
                print("Warning: Output folder is not empty. Old files will be deleted.")
                temp_folder_remover(svm_output)
                os.makedirs(svm_output, exist_ok=True)
            else:
                print("Error: Output folder is not empty.")
                print(
                    "If you want to overwrite the output folder, use --overwrite option.")
                sys.exit(1)
        else:
            os.makedirs(svm_output, exist_ok=True)

        ml_log_name = f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_SVM"

        with open(os.path.join(ml_output, f"{ml_log_name}_log_file.txt"), "w") as log_file:
            with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
                
                fia_file = combined_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.optimization_time_limit, "svm", args.feature_importance_analysis, args.save_model, resampling_strategy="cv", custom_scorer="MCC", fia_repeats=5, train=train_strains, test=test_strains, stratify=stratiy_random_split, feature_importance_analysis_strategy="permutation_importance", important_feature_limit=args.important_feature_limit)

    elif args.ml_algorithm == "gb":

        if os.path.exists(gb_output):
            if args.overwrite:
                print("Warning: Output folder is not empty. Old files will be deleted.")
                temp_folder_remover(gb_output)
                os.makedirs(gb_output, exist_ok=True)
            else:
                print("Error: Output folder is not empty.")
                print(
                    "If you want to overwrite the output folder, use --overwrite option.")
                sys.exit(1)
        else:
            os.makedirs(gb_output, exist_ok=True)

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

                    fia_file = combined_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.optimization_time_limit, "gb_auto_ml", args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, train=train_strains, test=test_strains, same_setup_run_count=same_setup_run_count, stratify=stratiy_random_split, feature_importance_analysis_strategy=args.feature_importance_analysis_strategy, important_feature_limit=args.important_feature_limit) 

        else:
            with open(os.path.join(ml_output, "log_file.txt"), "w") as log_file:
                with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):

                    fia_file = combined_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.optimization_time_limit, "gb", args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf, min_samples_split=args.min_samples_split, train=train_strains, test=test_strains, stratify=stratiy_random_split, feature_importance_analysis_strategy=args.feature_importance_analysis_strategy, important_feature_limit=args.important_feature_limit) 

    elif args.ml_algorithm == "histgb":

        if os.path.exists(hist_gb_output):
            if args.overwrite:
                print("Warning: Output folder is not empty. Old files will be deleted.")
                temp_folder_remover(hist_gb_output)
                os.makedirs(hist_gb_output, exist_ok=True)
            else:
                print("Error: Output folder is not empty.")
                print(
                    "If you want to overwrite the output folder, use --overwrite option.")
                sys.exit(1)
        else:
            os.makedirs(hist_gb_output, exist_ok=True)

        #TODO
        # Will be implemented later after tests
        # if args.parameter_optimization:

        #     same_setup_run_count = 1

        #     while True:

        #         if same_setup_run_count == 99:
        #             print("Error: Same setup run count reached 99.")
        #             print("Please change the output folder name.")
        #             sys.exit(1)

        #         if not os.path.exists(os.path.join(ml_output, f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_GB_AutoML_{same_setup_run_count}")):
        #             os.mkdir(os.path.join(
        #                 ml_output, f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_GB_AutoML_{same_setup_run_count}"))
        #             ml_output = os.path.join(
        #                 ml_output, f"seed_{args.random_state}_testsize_{args.test_train_split}_resampling_{args.resampling_strategy}_GB_AutoML_{same_setup_run_count}")
        #             break
        #         else:
        #             same_setup_run_count += 1

        #     with open(os.path.join(ml_output, "log_file.txt"), "w") as log_file:
        #         with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):

        #             fia_file = combined_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.optimization_time_limit, "gb_auto_ml", args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, train=train_strains, test=test_strains, same_setup_run_count=same_setup_run_count, stratify=stratiy_random_split, feature_importance_analysis_strategy=args.feature_importance_analysis_strategy, important_feature_limit=args.important_feature_limit) 

        # else:
        with open(os.path.join(ml_output, "log_file.txt"), "w") as log_file:
            with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file): 

                if args.min_samples_leaf == 1:
                    args.min_samples_leaf = 20

                fia_file = combined_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.optimization_time_limit, "histgb", args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf, min_samples_split=args.min_samples_split, train=train_strains, test=test_strains, stratify=stratiy_random_split, feature_importance_analysis_strategy=args.feature_importance_analysis_strategy, important_feature_limit=args.important_feature_limit) 


    if args.feature_importance_analysis:
        if args.annotation:
            if os.path.exists(fia_file):
                print(f"Annotating feature importance analysis file {fia_file}...")
                fia_file_annotation(fia_file, args.annotation)
                print(
                    f"Annotated feature importance analysis file is created, can be found in {fia_file}")
        elif os.path.exists(os.path.join(args.output, "mutations_annotations.tsv")):
            print(f"Annotating feature importance analysis file {fia_file}...")
            fia_file_annotation(fia_file, os.path.join(args.output, "mutations_annotations.tsv"))
            print(
                f"Annotated feature importance analysis file is created, can be found in {fia_file}")
                
        else:
            print("Warning: Annotations file does not exist.")
            print(
                "Feature importance analysis file will not be annotated.")

    if not args.keep_temp_files:
        print(f"Removing temp folder {ml_temp}...")
        temp_folder_remover(ml_temp)

    print(
        f"ML pipeline is finished, results can be found in the {os.path.abspath(ml_output)}")

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
        if len(os.listdir(binary_table_threshold_output)) > 0 and not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)

    # Create the output folder
    if not os.path.exists(binary_table_threshold_output):
        os.mkdir(binary_table_threshold_output)

    threshold_percentage_float = float(args.threshold_percentage)

    print("Creating binary table with threshold...")

    created_output_file = binary_table_threshold_with_percentage(
        args.input, binary_table_threshold_output, threshold_percentage_float)

    print(
        f"Binary table with threshold {threshold_percentage_float} is created. Can be found in: {created_output_file}")

    with open(args.input) as infile:
        line = infile.readline()
        splitted = line.split("\t")
        print(
            f"Number of mutations in the original binary table: {len(splitted) - 1}")
        original_table_mutations = len(splitted) - 1

    with open(created_output_file) as infile:
        line = infile.readline()
        splitted = line.split("\t")
        print(
            f"Number of mutations in the thresholded binary table: {len(splitted) - 1}")
        thresholded_table_mutations = len(splitted) - 1

    print(
        f"Percentage of mutations kept: {thresholded_table_mutations/original_table_mutations * 100:.2f}%")

    print(
        f"Thresholded binary table is created. Can be found in: {os.path.abspath(created_output_file)}")

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

    print("Phenotype dataframe is created. Can be found in: ", os.path.abspath(
        os.path.join(phenotype_output, "phenotype_table.tsv")))

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
    if os.path.exists(os.path.join(args.output, "phylogeny")) and os.path.isdir(os.path.join(args.output, "phylogeny")):
        if os.listdir(os.path.join(args.output, "phylogeny")) > 0 and not args.overwrite:
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

    print("Phylogenetic tree pipeline is finished, results can be found in the ", mash_output)

    if not args.keep_temp_files:
        print("Removing temp folder...")
        temp_folder_remover(mash_temp)

    end_time = time.time()

    print(time_function(start_time, end_time))


def fully_automated_pipeline(args):

    start_time = time.time()

    tool_list = ["PanACoTA", "snippy", "prokka", "panaroo", "pyseer"]

    for tool in tool_list:
        if not is_tool_installed(tool):
            print(f"Error: {tool} is not installed.")
            sys.exit(1)

    # Sanity checks

    if os.path.exists(args.input):
        if len(os.listdir(args.input)) == 0:
            print("Error: Input folder path is empty.")
            sys.exit(1)
    else:
        print("Error: Input folder path does not exist.")
        sys.exit(1)

    if os.path.exists(args.output) and os.path.isdir(args.output):
        if len(os.listdir(args.output)) > 0 and args.checkpoint:
            pass
        elif len(os.listdir(args.output)) > 0 and not args.overwrite:
            print("Error: Output folder is not empty.")
            print("Please provide an empty output folder or use the --overwrite option.")
            sys.exit(1)
    
    if args.ml_algorithm:
        for algorithm in args.ml_algorithm:
            if algorithm not in ["rf", "svm", "gb"]:
                print("Error: ML algorithm is not accepted.")
                sys.exit(1)

    automatix_runner(args)

    end_time = time.time()

    print(time_function(start_time, end_time))


def structman_pipeline(args):

    if args.input is None:
        print("Error: Input file or folder path is required.")
        sys.exit(1)

    if args.output is None:
        print("Error: Output folder path is required.")
        sys.exit(1)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    # Check if output folder empty
    if os.path.exists(os.path.join(args.output, "structman")) and os.path.isdir(os.path.join(args.output, "structman")):
        if not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")

    # Check if temp folder exists and create if not
    if not os.path.exists(args.temp):
        os.mkdir(args.temp)

    # Check if temp folder empty
    if os.path.exists(os.path.join(args.temp, "structman")) and os.path.isdir(os.path.join(args.temp, "structman")):
        if not args.overwrite:
            print("Error: Temp folder is not empty.")
            sys.exit(1)
    else:
        os.mkdir(os.path.join(args.temp, "structman"))

    structman_input_creator(args)

def prediction_pipeline(args):

    start_time = time.time()

    if args.input is None and args.prediction_table is None:
        print("Error: Input file or folder path is required.")
        sys.exit(1)

    if args.input and args.prediction_table:
        print("Error: Both input and prediction table cannot be used together.")
        sys.exit(1)
    
    # Check if output folder exists and is a directory
    if os.path.exists(os.path.join(args.output)) and os.path.isdir(os.path.join(args.output)):
        # Check if the directory is not empty
        if os.listdir(os.path.join(args.output)) and not args.overwrite:
            print("Error: Output folder is not empty.")
            sys.exit(1)

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    
    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")
    
    if not os.path.exists(args.temp):
        os.mkdir(args.temp)

    snippy_flag = True
    prokka_flag = True
    panaroo_flag = False
    gpa_flag = True
    prps_flag = False

    if args.no_gene_presence_absence and args.use_panaroo:
        print("Error: Gene presence absence and Panaroo cannot be used together.")
        sys.exit(1)

    if args.prps[0] is not None:
        prps_flag = True
            #args.prps[0] = os.path.join(args.output, "prps", "prps_scores.tsv")

    if args.no_gene_presence_absence:
        panaroo_flag = False
        gpa_flag = False
    if args.no_gene_annotation:
        prokka_flag = False
    if args.use_panaroo:
        panaroo_flag = True

    if args.verbosity > 3:
        print("Will run snippy: ", snippy_flag)
        print("Will run prokka: ", prokka_flag)
        print("Will run panaroo: ", panaroo_flag)
        print("Will run gene presence absence: ", gpa_flag)

    # Create the output folder
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(os.path.join(args.output, "snippy")):
        os.mkdir(os.path.join(args.output, "snippy"))
    if not os.path.exists(os.path.join(args.output, "prokka")):
        os.mkdir(os.path.join(args.output, "prokka"))
    if not os.path.exists(os.path.join(args.output, "panaroo")):
        os.mkdir(os.path.join(args.output, "panaroo"))
    if not os.path.exists(os.path.join(args.output, "cd-hit")):
        os.mkdir(os.path.join(args.output, "cd-hit"))

    snippy_output = os.path.join(args.output, "snippy")
    prokka_output = os.path.join(args.output, "prokka")
    panaroo_output = os.path.join(args.output, "panaroo")
    cd_hit_output = os.path.join(args.output, "cd-hit")

    # Create the temp folder for the panaroo input
    if not os.path.exists(os.path.join(args.temp, "panaroo")):
        os.mkdir(os.path.join(args.temp, "panaroo"))
    # Create the temp folder for the cdhit input
    if not os.path.exists(os.path.join(args.temp, "cdhit")):
        os.mkdir(os.path.join(args.temp, "cdhit"))


    if args.verbosity > 3:
        print(f"Output folder created: {args.output}")
        print(f"Snippy output folder created: {snippy_output}")
        print(f"Prokka output folder created: {prokka_output}")
        print(f"Panaroo output folder created: {panaroo_output}")

    do_not_remove_temp = False

    if args.keep_temp_files:
        do_not_remove_temp = True

    if args.custom_database and args.input:
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

    if args.input:

        if args.reference is None:
            print("Error: Reference genome is required.")
            sys.exit(1)

        input_files = process_data_for_prediction(args.input)

        with open(os.path.join(args.output, "strains.txt"), "w") as outfile:
            for strain in input_files:
                # Make sure path is same in both Windows and Linux
                strain_path = os.path.join(args.input, strain)
                strain_path = os.path.abspath(strain_path)
                strain_path = strain_path.replace("\\", "/")
                outfile.write(f"{strain_path}\n")

        input_file = os.path.join(args.output, "strains.txt")

        if args.verbosity > 3:
            print(f"Strains.txt file created: {input_file}")

        if input_file is not None:
            random_names = random_name_giver(
                input_file, os.path.join(args.output, "random_names.txt"))
            if args.verbosity > 3:
                print(f"Random names file created: {os.path.join(args.output, 'random_names.txt')}")

        with open(os.path.join(args.output, "random_names.txt")) as random_names_file:
            random_names = {}
            for random_names_file_line in random_names_file.readlines():
                random_names[random_names_file_line.split(
                    "\t")[0].strip()] = random_names_file_line.split("\t")[1].strip()

        input_file = os.path.join(args.output, "strains.txt")

        strain_list = []

        with open(input_file, "r") as infile:
            added_strains = []
            lines = infile.readlines()
            for line in lines:
                if os.path.splitext(line.split("/")[-1].strip())[0] not in added_strains:
                    added_strains.append(os.path.splitext(
                        line.split("/")[-1].strip())[0])
                    strain_list.append(line.strip())

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
            if prokka_flag and snippy_flag:
                if random_names[strain] in prokka_output_strains and random_names[strain] in snippy_output_strains:
                    strains_to_be_processed.append(random_names[strain])
                else:
                    strains_to_be_skiped.append(random_names[strain])
            elif prokka_flag and not snippy_flag:
                if random_names[strain] in prokka_output_strains:
                    strains_to_be_processed.append(random_names[strain])
                else:
                    strains_to_be_skiped.append(random_names[strain])
            elif snippy_flag and not prokka_flag:
                if random_names[strain] in snippy_output_strains:
                    strains_to_be_processed.append(random_names[strain])
                else:
                    strains_to_be_skiped.append(random_names[strain])

        print(f"Number of strains processed: {len(strains_to_be_processed)}")
        print(f"Number of strains skipped: {len(strains_to_be_skiped)}")

        snippy_processed_file_creator(snippy_output, os.path.join(
            args.output, "snippy_processed_strains.txt"))

        print("Creating binary mutation table...")
        # Create the binary table
        binary_table_creator(snippy_output, os.path.join(
            args.output, "binary_mutation_table.tsv"), args.threads, strains_to_be_processed)

        print("Creating annotation table...")

        annotation_file_from_snippy(snippy_output, args.output)

        if panaroo_flag:
            try:
                print("Creating panaroo input...")
                # Create the panaroo input
                panaroo_input_creator(os.path.join(args.output, "random_names.txt"), prokka_output, os.path.join(
                    args.temp, "panaroo"), strains_to_be_processed)

                print("Running panaroo...")
                # Run panaroo
                panaroo_runner(os.path.join(args.temp, "panaroo"), panaroo_output, os.path.join(
                    args.temp, "panaroo_log.txt"), args.threads)

                print("Adding gene presence absence information to the binary table...")
                # Add gene presence absence information to the binary table

                if not os.path.exists(os.path.join(panaroo_output, "gene_presence_absence.csv")):
                    print("Warning: Gene presence absence file does not exist.")
                    print(
                        "Gene presence absence information will not be added to the binary table.")
                    do_not_remove_temp = True

                else:
                    binary_mutation_table_gpa_information_adder_panaroo(os.path.join(args.output, "binary_mutation_table.tsv"), os.path.join(
                        panaroo_output, "gene_presence_absence.csv"), os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv"))
                    do_not_remove_temp = False

            except Exception as e:
                print("Error: Panaroo could not be run.")
                print(e)
                do_not_remove_temp = True

        if gpa_flag:

            try: 

                print("Creating gene presence absence input...")
                # Create the gene presence absence information
                print(f"CD-HIT preprecessor is running...")
                cdhit_preprocessor(os.path.join(args.output, "random_names.txt"), prokka_output, os.path.join(args.temp, "cdhit"), strains_to_be_processed)

                shutil.copy(os.path.join(args.temp, 'cdhit', 'protein_positions.csv'), os.path.join(args.output, 'cd-hit', 'protein_positions.csv'))

                print(f"CD-HIT is running...")
                cdhit_runner(os.path.join(args.temp, "cdhit", "combined_proteins.faa"), os.path.join(args.output, "cd-hit", "cdhit_output.txt"), n_cpu=args.threads)

                print(f"Gene presence-absence matrix is being created...")
                gene_presence_absence_file_creator(os.path.join(args.output, "cd-hit", "cdhit_output.txt.clstr"), strains_to_be_processed, os.path.join(args.temp, "cdhit"))

                print("Adding gene presence absence information to the binary table...")
                # Add gene presence absence information to the binary table

                if not os.path.exists(os.path.join(args.temp, "cdhit", "gene_presence_absence_matrix.csv")):
                    print("Warning: Gene presence absence file does not exist.")
                    print(
                        "Gene presence absence information will not be added to the binary table.")
                    do_not_remove_temp = True
                
                else:
                    binary_mutation_table_gpa_information_adder(os.path.join(args.output, "binary_mutation_table.tsv"), os.path.join(args.temp, "cdhit", "gene_presence_absence_matrix.csv"), os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv"))
                    do_not_remove_temp = False

            except Exception as e:
                print("Error: CD-HIT could not be run.")
                print(e)
                do_not_remove_temp = True

        if panaroo_flag or gpa_flag:
            if not os.path.exists(os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv")):
                print("Error: Gene presence absence information could not be added to the binary table.")
                print("Please check the logs for more information.")
                print("Using just mutations...")
                binary_table_to_be_used = os.path.join(args.output, "binary_mutation_table.tsv")
            
            else:
                binary_table_to_be_used = os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv")
        else:
            binary_table_to_be_used = os.path.join(args.output, "binary_mutation_table.tsv")

        print("Equalizing columns in the binary table...")

        if prps_flag:
            os.makedirs(os.path.join(args.temp, 'prps'), exist_ok=True)
            prps_ml_preprecessor(args.model_binary_table, args.prps[0], args.prps[1], os.path.join(args.temp, "prps"))
            equalize_columns(os.path.join(args.temp, 'prps', 'prps_filtered_table.tsv'), binary_table_to_be_used, os.path.join(
                args.output, "equalized_binary_mutation_table.tsv"))
        else:
            equalize_columns(args.model_binary_table, binary_table_to_be_used, os.path.join(
                args.output, "equalized_binary_mutation_table.tsv"))
        
        print("Running prediction...")
        
        predict(args.model, os.path.join(args.output, "equalized_binary_mutation_table.tsv"), args.output)

        print(f"Prediction is finished. Results can be found in the {os.path.join(args.output, 'predictions.csv')}")

    if args.prediction_table:
        prediction_table_to_be_used = args.prediction_table

        print("Equalizing columns in the binary table...")

        if prps_flag:
            os.makedirs(os.path.join(args.temp, 'prps'), exist_ok=True)
            prps_ml_preprecessor(prediction_table_to_be_used, args.prps[0], args.prps[1], os.path.join(args.temp, 'prps'))
            equalize_columns(os.path.join(args.temp, 'prps', 'prps_filtered_table.tsv'), prediction_table_to_be_used, os.path.join(
                args.output, "equalized_binary_mutation_table.tsv"))

        else:
            equalize_columns(args.model_binary_table, prediction_table_to_be_used, os.path.join(
                args.output, "equalized_binary_mutation_table.tsv"))
        
        print("Running prediction...")

        predict(args.model, os.path.join(args.output, "equalized_binary_mutation_table.tsv"), args.output)

        print(f"Prediction is finished. Results can be found in the {os.path.join(args.output, 'predictions.csv')}")

    if args.keep_temp_files:
        print("Warning, temp files will be kept this might take up space.")

    if not args.keep_temp_files:
        if not do_not_remove_temp:
            print("Removing temp folder...")
            temp_folder_remover(os.path.join(args.temp))
            temp_folder_remover(os.path.join(args.output, "snippy"))
            temp_folder_remover(os.path.join(args.output, "prokka"))
            temp_folder_remover(os.path.join(args.output, "panaroo"))

        else:
            print(
                "Warning: Temp folder will not be removed because gene presence absence file does not exist.")
            print(
                "Temp folder can be used for re-run panaroo again for gene presence absence information.")
            print("You can remove the temp folder manually. Temp folder path: ",
                  os.path.join(args.temp))
             
    end_time = time.time()
    print(time_function(start_time, end_time))

if __name__ == "__main__":
    main()
