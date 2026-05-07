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
import psutil

from sr_amr.utils import is_tool_installed, temp_folder_remover, time_function, copy_and_zip_file, ensure_conda_env, conda_env_wrapper
from sr_amr.version import __version__

import subprocess
import json

import warnings

warnings.filterwarnings("ignore")

accepted_ml_algorithms = ["rf", "svm", "gb", "histgb", "xgb", "lr"]

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

    parser_automatix.add_argument('--variant_calling_tool', type=str, help='variant calling tool to use, available selections: [snippy], default=snippy', default="snippy")
    parser_automatix.add_argument('--annotation_tool', type=str, help='annotation tool to use, available selections: [prokka, bakta], default=prokka', default="prokka")
    parser_automatix.add_argument('--gene_presence_absence_analysis_tool', type=str, help='gene presence absence analysis tool to use, available selections: [cd-hit, panaroo], default=cd-hit', default="cd-hit")

    parser_automatix.add_argument('--prokka_custom_database', type=str, nargs=2,
                                      help='creates and uses custom database for prokka, require path of the fasta file and genus name, default=None')
    
    parser_automatix.add_argument('--bakta_db', type=str, help='path to bakta database, if none given and bakta as annotation tool selected, will automatically download, default=None')

    parser_automatix.add_argument('--only_variants', action='store_true', help='only creates binary mutation table with mutations, without gene presence absence information, default=False')
    
    parser_automatix.add_argument('--no_feature_importance_analysis', action='store_true',help='do not run feature importance analysis on ML step, default=False')
    parser_automatix.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_automatix.add_argument(
        '--threads', type=int, help='number of threads to use, -1 for all available threads, default=1', default=1)
    parser_automatix.add_argument(
        '--ram', type=int, help='amount of ram to use in GB, default=8', default=8)
    parser_automatix.add_argument(
        '--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_automatix.add_argument('--overwrite', action='store_true',
                                  help='overwrite the output and temp folder if exists, default=False')
    parser_automatix.add_argument('--ml_algorithm', nargs='+',
                              help='classification algorithm to be used, available selections: [rf, svm, gb, histgb, lr, xgb], default=[rf, svm, lr, xgb]', default=["rf", "svm", "lr", "xgb"])
    parser_automatix.add_argument('--no_ml', action='store_true', help='do not run machine learning analysis, default=False')
    parser_automatix.add_argument('--fast', action='store_true', help='fast mode, does not run PanACoTA pipeline for phylogenetic tree analysis, default=False')
    parser_automatix.add_argument('--checkpoint', action='store_true',
                                  help='continues run from the checkpoint, default=False')

    parser_automatix.add_argument('--run_qc', action='store_true', help='run automated QC on input genomes, default=False')
    parser_automatix.add_argument('--qc_length_threshold', type=float, help='fraction of median length allowed, default=0.1', default=0.1)
    parser_automatix.add_argument('--qc_max_contigs', type=int, help='maximum allowed number of contigs, default=500', default=500)

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
    
    parser_main_pipeline.add_argument('--create_phenotype_from_folder', action='store_true',
                                      help='create phenotype file from the input folder that contains genomic files, default=True')
    
    parser_main_pipeline.add_argument('--variant_calling_tool', type=str, help='variant calling tool to use, available selections: [snippy], default=snippy', default="snippy")
    parser_main_pipeline.add_argument('--annotation_tool', type=str, help='annotation tool to use, available selections: [prokka, bakta], default=prokka', default="prokka")
    parser_main_pipeline.add_argument('--gene_presence_absence_analysis_tool', type=str, help='gene presence absence analysis tool to use, available selections: [cd-hit, panaroo], default=cd-hit', default="cd-hit")

    parser_main_pipeline.add_argument('--prokka_custom_database', type=str, nargs=2,
                                      help='creates and uses custom database for prokka, require path of the fasta file and genus name, default=None')
    
    parser_main_pipeline.add_argument('--bakta_db', type=str, help='path to bakta database, if none given and bakta as annotation tool selected, will automatically download, default=None')

    parser_main_pipeline.add_argument('--only_variants', action='store_true', help='only creates binary mutation table with mutations, without gene presence absence information, default=False')
    
    parser_main_pipeline.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_main_pipeline.add_argument('--overwrite', action='store_true',
                                      help='overwrite the output and temp folder if exists, default=False')
    parser_main_pipeline.add_argument(
        '--keep_temp_files', action='store_true', help='keep the temporary files, default=False')
    parser_main_pipeline.add_argument(
        '--threads', type=int, help='number of threads to use, -1 for all available threads, default=1', default=1)
    parser_main_pipeline.add_argument(
        '--ram', type=int, help='amount of ram to use in GB, -1 for all available ram, default=4', default=4)
    
    parser_main_pipeline.add_argument('--run_qc', action='store_true', help='run automated QC on input genomes, default=False')
    parser_main_pipeline.add_argument('--qc_length_threshold', type=float, help='fraction of median length allowed, default=0.1', default=0.1)
    parser_main_pipeline.add_argument('--qc_max_contigs', type=int, help='maximum allowed number of contigs, default=500', default=500)
    
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
        '--threads', type=int, help='number of threads to use, -1 for all available threads, default=1', default=1)
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
                             help='number of threads to use, -1 for all available threads, default=1', default=1)
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
        '--threads', type=int, help='number of threads to use, -1 for all available threads, default=1', default=1)
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
    parser_ml.add_argument('--sail_max_time', type=int,
                           help='maximum time in seconds for datasail to run, default=600', default=600)
    parser_ml.add_argument('--sail_stratify', action='store_true', help='whether to stratify the datasail split, default=False')
    parser_ml.add_argument('--sail_distance_matrix', type=str, help='distance matrix file path for datasail if it has been generated previously, default=None')
    parser_ml.add_argument('--train_strains_file', type=str,
                           help='train strains file path', default=None)
    parser_ml.add_argument('--test_strains_file', type=str,
                           help='test strains file path', default=None)
    parser_ml.add_argument('--validation_strains_file', type=str,
                            help='validation strains file path', default=None)
    parser_ml.add_argument('--overwrite', action='store_true',
                           help='overwrite the output folder if exists, default=False')
    parser_ml.add_argument(
        '--threads', type=int, help='number of threads to use, -1 for all available threads, default=1', default=1)
    parser_ml.add_argument(
        '--ram', type=int, help='amount of ram to use in GB, default=4', default=4)
    parser_ml.add_argument(
        '--temp', type=str, help='path of the temporary directory, default=output_folder/temp')
    parser_ml.add_argument('--keep_temp_files', action='store_true',
                           help='keep the temporary files, default=False')
    parser_ml.add_argument('--ml_algorithm', type=str,
                           help='classification algorithm to be used, available selections: [rf, svm, gb, histgb, xgb], default=rf', default="rf")
    parser_ml.add_argument('--test_train_split', type=float,
                           help='test train split ratio, default=0.20', default=0.20)
    parser_ml.add_argument('--random_state', type=int,
                           help='random state, default=42', default=42)
    parser_ml.add_argument('--n_estimators', type=int,
                           help='number of estimators for random forest, default=100', default=100)
    parser_ml.add_argument('--max_depth', type=int,
                           help='max depth for random forest & gradient boosting, default=10', default=10)
    parser_ml.add_argument('--min_samples_split', type=int,
                           help='min samples split for random forest, default=2', default=2)
    parser_ml.add_argument('--min_samples_leaf', type=int,
                           help='min samples leaf for random forest & gradient boosting, default=1', default=1)
    parser_ml.add_argument('--max_features', type=str,
                           help='max features for random forest, default=auto', default="auto")
    parser_ml.add_argument('--resampling_strategy', type=str,
                           help='resampling strategy for ml, available selections: [holdout, cv], default=holdout', default="holdout")
    parser_ml.add_argument(
        '--cv', type=int, help='applies Cross-Validation with given number of splits, default=4', default=4)
    parser_ml.add_argument(
        '--scoring', type=str, help='scoring method for cross-validation, available selections: [MCC,accuracy,f1,roc_auc], default=MCC', default="MCC")
    parser_ml.add_argument('--save_model', action='store_true',
                           help='save the ml model, default=False')
    parser_ml.add_argument('--feature_importance_analysis', action='store_true',
                           help='analyze feature importance, default=False')
    parser_ml.add_argument('--important_feature_limit', type=int,
                           help='number of reported maximum number of features in FIA file, if want all > 0.0, use -1, default=25', default=25)
    parser_ml.add_argument('--feature_importance_analysis_number_of_repeats', type=int,
                           help='number of repeats for feature importance analysis should be given with --feature_importance_analysis option, default=5', default=5)
    parser_ml.add_argument('--feature_importance_analysis_strategy', type=str,
                           help='strategy for feature importance analysis, available selections: [gini, permutation_importance], default=gini', default="gini")
    parser_ml.add_argument('--svm_kernel', type=str,
                           help='kernel for svm, available selections: [linear, poly, rbf, sigmoid], default=linear', default="linear")
    parser_ml.add_argument('--no_stratify_split', type=str,
                           help='if given, does not uses stratify in random split', default="False")
    parser_ml.add_argument('--param_grid_size', type=str,
                           help='size of the parameter grid for best parameter search,available selections: [small, medium, large], default=small', default="small")
    parser_ml.add_argument('--param_grid_low_memory_mode', action='store_true',
                           help='Only for XGB, if given, uses low memory mode for parameter grid search, if given memory is not 100 times more than datasize, will automatically activate, default=False')
    parser_ml.add_argument('--device', type=str,
                           help='Only for XGB, device to use for training, available selections: [cpu, cuda], default=cpu', default="cpu")
    parser_ml.add_argument('--parameter_search_strategy', type=str,
                           help='parameter search strategy, available selections: [grid_search, random_search], default=grid_search', default="grid_search")
    parser_ml.add_argument('--parameter_search_n_iter', type=int,
                           help='number of iterations for random parameter search, should be used with --parameter_search_strategy random_search option, default=20', default=20)
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
                                   help='number of threads to use, -1 for all available threads, default=1', default=1)
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

def run_variant_calling_and_annotation(strain, random_names, args):
    # if args.ram / args.threads < 8:
    #     print("Warning: Not enough ram for the processes. Minimum 8 GB of ram per thread is recommended.")
    from sr_amr.binary_tables import snippy_runner, prokka_runner, bakta_runner
    
    # Snippy creates issue with high memory usage, so it is limited to 100 GB
    snippy_ram = args.ram
    if args.ram > 100:
        #print(f"Warning: Due to restrictions of snippy, ram is limited to 100 GB for snippy run only.")
        snippy_ram = 100

    snippy_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[
                    0]], os.path.join(args.output, args.variant_calling_tool), args.reference, f"{args.temp}/snippy_log.txt", 1, snippy_ram, env_name=f"alpar-{args.variant_calling_tool}")
    
    if not args.only_variants:

        if args.annotation_tool == "bakta":
            bakta_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[
                        0]], os.path.join(args.output, args.annotation_tool), f"{args.temp}/bakta_log.txt", f"{args.temp}/bakta/", 1, args.bakta_db, env_name=f"alpar-{args.annotation_tool}")
            
        elif args.annotation_tool == "prokka":
            prokka_runner(strain, random_names[os.path.splitext(strain.split("/")[-1].strip())[
                        0]], os.path.join(args.output, args.annotation_tool), args.reference, f"{args.temp}/prokka_log.txt", 1, args.prokka_custom_database[1] if args.prokka_custom_database else None, env_name=f"alpar-{args.annotation_tool}")


def binary_table_pipeline(args):

    start_time = time.time()

    from sr_amr.binary_tables import check_and_download_bakta_db, random_name_giver, prokka_create_database, snippy_processed_file_creator, binary_table_creator, annotation_file_from_snippy, panaroo_input_creator, panaroo_runner, binary_mutation_table_gpa_information_adder_panaroo, cdhit_preprocessor, cdhit_runner, gene_presence_absence_file_creator, binary_mutation_table_gpa_information_adder
    from sr_amr.qc import run_qc_pipeline

    # Sanity checks

    status = -1

    # Annotation tools: Prokka, Bakta
    # Gene presence absence tools: Panaroo, CD-HIT
    # Variant calling tool: Snippy

    print("Starting the binary table creation pipeline...")
    print("Running sanity checks...")

    ensure_conda_env(f"alpar-{args.variant_calling_tool}")
    if not args.only_variants:
        ensure_conda_env(f"alpar-{args.annotation_tool}")
        ensure_conda_env(f"alpar-{args.gene_presence_absence_analysis_tool}")

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

    if args.prokka_custom_database:
        if len(args.prokka_custom_database) != 2:
            print("Error: Prokka custom database option should have two arguments. Required arguemnts: path of the fasta file and genus name")
            sys.exit(1)

        if not os.path.exists(args.prokka_custom_database[0]):
            print("Error: Prokka custom database fasta file does not exist.")
            sys.exit(1)

        if pathlib.Path(args.prokka_custom_database[0]).suffix != ".fasta":
            print("Error: Prokka custom database file extension is not accepted.")
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

    if args.temp is None:
        args.temp = os.path.join(args.output, "temp")

    if not os.path.exists(args.temp):
        os.makedirs(args.temp, exist_ok=True)
    else:
        print("Warning: Temp folder already exists. Will be overwritten by the run.")

    # Check if threads is positive
    if args.threads is not None:
        if args.threads < 0:
            # Use all available threads
            print(f"Using all available threads: {multiprocessing.cpu_count()}")
            args.threads = multiprocessing.cpu_count()
        elif args.threads == 0:
            print("Error: Number of threads should be positive.")
            sys.exit(1)
    else:
        args.threads = 1

    # Check if ram is positive
    if args.ram is not None:
        if args.ram <= 0:
            print("Using all available ram.")
            args.ram = int(psutil.virtual_memory().total / (1024 ** 3))
    else:
        args.ram = 4

    snippy_flag = True

    if args.only_variants:
        gene_annotation_flag = False
    else:
        gene_annotation_flag = True

    if args.verbosity > 3:
        print("Will run snippy: ", snippy_flag)
        print("Will run gene annotation: ", gene_annotation_flag)

    # Create the output folder
    os.makedirs(args.output, exist_ok=True)

    os.makedirs(os.path.join(args.output, "snippy"), exist_ok=True)
    
    if not args.only_variants:
        os.makedirs(os.path.join(args.output, args.annotation_tool), exist_ok=True)
        os.makedirs(os.path.join(args.output, args.gene_presence_absence_analysis_tool), exist_ok=True)
        if args.verbosity > 3:
            print(f"Output folder created: {args.output}")
        
        os.makedirs(os.path.join(args.temp, args.gene_presence_absence_analysis_tool), exist_ok=True)
        os.makedirs(os.path.join(args.temp, args.annotation_tool), exist_ok=True)

        if args.annotation_tool == "bakta":
            print("Checking Bakta database...")
            args.bakta_db = check_and_download_bakta_db(args.bakta_db, f"{args.temp}/bakta_db_init.log")
            if args.bakta_db is None:
                print("Error: Could not find or download Bakta database. Exiting.")
                sys.exit(1)

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

                STRAIN_FILE_MINIMUM_SIZE_LIMIT = 1000

                with open(os.path.join(args.output, "strains.txt"), "a") as outfile:
                    for strain in resistant_strains:
                        # Make sure path is same in both Windows and Linux
                        strain_path = os.path.join(resistant_path, strain)
                        strain_path = os.path.abspath(strain_path)
                        strain_path = strain_path.replace("\\", "/")
                        # Check if file size is larger than minimum size limit, teoratically bacterial wgs fasta files should be larger than 1 KB
                        if os.path.getsize(strain_path) < STRAIN_FILE_MINIMUM_SIZE_LIMIT:
                            print(f"Warning: {strain_path} file size is smaller than {STRAIN_FILE_MINIMUM_SIZE_LIMIT} bytes, skipping this file.")
                            continue
                        outfile.write(f"{strain_path}\n")
                    for strain in susceptible_strains:
                        strain_path = os.path.join(susceptible_path, strain)
                        strain_path = os.path.abspath(strain_path)
                        strain_path = strain_path.replace("\\", "/")
                        if os.path.getsize(strain_path) < STRAIN_FILE_MINIMUM_SIZE_LIMIT:
                            print(f"Warning: {strain_path} file size is smaller than {STRAIN_FILE_MINIMUM_SIZE_LIMIT} bytes, skipping this file.")
                            continue
                        outfile.write(f"{strain_path}\n")

            input_file = os.path.join(args.output, "strains.txt")

            if args.verbosity > 3:
                print(f"Strains.txt file created: {input_file}")

        if args.run_qc:
            strains_to_qc = []
            with open(os.path.join(args.output, "strains.txt"), "r") as infile:
                for line in infile:
                    strains_to_qc.append(line.strip())
            
            passed_strains = run_qc_pipeline(strains_to_qc, args.output, args.qc_length_threshold, args.qc_max_contigs)
            
            with open(os.path.join(args.output, "strains.txt"), "w") as outfile:
                for strain in passed_strains:
                    outfile.write(f"{strain}\n")

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

        if args.annotation_tool == "prokka" and args.prokka_custom_database:
            if len(args.prokka_custom_database) != 2:
                print("Error: Custom database option should have two arguments.")
                sys.exit(1)

            if not os.path.exists(args.prokka_custom_database[0]):
                print("Error: Custom database fasta file does not exist.")
                sys.exit(1)

            print("Creating custom database...")
            prokka_create_database(
                args.prokka_custom_database[0], args.prokka_custom_database[1], args.temp, args.threads, args.ram, env_name="alpar-cdhit")
            print("Custom database created.")

            with open(os.path.join(args.temp, "status.txt"), "w") as outfile:
                outfile.write(f"1")
                status = 1

    if status < 2:
        # Run snippy and annotation tool in parallel for each strain

        print(f"Number of strains to be processed: {len(strain_list)}")
        print("Running variant calling and annotation...")

        num_parallel_tasks = args.threads

        if args.ram / args.threads < 8 and args.annotation_tool == "bakta":
            print("Warning: Not enough ram for the processes. Minimum 4 GB of ram per thread is recommended.")
            print(f"Current ram per thread: {int(args.ram / args.threads)} GB")
            print(f"Setting number of parallel tasks to {int(args.ram / 8)} to avoid memory issues.")
            num_parallel_tasks = int(args.ram / 8)

        params = [(strain, random_names, args) for strain in strain_list]

        with multiprocessing.Pool(num_parallel_tasks) as pool:
            pool.starmap(run_variant_calling_and_annotation, params)

        # We will use a status file to indicate checkpoints

        with open(os.path.join(args.temp, "status.txt"), "w") as outfile:
            outfile.write(f"2")
            status = 2

    if status < 3:

        strains_to_be_processed = []

        variant_calling_output_strains = os.listdir(os.path.join(args.output, args.variant_calling_tool))
        if not args.only_variants:
            annotation_tool_output_strains = os.listdir(os.path.join(args.output, args.annotation_tool))

        strains_to_be_skiped = []

        for strain in random_names.keys():
            if random_names[strain] in variant_calling_output_strains:
                if not args.only_variants:
                    if random_names[strain] in annotation_tool_output_strains:
                        strains_to_be_processed.append(random_names[strain])
                    else:
                        print(f"Warning: {strain} is missing from annotation tool output, skipping this strain.")
                        strains_to_be_skiped.append(random_names[strain])
                else:
                    strains_to_be_processed.append(random_names[strain])
            else:
                print(f"Warning: {strain} is missing from variant calling output, skipping this strain.")
                strains_to_be_skiped.append(random_names[strain])

        print(f"Number of strains processed: {len(strains_to_be_processed)}")
        print(f"Number of strains skipped: {len(strains_to_be_skiped)}")

        snippy_processed_file_creator(os.path.join(args.output, args.variant_calling_tool), os.path.join(
            args.output, "snippy_processed_strains.txt"))

        print("Creating binary mutation table...")
        # Create the binary table
        binary_table_creator(os.path.join(args.output, args.variant_calling_tool), os.path.join(
            args.output, "binary_mutation_table.tsv"), args.threads, strains_to_be_processed, args.temp)

        print("Creating annotation table...")

        annotation_file_from_snippy(os.path.join(args.output, args.variant_calling_tool), args.output)

        with open(os.path.join(args.temp, "status.txt"), "w") as outfile:
            outfile.write(f"3")
            status = 3

    if status < 4:

        if not args.only_variants:
            
            try:
                print("Strains to be processed list is empty, checking file")
                if os.path.exists(os.path.join(args.output, "snippy_processed_strains.txt")):
                    with open(os.path.join(args.output, "snippy_processed_strains.txt"), "r") as infile:
                        strains_to_be_processed = [line.strip() for line in infile.readlines()]
                        print(f"Strains to be processed list is created from file, number of strains: {len(strains_to_be_processed)}")
                else:
                    print("Error: Strains to be processed list is empty and snippy_processed_strains.txt file does not exist.")
                    print("Gene presence absence analysis will be skipped.")
                    do_not_remove_temp = True
                    strains_to_be_processed = []
            except Exception as e:
                print("Error: Could not create strains to be processed list.")
                print(e)
                print("Gene presence absence analysis will be skipped.")
                do_not_remove_temp = True
                strains_to_be_processed = []

            os.makedirs(os.path.join(args.output, args.annotation_tool), exist_ok=True)

            if args.annotation_tool == "bakta" and args.gene_presence_absence_analysis_tool == "panaroo":
                    print("Warning: Panaroo is not fully compatible with Bakta annotations, automatically falling back to cd-hit for gene presence absence analysis instead. If you want to use panaroo, please use prokka for annotation.")
                    print("For more information: https://github.com/gtonkinhill/panaroo/issues/373")
                    args.gene_presence_absence_analysis_tool = "cd-hit"
            
            if args.gene_presence_absence_analysis_tool == "panaroo":

                try:
                    print("Creating panaroo input...")
                    # Create the panaroo input
                    panaroo_input_creator(os.path.join(args.output, "random_names.txt"), os.path.join(args.output, args.annotation_tool), os.path.join(
                        args.temp, args.gene_presence_absence_analysis_tool), strains_to_be_processed)

                    print("Running panaroo...")
                    # Run panaroo
                    panaroo_runner(os.path.join(args.temp, args.gene_presence_absence_analysis_tool), os.path.join(args.output, args.gene_presence_absence_analysis_tool), os.path.join(
                        args.temp, "panaroo_log.txt"), args.threads, env_name=f"alpar-{args.gene_presence_absence_analysis_tool}")

                    print("Adding gene presence absence information to the binary table...")
                    # Add gene presence absence information to the binary table

                    if not os.path.exists(os.path.join(os.path.join(args.output, args.gene_presence_absence_analysis_tool), "gene_presence_absence.csv")):
                        print("Warning: Gene presence absence file does not exist.")
                        print(
                            "Gene presence absence information will not be added to the binary table.")
                        do_not_remove_temp = True

                    else:
                        binary_mutation_table_gpa_information_adder_panaroo(os.path.join(args.output, "binary_mutation_table.tsv"), os.path.join(
                            os.path.join(args.output, args.gene_presence_absence_analysis_tool), "gene_presence_absence.csv"), os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv"))
                        do_not_remove_temp = False

                except Exception as e:
                    print("Error: Panaroo could not be run.")
                    print(e)
                    do_not_remove_temp = True

            if args.gene_presence_absence_analysis_tool == "cd-hit":

                try: 
                    print("Creating gene presence absence input...")
                    # Create the gene presence absence information
                    print(f"CD-HIT preprecessor is running...")
                    cdhit_preprocessor(os.path.join(args.output, "random_names.txt"), os.path.join(args.output, args.annotation_tool), os.path.join(args.temp, args.gene_presence_absence_analysis_tool), strains_to_be_processed)

                    shutil.copy(os.path.join(args.temp, args.gene_presence_absence_analysis_tool, 'protein_positions.csv'), os.path.join(args.output, args.gene_presence_absence_analysis_tool, 'protein_positions.csv'))

                    print(f"CD-HIT is running...")
                    cdhit_runner(os.path.join(args.temp, args.gene_presence_absence_analysis_tool, "combined_proteins.faa"), os.path.join(args.output, args.gene_presence_absence_analysis_tool, "cdhit_output.txt"), n_cpu=args.threads, env_name=f"alpar-{args.gene_presence_absence_analysis_tool}")

                    print(f"Gene presence-absence matrix is being created...")
                    gene_presence_absence_file_creator(os.path.join(args.output, args.gene_presence_absence_analysis_tool, "cdhit_output.txt.clstr"), strains_to_be_processed, os.path.join(args.temp, args.gene_presence_absence_analysis_tool))

                    print("Adding gene presence absence information to the binary table...")
                    # Add gene presence absence information to the binary table

                    if not os.path.exists(os.path.join(args.temp, args.gene_presence_absence_analysis_tool, "gene_presence_absence_matrix.csv")):
                        print("Warning: Gene presence absence file does not exist.")
                        print(
                            "Gene presence absence information will not be added to the binary table.")
                        do_not_remove_temp = True
                    
                    else:
                        binary_mutation_table_gpa_information_adder(os.path.join(args.output, "binary_mutation_table.tsv"), os.path.join(args.temp, args.gene_presence_absence_analysis_tool, "gene_presence_absence_matrix.csv"), os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv"))
                        do_not_remove_temp = False

                except Exception as e:
                    print("Error: CD-HIT could not be run.")
                    print(e)
                    do_not_remove_temp = True

        if args.create_phenotype_from_folder:
            print("Creating phenotype dataframe...")
            from sr_amr.binary_tables import phenotype_dataframe_creator, phenotype_dataframe_creator_post_processor
            # Create the phenotype dataframe
            phenotype_dataframe_creator(args.input, os.path.join(
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


@conda_env_wrapper("alpar-panacota")
def panacota_pipeline(args):

    start_time = time.time()

    from sr_amr.panacota import panacota_pre_processor, panacota_pipeline_runner, panacota_post_processor

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
                                type=args.data_type, min_seq_id=args.min_seq_id, mode=args.clustering_mode, core_genome_percentage=args.core_genome_percentage, env_name="alpar-panacota")

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


@conda_env_wrapper("alpar-pyseer")
def gwas_pipeline(args):

    start_time = time.time()

    from sr_amr.gwas import pyseer_genotype_matrix_creator, pyseer_phenotype_file_creator, pyseer_similarity_matrix_creator, pyseer_runner, pyseer_post_processor, pyseer_gwas_graph_creator, decision_tree_input_creator

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
            args.tree, os.path.join(gwas_output, "similarity_matrix.tsv"), env_name="alpar-pyseer")
        # pyseer_runner(genotype_file_path, phenotype_file_path, similarity_matrix, output_file_directory, threads):
        pyseer_runner(os.path.join(gwas_output, "genotype_matrix.tsv"), os.path.join(gwas_output, "pyseer_phenotypes"),
                    os.path.join(gwas_output, "similarity_matrix.tsv"), os.path.join(gwas_output, "gwas_results"), args.threads, env_name="alpar-pyseer")

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


@conda_env_wrapper("alpar-prps")
def prps_pipeline(args):

    start_time = time.time()

    from sr_amr.prps import PRPS_runner, PRPS_runner_continuous, PRPS_binary_check

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

        print("Checking if input file is binary or continuous...")
        if PRPS_binary_check(args.input):
            print("Input file is binary. Running PRPS for binary data...")
            PRPS_runner(args.tree, args.input, prps_output, prps_temp, env_name="alpar-prps")
        else:
            print("Input file is continuous. Running PRPS for continuous data...")
            PRPS_runner_continuous(args.tree, args.input, prps_output, prps_temp, env_name="alpar-prps")

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


def _run_datasail_pipeline(args):
    """
    Run DataSAIL pipeline for data splitting against information leakage.
    Returns: tuple of (train_strains, test_strains, datasail_output)
    Runs in the alpar-datasail conda environment via subprocess.
    """
    import json
    
    # First ensure the environment exists
    from sr_amr.utils import ensure_conda_env
    ensure_conda_env("alpar-datasail")
    
    train_strains = []
    test_strains = []
    
    datasail_temp = os.path.join(args.temp, "datasail")
    datasail_output = os.path.join(args.output, "datasail")

    if os.path.exists(os.path.join(datasail_output, "splits.tsv")):
        print("Warning: Split file already exists, it will be used for calculations. If you want to re-run the datasail, please remove the splits.tsv file from the output folder.")

    if os.path.exists(os.path.join(datasail_output, args.antibiotic, "splits.tsv")):
        print(f"Warning: Split file already exists at {os.path.join(datasail_output, args.antibiotic, 'splits.tsv')}, it will be used for calculations. If you want to re-run the datasail, please remove the splits.tsv file from the output folder.")
        datasail_output = os.path.join(datasail_output, args.antibiotic)
        
    else:
        # Check if output folder empty
        if os.path.exists(datasail_output) and os.path.isdir(datasail_output):
            if not args.overwrite:
                print("Error: Datasail output folder is not empty. If you want to overwrite, please use the --overwrite flag.")
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
        
        if args.sail_distance_matrix:
            print("Using provided distance matrix...")
            distance_matrix = args.sail_distance_matrix
        else:
            print("Creating distance matrix...")

            # Run datasail_pre_precessor in conda environment
            cmd = [
                "conda", "run", "-n", "alpar-datasail", "--no-capture-output",
                "python", "-c",
                f"""
import sys
sys.path.insert(0, '{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}')
from sr_amr.ds import datasail_pre_precessor
datasail_pre_precessor(
    '{args.sail}', 
    '{datasail_temp}', 
    {repr(random_names_dict)}, 
    '{datasail_output}', 
    {args.threads}, 
    env_name='alpar-datasail'
)
"""
            ]
            
            result = subprocess.run(cmd, capture_output=False)
            if result.returncode != 0:
                print("Error: Failed to run datasail_pre_precessor")
                sys.exit(1)
            
            distance_matrix = os.path.join(datasail_output, "distance_matrix.tsv")

            print(f"Distance matrix created: {distance_matrix}")
        
        print("Running datasail...")

        if not args.sail_epsilon:
            args.sail_epsilon = 0.1
        if not args.sail_delta:
            args.sail_delta = 0.1
        if not args.sail_solver:
            args.sail_solver = "SCIP"
        if args.sail_max_time:
            sail_max_time = args.sail_max_time
        else:
            sail_max_time = 600
        if args.sail_stratify:
            import pandas as pd
            phenotype_df = pd.read_csv(f'{args.phenotype}', sep='\t', index_col=0)
            phenotype_df_dict = phenotype_df.T.to_dict(orient='index')
        else:
            phenotype_df_dict = None
        
        # Run datasail_runner in conda environment
        cmd = [
            "conda", "run", "-n", "alpar-datasail", "--no-capture-output",
            "python", "-c",
            f"""
import sys
sys.path.insert(0, '{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}')
from sr_amr.ds import datasail_runner
import json

datasail_output_result = datasail_runner(
    '{distance_matrix}', 
    '{datasail_output}',
    splits={repr([float(1-args.test_train_split), float(args.test_train_split)])}, 
    cpus={args.threads}, 
    epsilon={args.sail_epsilon}, 
    delta={args.sail_delta}, 
    solver='{args.sail_solver}', 
    sail_max_time={sail_max_time}, 
    df_dict={repr(phenotype_df_dict)}, 
    antibiotic='{args.antibiotic}'
)
print(json.dumps({{'datasail_output': datasail_output_result}}))
"""
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("Error: Failed to run datasail_runner")
            print(result.stderr)
            sys.exit(1)
        
        # Parse output to get datasail_output path
        try:
            output_data = json.loads(result.stdout.split('\n')[-2])
            datasail_output = output_data['datasail_output']
        except:
            # If parsing fails, assume it's still the original path
            pass

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
    
    return train_strains, test_strains, datasail_output



@conda_env_wrapper("alpar-ml")
def ml_pipeline(args):

    start_time = time.time()

    if args.sail:
        ensure_conda_env("alpar-datasail")

    from sr_amr.ml import prps_ml_preprecessor, combined_ml
    from sr_amr.ml_common_files import fia_file_annotation

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
    algorithms_output = os.path.join(ml_output, args.ml_algorithm)

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
    validation_strains = []

    if args.sail:
        if args.train_strains_file or args.test_strains_file:
            print("Error: If you want to use DataSAIL as data split, train and test strains files should not be provided.")
            sys.exit(1)

        train_strains, test_strains, _ = _run_datasail_pipeline(args)

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

        if args.validation_strains_file:
            if os.path.exists(args.validation_strains_file):
                with open(args.validation_strains_file) as validation_file:
                    validation_strains_lines = validation_file.readlines()
                    for validation_strain_line in validation_strains_lines:
                        validation_strains.append(validation_strain_line.strip())
            else:
                print("Error: Validation strains file does not exist.")
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
    
    if args.validation_strains_file and not args.train_strains_file:
        print("Error: Validation strains file is provided but train strains file is missing.")
        sys.exit(1)
    
    if args.validation_strains_file and not args.test_strains_file:
        print("Error: Validation strains file is provided but test strains file is missing.")
        sys.exit(1)

    if args.no_stratify_split:
        stratiy_random_split = False
    else:
        stratiy_random_split = True

    print("Running ML pipeline...")

    with open(os.path.join(ml_output, "log_file.txt"), "w") as log_file:
        with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
            with open(os.path.join(ml_output, "parameters.txt"), "w") as parameters_file:
                parameters_file.write(f"ML algorithm: {args.ml_algorithm}\n")
                parameters_file.write(f"Scoring method: {args.scoring}\n")
                parameters_file.write(f"Test train split: {args.test_train_split}\n")
                parameters_file.write(f"Random state: {args.random_state}\n")
                parameters_file.write(f"Number of threads: {args.threads}\n")
                parameters_file.write(f"Amount of ram: {args.ram}\n")
                parameters_file.write(f"Feature importance analysis: {args.feature_importance_analysis}\n")
                if args.feature_importance_analysis:
                    parameters_file.write(f"Feature importance analysis strategy: {args.feature_importance_analysis_strategy}\n")
                    parameters_file.write(f"Number of repeats for feature importance analysis: {args.feature_importance_analysis_number_of_repeats}\n")
                if args.prps:
                    parameters_file.write(f"PRPS output used for ML pre-precessor: {args.prps}\n")
                    parameters_file.write(f"PRPS percentage used for ML pre-precessor: {PRPS_percentage}\n")
                if args.sail:
                    parameters_file.write(f"DataSAIL distance matrix used for data split: {args.sail_distance_matrix}\n")
                    parameters_file.write(f"DataSAIL epsilon used for data split: {args.sail_epsilon}\n")
                    parameters_file.write(f"DataSAIL delta used for data split: {args.sail_delta}\n")
                    parameters_file.write(f"DataSAIL solver used for data split: {args.sail_solver}\n")
                    parameters_file.write(f"DataSAIL max time used for data split: {args.sail_max_time}\n")
                    parameters_file.write(f"DataSAIL stratify used for data split: {args.sail_stratify}\n")
            fia_file = combined_ml(binary_mutation_table_path, args.phenotype, args.antibiotic, args.random_state, args.cv, args.test_train_split, ml_output, args.threads, ml_temp, args.ram, args.ml_algorithm, args.feature_importance_analysis, args.save_model, resampling_strategy=args.resampling_strategy, custom_scorer="MCC", fia_repeats=5, n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf, min_samples_split=args.min_samples_split, train=train_strains, test=test_strains, validation=validation_strains, stratify=stratiy_random_split, feature_importance_analysis_strategy=args.feature_importance_analysis_strategy, important_feature_limit=args.important_feature_limit, param_grid_size=args.param_grid_size, param_grid_low_memory_mode= args.param_grid_low_memory_mode, parameter_search_strategy=args.parameter_search_strategy, parameter_search_n_iter=args.parameter_search_n_iter, device=args.device) 

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

    from sr_amr.binary_table_threshold import binary_table_threshold_with_percentage

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

    from sr_amr.binary_tables import random_name_giver, phenotype_dataframe_creator

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


@conda_env_wrapper("alpar-mashtree")
def phylogenetic_tree_pipeline(args):

    start_time = time.time()

    from sr_amr.phylogeny_tree import mash_preprocessor, mash_distance_runner

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

    mash_distance_runner(mash_output, mash_temp, env_name="alpar-mashtree")

    print("Phylogenetic tree pipeline is finished, results can be found in the ", mash_output)

    if not args.keep_temp_files:
        print("Removing temp folder...")
        temp_folder_remover(mash_temp)

    end_time = time.time()

    print(time_function(start_time, end_time))


def fully_automated_pipeline(args):

    start_time = time.time()

    from sr_amr.full_automatix import automatix_runner

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
            if algorithm not in accepted_ml_algorithms:
                print("Error: ML algorithm is not accepted.")
                sys.exit(1)
    automatix_runner(args)

    end_time = time.time()

    print(time_function(start_time, end_time))


def structman_pipeline(args):

    from sr_amr.structman import structman_input_creator

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

    ensure_conda_env("alpar-snippy")
    if args.use_bakta:
        ensure_conda_env("alpar-bakta")
    else:
        ensure_conda_env("alpar-prokka")

    if not args.no_gene_presence_absence:
        if args.use_panaroo:
            ensure_conda_env("alpar-panaroo")
        else:
            ensure_conda_env("alpar-cdhit")

    from sr_amr.binary_tables import check_and_download_bakta_db, random_name_giver, prokka_create_database, snippy_processed_file_creator, binary_table_creator, annotation_file_from_snippy, panaroo_input_creator, panaroo_runner, binary_mutation_table_gpa_information_adder_panaroo, cdhit_preprocessor, cdhit_runner, gene_presence_absence_file_creator, binary_mutation_table_gpa_information_adder
    from sr_amr.ml import prps_ml_preprecessor
    from sr_amr.prediction import process_data_for_prediction, predict, equalize_columns

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
    gene_annotation_flag = True
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
        gene_annotation_flag = False
    if args.use_panaroo:
        panaroo_flag = True

    if args.verbosity > 3:
        print("Will run snippy: ", snippy_flag)
        print("Will run prokka: ", gene_annotation_flag)
        print("Will run panaroo: ", panaroo_flag)
        print("Will run gene presence absence: ", gpa_flag)

    # Create the output folder
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(os.path.join(args.output, "snippy")):
        os.makedirs(os.path.join(args.output, "snippy"), exist_ok=True)
    if args.use_panaroo and not os.path.exists(os.path.join(args.output, "panaroo")):
        os.makedirs(os.path.join(args.output, "panaroo"), exist_ok=True)
    if not os.path.exists(os.path.join(args.output, "cd-hit")):
        os.makedirs(os.path.join(args.output, "cd-hit"), exist_ok=True)
    if args.use_bakta and not os.path.exists(os.path.join(args.output, "bakta")):
        os.makedirs(os.path.join(args.output, "bakta"), exist_ok=True)
    if not args.use_bakta and not os.path.exists(os.path.join(args.output, "prokka")):
        os.makedirs(os.path.join(args.output, "prokka"), exist_ok=True)

    snippy_output = os.path.join(args.output, "snippy")
    prokka_output = os.path.join(args.output, "prokka")
    panaroo_output = os.path.join(args.output, "panaroo")
    bakta_output = os.path.join(args.output, "bakta")
    cd_hit_output = os.path.join(args.output, "cd-hit")

    # Create the temp folder for the panaroo input
    if args.use_panaroo and not os.path.exists(os.path.join(args.temp, "panaroo")):
        os.makedirs(os.path.join(args.temp, "panaroo"), exist_ok=True)
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
            args.custom_database[0], args.custom_database[1], args.temp, args.threads, args.ram, env_name="alpar-cdhit")
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

        annotation_output_to_use = bakta_output if args.use_bakta else prokka_output

        params = [(strain, random_names, snippy_output, annotation_output_to_use,
                   args, snippy_flag, gene_annotation_flag) for strain in strain_list]

        if args.custom_database:
            params = [(strain, random_names, snippy_output, annotation_output_to_use, args,
                       snippy_flag, gene_annotation_flag, args.custom_database[1]) for strain in strain_list]

        # with multiprocessing.Pool(num_parallel_tasks) as pool:
        #     pool.starmap(run_snippy_and_prokka, params)

        strains_to_be_processed = []

        prokka_output_strains = os.listdir(prokka_output)
        snippy_output_strains = os.listdir(snippy_output)

        strains_to_be_skiped = []

        for strain in random_names.keys():
            if gene_annotation_flag and snippy_flag:
                if random_names[strain] in prokka_output_strains and random_names[strain] in snippy_output_strains:
                    strains_to_be_processed.append(random_names[strain])
                else:
                    strains_to_be_skiped.append(random_names[strain])
            elif gene_annotation_flag and not snippy_flag:
                if random_names[strain] in prokka_output_strains:
                    strains_to_be_processed.append(random_names[strain])
                else:
                    strains_to_be_skiped.append(random_names[strain])
            elif snippy_flag and not gene_annotation_flag:
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
            args.output, "binary_mutation_table.tsv"), args.threads, strains_to_be_processed, args.temp)

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
                    args.temp, "panaroo_log.txt"), args.threads, env_name="alpar-panaroo")

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
                cdhit_runner(os.path.join(args.temp, "cdhit", "combined_proteins.faa"), os.path.join(args.output, "cd-hit", "cdhit_output.txt"), n_cpu=args.threads, env_name="alpar-cdhit")

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

