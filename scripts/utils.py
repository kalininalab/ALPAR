#!/usr/bin/env python3

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

# Get the path of the script
PATH_OF_SCRIPT = pathlib.Path(__file__).parent.resolve()


def is_tool_installed(tool_name):
    # 'which' for Unix-based systems, 'where' for Windows
    command = f"which {tool_name}" if os.name != 'nt' else f"where {tool_name}"
    return os.system(command) == 0


def check_contigs(input):
    """
    Check if the input file is a fasta file and if it has more than 1 contig
    """
    
    with open(input, 'r') as f:
        contigs = 0
        for line in f.readlines():
            if line.startswith('>'):
                contigs += 1
        if contigs > 1:
            return True
        else:
            return False

# Maybe can be improved by using parallel run
def snippy_runner(input, strain_random_name, output, reference, log_file, cpus = 1, memory = 4, parallel_run = False):

    """
    
    input (str): Path to the input file
    output (str): Path to the output directory
    reference (str): Path to the reference file
    temp (str): Path to the temporary directory
    cpus (int): Number of cpus to use
    memory (int): Amount of memory to use
    parallel_run (bool): If True, run snippy in parallel mode

    """

    #main_path = os.getcwd()

    # Create the output directory

    # if not os.path.exists(f"{output}/{strain_random_name}"):
    #     os.mkdir(f"{output}/{strain_random_name}")

    contigs = check_contigs(input)

    run_command = f"snippy --cpus {cpus} --ram {memory} --outdir {output}/{strain_random_name} --reference {reference} --force"

    if contigs == True:
        run_command = f"{run_command} --ctgs"
    

    run_command = f"{run_command} {input} >> {log_file} 2>&1"

    # Run snippy
    os.system(run_command)


# Maybe can be improved by using parallel run
def prokka_runner(input, strain_random_name, output, reference, log_file, cpus = 1, parallel_run = False):

    """

    input (str): Path to the input file
    output (str): Path to the output directory
    reference (str): Path to the reference file
    temp (str): Path to the temporary directory
    cpus (int): Number of cpus to use
    parallel_run (bool): If True, run prokka in parallel mode

    """

    main_path = os.getcwd()

    # Create the output directory

    # if not os.path.exists(f"{output}/{strain_random_name}"):
    #     os.mkdir(f"{output}/{strain_random_name}")

    # Check extension for reference file
    #reference_path = f"{PATH_OF_SCRIPT}/reference_files/{bacterium}.fasta"

    run_command = f"prokka --cpus {cpus} --outdir {output}/{strain_random_name} --proteins {reference} --force {input} >> {log_file} 2>&1"


    # Run prokka
    os.system(run_command)


def generate_random_key():
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))

random_names = {}

def random_name_giver(strains_text_file, random_name_file_output):

    with open(strains_text_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            strain = os.path.splitext(line.split("/")[-1].strip())[0]
            random_key = generate_random_key()
            # make sure it is unique
            while random_key in random_names.keys():
                random_key = generate_random_key()
            
            random_names[strain] = random_key

    with open(random_name_file_output, "w") as ofile:
        for key in random_names.keys():
            ofile.write(f"{key}\t{random_names[key]}\n")

    return random_names


def panaroo_input_creator(random_names_txt, prokka_output_folder, temp_folder, strains_to_be_processed):

    with open(random_names_txt, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            splitted = line.split("\t")
            given_random_name = splitted[1].strip()
            if given_random_name in strains_to_be_processed:
                for prokka_output_file in os.listdir(f"{prokka_output_folder}/{given_random_name}"):
                    if prokka_output_file.endswith(".gff"):
                        shutil.copyfile(f"{prokka_output_folder}/{given_random_name}/{prokka_output_file}", f"{temp_folder}/{given_random_name}.gff")


def panaroo_runner(panaroo_input_folder, panaroo_output_folder, log_file):

    run_command = f"panaroo -i {panaroo_input_folder}/*.gff -o {panaroo_output_folder}  --clean-mode strict >> {log_file} 2>&1"
    
    os.system(run_command)


def read_vcf_and_return_snp_class_list(vcf_path):

    with open(f"{vcf_path}") as infile:
        lines = infile.readlines()

    snp_list = []
    
    for line in lines:
        if not line.startswith("#"):
            splitted = line.split("\t")
            pos = splitted[1]
            ref = splitted[3]
            alt = splitted[4]
            splitted_info = splitted[7].split(';')
            for p in splitted_info:
                if p.startswith("TYPE"):
                    splitted_info2 = p.split('=')
                    mut_type = splitted_info2[1]

            temp_str = ref + ":" + alt
            temp_mutation_name = f"'{pos}','{temp_str}','{mut_type}'"
            #temp_tuple = (pos,temp_str,mut_type)

            snp_list.append(temp_mutation_name)

    return snp_list


def process_folder(snippy_out_path, folder):
    vcf_file_path = f"{snippy_out_path}/{folder}/snps.vcf"
    snp_list = read_vcf_and_return_snp_class_list(vcf_file_path)
    return folder, snp_list


def read_vcf_files_and_store_data(snippy_out_path, n_jobs, strains_to_be_processed):
    snp_combined_list = []
    mutation_dict = {}
    folders = [folder for folder in os.listdir(snippy_out_path) if os.path.isfile(f"{snippy_out_path}/{folder}/snps.vcf")]

    folders = [strain for strain in folders if strain in strains_to_be_processed]

    results = Parallel(n_jobs)(
        delayed(process_folder)(snippy_out_path, folder) for folder in folders
    )

    for folder, snp_list in results:
        mutation_dict[folder] = snp_list
        snp_combined_list = snp_combined_list + snp_list
    
    snp_combined_set = set(snp_combined_list)
    
    return mutation_dict, snp_combined_set


def temp_dict_creator(combined_set):
    temp_dict = {}

    for i in combined_set:
        temp_dict[i] = 0

    return temp_dict


def mutation_presence_absence_dict_creator(mut_dict, temp_dict, n_jobs):
    mutation_presence_absence_dict = {}

    results = Parallel(n_jobs)(
        delayed(strain_presence_absence)(mut_dict[strain], temp_dict) for strain in mut_dict.keys()
    )

    for idx, strain in enumerate(mut_dict.keys()):
        mutation_presence_absence_dict[strain] = results[idx]

    return mutation_presence_absence_dict


def strain_presence_absence(list_of_mutations, temp_dict):

    mut_presence_absence_dict = copy.deepcopy(temp_dict)

    for mutation in list_of_mutations:
        mut_presence_absence_dict[mutation] = 1
    
    return mut_presence_absence_dict


def binary_table_creator (input_folder, output_file, number_of_cores, strains_to_be_processed):

    number_of_cores = int(number_of_cores)

    mut_dict, snp_combined_set = read_vcf_files_and_store_data(f"{input_folder}", number_of_cores, strains_to_be_processed)

    temp_dict = temp_dict_creator(snp_combined_set)

    mutation_presence_absence_dict = mutation_presence_absence_dict_creator(mut_dict, temp_dict, number_of_cores)

    df = pd.DataFrame.from_dict(mutation_presence_absence_dict, orient='index')

    df.to_csv(f"{output_file}", sep="\t")


def binary_mutation_table_gpa_information_adder(binary_mutation_table, panaroo_output_gpa, binary_mutation_table_with_gpa_information):

    binary_mutation_table_df = pd.read_csv(binary_mutation_table, sep="\t", index_col=0)

    binary_mutation_table_df_dict = binary_mutation_table_df.to_dict(orient='index')

    binary_mutation_table_gpa_dict = copy.deepcopy(binary_mutation_table_df_dict)

    strain_index_dict = {}

    with open (panaroo_output_gpa) as infile:
        lines = infile.readlines()

        index_line = lines[0]

        index_line_split = index_line.split(",")

        for strain_idx in range(len(index_line_split)):
            
            if strain_idx < 3:
                continue

            strain_index_dict[strain_idx] = index_line_split[strain_idx].strip()


        for line in lines[1:]:
            splitted = line.split(",")

            gene = splitted[0].strip()

            cnt = 3

            for ocome in splitted[3:]:
                
                if ocome.strip() != "":
                    binary_mutation_table_gpa_dict[strain_index_dict[cnt]][gene] = 1
                    cnt += 1
                
                else:
                    binary_mutation_table_gpa_dict[strain_index_dict[cnt]][gene] = 0
                    cnt += 1


    binary_mutation_table_gpa_df = pd.DataFrame.from_dict(binary_mutation_table_gpa_dict, orient='index')

    binary_mutation_table_gpa_df.to_csv(binary_mutation_table_with_gpa_information, sep="\t")

        
def phenotype_dataframe_creator(data_folder_path, output_file, random_names_dict):
    list_of_antibiotics = os.listdir(data_folder_path)

    phenotype_dict = {}

    strain_phenotype_dict = {}

    for antibiotic in list_of_antibiotics:
        phenotype_dict[antibiotic] = 2

    for antibiotic in list_of_antibiotics:
        res_strains = os.listdir(f"{data_folder_path}/{antibiotic}/Resistant")

        sus_strains = os.listdir(f"{data_folder_path}/{antibiotic}/Susceptible")

        for strain in res_strains:
            if strain.endswith(".fna"):
                if not random_names_dict[strain[:-4]] in strain_phenotype_dict.keys():
                    strain_phenotype_dict[random_names_dict[strain[:-4]]] = copy.deepcopy(phenotype_dict)
                
                strain_phenotype_dict[random_names_dict[strain[:-4]]][antibiotic] = 1
        
        for strain in sus_strains:
            if strain.endswith(".fna"):
                if not random_names_dict[strain[:-4]] in strain_phenotype_dict.keys():
                    strain_phenotype_dict[random_names_dict[strain[:-4]]] = copy.deepcopy(phenotype_dict)
                
                strain_phenotype_dict[random_names_dict[strain[:-4]]][antibiotic] = 0

    df = pd.DataFrame.from_dict(strain_phenotype_dict, orient="index")

    df.to_csv(output_file, sep="\t")


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


def pyseer_similarity_matrix_creator(phylogenetic_tree, output_file):
    
    script_command = f"python {PATH_OF_SCRIPT}/phylogeny_distance.py --lmm {phylogenetic_tree} > {output_file}"

    os.system(script_command)


def pyseer_runner(genotype_file_path, phenotype_file_path, similarity_matrix, output_file_directory):

    phenotypes = os.listdir(f"{phenotype_file_path}")

    for phenotype in phenotypes:
        script_command = f"pyseer --lmm --phenotypes {phenotype_file_path}/{phenotype} --pres {genotype_file_path} --similarity {similarity_matrix} > {output_file_directory}/{phenotype}.tsv"
                
        if not os.path.exists(f"{output_file_directory}"):
            os.mkdir(f"{output_file_directory}")

        os.system(script_command)


def panacota_pre_processor(list_file, temp_folder, output_folder, random_names_dict):

    random_names_will_be_used = False

    if random_names_dict != None:
        random_names_will_be_used = True
        random_names = {}
        with open(random_names_dict, "r") as infile:
            lines = infile.readlines()
            for line in lines:
                splitted = line.split("\t")
                random_names[splitted[0].strip()] = splitted[1].strip()

    with open(list_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            strain_path = line.strip()
            if random_names_will_be_used:
                shutil.copy2(strain_path, f"{temp_folder}/{random_names[os.path.splitext(line.split('/')[-1].strip())[0]]}")
            else:
                shutil.copy2(strain_path, f"{temp_folder}/{line.split('/')[-1].strip()}")

    with open(os.path.join(output_folder, "panacota_input.lst"), "w") as ofile:
        for strain in os.listdir(temp_folder):
            ofile.write(f"{strain}\n")

def replace_values(line, replacements):
    for key, value in replacements.items():
        line = line.replace(key, value)
    return line

def panacota_post_processor(panacota_output_folder, run_name, type="nucl"):

    new_name_to_origin_name = {}

    with open(os.path.join(panacota_output_folder, "annotate_out", "LSTINFO-panacota_input.lst")) as infile:
        lines = infile.readlines()
        for line in lines:
            splitted = line.split("\t")
            new_name_to_origin_name[splitted[0].strip()] = splitted[1].strip()

    with open(os.path.join(panacota_output_folder, "tree", f"{run_name}.{type}.grp.aln.iqtree_tree.treefile")) as infile:
        line = infile.readline()
        with open(os.path.join(panacota_output_folder, "phylogenetic_tree.newick"), "w") as ofile:
            ofile.write(replace_values(line, new_name_to_origin_name))
    
            
# TODO PanACoTA needs all the files in one folder, need to process them before running PanACoTA
def panacota_pipeline_runner(list_file, dbpath, output_directory, run_name, n_cores, log_file, type="nucl", mode=1, min_seq_id=0.8):
    
    pc_annotate_command = f"PanACoTA annotate -l {list_file} -d {dbpath} -r {output_directory}/annotate_out -n {run_name} --threads {n_cores} >> {log_file} 2>&1"
    
    pc_pangenome_command = f"PanACoTA pangenome -l {output_directory}/annotate_out/LSTINFO-{list_file.split('/')[-1]} -d {output_directory}/annotate_out/Proteins/ -o {output_directory}/pangenome_out/ -n {run_name} --threads {n_cores}  >> {log_file} 2>&1"

    if n_cores > 1:
        pc_corepers_command = f"PanACoTA corepers -p {output_directory}/pangenome_out/PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}-th{n_cores}.lst -o {output_directory}/corepers_out/  >> {log_file} 2>&1"
        pc_align_command =  f"PanACoTA align -c {output_directory}/corepers_out/PersGenome_PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}-th{n_cores}.lst-all_1.lst -l {output_directory}/annotate_out/LSTINFO-{list_file.split('/')[-1]} -n {run_name} -d {output_directory}/annotate_out/ -o {output_directory}/align_out --threads {n_cores}  >> {log_file} 2>&1"
    
    else:
        pc_corepers_command = f"PanACoTA corepers -p {output_directory}/pangenome_out/PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}.lst -o {output_directory}/corepers_out/  >> {log_file} 2>&1"
        pc_align_command =  f"PanACoTA align -c {output_directory}/corepers_out/PersGenome_PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}.lst-all_1.lst -l {output_directory}/annotate_out/LSTINFO-{list_file.split('/')[-1]} -n {run_name} -d {output_directory}/annotate_out/ -o {output_directory}/align_out --threads {n_cores}  >> {log_file} 2>&1"

    pc_tree_command = f"PanACoTA tree -a {output_directory}/align_out/Phylo-{run_name}/{run_name}.{type}.grp.aln -o {output_directory}/tree/ --threads {n_cores}  >> {log_file} 2>&1"


    print(f"Running PanACoTA annotate...")
    os.system(pc_annotate_command)
    print(f"Running PanACoTA pangenome...")
    os.system(pc_pangenome_command)
    print(f"Running PanACoTA corepers...")
    os.system(pc_corepers_command)
    print(f"Running PanACoTA align...")
    # I don't know why but it needs to be run twice, otherwise it gives an error on some occasions, probably PanACoTA side issue
    for _ in range(2):
        os.system(pc_align_command)
    print(f"Running PanACoTA tree...")
    os.system(pc_tree_command)


def temp_folder_remover(temp_folder):
    shutil.rmtree(temp_folder)