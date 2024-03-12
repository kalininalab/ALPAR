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
import seaborn as sns
import matplotlib.pyplot as plt
import math
import logging

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


def prokka_create_database(faa_file, genus_name, temp_folder, cpus = 1, memory = 4):

    output_dir = os.path.join(temp_folder, genus_name)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    script_command = f"cd-hit -i {faa_file} -o {output_dir}/{genus_name} -T {cpus} -M {memory*1024} -g 1 -s 0.8 -c 0.9"

    print(script_command)
    
    p = subprocess.Popen(script_command, shell=True)

    p.wait()

    if os.path.exists(f"{output_dir}/{genus_name}.bak.clstr"):
        os.remove(f"{output_dir}/{genus_name}.bak.clstr")

    if os.path.exists(f"{output_dir}/{genus_name}.clstr"):
        os.remove(f"{output_dir}/{genus_name}.clstr")

    script_command = f"makeblastdb -dbtype prot -in {genus_name}"

    p = subprocess.Popen(script_command, cwd=output_dir, shell=True)

    p.wait()

    logging.basicConfig(filename=f"{os.path.join(temp_folder, 'db_path.txt')}", level=logging.INFO)

    command = "prokka --listdb"

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    for line in iter(process.stdout.readline, b''):
        line = line.decode('utf-8').strip()
        logging.info(line)

    process.communicate()

    process.wait()

    with open (f"{os.path.join(temp_folder, 'db_path.txt')}") as infile:
        
        lines = infile.readlines()

        for line in lines:
            if "Looking for databases in:" in line:
                splitted = line.split(":")
                db_path = splitted[-1].strip()
    
    if db_path == "":
        print("Database path could not be found, please check the log file for more information.")
        sys.exit(1)

    for file in os.listdir(output_dir):

        shutil.copy2(os.path.join(output_dir, file), os.path.join(db_path, "genus"))


def prokka_runner(input, strain_random_name, output, reference, log_file, cpus = 1, parallel_run = False, custom_db = None):

    """

    input (str): Path to the input file
    output (str): Path to the output directory
    reference (str): Path to the reference file
    temp (str): Path to the temporary directory
    cpus (int): Number of cpus to use
    parallel_run (bool): If True, run prokka in parallel mode

    """

    main_path = os.getcwd()

    if custom_db == None:
        run_command = f"prokka --cpus {cpus} --outdir {output}/{strain_random_name} --proteins {reference} --force {input} >> {log_file} 2>&1"
    else:
        run_command = f"prokka --cpus {cpus} --outdir {output}/{strain_random_name} --proteins {reference} --usegenus --genus {custom_db} --force {input} >> {log_file} 2>&1"

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


def panaroo_runner(panaroo_input_folder, panaroo_output_folder, log_file, cpus):

    # Panaroo run fails if given threads are more than 16, issue has been raised their github
    
    # if cpus <= 16:
    #     cpus = cpus

    # else:
    #     cpus = 16 

    run_command = f"panaroo -i {panaroo_input_folder}/*.gff -o {panaroo_output_folder} --clean-mode strict -t {cpus} >> {log_file} 2>&1"
    
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
                    try:
                        binary_mutation_table_gpa_dict[strain_index_dict[cnt]][gene] = 1
                        cnt += 1
                    except:
                        cnt += 1
                    
                else:
                    try:
                        binary_mutation_table_gpa_dict[strain_index_dict[cnt]][gene] = 0
                        cnt += 1
                    except:
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
        if not antibiotic.startswith("."):
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


def panacota_pre_processor(list_file, temp_folder, output_folder, random_names_dict):

    random_names_will_be_used = False

    if random_names_dict != None:
        random_names_will_be_used = True
        random_names = {}
        random_names_to_strains = {}
        with open(random_names_dict, "r") as infile:
            lines = infile.readlines()
            for line in lines:
                splitted = line.split("\t")
                random_names[splitted[0].strip()] = splitted[1].strip()
                random_names_to_strains[splitted[1].strip()] = splitted[0].strip()
    
    already_copied = []

    with open(list_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            strain_path = line.strip()
            if random_names_will_be_used:
                if not random_names[os.path.splitext(line.split('/')[-1].strip())[0]] in already_copied:
                    shutil.copy2(strain_path, f"{temp_folder}/{random_names[os.path.splitext(line.split('/')[-1].strip())[0]]}")
                    already_copied.append(random_names[os.path.splitext(line.split('/')[-1].strip())[0]])
            else:
                if not line.split('/')[-1].strip() in already_copied:
                    shutil.copy2(strain_path, f"{temp_folder}/{line.split('/')[-1].strip()}")
                    already_copied.append(line.split('/')[-1].strip())
                
    the_names_will_be_skipped = []
    
    if os.path.exists(f"{os.path.dirname(output_folder)}/snippy"):

        snippy_dir = os.listdir(f"{os.path.dirname(output_folder)}/snippy")

        for strain in os.listdir(temp_folder):
            if strain not in snippy_dir:
                the_names_will_be_skipped.append(strain)
    
    print(f"The following strains will be skipped: {the_names_will_be_skipped}")
   
    with open(os.path.join(output_folder, "panacota_input.lst"), "w") as ofile:
        for strain in os.listdir(temp_folder):
            if not strain in the_names_will_be_skipped:
                ofile.write(f"{strain}\n")

def replace_values(line, replacements):
    for key, value in replacements.items():
        line = line.replace(key, value)
    return line


def panacota_post_processor(panacota_output_folder, run_name, output_folder, type="nucl"):

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

    shutil.copyfile(f"{panacota_output_folder}/phylogenetic_tree.newick", f"{output_folder}/phylogenetic_tree.newick")
    
            
# TODO PanACoTA needs all the files in one folder, need to process them before running PanACoTA
def panacota_pipeline_runner(list_file, dbpath, output_directory, run_name, n_cores, log_file, type="nucl", mode=1, min_seq_id=0.8):
    
    pc_annotate_command = f"PanACoTA annotate -l {list_file} -d {dbpath} -r {output_directory}/annotate_out -n {run_name} --threads {n_cores} >> {log_file} 2>&1"
    
    pc_pangenome_command = f"PanACoTA pangenome -l {output_directory}/annotate_out/LSTINFO-{list_file.split('/')[-1]} -d {output_directory}/annotate_out/Proteins/ -o {output_directory}/pangenome_out/ -n {run_name} --threads {n_cores} >> {log_file} 2>&1"

    if n_cores > 1:
        pc_corepers_command = f"PanACoTA corepers -p {output_directory}/pangenome_out/PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}-th{n_cores}.lst -o {output_directory}/corepers_out/  >> {log_file} 2>&1"
        pc_align_command =  f"PanACoTA align -c {output_directory}/corepers_out/PersGenome_PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}-th{n_cores}.lst-all_1.lst -l {output_directory}/annotate_out/LSTINFO-{list_file.split('/')[-1]} -n {run_name} -d {output_directory}/annotate_out/ -o {output_directory}/align_out --threads {n_cores} >> {log_file} 2>&1"
    
    else:
        pc_corepers_command = f"PanACoTA corepers -p {output_directory}/pangenome_out/PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}.lst -o {output_directory}/corepers_out/  >> {log_file} 2>&1"
        pc_align_command =  f"PanACoTA align -c {output_directory}/corepers_out/PersGenome_PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}.lst-all_1.lst -l {output_directory}/annotate_out/LSTINFO-{list_file.split('/')[-1]} -n {run_name} -d {output_directory}/annotate_out/ -o {output_directory}/align_out --threads {n_cores} >> {log_file} 2>&1"

    pc_tree_command = f"PanACoTA tree -a {output_directory}/align_out/Phylo-{run_name}/{run_name}.{type}.grp.aln -o {output_directory}/tree/ --threads {n_cores} >> {log_file} 2>&1"

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


def temp_folder_remover(path):
    for i in range(5):
        try:
            shutil.rmtree(path, ignore_errors=True)
            break 
        except OSError:
            time.sleep(1)
    else:
        print(f"Failed to delete {path}")


def binary_table_threshold_with_percentage(binary_table, output_folder, threshold_percentage):

    df = pd.read_csv(binary_table, sep="\t")

    print(f"Number of mutations in the table: {len(df.columns[1:])}")

    cols_to_be_dropped = []

    amount_of_strains = int(df.shape[0])

    threshold_value = (threshold_percentage / 100) * amount_of_strains

    for col in df.columns[1:]:
        if int(df[col].sum()) <= int(threshold_value):
            cols_to_be_dropped.append(col)

    print(f"Number of mutations to be dropped: {len(cols_to_be_dropped)}")

    dropped_df = df.drop(cols_to_be_dropped, axis=1)

    print(f"Number of mutations in the table after dropping: {len(dropped_df.columns[1:])}")
            
    dropped_df.to_csv(os.path.join(output_folder, f"binary_mutation_table_threshold_{threshold_percentage}_percent.tsv"), sep="\t", index=False)


def mash_preprocessor(strains_text_file, temp_folder, random_names_dict):

    os.mkdir(os.path.join(temp_folder, "fasta_files"))

    fasta_files_folder = os.path.join(temp_folder, "fasta_files")

    if random_names_dict != None:
        random_names_will_be_used = True
        random_names = {}
        with open(random_names_dict, "r") as infile:
            lines = infile.readlines()
            for line in lines:
                splitted = line.split("\t")
                random_names[splitted[0].strip()] = splitted[1].strip()

    with open (strains_text_file) as infile:
        lines = infile.readlines()
        for line in lines:
            strain_path = line.strip()
            if random_names_will_be_used:
                shutil.copy2(strain_path, f"{fasta_files_folder}/{random_names[os.path.splitext(line.split('/')[-1].strip())[0]]}.fna")
            else:
                shutil.copy2(strain_path, f"{fasta_files_folder}/{line.split('/')[-1].strip()}")


def mash_distance_runner(output_folder, temp_folder):

    # Get the absolute path of the current file
    temp_folder_abs_path = os.path.abspath(temp_folder)

    print(temp_folder_abs_path)

    mash_tree_command = f"mashtree {temp_folder_abs_path}/fasta_files/* > {output_folder}/phylogenetic_tree.tree"
    
    os.system(mash_tree_command)


def time_function(start_time, end_time):

    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"Elapsed time: {int(hours)} hours, {int(minutes)} minutes, {seconds:.2f} seconds"
