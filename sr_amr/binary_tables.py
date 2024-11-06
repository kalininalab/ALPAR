#!/usr/bin/env python3

import os
import sys
import subprocess
import random
import string
import shutil
import pandas as pd
from joblib import Parallel, delayed
import copy
import logging
import csv
from collections import defaultdict
import re
from sr_amr.gbk_parser import gbk_parser

import warnings

warnings.filterwarnings("ignore")

csv.field_size_limit(sys.maxsize)

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


def snippy_runner(input, strain_random_name, output, reference, log_file, cpus=1, memory=8, parallel_run=False):
    """
    input (str): Path to the input file
    output (str): Path to the output directory
    reference (str): Path to the reference file
    temp (str): Path to the temporary directory
    cpus (int): Number of cpus to use
    memory (int): Amount of memory to use
    parallel_run (bool): If True, run snippy in parallel mode

    """

    #contigs = check_contigs(input)

    contigs = True

    run_command = f"snippy --cpus {cpus} --ram {memory} --outdir {output}/{strain_random_name} --reference {reference} --force"

    if contigs == True:
        run_command = f"{run_command} --ctgs"

    run_command = f"{run_command} {input} >> {log_file} 2>&1"

    os.system(run_command)  


def prokka_create_database(faa_file, genus_name, temp_folder, cpus=1, memory=8):

    output_dir = os.path.join(temp_folder, genus_name)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    log_file = os.path.join(temp_folder, "prokka_db.log")

    script_command = f"cd-hit -i {faa_file} -o {output_dir}/{genus_name} -T {cpus} -M {memory*1024} -g 1 -s 0.8 -c 0.9 >> {log_file} 2>&1"

    p = subprocess.Popen(script_command, shell=True)

    p.wait()

    if os.path.exists(f"{output_dir}/{genus_name}.bak.clstr"):
        os.remove(f"{output_dir}/{genus_name}.bak.clstr")

    if os.path.exists(f"{output_dir}/{genus_name}.clstr"):
        os.remove(f"{output_dir}/{genus_name}.clstr")

    script_command = f"makeblastdb -dbtype prot -in {genus_name}"

    p = subprocess.Popen(script_command, cwd=output_dir, shell=True)

    p.wait()

    logging.basicConfig(
        filename=f"{os.path.join(temp_folder, 'db_path.txt')}", level=logging.INFO)

    command = "prokka --listdb"

    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    for line in iter(process.stdout.readline, b''):
        line = line.decode('utf-8').strip()
        logging.info(line)

    process.communicate()

    process.wait()

    with open(f"{os.path.join(temp_folder, 'db_path.txt')}") as infile:

        lines = infile.readlines()

        for line in lines:
            if "Looking for databases in:" in line:
                splitted = line.split(":")
                db_path = splitted[-1].strip()

    if db_path == "":
        print("Database path could not be found, please check the log file for more information.")
        sys.exit(1)

    for file in os.listdir(output_dir):

        shutil.copy2(os.path.join(output_dir, file),
                     os.path.join(db_path, "genus"))


def prokka_runner(input, strain_random_name, output, reference, log_file, cpus=1, parallel_run=False, custom_db=None):
    """

    input (str): Path to the input file
    output (str): Path to the output directory
    reference (str): Path to the reference file
    temp (str): Path to the temporary directory
    cpus (int): Number of cpus to use
    parallel_run (bool): If True, run prokka in parallel mode

    """

    if custom_db == None:
        run_command = f"prokka --cpus {cpus} --outdir {output}/{strain_random_name} --proteins {reference} --compliant --force {input} >> {log_file} 2>&1"
    else:
        run_command = f"prokka --cpus {cpus} --outdir {output}/{strain_random_name} --proteins {reference} --usegenus --genus {custom_db} --compliant --force {input} >> {log_file} 2>&1"

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
                        shutil.copyfile(
                            f"{prokka_output_folder}/{given_random_name}/{prokka_output_file}", f"{temp_folder}/{given_random_name}.gff")


def panaroo_runner(panaroo_input_folder, panaroo_output_folder, log_file, cpus):

    # Panaroo fails with more than 30 cpus
    if cpus >= 32:
        cpus = 30

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

            snp_list.append(temp_mutation_name)

    return snp_list


def process_folder(snippy_out_path, folder):
    vcf_file_path = f"{snippy_out_path}/{folder}/snps.vcf"
    snp_list = read_vcf_and_return_snp_class_list(vcf_file_path)
    return folder, snp_list


def read_vcf_files_and_store_data(snippy_out_path, n_jobs, strains_to_be_processed):
    snp_combined_list = []
    mutation_dict = {}
    folders = [folder for folder in os.listdir(snippy_out_path) if os.path.isfile(
        f"{snippy_out_path}/{folder}/snps.vcf")]

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


def binary_table_creator(input_folder, output_file, number_of_cores, strains_to_be_processed):
    number_of_cores = int(number_of_cores)

    mut_dict, snp_combined_set = read_vcf_files_and_store_data(
        f"{input_folder}", number_of_cores, strains_to_be_processed)

    temp_dict = temp_dict_creator(snp_combined_set)

    mutation_presence_absence_dict = mutation_presence_absence_dict_creator(
        mut_dict, temp_dict, number_of_cores)

    headers = ['Strain'] + list(snp_combined_set) 

    data_rows = []
    for strain, mutations in mutation_presence_absence_dict.items():
        row = [strain] + [mutations.get(mutation, 0) for mutation in snp_combined_set]
        data_rows.append(row)

    with open(f"{output_file}", "w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(data_rows)


def binary_mutation_table_gpa_information_adder_panaroo(binary_mutation_table, panaroo_output_gpa, binary_mutation_table_with_gpa_information):

    with open(binary_mutation_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        # Create a mapping of header names to their indices
        header_indices = {name: index for index, name in enumerate(headers[1:], start=1)}
        binary_table_dict = {}
        for row in reader:
            strain = row[0]
            mutations = row[1:]
            binary_table_dict[strain] = {mutation_name: mutations[header_indices[mutation_name]-1] for mutation_name in headers[1:]}

    binary_mutation_table_gpa_dict = copy.deepcopy(binary_table_dict)

    strain_index_dict = {}

    with open(panaroo_output_gpa) as infile:
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
                        binary_mutation_table_gpa_dict[strain_index_dict[cnt]][gene] = "1"
                        cnt += 1
                    except:
                        cnt += 1
                else:
                    try:
                        binary_mutation_table_gpa_dict[strain_index_dict[cnt]][gene] = "0"
                        cnt += 1
                    except:
                        cnt += 1

    with open(binary_mutation_table_with_gpa_information, 'w') as ofile:
        headers = ['Strain'] + list(next(iter(binary_mutation_table_gpa_dict.values())).keys())
        ofile.write('\t'.join(headers) + '\n') 
        
        for strain, mutations in binary_mutation_table_gpa_dict.items():
            row = [strain] + [mutations[mutation] for mutation in headers[1:]] 
            ofile.write('\t'.join(row) + '\n')

    table_binary_maker(binary_mutation_table_with_gpa_information)


def binary_mutation_table_gpa_information_adder(binary_mutation_table, gpa_file, binary_mutation_table_with_gpa_information):

    with open(binary_mutation_table, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        # Create a mapping of header names to their indices
        header_indices = {name: index for index, name in enumerate(headers[1:], start=1)}
        binary_table_dict = {}
        for row in reader:
            strain = row[0]
            mutations = row[1:]
            binary_table_dict[strain] = {mutation_name: mutations[header_indices[mutation_name]-1] for mutation_name in headers[1:]}

    binary_mutation_table_gpa_dict = copy.deepcopy(binary_table_dict)

    strain_index_dict = {}

    with open(gpa_file) as infile:
        lines = infile.readlines()
        index_line = lines[0]
        index_line_split = index_line.split(",")
        for strain_idx in range(len(index_line_split)):
            if strain_idx < 1:
                continue
            strain_index_dict[strain_idx] = index_line_split[strain_idx].strip()
        for line in lines[1:]:
            splitted = line.split(",")
            gene = splitted[0].strip()
            cnt = 1
            for ocome in splitted[1:]:
                if ocome.strip() == 1 or ocome.strip() == '1':
                    try:
                        binary_mutation_table_gpa_dict[strain_index_dict[cnt]][gene] = "1"
                        cnt += 1
                    except:
                        cnt += 1
                elif ocome.strip() == 0 or ocome.strip() == '0':
                    try:
                        binary_mutation_table_gpa_dict[strain_index_dict[cnt]][gene] = "0"
                        cnt += 1
                    except:
                        cnt += 1

    with open(binary_mutation_table_with_gpa_information, 'w') as ofile:
        headers = ['Strain'] + list(next(iter(binary_mutation_table_gpa_dict.values())).keys())
        ofile.write('\t'.join(headers) + '\n') 
        
        for strain, mutations in binary_mutation_table_gpa_dict.items():
            row = [strain] + [mutations[mutation] for mutation in headers[1:]] 
            ofile.write('\t'.join(row) + '\n')

    table_binary_maker(binary_mutation_table_with_gpa_information)
    

def phenotype_dataframe_creator(data_folder_path, output_file, random_names_dict):
    list_of_antibiotics = os.listdir(data_folder_path)

    phenotype_dict = {}

    strain_phenotype_dict = {}

    for antibiotic in list_of_antibiotics:
        phenotype_dict[antibiotic] = 2

    for antibiotic in list_of_antibiotics:
        if not antibiotic.startswith("."):
            res_strains = os.listdir(
                f"{data_folder_path}/{antibiotic}/Resistant")
            sus_strains = os.listdir(
                f"{data_folder_path}/{antibiotic}/Susceptible")
            for strain in res_strains:
                if strain.endswith(".fna"):
                    if not random_names_dict[strain[:-4]] in strain_phenotype_dict.keys():
                        strain_phenotype_dict[random_names_dict[strain[:-4]]
                                              ] = copy.deepcopy(phenotype_dict)
                    strain_phenotype_dict[random_names_dict[strain[:-4]]
                                          ][antibiotic] = 1
            for strain in sus_strains:
                if strain.endswith(".fna"):
                    if not random_names_dict[strain[:-4]] in strain_phenotype_dict.keys():
                        strain_phenotype_dict[random_names_dict[strain[:-4]]
                                              ] = copy.deepcopy(phenotype_dict)
                    strain_phenotype_dict[random_names_dict[strain[:-4]]
                                          ][antibiotic] = 0

    df = pd.DataFrame.from_dict(strain_phenotype_dict, orient="index")
    df.to_csv(output_file, sep="\t")


import csv

def phenotype_dataframe_creator_post_processor(genotype_data, phenotype_data):

    with open(genotype_data, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        genotype_header = next(reader)
        genotype_dict = {rows[0]: rows[1:] for rows in reader}

    with open(phenotype_data, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        phenotype_header = next(reader)
        phenotype_dict = {rows[0]: rows[1:] for rows in reader}

    strains_to_be_dropped = [strain for strain in phenotype_dict if strain not in genotype_dict]

    for strain in strains_to_be_dropped:
        del phenotype_dict[strain]

    with open(phenotype_data, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(phenotype_header)
        for key, value in phenotype_dict.items():
            writer.writerow([key] + value)


def snippy_processed_file_creator(snippy_output_folder, snippy_processed_text_file):

    processed_by_snippy = os.listdir(snippy_output_folder)

    with open(snippy_processed_text_file, "w") as ofile:
        for strain in processed_by_snippy:
            ofile.write(f"{strain}\n")


def annotation_file_from_snippy(snippy_output_folder, output_folder):

    annotation_dict = {}

    snippy_output_strains = os.listdir(snippy_output_folder)

    for strain in snippy_output_strains:
        if "snps.tab" in os.listdir(os.path.join(snippy_output_folder, strain)):
            with open(os.path.join(snippy_output_folder, strain, "snps.tab")) as infile:
                lines = infile.readlines()
            
            for line in lines[1:]:
                splitted = line.split("\t")
                mutation = f"'{splitted[1]}','{splitted[3]}:{splitted[4]}','{splitted[2]}'"
                if not mutation in annotation_dict.keys():
                    annotation_dict[mutation] = [splitted[10], splitted[12], splitted[13]]

    with open(os.path.join(output_folder, "mutations_annotations.tsv"), "w") as ofile:
        ofile.write(f"Mutation\tEFFECT\tGENE\tPRODUCT\n")
        for key in annotation_dict.keys():
            ofile.write(f"{key}\t{annotation_dict[key][0]}\t{annotation_dict[key][1]}\t{annotation_dict[key][2]}")

# After addition of gene presence absence with panaroo, some values become "0.0" or "1.0" instead of "0" or "1". This function changes them to "0" or "1"
def table_binary_maker(binary_mutation_table):
    with open(binary_mutation_table, "r") as infile:
        lines = infile.readlines()

    with open(binary_mutation_table, "w") as ofile:
        ofile.write(lines[0])
        for line in lines[1:]:
            line = line.replace("0.0", "0")
            line = line.replace("1.0", "1")

            ofile.write(line)

def replace_non_alphanumeric(string):
    # Replace every character that is not a letter or number with '_'
    return re.sub(r'\W', '_', string)

def combine_faa_files(temp_folder):
    command = f'cat {temp_folder}/*.faa > {temp_folder}/combined_proteins.faa'
    subprocess.run(command, shell=True, check=True)

def delete_unwanted_faa_files(temp_folder):
    for file_name in os.listdir(temp_folder):
        if file_name.endswith(".faa") and file_name != "combined_proteins.faa":
            os.remove(os.path.join(temp_folder, file_name))

def cdhit_protein_dictionary_creator(fasta_file):
    protein_dict = {}
    with open(fasta_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            if line.startswith(">"):
                protein_id = line.strip().split(" ")[0][1:]
                protein_names = line.strip().split(" ")[1:]
                protein_name = "_".join(protein_names)
                if protein_name.endswith("_"):
                    protein_name = protein_name[:-1]
                protein_dict[protein_id] = protein_name
    return protein_dict

def cdhit_protein_name_corrector(input_file, output_file):
    fasta_file_name = input_file.split("/")[-1].split(".")[0]
    with open(output_file, "w") as outfile:
        with open(input_file, "r") as infile:
            lines = infile.readlines()
            for line in lines:
                if line.startswith(">"):
                    protein_original_name = line[1:].strip()
                    protein_original_name = replace_non_alphanumeric(protein_original_name)
                    if protein_original_name.endswith("_"):
                        protein_original_name = protein_original_name[:-1]
                    outfile.write(f">{fasta_file_name}_{protein_original_name}\n")
                else:
                    outfile.write(line)

def cdhit_preprocessor(random_names_txt, prokka_output_folder, temp_folder, strains_to_be_processed):
    protein_positions = {}
    with open(random_names_txt, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            splitted = line.split("\t")
            given_random_name = splitted[1].strip()
            if given_random_name in strains_to_be_processed:
                for prokka_output_file in os.listdir(f"{os.path.join(prokka_output_folder, given_random_name)}"):
                    if prokka_output_file.endswith(".faa"):
                        if prokka_output_file.startswith("PROKKA"):
                            shutil.copyfile(os.path.join(f"{prokka_output_folder}", f"{given_random_name}", f"{prokka_output_file}"), os.path.join(f"{temp_folder}", f"{given_random_name}.faa"))
                            cdhit_protein_name_corrector(os.path.join(f"{temp_folder}", f"{given_random_name}.faa"), os.path.join(f"{temp_folder}", f"{given_random_name}_corrected.faa"))
                            os.remove(os.path.join(f"{temp_folder}", f"{given_random_name}.faa"))
                    if prokka_output_file.endswith(".gbk"):
                        if prokka_output_file.startswith("PROKKA"):
                            current_records = gbk_parser(given_random_name, os.path.join(f"{prokka_output_folder}", f"{given_random_name}", f"{prokka_output_file}"))
                            for record in current_records:
                                if record.protein_gpa_name:
                                    curr_gpa_name = replace_non_alphanumeric(record.protein_gpa_name)
                                    if curr_gpa_name.endswith("_"):
                                        curr_gpa_name = curr_gpa_name[:-1]
                                    protein_middle_point = (record.start + record.end) // 2
                                    protein_positions[curr_gpa_name] = protein_middle_point
    
    with open(os.path.join(temp_folder, "protein_positions.csv"), "w") as outfile:
        for protein_name, protein_position in protein_positions.items():
            outfile.write(f"{protein_name},{protein_position}\n")
                      
    combine_faa_files(temp_folder)
    delete_unwanted_faa_files(temp_folder)
                                                                         
def cdhit_runner(input_file,
    output_file,
    id=0.95,
    n_cpu=1,
    s=0.0,  # length difference cutoff (%), default 0.0
    aL=0.0,  # alignment coverage for the longer sequence
    AL=99999999,  # alignment coverage control for the longer sequence
    aS=0.0,  # alignment coverage for the shorter sequence
    AS=99999999,  # alignment coverage control for the shorter sequence
    memory=0,  # memory limit (in MB) for the program, default 8000
    accurate=True,  # use the slower but more accurate options
    use_local=False,  #whether to use local or global sequence alignment
    word_length=None,
    min_length=None,
    quiet=False):

    cdhit_script = "cd-hit"
    cdhit_script += " -T " + str(n_cpu)
    cdhit_script += " -i " + input_file
    cdhit_script += " -o " + output_file
    cdhit_script += " -c " + str(id)
    cdhit_script += " -s " + str(s)
    cdhit_script += " -aL " + str(aL)
    cdhit_script += " -AL " + str(AL)
    cdhit_script += " -aS " + str(aS)
    cdhit_script += " -AS " + str(AS)
    cdhit_script += " -M " + str(memory)
    cdhit_script += " -d 0"

    if use_local:
        cdhit_script += " -G 0"

    if accurate:
        cdhit_script += " -g 1 -n 2"

    if (word_length is not None) and (not accurate):
        cdhit_script += " -n " + str(word_length)

    if min_length is not None:
        cdhit_script += " -l " + str(min_length)

    if not quiet:
        print("running cdhit_script: " + cdhit_script)
    else:
        cdhit_script += " > /dev/null"

    subprocess.run(cdhit_script, shell=True, check=True)

    return

def save_matrix_to_csv(matrix, file_path):
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        headers = ['Protein'] + sorted(matrix[next(iter(matrix))].keys())
        writer.writerow(headers)
        for cluster_id, genomes in matrix.items():
            clean_cluster_id = str(cluster_id).strip()
            row = [clean_cluster_id] + [genomes.get(genome, 0) for genome in headers[1:]]
            writer.writerow(row)

def gene_presence_absence_file_creator(clstr_file, strains_to_be_processed, output_path):
    # Dictionary to store clusters
    clusters = defaultdict(set)
    cluster_represantative = {}

    # Parse the .clstr file
    with open(clstr_file, 'r') as f:
        cluster_id = None
        for line in f:
            if line.startswith('>Cluster'):
                cluster_id = line.strip().split()[-1]
            else:
                protein_id = line.strip().split('>')[1].split('...')[0]
                clusters[cluster_id].add(protein_id)
                if "*" in line:
                    cluster_represantative[cluster_id] = protein_id

    temp_dict = {}

    for strain in strains_to_be_processed:
        temp_dict[strain] = 0

    # Prepare the presence-absence matrix
    matrix = {}
    for cluster_id, proteins in clusters.items():
        matrix[f"{cluster_represantative[cluster_id]}"] = copy.deepcopy(temp_dict)
        for protein in proteins:
            genome = protein.split("_")[0]
            if genome in strains_to_be_processed:
                matrix[f"{cluster_represantative[cluster_id]}"][genome] = 1
            else:
                print(f"Genome {genome} is not in the list of strains to be processed")
    
    print(len(matrix))

    # Save the presence-absence matrix to a CSV file
    save_matrix_to_csv(matrix, os.path.join(output_path, 'gene_presence_absence_matrix.csv'))

    print("Gene presence-absence matrix created successfully")

