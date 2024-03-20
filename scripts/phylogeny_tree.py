#!/usr/bin/env python3

import os
import shutil

def mash_preprocessor(strains_text_file, output_folder, temp_folder, random_names_dict):

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

    if os.path.exists(f"{os.path.dirname(output_folder)}/snippy"):

        snippy_dir = os.listdir(f"{os.path.dirname(output_folder)}/snippy")

        the_names_will_be_skipped = []

        copied_files = os.listdir(f"{fasta_files_folder}")

        for file in copied_files:
            if file not in snippy_dir:
                the_names_will_be_skipped.append(file)

        for file in the_names_will_be_skipped:
            os.remove(f"{fasta_files_folder}/{file}")

def mash_distance_runner(output_folder, temp_folder):

    # Get the absolute path of the current file
    temp_folder_abs_path = os.path.abspath(temp_folder)

    print(temp_folder_abs_path)

    mash_tree_command = f"mashtree {temp_folder_abs_path}/fasta_files/* > {output_folder}/phylogenetic_tree.tree"
    
    os.system(mash_tree_command)
