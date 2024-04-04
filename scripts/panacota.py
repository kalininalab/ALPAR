#!/usr/bin/env python3

import os
import shutil
import pathlib

# Get the path of the script
PATH_OF_SCRIPT = pathlib.Path(__file__).parent.resolve()


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
                random_names_to_strains[splitted[1].strip(
                )] = splitted[0].strip()

    already_copied = []

    with open(list_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            strain_path = line.strip()
            if random_names_will_be_used:
                if not random_names[os.path.splitext(line.split('/')[-1].strip())[0]] in already_copied:
                    shutil.copy2(
                        strain_path, f"{temp_folder}/{random_names[os.path.splitext(line.split('/')[-1].strip())[0]]}")
                    already_copied.append(
                        random_names[os.path.splitext(line.split('/')[-1].strip())[0]])
            else:
                if not line.split('/')[-1].strip() in already_copied:
                    shutil.copy2(
                        strain_path, f"{temp_folder}/{line.split('/')[-1].strip()}")
                    already_copied.append(line.split('/')[-1].strip())

    the_names_will_be_skipped = []

    if os.path.exists(f"{os.path.dirname(output_folder)}/snippy_processed_strains.txt"):

        snippy_dir = []

        with open(f"{os.path.dirname(output_folder)}/snippy_processed_strains.txt", "r") as snippy_processed_strains_file:
            snippy_processed_strains_file_lines = snippy_processed_strains_file.readlines()
            for line in snippy_processed_strains_file_lines:
                snippy_dir.append(line.strip())

        for strain in os.listdir(temp_folder):
            if strain not in snippy_dir:
                the_names_will_be_skipped.append(strain)

    if len(the_names_will_be_skipped) > 0:
        print(
            f"The following strains will be skipped for panacota: {the_names_will_be_skipped}")

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

    shutil.copyfile(f"{panacota_output_folder}/phylogenetic_tree.newick",
                    f"{output_folder}/phylogenetic_tree.newick")


def panacota_pipeline_runner(list_file, dbpath, output_directory, run_name, n_cores, log_file, type="nucl", mode=1, min_seq_id=0.8, core_genome_percentage=1):

    pc_annotate_command = f"PanACoTA annotate -l {list_file} -d {dbpath} -r {output_directory}/annotate_out -n {run_name} --threads {n_cores} >> {log_file} 2>&1"

    pc_pangenome_command = f"PanACoTA pangenome -l {output_directory}/annotate_out/LSTINFO-{list_file.split('/')[-1]} -d {output_directory}/annotate_out/Proteins/ -o {output_directory}/pangenome_out/ -n {run_name} --threads {n_cores} >> {log_file} 2>&1"

    if n_cores > 1:
        pc_corepers_command = f"PanACoTA corepers -p {output_directory}/pangenome_out/PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}-th{n_cores}.lst -o {output_directory}/corepers_out/ -t {core_genome_percentage} >> {log_file} 2>&1"
        pc_align_command = f"PanACoTA align -c {output_directory}/corepers_out/PersGenome_PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}-th{n_cores}.lst-all_1.lst -l {output_directory}/annotate_out/LSTINFO-{list_file.split('/')[-1]} -n {run_name} -d {output_directory}/annotate_out/ -o {output_directory}/align_out --threads {n_cores} >> {log_file} 2>&1"

    else:
        pc_corepers_command = f"PanACoTA corepers -p {output_directory}/pangenome_out/PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}.lst -o {output_directory}/corepers_out/ -t {core_genome_percentage} >> {log_file} 2>&1"
        pc_align_command = f"PanACoTA align -c {output_directory}/corepers_out/PersGenome_PanGenome-{run_name}.All.prt-clust-{min_seq_id}-mode{mode}.lst-all_1.lst -l {output_directory}/annotate_out/LSTINFO-{list_file.split('/')[-1]} -n {run_name} -d {output_directory}/annotate_out/ -o {output_directory}/align_out --threads {n_cores} >> {log_file} 2>&1"

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
