from datasail.routine import datasail_main
import pandas as pd
from pathlib import Path
import os
import shutil


def datasail_pre_precessor(strains_text_file, temp_folder, random_names_dict, output_folder, cpus=1):

    os.mkdir(os.path.join(temp_folder, "fasta_files"))

    fasta_files_folder = os.path.join(temp_folder, "fasta_files")

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

    with open(strains_text_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            strain_path = line.strip()
            if random_names_will_be_used:
                if not random_names[os.path.splitext(line.split('/')[-1].strip())[0]] in already_copied:
                    shutil.copy2(
                        strain_path, f"{fasta_files_folder}/{random_names[os.path.splitext(line.split('/')[-1].strip())[0]]}")
                    already_copied.append(
                        random_names[os.path.splitext(line.split('/')[-1].strip())[0]])
            else:
                if not line.split('/')[-1].strip() in already_copied:
                    shutil.copy2(
                        strain_path, f"{fasta_files_folder}/{line.split('/')[-1].strip()}")
                    already_copied.append(line.split('/')[-1].strip())

    if os.path.exists(f"{os.path.dirname(output_folder)}/snippy_processed_strains.txt"):

        snippy_dir = []

        with open(f"{os.path.dirname(output_folder)}/snippy_processed_strains.txt", "r") as snippy_processed_strains_file:
            snippy_processed_strains_file_lines = snippy_processed_strains_file.readlines()
            for line in snippy_processed_strains_file_lines:
                snippy_dir.append(line.strip())

        the_names_will_be_skipped = []

        copied_files = os.listdir(f"{fasta_files_folder}")

        for file in copied_files:
            if file not in snippy_dir:
                the_names_will_be_skipped.append(file)

        for file in the_names_will_be_skipped:
            os.remove(f"{fasta_files_folder}/{file}")

    mash_sketch_command = f"mash sketch -p {cpus} -o {temp_folder}/mash_sketch {fasta_files_folder}/*"

    os.system(mash_sketch_command)

    mash_dist_command = f"mash dist -p {cpus} -t {temp_folder}/mash_sketch.msh {temp_folder}/mash_sketch.msh > {temp_folder}/distance_matrix.tsv"

    os.system(mash_dist_command)

    with open(f"{output_folder}/distance_matrix.tsv", "w") as ofile:
        with open(f"{temp_folder}/distance_matrix.tsv") as dist_matr_file:
            lines = dist_matr_file.readlines()

        first_line_splitted = lines[0].split("\t")

        first_line = ""
        for split in first_line_splitted:
            if "/" in split:
                splitted2 = split.split("/")
                first_line += f"{splitted2[-1].strip()}\t"
            else:
                first_line += f"{split.strip()}\t"

        ofile.write(f"{first_line.strip()}\n")

        for line in lines[1:]:

            print_line = ""

            splitted = line.split("\t")

            for split in splitted:
                if "/" in split:
                    splitted2 = split.split("/")
                    print_line += f"{splitted2[-1].strip()}\t"
                else:
                    print_line += f"{split.strip()}\t"

            ofile.write(f"{print_line.strip()}\n")


def datasail_runner(distance_matrix, output_folder, splits=[0.8, 0.2], cpus=1, max_time=600):
    """
    Runs the datasail algorithm on the distance matrix.

    Parameters
    ----------
    distance_matrix : pandas.DataFrame
        Distance matrix between all pairs of nodes in the tree.
    output_folder : str
        Path to the output folder.
    temp_folder : str
        Path to the temporary folder.

    Returns
    -------
    datasail_output 
    """

    dm = pd.read_csv(distance_matrix, sep="\t", index_col=0, header=0)

    splits, _, _ = datasail_main(output=None, techniques=["C1e"], splits=splits, names=["train", "test"], e_type="P", e_data=[(n, "a" * i) for i, n in enumerate(dm.columns)], e_dist=Path(distance_matrix), max_sec=max_time, threads=cpus, inter=None, max_sol=10, verbosity="I", delta=0.1,
                                 epsilon=0.1, runs=1, solver="SCIP", cache=False, cache_dir=None, linkage="average", e_weights=None, e_strat=None, e_sim=None, e_args="", e_clusters=50, f_type=None, f_data=None, f_weights=None, f_strat=None, f_sim=None, f_dist=None, f_args="", f_clusters=50, cli=False, logdir=None)

    with open(f"{output_folder}/splits.tsv", "w") as ofile:
        for key in splits["C1e"][0]:
            ofile.write(f"{key}\t{splits['C1e'][0][key]}\n")

    return None
