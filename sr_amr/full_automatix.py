import os
import sys
from datetime import datetime
import random
import string

import warnings

warnings.filterwarnings("ignore")

def generate_random_key():
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))

def run_command(command, error_message):
    exit_code = os.system(command)
    if exit_code != 0:
        print(f"Error: {error_message}")
        sys.exit(1)

def automatix_runner(args):    

    if args.temp == None:
        args.temp = f"{args.output}/temp"
    
    create_binary_tables_script = f"alpar create_binary_tables -i '{args.input}' -o '{args.output}' --reference '{args.reference}' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --create_phenotype_from_folder"

    create_binary_tables_script += f" --variant_calling_tool {args.variant_calling_tool}"
    create_binary_tables_script += f" --annotation_tool {args.annotation_tool}"
    create_binary_tables_script += f" --gene_presence_absence_analysis_tool {args.gene_presence_absence_analysis_tool}"

    if args.prokka_custom_database != None:
        create_binary_tables_script += f" --prokka_custom_database '{args.prokka_custom_database[0]}' '{args.prokka_custom_database[1]}'"
    
    if args.keep_temp_files:
        create_binary_tables_script += " --keep_temp_files"

    if args.only_variants:
        create_binary_tables_script += " --only_variants"
    
    if args.checkpoint:
        create_binary_tables_script += " --checkpoint"

    if args.verbosity:
        create_binary_tables_script += f" --verbosity {args.verbosity}"

    if args.overwrite:
        create_binary_tables_script += " --overwrite"

    if args.bakta_db:
        create_binary_tables_script += f" --bakta_db '{args.bakta_db}'"

    if args.run_qc:
        create_binary_tables_script += " --run_qc"
        create_binary_tables_script += f" --qc_length_threshold {args.qc_length_threshold}"
        create_binary_tables_script += f" --qc_max_contigs {args.qc_max_contigs}"

    # Run all steps here before ml:

    print("Creating binary tables...")
    run_command(create_binary_tables_script, "Creating binary tables failed!")

    if args.only_variants:
        files_to_be_checked_list = ["binary_mutation_table.tsv", "mutations_annotations.tsv"]
    
    else:
        files_to_be_checked_list = ["binary_mutation_table.tsv", "binary_mutation_table_with_gene_presence_absence.tsv", "mutations_annotations.tsv"]

    # if not files_to_be_checked(files_to_be_checked_list, args.output):
    #     print("Error in creating binary tables!")
    #     print("Please check the logs and try again!")
    #     sys.exit(1)

    run_status_writer(f"{args.output}/status.txt", "Binary tables created")

    print("Creating binary tables done!") 

    print("Thresholding binary tables...")

    binary_table_threshold_script = ""
    if "binary_mutation_table_with_gene_presence_absence.tsv" not in os.listdir(args.output):
        if "binary_mutation_table.tsv" in os.listdir(args.output):
            binary_table_threshold_script = f"alpar binary_table_threshold -i '{args.output}/binary_mutation_table.tsv' -o '{args.output}'"
            
    else:
        binary_table_threshold_script = f"alpar binary_table_threshold -i '{args.output}/binary_mutation_table_with_gene_presence_absence.tsv' -o '{args.output}'"

    if binary_table_threshold_script == "":
        print("Error: Could not find binary mutation tables for thresholding!")
        sys.exit(1)

    if args.keep_temp_files:
                binary_table_threshold_script += " --keep_temp_files"
    
    if args.overwrite:
        binary_table_threshold_script += " --overwrite"

    run_command(binary_table_threshold_script, "Thresholding binary tables failed!")

    run_status_writer(f"{args.output}/status.txt", "Binary tables thresholded")

    print("Thresholding binary tables done!")

    if not args.fast:
        print("Running PanACoTA...")

        panacota_script = f"alpar panacota -i '{args.output}/strains.txt' -o '{args.output}' --temp '{args.temp}' --threads {args.threads} --random_names_dict '{args.output}/random_names.txt'"

        if args.keep_temp_files:
            panacota_script += " --keep_temp_files"

        if args.overwrite:
            panacota_script += " --overwrite"

        run_command(panacota_script, "PanACoTA failed!")

    files_to_be_checked_list = ["phylogenetic_tree.newick"]

    if not files_to_be_checked(files_to_be_checked_list, os.path.join(args.output, "panacota")):
        print("Phylogenetic tree has not been created by PanACoTA!")
        print("Phylogenetic tree will be created with MashTree...")
    else:
        run_status_writer(f"{args.output}/status.txt", "PanACoTA done")
        print("PanACoTA done!")

    print("Running phylogenetic tree...")

    phylogenetic_tree_script = f"alpar phylogenetic_tree -i '{args.output}/strains.txt' -o '{args.output}' --temp '{args.temp}' --random_names_dict '{args.output}/random_names.txt'"
    
    if args.keep_temp_files:
        phylogenetic_tree_script += " --keep_temp_files"

    if args.overwrite:
        phylogenetic_tree_script += " --overwrite"

    run_command(phylogenetic_tree_script, "Creating phylogenetic tree failed!")

    files_to_be_checked_list = ["phylogenetic_tree.tree"]

    if not files_to_be_checked(files_to_be_checked_list, os.path.join(args.output, "phylogeny")):
        print("Error in creating phylogenetic tree!")
        print("Please check the logs and try again!")
        sys.exit(1)

    run_status_writer(f"{args.output}/status.txt", "Phylogenetic tree done")
    print("Phylogenetic tree done!")

    print("Running PRPS...")

    binary_mutation_table_for_prps = None
    for file in os.listdir(f"{args.output}/binary_table_threshold/"):
        if file.endswith(".tsv"):
            binary_mutation_table_for_prps =  os.path.join(f"{args.output}/binary_table_threshold/", file)
            break

    if binary_mutation_table_for_prps == None:
        print("Error in finding thresholded binary mutation table for PRPS!")
        print("Falling back to the default binary mutation table...")
        if os.path.exists(os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv")):
            binary_mutation_table_for_prps = os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv")
            print("Using binary mutation table with gene presence absence...")
        else:
            binary_mutation_table_for_prps = os.path.join(args.output, "binary_mutation_table.tsv")
            print("Using binary mutation table...")

    prps_script = f"alpar prps -i '{binary_mutation_table_for_prps}' -o '{args.output}' --temp '{args.temp}' --threads {args.threads}"
    
    if args.keep_temp_files:
        prps_script += " --keep_temp_files"

    if os.path.exists(os.path.join(args.output, "panacota", "phylogenetic_tree.newick")):
        print("Using phylogenetic tree has been created by PanACoTA!")
        prps_script += f" -t '{args.output}/panacota/phylogenetic_tree.newick'"
        
    else:
        print("Phylogenetic tree has been created by the MashTree!")
        prps_script += f" -t '{args.output}/phylogeny/phylogenetic_tree.tree'"

    if args.overwrite:
        prps_script += " --overwrite"

    run_command(prps_script, "PRPS calculation failed!")

    files_to_be_checked_list = ["prps_score.tsv"]

    if not files_to_be_checked(files_to_be_checked_list, os.path.join(args.output, "prps")):
        print("Error in PRPS calculation!")
        print("Please check the logs and try again!")
        sys.exit(1)

    run_status_writer(f"{args.output}/status.txt", "PRPS done")
    print("PRPS done!")

    print("Running GWAS...")

    binary_mutation_table_for_gwas = None
    for file in os.listdir(f"{args.output}/binary_table_threshold/"):
        if file.endswith(".tsv"):
            binary_mutation_table_for_gwas =  os.path.join(f"{args.output}/binary_table_threshold/", file)
            break

    if binary_mutation_table_for_gwas == None:
        print("Error in finding thresholded binary mutation table for PRPS!")
        print("Falling back to the default binary mutation table...")
        if os.path.exists(os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv")):
            binary_mutation_table_for_gwas = os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv")
            print("Using binary mutation table with gene presence absence...")
        else:
            binary_mutation_table_for_gwas = os.path.join(args.output, "binary_mutation_table.tsv")
            print("Using binary mutation table...")

    gwas_script = f"alpar gwas -i '{binary_mutation_table_for_gwas}' -o '{args.output}' -p '{args.output}/phenotype_table.tsv' --threads {args.threads}"

    if os.path.exists(os.path.join(args.output, "panacota", "phylogenetic_tree.newick")):
        print("Using phylogenetic tree has been created by PanACoTA!")
        gwas_script += f" -t '{args.output}/panacota/phylogenetic_tree.newick'"
    
    else:
        print("Phylogenetic tree has been created by the MashTree!")
        gwas_script += f" -t '{args.output}/phylogeny/phylogenetic_tree.tree'"

    if args.overwrite:
        gwas_script += " --overwrite"

    run_command(gwas_script, "GWAS failed!")

    files_to_be_checked_list = ["genotype_matrix.tsv", "similarity_matrix.tsv"]

    if not files_to_be_checked(files_to_be_checked_list, os.path.join(args.output, "gwas")):
        print("Error in GWAS!")
        print("Please check the logs and try again!")
        sys.exit(1)

    run_status_writer(f"{args.output}/status.txt", "GWAS done")
    print("GWAS done!")

    # ML calculations starts here:

    if not args.no_ml:

        antibiotics_list = []

        with open(f"{args.output}/phenotype_table.tsv", "r") as phenotype_file:
            header_line = phenotype_file.readline()
            splitted = header_line.split("\t")
            for i in splitted[1:]:
                if not i.startswith("."):
                    antibiotics_list.append(i.strip())

        binary_mutation_table_for_ml = None
        for file in os.listdir(f"{args.output}/binary_table_threshold/"):
            if file.endswith(".tsv"):
                binary_mutation_table_for_ml =  os.path.join(f"{args.output}/binary_table_threshold/", file)
                break

        if binary_mutation_table_for_ml == None:
            print("Error in finding thresholded binary mutation table for ML!")
            print("Falling back to the default binary mutation table...")
            if os.path.exists(os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv")):
                binary_mutation_table_for_ml = os.path.join(args.output, "binary_mutation_table_with_gene_presence_absence.tsv")
                print("Using binary mutation table with gene presence absence...")
            else:
                binary_mutation_table_for_ml = os.path.join(args.output, "binary_mutation_table.tsv")
                print("Using binary mutation table...")

        for abiotic in antibiotics_list:
            for algorithm in args.ml_algorithm:
                algo_output_path = os.path.join(args.output, "ml_results", algorithm, abiotic)
                os.makedirs(algo_output_path, exist_ok=True)

                ml_script = f"alpar ml -i '{binary_mutation_table_for_ml}' -o '{algo_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} -a '{abiotic}' --save_model --prps '{args.output}/prps/prps_score.tsv' --annotation '{args.output}/mutations_annotations.tsv' --ml_algorithm '{algorithm}' --overwrite"

                if not args.no_feature_importance_analysis:
                    ml_script += " --feature_importance_analysis"

                if not args.no_datasail:
                    ml_script += f" --sail '{args.output}/strains.txt'"

                print(f"Running {algorithm} for {abiotic}...")
                run_command(ml_script, f"Running {algorithm} for {abiotic} failed!")
                run_status_writer(f"{args.output}/status.txt", f"{algorithm} done for {abiotic}")

        # After ML and GWAS are done, compare them
        if not args.no_feature_importance_analysis:
            print("Comparing GWAS and ML results...")
            for abiotic in antibiotics_list:
                pyseer_results = os.path.join(args.output, "gwas", "gwas_results", f"{abiotic}.tsv")
                
                for algorithm in args.ml_algorithm:
                    ml_importance = os.path.join(args.output, "ml_results", algorithm, abiotic, f"{abiotic}_feature_importances.tsv")
                    if os.path.exists(pyseer_results) and os.path.exists(ml_importance):
                        comparison_output = os.path.join(args.output, "ml_results", algorithm, abiotic, f"{abiotic}_gwas_ml_comparison.txt")
                        from sr_amr.gwas import compare_gwas_ml
                        compare_gwas_ml(pyseer_results, ml_importance, comparison_output)
            
    run_status_writer(f"{args.output}/status.txt", "Full automatix pipeline is done!")
    print("Full automatix pipeline is done!")


def files_to_be_checked(files_list, output_folder):
    for file in files_list:
        if not os.path.exists(os.path.join(output_folder, file)):
            return False
    return True

def run_status_writer(status_file, status):
    with open(status_file, "a") as status_file:
        status_file.write(f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}\t{status}\n")
