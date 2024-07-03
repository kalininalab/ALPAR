import os
import sys
from datetime import datetime
import random
import string

def generate_random_key():
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))

def automatix_runner(args):    

    if args.temp == None:
        args.temp = f"{args.output}/temp"
    
    create_binary_tables_script = f"alpar create_binary_tables -i '{args.input}' -o '{args.output}' --reference '{args.reference}' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --create_phenotype_from_folder {args.input}"

    if args.custom_database != None:
        custom_db_name = generate_random_key()
        create_binary_tables_script += f" --custom_database '{args.custom_database}' {custom_db_name}"
    
    if args.keep_temp_files:
        create_binary_tables_script += " --keep_temp_files"

    if args.just_mutations:
        create_binary_tables_script += " --no_gene_presence_absence"
    
    if args.checkpoint:
        create_binary_tables_script += " --checkpoint"

    if args.verbosity:
        create_binary_tables_script += f" --verbosity {args.verbosity}"

    if args.just_mutations:
        binary_table_threshold_script = f"alpar binary_table_threshold -i '{args.output}/binary_mutation_table.tsv' -o '{args.output}'"
        if args.keep_temp_files:
            binary_table_threshold_script += " --keep_temp_files"

    else:
        binary_table_threshold_script = f"alpar binary_table_threshold -i '{args.output}/binary_mutation_table_with_gene_presence_absence.tsv' -o '{args.output}'"
        if args.keep_temp_files:
            binary_table_threshold_script += " --keep_temp_files"

    panacota_script = f"alpar panacota -i '{args.output}/strains.txt' -o '{args.output}' --temp '{args.temp}' --threads {args.threads} --random_names_dict '{args.output}/random_names.txt'"

    if args.keep_temp_files:
        panacota_script += " --keep_temp_files"

    phylogenetic_tree_script = f"alpar phylogenetic_tree -i '{args.output}/strains.txt' -o '{args.output}' --temp '{args.temp}' --random_names_dict '{args.output}/random_names.txt'"
    
    if args.keep_temp_files:
        phylogenetic_tree_script += " --keep_temp_files"

    if args.just_mutations:
        prps_script = f"alpar prps -i '{args.output}/binary_table_threshold/binary_mutation_table_threshold_0.2_percent.tsv' -o '{args.output}' --temp '{args.temp}' --threads {args.threads}"

        if args.keep_temp_files:
            prps_script += " --keep_temp_files"

    else:
        prps_script = f"alpar prps -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{args.output}' --temp '{args.temp}' --threads {args.threads}"
        
        if args.keep_temp_files:
            prps_script += " --keep_temp_files"

    if args.just_mutations:
        gwas_script = f"alpar gwas -i '{args.output}/binary_table_threshold/binary_mutation_table_threshold_0.2_percent.tsv' -o '{args.output}' -p '{args.output}/phenotype_table.tsv' --threads {args.threads}"
    else:
        gwas_script = f"alpar gwas -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{args.output}' -p '{args.output}/phenotype_table.tsv' --threads {args.threads}"

    if args.overwrite:
        create_binary_tables_script += " --overwrite"
        binary_table_threshold_script += " --overwrite"
        panacota_script += " --overwrite"
        phylogenetic_tree_script += " --overwrite"
        prps_script += " --overwrite"
        gwas_script += " --overwrite"

    # Run all steps here before ml:

    print("Creating binary tables...")
    os.system(create_binary_tables_script)

    if args.just_mutations:
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

    if "binary_mutation_table_with_gene_presence_absence.tsv" not in os.listdir(args.output):
        if "binary_mutation_table.tsv" in os.listdir(args.output):
            binary_table_threshold_script = f"alpar binary_table_threshold -i '{args.output}/binary_mutation_table.tsv' -o '{args.output}'"
            if args.keep_temp_files:
                binary_table_threshold_script += " --keep_temp_files"

    os.system(binary_table_threshold_script)
    

    if args.just_mutations:
        files_to_be_checked_list = ["binary_mutation_table_threshold_0.2_percent.tsv"]
    
    else:
        files_to_be_checked_list = ["binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv"]

    if not files_to_be_checked(files_to_be_checked_list, os.path.join(args.output, "binary_table_threshold")):
        print("Error in thresholding binary tables!")
        print("Please check the logs and try again!")
        sys.exit(1)

    run_status_writer(f"{args.output}/status.txt", "Binary tables thresholded")

    print("Thresholding binary tables done!")

    if not args.fast:
        print("Running PanACoTA...")
        os.system(panacota_script)

    files_to_be_checked_list = ["phylogenetic_tree.newick"]

    if not files_to_be_checked(files_to_be_checked_list, os.path.join(args.output, "panacota")):
        print("Phylogenetic tree has not been created by PanACoTA!")
        print("Phylogenetic tree will be created with MashTree...")
    else:
        run_status_writer(f"{args.output}/status.txt", "PanACoTA done")
        print("PanACoTA done!")

    print("Running phylogenetic tree...")

    os.system(phylogenetic_tree_script)

    files_to_be_checked_list = ["phylogenetic_tree.tree"]

    if not files_to_be_checked(files_to_be_checked_list, os.path.join(args.output, "phylogeny")):
        print("Error in creating phylogenetic tree!")
        print("Please check the logs and try again!")
        sys.exit(1)

    run_status_writer(f"{args.output}/status.txt", "Phylogenetic tree done")
    print("Phylogenetic tree done!")

    if os.path.exists(os.path.join(args.output, "panacota", "phylogenetic_tree.newick")):
        print("Using phylogenetic tree has been created by PanACoTA!")
        prps_script += f" -t '{args.output}/panacota/phylogenetic_tree.newick'"
        gwas_script += f" -t '{args.output}/panacota/phylogenetic_tree.newick'"
    
    else:
        print("Phylogenetic tree has been created by the MashTree!")
        prps_script += f" -t '{args.output}/phylogeny/phylogenetic_tree.tree'"
        gwas_script += f" -t '{args.output}/phylogeny/phylogenetic_tree.tree'"

    print("Running PRPS...")
    os.system(prps_script)

    files_to_be_checked_list = ["prps_score.tsv"]

    if not files_to_be_checked(files_to_be_checked_list, os.path.join(args.output, "prps")):
        print("Error in PRPS calculation!")
        print("Please check the logs and try again!")
        sys.exit(1)

    run_status_writer(f"{args.output}/status.txt", "PRPS done")
    print("PRPS done!")

    print("Running GWAS...")
    os.system(gwas_script)

    files_to_be_checked_list = ["genotype_matrix.tsv", "similarity_matrix.tsv"]

    if not files_to_be_checked(files_to_be_checked_list, os.path.join(args.output, "gwas")):
        print("Error in GWAS!")
        print("Please check the logs and try again!")
        sys.exit(1)

    run_status_writer(f"{args.output}/status.txt", "GWAS done")
    print("GWAS done!")

    # ML calculations starts here:

    if not args.no_ml:

        if "rf" in args.ml_algorithm:
            rf_output_path = os.path.join(args.output, "rf_output")
            os.makedirs(rf_output_path, exist_ok=True)
        
        if "svm" in args.ml_algorithm:
            svm_output_path = os.path.join(args.output, "svm_output")
            os.makedirs(svm_output_path, exist_ok=True)
        
        if "gb" in args.ml_algorithm:
            gb_output_path = os.path.join(args.output, "gb_output")
            os.makedirs(gb_output_path, exist_ok=True)

        antibiotics_list = []

        with open(f"{args.output}/phenotype_table.tsv", "r") as phenotype_file:
            header_line = phenotype_file.readline()
            splitted = header_line.split("\t")
            for i in splitted[1:]:
                if not i.startswith("."):
                    antibiotics_list.append(i.strip())

        for abiotic in antibiotics_list:
            
            if "rf" in args.ml_algorithm:
                abiotic_rf_output_path = os.path.join(rf_output_path, f"{abiotic}")
                os.makedirs(abiotic_rf_output_path, exist_ok=True)

            if "svm" in args.ml_algorithm:
                abiotic_svm_output_path = os.path.join(svm_output_path, f"{abiotic}")
                os.makedirs(abiotic_svm_output_path, exist_ok=True)

            if "gb" in args.ml_algorithm:
                abiotic_gb_output_path = os.path.join(gb_output_path, f"{abiotic}")
                os.makedirs(abiotic_gb_output_path, exist_ok=True)

            if args.no_feature_importance_analysis:
                
                if "rf" in args.ml_algorithm:
                    ml_script_rf = f"alpar ml -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{abiotic_rf_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --sail {args.output}/strains.txt -a '{abiotic}' --save_model --prps '{args.output}/prps/prps_score.tsv' --parameter_optimization --annotation '{args.output}/mutations_annotations.tsv' --ml_algorithm 'rf' --overwrite"

                if "svm" in args.ml_algorithm:
                    ml_script_svm = f"alpar ml -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{abiotic_svm_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --sail {args.output}/strains.txt -a '{abiotic}' --save_model --prps '{args.output}/prps/prps_score.tsv' --parameter_optimization --annotation '{args.output}/mutations_annotations.tsv' --ml_algorithm 'svm' --overwrite"

                if "gb" in args.ml_algorithm:
                    ml_script_gb = f"alpar ml -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{abiotic_gb_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --sail {args.output}/strains.txt -a '{abiotic}' --save_model --prps '{args.output}/prps/prps_score.tsv' --parameter_optimization --annotation '{args.output}/mutations_annotations.tsv' --ml_algorithm 'gb' --overwrite"

            else:
                
                if "rf" in args.ml_algorithm:
                    ml_script_rf = f"alpar ml -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{abiotic_rf_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --sail {args.output}/strains.txt -a '{abiotic}' --save_model --feature_importance_analysis --prps '{args.output}/prps/prps_score.tsv' --parameter_optimization --annotation '{args.output}/mutations_annotations.tsv' --ml_algorithm 'rf' --overwrite"

                if "svm" in args.ml_algorithm:
                    ml_script_svm = f"alpar ml -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{abiotic_svm_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --sail {args.output}/strains.txt -a '{abiotic}' --save_model --feature_importance_analysis --prps '{args.output}/prps/prps_score.tsv' --parameter_optimization --annotation '{args.output}/mutations_annotations.tsv' --ml_algorithm 'svm' --overwrite"

                if "gb" in args.ml_algorithm:
                    ml_script_gb = f"alpar ml -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{abiotic_gb_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --sail {args.output}/strains.txt -a '{abiotic}' --save_model --feature_importance_analysis --prps '{args.output}/prps/prps_score.tsv' --parameter_optimization --annotation '{args.output}/mutations_annotations.tsv' --ml_algorithm 'gb' --overwrite"

            if "rf" in args.ml_algorithm:
                print(f"Running Random Forest for {abiotic}...")
                os.system(ml_script_rf)
                run_status_writer(f"{args.output}/status.txt", f"RF done for {abiotic}")

            if "svm" in args.ml_algorithm:
                print(f"Running Support Vector Machine for {abiotic}...")
                os.system(ml_script_svm)
                run_status_writer(f"{args.output}/status.txt", f"SVM done for {abiotic}")

            if "gb" in args.ml_algorithm:
                print(f"Running Gradient Boosting for {abiotic}...")
                os.system(ml_script_gb)
                run_status_writer(f"{args.output}/status.txt", f"GB done for {abiotic}")
            
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