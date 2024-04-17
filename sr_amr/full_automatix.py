import os

def automatix_runner(args):    

    # Get the absolute path of the current script
    current_script_path = os.path.abspath(__file__)

    # Get the directory containing the current script
    current_directory = os.path.dirname(current_script_path)

    if args.temp == None:
        args.temp = f"{args.output}/temp"

    # Create the temporary directory
    os.makedirs(args.temp, exist_ok=True)
    
    create_binary_tables_script = f"sr-amr create_binary_tables -i '{args.input}' -o '{args.output}' --reference '{args.reference}' --temp '{args.temp}' --threads {args.threads} --ram {args.ram}"

    if args.custom_database != None:
        create_binary_tables_script += f" --custom_database '{args.custom_database}'"
    
    if args.keep_temp_files:
        create_binary_tables_script += " --keep_temp_files"

    if args.just_mutations:
        binary_table_threshold_script = f"sr-amr binary_table_threshold -i '{args.output}/binary_mutation_table.tsv' -o '{args.output}'"
        if args.keep_temp_files:
            binary_table_threshold_script += " --keep_temp_files"

    else:
        binary_table_threshold_script = f"sr-amr binary_table_threshold -i '{args.output}/binary_mutation_table_with_gene_presence_absence.tsv' -o '{args.output}'"
        if args.keep_temp_files:
            binary_table_threshold_script += " --keep_temp_files"

    panacota_script = f"sr-amr panacota -i '{args.output}/strains.txt' -o '{args.output}' --temp '{args.temp}' --threads {args.threads} --random_names_dict '{args.output}/random_names.txt'"

    if args.keep_temp_files:
        panacota_script += " --keep_temp_files"

    phylogenetic_tree_script = f"sr-amr phylogenetic_tree -i '{args.output}/strains.txt' -o '{args.output}' --temp '{args.temp}' --random_names_dict '{args.output}/random_names.txt'"
    
    if args.keep_temp_files:
        phylogenetic_tree_script += " --keep_temp_files"

    if args.just_mutations:
        prps_script = f"sr-amr prps -i '{args.output}/binary_table_threshold/binary_mutation_table_threshold_0.2_percent.tsv' -o '{args.output}' --temp '{args.temp}' --threads {args.threads}"

        if args.keep_temp_files:
            prps_script += " --keep_temp_files"

    else:
        prps_script = f"sr-amr prps -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{args.output}' --temp '{args.temp}' --threads {args.threads}"
        
        if args.keep_temp_files:
            prps_script += " --keep_temp_files"

    if args.just_mutations:
        gwas_script = f"sr-amr gwas -i '{args.output}/binary_table_threshold/binary_mutation_table_threshold_0.2_percent.tsv' -o '{args.output}' -p '{args.output}/phenotype_table.tsv'"
        if args.keep_temp_files:
            gwas_script += " --keep_temp_files"
    else:
        gwas_script = f"sr-amr gwas -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{args.output}' -p '{args.output}/phenotype_table.tsv'"
        if args.keep_temp_files:
            gwas_script += " --keep_temp_files"

    # Run all steps here before ml:

    print("Creating binary tables...")
    os.system(create_binary_tables_script)
    print("Creating binary tables done!")

    print("Thresholding binary tables...")
    os.system(binary_table_threshold_script)
    print("Thresholding binary tables done!")

    print("Running PanACoTA...")
    os.system(panacota_script)
    print("PanACoTA done!")

    print("Running phylogenetic tree...")
    os.system(phylogenetic_tree_script)
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
    print("PRPS done!")

    print("Running GWAS...")
    os.system(gwas_script)
    print("GWAS done!")

    # ML calculations starts here:

    rf_output_path = os.path.join(args.output, "rf_output")
    svm_output_path = os.path.join(args.output, "svm_output")
    gb_output_path = os.path.join(args.output, "gb_output")

    os.makedirs(rf_output_path, exist_ok=True)
    os.makedirs(svm_output_path, exist_ok=True)
    os.makedirs(gb_output_path, exist_ok=True)

    antibiotics_list = []

    with open(f"{args.output}/phenotype_table.tsv", "r") as phenotype_file:
        header_line = phenotype_file.readline()
        splitted = header_line.split("\t")
        for i in splitted[1:]:
            if not i.startswith("."):
                antibiotics_list.append(i)

    for abiotic in antibiotics_list:
        abiotic_rf_output_path = os.path.join(rf_output_path, f"{abiotic}")
        abiotic_svm_output_path = os.path.join(svm_output_path, f"{abiotic}")
        abiotic_gb_output_path = os.path.join(gb_output_path, f"{abiotic}")

        os.makedirs(abiotic_rf_output_path, exist_ok=True)
        os.makedirs(abiotic_svm_output_path, exist_ok=True)
        os.makedirs(abiotic_gb_output_path, exist_ok=True)

        ml_script_rf = f"sr-amr ml -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{abiotic_rf_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --sail {args.output}/strains.txt -a '{abiotic}' --save_model --feature_importance_analysis --prps '{args.output}/prps/prps_score.tsv --parameter_optimization --ml_algorithm 'rf'"

        ml_script_svm = f"sr-amr ml -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{abiotic_rf_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --sail {args.output}/strains.txt -a '{abiotic}' --save_model --feature_importance_analysis --prps '{args.output}/prps/prps_score.tsv --parameter_optimization --ml_algorithm 'svm'"

        ml_script_gb = f"sr-amr ml -i '{args.output}/binary_table_threshold/binary_mutation_table_with_gene_presence_absence_threshold_0.2_percent.tsv' -o '{abiotic_rf_output_path}' -p '{args.output}/phenotype_table.tsv' --temp '{args.temp}' --threads {args.threads} --ram {args.ram} --sail {args.output}/strains.txt -a '{abiotic}' --save_model --feature_importance_analysis --prps '{args.output}/prps/prps_score.tsv --parameter_optimization --ml_algorithm 'gb'"

        print(f"Running Random Forest for {abiotic}...")
        os.system(ml_script_rf)

        print(f"Running Support Vector Machine for {abiotic}...")
        os.system(ml_script_svm)

        print(f"Running Gradient Boosting for {abiotic}...")
        os.system(ml_script_gb)

    print("Full automatix pipeline is done!")