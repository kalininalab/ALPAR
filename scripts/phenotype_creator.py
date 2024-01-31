import os
import sys
import copy
import pandas as pd


random_names_dict = {}

with open(f"./NEW_strains_random_names_dict.tsv", "r") as infile:
    lines = infile.readlines()
    for line in lines:
        random_names_dict[line.split("\t")[0]] = line.split("\t")[1].strip()

def phenotype_dataframe_creator(data_folder_path, output_file, random_names_dict=random_names_dict):
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

    print(df)

    df.to_csv(output_file, sep="\t")

