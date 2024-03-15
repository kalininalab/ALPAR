import pandas as pd
import os


def phenotype_dataframe_creator_post_processor(genotype_data, phenotype_data):
    
    genotype_df = pd.read_csv(genotype_data, sep="\t", index_col=0, header=0)

    phenotype_df = pd.read_csv(phenotype_data, sep="\t", index_col=0, header=0)

    strains_to_be_dropped = []

    phenotype_strains = phenotype_df.index.to_list()

    genotype_strains = genotype_df.index.to_list()

    for strain in phenotype_strains:
        if not strain in genotype_strains:
            strains_to_be_dropped.append(strain)

    phenotype_df_dropped = phenotype_df.drop(strains_to_be_dropped, axis=0)

    print(strains_to_be_dropped)

    phenotype_df_dropped.to_csv(phenotype_data, sep="\t")


phenotype_dataframe_creator_post_processor("/scratch/SCRATCH_SAS/alper/SR-AMR/example/example_output/binary_mutation_table.tsv", "/scratch/SCRATCH_SAS/alper/SR-AMR/example/example_output/phenotype_table.tsv")