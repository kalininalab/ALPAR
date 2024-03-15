#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import pandas as pd
import os
import numpy as np
import sklearn

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


#phenotype_dataframe_creator_post_processor("/scratch/SCRATCH_SAS/alper/SR-AMR/example/example_output/binary_mutation_table.tsv", "/scratch/SCRATCH_SAS/alper/SR-AMR/example/example_output/phenotype_table.tsv")

train_strains = []
test_strains = []
    
with open(f"/home/alper/git/SR-AMR/example/example_output/datasail/splits.tsv") as splits_file:
        lines = splits_file.readlines()
        for line in lines:
            splitted = line.split("\t")
            if splitted[1].strip() == "train":
                train_strains.append(splitted[0].strip())
            elif splitted[1].strip() == "test":
                test_strains.append(splitted[0].strip())

def rf_auto_ml(phenotype_table, antibiotic, resampling_strategy="holdout", custom_scorer="MCC", fia_repeats=5, train=[], test=[], same_setup_run_count=1):

    #output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_RF_AutoML"

    #genotype_df = pd.read_csv(binary_mutation_table, sep="\t", index_col=0, header=0)
    phenotype_df = pd.read_csv(phenotype_table, sep="\t", index_col=0, header=0)

    strains_to_be_skipped_phenotype = []
    # for strain in phenotype_df.index.to_list():
    #     if strain not in genotype_df.index.to_list():
    #         strains_to_be_skipped_phenotype.append(strain)
    
    phenotype_df = phenotype_df.drop(strains_to_be_skipped_phenotype, axis=0)

    # Make sure rows are matching
    #phenotype_df = phenotype_df.reindex(genotype_df.index)

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    # Get rid of uninformative strains for given antibiotic
    strains_to_be_skipped = []

    for strain in phenotype_df.index.to_list():
        if phenotype_df.loc[strain, antibiotic] == "2" or phenotype_df.loc[strain, antibiotic] == 2:
            strains_to_be_skipped.append(strain)

    print(strains_to_be_skipped)

    print(phenotype_df.shape)
    
    #genotype_df = genotype_df.drop(strains_to_be_skipped, axis=0)
    phenotype_df = phenotype_df.drop(strains_to_be_skipped, axis=0)

    print(phenotype_df.shape)
    print(phenotype_df["ciprofloxacin"])
    print(phenotype_df["ciprofloxacin"].value_counts())

    #genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    if len(train)==0 and len(test)==0:
        #X = genotype_array[:, :].astype(int)
        y = phenotype_array[:, index_of_antibiotic].astype(int)
        #X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))
    
    else:
        #X_train = []
        y_train = []
        #X_test = []
        y_test = []

        for train_strain in train:
            if train_strain in phenotype_df.index.to_list():
                #X_train.append(genotype_df[:].loc[train_strain].astype(int))
                y_train.append(phenotype_df[antibiotic].loc[train_strain].astype(int))
        for test_strain in test:
            if test_strain in phenotype_df.index.to_list():
                #X_test.append(genotype_df[:].loc[test_strain].astype(int))
                y_test.append(phenotype_df[antibiotic].loc[test_strain].astype(int))

        #X_train = np.array(X_train)
        y_train = np.array(y_train)
        #X_test = np.array(X_test)
        y_test = np.array(y_test)
        #print("X_train")
        #print(X_train)
        print("y_train")
        print(y_train)
        print(len(y_train))
        #print("X_test")
        #print(X_test)
        print("y_test")
        print(y_test)
        print(len(y_test))

rf_auto_ml(f"/home/alper/git/SR-AMR/example/example_output/datasail/phenotype_table.tsv", "ciprofloxacin", train=train_strains, test=test_strains)
