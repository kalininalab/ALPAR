import sklearn.model_selection
import sklearn.datasets
import sklearn.metrics
from sklearn.inspection import permutation_importance
import pandas as pd
import numpy as np
import pickle
import os
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier


def output_file_writer(outfile, y_test, y_hat, cls = None, best_c = None):

    with open(outfile, "w") as ofile:

        # New version of sklearn is not working with auto sk-learn ensemble models
        
        # if cls:
        #     ofile.write(str(cls.show_models()))
        #     ofile.write("\n")
        if best_c:
            ofile.write("C: " + str(best_c))
            ofile.write("\n")
        ofile.write("Accuracy score: " + str(sklearn.metrics.accuracy_score(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Balanced Accuracy score: " + str(sklearn.metrics.balanced_accuracy_score(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Brier score loss: " + str(sklearn.metrics.brier_score_loss(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("F1 score macro: " + str(sklearn.metrics.f1_score(y_test, y_hat, average='macro')))
        ofile.write("\n")
        ofile.write("F1 score micro: " + str(sklearn.metrics.f1_score(y_test, y_hat, average='micro')))
        ofile.write("\n")
        ofile.write("F1 score weighted: " + str(sklearn.metrics.f1_score(y_test, y_hat, average='weighted')))
        ofile.write("\n")
        ofile.write("F1 score binary: " + str(sklearn.metrics.f1_score(y_test, y_hat, average='binary')))
        ofile.write("\n")
        ofile.write("Precision score: " + str(sklearn.metrics.precision_score(y_test, y_hat, average='binary')))
        ofile.write("\n")
        ofile.write("Recall score: " + str(sklearn.metrics.recall_score(y_test, y_hat, average='binary')))
        ofile.write("\n")
        ofile.write("Confussion matrix: " + str(sklearn.metrics.confusion_matrix(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("ROC Curve: " + str(sklearn.metrics.roc_curve(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("ROC AUC Score: " + str(sklearn.metrics.roc_auc_score(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Jaccard score: " + str(sklearn.metrics.jaccard_score(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Hinge loss: " + str(sklearn.metrics.hinge_loss(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Hamming loss: " + str(sklearn.metrics.hamming_loss(y_test, y_hat)))
        ofile.write("\n")
        ofile.write(
            "Fbeta score macro: " + str(sklearn.metrics.fbeta_score(y_test, y_hat, average='macro', beta=0.5)))
        ofile.write("\n")
        ofile.write(
            "Fbeta score micro: " + str(sklearn.metrics.fbeta_score(y_test, y_hat, average='micro', beta=0.5)))
        ofile.write("\n")
        ofile.write("Fbeta score weighted: " + str(
            sklearn.metrics.fbeta_score(y_test, y_hat, average='weighted', beta=0.5)))
        ofile.write("\n")
        ofile.write("Log loss: " + str(sklearn.metrics.log_loss(y_test, y_hat)))
        ofile.write("\n")
        ofile.write("Matthews correlation coefficient: " + str(sklearn.metrics.matthews_corrcoef(y_test, y_hat)))


def rf(binary_mutation_table, phenotype_table, antibiotic, random_seed, cv_split, test_size, output_folder, n_jobs, feature_importance_analysis = False, save_model = False, resampling_strategy="holdout", fia_repeats=5, custom_scorer="MCC", n_estimators=100, max_depth=2, min_samples_leaf=1, min_samples_split=2, train=[], test=[], same_setup_run_count=1):
    
    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_RF"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t", index_col=0, header=0)
    phenotype_df = pd.read_csv(phenotype_table, sep="\t", index_col=0, header=0)
    
    strains_to_be_skipped_phenotype = []
    for strain in phenotype_df.index.to_list():
        if strain not in genotype_df.index.to_list():
            strains_to_be_skipped_phenotype.append(strain)
    
    phenotype_df = phenotype_df.drop(strains_to_be_skipped_phenotype, axis=0)

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    # Get rid of uninformative strains for given antibiotic
    strains_to_be_skipped = []

    for strain in phenotype_df.index.to_list():
        if phenotype_df.loc[strain, antibiotic] == "2" or phenotype_df.loc[strain, antibiotic] == 2:
            strains_to_be_skipped.append(strain)
    
    genotype_df = genotype_df.drop(strains_to_be_skipped, axis=0)
    phenotype_df = phenotype_df.drop(strains_to_be_skipped, axis=0)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    if len(train)==0 and len(test)==0:
        X = genotype_array[:, :].astype(int)
        y = phenotype_array[:, index_of_antibiotic].astype(int)
        X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))
    
    else:
        X_train = []
        y_train = []
        X_test = []
        y_test = []

        for train_strain in train:
            if train_strain in genotype_df.index.to_list():
                X_train.append(genotype_df[:].loc[train_strain].astype(int))
                y_train.append(phenotype_df[antibiotic].loc[train_strain].astype(int))
        for test_strain in test:
            if test_strain in genotype_df.index.to_list():
                X_test.append(genotype_df[:].loc[test_strain].astype(int))
                y_test.append(phenotype_df[antibiotic].loc[test_strain].astype(int))

        X_train = np.array(X_train)
        y_train = np.array(y_train)
        X_test = np.array(X_test)
        y_test = np.array(y_test)

    rf_cls = RandomForestClassifier(class_weight={0: sum(y_train), 1: len(y_train) - sum(y_train)}, n_estimators=n_estimators, max_depth=max_depth, min_samples_leaf=min_samples_leaf, min_samples_split=min_samples_split)

    param_grid = {
    'n_estimators': [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100],
    'max_depth': [2, 5, 7, 9]
    }

    if resampling_strategy == "cv":

        if custom_scorer == "MCC":
            scorer = "matthews_corrcoef"

        # Create a GridSearchCV instance to find the best hyperparameters
            
        grid_search = GridSearchCV(rf_cls, param_grid, cv=cv_split, scoring=scorer)
        grid_search.fit(X_train, y_train)

        y_hat = grid_search.predict(X_test)

    else:

        rf_cls.fit(X_train, y_train)
        y_hat = rf_cls.predict(X_test)

    outfile = os.path.join(output_folder, f"{output_file_template}_Result") 

    output_file_writer(outfile, y_test, y_hat)

    if feature_importance_analysis:

        r = permutation_importance(rf_cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA"), "w") as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:

        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(rf_cls, open(model_file, 'wb'))


def svm(binary_mutation_table, phenotype_table, antibiotic, random_seed, test_size, output_folder, n_jobs, feature_importance_analysis = False, save_model = False, resampling_strategy="holdout", fia_repeats=5, optimization=False, kernel="linear", train=[], test=[]):

    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_SVM"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t", index_col=0, header=0)
    phenotype_df = pd.read_csv(phenotype_table, sep="\t", index_col=0, header=0)

    strains_to_be_skipped_phenotype = []
    for strain in phenotype_df.index.to_list():
        if strain not in genotype_df.index.to_list():
            strains_to_be_skipped_phenotype.append(strain)
    
    phenotype_df = phenotype_df.drop(strains_to_be_skipped_phenotype, axis=0)

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    # Get rid of uninformative strains for given antibiotic
    strains_to_be_skipped = []

    for strain in phenotype_df.index.to_list():
        if phenotype_df.loc[strain, antibiotic] == "2" or phenotype_df.loc[strain, antibiotic] == 2:
            strains_to_be_skipped.append(strain)
    
    genotype_df = genotype_df.drop(strains_to_be_skipped, axis=0)
    phenotype_df = phenotype_df.drop(strains_to_be_skipped, axis=0)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    if len(train)==0 and len(test)==0:
        X = genotype_array[:, :].astype(int)
        y = phenotype_array[:, index_of_antibiotic].astype(int)
        X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))
    
    else:
        X_train = []
        y_train = []
        X_test = []
        y_test = []

        for train_strain in train:
            if train_strain in genotype_df.index.to_list():
                X_train.append(genotype_df[:].loc[train_strain].astype(int))
                y_train.append(phenotype_df[antibiotic].loc[train_strain].astype(int))
        for test_strain in test:
            if test_strain in genotype_df.index.to_list():
                X_test.append(genotype_df[:].loc[test_strain].astype(int))
                y_test.append(phenotype_df[antibiotic].loc[test_strain].astype(int))

        X_train = np.array(X_train)
        y_train = np.array(y_train)
        X_test = np.array(X_test)
        y_test = np.array(y_test)

    best_model_mcc = -1.0
    bm_c = 0

    max_c_range = 2

    if optimization:
        max_c_range = 11

    for c_val in np.arange(1, max_c_range, 1):
        
        svm_cls = SVC(class_weight={0: sum(y_train), 1: len(y_train) - sum(y_train)}, kernel=kernel, C=c_val)
        svm_cls.fit(X_train, y_train)

        y_hat = svm_cls.predict(X_test)

        cur_mcc_val = sklearn.metrics.matthews_corrcoef(y_test, y_hat)
        if cur_mcc_val > best_model_mcc:
            best_model_mcc = cur_mcc_val
            best_model = svm_cls
            bm_c = c_val
    
    outfile = os.path.join(output_folder, f"{output_file_template}_Result") 

    output_file_writer(outfile, y_test, y_hat, best_c=bm_c)
    
    if feature_importance_analysis:

        r = permutation_importance(best_model, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA"), "w") as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:

        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(best_model, open(model_file, 'wb'))

    
def svm_cv(binary_mutation_table, phenotype_table, antibiotic, random_seed, test_size, output_folder, n_jobs, cv_split, feature_importance_analysis = False, save_model = False, resampling_strategy="cv", fia_repeats=5, optimization=False, custom_scorer="MCC", kernel="linear", train=[], test=[]):
    
    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_SVM"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t", index_col=0, header=0)
    phenotype_df = pd.read_csv(phenotype_table, sep="\t", index_col=0, header=0)

    strains_to_be_skipped_phenotype = []
    for strain in phenotype_df.index.to_list():
        if strain not in genotype_df.index.to_list():
            strains_to_be_skipped_phenotype.append(strain)
    
    phenotype_df = phenotype_df.drop(strains_to_be_skipped_phenotype, axis=0)

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    # Get rid of uninformative strains for given antibiotic
    strains_to_be_skipped = []

    for strain in phenotype_df.index.to_list():
        if phenotype_df.loc[strain, antibiotic] == "2" or phenotype_df.loc[strain, antibiotic] == 2:
            strains_to_be_skipped.append(strain)
    
    genotype_df = genotype_df.drop(strains_to_be_skipped, axis=0)
    phenotype_df = phenotype_df.drop(strains_to_be_skipped, axis=0)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    if len(train)==0 and len(test)==0:
        X = genotype_array[:, :].astype(int)
        y = phenotype_array[:, index_of_antibiotic].astype(int)
        X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))
    
    else:
        X_train = []
        y_train = []
        X_test = []
        y_test = []

        for train_strain in train:
            if train_strain in genotype_df.index.to_list():
                X_train.append(genotype_df[:].loc[train_strain].astype(int))
                y_train.append(phenotype_df[antibiotic].loc[train_strain].astype(int))
        for test_strain in test:
            if test_strain in genotype_df.index.to_list():
                X_test.append(genotype_df[:].loc[test_strain].astype(int))
                y_test.append(phenotype_df[antibiotic].loc[test_strain].astype(int))

        X_train = np.array(X_train)
        y_train = np.array(y_train)
        X_test = np.array(X_test)
        y_test = np.array(y_test)

    # Define the hyperparameter grid to search for the best 'C' value
    param_grid = {'C': [1, 10, 100]}

    if not optimization:
        param_grid = {'C': [1]}

    # Create an SVM model
    svm_cls = SVC(class_weight={0: sum(y_train), 1: len(y_train) - sum(y_train)}, kernel=kernel)

    if custom_scorer == "MCC":
        scorer = "matthews_corrcoef"
    # Create a GridSearchCV instance to find the best hyperparameters
    grid_search = GridSearchCV(svm_cls, param_grid, cv=cv_split, scoring=scorer)
    grid_search.fit(X_train, y_train)

    # Get the best hyperparameters
    best_params = grid_search.best_params_

    # Create the final model with the best hyperparameters
    final_svm_model = SVC(class_weight={0: sum(y_train), 1: len(y_train) - sum(y_train)}, kernel=kernel, C=best_params['C'])
    final_svm_model.fit(X_train, y_train)
            
    y_hat = final_svm_model.predict(X_test)

    outfile = os.path.join(output_folder, f"{output_file_template}_Result") 

    output_file_writer(outfile, y_test, y_hat, best_c=str(best_params['C']))
    
    if feature_importance_analysis:

        r = permutation_importance(final_svm_model, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA"), "w") as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:

        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(final_svm_model, open(model_file, 'wb'))

    
def prps_ml_preprecessor(binary_mutation_table, prps_score_file, prps_percentage, temp_path):

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t")

    prps_scores = {}

    with open(prps_score_file, "r") as prps_file:
        prps_score_lines = prps_file.readlines()

    for line in prps_score_lines:
        splitted = line.split("\t")
        prps_scores[splitted[0].strip()] = float(splitted[1].strip())
    
    sorted_prps_scores = {k: v for k, v in sorted(prps_scores.items(), key=lambda item: item[1], reverse=True)}

    length_of_prps = len(sorted_prps_scores.keys())

    prps_percentage = float(prps_percentage)

    amount_of_cols_to_be_kept = (prps_percentage / 100) * length_of_prps

    cols_to_be_dropped = []

    genotype_df_columns = genotype_df.columns

    count = amount_of_cols_to_be_kept
    for key in sorted_prps_scores.keys():
        if count < 0:     
            if key in genotype_df_columns:
                cols_to_be_dropped.append(key)
            else:
                print(f"Warning: {key} is not found in the genotype table. It will be ignored.")
        count -= 1

    genotype_df_dropped = genotype_df.drop(columns=cols_to_be_dropped, axis=1)

    genotype_df_dropped.to_csv(os.path.join(temp_path, "prps_filtered_table.tsv"), sep="\t", index=False)


def gb(binary_mutation_table, phenotype_table, antibiotic, random_seed, cv_split, test_size, output_folder, n_jobs, feature_importance_analysis = False, save_model = False, resampling_strategy="holdout", fia_repeats=5, custom_scorer="MCC", n_estimators=100, max_depth=2, min_samples_leaf=1, min_samples_split=2, train=[], test=[], same_setup_run_count=1):
    
    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_GB"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t", index_col=0, header=0)
    phenotype_df = pd.read_csv(phenotype_table, sep="\t", index_col=0, header=0)

    strains_to_be_skipped_phenotype = []
    for strain in phenotype_df.index.to_list():
        if strain not in genotype_df.index.to_list():
            strains_to_be_skipped_phenotype.append(strain)
    
    phenotype_df = phenotype_df.drop(strains_to_be_skipped_phenotype, axis=0)

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    # Get rid of uninformative strains for given antibiotic
    strains_to_be_skipped = []

    for strain in phenotype_df.index.to_list():
        if phenotype_df.loc[strain, antibiotic] == "2" or phenotype_df.loc[strain, antibiotic] == 2:
            strains_to_be_skipped.append(strain)
    
    genotype_df = genotype_df.drop(strains_to_be_skipped, axis=0)
    phenotype_df = phenotype_df.drop(strains_to_be_skipped, axis=0)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    if len(train)==0 and len(test)==0:
        X = genotype_array[:, :].astype(int)
        y = phenotype_array[:, index_of_antibiotic].astype(int)
        X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))
    
    else:
        X_train = []
        y_train = []
        X_test = []
        y_test = []

        for train_strain in train:
            if train_strain in genotype_df.index.to_list():
                X_train.append(genotype_df[:].loc[train_strain].astype(int))
                y_train.append(phenotype_df[antibiotic].loc[train_strain].astype(int))
        for test_strain in test:
            if test_strain in genotype_df.index.to_list():
                X_test.append(genotype_df[:].loc[test_strain].astype(int))
                y_test.append(phenotype_df[antibiotic].loc[test_strain].astype(int))

        X_train = np.array(X_train)
        y_train = np.array(y_train)
        X_test = np.array(X_test)
        y_test = np.array(y_test)

    gb_cls = GradientBoostingClassifier(class_weight={0: sum(y_train), 1: len(y_train) - sum(y_train)}, n_estimators=n_estimators, max_depth=max_depth, min_samples_leaf=min_samples_leaf, min_samples_split=min_samples_split)

    param_grid = {
    'n_estimators': [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100],
    'max_depth': [2, 5, 7, 9]
    }

    if resampling_strategy == "cv":

        if custom_scorer == "MCC":
            scorer = "matthews_corrcoef"

        # Create a GridSearchCV instance to find the best hyperparameters
            
        grid_search = GridSearchCV(gb_cls, param_grid, cv=cv_split, scoring=scorer)
        grid_search.fit(X_train, y_train)

        y_hat = grid_search.predict(X_test)

    else:

        gb_cls.fit(X_train, y_train)
        y_hat = gb_cls.predict(X_test)

    outfile = os.path.join(output_folder, f"{output_file_template}_Result") 

    output_file_writer(outfile, y_test, y_hat)

    if feature_importance_analysis:

        r = permutation_importance(gb_cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA"), "w") as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:

        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(gb_cls, open(model_file, 'wb'))
