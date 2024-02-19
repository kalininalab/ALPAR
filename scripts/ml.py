import sys
import autosklearn
import autosklearn.classification
import sklearn.model_selection
import sklearn.datasets
import sklearn.metrics
from sklearn.inspection import permutation_importance
from autosklearn.experimental.askl2 import AutoSklearn2Classifier
import pandas as pd
from pprint import pprint
import numpy as np
import pickle
from random import randint
import os
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

mcc_scorer = autosklearn.metrics.make_scorer(
    "mcc",
    sklearn.metrics.matthews_corrcoef,
    greater_is_better=True
)


def output_file_writer(outfile, y_test, y_hat, cls = None, best_c = None):

    with open(outfile, "w") as ofile:
        if cls:
            ofile.write(str(cls.show_models()))
            ofile.write("\n")
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


def rf_auto_ml(binary_mutation_table, phenotype_table, antibiotic, random_seed, cv_split, test_size, output_folder, n_jobs, temp_folder, ram, optimization_time_limit, feature_importance_analysis = False, save_model = False, resampling_strategy="holdout", custom_scorer="MCC", fia_repeats=5):

    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_RF_AutoML"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t")
    phenotype_df = pd.read_csv(phenotype_table, sep="\t")

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    X = genotype_array[:, 1:].astype(int)
    y = phenotype_array[:, index_of_antibiotic].astype(int)

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))

    if custom_scorer == "MCC":
        scorer = mcc_scorer
    else:
        scorer = custom_scorer

    classifier = "random_forest"

    if resampling_strategy == "cv":
        resampling_strategy_arguments = {"folds": cv_split, 'train_size': float(1.00) - float(test_size)}

    elif resampling_strategy == "holdout":
        resampling_strategy_arguments = {"train_size": float(1.00) - float(test_size)}

    cls = autosklearn.classification.AutoSklearnClassifier(
        memory_limit=float(ram) * 1024,
        time_left_for_this_task=optimization_time_limit,
        include={'classifier': [classifier]},
        resampling_strategy=f'{resampling_strategy}',
        resampling_strategy_arguments=resampling_strategy_arguments,
        delete_tmp_folder_after_terminate=False,
        ensemble_size=1,
        metric=scorer,
        tmp_folder=os.path.join(temp_folder, f"{output_file_template}_temp"),
        n_jobs=n_jobs
    )
    cls.fit(X_train, y_train)
    y_hat = cls.predict(X_test)

    outfile = os.path.join(output_folder, f"{output_file_template}_Result") 

    if feature_importance_analysis:

        r = permutation_importance(cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA", "w")) as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:
        
        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(cls, open(model_file, 'wb'))

    output_file_writer(outfile, y_test, y_hat, cls)


def rf(binary_mutation_table, phenotype_table, antibiotic, random_seed, cv_split, test_size, output_folder, n_jobs, feature_importance_analysis = False, save_model = False, resampling_strategy="holdout", fia_repeats=5, custom_scorer="MCC", n_estimators=100, max_depth=2, min_samples_leaf=1, min_samples_split=2):
    
    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_RF"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t")
    phenotype_df = pd.read_csv(phenotype_table, sep="\t")

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    X = genotype_array[:, 1:].astype(int)
    y = phenotype_array[:, index_of_antibiotic].astype(int)

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))

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

    if feature_importance_analysis:

        r = permutation_importance(rf_cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA", "w")) as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:

        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(rf_cls, open(model_file, 'wb'))

    output_file_writer(outfile, y_test, y_hat)



def svm(binary_mutation_table, phenotype_table, antibiotic, random_seed, test_size, output_folder, n_jobs, feature_importance_analysis = False, save_model = False, resampling_strategy="holdout", fia_repeats=5, optimization=False, kernel="linear"):

    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_SVM"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t")
    phenotype_df = pd.read_csv(phenotype_table, sep="\t")

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    X = genotype_array[:, 1:].astype(int)
    y = phenotype_array[:, index_of_antibiotic].astype(int)

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))

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
    
    if feature_importance_analysis:

        r = permutation_importance(best_model, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA", "w")) as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:

        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(best_model, open(model_file, 'wb'))

    output_file_writer(outfile, y_test, y_hat, best_c=bm_c)
    

def svm_cv(binary_mutation_table, phenotype_table, antibiotic, random_seed, test_size, output_folder, n_jobs, cv_split, feature_importance_analysis = False, save_model = False, resampling_strategy="cv", fia_repeats=5, optimization=False, custom_scorer="MCC", kernel="linear"):
    
    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_SVM"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t")
    phenotype_df = pd.read_csv(phenotype_table, sep="\t")

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    X = genotype_array[:, 1:].astype(int)
    y = phenotype_array[:, index_of_antibiotic].astype(int)

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))

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
    
    if feature_importance_analysis:

        r = permutation_importance(final_svm_model, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA", "w")) as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:

        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(final_svm_model, open(model_file, 'wb'))

    output_file_writer(outfile, y_test, y_hat, best_c=str(best_params['C']))

    
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


def gb_auto_ml(binary_mutation_table, phenotype_table, antibiotic, random_seed, cv_split, test_size, output_folder, n_jobs, temp_folder, ram, optimization_time_limit, feature_importance_analysis = False, save_model = False, resampling_strategy="holdout", custom_scorer="MCC", fia_repeats=5):

    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_GB_AutoML"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t")
    phenotype_df = pd.read_csv(phenotype_table, sep="\t")

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    X = genotype_array[:, 1:].astype(int)
    y = phenotype_array[:, index_of_antibiotic].astype(int)

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))

    if custom_scorer == "MCC":
        scorer = mcc_scorer
    else:
        scorer = custom_scorer

    classifier = "gradient_boosting"

    if resampling_strategy == "cv":
        resampling_strategy_arguments = {"folds": cv_split, 'train_size': float(1.00) - float(test_size)}

    elif resampling_strategy == "holdout":
        resampling_strategy_arguments = {"train_size": float(1.00) - float(test_size)}

    cls = autosklearn.classification.AutoSklearnClassifier(
        memory_limit=float(ram) * 1024,
        time_left_for_this_task=optimization_time_limit,
        include={'classifier': [classifier]},
        resampling_strategy=f'{resampling_strategy}',
        resampling_strategy_arguments=resampling_strategy_arguments,
        delete_tmp_folder_after_terminate=False,
        ensemble_size=1,
        metric=scorer,
        tmp_folder=os.path.join(temp_folder, f"{output_file_template}_temp"),
        n_jobs=n_jobs
    )
    cls.fit(X_train, y_train)
    y_hat = cls.predict(X_test)

    outfile = os.path.join(output_folder, f"{output_file_template}_Result") 

    if feature_importance_analysis:

        r = permutation_importance(cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA", "w")) as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:
        
        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(cls, open(model_file, 'wb'))

    output_file_writer(outfile, y_test, y_hat, cls)


def gb(binary_mutation_table, phenotype_table, antibiotic, random_seed, cv_split, test_size, output_folder, n_jobs, feature_importance_analysis = False, save_model = False, resampling_strategy="holdout", fia_repeats=5, custom_scorer="MCC", n_estimators=100, max_depth=2, min_samples_leaf=1, min_samples_split=2):
    
    output_file_template = f"seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_GB"

    genotype_df = pd.read_csv(binary_mutation_table, sep="\t")
    phenotype_df = pd.read_csv(phenotype_table, sep="\t")

    # Make sure rows are matching
    phenotype_df = phenotype_df.reindex(genotype_df.index)

    genotype_array = genotype_df.to_numpy()
    phenotype_array = phenotype_df.to_numpy()

    index_of_antibiotic = phenotype_df.columns.get_loc(antibiotic)

    X = genotype_array[:, 1:].astype(int)
    y = phenotype_array[:, index_of_antibiotic].astype(int)

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=random_seed, test_size=float(test_size))

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

    if feature_importance_analysis:

        r = permutation_importance(gb_cls, X_test, y_test, n_repeats=fia_repeats, random_state=random_seed, n_jobs=n_jobs)

        with open(os.path.join(output_folder, f"{output_file_template}_FIA", "w")) as ofile:
            for i in r.importances_mean.argsort()[::-1]:
                if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
                    ofile.write(f"{genotype_df.columns[i]:<8};{r.importances_mean[i]:.3f};+/-{r.importances_std[i]:.3f}\n")

    if save_model:

        model_file = os.path.join(output_folder, f"{output_file_template}_model.sav")
        pickle.dump(gb_cls, open(model_file, 'wb'))

    output_file_writer(outfile, y_test, y_hat)
