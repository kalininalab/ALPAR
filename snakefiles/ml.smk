from pathlib import Path

rule prps_ml_preprocessor:
    input:
        binary_mutation_table = rules.pivot_merged_features_miller.output[0],
        prps_score_file = rules.prps_runner.output[0],
    output: TEMP_DIR / "prps_filtered_table.tsv",
    benchmark: BENCHMARKS_DIR / "prps_ml_preprocessor.tsv",
    log: LOGS_DIR / "ml" / "prps_ml_preprocessor.log",
    params:
        prps_percentage = 30
    conda: ENVS_DIR.format("python313")
    threads: 1
    script:
        SCRIPTS_DIR / "prps_ml_preprocessor.py"

rule copy_and_zip_file:
    input: rules.pivot_merged_features_miller.output[0] #TODO:prps_ml_preprocessor
    output: TEMP_DIR / "ml" / "model_binary_mutation_table.tar.gz"
    benchmark: BENCHMARKS_DIR / "copy_and_zip_file.tsv"
    log: LOGS_DIR / "ml" / "copy_and_zip_file.log"
    threads: 1
    shell:
        """
        gzip -c {input} > {output} 2> {log}
        """

rule combined_ml:
    input:
        binary_mutation_table = rules.pivot_merged_features_miller.output[0],
        phenotype_dataframe   = rules.phenotype_dataframe_creator.output[0],
    output:
        best_params     = TEMP_DIR / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_best_params.txt",
        model_file      = TEMP_DIR / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_model.sav",
        result          = TEMP_DIR / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_Result.txt",
        fia_permutation = TEMP_DIR / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_permutation_importance.txt",
        fia_weights     = TEMP_DIR / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_absolute_feature_weights.txt",
        fia_gini        = TEMP_DIR / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_gini.txt",
    benchmark: BENCHMARKS_DIR / "combined_ml_{antibiotic}_seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}.tsv"
    log:       LOGS_DIR / "ml" / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}.log"
    params:
        feature_importance_analysis = True,
        save_model = True,
        feature_importance_analysis_strategy = "gini",
    conda: ENVS_DIR.format("python313")
    threads: 1
    resources:
        mem_gb = lambda wildcards: workflow.global_resources.get("mem_gb", 4),
    script:
        SCRIPTS_DIR / "combined_ml.py"

rule ml:
    input:
        expand(
            rules.combined_ml.output,
            antibiotic=ANTIBIOTICS,
            random_seed=[42],
            test_size=[0.2],
            resampling_strategy=["cv"],
            model_type=["xgb"],
        )
    output: touch(TEMP_DIR / "flags" / "ml_runner.done")
