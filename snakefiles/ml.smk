from pathlib import Path

ML_OUT_DIR = OUT_DIR / "ml"
ML_LOGS_DIR = ML_OUT_DIR / "logs"

rule merge_features:
    input:
        rules.merge_binary_features.output,
        rules.gather_bubble_features.output,
    output: ML_OUT_DIR / "merged_table_{antibiotic}.tsv",
    threads: 1,
    shell:
        r"""
        cat {input} > {output}
        """

rule pivot_merged_features_miller:
    input: rules.merge_features.output,
    output: ML_OUT_DIR / "merged_table_pivot_{antibiotic}.tsv",
    log: ML_LOGS_DIR / "pivot_merged_features_miller_{antibiotic}.log",
    benchmark: BENCHMARKS_DIR / "pivot_merged_features_miller_{antibiotic}.tsv",
    conda: ENVS_DIR.format("miller"),
    threads: workflow.cores,
    shell:
        r"""
        mlr --tsv --implicit-tsv-header \
            label hash,feature,value \
            then reshape -s feature,value \
            then unsparsify --fill-with '' \
            {input} > {output} 2> {log}
        """

rule prps_ml_preprocessor:
    input:
        binary_mutation_table = rules.pivot_merged_features_miller.output[0],
        prps_score_file = rules.prps_runner.output[0],
    output: ML_OUT_DIR / "prps_filtered_table.tsv",
    benchmark: BENCHMARKS_DIR / "prps_ml_preprocessor.tsv",
    log: ML_LOGS_DIR / "prps_ml_preprocessor.log",
    params:
        prps_percentage = 30
    conda: ENVS_DIR.format("python313")
    threads: 1
    script:
        SCRIPTS_DIR / "prps_ml_preprocessor.py"

rule copy_and_zip_file:
    input: rules.pivot_merged_features_miller.output[0] #TODO:prps_ml_preprocessor
    output: ML_OUT_DIR / "model_binary_mutation_table.tar.gz"
    benchmark: BENCHMARKS_DIR / "copy_and_zip_file.tsv"
    log: ML_LOGS_DIR / "copy_and_zip_file.log"
    threads: 1
    shell:
        """
        gzip -c {input} > {output} 2> {log}
        """

rule combined_ml:
    input:
        binary_mutation_table = rules.pivot_merged_features_miller.output[0],
        phenotype_table = rules.phenotype_dataframe_creator.output[0],
        train = lambda wildcards: expand(rules.split_train_test.output, split_category=["train"], **wildcards),
        test = lambda wildcards: expand(rules.split_train_test.output, split_category=["test"], **wildcards),
    output:
        best_params = ML_OUT_DIR / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}_best_params.txt",
        model_file  = ML_OUT_DIR / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}_model.sav",
        result      = ML_OUT_DIR / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}_Result.txt",
        fia         = ML_OUT_DIR / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}.txt",
    benchmark: BENCHMARKS_DIR / "combined_ml_{antibiotic}_seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}.tsv"
    log: ML_LOGS_DIR / "{antibiotic}" / "seed_{random_seed}_testsize_{test_size}_resampling_{resampling_strategy}_{model_type}_FIA_{feature_importance_analysis_strategy}.log"
    params:
        feature_importance_analysis = True,
        save_model = True,
    conda: ENVS_DIR.format("ml")
    threads: lambda wildcards: workflow.cores // len(ANTIBIOTICS)
    resources:
        mem_gb = lambda wildcards: workflow.global_resources.get("mem_gb", 4) // len(ANTIBIOTICS),
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
            feature_importance_analysis_strategy=["gini"]
        )
    output: touch(OUT_DIR / "flags" / "ml_runner.done")
