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
    input: rules.pivot_merged_features_miller.output[0]
    output: TEMP_DIR / "ml" / "model_binary_mutation_table.tar.gz"
    benchmark: BENCHMARKS_DIR / "copy_and_zip_file.tsv"
    log: LOGS_DIR / "ml" / "copy_and_zip_file.log"
    threads: 1
    shell:
        """
        gzip -c {input} > {output} 2> {log}
        """
