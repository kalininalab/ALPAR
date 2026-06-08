from pathlib import Path

rule prps_runner:
    input:
        phylogeny_tree = rules.mashtree_runner.output[0],
        feature_matrix = rules.pivot_merged_features_miller.output[0],
    output: TEMP_DIR / "prps_scores.tsv",
    benchmark: BENCHMARKS_DIR / "prps_runner.tsv",
    log: LOGS_DIR / "prps_runner.log"
    conda: ENVS_DIR.format("prps")
    threads: 1
    script:
        SCRIPTS_DIR / "prps.py"

rule prps:
    input: rules.prps_runner.output
    output: touch(TEMP_DIR / "flags" / "prps.done")
