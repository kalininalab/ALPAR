from pathlib import Path

PRPS_OUT_DIR = OUT_DIR / "prps"
PRPS_LOGS_DIR = PRPS_OUT_DIR / "logs"

rule prps_runner:
    input:
        phylogeny_tree = rules.mashtree_runner.output[0],
        feature_matrix = rules.merge_binary_features.output[0],
    output: PRPS_OUT_DIR / "prps_scores.tsv",
    benchmark: BENCHMARKS_DIR / "prps_runner.tsv",
    log: PRPS_LOGS_DIR / "prps_runner.log",
    conda: ENVS_DIR.format("prps"),
    threads: 1
    script:
        SCRIPTS_DIR / "prps.py"

rule prps:
    input: rules.prps_runner.output
    output: touch(OUT_DIR / "flags" / "prps.done")
