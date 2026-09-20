from pathlib import Path

PRPS_OUT_DIR = OUT_DIR / "prps"
PRPS_LOGS_DIR = PRPS_OUT_DIR / "logs"

rule prps_runner:
    group: "prps_runner_batch"
    input:
        phylogeny_tree = rules.mashtree_runner.output[0],
        feature_matrix = rules.pivot_merged_features_miller.output[0],
    output: PRPS_OUT_DIR / "prps_scores_{antibiotic}.tsv",
    benchmark: BENCHMARKS_DIR / "prps_runner_{antibiotic}.tsv",
    log: PRPS_LOGS_DIR / "prps_runner_{antibiotic}.log",
    conda: ENVS_DIR.format("prps"),
    container: CONTAINERS.format("prps:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "prps.py"

rule prps:
    localrule: True
    input: lambda wc: expand(rules.prps_runner.output, sample = get_sample_names(wc), antibiotic = ANTIBIOTICS)
    output: touch(OUT_DIR / "flags" / "prps.done")
