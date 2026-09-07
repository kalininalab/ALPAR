from pathlib import Path

FEATURE_TABLE_OUT_DIR = OUT_DIR / "feature_table"
FEATURE_TABLES_LOGS_DIR = FEATURE_TABLE_OUT_DIR / "logs"


rule merge_features:
    input:
        snp = rules.binary_mutation_table.output,
        gpa = rules.binary_gpa.output,
        pangenome_train = rules.gather_bubble_features.output.features,
        pangenome_test = rules.gather_bubble_features_test.output,
    output: FEATURE_TABLE_OUT_DIR / "merged_table_{antibiotic}.tsv",
    threads: 1,
    shell:
        r"""
        cat \
            {input.snp} \
            {input.gpa} \
            {input.pangenome_train} \
            {input.pangenome_test} \
            > {output}
        """

rule pivot_merged_features_miller:
    input: rules.merge_features.output,
    output: FEATURE_TABLE_OUT_DIR / "merged_table_pivot_{antibiotic}.tsv",
    log: FEATURE_TABLES_LOGS_DIR / "pivot_merged_features_miller_{antibiotic}.log",
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

rule feature_table:
    input:
        expand(rules.pivot_merged_features_miller.output, antibiotic = ANTIBIOTICS)
    output: touch(OUT_DIR / "flags" / "feature_table.done")
