from pathlib import Path

FEATURE_TABLE_OUT_DIR = OUT_DIR / "feature_table"
FEATURE_TABLES_LOGS_DIR = FEATURE_TABLE_OUT_DIR / "logs"


rule merge_features:
    localrule: True
    input:
        snp = rules.binary_mutation_table.output,
        gpa = rules.binary_gpa.output,
        pangenome_train = rules.merge_bubble_features.output,
        pangenome_test = rules.merge_bubble_features_test.output,
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
    group: "pivot_merged_features_miller_batch"
    input: rules.merge_features.output,
    output: FEATURE_TABLE_OUT_DIR / "merged_table_pivot_{antibiotic}.tsv",
    log: FEATURE_TABLES_LOGS_DIR / "pivot_merged_features_miller_{antibiotic}.log",
    benchmark: BENCHMARKS_DIR / "pivot_merged_features_miller_{antibiotic}.tsv",
    conda: ENVS_DIR.format("miller"),
    container: CONTAINERS.format("miller:1.0.0")
    threads: workflow.cores,
    params:
        # The rule keeps its historical name. Miller's in-memory reshape cannot
        # hold the multi-GB merged table, so the script streams instead.
        script = SCRIPTS_DIR / "long_to_wide.sh",
        sort_memory = lambda wc, resources: f"{max(resources.mem_mb // 2, 500)}M",
        # Sort spills about the input size; keep it on shared storage beside the output.
        sort_tmp = subpath(output[0], parent=True),
    shell:
        r"""
        # Rows are samples (field 1), columns are features (field 2), empty fill.
        bash {params.script:q} {input:q} {output:q} {params.sort_tmp:q} {threads} {params.sort_memory} \
            1 2 hash '' 2> {log:q}
        """

rule feature_table:
    localrule: True
    input:
        expand(rules.pivot_merged_features_miller.output, antibiotic = ANTIBIOTICS)
    output: touch(OUT_DIR / "flags" / "feature_table.done")
