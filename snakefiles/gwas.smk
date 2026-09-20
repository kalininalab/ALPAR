from pathlib import Path

GWAS_OUT_DIR = OUT_DIR / "gwas"
GWAS_LOGS_DIR = GWAS_OUT_DIR / "logs"

rule pyseer_genotype_matrix_creator:
    input: rules.binary_mutation_table.output # GWAS only on SNP
    output: GWAS_OUT_DIR / "genotype_matrix.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_genotype_matrix_creator.tsv"
    log: GWAS_LOGS_DIR / "pyseer_genotype_matrix_creator.log"
    conda: ENVS_DIR.format("miller")
    container: CONTAINERS.format("miller:1.0.0")
    threads: workflow.cores
    shell:
        r"""
        mlr --tsv --implicit-tsv-header \
            label hash,feature,value \
            then reshape -s hash,value \
            then unsparsify --fill-with 0 \
            {input} > {output} 2> {log}
        """

rule pyseer_phenotype_file_creator:
    group: "pyseer_phenotype_file_creator_batch"
    input: rules.phenotype_dataframe_creator.output
    output: GWAS_OUT_DIR / "pyseer_phenotype_file_{antibiotic}.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_phenotype_file_creator_{antibiotic}.tsv"
    log: GWAS_LOGS_DIR / "pyseer_phenotype_file_creator_{antibiotic}.log"
    conda: ENVS_DIR.format("miller")
    container: CONTAINERS.format("miller:1.0.0")
    threads: 1
    shell:
        r"""
        first_col=$(head -1 {input} | cut -f1)
        mlr --tsv \
            rename "$first_col,samples" \
            then cut -f samples,{wildcards.antibiotic} \
            then filter '$["{wildcards.antibiotic}"] == 0 || $["{wildcards.antibiotic}"] == 1' \
            {input} > {output} 2> {log}
        """

rule pyseer_similarity_matrix_creator:
    input:
        phylogeny = rules.mashtree_runner.output[0],
    output: GWAS_OUT_DIR / "similarity_matrix.tsv"
    params:
        output_format = "newick",
        midpoint = False,
        method = "lmm", # topology
    benchmark: BENCHMARKS_DIR / "pyseer_similarity_matrix_creator.tsv"
    log: GWAS_LOGS_DIR / "pyseer_similarity_matrix_creator.log"
    conda: ENVS_DIR.format("gwas")
    container: CONTAINERS.format("gwas:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "phylogeny_distance.py"

rule pyseer_runner:
    group: "pyseer_runner_batch"
    input:
        phenotype = rules.pyseer_phenotype_file_creator.output,
        genotype = rules.pyseer_genotype_matrix_creator.output,
        similarity_matrix = rules.pyseer_similarity_matrix_creator.output,
    output: GWAS_OUT_DIR / "pyseer_results" / "{antibiotic}.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_runner_{antibiotic}.tsv"
    log: GWAS_LOGS_DIR / "pyseer_runner_{antibiotic}.log"
    conda: ENVS_DIR.format("pyseer")
    container: CONTAINERS.format("pyseer:1.0.0")
    threads: workflow.cores
    shell:
        r"""
        pyseer \
            --lmm \
            --phenotypes {input.phenotype} \
            --pres {input.genotype} \
            --similarity {input.similarity_matrix} \
            --cpu {threads} \
            > {output} 2> {log}
        """

rule pyseer_post_processor_sort:
    group: "pyseer_post_processor_sort_batch"
    input: rules.pyseer_runner.output,
    output: GWAS_OUT_DIR / "pyseer_results_sorted" / "{antibiotic}.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_post_processor_{antibiotic}.tsv"
    log: GWAS_LOGS_DIR / "pyseer_post_processor_{antibiotic}.log"
    conda: ENVS_DIR.format("miller")
    container: CONTAINERS.format("miller:1.0.0")
    threads: 1
    shell:
        r"""
        mlr --tsv \
            sort -n lrt-pvalue \
            {input} > {output} 2> {log}
        """

rule pyseer_post_processor_clean:
    group: "pyseer_post_processor_clean_batch"
    input: rules.pyseer_post_processor_sort.output,
    output: GWAS_OUT_DIR / "pyseer_results_sorted_cleaned" / "{antibiotic}.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_post_processor_clean_{antibiotic}.tsv"
    log: GWAS_LOGS_DIR / "pyseer_post_processor_clean_{antibiotic}.log"
    params:
        mock = lookup(dpath="mock", within=config, default=False),
    conda: ENVS_DIR.format("miller")
    container: CONTAINERS.format("miller:1.0.0")
    threads: 1
    shell:
        r"""
        mlr --tsv \
            then filter '$[NF] != "bad-chisq"' \
            {input} > {output} 2> {log}
        
        if [ $(wc -l < {output}) -eq 0 ]; then
            echo "No valid results after cleaning {input}" >> {log}
            head -1 {input} > {output}
        fi

        if [ "{params.mock}" = "True" ]; then
            echo "Mock mode enabled" >> {log}
            cat {input} > {output}
        fi
        """

rule pyseer_gwas_graph_creator:
    group: "pyseer_gwas_graph_creator_batch"
    input:
        gwas_results = rules.pyseer_post_processor_clean.output[0],
        gwas_postprocessed = rules.pyseer_post_processor_sort.output[0],
    output: GWAS_OUT_DIR / "graphs" / "{antibiotic}.jpg"
    log: GWAS_LOGS_DIR / "pyseer_gwas_graph_creator_{antibiotic}.log"
    conda: ENVS_DIR.format("gwas")
    container: CONTAINERS.format("gwas:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "pyseer_gwas_graph_creator.py"

rule decision_tree_input_creator:
    group: "decision_tree_input_creator_batch"
    input:
        binary_table = rules.pyseer_genotype_matrix_creator.input[0],
        phenotype_file = rules.pyseer_phenotype_file_creator.input[0],
        pyseer_output_raw = rules.pyseer_post_processor_clean.output[0],
        pyseer_output_sorted_cleaned = rules.pyseer_post_processor_clean.output[0],
    output:
        tree_result = GWAS_OUT_DIR / "decision_tree" / "{antibiotic}_result.txt",
        tree_model = GWAS_OUT_DIR / "decision_tree" / "{antibiotic}_model.pkl",
    params:
        antibiotic = lambda wildcards: wildcards.antibiotic,
    benchmark: BENCHMARKS_DIR / "decision_tree_input_creator_{antibiotic}.tsv"
    log: GWAS_LOGS_DIR / "decision_tree_input_creator_{antibiotic}.log"
    conda: ENVS_DIR.format("gwas")
    container: CONTAINERS.format("gwas:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "decision_tree_input_creator.py"

rule gwas:
    localrule: True
    input:
        expand(rules.pyseer_gwas_graph_creator.output, antibiotic=ANTIBIOTICS),
        expand(rules.decision_tree_input_creator.output.tree_result, antibiotic=ANTIBIOTICS),
        expand(rules.decision_tree_input_creator.output.tree_model, antibiotic=ANTIBIOTICS),
    output: touch(OUT_DIR / "flags" / "gwas.done")
