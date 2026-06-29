from pathlib import Path

rule pyseer_genotype_matrix_creator:
    input: rules.binary_mutation_table.output # GWAS only on SNP
    output: TEMP_DIR/ "gwas" / "genotype_matrix.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_genotype_matrix_creator.tsv"
    log: LOGS_DIR / "gwas" / "pyseer_genotype_matrix_creator.log"
    conda: ENVS_DIR.format("miller")
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
    input: rules.phenotype_dataframe_creator.output
    output: TEMP_DIR / "gwas" / "pyseer_phenotype_file_{antibiotic}.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_phenotype_file_creator_{antibiotic}.tsv"
    log: LOGS_DIR / "gwas" / "pyseer_phenotype_file_creator_{antibiotic}.log"
    conda: ENVS_DIR.format("miller")
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
    output: TEMP_DIR / "gwas" / "similarity_matrix.tsv"
    params:
        output_format = "newick",
        midpoint = False,
        method = "lmm", # topology
    benchmark: BENCHMARKS_DIR / "pyseer_similarity_matrix_creator.tsv"
    log: LOGS_DIR / "gwas" / "pyseer_similarity_matrix_creator.log"
    conda: ENVS_DIR.format("gwas")
    threads: 1
    script:
        SCRIPTS_DIR / "phylogeny_distance.py"

rule pyseer_runner:
    input:
        phenotype = rules.pyseer_phenotype_file_creator.output,
        genotype = rules.pyseer_genotype_matrix_creator.output,
        similarity_matrix = rules.pyseer_similarity_matrix_creator.output,
    output: TEMP_DIR / "gwas" / "pyseer_results" / "{antibiotic}.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_runner_{antibiotic}.tsv"
    log: LOGS_DIR / "gwas" / "pyseer_runner_{antibiotic}.log"
    conda: ENVS_DIR.format("pyseer")
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
    input: rules.pyseer_runner.output,
    output: TEMP_DIR / "gwas" / "pyseer_results_sorted" / "{antibiotic}.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_post_processor_sort_{antibiotic}.tsv"
    log: LOGS_DIR / "gwas" / "pyseer_post_processor_sort_{antibiotic}.log"
    conda: ENVS_DIR.format("miller")
    threads: 1
    shell:
        r"""
        mlr --tsv \
            sort -n lrt-pvalue \
            {input} > {output} 2> {log}
        
        if [ ! -f {output} ]; then
            head -1 {input} > {output}
        fi
        """

rule pyseer_post_processor_clean:
    input: rules.pyseer_runner.output,
    output: TEMP_DIR / "gwas" / "pyseer_results_sorted_cleaned" / "{antibiotic}.tsv"
    benchmark: BENCHMARKS_DIR / "pyseer_post_processor_clean_{antibiotic}.tsv"
    log: LOGS_DIR / "gwas" / "pyseer_post_processor_clean_{antibiotic}.log"
    conda: ENVS_DIR.format("miller")
    threads: 1
    shell:
        r"""
        mlr --tsv \
            then filter '$[NF] != "bad-chisq"' \
            {input} > {output} 2> {log}
        
        if [ ! -f {output} ]; then
            head -1 {input} > {output}
        fi
        """

rule pyseer_gwas_graph_creator:
    input:
        gwas_results = rules.pyseer_runner.output[0],
        gwas_postprocessed = rules.pyseer_post_processor_sort.output[0],
    output: TEMP_DIR / "gwas" / "graphs" / "{antibiotic}.jpg"
    log: LOGS_DIR / "gwas" / "pyseer_gwas_graph_creator_{antibiotic}.log"
    conda: ENVS_DIR.format("gwas")
    threads: 1
    script:
        SCRIPTS_DIR / "pyseer_gwas_graph_creator.py"

rule decision_tree_input_creator:
    input:
        binary_table = rules.pyseer_genotype_matrix_creator.input[0],
        phenotype_file = rules.pyseer_phenotype_file_creator.input[0],
        pyseer_output_raw = rules.pyseer_runner.output[0],
        pyseer_output_sorted_cleaned = rules.pyseer_post_processor_clean.output[0],
    output:
        tree_result = TEMP_DIR / "gwas" / "decision_tree" / "{antibiotic}_result.txt",
        tree_model = TEMP_DIR / "gwas" / "decision_tree" / "{antibiotic}_model.pkl",
    params:
        antibiotic = lambda wildcards: wildcards.antibiotic,
    benchmark: BENCHMARKS_DIR / "decision_tree_input_creator_{antibiotic}.tsv"
    log: LOGS_DIR / "gwas" / "decision_tree_input_creator_{antibiotic}.log"
    conda: ENVS_DIR.format("gwas")
    threads: 1
    script:
        SCRIPTS_DIR / "decision_tree_input_creator.py"

rule gwas:
    input:
        expand(rules.pyseer_gwas_graph_creator.output, antibiotic=ANTIBIOTICS),
        expand(rules.decision_tree_input_creator.output.tree_result, antibiotic=ANTIBIOTICS),
        expand(rules.decision_tree_input_creator.output.tree_model, antibiotic=ANTIBIOTICS),
    output: touch(TEMP_DIR / "flags" / "gwas.done")
