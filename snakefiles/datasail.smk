from pathlib import Path


rule mash_sketch:
    input:
        sample_store = rules.rename_files.output.store,
        samples = lambda wildcards: expand(
            Path(rules.rename_files.output.store) / "{sample}",
            sample = get_sample_names(wildcards)
        ),
    output: TEMP_DIR / "datasail" / "mash_sketch.msh",
    log: LOGS_DIR / "mash_sketch.log",
    benchmark: BENCHMARKS_DIR / "mash_sketch.tsv",
    conda: ENVS_DIR.format("datasail"),
    threads: 8,
    shell:
        r"""
        mash \
            sketch \
            -p {threads} \
            -o {output} \
            {input.samples} \
            >> {log} 2>&1
        """

rule mash_dist:
    input: rules.mash_sketch.output,
    output: TEMP_DIR / "datasail" / "distance_matrix.tsv",
    log: LOGS_DIR / "mash_dist.log",
    benchmark: BENCHMARKS_DIR / "mash_dist.tsv",
    conda: ENVS_DIR.format("datasail"),
    threads: 8,
    shell:
        r"""
        mash \
            dist \
            -t \
            -p {threads} \
            {input} \
            {input} \
            > {output} \
            2>> {log}
        """

rule datasail_pre_processor:
    input: rules.mash_dist.output,
    output: TEMP_DIR / "datasail" / "distance_matrix_preprocessed.tsv",
    log: LOGS_DIR / "datasail_preprocessor.log",
    benchmark: BENCHMARKS_DIR / "datasail_preprocessor.tsv",
    conda: ENVS_DIR.format("miller"),
    threads: 1,
    shell:
        r"""
        mlr --tsv \
            rename -r '^.*/([^/]+)$,\1' \
            then put '${{#query}} = sub(${{#query}}, "^.*/", "")' \
            {input} > {output} 2> {log}
        """

rule datasail_runner:
    input:
        distance_matrix = rules.datasail_pre_processor.output[0],
        phenotype_dataframe = rules.phenotype_dataframe_creator.output[0],
    output: TEMP_DIR / "datasail" / "{antibiotic}" / "splits.tsv",
    log: LOGS_DIR / "datasail_runner_{antibiotic}.log",
    benchmark: BENCHMARKS_DIR / "datasail_runner_{antibiotic}.tsv",
    conda: ENVS_DIR.format("datasail"),
    params:
        techniques = "C1e",
        splits = [0.8, 0.2],
        names = ["train", "test"],
        e_type = "P",
        max_sec = 600,
        verbose = "D",
        delta = 0.1,
        epsilon = 0.1,
        runs = 1,
        solver = "SCIP",
        linkage = "average",
        e_clusters = 50,
        mock = lookup(dpath="mock", within=config, default=False),
    threads: 1,
    script:
        SCRIPTS_DIR / "datasail_runner.py"

rule split_train_test:
    input: rules.datasail_runner.output,
    output: TEMP_DIR / "datasail" / "{antibiotic}" / "{split_category}.txt",
    log: LOGS_DIR / "split_train_test_{antibiotic}_{split_category}.log",
    benchmark: BENCHMARKS_DIR / "split_train_test_{antibiotic}_{split_category}.tsv",
    threads: 1,
    shell:
        r"""
        grep -P '\t{wildcards.split_category}$' {input} | cut -f1 > {output} 2>> {log}
        """

rule split_phenotype_dataframe:
    input:
        phenotype_dataframe = rules.phenotype_dataframe_creator.output,
        split_category = rules.split_train_test.output,
    output: TEMP_DIR / "datasail" / "{antibiotic}" / "{split_category}_phenotype_dataframe.tsv",
    log: LOGS_DIR / "split_phenotype_dataframe_{antibiotic}_{split_category}.log",
    benchmark: BENCHMARKS_DIR / "split_phenotype_dataframe_{antibiotic}_{split_category}.tsv",
    conda: ENVS_DIR.format("miller"),
    threads: 1,
    shell:
        r"""
	mlr --tsv \
        	join \
	        --implicit-tsv-header \
	        -j checksum \
	        -l 1 \
	        -r checksum \
	        -f {input.split_category} \
	        {input.phenotype_dataframe} \
        	> {output} 2> {log}        
	"""

rule datasail:
    input:
        lambda wildcards: expand(
            rules.datasail_runner.output,
            antibiotic = ANTIBIOTICS
        )
    output: touch(TEMP_DIR / "flags" / "datasail.done")
