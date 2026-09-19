from pathlib import Path

SNP_OUT_DIR = OUT_DIR / "snp"
SNP_LOGS_DIR = SNP_OUT_DIR / "logs"


# -----------------------
# Mutation Analysis
# -----------------------

rule snippy_runner:
    group: "snippy_batch"
    input:
        sample_store = rules.rename_files.output.store,
        sample = Path(rules.rename_files.output.store) / "{sample}",
        reference = GBFF_FILE,
    output:
        vcf = SNP_OUT_DIR / "snippy" / "{sample}" / "snps.vcf",
        tab = SNP_OUT_DIR / "snippy" / "{sample}" / "snps.tab",
    log: SNP_LOGS_DIR / "snippy_runner" / "{sample}.log"
    benchmark: BENCHMARKS_DIR / "snippy_{sample}.tsv"
    params:
        out_dir = subpath(output.vcf, parent=True),
        ram_gb = lambda wildcards, resources: resources.mem_mb // 1000
    conda: ENVS_DIR.format("snippy")
    container: CONTAINERS.format("snippy:1.0.0")
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        r"""
        snippy \
            --ctgs {input.sample} \
            --outdir {params.out_dir} \
            --reference {input.reference} \
            --cpus {threads} \
            --ram {params.ram_gb} \
            --force \
            >> {log} 2>&1
        """


# -----------------------
# Feature Importance Analysis
# -----------------------

rule annotation_file_from_snippy:
    group: "annotation_file_from_snippy_batch"
    input: rules.snippy_runner.output.tab
    output: SNP_OUT_DIR / "annotation_file_from_snippy" / "{sample}.tsv"
    benchmark: BENCHMARKS_DIR / "annotation_file_from_snippy_{sample}.tsv"
    log: SNP_LOGS_DIR / "annotation_file_from_snippy" / "{sample}.log"
    conda: ENVS_DIR.format("miller")
    container: CONTAINERS.format("miller:1.0.0")
    threads: 1
    shell:
        r"""
        mlr --tsv \
            put '$Mutation = $POS . "," . $REF . ":" . $ALT . "," . $TYPE' \
            then cut -o -f Mutation,EFFECT,GENE,PRODUCT \
            {input} > {output} 2> {log}
        """


rule gather_annotation_file_from_snippy:
    input:
        lambda wc: expand(
            rules.annotation_file_from_snippy.output,
            sample = get_sample_names(wc)
        )
    output: SNP_OUT_DIR / "mutations_annotations.tsv"
    benchmark: BENCHMARKS_DIR / "gather_annotation_file_from_snippy.tsv"
    log: SNP_LOGS_DIR / "gather_annotation_file_from_snippy.log"
    conda: ENVS_DIR.format("miller")
    container: CONTAINERS.format("miller:1.0.0")
    threads: 1
    shell:
        r"""
        mlr --tsv \
            cat \
            then head -n 1 -g Mutation \
            {input} > {output} 2> {log}
        """


rule binary_mutation_table:
    input:
        lambda wc: expand(
            rules.snippy_runner.output.vcf,
            sample = get_sample_names(wc)
        )
    output: SNP_OUT_DIR / "binary_mutation_table.tsv"
    benchmark: BENCHMARKS_DIR / "binary_mutation_table.tsv"
    log: SNP_LOGS_DIR / "binary_mutation_table.log"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "binary_mutation_table.py"


rule snp:
    localrule: True
    input:
        rules.gather_annotation_file_from_snippy.output,
        rules.binary_mutation_table.output
    output: touch(OUT_DIR / "flags" / "snp.done")
