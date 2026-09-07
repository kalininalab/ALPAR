from pathlib import Path

SNP_OUT_DIR = OUT_DIR / "snp"
SNP_LOGS_DIR = SNP_OUT_DIR / "logs"


# -----------------------
# Mutation Analysis
# -----------------------

rule snippy_runner:
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
        out_dir = subpath(output.vcf, parent=True)
    conda: ENVS_DIR.format("snippy")
    threads: 1
    resources:
        mem_gb = 1
    shell:
        r"""
        snippy \
            --ctgs {input.sample} \
            --outdir {params.out_dir} \
            --reference {input.reference} \
            --cpus {threads} \
            --ram {resources.mem_gb} \
            --force \
            >> {log} 2>&1
        """


# -----------------------
# Feature Importance Analysis
# -----------------------

rule annotation_file_from_snippy:
    input:
        lambda wc: expand(
            rules.snippy_runner.output.tab,
            sample = get_sample_names(wc)
        )
    output: SNP_OUT_DIR / "mutations_annotations.tsv"
    benchmark: BENCHMARKS_DIR / "annotation_file_from_snippy.tsv"
    log: SNP_LOGS_DIR / "annotation_file_from_snippy.log"
    conda: ENVS_DIR.format("python313")
    threads: MAX_PYTHON_THREADS
    script:
        SCRIPTS_DIR / "annotation_file_from_snippy.py"


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
    threads: 1
    script:
        SCRIPTS_DIR / "binary_mutation_table.py"


rule snp:
    input:
        rules.annotation_file_from_snippy.output,
        rules.binary_mutation_table.output
    output: touch(OUT_DIR / "flags" / "snp.done")
