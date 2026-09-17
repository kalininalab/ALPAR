from pathlib import Path

GENOME_ANNOTATION_OUT_DIR = OUT_DIR / "genome_annotation"
GENOME_ANNOTATION_LOGS_DIR = GENOME_ANNOTATION_OUT_DIR / "logs"


# -----------------------
# Genome annotation
# -----------------------

rule cd_hit_create_db:
    input: FASTA_FILE
    output: GENOME_ANNOTATION_OUT_DIR / GENUS / GENUS
    log: GENOME_ANNOTATION_LOGS_DIR / "cd_hit_create_db.log"
    benchmark: BENCHMARKS_DIR / "cdhit_create_db.tsv"
    conda: ENVS_DIR.format("cd-hit")
    container: CONTAINERS.format("cd-hit:1.0.0")
    threads: workflow.cores
    shell:
        r"""
        cd-hit \
            -i {input} \
            -o {output} \
            -T {threads} \
            -M 0 \
            -g 1 \
            -s 0.8 \
            -c 0.9 \
            >> {log} 2>&1

        # Dump the cluster files
        rm -f {output}.clstr
        rm -f {output}.bak.clstr
        """


rule makeblastdb:
    input: rules.cd_hit_create_db.output
    output: touch(GENOME_ANNOTATION_OUT_DIR / "makeblastdb.done")
    log: GENOME_ANNOTATION_LOGS_DIR / "makeblastdb.log"
    benchmark: BENCHMARKS_DIR / "makeblastdb.tsv"
    conda: ENVS_DIR.format("makeblastdb")
    container: CONTAINERS.format("makeblastdb:1.0.0")
    shell:
        r"""
        makeblastdb \
            -in {input} \
            -out {input} \
            -dbtype prot \
            >> {log} 2>&1
        """


rule prokka_listdb:
    input:
        rules.makeblastdb.output,
        db_dir = rules.cd_hit_create_db.output
    output: touch(GENOME_ANNOTATION_OUT_DIR / "prokka_listdb.done"),
    log: GENOME_ANNOTATION_LOGS_DIR / "prokka_listdb.log"
    benchmark: BENCHMARKS_DIR / "prokka_listdb.tsv"
    conda: ENVS_DIR.format("prokka")
    container: CONTAINERS.format("prokka:1.0.0")
    shell:
        r"""
        DB_DIR=$(dirname {input.db_dir})
        PROKKA_DB_DIR="$(dirname $(dirname $(which prokka)))/db/genus"
        cp -a "$DB_DIR/." $PROKKA_DB_DIR

        echo $PROKKA_DB_DIR > {log}
        prokka --listdb >> {log} 2>&1
        """


#TODO Implement branching logic if no reference is given
rule prokka_runner:
    input:
        rules.prokka_listdb.output,
        sample_store = rules.rename_files.output.store,
        sample = Path(rules.rename_files.output.store) / "{sample}",
        reference = GBFF_FILE,
    output:
        gff = GENOME_ANNOTATION_OUT_DIR / "prokka" / "{sample}" / "{sample}.gff",
        faa = GENOME_ANNOTATION_OUT_DIR / "prokka" / "{sample}" / "{sample}.faa",
        gbk = GENOME_ANNOTATION_OUT_DIR / "prokka" / "{sample}" / "{sample}.gbk",
    log: GENOME_ANNOTATION_LOGS_DIR / "prokka_runner" / "{sample}.log"
    benchmark: BENCHMARKS_DIR / "prokka_{sample}.tsv"
    params:
        genus = GENUS,
        outdir = subpath(output.gff, parent=True),
    threads: 1
    resources:
        mem_mb = 600
    conda: ENVS_DIR.format("prokka")
    container: CONTAINERS.format("prokka:1.0.0")
    shell:
        r"""
        input_file=$(readlink -f {input.sample})
        prokka $input_file \
            --outdir {params.outdir} \
            --prefix {wildcards.sample} \
            --proteins {input.reference} \
            --usegenus \
            --genus {params.genus} \
            --cpus {threads} \
            --compliant \
            --force \
            >> {log} 2>&1
        """

rule genome_annotation:
    localrule: True
    input:
        lambda wc: expand(rules.prokka_runner.output, sample = get_sample_names(wc))
    output: touch(OUT_DIR / "flags" / "genome_annotation.done")
