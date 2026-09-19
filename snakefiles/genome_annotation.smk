from pathlib import Path

GENOME_ANNOTATION_OUT_DIR = OUT_DIR / "genome_annotation"
GENOME_ANNOTATION_LOGS_DIR = GENOME_ANNOTATION_OUT_DIR / "logs"


# -----------------------
# Genome annotation
# -----------------------

rule cd_hit_create_db:
    input: FASTA_FILE
    # Prokka normalizes --genus with ucfirst(lc(...)) before looking up the DB.
    output: GENOME_ANNOTATION_OUT_DIR / "prokka_db" / "genus" / GENUS.lower().capitalize()
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
    output: touch(GENOME_ANNOTATION_OUT_DIR / "prokka_db" / "makeblastdb.done")
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


#TODO Implement branching logic if no reference is given
rule prokka_runner:
    group: "prokka_batch"
    input:
        database = rules.makeblastdb.output,
        genus_database = rules.cd_hit_create_db.output,
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
        genus = GENUS.lower().capitalize(),
        genus_db_dir = GENOME_ANNOTATION_OUT_DIR / "prokka_db" / "genus",
        outdir = subpath(output.gff, parent=True),
    threads: 1
    resources:
        mem_mb = 600
    conda: ENVS_DIR.format("prokka")
    container: CONTAINERS.format("prokka:1.0.0")
    shell:
        r"""
        input_file=$(readlink -f {input.sample})

        # --dbdir replaces Prokka's complete database root. Build a writable,
        # job-local view containing the bundled DBs and our shared genus DB.
        bundled_db_dir="$(dirname "$(dirname "$(command -v prokka)")")/db"
        job_db_dir="$(mktemp -d)"
        trap 'rm -rf "$job_db_dir"' EXIT

        for bundled_db in "$bundled_db_dir"/*; do
            if [[ "$(basename "$bundled_db")" != "genus" ]]; then
                ln -s "$bundled_db" "$job_db_dir/$(basename "$bundled_db")"
            fi
        done
        ln -s {params.genus_db_dir:q} "$job_db_dir/genus"

        prokka "$input_file" \
            --dbdir "$job_db_dir" \
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
