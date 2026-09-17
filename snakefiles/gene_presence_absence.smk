from pathlib import Path

GENE_PRESENCE_ABSENCE_OUT_DIR = OUT_DIR / "gene_presence_absence"
GENE_PRESENCE_ABSENCE_LOGS_DIR = GENE_PRESENCE_ABSENCE_OUT_DIR / "logs"


# -----------------------
# Gene Presence Absence: Panaroo
# -----------------------

rule panaroo_runner:
    input: 
        lambda wc: expand(
            rules.prokka_runner.output.gff,
            sample = get_sample_names(wc)
        ),
    output:
        gpa = GENE_PRESENCE_ABSENCE_OUT_DIR / "panaroo" / "gene_presence_absence.csv",
        gene_data = GENE_PRESENCE_ABSENCE_OUT_DIR / "panaroo" / "gene_data.csv",
    log: GENE_PRESENCE_ABSENCE_LOGS_DIR / "panaroo_runner.log"
    benchmark: BENCHMARKS_DIR / "panaroo.tsv"
    params:
        outdir = subpath(output.gpa, parent=True),
        seq_identity_threshold = 0.8,
        seq_len_diff_cutoff = 0.8,
    conda: ENVS_DIR.format("panaroo")
    container: CONTAINERS.format("panaroo:1.0.0")
    threads: 30
    shell:
        r"""
        panaroo \
            --input {input} \
            --out_dir {params.outdir} \
            --threshold {params.seq_identity_threshold} \
            --len_dif_percent {params.seq_len_diff_cutoff} \
            --codons \
            --clean-mode strict \
            --threads {threads} \
            >> {log} 2>&1
        """


rule binary_gpa_panaroo:
    input: rules.panaroo_runner.output.gpa,
    output: GENE_PRESENCE_ABSENCE_OUT_DIR / "binary_gpa_panaroo.tsv"
    benchmark: BENCHMARKS_DIR / "binary_gpa_panaroo.py.tsv"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "binary_gpa_panaroo.py"


# -----------------------
# Gene Presence Absence: CD-HIT
# -----------------------

rule cdhit_protein_name_corrector:
    localrule: True
    input: rules.prokka_runner.output.faa,
    output: GENE_PRESENCE_ABSENCE_OUT_DIR / "cd-hit" / "{sample}.faa",
    shell:
        r"""
        awk 'BEGIN {{OFS=""}} \
            !/^>/ {{print; next}} \
            /^>/ {{
                header = substr($0, 2);
                gsub(/[^a-zA-Z0-9]/, "_", header);
                sub(/_+$/, "", header);
                print ">{wildcards.sample}_" header
            }}' {input} > {output}
        """


rule combine_faa_files:
    localrule: True
    input: 
        lambda wc: expand(
            rules.cdhit_protein_name_corrector.output,
            sample = get_sample_names(wc)
        )
    output: GENE_PRESENCE_ABSENCE_OUT_DIR / "cd-hit" / "combined_proteins.faa",
    shell:
        r"""
        cat {input} > {output}
        """


rule cdhit_protein_positions:
    input: 
        lambda wc: expand(
            rules.prokka_runner.output.gbk,
            sample = get_sample_names(wc)
        )
    output: GENE_PRESENCE_ABSENCE_OUT_DIR / "cd-hit" / "protein_positions.csv",
    benchmark: BENCHMARKS_DIR / "cdhit_protein_positions.py.tsv"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "cdhit_protein_positions.py"


rule cdhit_runner:
    input: rules.combine_faa_files.output
    output: 
        faa = GENE_PRESENCE_ABSENCE_OUT_DIR / "cd-hit" / "cdhit_output.faa",
        clstr = GENE_PRESENCE_ABSENCE_OUT_DIR / "cd-hit" / "cdhit_output.faa.clstr",
    log: GENE_PRESENCE_ABSENCE_LOGS_DIR / "cdhit_runner.log"
    benchmark: BENCHMARKS_DIR / "cdhit.tsv"
    params:
        seq_identity_threshold = 0.7,
        length_difference_cutoff = 0.0, # (%)
        aln_cov_longer_seq = 0.0, # alingment coverage for the longer sequence
        aln_cov_control_longer_seq = 99_999_999, # alignment coverage control for the longer sequence
        aln_cov_shorter_seq = 0.0, # alignment coverage for the shorter sequence
        aln_cov_control_shorter_seq = 99_999_999, # alignment coverage control for the shorter sequence
        unlimited_memory = 0, # memory limit (in MB) for the program; 0 for unlimited;
    threads: workflow.cores
    conda: ENVS_DIR.format("cd-hit")
    container: CONTAINERS.format("cd-hit:1.0.0")
    shell:
        r"""
        cd-hit \
            -i {input} \
            -o {output.faa} \
            -T {threads} \
            -M {params.unlimited_memory} \
            -c {params.seq_identity_threshold} \
            -s {params.length_difference_cutoff} \
            -aL {params.aln_cov_longer_seq} \
            -AL {params.aln_cov_control_longer_seq} \
            -aS {params.aln_cov_shorter_seq} \
            -AS {params.aln_cov_control_shorter_seq} \
            -d 0 \
            > {log} 2>&1
        """


rule binary_gpa_cdhit:
    input:
        rules.cdhit_runner.output.clstr,
    output: GENE_PRESENCE_ABSENCE_OUT_DIR / "binary_gpa_cdhit.tsv"
    benchmark: BENCHMARKS_DIR / "binary_gpa_cdhit.tsv"
    log: GENE_PRESENCE_ABSENCE_LOGS_DIR / "binary_gpa_cdhit.log"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "binary_gpa_cdhit.py"


# -----------------------
# Feature Tables
# -----------------------

rule binary_gpa:
    localrule: True
    input:
        branch(
            condition=lookup(dpath="gpa_method", within=config, default="cd-hit"),
            cases={
                "cd-hit": rules.binary_gpa_cdhit.output,
                "panaroo": rules.binary_gpa_panaroo.output,
            }
        )
    output: GENE_PRESENCE_ABSENCE_OUT_DIR / "binary_gpa.tsv"
    log: GENE_PRESENCE_ABSENCE_LOGS_DIR / "binary_gpa.log"
    shell:
        r"""
        ln -srv {input} {output} >> {log} 2>&1
        """

rule gene_presence_absence:
    localrule: True
    input:
        rules.binary_gpa.output
    output: touch(OUT_DIR / "flags" / "gene_presence_absence.done")
