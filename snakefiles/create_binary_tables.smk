# snakemake --use-conda --jobs 8 --resources mem_gb=4 --benchmark-extended --rerun-triggers mtime
from pathlib import Path

configfile: "config/config.yaml"

# -----------------------
# Global Variables
# -----------------------

GENUS = config.get("genus")

# -----------------------
# Directories and Files
# -----------------------

IN_DIR = Path(config.get("input_dir", "data"))
OUT_DIR = Path(config.get("output_dir", "out"))
TEMP_DIR = Path(config.get("temp_dir", "temp"))
GBFF_FILE = Path(config.get("gbff_file"))
FASTA_FILE = Path(config.get("fasta_file"))
CHECKSUM_DIR = TEMP_DIR / "data_checksum"


# -----------------------
# Create link to files in singleton directory
# -----------------------

checkpoint rename_files:
    input: IN_DIR
    output:
        store = directory(CHECKSUM_DIR),
        mapping = OUT_DIR / "all_files.tsv",
    shell:
        r"""
        mkdir -p {output.store}
        echo -e "checksum\tfilepath" > {output.mapping}
        
        find {input} -type f -name "*.fna" -print0 | while read -r -d '' file; do

            checksum=$(shasum "$file" -a 1 | cut -d ' ' -f 1)
            ln -sr "$file" "{output.store}/$checksum"
            echo -e "$checksum\t$file" >> {output.mapping}

        done
        """


def get_sample_names() -> list[str]:
    rename_checkpoint = checkpoints.rename_files.get()
    rename_folder = Path(rename_checkpoint.output.store)
    sample_names = [
        f.name 
        for f in rename_folder.iterdir()
        if f.name != ".snakemake_timestamp"
    ]
    return sample_names


# -----------------------
# Genome annotation
# -----------------------

rule cd_hit_create_db:
    input: FASTA_FILE
    output: TEMP_DIR / GENUS / GENUS
    log: TEMP_DIR / "prokka_db.log"
    conda: "envs/cd-hit.yaml"
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
    output: touch(TEMP_DIR / "flags" / "makeblastdb")
    log: TEMP_DIR / "makeblastdb.log"
    conda: "envs/makeblastdb.yaml"
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
    output: touch(TEMP_DIR / "flags" / "prokka_listdb"),
    params:
        prokka_listdb = TEMP_DIR / "db_path.txt"
    conda: "envs/alpar.yaml"
    shell:
        r"""
        DB_DIR=$(dirname {input.db_dir})
        PROKKA_DB_DIR="$(dirname $(dirname $(which prokka)))/db/genus"
        cp -a "$DB_DIR/." $PROKKA_DB_DIR
        prokka --listdb > {params.prokka_listdb} 2>&1
        """


#TODO Implement branching logic if no reference is given
rule prokka_runner:
    input:
        rules.prokka_listdb.output,
        sample = CHECKSUM_DIR / "{sample}",
        reference = GBFF_FILE,
    output: 
        gff = OUT_DIR / "prokka" / "{sample}" / "{sample}.gff",
        faa = OUT_DIR / "prokka" / "{sample}" / "{sample}.faa",
        gbk = OUT_DIR / "prokka" / "{sample}" / "{sample}.gbk",
    log: TEMP_DIR / "prokka" / "{sample}.log"
    benchmark: TEMP_DIR / "benchmarks" / "prokka_{sample}.tsv"
    params:
        genus = GENUS,
        outdir = subpath(output.gff, parent=True),
    threads: 4
    resources:
        mem_mb = 500
    conda: "envs/prokka.yaml"
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


# -----------------------
# Gene Presence Absence: Panaroo
# -----------------------

rule panaroo_runner:
    input: 
        lambda wc: expand(
            rules.prokka_runner.output.gff,
            sample = get_sample_names()
        ),
    output: OUT_DIR / "panaroo" / "gene_presence_absence.csv",
    log: TEMP_DIR / "panaroo.log"
    benchmark: TEMP_DIR / "benchmarks" / "panaroo.tsv"
    params:
        outdir = subpath(output[0], parent=True),
    conda: "envs/panaroo.yaml"
    threads: 30
    shell:
        r"""
        panaroo \
            --input {input} \
            --out_dir {params.outdir} \
            --clean-mode strict \
            --threads {threads} \
            >> {log} 2>&1
        """


rule binary_gpa_panaroo:
    input: rules.panaroo_runner.output,
    output: OUT_DIR / "binary_gpa_panaroo.tsv"
    conda: "envs/python312.yaml"
    threads: 1
    script:
        "scripts/binary_gpa_panaroo.py"


# -----------------------
# Gene Presence Absence: CD-HIT
# -----------------------

rule cdhit_protein_name_corrector:
    input: rules.prokka_runner.output.faa,
    output: TEMP_DIR / "cd-hit" / "{sample}.faa",
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
    input: 
        lambda wc: expand(
            rules.cdhit_protein_name_corrector.output,
            sample = get_sample_names()
        )
    output: OUT_DIR / "cd-hit" / "combined_proteins.faa",
    shell:
        r"""
        cat {input} > {output}
        """


rule cdhit_protein_positions:
    input: 
        lambda wc: expand(
            rules.prokka_runner.output.gbk,
            sample = get_sample_names()
        )
    output: OUT_DIR / "cd-hit" / "protein_positions.csv",
    conda: "envs/python312.yaml"
    threads: 4
    script:
        "scripts/cdhit_protein_positions.py"


rule cdhit_runner:
    input: rules.combine_faa_files.output
    output: 
        faa = OUT_DIR / "cd-hit" / "cdhit_output.faa",
        clstr = OUT_DIR / "cd-hit" / "cdhit_output.faa.clstr",
    benchmark: TEMP_DIR / "benchmarks" / "cdhit.tsv"
    params:
        seq_identity_threshold = 0.95,
        length_difference_cutoff = 0.0, # (%)
        aln_cov_longer_seq = 0.0, # alingment coverage for the longer sequence
        aln_cov_control_longer_seq = 99_999_999, # alignment coverage control for the longer sequence
        aln_cov_shorter_seq = 0.0, # alignment coverage for the shorter sequence
        aln_cov_control_shorter_seq = 99_999_999, # alignment coverage control for the shorter sequence
        unlimited_memory = 0, # memory limit (in MB) for the program; 0 for unlimitted;
    threads: workflow.cores
    conda: "envs/cd-hit.yaml"
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
            > /dev/null
        """


rule binary_gpa_cdhit:
    input: 
        clstr = rules.cdhit_runner.output.clstr,
        positions = rules.cdhit_protein_positions.output,
    output: OUT_DIR / "binary_gpa_cdhit.tsv"
    conda: "envs/python312.yaml"
    script:
        "scripts/binary_gpa_cdhit.py"


# -----------------------
# Binary Gene Presence Absence
# -----------------------

rule binary_gpa:
    input:
        branch(
            lookup(dpath="gpa_method", within=config, default="cd-hit"),
            then=rules.binary_gpa_cdhit.output,
            otherwise=rules.binary_gpa_panaroo.output,
        )
    output: OUT_DIR / "binary_gpa.tsv"
    shell:
        r"""
        mv {input} {output}
        """


# -----------------------
# Mutation Analysis
# -----------------------

rule snippy_runner:
    input:
        sample = CHECKSUM_DIR / "{sample}",
        reference = GBFF_FILE,
    output:
        vcf = OUT_DIR / "snippy" / "{sample}" / "snps.vcf",
        tab = OUT_DIR / "snippy" / "{sample}" / "snps.tab",
    log: TEMP_DIR / "snippy" / "{sample}.log"
    benchmark: TEMP_DIR / "benchmarks" / "snippy_{sample}.tsv"
    params:
        out_dir = subpath(output.vcf, parent=True)
    conda: "envs/snippy.yaml"
    threads: workflow.cores
    resources:
        mem_gb = workflow.global_resources.get("mem_gb", 4)
    shell:
        r"""
        input_file=$(readlink -f {input.sample})
        snippy \
            --ctgs $input_file \
            --outdir {params.out_dir} \
            --reference {input.reference} \
            --cpus {threads} \
            --ram {resources.mem_gb} \
            --force \
            >> {log} 2>&1
        """


rule binary_mutation_table:
    input: 
        lambda wc: expand(
            rules.snippy_runner.output.vcf,
            sample = get_sample_names()
        )
    output: OUT_DIR / "binary_mutation_table.tsv"
    conda: "envs/python312.yaml"
    threads: 4
    script:
        "scripts/binary_mutation_table.py"


# -----------------------
# Featre Importance Analysis
# -----------------------

rule annotation_file_from_snippy:
    input:
        lambda wc: expand(
            rules.snippy_runner.output.tab,
            sample = get_sample_names()
        )
    output: OUT_DIR / "mutations_annotations.tsv"
    conda: "envs/python312.yaml"
    threads: 4
    script:
        "scripts/annotation_file_from_snippy.py"


# -----------------------
# Binary Tables
# -----------------------

rule phenotype_dataframe_creator:
    input: rules.rename_files.output.mapping
    output: OUT_DIR / "phenotype_table.tsv"
    conda: "envs/python312.yaml"
    threads: 1
    script:
        "scripts/phenotype_dataframe_creator.py"


rule merge_binary_tables:
    input:
        rules.binary_mutation_table.output,
        rules.binary_gpa.output,
    output: OUT_DIR / "merged_binary_table.tsv"
    conda: "envs/python312.yaml"
    threads: 1
    script:
        "scripts/merge_binary_tables.py"


rule create_binary_tables:
    input:
        rules.merge_binary_tables.output,
        rules.phenotype_dataframe_creator.output,
        rules.annotation_file_from_snippy.output,
    default_target: True

