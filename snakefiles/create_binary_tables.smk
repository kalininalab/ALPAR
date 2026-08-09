from pathlib import Path


# -----------------------
# Create link to files in singleton directory
# -----------------------

checkpoint rename_files:
    input: IN_DIR
    output:
        store = directory(TEMP_DIR / "data_checksum"),
        mapping = OUT_DIR / "all_files.tsv",
    benchmark: BENCHMARKS_DIR / "rename_files.tsv"
    log: LOGS_DIR / "rename_files.log"
    threads: 1
    shell:
        r"""
        mkdir -p {output.store}
        echo -e "checksum\tfilepath" > {output.mapping}
        
        find -L {input} -type f \( -name "*.fna" -o -name "*.fasta" -o -name "*.faa" \) -print0 | \
        while read -r -d '' file; do
        
            resolved_path=$(readlink -f "$file")
            checksum=$(shasum "$file" -a 1 | cut -d ' ' -f 1)
            echo "Resolving $file to $resolved_path with checksum $checksum" >> {log}

            if [ ! -e "{output.store}/$checksum" ]; then
                ln -srv "$resolved_path" "{output.store}/$checksum" >> {log} 2>&1
            else
                echo "Checksum file $checksum already exists, skipping" >> {log}
            fi

            echo -e "$checksum\t$file" >> {output.mapping}
        
        done
        """

def get_sample_names(wildcards) -> list[str]:
    rename_checkpoint = checkpoints.rename_files.get(**wildcards)
    rename_folder = Path(rename_checkpoint.output.store)
    sample_names = [
        f.name 
        for f in rename_folder.iterdir()
        if f.name != ".snakemake_timestamp"
    ]
    return sample_names


# -----------------------
# Phenotype
# -----------------------

rule phenotype_dataframe_creator:
    input: rules.rename_files.output.mapping
    output: OUT_DIR / "phenotype_table.tsv"
    benchmark: TEMP_DIR / "benchmarks" / "phenotype_dataframe_creator.tsv"
    log: LOGS_DIR / "phenotype_dataframe_creator.log"
    conda: ENVS_DIR.format("python313")
    params:
        resistance_status_mapping = RESISTANCE_STATUS_MAPPING,
        antibiotics = ANTIBIOTICS,
    threads: 1
    script:
        SCRIPTS_DIR / "phenotype_dataframe_creator.py"


# -----------------------
# Genome annotation
# -----------------------

rule cd_hit_create_db:
    input: FASTA_FILE
    output: TEMP_DIR / GENUS / GENUS
    log: LOGS_DIR / "cd_hit_create_db.log"
    benchmark: BENCHMARKS_DIR / "cdhit_create_db.tsv"
    conda: ENVS_DIR.format("cd-hit")
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
    output: touch(TEMP_DIR / "flags" / "makeblastdb.done")
    log: LOGS_DIR / "makeblastdb.log"
    benchmark: BENCHMARKS_DIR / "makeblastdb.tsv"
    conda: ENVS_DIR.format("makeblastdb")
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
    output: touch(TEMP_DIR / "flags" / "prokka_listdb.done"),
    log: LOGS_DIR / "prokka_listdb.log"
    benchmark: BENCHMARKS_DIR / "prokka_listdb.tsv"
    conda: ENVS_DIR.format("prokka")
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
        gff = TEMP_DIR / "prokka" / "{sample}" / "{sample}.gff",
        faa = TEMP_DIR / "prokka" / "{sample}" / "{sample}.faa",
        gbk = TEMP_DIR / "prokka" / "{sample}" / "{sample}.gbk",
    log: LOGS_DIR / "prokka_runner" / "{sample}.log"
    benchmark: BENCHMARKS_DIR / "prokka_{sample}.tsv"
    params:
        genus = GENUS,
        outdir = subpath(output.gff, parent=True),
    threads: 1
    resources:
        mem_mb = 600
    conda: ENVS_DIR.format("prokka")
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
            sample = get_sample_names(wc)
        ),
    output:
        gpa = TEMP_DIR / "panaroo" / "gene_presence_absence.csv",
        gene_data = TEMP_DIR / "panaroo" / "gene_data.csv",
    log: LOGS_DIR / "panaroo_runner.log"
    benchmark: BENCHMARKS_DIR / "panaroo.tsv"
    params:
        outdir = subpath(output.gpa, parent=True),
        seq_identity_threshold = 0.8,
        seq_len_diff_cutoff = 0.8,
    conda: ENVS_DIR.format("panaroo")
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
    output: TEMP_DIR / "binary_gpa_panaroo.tsv"
    benchmark: BENCHMARKS_DIR / "binary_gpa_panaroo.py.tsv"
    conda: ENVS_DIR.format("python313")
    threads: 1
    script:
        SCRIPTS_DIR / "binary_gpa_panaroo.py"


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
            sample = get_sample_names(wc)
        )
    output: TEMP_DIR / "cd-hit" / "combined_proteins.faa",
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
    output: OUT_DIR / "cd-hit" / "protein_positions.csv",
    benchmark: BENCHMARKS_DIR / "cdhit_protein_positions.py.tsv"
    conda: ENVS_DIR.format("python313")
    threads: 1
    script:
        SCRIPTS_DIR / "cdhit_protein_positions.py"


rule cdhit_runner:
    input: rules.combine_faa_files.output
    output: 
        faa = TEMP_DIR / "cd-hit" / "cdhit_output.faa",
        clstr = TEMP_DIR / "cd-hit" / "cdhit_output.faa.clstr",
    log: LOGS_DIR / "cdhit_runner.log"
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
    output: TEMP_DIR / "binary_gpa_cdhit.tsv"
    benchmark: BENCHMARKS_DIR / "binary_gpa_cdhit.tsv"
    log: LOGS_DIR / "binary_gpa_cdhit.log"
    conda: ENVS_DIR.format("python313")
    threads: 1
    script:
        SCRIPTS_DIR / "binary_gpa_cdhit.py"

# -----------------------
# Mutation Analysis
# -----------------------

rule snippy_runner:
    input:
        sample_store = rules.rename_files.output.store,
        sample = Path(rules.rename_files.output.store) / "{sample}",
        reference = GBFF_FILE,
    output:
        vcf = TEMP_DIR / "snippy" / "{sample}" / "snps.vcf",
        tab = TEMP_DIR / "snippy" / "{sample}" / "snps.tab",
    log: LOGS_DIR / "snippy_runner" / "{sample}.log"
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
    output: TEMP_DIR / "mutations_annotations.tsv"
    benchmark: BENCHMARKS_DIR / "annotation_file_from_snippy.tsv"
    log: LOGS_DIR / "annotation_file_from_snippy.log"
    conda: ENVS_DIR.format("python313")
    threads: MAX_PYTHON_THREADS
    script:
        SCRIPTS_DIR / "annotation_file_from_snippy.py"


# -----------------------
# Feature Tables
# -----------------------

rule binary_gpa:
    input:
        branch(
            condition=lookup(dpath="gpa_method", within=config, default="cd-hit"),
            cases={
                "cd-hit": rules.binary_gpa_cdhit.output,
                "panaroo": rules.binary_gpa_panaroo.output,
            }
        )
    output: TEMP_DIR / "binary_gpa.tsv"
    log: LOGS_DIR / "binary_gpa.log"
    shell:
        r"""
        ln -srv {input} {output} >> {log} 2>&1
        """

rule binary_mutation_table:
    input:
        lambda wc: expand(
            rules.snippy_runner.output.vcf,
            sample = get_sample_names(wc)
        )
    output: TEMP_DIR / "binary_mutation_table.tsv"
    benchmark: BENCHMARKS_DIR / "binary_mutation_table.tsv"
    log: LOGS_DIR / "binary_mutation_table.log"
    conda: ENVS_DIR.format("python313")
    threads: 1
    script:
        SCRIPTS_DIR / "binary_mutation_table.py"

# -----------------------
# Snakefile Target
# -----------------------

rule merge_binary_features:
    input:
        rules.binary_mutation_table.output,
        rules.binary_gpa.output,
    output: OUT_DIR / "merged_binary_table.tsv"
    threads: 1
    shell:
        r"""
        cat {input} > {output}
        """

rule create_binary_tables:
    input:
        rules.merge_binary_features.output,
        rules.annotation_file_from_snippy.output,
        rules.cdhit_protein_positions.output,
    output: touch(TEMP_DIR / "flags" / "create_binary_tables.done")
