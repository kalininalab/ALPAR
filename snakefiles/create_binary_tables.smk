import functools
from pathlib import Path

workflow_dir = Path(workflow.basedir)
configfile: workflow_dir / "snakefiles" / "config" / "config.yaml"

# -----------------------
# Global Variables
# -----------------------

MAX_PYTHON_THREADS = min(workflow.cores, 32)
JOB_BATCH_SIZE = 100

GENUS = config.get("genus")
RESISTANCE_STATUS_MAPPING = {
    'Resistant': 1,
    'Susceptible': 0,
}

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
        
        find {input} -type f \( -name "*.fna" -o -name "*.fasta" -o -name "*.faa" \) -print0 | \
        while read -r -d '' file; do
        
            checksum=$(shasum "$file" -a 1 | cut -d ' ' -f 1)
            ln -sr "$file" "{output.store}/$checksum"
            echo -e "$checksum\t$file" >> {output.mapping}

        done
        """

@functools.lru_cache
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

# TODO: Should strains not present in an antibiotic, be assumed susceptible?
rule phenotype_dataframe_creator:
    input: rules.rename_files.output.mapping
    output: OUT_DIR / "phenotype_table.tsv"
    benchmark: TEMP_DIR / "benchmarks" / "phenotype_dataframe_creator.py.tsv"
    conda: "envs/python312.yaml"
    params:
        resistance_status_mapping = RESISTANCE_STATUS_MAPPING
    threads: 1
    script:
        "scripts/phenotype_dataframe_creator.py"


# -----------------------
# Genome annotation
# -----------------------

rule cd_hit_create_db:
    input: FASTA_FILE
    output: TEMP_DIR / GENUS / GENUS
    log: TEMP_DIR / "logs" / "cd_hit_create_db.log"
    benchmark: TEMP_DIR / "benchmarks" / "cdhit_create_db.tsv"
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
    log: TEMP_DIR / "logs" / "makeblastdb.log"
    benchmark: TEMP_DIR / "benchmarks" / "makeblastdb.tsv"
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
    conda: "envs/prokka.yaml"
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
        sample_store = rules.rename_files.output.store,
        sample = Path(rules.rename_files.output.store) / "{sample}",
        reference = GBFF_FILE,
    output:
        gff = OUT_DIR / "prokka" / "{sample}" / "{sample}.gff",
        faa = OUT_DIR / "prokka" / "{sample}" / "{sample}.faa",
        gbk = OUT_DIR / "prokka" / "{sample}" / "{sample}.gbk",
    log: TEMP_DIR / "logs" / "prokka_runner" / "{sample}.log"
    benchmark: TEMP_DIR / "benchmarks" / "prokka_{sample}.tsv"
    params:
        genus = GENUS,
        outdir = subpath(output.gff, parent=True),
    threads: 1
    resources:
        mem_mb = 600
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
            sample = get_sample_names(wc)
        ),
    output:
        gpa = OUT_DIR / "panaroo" / "gene_presence_absence.csv",
        gene_data = OUT_DIR / "panaroo" / "gene_data.csv",
    log: TEMP_DIR / "logs" / "panaroo_runner.log"
    benchmark: TEMP_DIR / "benchmarks" / "panaroo.tsv"
    params:
        outdir = subpath(output.gpa, parent=True),
        seq_identity_threshold = 0.8,
        seq_len_diff_cutoff = 0.8,
    conda: "envs/panaroo.yaml"
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
    output: OUT_DIR / "binary_gpa_panaroo.tsv"
    benchmark: TEMP_DIR / "benchmarks" / "binary_gpa_panaroo.py.tsv"
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
            sample = get_sample_names(wc)
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
            sample = get_sample_names(wc)
        )
    output: OUT_DIR / "cd-hit" / "protein_positions.csv",
    benchmark: TEMP_DIR / "benchmarks" / "cdhit_protein_positions.py.tsv"
    conda: "envs/python312.yaml"
    threads: MAX_PYTHON_THREADS
    script:
        "scripts/cdhit_protein_positions.py"


rule cdhit_runner:
    input: rules.combine_faa_files.output
    output: 
        faa = OUT_DIR / "cd-hit" / "cdhit_output.faa",
        clstr = OUT_DIR / "cd-hit" / "cdhit_output.faa.clstr",
    log: TEMP_DIR / "logs" / "cdhit_runner.log"
    benchmark: TEMP_DIR / "benchmarks" / "cdhit.tsv"
    params:
        seq_identity_threshold = 0.7,
        length_difference_cutoff = 0.0, # (%)
        aln_cov_longer_seq = 0.0, # alingment coverage for the longer sequence
        aln_cov_control_longer_seq = 99_999_999, # alignment coverage control for the longer sequence
        aln_cov_shorter_seq = 0.0, # alignment coverage for the shorter sequence
        aln_cov_control_shorter_seq = 99_999_999, # alignment coverage control for the shorter sequence
        unlimited_memory = 0, # memory limit (in MB) for the program; 0 for unlimited;
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
            > {log} 2>&1
        """


rule binary_gpa_cdhit:
    input: 
        rules.cdhit_runner.output.clstr,
    output: OUT_DIR / "binary_gpa_cdhit.tsv"
    benchmark: TEMP_DIR / "benchmarks" / "binary_gpa_cdhit.py.tsv"
    conda: "envs/python312.yaml"
    threads: 1
    script:
        "scripts/binary_gpa_cdhit.py"


# ------------------------
# Split Clusters
# ------------------------

checkpoint split_cluster_fasta:
    input:
        cdhit_clstr = rules.cdhit_runner.output.clstr,
        combined_proteins = rules.combine_faa_files.output[0],
    output: directory(OUT_DIR / "cluster_sequences"),
    log: TEMP_DIR / "logs" / "split_cluster_fasta.log"
    benchmark: TEMP_DIR / "benchmarks" / "split_cluster_fasta.py.tsv"
    params:
        file_ext = ".fasta"
    conda: "envs/python312.yaml"
    threads: 8
    script:
        "scripts/split_cluster_fasta.py"


# -----------------------
# Multiple Sequence Alignment: MAFFT
# -----------------------

@functools.lru_cache
def get_cluster_files(wildcards) -> list[Path]:
    cluster_checkpoint = checkpoints.split_cluster_fasta.get(**wildcards)
    cluster_folder = Path(cluster_checkpoint.output[0])
    file_ext = rules.split_cluster_fasta.params.file_ext
    return sorted(cluster_folder.glob(f"*{file_ext}"))

def batched_clusters(wildcards) -> list[Path]:
    cluster_files = get_cluster_files(wildcards)
    batch_num = int(wildcards.batch_num)
    start = batch_num * JOB_BATCH_SIZE
    end = min(start + JOB_BATCH_SIZE, len(cluster_files))
    return cluster_files[start:end]

rule batch_align_clusters:
    input:
        cluster_store = rules.split_cluster_fasta.output,
        batched_clusters = batched_clusters
    output: directory(TEMP_DIR / "batch_align_clusters" / "batch_{batch_num}")
    log: TEMP_DIR / "logs" / "batch_align_clusters" / "batch_{batch_num}.log"
    benchmark: TEMP_DIR / "benchmarks" / "batch_align_clusters_batch_{batch_num}.tsv"
    conda: "envs/mafft.yaml"
    threads: 1
    shell:
        r"""
        mkdir -p {output}

        for i in {input.batched_clusters}; do
            OUT_FILE="{output}/$(basename $i)"
            FASTA_COUNT=$(grep -c "^>" $i)

            if [ $FASTA_COUNT -eq 1 ]; then
                echo ">>$i with $FASTA_COUNT sequence, skipping alignment" >> {log}
                ln -sr $i $OUT_FILE
            else
                echo ">>$i with $FASTA_COUNT sequences" >> {log}
                mafft --auto --thread {threads} $i > $OUT_FILE 2>> {log}
            fi
        done
        """

rule gather_align_clusters:
    input:
        lambda wc: expand(
            rules.batch_align_clusters.output,
            batch_num = range((len(get_cluster_files(wc)) - 1) // JOB_BATCH_SIZE + 1)
        )
    output: directory(OUT_DIR / "cluster_alignments")
    shell:
        r"""
        mkdir -p {output}
        for batch_dir in {input}; do
            for aln_file in $batch_dir/*.fasta; do
                [ -f "$aln_file" ] || continue
                ln -sr $aln_file {output}/$(basename $aln_file)
            done
        done
        """


# -----------------------
# Cluster True Variants: MAFFT
# -----------------------


checkpoint split_cluster_by_phenotype:
    input:
        clstr_store = rules.split_cluster_fasta.output,
        clstr_sequences = get_cluster_files,
        phenotypes = rules.phenotype_dataframe_creator.output[0],
    output:
        store = directory(TEMP_DIR / "clusters_by_phenotype")
    log: TEMP_DIR / "logs" / "split_cluster_by_phenotype.log"
    benchmark: TEMP_DIR / "benchmarks" / "split_cluster_by_phenotype.py.tsv"
    params:
        resistance_status_mapping = RESISTANCE_STATUS_MAPPING
    conda: "envs/python312.yaml"
    threads: MAX_PYTHON_THREADS
    script:
        "scripts/split_cluster_by_phenotype.py"

def get_all_clustered_antibiotics(wildcards) -> tuple[Path]:
    cluster_checkpoint = checkpoints.split_cluster_by_phenotype.get(**wildcards)
    cluster_folder = Path(cluster_checkpoint.output.store)
    antibiotics = tuple(
        antibiotic_folder.name
        for antibiotic_folder in cluster_folder.iterdir()
        if antibiotic_folder.is_dir()
    )
    return antibiotics

@functools.lru_cache
def get_cluster_by_phenotype_files(resistance_status: str, wildcards) -> list[Path]:
    cluster_checkpoint = checkpoints.split_cluster_by_phenotype.get(**wildcards)
    folder = Path(cluster_checkpoint.output.store) / wildcards.antibiotic / resistance_status
    file_ext = rules.split_cluster_fasta.params.file_ext
    return sorted(folder.glob(f"*{file_ext}"))

def input_batched_cluster_by_phenotype(resistance_status: str, wildcards) -> list[Path]:
    batch_files = get_cluster_by_phenotype_files(resistance_status, wildcards)
    batch_num = int(wildcards.batch_num)
    start = batch_num * JOB_BATCH_SIZE
    end = min(start + JOB_BATCH_SIZE, len(batch_files))
    return batch_files[start:end]

def output_batched_cluster_by_phenotype(suffix: str, wildcards, input) -> list[Path]:
    batch_num = int(wildcards.batch_num)
    clstr_out_dir = TEMP_DIR / "cluster_antibiotic_alignment" / wildcards.antibiotic / f"batch_{batch_num}"
    batch_clstr = (Path(f).stem for f in input.clstr_susceptible)
    return [clstr_out_dir / f"{clstr_name}_{suffix}" for clstr_name in batch_clstr]

rule align_clusters_by_phenotype_and_map_variants:
    input:
        clstr_store = rules.split_cluster_by_phenotype.output,
        clstr_susceptible = lambda wc: input_batched_cluster_by_phenotype("Susceptible", wc),
        clstr_resistant = lambda wc: input_batched_cluster_by_phenotype("Resistant", wc),
    output: directory(TEMP_DIR / "cluster_antibiotic_alignment" / "{antibiotic}" / "batch_{batch_num}")
    log: TEMP_DIR / "logs" / "align_clusters_by_phenotype_and_map_variants" / "{antibiotic}" / "batch_{batch_num}.log"
    benchmark: TEMP_DIR / "benchmarks" / "align_clusters_by_phenotype_and_map_variants_{antibiotic}_batch_{batch_num}.tsv"
    params:
        aln_susceptible = lambda wc, input: output_batched_cluster_by_phenotype("only_susceptible.fasta", wc, input),
        aln_resistant = lambda wc, input: output_batched_cluster_by_phenotype("all.fasta", wc, input),
        map_files = lambda wc, input: output_batched_cluster_by_phenotype("all.fasta.map", wc, input),
    conda: "envs/mafft.yaml"
    threads: 1
    shell:
        r"""
        # Create output directory if it doesn't exist
        mkdir -p {output}
        mkdir -p $(dirname {log})

        susceptible_array=({input.clstr_susceptible})
        resistant_array=({input.clstr_resistant})
        base_aln_array=({params.aln_susceptible})
        all_aln_array=({params.aln_resistant})
        all_map_array=({params.map_files})

        for i in "${{!susceptible_array[@]}}"; do
            SUSCEPTIBLE_FILE="${{susceptible_array[$i]}}"
            RESISTANT_FILE="${{resistant_array[$i]}}"
            BASE_ALN="${{base_aln_array[$i]}}"
            ALL_ALN="${{all_aln_array[$i]}}"
            ALL_MAP="${{all_map_array[$i]}}"

            CLSTR_NAME="$(basename ${{SUSCEPTIBLE_FILE%.faa}})"

            echo ">>Aligning Susceptibles: $CLSTR_NAME" >> {log}
            mafft \
                --auto \
                --thread {threads} \
                $SUSCEPTIBLE_FILE \
                > $BASE_ALN \
                2>> {log}

            echo ">>Aligning Resistant: $CLSTR_NAME" >> {log}
            mafft \
                --auto \
                --thread {threads} \
                --add $RESISTANT_FILE \
                --mapout \
                --keeplength \
                $BASE_ALN \
                > $ALL_ALN \
                2>> {log}
            mv $RESISTANT_FILE.map $ALL_MAP
        done
        """

rule gather_align_and_map_clusters_by_batch:
    input:
        lambda wc: expand(
            rules.align_clusters_by_phenotype_and_map_variants.output,
            antibiotic = wc.antibiotic,
            batch_num = range((len(get_cluster_by_phenotype_files("Susceptible", wc)) - 1) // JOB_BATCH_SIZE + 1)
        )
    output:
        susceptible_alignments = directory(OUT_DIR / "susceptible_alignments" / "{antibiotic}"),
        all_alignments = directory(OUT_DIR / "all_alignments" / "{antibiotic}"),
        map_files = directory(OUT_DIR / "map_files" / "{antibiotic}"),
    threads: 1
    shell:
        r"""
        mkdir -p {output.susceptible_alignments} {output.all_alignments} {output.map_files}

        for batch_dir in {input}; do
            for file in "$batch_dir"/*; do
                [ -f "$file" ] || continue
                case "$file" in
                    *_only_susceptible.fasta)
                        ln -sr "$file" {output.susceptible_alignments}/$(basename "$file") ;;
                    *_all.fasta)
                        ln -sr "$file" {output.all_alignments}/$(basename "$file") ;;
                    *.map)
                        ln -sr "$file" {output.map_files}/$(basename "$file") ;;
                esac
            done
        done
        """


# -----------------------
# Cluster Variants: PanPA
# -----------------------


rule panpa_build_index_by_phenotype:
    input: rules.gather_align_and_map_clusters_by_batch.output.susceptible_alignments
    output: OUT_DIR / "cluster_panpa" / "{antibiotic}" / "index.pickle"
    log: TEMP_DIR / "logs" / "cluster_panpa_build_index_by_phenotype" / "{antibiotic}.log"
    benchmark: TEMP_DIR / "benchmarks" / "cluster_panpa_build_index_by_phenotype_{antibiotic}.tsv"
    params:
        kmer_size = 10,
        window_size = 15,
        seed_limit = 0,
    conda: "envs/panpa.yaml"
    threads: 1
    shell:
        r"""
        PanPA \
            --log_file {log} \
            build_index \
            --in_dir {input} \
            --out_index {output} \
            --seeding_alg wk_min \
            --kmer_size {params.kmer_size} \
            --window {params.window_size} \
            --seed_limit {params.seed_limit}
        """

rule panpa_build_gfa_by_phenotype:
    input: rules.gather_align_and_map_clusters_by_batch.output.susceptible_alignments
    output: directory(OUT_DIR / "cluster_panpa" / "{antibiotic}" / "gfa")
    log: TEMP_DIR / "logs" / "cluster_panpa_build_gfa" / "{antibiotic}.log"
    benchmark: TEMP_DIR / "benchmarks" / "cluster_panpa_build_gfa_{antibiotic}.tsv"
    conda: "envs/panpa.yaml"
    threads: workflow.cores
    shell:
        r"""
        PanPA \
            --log_file {log} \
            build_gfa \
            --in_dir {input} \
            --out_dir {output} \
            --cores {threads}
        """

rule panpa_align_cluster_single_target:
    input:
        clstr_resistant = lambda wc: get_cluster_by_phenotype_files("Resistant", wc),
        gfa_folder = rules.panpa_build_gfa_by_phenotype.output[0],
    output: directory(TEMP_DIR / "cluster_panpa_alignments" / "{antibiotic}")
    log: TEMP_DIR / "logs" / "panpa_align_cluster_single_target" / "{antibiotic}.log"
    benchmark: TEMP_DIR / "benchmarks" / "panpa_align_cluster_single_target_{antibiotic}.log"
    conda: "envs/panpa.yaml"
    threads: workflow.cores
    shell:
        r"""
        # Create output directory if it doesn't exist
        mkdir -p {output}
        mkdir -p $(dirname {log})

        for i in {input.clstr_resistant}; do
            CLSTR_NAME="$(basename ${{i%.faa}})"
            TEMP_RENAME_INPUT={output}/$CLSTR_NAME.fasta
            TEMP_LOG={log}.temp

            ln -sr $i $TEMP_RENAME_INPUT

            PanPA \
                --log_file $TEMP_LOG \
                align_single \
                --gfa_files "{input.gfa_folder}/${{CLSTR_NAME}}_only_susceptible.fasta.gfa" \
                --seqs $TEMP_RENAME_INPUT \
                --cores {threads} \
                --out_gaf {output}/$CLSTR_NAME.gaf
            
            cat $TEMP_LOG >> {log}
            rm $TEMP_LOG
            rm $TEMP_RENAME_INPUT
        done
        """

rule gather_panpa_all_clusters:
    input:
        lambda wc: expand(
            rules.panpa_build_index_by_phenotype.output,
            antibiotic = get_all_clustered_antibiotics(wc)
        ),
        lambda wc: expand(
            rules.panpa_build_gfa_by_phenotype.output,
            antibiotic = get_all_clustered_antibiotics(wc)
        ),
        lambda wc: expand(
            rules.panpa_align_cluster_single_target.output,
            antibiotic = get_all_clustered_antibiotics(wc)
        ),
    output: touch(TEMP_DIR / "flags" / "gather_panpa_all_clusters.done")

# -----------------------
# Panproteome Graph: PanPA
# -----------------------

rule panpa_build_index:
    input: rules.gather_align_clusters.output
    output: OUT_DIR / "panpa" / "index.pickle"
    log: TEMP_DIR / "logs" / "panpa_build_index.log"
    benchmark: TEMP_DIR / "benchmarks" / "panpa_build_index.tsv"
    params:
        kmer_size = 10,
        window_size = 15,
        seed_limit = 0,
    conda: "envs/panpa.yaml"
    threads: 1
    shell:
        r"""
        PanPA build_index \
            --in_dir {input} \
            --out_index {output} \
            --seeding_alg wk_min \
            --kmer_size {params.kmer_size} \
            --window {params.window_size} \
            --seed_limit {params.seed_limit} \
            >> {log} 2>&1
        """

rule panpa_build_gfa:
    input: rules.gather_align_clusters.output
    output: directory(OUT_DIR / "panpa" / "gfa")
    log: TEMP_DIR / "logs" / "panpa_build_gfa.log"
    benchmark: TEMP_DIR / "benchmarks" / "panpa_build_gfa.tsv"
    conda: "envs/panpa.yaml"
    threads: workflow.cores
    shell:
        r"""
        PanPA build_gfa \
            --in_dir {input} \
            --out_dir {output} \
            --cores {threads} \
            >> {log} 2>&1
        """


# -----------------------
# Binary Gene Presence Absence
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
    output: OUT_DIR / "binary_gpa.tsv"
    shell:
        r"""
        ln -sr {input} {output}
        """


# -----------------------
# Mutation Analysis
# -----------------------

rule snippy_runner:
    input:
        sample_store = rules.rename_files.output.store,
        sample = Path(rules.rename_files.output.store) / "{sample}",
        reference = GBFF_FILE,
    output:
        vcf = OUT_DIR / "snippy" / "{sample}" / "snps.vcf",
        tab = OUT_DIR / "snippy" / "{sample}" / "snps.tab",
    log: TEMP_DIR / "logs" / "snippy_runner" / "{sample}.log"
    benchmark: TEMP_DIR / "benchmarks" / "snippy_{sample}.tsv"
    params:
        out_dir = subpath(output.vcf, parent=True)
    conda: "envs/snippy.yaml"
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


rule binary_mutation_table:
    input:
        lambda wc: expand(
            rules.snippy_runner.output.vcf,
            sample = get_sample_names(wc)
        )
    output: OUT_DIR / "binary_mutation_table.tsv"
    benchmark: TEMP_DIR / "benchmarks" / "binary_mutation_table.py.tsv"
    conda: "envs/python312.yaml"
    threads: MAX_PYTHON_THREADS
    script:
        "scripts/binary_mutation_table.py"


# -----------------------
# Feature Importance Analysis
# -----------------------

rule annotation_file_from_snippy:
    input:
        lambda wc: expand(
            rules.snippy_runner.output.tab,
            sample = get_sample_names(wc)
        )
    output: OUT_DIR / "mutations_annotations.tsv"
    benchmark: TEMP_DIR / "benchmarks" / "annotation_file_from_snippy.py.tsv"
    conda: "envs/python312.yaml"
    threads: MAX_PYTHON_THREADS
    script:
        "scripts/annotation_file_from_snippy.py"


# -----------------------
# Binary Tables
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
        rules.phenotype_dataframe_creator.output,
        rules.annotation_file_from_snippy.output,
        rules.cdhit_protein_positions.output,
        rules.panpa_build_index.output,
        rules.panpa_build_gfa.output,
        rules.gather_panpa_all_clusters.output,
    default_target: True

