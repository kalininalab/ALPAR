from pathlib import Path


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
    log: LOGS_DIR / "split_cluster_by_phenotype.log"
    benchmark: BENCHMARKS_DIR / "split_cluster_by_phenotype.tsv"
    params:
        resistance_status_mapping = RESISTANCE_STATUS_MAPPING
    conda: "envs/python313.yaml"
    threads: MAX_PYTHON_THREADS
    script:
        "scripts/split_cluster_by_phenotype.py"

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
    log: LOGS_DIR / "align_clusters_by_phenotype_and_map_variants" / "{antibiotic}" / "batch_{batch_num}.log"
    benchmark: BENCHMARKS_DIR / "align_clusters_by_phenotype_and_map_variants_{antibiotic}_batch_{batch_num}.tsv"
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
        susceptible_alignments = directory(TEMP_DIR / "susceptible_alignments" / "{antibiotic}"),
        all_alignments = directory(TEMP_DIR / "all_alignments" / "{antibiotic}"),
        map_files = directory(OUT_DIR / "map_files" / "{antibiotic}"),
    log: LOGS_DIR / "gather_align_and_map_clusters_by_batch" / "{antibiotic}.log"
    threads: 1
    shell:
        r"""
        mkdir -p {output.susceptible_alignments} {output.all_alignments} {output.map_files}

        for batch_dir in {input}; do
            for file in "$batch_dir"/*; do
                [ -f "$file" ] || continue
                case "$file" in
                    *_only_susceptible.fasta)
                        ln -srv "$file" {output.susceptible_alignments}/$(basename "$file") >> {log} 2>&1 ;;
                    *_all.fasta)
                        ln -srv "$file" {output.all_alignments}/$(basename "$file") >> {log} 2>&1 ;;
                    *.map)
                        cp -v "$file" {output.map_files}/$(basename "$file") >> {log} 2>&1 ;;
                        # ln -srv "$file" {output.map_files}/$(basename "$file") >> {log} 2>&1 ;;
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
    log: LOGS_DIR / "cluster_panpa_build_index_by_phenotype" / "{antibiotic}.log"
    benchmark: BENCHMARKS_DIR / "cluster_panpa_build_index_by_phenotype_{antibiotic}.tsv"
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
    log: LOGS_DIR / "cluster_panpa_build_gfa" / "{antibiotic}.log"
    benchmark: BENCHMARKS_DIR / "cluster_panpa_build_gfa_{antibiotic}.tsv"
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
    output: directory(OUT_DIR / "cluster_panpa_alignments" / "{antibiotic}")
    log: LOGS_DIR / "panpa_align_cluster_single_target" / "{antibiotic}.log"
    benchmark: BENCHMARKS_DIR / "panpa_align_cluster_single_target_{antibiotic}.log"
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
            antibiotic = ANTIBIOTICS
        ),
        lambda wc: expand(
            rules.panpa_build_gfa_by_phenotype.output,
            antibiotic = ANTIBIOTICS
        ),
        lambda wc: expand(
            rules.panpa_align_cluster_single_target.output,
            antibiotic = ANTIBIOTICS
        ),
    output: touch(TEMP_DIR / "flags" / "gather_panpa_all_clusters.done")
