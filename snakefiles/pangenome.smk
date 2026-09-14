from pathlib import Path

PANGENOME_OUT_DIR = OUT_DIR / "pangenome"
PANGENOME_LOGS_DIR = PANGENOME_OUT_DIR / "logs"

# ------------------------
# Split Clusters
# ------------------------

checkpoint split_cluster_fasta:
    input:
        cdhit_clstr = rules.cdhit_runner.output.clstr,
        combined_proteins = rules.combine_faa_files.output[0],
    output: directory(PANGENOME_OUT_DIR / "cluster_sequences"),
    log: PANGENOME_LOGS_DIR / "split_cluster_fasta.log",
    benchmark: BENCHMARKS_DIR / "split_cluster_fasta.tsv",
    params:
        file_ext = ".fasta"
    conda: ENVS_DIR.format("python313"),
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1,
    script:
        SCRIPTS_DIR / "split_cluster_fasta.py"

rule cluster_fasta_splits:
    input:
        cluster_store = rules.split_cluster_fasta.output[0],
        datasail_splits = rules.datasail_runner.output[0],
    output: directory(PANGENOME_OUT_DIR / "cluster_fasta_splits" / "{antibiotic}" / "{split_category}"),
    log: PANGENOME_LOGS_DIR / "cluster_fasta_splits" / "{antibiotic}_{split_category}.log",
    benchmark: BENCHMARKS_DIR / "cluster_fasta_splits_{antibiotic}_{split_category}.tsv",
    wildcard_constraints:
        split_category = "train|test",
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "cluster_fasta_splits.py"

# -----------------------
# Multiple Sequence Alignment: MAFFT
# -----------------------

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
    output: directory(PANGENOME_OUT_DIR / "batch_align_clusters" / "batch_{batch_num}")
    log: PANGENOME_LOGS_DIR / "batch_align_clusters" / "batch_{batch_num}.log"
    benchmark: BENCHMARKS_DIR / "batch_align_clusters_batch_{batch_num}.tsv"
    conda: ENVS_DIR.format("mafft")
    container: CONTAINERS.format("mafft:1.0.0")
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
    output: directory(PANGENOME_OUT_DIR / "cluster_alignments")
    log: PANGENOME_LOGS_DIR / "gather_align_clusters.log"
    shell:
        r"""
        mkdir -p {output}
        for batch_dir in {input}; do
            for aln_file in $batch_dir/*.fasta; do
                [ -f "$aln_file" ] || continue
                ln -srv $aln_file {output}/$(basename $aln_file) >> {log} 2>&1
            done
        done
        """


# -----------------------
# Panproteome Graph: PanPA
# -----------------------

rule panpa_build_index:
    input: rules.gather_align_clusters.output
    output: PANGENOME_OUT_DIR / "panpa" / "index.pickle"
    log: PANGENOME_LOGS_DIR / "panpa_build_index.log"
    benchmark: BENCHMARKS_DIR / "panpa_build_index.tsv"
    params:
        kmer_size = 10,
        window_size = 15,
        seed_limit = 0,
    conda: ENVS_DIR.format("panpa-vcf")
    container: CONTAINERS.format("panpa-vcf:1.0.0")
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

checkpoint panpa_build_gfa:
    input: rules.gather_align_clusters.output
    output: directory(PANGENOME_OUT_DIR / "panpa" / "gfa")
    log: PANGENOME_LOGS_DIR / "panpa_build_gfa.log"
    benchmark: BENCHMARKS_DIR / "panpa_build_gfa.tsv"
    conda: ENVS_DIR.format("panpa-vcf")
    container: CONTAINERS.format("panpa-vcf:1.0.0")
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


# -----------------------
# Bubble Annotation
# -----------------------

def get_panpa_graphs(wildcards) -> list[Path]:
    panpa_checkpoint = checkpoints.panpa_build_gfa.get(**wildcards)
    panpa_folder = Path(panpa_checkpoint.output[0])
    return sorted(panpa_folder.glob(f"*.gfa"))

def batched_bubblegun(wildcards) -> list[Path]:
    panpa_graphs = get_panpa_graphs(wildcards)
    batch_num = int(wildcards.batch_num)
    start = batch_num * JOB_BATCH_SIZE
    end = min(start + JOB_BATCH_SIZE, len(panpa_graphs))
    return panpa_graphs[start:end]

rule batched_bubblegun_runner:
    input:
        panpa_graph_folder = rules.panpa_build_gfa.output[0],
        batch_graphs = batched_bubblegun
    output: directory(PANGENOME_OUT_DIR / "bubblegun_batches" / "batch_{batch_num}")
    log: PANGENOME_LOGS_DIR / "batched_bubblegun_runner" / "batch_{batch_num}.log"
    benchmark: BENCHMARKS_DIR / "batched_bubblegun_runner_batch_{batch_num}.tsv"
    conda: ENVS_DIR.format("bubblegun")
    container: CONTAINERS.format("bubblegun:1.0.0")
    threads: 1
    shell:
        r"""
        mkdir -p {output}
        mkdir -p $(dirname {log})

        for i in {input.batch_graphs}; do
            OUT_FILE="{output}/$(basename ${{i%.gfa}}).json"
            TEMP_LOG=$(mktemp --suffix .log)
            echo ">> Processing $i" >> $TEMP_LOG

            BubbleGun \
                --log_file $TEMP_LOG \
                --in_graph $i \
                bchains \
                --bubble_json $OUT_FILE \
                >> $TEMP_LOG 2>&1

            if [ ! -f "$OUT_FILE" ]; then
                echo "No bubbles found in $i." >> $TEMP_LOG
                echo '{{}}' > "$OUT_FILE"
            fi

            cat $TEMP_LOG >> {log}
            rm $TEMP_LOG
        done
        """


# -----------------------
# Bubble Features: Train
# -----------------------

def batched_bubble_features(wildcards) -> list[Path]:
    panpa_graphs = get_panpa_graphs(wildcards)
    batch_num = int(wildcards.batch_num)
    start = batch_num * JOB_BATCH_SIZE
    end = min(start + JOB_BATCH_SIZE, len(panpa_graphs))
    return panpa_graphs[start:end]

rule batch_bubble_features:
    input:
        script                = SCRIPTS_DIR / "bubble_features.py",
        panpa_graph_folder    = rules.panpa_build_gfa.output,
        bubblegun_folder      = rules.batched_bubblegun_runner.output,
        batch_graphs          = batched_bubble_features,
        split_phenotype_table = lambda wildcards: rules.split_phenotype_dataframe.output[0].format(split_category="train", **wildcards),
    output:
        outdir     = directory(PANGENOME_OUT_DIR / "bubble_features_batches" / "batch_{batch_num}" / "{antibiotic}"),
        lor_lookup = directory(PANGENOME_OUT_DIR / "bubble_features_batches" / "batch_{batch_num}" / "{antibiotic}_lor_lookup"),
    log: PANGENOME_LOGS_DIR / "batch_bubble_features" / "batch_{batch_num}_{antibiotic}.log"
    benchmark: BENCHMARKS_DIR / "batch_bubble_features_batch_{batch_num}_{antibiotic}.tsv"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    shell:
        r"""
        mkdir -p {output.outdir}
        mkdir -p {output.lor_lookup}
        mkdir -p $(dirname {log})

        for gfa_file in {input.batch_graphs}; do
            CLUSTER=$(basename "${{gfa_file%.gfa}}")
            echo ">> Processing cluster $CLUSTER" >> {log}

            BUBBLE_JSON={input.bubblegun_folder}/$CLUSTER.json
            OUTPUT_FILE="{output.outdir}/${{CLUSTER}}.tsv"
            LOR_LOOKUP_FILE="{output.lor_lookup}/$CLUSTER.tsv"
            TEMP_LOG=$(mktemp --suffix .log)

            python {input.script} \
                --gfa-file $gfa_file \
                --bubble-gun $BUBBLE_JSON \
                --phenotype-table {input.split_phenotype_table} \
                --log-file $TEMP_LOG \
                --antibiotic {wildcards.antibiotic} \
                --output-file $OUTPUT_FILE \
                --lor-lookup-file $LOR_LOOKUP_FILE \
                >> {log} 2>&1

            cat $TEMP_LOG >> {log}
            rm $TEMP_LOG
        done
        """


rule gather_bubble_features:
    input:
        features = lambda wildcards: expand(
            rules.batch_bubble_features.output.outdir,
            batch_num = range((len(get_panpa_graphs(wildcards)) - 1) // JOB_BATCH_SIZE + 1),
            **wildcards
        ),
        lor_lookup = lambda wildcards: expand(
            rules.batch_bubble_features.output.lor_lookup,
            batch_num = range((len(get_panpa_graphs(wildcards)) - 1) // JOB_BATCH_SIZE + 1),
            **wildcards
        )
    output:
        features = PANGENOME_OUT_DIR / "bubble_features_{antibiotic}.tsv",
        lor_lookup = directory(PANGENOME_OUT_DIR / "bubble_features_{antibiotic}_lor_lookup"),
    log: PANGENOME_LOGS_DIR / "gather_bubble_features_{antibiotic}.log"
    shell:
        r"""
        for batch_dir in {input.features}; do
            cat $batch_dir/*.tsv >> {output.features}
        done

        mkdir -p {output.lor_lookup}
        for batch_dir in {input.lor_lookup}; do
            ln -srv $batch_dir/*.tsv {output.lor_lookup} >> {log} 2>&1
        done
        """


rule batch_bubble_features_complete:
    input:
        features = lambda wildcards: expand(
            rules.batch_bubble_features.output.outdir,
            batch_num=range((len(get_panpa_graphs(wildcards)) - 1) // JOB_BATCH_SIZE + 1),
            **wildcards,
        ),
        lor_lookup = lambda wildcards: expand(
            rules.batch_bubble_features.output.lor_lookup,
            batch_num=range((len(get_panpa_graphs(wildcards)) - 1) // JOB_BATCH_SIZE + 1),
            **wildcards,
        ),
    output: touch(PANGENOME_OUT_DIR / "bubble_features_train_{antibiotic}.done")


# -----------------------
# Bubble Features: Test
# -----------------------

def batched_panpa_graphs(wildcards) -> list[Path]:
    panpa_graphs = get_panpa_graphs(wildcards)
    batch_num = int(wildcards.batch_num)
    start = batch_num * JOB_BATCH_SIZE
    end = min(start + JOB_BATCH_SIZE, len(panpa_graphs))
    return panpa_graphs[start:end]

rule batch_panpa_align:
    input:
        panpa_graph_folder = rules.panpa_build_gfa.output,
        test_sequence_folder = lambda wildcards: rules.cluster_fasta_splits.output[0].format(
                split_category="test",
                **wildcards,
            ),
        batch_graphs = batched_panpa_graphs,
    output: directory(PANGENOME_OUT_DIR / "panpa_alignment_batches" / "{antibiotic}" / "batch_{batch_num}"),
    log: PANGENOME_LOGS_DIR / "batch_panpa_align" / "{antibiotic}_batch_{batch_num}.log",
    benchmark: BENCHMARKS_DIR / "batch_panpa_align_{antibiotic}_{batch_num}.tsv",
    conda: ENVS_DIR.format("panpa-vcf")
    container: CONTAINERS.format("panpa-vcf:1.0.0")
    threads: 1
    shell:
        r"""
        mkdir -p "{output}" "$(dirname {log})"
        : > "{log}"

        for gfa_file in {input.batch_graphs}; do
            CLUSTER=$(basename "${{gfa_file%.gfa}}")
            QUERY_FASTA="{input.test_sequence_folder}/${{CLUSTER}}"
            OUTPUT_GAF="{output}/${{CLUSTER}}.gaf"

            echo ">> Processing $CLUSTER" >> "{log}"

            if [ ! -s "$QUERY_FASTA" ]; then
                echo "No test sequences; creating an empty GAF." >> "{log}"
                : > "$OUTPUT_GAF"
                continue
            fi

            TEMP_LOG=$(mktemp --suffix=.log)

            if PanPA \
                --log_file "$TEMP_LOG" \
                align_single \
                --gfa_files "$gfa_file" \
                --seqs "$QUERY_FASTA" \
                --cores {threads} \
                --out_gaf "$OUTPUT_GAF" \
                >> "$TEMP_LOG" 2>&1
            then
                STATUS=0
            else
                STATUS=$?
            fi

            cat "$TEMP_LOG" >> "{log}"
            rm "$TEMP_LOG"

            if [ "$STATUS" -ne 0 ]; then
                echo "PanPA failed for $CLUSTER." >> "{log}"
                exit "$STATUS"
            fi

            # Ensure that every cluster has a GAF, even if PanPA emitted none.
            if [ ! -f "$OUTPUT_GAF" ]; then
                : > "$OUTPUT_GAF"
            fi
        done
        """

rule batch_gaf_lor_features:
    input:
        gaf_dir = rules.batch_panpa_align.output,
        lor_lookup_dir = rules.batch_bubble_features.output.lor_lookup,
    output: directory(PANGENOME_OUT_DIR / "bubble_features_test_batches" / "batch_{batch_num}" / "{antibiotic}"),
    log: PANGENOME_LOGS_DIR / "batch_gaf_lor_features" / "{antibiotic}_batch_{batch_num}.log",
    benchmark: BENCHMARKS_DIR / "batch_gaf_lor_features_{antibiotic}_{batch_num}.tsv",
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    params:
        script = SCRIPTS_DIR / "gaf_lor_features.py"
    threads: 1
    shell:
        r"""
        mkdir -p "{output}" "$(dirname {log})"
        : > "{log}"

        for gaf_file in "{input.gaf_dir}"/*.gaf; do
            [ -f "$gaf_file" ] || continue

            CLUSTER=$(basename "${{gaf_file%.gaf}}")
            LOOKUP_FILE="{input.lor_lookup_dir}/${{CLUSTER}}.tsv"
            OUTPUT_FILE="{output}/${{CLUSTER}}.tsv"

            echo ">> Mapping $CLUSTER" >> "{log}"
            if [ ! -f "$LOOKUP_FILE" ]; then
                echo "Missing LOR lookup: $LOOKUP_FILE" >> "{log}"
                exit 1
            fi

            TEMP_LOG=$(mktemp --suffix=.log)

            if python "{params.script}" \
                --gaf-file "$gaf_file" \
                --lor-lookup-file "$LOOKUP_FILE" \
                --output-file "$OUTPUT_FILE" \
                --log-file "$TEMP_LOG" \
                >> "{log}" 2>&1
            then
                STATUS=0
            else
                STATUS=$?
            fi

            cat "$TEMP_LOG" >> "{log}"
            rm "$TEMP_LOG"

            if [ "$STATUS" -ne 0 ]; then
                echo "GAF-to-LOR mapping failed for $CLUSTER." >> "{log}"
                exit "$STATUS"
            fi
        done
        """


rule batch_gaf_lor_features_complete:
    input:
        lambda wildcards: expand(
            rules.batch_gaf_lor_features.output,
            batch_num=range((len(get_panpa_graphs(wildcards)) - 1) // JOB_BATCH_SIZE + 1),
            **wildcards,
        )
    output: touch(PANGENOME_OUT_DIR / "bubble_features_test_{antibiotic}.done")


rule gather_bubble_features_test:
    input:
        lambda wildcards: expand(
            rules.batch_gaf_lor_features.output,
            batch_num=range((len(get_panpa_graphs(wildcards)) - 1) // JOB_BATCH_SIZE + 1),
            **wildcards,
        )
    output: PANGENOME_OUT_DIR / "bubble_features_test_{antibiotic}.tsv"
    shell:
        r"""
        : > "{output}"
        for batch_dir in {input}; do
            for tsv in "$batch_dir"/*.tsv; do
                [ -f "$tsv" ] || continue
                cat "$tsv" >> "{output}"
            done
        done
        """


rule gather_panpa_alignments:
    input:
        lambda wildcards: expand(
            rules.batch_panpa_align.output,
            batch_num=range((len(get_panpa_graphs(wildcards)) - 1) // JOB_BATCH_SIZE + 1),
            **wildcards,
        )
    output: directory(PANGENOME_OUT_DIR / "panpa" / "alignments" / "{antibiotic}"),
    log: PANGENOME_LOGS_DIR / "gather_panpa_alignments_{antibiotic}.log",
    shell:
        r"""
        mkdir -p "{output}"
        : > "{log}"

        for batch_dir in {input}; do
            for gaf_file in "$batch_dir"/*.gaf; do
                [ -f "$gaf_file" ] || continue

                ln -srv \
                    "$gaf_file" \
                    "{output}/$(basename "$gaf_file")" \
                    >> "{log}" 2>&1
            done
        done
        """

# -----------------------
# Snakefile Target
# -----------------------

rule pangenome:
    input:
        bubble_features_train = expand(rules.batch_bubble_features_complete.output, antibiotic=ANTIBIOTICS),
        panpa_index = rules.panpa_build_index.output,
        bubble_features_test = expand(rules.batch_gaf_lor_features_complete.output, antibiotic=ANTIBIOTICS),
    output: touch(OUT_DIR / "flags" / "pangenome.done")
