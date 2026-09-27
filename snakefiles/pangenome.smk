from pathlib import Path

PANGENOME_OUT_DIR = OUT_DIR / "pangenome"
PANGENOME_LOGS_DIR = PANGENOME_OUT_DIR / "logs"

# Clusters are processed in a fixed number of shards rather than one job per
# cluster. Tens of thousands of per-cluster jobs made Snakemake's checkpoint
# re-evaluation and group validation take hours and many GiB in the controller.
# A static shard count keeps the DAG known up front and needs no checkpoint.
PANGENOME_SHARD_COUNT = int(config.get("pangenome_shards", 256))
PANGENOME_SHARDS = tuple(f"{index:04d}" for index in range(PANGENOME_SHARD_COUNT))

wildcard_constraints:
    shard = r"\d{4}"

# ------------------------
# Split Clusters
# ------------------------

rule split_cluster_fasta:
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


rule shard_clusters:
    localrule: True
    input: rules.split_cluster_fasta.output[0],
    output: expand(PANGENOME_OUT_DIR / "cluster_shards" / "{shard}.txt", shard=PANGENOME_SHARDS),
    script:
        SCRIPTS_DIR / "shard_clusters.py"


rule cluster_fasta_splits:
    group: "cluster_fasta_splits_batch"
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

# Per-cluster files retain the FASTA filename (e.g. Cluster_0.fasta).
# PanPA appends .gfa to that filename; retaining it also preserves feature IDs.
SHARD_MANIFEST = PANGENOME_OUT_DIR / "cluster_shards" / "{shard}.txt"


# -----------------------
# Multiple Sequence Alignment: MAFFT
# -----------------------

rule align_clusters:
    input:
        manifest = SHARD_MANIFEST,
        cluster_store = rules.split_cluster_fasta.output[0],
    output: directory(PANGENOME_OUT_DIR / "cluster_alignments" / "{shard}")
    log: PANGENOME_LOGS_DIR / "align_clusters" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "align_clusters_{shard}.tsv"
    conda: ENVS_DIR.format("mafft")
    container: CONTAINERS.format("mafft:1.0.0")
    threads: 1
    shell:
        r"""
        mkdir -p {output:q}
        : > {log:q}
        while IFS= read -r cluster; do
            fasta={input.cluster_store:q}/"$cluster"
            echo ">> $cluster" >> {log:q}
            if [ "$(grep -c '^>' "$fasta")" -eq 1 ]; then
                echo "Single sequence; skipping alignment." >> {log:q}
                cp "$fasta" {output:q}/"$cluster"
            else
                mafft --auto --thread {threads} "$fasta" > {output:q}/"$cluster" 2>> {log:q}
            fi
        done < {input.manifest:q}
        """


# -----------------------
# Panproteome Graph: PanPA
# -----------------------

rule panpa_alignment_list:
    localrule: True
    input: expand(rules.align_clusters.output, shard=PANGENOME_SHARDS)
    output: PANGENOME_OUT_DIR / "panpa" / "alignments.txt"
    script:
        SCRIPTS_DIR / "write_input_paths.py"


rule panpa_build_index:
    input:
        alignment_list = rules.panpa_alignment_list.output,
        alignments = expand(rules.align_clusters.output, shard=PANGENOME_SHARDS),
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
            --log_file {log:q} \
            build_index \
            --fasta_list {input.alignment_list:q} \
            --out_index {output:q} \
            --seeding_alg wk_min \
            --kmer_size {params.kmer_size} \
            --window {params.window_size} \
            --seed_limit {params.seed_limit}
        """


rule panpa_build_gfa:
    input:
        manifest = SHARD_MANIFEST,
        alignments = rules.align_clusters.output[0],
    output: directory(PANGENOME_OUT_DIR / "panpa" / "gfa" / "{shard}")
    log: PANGENOME_LOGS_DIR / "panpa_build_gfa" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "panpa_build_gfa_{shard}.tsv"
    conda: ENVS_DIR.format("panpa-vcf")
    container: CONTAINERS.format("panpa-vcf:1.0.0")
    threads: 1
    shell:
        r"""
        mkdir -p {output:q}
        : > {log:q}
        cluster_log=$(mktemp)
        trap 'rm -f "$cluster_log"' EXIT
        while IFS= read -r cluster; do
            echo ">> $cluster" >> {log:q}
            PanPA \
                --log_file "$cluster_log" \
                build_gfa \
                --fasta_files {input.alignments:q}/"$cluster" \
                --out_dir {output:q} \
                --cores {threads}
            cat "$cluster_log" >> {log:q}
        done < {input.manifest:q}
        """


# -----------------------
# Bubble Annotation
# -----------------------

rule bubblegun_runner:
    input:
        manifest = SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa.output[0],
    output: directory(PANGENOME_OUT_DIR / "bubblegun" / "{shard}")
    log: PANGENOME_LOGS_DIR / "bubblegun_runner" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "bubblegun_runner_{shard}.tsv"
    conda: ENVS_DIR.format("bubblegun")
    container: CONTAINERS.format("bubblegun:1.0.0")
    threads: 1
    shell:
        r"""
        mkdir -p {output:q}
        : > {log:q}
        cluster_log=$(mktemp)
        trap 'rm -f "$cluster_log"' EXIT
        while IFS= read -r cluster; do
            json={output:q}/"$cluster.json"
            echo ">> $cluster" >> {log:q}
            BubbleGun \
                --log_file "$cluster_log" \
                --in_graph {input.gfa_dir:q}/"$cluster.gfa" \
                bchains \
                --bubble_json "$json" \
                >> {log:q} 2>&1
            cat "$cluster_log" >> {log:q}

            if [ ! -f "$json" ]; then
                echo "No bubbles found." >> {log:q}
                echo '{{}}' > "$json"
            fi
        done < {input.manifest:q}
        """


# -----------------------
# Bubble Features: Train
# -----------------------

rule bubble_features:
    input:
        manifest = SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa.output[0],
        bubblegun_dir = rules.bubblegun_runner.output[0],
        phenotype_table = lambda wc: rules.split_phenotype_dataframe.output[0].format(
            antibiotic=wc.antibiotic, split_category="train"
        ),
    output:
        output_dir = directory(PANGENOME_OUT_DIR / "bubble_features" / "train" / "{antibiotic}" / "{shard}"),
        lor_lookup_dir = directory(PANGENOME_OUT_DIR / "bubble_features_{antibiotic}_lor_lookup" / "{shard}"),
    log: PANGENOME_LOGS_DIR / "bubble_features" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "bubble_features_{antibiotic}_{shard}.tsv"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "bubble_features_shard.py"


rule merge_bubble_features:
    localrule: True
    input: expand(rules.bubble_features.output.output_dir, shard=PANGENOME_SHARDS, allow_missing=True)
    output: PANGENOME_OUT_DIR / "bubble_features_{antibiotic}.tsv"
    shell:
        "find {input:q} -maxdepth 1 -type f -name '*.tsv' -print0 | sort -z | xargs -0 -r cat > {output:q}"


rule bubble_features_complete:
    localrule: True
    input: expand(rules.bubble_features.output, shard=PANGENOME_SHARDS, allow_missing=True)
    output: touch(PANGENOME_OUT_DIR / "bubble_features_train_{antibiotic}.done")


# -----------------------
# Bubble Features: Test
# -----------------------

rule panpa_align:
    input:
        manifest = SHARD_MANIFEST,
        gfa_dir = rules.panpa_build_gfa.output[0],
        test_sequence_folder = lambda wc: rules.cluster_fasta_splits.output[0].format(
            antibiotic=wc.antibiotic, split_category="test",
        ),
    output: directory(PANGENOME_OUT_DIR / "panpa" / "alignments" / "{antibiotic}" / "{shard}")
    log: PANGENOME_LOGS_DIR / "panpa_align" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "panpa_align_{antibiotic}_{shard}.tsv"
    conda: ENVS_DIR.format("panpa-vcf")
    container: CONTAINERS.format("panpa-vcf:1.0.0")
    threads: 1
    shell:
        r"""
        mkdir -p {output:q}
        : > {log:q}
        cluster_log=$(mktemp)
        trap 'rm -f "$cluster_log"' EXIT
        while IFS= read -r cluster; do
            query_fasta={input.test_sequence_folder:q}/"$cluster"
            gaf={output:q}/"$cluster.gaf"
            echo ">> $cluster" >> {log:q}
            if [ ! -f "$query_fasta" ]; then
                printf 'Missing test FASTA: %s\n' "$query_fasta" >> {log:q}
                exit 1
            fi
            if [ ! -s "$query_fasta" ]; then
                echo "No test sequences; creating an empty GAF." >> {log:q}
                : > "$gaf"
                continue
            fi
            PanPA \
                --log_file "$cluster_log" \
                align_single \
                --gfa_files {input.gfa_dir:q}/"$cluster.gfa" \
                --seqs "$query_fasta" \
                --cores {threads} \
                --out_gaf "$gaf" \
                >> {log:q} 2>&1
            cat "$cluster_log" >> {log:q}

            # PanPA can finish successfully without emitting any alignments.
            if [ ! -f "$gaf" ]; then
                : > "$gaf"
            fi
        done < {input.manifest:q}
        """


rule gaf_lor_features:
    input:
        manifest = SHARD_MANIFEST,
        gaf_dir = rules.panpa_align.output[0],
        lor_lookup_dir = rules.bubble_features.output.lor_lookup_dir,
    output:
        output_dir = directory(PANGENOME_OUT_DIR / "bubble_features" / "test" / "{antibiotic}" / "{shard}"),
    log: PANGENOME_LOGS_DIR / "gaf_lor_features" / "{antibiotic}" / "{shard}.log"
    benchmark: BENCHMARKS_DIR / "gaf_lor_features_{antibiotic}_{shard}.tsv"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "gaf_lor_features_shard.py"


rule gaf_lor_features_complete:
    localrule: True
    input: expand(rules.gaf_lor_features.output, shard=PANGENOME_SHARDS, allow_missing=True)
    output: touch(PANGENOME_OUT_DIR / "bubble_features_test_{antibiotic}.done")


rule merge_bubble_features_test:
    localrule: True
    input: expand(rules.gaf_lor_features.output.output_dir, shard=PANGENOME_SHARDS, allow_missing=True)
    output: PANGENOME_OUT_DIR / "bubble_features_test_{antibiotic}.tsv"
    shell:
        "find {input:q} -maxdepth 1 -type f -name '*.tsv' -print0 | sort -z | xargs -0 -r cat > {output:q}"


# -----------------------
# Snakefile Target
# -----------------------

rule pangenome:
    localrule: True
    input:
        bubble_features_train = expand(rules.bubble_features_complete.output, antibiotic=ANTIBIOTICS),
        panpa_index = rules.panpa_build_index.output,
        bubble_features_test = expand(rules.gaf_lor_features_complete.output, antibiotic=ANTIBIOTICS),
    output: touch(OUT_DIR / "flags" / "pangenome.done")
