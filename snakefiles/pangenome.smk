from pathlib import Path

PANGENOME_OUT_DIR = OUT_DIR / "pangenome"
PANGENOME_LOGS_DIR = PANGENOME_OUT_DIR / "logs"

# Clusters are processed in a fixed number of shards rather than one job per
# cluster. Tens of thousands of per-cluster jobs made Snakemake's checkpoint
# re-evaluation and group validation take hours and many GiB in the controller.
# A static shard count keeps the DAG known up front and needs no checkpoint.
#
# Each shard is one multithreaded job that processes its clusters in parallel
# (for_each_cluster.sh), giving a few large HTCondor jobs. Snakemake job groups
# are not used here: grouping 256 small shards took 22 minutes to plan and then
# deadlocked ("Out of jobs ready to be started", snakemake issue #823).
PANGENOME_SHARD_COUNT = int(config.get("pangenome_shards", 8))
PANGENOME_SHARDS = tuple(f"{index:04d}" for index in range(PANGENOME_SHARD_COUNT))
FOR_EACH_CLUSTER = SCRIPTS_DIR / "for_each_cluster.sh"

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
# In each for_each_cluster.sh snippet, $1 is the cluster and $2 a tool log path.
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
    params:
        runner = FOR_EACH_CLUSTER,
    threads: workflow.cores
    shell:
        r"""
        mkdir -p {output:q}
        export STORE={input.cluster_store:q} OUT={output:q}
        bash {params.runner:q} {input.manifest:q} {threads} {log:q} '
            if [ "$(grep -c "^>" "$STORE/$1")" -eq 1 ]; then
                echo "Single sequence; skipping alignment."
                cp "$STORE/$1" "$OUT/$1"
            else
                mafft --auto --thread 1 "$STORE/$1" > "$OUT/$1"
            fi'
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
    params:
        runner = FOR_EACH_CLUSTER,
    threads: workflow.cores
    shell:
        r"""
        mkdir -p {output:q}
        export ALIGNMENTS={input.alignments:q} OUT={output:q}
        bash {params.runner:q} {input.manifest:q} {threads} {log:q} '
            PanPA --log_file "$2" build_gfa --fasta_files "$ALIGNMENTS/$1" --out_dir "$OUT" --cores 1'
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
    params:
        runner = FOR_EACH_CLUSTER,
    threads: workflow.cores
    shell:
        r"""
        mkdir -p {output:q}
        export GFA={input.gfa_dir:q} OUT={output:q}
        bash {params.runner:q} {input.manifest:q} {threads} {log:q} '
            BubbleGun --log_file "$2" --in_graph "$GFA/$1.gfa" bchains --bubble_json "$OUT/$1.json"
            if [ ! -f "$OUT/$1.json" ]; then
                echo "No bubbles found."
                echo "{{}}" > "$OUT/$1.json"
            fi'
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
    params:
        script = SCRIPTS_DIR / "bubble_features.py",
        runner = FOR_EACH_CLUSTER,
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: workflow.cores
    shell:
        r"""
        mkdir -p {output.output_dir:q} {output.lor_lookup_dir:q}
        export SCRIPT={params.script:q} GFA={input.gfa_dir:q} BUBBLES={input.bubblegun_dir:q} \
            PHENOTYPES={input.phenotype_table:q} ANTIBIOTIC={wildcards.antibiotic:q} \
            OUT={output.output_dir:q} LOOKUP={output.lor_lookup_dir:q}
        # main() logs and swallows exceptions, so a missing output marks a failed cluster.
        bash {params.runner:q} {input.manifest:q} {threads} {log:q} '
            python "$SCRIPT" \
                --gfa-file "$GFA/$1.gfa" --bubble-gun "$BUBBLES/$1.json" \
                --phenotype-table "$PHENOTYPES" --antibiotic "$ANTIBIOTIC" --log-file "$2" \
                --output-file "$OUT/$1.tsv" --lor-lookup-file "$LOOKUP/$1.tsv"
            test -f "$OUT/$1.tsv" && test -f "$LOOKUP/$1.tsv"'
        """


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
    params:
        runner = FOR_EACH_CLUSTER,
    threads: workflow.cores
    shell:
        r"""
        mkdir -p {output:q}
        export GFA={input.gfa_dir:q} QUERIES={input.test_sequence_folder:q} OUT={output:q}
        bash {params.runner:q} {input.manifest:q} {threads} {log:q} '
            query="$QUERIES/$1" gaf="$OUT/$1.gaf"
            if [ ! -f "$query" ]; then
                echo "Missing test FASTA: $query"
                exit 1
            fi
            if [ ! -s "$query" ]; then
                echo "No test sequences; creating an empty GAF."
                : > "$gaf"
                exit 0
            fi
            PanPA --log_file "$2" align_single \
                --gfa_files "$GFA/$1.gfa" --seqs "$query" --cores 1 --out_gaf "$gaf"
            # PanPA can finish successfully without emitting any alignments.
            if [ ! -f "$gaf" ]; then
                : > "$gaf"
            fi'
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
    params:
        script = SCRIPTS_DIR / "gaf_lor_features.py",
        runner = FOR_EACH_CLUSTER,
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: workflow.cores
    shell:
        r"""
        mkdir -p {output.output_dir:q}
        export SCRIPT={params.script:q} GAF={input.gaf_dir:q} LOOKUP={input.lor_lookup_dir:q} \
            OUT={output.output_dir:q}
        bash {params.runner:q} {input.manifest:q} {threads} {log:q} '
            python "$SCRIPT" --gaf-file "$GAF/$1.gaf" --lor-lookup-file "$LOOKUP/$1.tsv" \
                --output-file "$OUT/$1.tsv" --log-file "$2"'
        """


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
