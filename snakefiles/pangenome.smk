from functools import lru_cache
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

# Cluster wildcards retain the FASTA filename (e.g. Cluster_0.fasta).
# PanPA appends .gfa to that filename; retaining it also preserves feature IDs.
wildcard_constraints:
    cluster = r"[^/]+\.fasta"


@lru_cache(maxsize=8)
def _pangenome_cluster_names(folder, directory_mtime_ns):
    """Scan each completed checkpoint directory once per generation."""
    return tuple(sorted(path.name for path in folder.glob("*.fasta")))


def get_pangenome_clusters(wildcards):
    # Always consult the checkpoint before the cache so Snakemake can defer DAG
    # expansion. A rebuilt directory invalidates the cached names.
    cluster_checkpoint = checkpoints.split_cluster_fasta.get()
    folder = Path(cluster_checkpoint.output[0])
    return _pangenome_cluster_names(folder, folder.stat().st_mtime_ns)


# -----------------------
# Multiple Sequence Alignment: MAFFT
# -----------------------

rule align_clusters:
    group: "align_clusters"
    input:
        cluster_store = rules.split_cluster_fasta.output[0],
        fasta = PANGENOME_OUT_DIR / "cluster_sequences" / "{cluster}",
    output: PANGENOME_OUT_DIR / "cluster_alignments" / "{cluster}"
    log: PANGENOME_LOGS_DIR / "align_clusters" / "{cluster}.log"
    benchmark: BENCHMARKS_DIR / "align_clusters_{cluster}.tsv"
    conda: ENVS_DIR.format("mafft")
    container: CONTAINERS.format("mafft:1.0.0")
    threads: 1
    shell:
        r"""
        FASTA_COUNT=$(grep -c "^>" {input.fasta:q})
        if [ "$FASTA_COUNT" -eq 1 ]; then
            echo "Single sequence; skipping alignment." > {log:q}
            cp {input.fasta:q} {output:q}
        else
            mafft --auto --thread {threads} {input.fasta:q} > {output:q} 2> {log:q}
        fi
        """


# -----------------------
# Panproteome Graph: PanPA
# -----------------------

rule panpa_alignment_list:
    localrule: True
    input:
        lambda wc: expand(rules.align_clusters.output, cluster=get_pangenome_clusters(wc))
    output: PANGENOME_OUT_DIR / "panpa" / "alignments.txt"
    script:
        SCRIPTS_DIR / "write_input_paths.py"


rule panpa_build_index:
    input:
        alignment_list = rules.panpa_alignment_list.output,
        alignments = lambda wc: expand(rules.align_clusters.output, cluster=get_pangenome_clusters(wc)),
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
    group: "panpa_build_gfa"
    input: rules.align_clusters.output
    output: PANGENOME_OUT_DIR / "panpa" / "gfa" / "{cluster}.gfa"
    log: PANGENOME_LOGS_DIR / "panpa_build_gfa" / "{cluster}.log"
    benchmark: BENCHMARKS_DIR / "panpa_build_gfa_{cluster}.tsv"
    params:
        out_dir = subpath(output[0], parent=True),
    conda: ENVS_DIR.format("panpa-vcf")
    container: CONTAINERS.format("panpa-vcf:1.0.0")
    threads: 1
    shell:
        r"""
        PanPA \
            --log_file {log:q} \
            build_gfa \
            --fasta_files {input:q} \
            --out_dir {params.out_dir:q} \
            --cores {threads}
        """


# -----------------------
# Bubble Annotation
# -----------------------

rule bubblegun_runner:
    group: "bubblegun_runner"
    input: rules.panpa_build_gfa.output
    output: PANGENOME_OUT_DIR / "bubblegun" / "{cluster}.json"
    log: PANGENOME_LOGS_DIR / "bubblegun_runner" / "{cluster}.log"
    benchmark: BENCHMARKS_DIR / "bubblegun_runner_{cluster}.tsv"
    conda: ENVS_DIR.format("bubblegun")
    container: CONTAINERS.format("bubblegun:1.0.0")
    threads: 1
    shell:
        r"""
        BubbleGun \
            --log_file {log:q} \
            --in_graph {input:q} \
            bchains \
            --bubble_json {output:q} \
            >> {log:q} 2>&1

        if [ ! -f {output:q} ]; then
            echo "No bubbles found." >> {log:q}
            echo '{{}}' > {output:q}
        fi
        """


# -----------------------
# Bubble Features: Train
# -----------------------

rule bubble_features:
    group: "bubble_features"
    input:
        gfa_file = rules.panpa_build_gfa.output[0],
        bubble_gun = rules.bubblegun_runner.output[0],
        phenotype_table = lambda wc: rules.split_phenotype_dataframe.output[0].format(
            antibiotic=wc.antibiotic, split_category="train"
        ),
    output:
        output_file = PANGENOME_OUT_DIR / "bubble_features" / "train" / "{antibiotic}" / "{cluster}.tsv",
        lor_lookup_file = PANGENOME_OUT_DIR / "bubble_features_{antibiotic}_lor_lookup" / "{cluster}.tsv",
    log: PANGENOME_LOGS_DIR / "bubble_features" / "{antibiotic}" / "{cluster}.log"
    benchmark: BENCHMARKS_DIR / "bubble_features_{antibiotic}_{cluster}.tsv"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "bubble_features.py"


rule merge_bubble_features:
    localrule: True
    input:
        lambda wc: expand(
            rules.bubble_features.output.output_file,
            cluster=get_pangenome_clusters(wc), antibiotic=wc.antibiotic,
        )
    output: PANGENOME_OUT_DIR / "bubble_features_{antibiotic}.tsv"
    script:
        SCRIPTS_DIR / "concatenate_files.py"


rule bubble_features_complete:
    localrule: True
    input:
        lambda wc: expand(
            rules.bubble_features.output,
            cluster=get_pangenome_clusters(wc), antibiotic=wc.antibiotic,
        )
    output: touch(PANGENOME_OUT_DIR / "bubble_features_train_{antibiotic}.done")


# -----------------------
# Bubble Features: Test
# -----------------------

rule panpa_align:
    group: "panpa_align"
    input:
        gfa_file = rules.panpa_build_gfa.output[0],
        test_sequence_folder = lambda wc: rules.cluster_fasta_splits.output[0].format(
            antibiotic=wc.antibiotic, split_category="test",
        ),
    output: PANGENOME_OUT_DIR / "panpa" / "alignments" / "{antibiotic}" / "{cluster}.gaf"
    log: PANGENOME_LOGS_DIR / "panpa_align" / "{antibiotic}" / "{cluster}.log"
    benchmark: BENCHMARKS_DIR / "panpa_align_{antibiotic}_{cluster}.tsv"
    params:
        query_fasta = lambda wc, input: Path(input.test_sequence_folder) / wc.cluster,
    conda: ENVS_DIR.format("panpa-vcf")
    container: CONTAINERS.format("panpa-vcf:1.0.0")
    threads: 1
    shell:
        r"""
        if [ ! -f {params.query_fasta:q} ]; then
            printf 'Missing test FASTA: %s\n' {params.query_fasta:q} > {log:q}
            exit 1
        fi
        if [ ! -s {params.query_fasta:q} ]; then
            echo "No test sequences; creating an empty GAF." > {log:q}
            : > {output:q}
        else
            PanPA \
                --log_file {log:q} \
                align_single \
                --gfa_files {input.gfa_file:q} \
                --seqs {params.query_fasta:q} \
                --cores {threads} \
                --out_gaf {output:q} \
                >> {log:q} 2>&1

            # PanPA can finish successfully without emitting any alignments.
            if [ ! -f {output:q} ]; then
                : > {output:q}
            fi
        fi
        """


rule gaf_lor_features:
    group: "gaf_lor_features"
    input:
        gaf_file = rules.panpa_align.output[0],
        lor_lookup_file = rules.bubble_features.output.lor_lookup_file,
    output:
        output_file = PANGENOME_OUT_DIR / "bubble_features" / "test" / "{antibiotic}" / "{cluster}.tsv",
    log: PANGENOME_LOGS_DIR / "gaf_lor_features" / "{antibiotic}" / "{cluster}.log"
    benchmark: BENCHMARKS_DIR / "gaf_lor_features_{antibiotic}_{cluster}.tsv"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    threads: 1
    script:
        SCRIPTS_DIR / "gaf_lor_features.py"


rule gaf_lor_features_complete:
    localrule: True
    input:
        lambda wc: expand(
            rules.gaf_lor_features.output,
            cluster=get_pangenome_clusters(wc), antibiotic=wc.antibiotic,
        )
    output: touch(PANGENOME_OUT_DIR / "bubble_features_test_{antibiotic}.done")


rule merge_bubble_features_test:
    localrule: True
    input:
        lambda wc: expand(
            rules.gaf_lor_features.output,
            cluster=get_pangenome_clusters(wc), antibiotic=wc.antibiotic,
        )
    output: PANGENOME_OUT_DIR / "bubble_features_test_{antibiotic}.tsv"
    script:
        SCRIPTS_DIR / "concatenate_files.py"


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
