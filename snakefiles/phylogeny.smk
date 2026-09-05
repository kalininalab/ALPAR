from pathlib import Path

PHYLOGENY_OUT_DIR = OUT_DIR / "phylogeny"
PHYLOGENY_LOGS_DIR = PHYLOGENY_OUT_DIR / "logs"


# -----------------------
# Mashtree
# -----------------------

rule mashtree_preprocessor:
    input: Path(rules.rename_files.output.store) / "{sample}",
    output: PHYLOGENY_OUT_DIR / "mashtree_preprocessor" / "{sample}.fasta",
    log: PHYLOGENY_LOGS_DIR / "mashtree_preprocessor" / "{sample}.log"
    threads: 1
    shell:
        r"""
        ln -srv {input} {output} >> {log} 2>&1
        """

rule mashtree_runner:
    input:
        lambda wc: expand(
            rules.mashtree_preprocessor.output[0],
            sample = get_sample_names(wc)
        ),
    output: PHYLOGENY_OUT_DIR / "phylogenetic_tree.dnd",
    benchmark: BENCHMARKS_DIR / "mashtree_runner.tsv",
    log: PHYLOGENY_LOGS_DIR / "mashtree_runner.log"
    conda: ENVS_DIR.format("mashtree")
    threads: workflow.cores
    shell:
        r"""
        mashtree \
            {input} \
            --numcpus {threads} \
            --outtree {output} \
            >> {log} 2>&1
        """

rule phylogeny:
    input: rules.mashtree_runner.output
    output: touch(OUT_DIR / "flags" / "phylogeny.done")
