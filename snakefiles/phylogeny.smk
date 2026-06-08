from pathlib import Path


# -----------------------
# Mashtree
# -----------------------

rule mashtree_preprocessor:
    input: Path(rules.rename_files.output.store) / "{sample}",
    output: TEMP_DIR / "mashtree_preprocessor" / "{sample}.fasta",
    log: LOGS_DIR / "mashtree_preprocessor" / "{sample}.log"
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
    output: OUT_DIR / "phylogenetic_tree.dnd",
    benchmark: BENCHMARKS_DIR / "mashtree_runner.tsv",
    log: LOGS_DIR / "mashtree_runner.log"
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
    output: touch(TEMP_DIR / "flags" / "phylogeny.done")
