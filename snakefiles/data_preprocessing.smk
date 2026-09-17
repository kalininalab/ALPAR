from pathlib import Path

DATA_PREPROCESSING_OUT_DIR = OUT_DIR / "data_preprocessing"
DATA_PREPROCESSING_LOGS_DIR = DATA_PREPROCESSING_OUT_DIR / "logs"


# -----------------------
# Create link to files in singleton directory
# -----------------------

checkpoint rename_files:
    localrule: True
    input: IN_DIR
    output:
        store = directory(DATA_PREPROCESSING_OUT_DIR / "data_checksum"),
        mapping = DATA_PREPROCESSING_OUT_DIR / "all_files.tsv",
    benchmark: BENCHMARKS_DIR / "rename_files.tsv"
    log: DATA_PREPROCESSING_LOGS_DIR / "rename_files.log"
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
    output: DATA_PREPROCESSING_OUT_DIR / "phenotype_table.tsv"
    benchmark: BENCHMARKS_DIR / "phenotype_dataframe_creator.tsv"
    log: DATA_PREPROCESSING_LOGS_DIR / "phenotype_dataframe_creator.log"
    conda: ENVS_DIR.format("python313")
    container: CONTAINERS.format("python313:1.0.0")
    params:
        resistance_status_mapping = RESISTANCE_STATUS_MAPPING,
        antibiotics = ANTIBIOTICS,
    threads: 1
    script:
        SCRIPTS_DIR / "phenotype_dataframe_creator.py"

rule data_preprocessing:
    localrule: True
    input:
        rules.phenotype_dataframe_creator.output
    output: touch(OUT_DIR / "flags" / "data_preprocessing.done")
