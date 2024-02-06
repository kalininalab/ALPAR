# SR-AMR

Pipeline for generating binary matrices suitable for machine learning from genomic fasta data.

## Installation

```console
git clone https://github.com/kalininalab/SR-AMR.git

cd SR-AMR

conda env create -f SR-AMR_env.yml

conda install snippy vt=0.57721
```

## Create Binary Tables
- -i :
    - For genotype information:
        - File that contains path of genomic fasta files per line
        or
        - Folder that have structure: input_folder -> antibiotic -> [Resistant, Susceptible]
        ```
        input_folder
        └───antibiotic1
        │   └───Resistant
        │   │   │   fasta1.fna
        │   │   │   fasta2.fna
        │   │   │   ...
        │   │
        │   └───Susceptible
        │       │   fasta3.fna
        │       │   fasta4.fna
        │       │   ...
        │   
        └───antibiotic2
        │   └───Resistant
        │   │   │   fasta5.fna
        │   │   │   fasta6.fna
        │   │   │   ...
        │   │
        │   └───Susceptible
        │       │   fasta2.fna
        │       │   fasta1.fna
        │       │   ...
        ```
    - For phenotype information:
        - `--create_phenotype_from_folder` and `--phenotype_folder {Genomes_folder_path}` should used
        - Genomes_folder_path should have structure: input_folder -> antibiotic -> [Resistant, Susceptible] -> genomic fasta files
        ```
        input_folder
        └───antibiotic1
        │   └───Resistant
        │   │   │   fasta1.fna
        │   │   │   fasta2.fna
        │   │   │   ...
        │   │
        │   └───Susceptible
        │       │   fasta3.fna
        │       │   fasta4.fna
        │       │   ...
        │   
        └───antibiotic2
        │   └───Resistant
        │   │   │   fasta5.fna
        │   │   │   fasta6.fna
        │   │   │   ...
        │   │
        │   └───Susceptible
        │       │   fasta2.fna
        │       │   fasta1.fna
        │       │   ...  
        ```
- -o : Output folder path
- --reference : Path of reference file, accepted file formats: '.gbk', '.gbff'

Basic usage:
`python scripts/sr_pipeline.py create_binary_tables -i example/example_files/ -o example_output/ --reference example/reference.gbff`

## Panacota
- -i : txt file that contains path of each strain per line or input folder path, can be found in create_binary_tables output path as `strains.txt`
- -o : Output folder path

Basic usage:
`python scripts/sr_pipeline.py panacota -i example_output/strains.txt -o example_output/`

## GWAS
- -i : binary mutation table path that is created via create_binary_tables command, can be found in create_binary_tables output path as `binary_mutation_table_with_gene_presence_absence.tsv` or `binary_mutation_table.tsv`
- -p : phenotype table path,  can be found in create_binary_tables output path as `phenotype_table.tsv` if `--create_phenotype_from_folder` is used
- -t : phylogenetic tree path, can be found in panacota output path as `tree/ .treefile`
- -o : path of the output folder

Basic usage:
`python scripts/sr_pipeline.py gwas -i example/example_output/binary_mutation_table_with_gene_presence_absence.tsv -p example/example_output/phenotype_table.tsv -t example/example_output/panacota/tree/ .treefile -o example_output/`