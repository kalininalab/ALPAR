# SR-AMR

Pipeline for generating binary matrices suitable for machine learning from genomic fasta data.

## Create Binary Tables
- Required input:
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
        │   
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
        │   
        ```
- Required output: Empty folder path
- Required reference: Path of reference file, accepted file formats: '.gbk', '.gbff'

Basic usage:
`python scripts/sr_pipeline.py create_binary_tables -i example/example_files/ -o sr_test_out/ --reference example/reference.gbff`

