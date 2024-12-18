# ALPAR - Automated Learning Pipeline for Antimicrobial Resistance

Automated Learning Pipeline for Antimicrobial Resistance

![Pipeline](./flowcharts/ALPAR_Flowchart_new.jpg?raw=true "Pipeline")

## Installation

Single-Reference AMR is installable from [conda](https://anaconda.org/kalininalab/ALPAR) using [mamba](https://mamba.readthedocs.io/en/latest/installation.html#existing-conda-install>):

To install it into the new environment:

`````shell
mamba create -n alpar -c conda-forge -c kalininalab -c bioconda -c etetoolkit alpar
conda activate alpar
pip install panacota
`````

Or to install it into the already existing environment:

`````shell
mamba install -c conda-forge -c kalininalab -c bioconda -c etetoolkit alpar
pip install panacota
`````

## Example Files

Example files can be downloaded from:

[Example files](https://github.com/kalininalab/ALPAR-Example/blob/main/example.zip)

## Automatic Pipeline

From genomic files, creates binary mutation and phenotype tables, applies thresholds, creates phylogenetic tree, conducts GWAS analysis, calculates PRPS score and trains machine learning models with conducting feature importance analysis and splitting data aginst information leakage with [DataSAIL](https://github.com/kalininalab/DataSAIL) against all the given antibiotics.

- Input, `-i`: Path of folder that have structure: input_folder -> antibiotic -> [Resistant, Susceptible]

    `````shell
    📦input_folder
    ┣ 📂antibiotic1
    ┃ ┣ 📂Resistant
    ┃ ┃ ┣ 📜fasta1.fna
    ┃ ┃ ┗ 📜fasta2.fna
    ┃ ┃ ┗ ...
    ┃ ┗ 📂Susceptible
    ┃ ┃ ┣ 📜fasta3.fna
    ┃ ┃ ┗ 📜fasta4.fna
    ┃ ┃ ┗ ...
    ┗ 📂antibiotic2
    ┃ ┣ 📂Resistant
    ┃ ┃ ┣ 📜fasta2.fna
    ┃ ┃ ┗ 📜fasta5.fna
    ┃ ┃ ┗ ...
    ┃ ┗ 📂Susceptible
    ┃ ┃ ┣ 📜fasta2.fna
    ┃ ┃ ┗ 📜fasta3.fna
    ┃ ┃ ┗ ...
    ┗ 📂...
    `````

- Output, `-o`: Output folder path, where the output will be stored. If path exist, `--overwrite` option can be used to overwrite existing output.

- Reference, `--reference`: Reference file path, accepted file formats are: `.gbk .gbff`

- Custom database (Optional), `--custom_database`: Fasta file path for protein database creation, can be downloaded from [UniProt](https://www.uniprot.org/) accepted file formats are: `.fasta`

Basic usage:
`````shell
alpar automatix -i example/example_files/ -o example/example_output/ --reference example/reference.gbff
`````

## Create Binary Tables

From genomic files, creates binary mutation and phenotype tables

- Input, `-i`: Path of file that contains path of genomic fasta files per line or path of folder that have structure: input_folder -> antibiotic -> [Resistant, Susceptible]

    `````shell
    📦input_folder
    ┣ 📂antibiotic1
    ┃ ┣ 📂Resistant
    ┃ ┃ ┣ 📜fasta1.fna
    ┃ ┃ ┗ 📜fasta2.fna
    ┃ ┃ ┗ ...
    ┃ ┗ 📂Susceptible
    ┃ ┃ ┣ 📜fasta3.fna
    ┃ ┃ ┗ 📜fasta4.fna
    ┃ ┃ ┗ ...
    ┗ 📂antibiotic2
    ┃ ┣ 📂Resistant
    ┃ ┃ ┣ 📜fasta2.fna
    ┃ ┃ ┗ 📜fasta5.fna
    ┃ ┃ ┗ ...
    ┃ ┗ 📂Susceptible
    ┃ ┃ ┣ 📜fasta2.fna
    ┃ ┃ ┗ 📜fasta3.fna
    ┃ ┃ ┗ ...
    ┗ 📂...
    `````

- Output, `-o`: Output folder path, where the output will be stored. If path exist, `--overwrite` option can be used to overwrite existing output.

- Reference, `--reference`: Reference file path, accepted file formats are: `.gbk .gbff`

- Custom database (Optional), `--custom_database`: Fasta file path for protein database creation, can be downloaded from [UniProt](https://www.uniprot.org/) accepted file formats are: `.fasta`

- Creation of phenotype table (Optional):
    - `--create_phenotype_from_folder` should be used
    - Genomes_folder_path should have a structure: input_folder -> antibiotic -> [Resistant, Susceptible] -> genomic fasta files

    `````shell
    📦input_folder
    ┣ 📂antibiotic1
    ┃ ┣ 📂Resistant
    ┃ ┃ ┣ 📜fasta1.fna
    ┃ ┃ ┗ 📜fasta2.fna
    ┃ ┃ ┗ ...
    ┃ ┗ 📂Susceptible
    ┃ ┃ ┣ 📜fasta3.fna
    ┃ ┃ ┗ 📜fasta4.fna
    ┃ ┃ ┗ ...
    ┗ 📂antibiotic2
    ┃ ┣ 📂Resistant
    ┃ ┃ ┣ 📜fasta2.fna
    ┃ ┃ ┗ 📜fasta5.fna
    ┃ ┃ ┗ ...
    ┃ ┗ 📂Susceptible
    ┃ ┃ ┣ 📜fasta2.fna
    ┃ ┃ ┗ 📜fasta3.fna
    ┃ ┃ ┗ ...
    ┗ 📂...
    `````

Basic usage:
`````shell
alpar create_binary_tables -i example/example_files/ -o example/example_output/ --reference example/reference.gbff
`````

## Binary Table Threshold

Applies threshold to binary mutation table, and drops columns that has less than threshold percentage, useful to reduce sequencing errors in the data.

- Input, `-i`: Binary mutation table path

- Output, `-o`: Output folder path, where the output will be stored. If path exist, `--overwrite` option can be used to overwrite existing output.

- Threshold percentage, `--threshold_percentage`: Threshold percentage value to be used to drop columns. If column sum is less than this value, columns will be deleted from table

Basic usage:
`````shell
alpar binary_tables_threshold -i example/example_output/binary_mutation_table.tsv -o example/example_output/ 
`````

## Phylogenetic Tree

Runs Phylogeny pipeline to create phylogenetic tree. (Alignment free)

- Input, `-i`: Text file that contains path of each strain per line. It can be found in create_binary_tables output path as `strains.txt`

- Output, `-o`: Output folder path, where the output will be stored. If path exist, `--overwrite` option can be used to overwrite existing output.

- Random names dictionary path, `--random_names_dict`: Random names text file path. If not provided, strain's original names will be used for phylogenetic tree

Basic usage:
`````shell
alpar phylogenetic_tree -i example/example_output/strains.txt -o example/example_output/ --random_names_dict example/example_output/random_names.txt 
`````

## Panacota

Runs PanACoTA pipeline to create phylogenetic tree. (Alignment based)

- Input, `-i`: Text file that contains path of each strain per line. It can be found in create_binary_tables output path as `strains.txt`

- Output, `-o`: Output folder path, where the output will be stored. If path exist, `--overwrite` option can be used to overwrite existing output.

- Random names dictionary path, `--random_names_dict`: Random names text file path. If not provided, strain's original names will be used for phylogenetic tree

Basic usage:
`````shell
alpar panacota -i example/example_output/strains.txt -o example/example_output/
`````

## GWAS

Runs GWAS analysis to detect important mutations in the data.

- Input, `-i`:  Binary mutation table path that is created via create_binary_tables command, can be found in create_binary_tables output path as `binary_mutation_table_with_gene_presence_absence.tsv` or `binary_mutation_table.tsv` or if threshold applied, can be found in binary_table_threshold output path as `binary_mutation_table_threshold_*_percent.tsv`

- Phenotype, `-p`:  Binary phenotype table path,  can be found in create_binary_tables output path as `phenotype_table.tsv` if `--create_phenotype_from_folder` is used. Can also created manually and used.

- Tree, `-t` : Phylogenetic tree path, can be found in panacota output path as `phylogenetic_tree.newick` or phylogeny output path as `phylogenetic_tree.tree`

- Output, `-o`: Output folder path, where the output will be stored. If path exist, `--overwrite` option can be used to overwrite existing output.

Basic usage:
`````shell
alpar gwas -i example/example_output/binary_mutation_table_with_gene_presence_absence.tsv -p example/example_output/phenotype_table.tsv -t example/example_output/phylogeny/phylogenetic_tree.tree -o example_output/
`````

## PRPS

Runs PRPS (Phylogeny-Related Parallelism Score) to detect the mutations are more likely associated with phylogeny rather than antimicrobial resistance.

- Input, `-i`:  Binary mutation table path that is created via create_binary_tables command, can be found in create_binary_tables output path as `binary_mutation_table.tsv` or if threshold applied, can be found in binary_table_threshold output path as `binary_mutation_table_threshold_*_percent.tsv`

- Tree, `-t` : Phylogenetic tree path, can be found in panacota output path as `phylogenetic_tree.newick` or phylogeny output path as `phylogenetic_tree.tree`

- Output, `-o`: Output folder path, where the output will be stored. If path exist, `--overwrite` option can be used to overwrite existing output.

Basic usage:
`````shell
alpar prps -i example/example_output/binary_mutation_table.tsv -t example/example_output/phylogeny/phylogenetic_tree.tree -o example_output/
`````

## ML

Trains machine learning models with classification algorithms on the data and optimizes.
<br>
Available Classification algorithms: Random Forest, Support Vector Machine and Gradient Boosting

- Input, `-i`:  Binary mutation table path that is created via create_binary_tables command, can be found in create_binary_tables output path as `binary_mutation_table_with_gene_presence_absence.tsv` or `binary_mutation_table.tsv`

- Phenotype, `-p`:  Binary phenotype table path,  can be found in create_binary_tables output path as `phenotype_table.tsv` if `--create_phenotype_from_folder` is used. Can also created manually and used.

- Output, `-o`: Output folder path, where the output will be stored. If path exist, `--overwrite` option can be used to overwrite existing output.

- Antibiotic, `-a`: Antibiotic name that model will be trained. Should match the name with column that represents phenotype in binary phenotype table

- Optional arguments:
    - Machine learning algorithm, `--ml_algorithm`: Classification algorithm to be used, available selections: [rf, svm, gb]
    - Resampling strategy, `--resampling_strategy`: Resampling strategy to be used, available selections: [holdout, cv]
    - Parameter optimization, `--parameter_optimization`: Parameter optimization for model with autosklearn (https://automl.github.io/auto-sklearn/master/index.html)
    - Save model, `--save_model`: Save model
    - Feature importance analysis, `--feature_importance_analysis`: Analyze important features in the model with gini importance (for RF & GB) or permutation importance (for SVM, RF and GB)
    - Datasail, `--sail`: Splits data into training and test sets against information leakage to train better models. Requires text file that contains path of each strain per line. It can be found in create_binary_tables output path as `strains.txt` 

    More optional arguments can be found in help page: 
    `````shell
    python alpar/sr_amr.py ml -h
    `````

Basic usage:
`````shell
alpar ml -i example/example_output/binary_mutation_table.tsv -p example/example_output/phenotype_table.tsv -o example_output/ -a amikacin
`````
