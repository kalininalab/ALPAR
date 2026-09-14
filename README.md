# ALPAR - Automated Learning Pipeline for Antimicrobial Resistance

Automated Learning Pipeline for Antimicrobial Resistance

![Pipeline](docs/imgs/ALPAR_Flowchart_new.png?raw=true "Pipeline")

`Antimicrobial resistance (AMR)` happens when bacteria evolve so that antibiotics no longer work against them, making infections harder to treat.

This paper introduces ALPAR, a computer tool that combines genetics, bioinformatics, and machine learning to predict whether bacteria will be resistant to antibiotics.

Normally, studying antibiotic resistance requires many different software programs and technical steps, which can be difficult and time-consuming. ALPAR automates most of the process in a single pipeline.

Researchers can give ALPAR bacterial genome sequences, and it automatically finds genetic mutations, builds data tables, trains machine learning models, and identifies genes linked to resistance.

The system also includes safeguards to reduce common mistakes, such as models learning from closely related bacteria instead of finding the true causes of resistance.

To test ALPAR, the researchers analyzed about 1,000 E. coli genomes and used machine learning to predict resistance to the antibiotic ciprofloxacin.

The tool successfully identified well-known resistance mutations in genes such as gyrA and parC, showing that it was finding biologically meaningful patterns.

ALPAR performed much better than a traditional genome-wide association study (GWAS) approach on this dataset.

The tool was also tested in international antimicrobial resistance prediction competitions (CAMDA), where it won the 2024 challenge and placed third in the 2025 challenge.

When compared with an established rule-based resistance prediction tool, ALPAR achieved better overall prediction performance across many bacterial species and antibiotics.

The authors conclude that ALPAR can help scientists more quickly discover resistance-related mutations and make accurate predictions from bacterial DNA data.

In simple terms, ALPAR is like a smart assistant that reads bacterial DNA and helps researchers predict which antibiotics are likely to fail, potentially supporting better treatment decisions and the fight against antibiotic-resistant infections

## Folder Strucutre

ALPAR/
├── sr_amr/                 # The actual program (Python package)
│   ├── amr.py              # CLI entry point → `alpar` command
│   ├── full_automatix.py   # End-to-end pipeline
│   ├── binary_tables.py    # FASTA → mutation/phenotype tables
│   ├── gwas.py, ml.py, prps.py, panacota.py, ...
│   ├── envs/               # One conda env per tool (snippy, prokka, ml, …)
│   └── card_data/          # CARD antibiotic/pathogen lookup tables
├── docs/                   # Sphinx docs (install + subcommands)
├── recipe/                 # Conda package recipe
├── tests/
├── tool_results/           # Precomputed phenotype TSVs (CAMDA, CABBAGE, BV-BRC)
├── flowcharts/
├── setup.py                # pip install → creates `alpar` command
├── environment.yml         # Lightweight Python deps
└── pixi.toml               # Local-dev env (Python only; not the Linux tools)

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

### Installing with pixi

Alternatively, this repository includes a [pixi.toml](pixi.toml) manifest that installs the dependencies listed in [environment.yml](environment.yml) via [pixi](https://pixi.sh/):

`````shell
pixi install
`````

With pixi, there is no separate `conda activate` step:

- Run one-off commands with `pixi run <command>`, e.g. `pixi run alpar --help`. Pixi resolves/activates the environment automatically for that command.
- Or run `pixi shell` once per terminal session to drop into a shell with the environment already active, similar to `conda activate`.

### Windows: Linux/macOS-only tools require WSL

ALPAR's `automatix`/`create_binary_tables` pipeline creates additional conda environments on demand for tools such as `snippy`, `prokka`, `cd-hit`, `panaroo`, `mashtree`, `bakta`, `pyseer`, and `panacota` (see [sr_amr/envs](sr_amr/envs)). Several of these packages, and their dependencies (e.g. `bcftools`, `aragorn`), are **only published for Linux/macOS on bioconda and have no `win-64` build**. As a result, the pipeline cannot run natively on Windows — environment creation for those tools will fail with errors like `nothing provides bcftools` or `PackagesNotFoundError`.

To run ALPAR on a Windows machine, use **WSL2 (Windows Subsystem for Linux)**:

1. Install WSL2 with a Linux distribution (run in Windows PowerShell as Administrator):

    `````powershell
    wsl --install -d Ubuntu
    `````

    Restart if prompted, then launch "Ubuntu" from the Start menu and finish the first-run setup (create a Linux username/password).

2. Inside the WSL Ubuntu terminal, install pixi:

    `````shell
    curl -fsSL https://pixi.sh/install.sh | sh
    exec $SHELL
    `````

3. Open this project **from within WSL** (either clone it inside the Linux filesystem, e.g. `~/ALPAR`, or open the existing Windows checkout via its `/mnt/c/...` path), then install and run as usual:

    `````shell
    cd ~/ALPAR   # or: cd "/mnt/c/Users/<you>/OneDrive - Danaher/Documents/GitHub/ALPAR"
    pixi install
    pixi run alpar automatix -i example/example_files/ -o example/example_output/ --reference example/reference.gbff
    `````

    Cloning/copying the project into the native Linux filesystem (e.g. `~/ALPAR`) instead of `/mnt/c/...` is recommended for better performance.

## Example Files

Example files can be downloaded from:

[Example files](https://www.bv-brc.org/)

This repository also ships a small [example/](example/) folder with **synthetic placeholder genomes** (randomly generated DNA, not real bacterial sequences) so you can test that the pipeline and CLI commands run end-to-end. Strain IDs and Resistant/Susceptible labels are copied from the real [tool_results/CAMDA2025/phenotype_Neisseria_gonorrhoeae.tsv](tool_results/CAMDA2025/phenotype_Neisseria_gonorrhoeae.tsv) table, but the sequence content itself is fake, so results from this example are not scientifically meaningful. See [example/README.md](example/README.md) for details, or run it directly:

`````shell
pixi run alpar automatix -i example/example_files/ -o example/example_output/ --reference example/reference.gbff
`````

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

## Citation

If you use ALPAR, please cite the preprint:

Yurtseven et al. ALPAR: Automated Learning Pipeline for Antimicrobial Resistance. bioRxiv (2025). doi:10.1101/2025.07.08.663126  
Preprint: https://www.biorxiv.org/content/10.1101/2025.07.08.663126v1

### Plain text
Yurtseven et al (2025). ALPAR: Automated Learning Pipeline for Antimicrobial Resistance. bioRxiv. doi:10.1101/2025.07.08.663126

### BibTeX
```bibtex
@article{ALPAR_2025_preprint,
  title   = {ALPAR: Automated Learning Pipeline for Antimicrobial Resistance},
  author  = {Yurtseven, Alper and Joeres, Roman and Kalinina, Olga V.},
  year    = {2025},
  journal = {bioRxiv},
  publisher = {Cold Spring Harbor Laboratory},
  doi     = {10.1101/2025.07.08.663126},
  url     = {https://www.biorxiv.org/content/10.1101/2025.07.08.663126v1},
  note    = {Preprint. Not peer reviewed.}
}
```
