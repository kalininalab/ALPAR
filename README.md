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

In simple terms, ALPAR is like a smart assistant that reads bacterial DNA and helps researchers predict which antibiotics are likely to fail, potentially supporting better treatment decisions and the fight against antibiotic-resistant infections.

The installable CLI is `alpar`. The Python package on disk is still named `sr_amr` (from the earlier "Single-Reference AMR" name).

## Repository layout

```
ALPAR/
├── sr_amr/                 # Python package and CLI implementation
│   ├── amr.py              # CLI entry: `alpar <subcommand>`
│   ├── full_automatix.py   # End-to-end orchestrator
│   ├── binary_tables.py    # Variant calling → mutation/phenotype tables
│   ├── gwas.py, prps.py, ml.py, prediction.py, ...
│   ├── envs/*.yaml         # Extra conda envs created on demand (snippy, prokka, …)
│   └── card_data/          # CARD antibiotic/pathogen lookup tables
├── setup.py / pixi.toml / environment.yml   # Install metadata
├── recipe/                 # Conda-forge/bioconda packaging
├── docs/ + flowcharts/     # Docs and pipeline diagrams
├── tests/                  # Unit tests
├── tool_results/           # Published CAMDA / BV-BRC / CABBAGE outputs (not runtime input)
└── example/                # Optional synthetic genomes for a smoke test (see Example Files)
```

`pixi.toml` lists `win-64` as a platform. That applies only to the thin Python layer (pandas, numpy, biopython). The full pipeline still requires Linux or macOS.

## How to use

Install on Linux, macOS, or **WSL2 on Windows** (see [Windows](#windows-linuxmacos-only-tools-require-wsl)). Then either run the full pipeline or individual steps.

Typical first run:

```shell
pixi run alpar automatix \
  -i path/to/input_folder \
  -o path/to/output \
  --reference path/to/reference.gbff
```

Arrange genomes like this (`input_folder` → antibiotic → `Resistant` / `Susceptible` → FASTA files). You also need a reference genome in `.gbk` / `.gbff`.

```
input_folder/
  ciprofloxacin/
    Resistant/     strain1.fna, strain2.fna, ...
    Susceptible/   strain3.fna, strain4.fna, ...
  amikacin/
    Resistant/     ...
    Susceptible/   ...
```

| Command | What it does |
|---|---|
| `alpar automatix` | Full pipeline: tables → tree → GWAS → PRPS → ML |
| `alpar create_binary_tables` | Snippy + Prokka/Bakta + CD-HIT/Panaroo → mutation/phenotype TSVs |
| `alpar binary_table_threshold` | Drop rare columns (sequencing-error filter) |
| `alpar phylogenetic_tree` / `alpar panacota` | Mash tree vs alignment-based tree |
| `alpar gwas` / `alpar prps` | Association + phylogeny-confounder score |
| `alpar ml` | Train RF / SVM / GB / XGB, optional DataSAIL split |
| `alpar prediction` | Score new strains with a saved model |

`tool_results/` is a reference for expected phenotype TSV and model-output layout, not something you pass as `-i`. Detailed flags for each command are in the sections below.

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

#### Corporate TLS / custom CA certificates

If `pixi install` fails with a TLS or certificate error (common behind a company proxy or custom CA), tell Pixi to use the **OS certificate store** instead of the bundled Mozilla roots (`webpki`). These two forms do the same thing; they differ only in where the setting lives and how long it lasts.

| | `tls-root-certs = "system"` | `export PIXI_TLS_ROOT_CERTS=system` |
|---|---|---|
| What it is | Line in a Pixi **config file** | Environment variable for that shell/process |
| Effect | Use the system CA store | Same |
| Lifetime | Persists until you edit/remove the file | Lasts for that terminal (or until you add it to your profile) |
| Scope | The config file you put it in (project, user, or machine) | Every `pixi` command in that environment |
| Precedence | Lower than CLI / env | Overrides the config file |

Put the config-file form in a Pixi config, **not** in `pixi.toml`. Usual locations:

- Project: `.pixi/config.toml`
- User (Linux/WSL/macOS): `~/.pixi/config.toml`
- User (Windows): `%USERPROFILE%\.pixi\config.toml`

```toml
tls-root-certs = "system"
```

Temporary (or in your shell profile):

```shell
# Linux, macOS, WSL
export PIXI_TLS_ROOT_CERTS=system
pixi install
```

```powershell
# Windows PowerShell
$env:PIXI_TLS_ROOT_CERTS = "system"
pixi install
```

Same meaning as `pixi install --tls-root-certs system`. If both the file and the env var are set, the env var / CLI flag wins. `SSL_CERT_FILE` / `SSL_CERT_DIR` win over either.

This only matters for the standalone Pixi binary (GitHub / install script). Conda-forge Pixi already uses the system store, so the setting is accepted but ignored.

`PIXI_TLS_ROOT_CERTS` is read only by **`pixi`**, after Pixi is already installed. It does **not** fix `curl` errors while downloading `https://pixi.sh/install.sh` or the GitHub `.tar.gz` inside that script. Those use Ubuntu's CA store. On WSL, `system` means **Linux** certificates, not the Windows/company store that your browser already trusts.

If `curl` fails with `SSL certificate ... self-signed certificate in certificate chain (60)` (typical on a corporate network):

1. `export PIXI_TLS_ROOT_CERTS=system` will not help yet.
2. The official installer may get past `pixi.sh` with `curl -k` and then fail again on the GitHub binary download. Install the binary yourself:

    ```shell
    cd ~
    mkdir -p ~/.pixi/bin
    curl -kL -o /tmp/pixi.tar.gz \
      https://github.com/prefix-dev/pixi/releases/latest/download/pixi-x86_64-unknown-linux-musl.tar.gz
    tar -xzf /tmp/pixi.tar.gz -C ~/.pixi/bin
    chmod +x ~/.pixi/bin/pixi
    export PATH="$HOME/.pixi/bin:$PATH"
    grep -q '.pixi/bin' ~/.bashrc || echo 'export PATH="$HOME/.pixi/bin:$PATH"' >> ~/.bashrc
    pixi --version
    ```

3. Then use `export PIXI_TLS_ROOT_CERTS=system` for `pixi install`. If that still fails with TLS, Ubuntu still lacks the company CA — adding that CA to Linux is the lasting fix (`-k` is only for the download).

If `pixi install` itself fails with `invalid peer certificate: UnknownIssuer` (often while fetching from `conda.anaconda.org`), Pixi is installed but Ubuntu still does not trust the interceptor CA. A warning that the lock file uses an older format (`v6` → run `pixi lock` for `v7`) is harmless and is not the cause.

On many corporate networks the interceptor is **Cato Networks** (`CN=Cato Networks Root CA`). Windows already trusts it; WSL does not. `PIXI_TLS_ROOT_CERTS=system` only works **after** that CA is in the Linux trust store.

In Ubuntu (`~$`):

```shell
# 1. Save the proxy root from the live chain
echo | openssl s_client -showcerts -servername conda.anaconda.org -connect conda.anaconda.org:443 2>/dev/null \
  | awk '/BEGIN CERTIFICATE/{n++} {print > ("/tmp/cato-cert" n ".crt")}'

# Confirm the last file is the proxy root (Issuer and Subject match)
openssl x509 -in /tmp/cato-cert3.crt -noout -subject -issuer

# 2. Trust it in Ubuntu
sudo cp /tmp/cato-cert3.crt /usr/local/share/ca-certificates/cato-networks-root-ca.crt
sudo update-ca-certificates

# 3. Tell Pixi to use those system certs, then retry
export PIXI_TLS_ROOT_CERTS=system
grep -q 'PIXI_TLS_ROOT_CERTS' ~/.bashrc || echo 'export PIXI_TLS_ROOT_CERTS=system' >> ~/.bashrc

cd ~/ALPAR
pixi install
```

If `/tmp/cato-cert3.crt` is missing, list `/tmp/cato-cert*.crt` and use the file whose `subject` and `issuer` both name the company/proxy root (for Cato, `Cato Networks Root CA`).

Optional, so you do not have to export every time:

```shell
mkdir -p ~/.pixi
printf 'tls-root-certs = "system"\n' >> ~/.pixi/config.toml
```

After `pixi install` works, ignore the lock-format warning or run `pixi lock` later to upgrade v6 → v7.

### Windows: Linux/macOS-only tools require WSL

**You need WSL to run this pipeline on Windows.** Native Windows can browse the code, read `tool_results/`, and maybe install the Python-only deps via pixi. It cannot run variant calling, annotation, trees, or the full pipeline.

`automatix` / `create_binary_tables` create extra conda environments on demand for tools such as `snippy`, `prokka`, `cd-hit`, `panaroo`, `mashtree`, `bakta`, `pyseer`, and `panacota` (see [sr_amr/envs](sr_amr/envs)). Those packages and their dependencies (e.g. `bcftools`, `aragorn`) are **only published for Linux/macOS on bioconda and have no `win-64` build**. Environment creation on native Windows fails with errors like `nothing provides bcftools` or `PackagesNotFoundError`.

To run ALPAR on a Windows machine, use **WSL2 (Windows Subsystem for Linux)**.

**Where to run each command:** PowerShell is only for installing WSL/Ubuntu. The Pixi installer (`curl ... | sh`), `pixi install`, and `alpar` must run in an **Ubuntu** terminal, not PowerShell and not a Cursor PowerShell terminal. That `curl` script is a Linux installer.

1. Install WSL2 with a Linux distribution (run in Windows PowerShell as Administrator):

    `````powershell
    wsl --install -d Ubuntu
    `````

    Let it finish. Restart if Windows asks. After reboot, Ubuntu often opens by itself and asks you to create a Linux username and password. Do that first. Confirm a distro is present with `wsl -l -v`.

2. Open the Ubuntu terminal (you want a prompt like `user@pc:~$`, not `PS C:\...>`):

    | Method | What to do |
    |---|---|
    | Start menu | Search **Ubuntu** and open it |
    | PowerShell | `wsl` or `wsl -d Ubuntu` |
    | Cursor | Terminal dropdown → **New Terminal** → pick a **WSL / Ubuntu** profile (not PowerShell) |

    Cursor is optional. Use it only after Ubuntu exists, and only if the terminal is a WSL profile. If the prompt still starts with `PS`, you are in PowerShell — switch profiles.

3. Inside that Ubuntu terminal, install pixi (or mamba/conda). **Do not** run this from PowerShell:

    `````shell
    curl -fsSL https://pixi.sh/install.sh | sh
    exec $SHELL
    `````

4. Clone the project **into the Linux filesystem** and install from there:

    `````shell
    git clone <this-repo-url> ~/ALPAR
    cd ~/ALPAR
    pixi install
    pixi run alpar --help
    `````

    Prefer `~/ALPAR` over `/mnt/c/...`. The Windows filesystem under WSL is slower, and OneDrive-synced folders in particular can break long bioinformatics jobs. If you must use the Windows checkout, the path looks like `/mnt/c/Users/<you>/.../ALPAR`.

#### Reading the Ubuntu prompt

A prompt like `billchung@CEPSNYLPENG1754:~$` means Ubuntu is running (not PowerShell):

```
user@computer:folder$
│    │         │
│    │         └─ current folder (`~` is Linux home)
│    └─ Windows computer name
└─ Linux username
```

| Prompt | Meaning |
|---|---|
| `user@pc:~$` | Ubuntu, in Linux home (`/home/user`) |
| `user@pc:/mnt/c/Users/...$` | Ubuntu, but standing on the Windows disk |
| `PS C:\Users\...>` | PowerShell — do not run `curl ... install.sh` here |

`~` is `/home/<linux-user>`. That is **not** your Windows user folder. Ubuntu can still *reach* Windows files here:

```
/mnt/c/Users/<windows-user>
```

Examples:

```shell
ls ~
ls /mnt/c/Users/<windows-user>
ls "/mnt/c/Users/<windows-user>/OneDrive - Danaher/Documents/GitHub/ALPAR"
```

Cursor can open a WSL terminal already sitting in the OneDrive path. That is still Ubuntu, so you can install Pixi there, but `cd ~` before you clone or run the pipeline.

Windows can see Linux files in File Explorer at `\\wsl$\Ubuntu\home\<linux-user>`.

#### Windows copy vs Ubuntu copy

A clone or copy in Ubuntu (`~/ALPAR`) is a **second, separate folder** from the ALPAR that Cursor has open on Windows (for example `C:\Users\...\Documents\GitHub\ALPAR`). They are not the same files.

| | Windows (Cursor / OneDrive) | Ubuntu |
|---|---|---|
| Path | `C:\Users\...\GitHub\ALPAR` | `/home/<user>/ALPAR` (`~/ALPAR`) |
| Seen in Ubuntu as | `/mnt/c/Users/.../GitHub/ALPAR` | `~/ALPAR` |
| Same files? | No — two copies | No |

Edits in one do **not** appear in the other unless you `git push` / `git pull` or copy files yourself.

If you already have the Windows repo, you do not have to clone again. Copying still creates a second copy:

```shell
mkdir -p ~/ALPAR
cp -a "/mnt/c/Users/<you>/OneDrive - Danaher/Documents/GitHub/ALPAR/." ~/ALPAR/
```

Pick one folder to edit day to day so the two copies do not drift.

#### Pushing the Ubuntu copy to GitHub

Yes. The Ubuntu copy can push to GitHub the same way the Windows copy can, as long as it is a git repo with a remote.

- **`git clone` into `~/ALPAR`:** `origin` is already set. After commits: `git push -u origin HEAD`
- **`cp -a` from OneDrive including `.git`:** remotes and history come along. You can push from either folder.
- **No `.git` folder:** it is just files. Clone instead, or `git init` and `git remote add origin <url>` — that is a new repo, not a continuation of the Windows one.

You need GitHub login **inside Ubuntu** (HTTPS prompt, `gh auth login`, or an SSH key in WSL). Windows credentials do not always carry over.

Do not push both copies to the same branch with different un-pulled commits. Pick one folder to commit from, or always `git pull` before `git push`.

## Example Files

Example files can be downloaded from:

[Example files](https://www.bv-brc.org/)

This repository can also include a small [example/](example/) folder with **synthetic placeholder genomes** (randomly generated DNA, not real bacterial sequences) so you can test that the pipeline and CLI commands run end-to-end. That folder is listed in `.gitignore`, so it may be missing from a fresh clone. Strain IDs and Resistant/Susceptible labels are copied from the real [tool_results/CAMDA2025/phenotype_Neisseria_gonorrhoeae.tsv](tool_results/CAMDA2025/phenotype_Neisseria_gonorrhoeae.tsv) table, but the sequence content itself is fake, so results from this example are not scientifically meaningful. See [example/README.md](example/README.md) for details if the folder is present, or run it directly:

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
