# ALPAR command inputs

Every command is `alpar <subcommand> ...`. This page describes **what data each subcommand reads**. Flags that only control runtime (threads, RAM, overwrite, verbosity, temp) are listed briefly; they are not input data.

`alpar automatix` runs the later steps for you. The other commands are the same stages, called one at a time. Outputs of an earlier command become inputs of the next.

```
genomes + reference
        │
        ▼
create_binary_tables  →  mutation TSV, phenotype TSV, strains.txt, annotations
        │
        ├─► binary_table_threshold   (optional filter)
        ├─► phylogenetic_tree / panacota  →  tree
        ├─► gwas   (mutation + phenotype + tree)
        ├─► prps   (mutation + tree)
        └─► ml     (mutation + phenotype [+ PRPS] [+ strains.txt])
                    │
                    ├─► structman   (feature-importance file + annotations + reference)
                    └─► prediction  (new genomes or table + saved model)
```

`create_phenotype_table` is only needed if you did not already build a phenotype table from the folder layout.

---

## Shared genome folder (`-i` for automatix and related commands)

Used by **`automatix`** (required) and optionally by **`create_binary_tables`**, **`create_phenotype_table`**, **`panacota`**, **`phylogenetic_tree`**, **`ml --sail`**, and **`prediction`**.

```
input_folder/
  <antibiotic_name>/
    Resistant/
      strain1.fna
      strain2.fna
    Susceptible/
      strain3.fna
      strain4.fna
  <another_antibiotic>/
    Resistant/
      ...
    Susceptible/
      ...
```

Rules:

- Top-level folder names are **antibiotic names**. They must match later `-a` / phenotype column names (for example `ciprofloxacin`, `amikacin`).
- Each antibiotic folder must contain **`Resistant`** and **`Susceptible`**.
- Genome files are FASTA. Phenotype creation looks for **`.fna`**. The same strain can appear under more than one antibiotic.
- A strain in `Resistant` is coded **1**; in `Susceptible` **0**. If a strain is not listed under an antibiotic, that cell is **2** (not tested / missing).

`create_binary_tables` also accepts a **text file** instead of this folder: one genome path per line. That path list does not encode resistance labels. Use `--create_phenotype_from_folder` plus the folder layout (or `create_phenotype_table` / a hand-made TSV) if you need phenotypes.

### Reference genome (`--reference`)

Required for **`automatix`**, **`create_binary_tables`**, and **`structman`**. Optional for **`prediction`** when new genomes still need variant calling.

Accepted formats: **`.gbk`** or **`.gbff`** (GenBank). This is the single reference used for Snippy and annotation.

### Strain list (`strains.txt`)

Written by **`create_binary_tables`**. One genome path per line. Later used as `-i` for **`panacota`** and **`phylogenetic_tree`**, and as `--sail` for **`ml`**.

### Random names (`random_names.txt`)

Written by **`create_binary_tables`**. Tab-separated: original strain ID → short random ID used inside some tools. Pass it as `--random_names_dict` when you want trees or phenotype tables to use the same IDs.

### Phenotype table (`phenotype_table.tsv`)

Tab-separated. First column is strain ID (or the random ID if names were remapped). Other columns are antibiotics. Values: **1** resistant, **0** susceptible, **2** missing.

You can build this from the folder layout or write it by hand. Examples of the layout (real published labels) live under `tool_results/` (for example `tool_results/CAMDA2025/phenotype_Neisseria_gonorrhoeae.tsv`).

### Binary mutation table

TSV from **`create_binary_tables`**:

- `binary_mutation_table.tsv` — variants only
- `binary_mutation_table_with_gene_presence_absence.tsv` — variants plus gene presence/absence (unless `--only_variants`)

Rows are strains; columns are features (mutations / genes); values are **0/1**. After thresholding, the file is named like `binary_mutation_table_threshold_*_percent.tsv`.

---

## `alpar automatix`

End-to-end: binary tables → threshold → tree → GWAS → PRPS → ML (unless `--no_ml`).

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | Genome folder in the layout above |
| `-o` / `--output` | yes | Output directory |
| `--reference` | yes | Reference `.gbk` / `.gbff` |
| `--prokka_custom_database` | no | Two values: protein FASTA (e.g. UniProt) and genus name |
| `--bakta_db` | no | Existing Bakta database directory (downloaded if missing and Bakta is selected) |

Example:

```shell
alpar automatix -i path/to/input_folder -o path/to/output --reference path/to/reference.gbff
```

Useful switches (not extra data files): `--only_variants`, `--fast` (skip PanACoTA; Mash tree only), `--no_ml`, `--no_datasail`, `--run_qc`, `--ml_algorithm rf svm lr xgb`, `--checkpoint`.

---

## `alpar create_binary_tables`

From genomes, write mutation/phenotype tables and strain lists.

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | Genome folder **or** text file of genome paths (one per line) |
| `-o` / `--output` | yes | Output directory |
| `--reference` | yes | Reference `.gbk` / `.gbff` |
| `--create_phenotype_from_folder` | no | Build `phenotype_table.tsv` from the Resistant/Susceptible folders |
| `--prokka_custom_database` | no | Protein FASTA + genus name |
| `--bakta_db` | no | Bakta database directory |

Main outputs you will feed into later commands:

- `binary_mutation_table.tsv`
- `binary_mutation_table_with_gene_presence_absence.tsv`
- `phenotype_table.tsv` (if created from the folder)
- `strains.txt`
- `random_names.txt`
- `mutations_annotations.tsv`

Example:

```shell
alpar create_binary_tables -i path/to/input_folder -o path/to/output --reference path/to/reference.gbff --create_phenotype_from_folder
```

---

## `alpar create_phenotype_table`

Builds only the phenotype TSV from the genome folder. Skip this if you already used `--create_phenotype_from_folder`.

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | Same antibiotic / Resistant / Susceptible folder |
| `-o` / `--output` | yes | Output directory |
| `--random_names_dict` | no | `random_names.txt` so IDs match the mutation table |

---

## `alpar binary_table_threshold`

Drops rare mutation columns (default threshold 0.2).

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | A binary mutation TSV from `create_binary_tables` |
| `-o` / `--output` | yes | Output directory |
| `--threshold_percentage` | no | Fraction; columns below this are dropped (default `0.2`) |

The filtered TSV can replace the original table as `-i` for **gwas**, **prps**, and **ml**.

---

## `alpar phylogenetic_tree`

Alignment-free tree (Mashtree).

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | `strains.txt` **or** the genome folder |
| `-o` / `--output` | yes | Output directory |
| `--random_names_dict` | no | `random_names.txt` |

Typical tree file for later steps: `phylogenetic_tree.tree` under the phylogeny output folder.

---

## `alpar panacota`

Alignment-based tree (PanACoTA). Slower than Mashtree; `automatix --fast` skips this.

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | `strains.txt` **or** the genome folder |
| `-o` / `--output` | yes | Output directory |
| `--random_names_dict` | no | `random_names.txt` |
| `--data_type` | no | `nucl` (default) or `prot` |

Typical tree file: `phylogenetic_tree.newick`.

---

## `alpar gwas`

Association test (Pyseer). Needs genotypes, labels, and a tree.

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | Binary mutation TSV (raw, with GPA, or thresholded) |
| `-p` / `--phenotype` | yes | `phenotype_table.tsv` |
| `-t` / `--tree` | yes | Tree from `phylogenetic_tree` or `panacota` |
| `-o` / `--output` | yes | Output directory |

---

## `alpar prps`

Phylogeny-Related Parallelism Score: which mutations track the tree more than AMR.

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | Binary mutation TSV |
| `-t` / `--tree` | yes | Same tree as GWAS |
| `-o` / `--output` | yes | Output directory |

The PRPS score file can be passed to **`ml --prps`**.

---

## `alpar ml`

Train a classifier for **one** antibiotic column.

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | Binary mutation TSV |
| `-p` / `--phenotype` | yes | Phenotype TSV |
| `-a` / `--antibiotic` | yes | Column name in the phenotype table (must match the folder/antibiotic name) |
| `-o` / `--output` | yes | Output directory |
| `--annotation` | no | `mutations_annotations.tsv` |
| `--prps` | no | PRPS score file from `alpar prps` |
| `--prps_percentage` | no | Keep top N% of PRPS scores (default `30`; use with `--prps`) |
| `--sail` | no | `strains.txt` or genome folder — DataSAIL split to reduce leakage |
| `--sail_distance_matrix` | no | Precomputed DataSAIL distance matrix |
| `--train_strains_file` / `--test_strains_file` / `--validation_strains_file` | no | Optional explicit strain-ID lists |

`--ml_algorithm` is `rf`, `svm`, `gb`, `histgb`, `xgb`, or `lr` (default `rf`). Use `--save_model` if you will run **`prediction`**. `--feature_importance_analysis` writes the file **`structman`** needs.

Example:

```shell
alpar ml -i output/binary_mutation_table.tsv -p output/phenotype_table.tsv -o ml_out/ -a amikacin --save_model
```

---

## `alpar structman`

Turns feature-importance results into a StructMAn input (not a genome command).

| Flag | Required | What it is |
|---|---|---|
| `-i` / `--input` | yes | One feature-importance file, **or** a text file listing those paths (one per line) |
| `--annotation` | yes | `mutations_annotations.tsv` from `create_binary_tables` |
| `--reference` | yes | Same reference `.gbk` / `.gbff` used for variant calling |
| `-o` / `--output` | yes | Output directory |

---

## `alpar prediction`

Score new strains with saved models.

| Flag | Required | What it is |
|---|---|---|
| `-b` / `--model_binary_table` | yes | Mutation TSV **used to train** the model (feature columns must match) |
| `--model` | yes | One or more saved model files |
| `--model_names` | yes | Labels for those models, same order as `--model` |
| `-o` / `--output` | yes | Output directory |
| `-i` / `--input` | one of `-i` or `-p` | New genomes: `strains.txt` **or** genome folder |
| `-p` / `--prediction_table` | one of `-i` or `-p` | Already-built mutation table for the new strains |
| `--reference` | if using `-i` genomes | Reference `.gbk` / `.gbff` for variant calling |
| `--prps` | no | Two values: PRPS table path and percentage (default percentage `30`) |

If you pass new FASTA via `-i`, the same variant-calling / annotation tools as `create_binary_tables` are used so features line up with `-b`.

---

## Which file goes where

| You have | Pass it to |
|---|---|
| Genome folder (antibiotic / Resistant / Susceptible) | `automatix -i`, `create_binary_tables -i`, `create_phenotype_table -i`, optionally tree / sail / prediction |
| List of FASTA paths | `create_binary_tables -i`, `panacota -i`, `phylogenetic_tree -i`, `ml --sail`, `prediction -i` |
| Reference `.gbk` / `.gbff` | `automatix --reference`, `create_binary_tables --reference`, `structman --reference`, `prediction --reference` |
| `binary_mutation_table*.tsv` | `binary_table_threshold -i`, `gwas -i`, `prps -i`, `ml -i`, `prediction -b` |
| `phenotype_table.tsv` | `gwas -p`, `ml -p` |
| Tree (`.tree` / `.newick`) | `gwas -t`, `prps -t` |
| PRPS scores | `ml --prps`, `prediction --prps` |
| `mutations_annotations.tsv` | `ml --annotation`, `structman --annotation` |
| Feature-importance file | `structman -i` |
| Saved model | `prediction --model` |

For install and WSL notes, see [README.md](README.md). For every flag, run `alpar <subcommand> -h`.

---

## Public data sources

These are the main public genotype–phenotype resources used in CAMDA AMR work (including the BIOTIA-DX RESISTANCE 2026 preprint) and mirrored as result tables under `tool_results/`. They are **not** ALPAR runtime inputs; download genomes and labels from the sites below, then arrange them as described above.

| Source | What it is | Links |
|---|---|---|
| CAMDA | Annual AMR prediction challenge: assembled genomes plus sequestered test labels. 2026 used 800 train / 250 test isolates per species–drug pair. | [CAMDA 2026](https://bipress.boku.ac.at/camda2026/) · [CAMDA 2025](https://bipress.boku.ac.at/camda2025/) · [CAMDA 2024](https://bipress.boku.ac.at/camda-play/the-camda-contest-challenges) |
| BV-BRC | Bacterial and Viral Bioinformatics Resource Center (successor to PATRIC). Public isolate genomes with AMR phenotypes. | [BV-BRC](https://www.bv-brc.org/) · [Olson et al. 2023, *NAR*](https://doi.org/10.1093/nar/gkac1003) |
| CARD | Comprehensive Antibiotic Resistance Database: curated AMR genes, mutations, and ontology (ARO). | [CARD](https://card.mcmaster.ca/) · [Alcock et al. 2023, *NAR*](https://doi.org/10.1093/nar/gkac920) |
| CABBAGE | Comprehensive Assessment of Bacterial-Based AMR prediction from Genotypes. Unified genotype–phenotype database (~170k isolates, ~1.7M pairs) served through EMBL-EBI. | [EMBL-EBI AMR Portal](https://www.ebi.ac.uk/amr/) · [Downloads](https://www.ebi.ac.uk/amr/developers/) · [Kovaka et al., *NAR*](https://doi.org/10.1093/nar/gkag780) · [Preprint](https://doi.org/10.1101/2025.11.12.688105) |
