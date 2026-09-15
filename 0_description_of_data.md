# Description of local `data/` folders

This file describes the three datasets under `data/`. That folder is gitignored (see `.gitignore`); nothing here is meant to be pushed to GitHub. Command-line inputs for ALPAR are documented in [Data.md](Data.md).

| Folder | Role |
|---|---|
| `data/cabbage/` | CABBAGE / EMBL-EBI AMR Portal tables (isolate phenotypes + in silico genotypes) |
| `data/camda_amr_challenge/` | CAMDA 2026 assembled genomes, laid out for ALPAR |
| `data/card/` | CARD download: ARO ontology, detection models, reference FASTA, SNPs, CARD:Epi |

---

## cabbage

Source: [CABBAGE](https://doi.org/10.1093/nar/gkag780) via the [EMBL-EBI AMR Portal](https://www.ebi.ac.uk/amr/). Three Apache Parquet tables from the portal “views” (same join keys as [EBI developer docs](https://www.ebi.ac.uk/amr/developers/)).

| File | Rows (this copy) | Columns | What each row is |
|---|---|---|---|
| `phenotype.parquet` | 1,714,486 | 34 | One isolate × one antibiotic AST result |
| `genotype.parquet` | 1,370,205 | 33 | One isolate × one predicted AMR gene/element (AMRFinderPlus / mettannotator) |
| `phenotype_genotype_merged.parquet` | 109,338 | 55 | Phenotype row joined to a matching genotype row |

### How the tables relate

```
phenotype.parquet          genotype.parquet
(BioSample_ID,             (BioSample_ID,
 assembly_ID,               assembly_ID,
 antibiotic_ontology)       antibiotic_ontology)
            \                 /
             \               /
              \             /
     phenotype_genotype_merged.parquet
     (inner join on those three keys)
```

- **`BioSample_ID`** — NCBI/ENA sample (e.g. `SAMEA1028830`).
- **`assembly_ID`** — genome assembly (`GCA_…` or `ERZ…`).
- **`antibiotic_ontology`** — ARO-style antibiotic ID (links toward CARD ARO terms).

The merged table is **much smaller** than either parent: a phenotype row is kept only when the same sample, assembly, and antibiotic also have a genotype call. One isolate can appear many times (many drugs in phenotype; many genes in genotype; many gene×drug pairs in the merge).

### Sample: `phenotype.parquet` (first 5 rows, selected columns)

| BioSample_ID | assembly_ID | organism | antibiotic_name | resistance_phenotype | measurement |
|---|---|---|---|---|---|
| SAMEA1028830 | GCA_001096525.1 | Streptococcus pneumoniae | trimethoprim-sulfamethoxazole | susceptible | (empty) |
| SAMEA1024082 | GCA_001164245.1 | Streptococcus pneumoniae | trimethoprim-sulfamethoxazole | resistant | (empty) |
| SAMEA1024159 | GCA_001132945.1 | Streptococcus pneumoniae | trimethoprim-sulfamethoxazole | resistant | (empty) |
| SAMEA1024568 | GCA_001133945.1 | Streptococcus pneumoniae | trimethoprim-sulfamethoxazole | resistant | (empty) |
| SAMEA1024537 | GCA_001098705.1 | Streptococcus pneumoniae | trimethoprim-sulfamethoxazole | resistant | (empty) |

Other useful columns: `SRA_accession`, `collection_year`, `country`, `host`, `ast_standard`, `laboratory_typing_method`, `measurement_sign`, `measurement_units`.

### Sample: `genotype.parquet` (first 5 rows, selected columns)

| BioSample_ID | assembly_ID | organism | gene_symbol | amr_element_symbol | antibiotic_name |
|---|---|---|---|---|---|
| SAMEA5075807 | ERZ26384870 | Escherichia coli | neo | aph(3')-Ia | kanamycin A |
| SAMEA5075807 | ERZ26384870 | Escherichia coli | tetA | tet(A) | tetracycline |
| SAMEA5075807 | ERZ26384870 | Escherichia coli | dhfrI | dfrA5 | trimethoprim |
| SAMEA5075807 | ERZ26384870 | Escherichia coli | ampC | blaEC | beta-lactam antibiotic |
| SAMEA5075807 | ERZ26384870 | Escherichia coli | dfrD | dfrA51 | trimethoprim |

The first five rows are **five genes in the same isolate**. Other useful columns: `region_start`, `region_end`, `class`, `subclass`, `amrfinderplus_method`, `reference_accession`.

### Sample: `phenotype_genotype_merged.parquet` (first 5 rows, selected columns)

| BioSample_ID | assembly_ID | organism | antibiotic_name | resistance_phenotype | gene_symbol |
|---|---|---|---|---|---|
| SAMN13836999 | GCA_018263815.1 | Escherichia coli | tetracycline | resistant | tetA_1 |
| SAMN13836999 | GCA_018263815.1 | Escherichia coli | tetracycline | resistant | tetA_2 |
| SAMN13836995 | GCA_018263235.1 | Escherichia coli | tetracycline | resistant | tetA_1 |
| SAMN13836995 | GCA_018263235.1 | Escherichia coli | tetracycline | resistant | tetA_2 |
| SAMN18988302 | GCA_018263255.1 | Salmonella enterica | tetracycline | resistant | tetA_2 |

Same isolate can have more than one gene (`tetA_1` and `tetA_2`) for the same drug phenotype.

These tables are **not** ALPAR `-i` input. To train ALPAR you need FASTA plus Resistant/Susceptible folders (see CAMDA below) or a hand-built phenotype TSV.

---

## camda_amr_challenge

Source: [CAMDA 2026 AMR challenge](https://bipress.boku.ac.at/camda2026/the-camda-contest-challenges/#amr). Training isolates are drawn from CABBAGE (see `data/camda_amr_challenge/README.txt`). Testing isolates are a held-out set with **no public labels**.

This copy was renamed so it matches ALPAR: `Resistant` / `Susceptible` (capitalized) and `.fna` (decompressed from `.fa` / `.fa.gz`). Sequences were not edited.

### Layout and how the parts relate

```
camda_amr_challenge/
  README.txt                          # terms, provenance, folder rules
  training_dataset/
    <species>__<antibiotic>/
      Resistant/*.fna                 # phenotype = 1 for that drug
      Susceptible/*.fna               # phenotype = 0 for that drug
  testing_dataset/
    <species>__<antibiotic>/*.fna     # unlabeled; for prediction only
```

- Six species–drug **pairs**. Folder name is `species__antibiotic`.
- **Training** encodes the label in the folder (`Resistant` vs `Susceptible`). Pass `training_dataset/` to `alpar automatix -i` or `alpar create_binary_tables -i`.
- **Testing** has the same six pair folders but **no** Resistant/Susceptible split. Use only as `alpar prediction -i`. CAMDA sequesters test phenotypes for scoring.
- Training and testing genomes for a pair are the same species/drug problem; they are **not** joined by a shared ID table in this folder.

### Counts (this copy)

| Pair folder | Train Resistant | Train Susceptible | Test (unlabeled) |
|---|---|---|---|
| `acinetobacter_baumannii__imipenem` | 400 | 400 | 250 |
| `campylobacter_coli__tetracycline` | 400 | 400 | 250 |
| `campylobacter_jejuni__nalidixic_acid` | 400 | 400 | 250 |
| `escherichia_coli__piperacillin` | 400 | 400 | 250 |
| `klebsiella_pneumoniae__imipenem` | 400 | 400 | 250 |
| `streptococcus_pneumoniae__penicillin` | 400 | 400 | 250 |
| **Total** | **2400** | **2400** | **1500** |

That matches the 2026 challenge design: 800 train (400/400) and 250 test per pair.

### Sample files (first five training names, one pair)

`training_dataset/acinetobacter_baumannii__imipenem/Resistant/`:

| File | FASTA header (first contig) |
|---|---|
| `Aba_IMI_R_0.fna` | `>Aba_IMI_R_0.0` |
| `Aba_IMI_R_1.fna` | (same naming: `>Aba_IMI_R_1.0`) |
| `Aba_IMI_R_2.fna` | … |
| `Aba_IMI_R_3.fna` | … |
| `Aba_IMI_R_4.fna` | … |

First bases of `Aba_IMI_R_0.fna` (assembled contig, not a table):

```
>Aba_IMI_R_0.0
GTCACTTAAATTTGAGTAGATATGAGAAGTCTTAACTTTTAAATCTTTGTGCAACAAAGCCTAATTATACTAAAGAAATC
```

Test files use the same FASTA style without an R/S suffix in the folder, e.g. `testing_dataset/escherichia_coli__piperacillin/Eco_PIP_9.fna` starts `>Eco_PIP_9.0`.

There is no phenotype TSV in this folder; ALPAR builds one from the Resistant/Susceptible directories.

**Use agreement:** CAMDA requires acknowledging CABBAGE/CAMDA and submitting results for the CAMDA proceedings if you analyze this download. Do not redistribute the genomes without checking current CAMDA terms.

---

## card

Source: [CARD](https://card.mcmaster.ca/) (Alcock et al. 2023, *NAR*). This folder is a CARD data download plus CARD:Epi relationship tables. CARD materials have a **commercial-use restriction** (see `CARD-Download-README.txt`); ontologies are CC-BY 4.0 (`Ontology-Download-README.txt`).

CARD is a **reference catalog** (genes, mutations, ontology), not a per-isolate AST table like CABBAGE. CABBAGE `antibiotic_ontology` values are meant to line up with CARD ARO accessions.

### How the files relate

```
aro.tsv / aro.json / aro.obo / aro.owl     ← ARO terms (genes, drugs, mechanisms)
        │
        ├─ aro_index.tsv                   ← each detection model + GenBank accessions + ARO
        ├─ aro_categories.tsv              ← ARO terms grouped (family / drug class / mechanism)
        ├─ aro_categories_index.tsv        ← GenBank accessions × those categories
        ├─ PMID.tsv                        ← literature for each ARO term
        ├─ snps.txt                        ← curated resistance mutations on variant models
        ├─ card.json                       ← full RGI detection models (sequences, cutoffs, SNPs)
        ├─ nucleotide_*.fasta              ← DNA reference sequences by model type
        └─ protein_*.fasta                 ← protein reference sequences by model type

mo.tsv / ro.tsv / ncbi_taxonomy.tsv        ← other ontologies (model types, relations, taxa)
shortname_antibiotics.tsv
shortname_pathogens.tsv                    ← abbreviations used in CARD Short Name

card_epi_relationship.tsv                  ← ARO term ↔ epidemiology term
        │
        └─ card_epi_relationship_pmid.tsv  ← papers supporting each relationship (re_id)
```

Join keys that show up in several files:

- **`ARO Accession`** (`ARO:3000005`) — primary CARD term.
- **`CVTERM ID`** — internal CARD integer; also on PMID and CARD:Epi tables.
- **`CARD Short Name`** — compact gene ID used by RGI (e.g. `vanA`, `Abau_gyrA_FLO`).
- **`Protein Accession` / `DNA Accession`** — GenBank IDs shared by `aro_index.tsv`, `aro_categories_index.tsv`, and FASTA headers.
- **`re_id`** — links `card_epi_relationship.tsv` to `card_epi_relationship_pmid.tsv`.

FASTA files split by **model type** (see `mo.tsv`): homolog models are presence/absence genes (e.g. beta-lactamases); variant models are wild-type references used to map resistance SNPs. Using variant FASTA without SNP screening produces false positives (`CARD-Download-README.txt`).

### Sample: `aro.tsv` (first 5 rows, selected columns)

| Accession | Name | CARD Short Name | Description (truncated) |
|---|---|---|---|
| ARO:3000005 | vanD | vanD | D-Ala-D-Ala ligase homolog; vancomycin and teicoplanin resistance |
| ARO:3000010 | vanA | vanA | D-Ala-D-Lac ligase; isolated from VRE; vancomycin and teicoplanin |
| ARO:3000013 | vanB | vanB | Similar to VanA; vancomycin resistance, not teicoplanin |
| ARO:3000024 | patA | patA | ABC transporter (S. pneumoniae) with PatB; fluoroquinolone resistance |
| ARO:3000025 | patB | patB | ABC transporter (S. pneumoniae) with PatA; fluoroquinolone resistance |

### Sample: `aro_index.tsv` (first 5 rows, selected columns)

| ARO Accession | ARO Name | DNA Accession | AMR Gene Family | Drug Class | Resistance Mechanism | CARD Short Name |
|---|---|---|---|---|---|---|
| ARO:3005099 | 23S rRNA methyltransferase Erm(A) | AF002716.1 | Erm-like 23S rRNA methyltransferase | lincosamide; macrolide; streptogramin | antibiotic target alteration | Spyo_ErmA_MLSb |
| ARO:3002523 | AAC(2')-Ia | L06156.2 | AAC(2') | aminoglycoside antibiotic | antibiotic inactivation | AAC(2')-Ia |
| ARO:3002524 | AAC(2')-Ib | U41471.1 | AAC(2') | aminoglycoside antibiotic | antibiotic inactivation | AAC(2')-Ib |
| ARO:3002525 | AAC(2')-Ic | AL123456.3 | AAC(2') | aminoglycoside antibiotic | antibiotic inactivation | AAC(2')-Ic |
| ARO:3002526 | AAC(2')-Id | U72743.1 | AAC(2') | aminoglycoside antibiotic | antibiotic inactivation | AAC(2')-Id |

`aro_categories_index.tsv` repeats Protein/DNA accession plus family, drug class, and mechanism (no ARO column). `aro_categories.tsv` is the category vocabulary (`AMR Gene Family`, `ARO:3004272`, `16S rRNA methyltransferase (A1408)`, …).

### Sample: `snps.txt` (first 5 rows, selected columns)

| Accession | Name | Model Type | Mutations | CARD Short Name |
|---|---|---|---|---|
| 3003817 | A. baumannii gyrA conferring resistance to fluoroquinolones | protein variant model | S81L | Abau_gyrA_FLO |
| 3003817 | A. baumannii gyrA conferring resistance to fluoroquinolones | protein variant model | G79C | Abau_gyrA_FLO |
| 3003818 | A. baumannii parC conferring resistance to fluoroquinolones | protein variant model | S84L | Abau_parC_FLO |
| 3003818 | A. baumannii parC conferring resistance to fluoroquinolones | protein variant model | V104I | Abau_parC_FLO |
| 3003818 | A. baumannii parC conferring resistance to fluoroquinolones | protein variant model | D105E | Abau_parC_FLO |

One ARO/model can have many mutation rows.

### Sample: `PMID.tsv` (first 5 rows, selected columns)

| ARO Accession | ARO Name | PMID |
|---|---|---|
| ARO:3003389 | M. leprae folP with mutation conferring resistance to dapsone | 11709358;15603834;21115799 |
| ARO:3003390 | E. coli ompF with mutation conferring resistance to beta-lactams | 10639355 |
| ARO:3003057 | smeF | 11709330 |
| ARO:3003392 | M. tuberculosis katG mutations conferring resistance to isoniazid | 35944069;36635309;… |
| ARO:3003393 | M. tuberculosis inhA mutations conferring resistance to isoniazid | 35944069;30337678;… |

### Sample: CARD:Epi (`card_epi_relationship.tsv`, first 5 rows)

| re_id | aro_accession | cvterm_name | epi_lexicon | epi_normalized_term_name | paper_count |
|---|---|---|---|---|---|
| 1 | ARO:3003976 | 16S rRNA with mutation conferring resistance to pactamycin | NCBI_TAXON | Dichelobacter nodosus | 1 |
| 2 | ARO:3003976 | 16S rRNA with mutation conferring resistance to pactamycin | NCBI_TAXON | Mycobacterium tuberculosis | 1 |
| 3 | ARO:3003976 | 16S rRNA with mutation conferring resistance to pactamycin | NCBI_TAXON | Mycoplasma bovis | 1 |
| 4 | ARO:3003976 | 16S rRNA with mutation conferring resistance to pactamycin | SO | gene | 1 |
| 5 | ARO:3003976 | 16S rRNA with mutation conferring resistance to pactamycin | SO | insertion_sequence | 1 |

`card_epi_relationship_pmid.tsv` adds the papers for each `re_id` (`pmid`, `pub_year`, `sentence_index`).

### Sample: short-name tables (first 5 rows)

`shortname_antibiotics.tsv`:

| AAC Abbreviation | Molecule |
|---|---|
| AMG | Aminoglycosides |
| AMK | Amikacin |
| AMU | Aminocoumarin |
| AMX | Amoxicillin |
| ATM | Aztreonam |

`shortname_pathogens.tsv`:

| Abbreviation | Pathogen |
|---|---|
| Abau | Acinetobacter baumannii |
| Acla | Alkalihalobacillus clausii |
| Afab | Agrobacterium fabrum |
| Afum | Aspergillus fumigatus |
| Bado | Bifidobacterium adolescentis |

### FASTA headers (protein homolog DNA file)

`nucleotide_fasta_protein_homolog_model.fasta` headers embed GenBank, strand, coordinates, ARO, and gene name — the same ARO/accession space as `aro_index.tsv`:

```
>gb|GQ343019.1|+|132-1023|ARO:3002999|CblA-1 [mixed culture bacterium AX_gF3SD01_15]
>gb|HQ845196.1|+|0-861|ARO:3001109|SHV-52 [Klebsiella pneumoniae]
>gb|AF028812.1|+|392-887|ARO:3002867|dfrF [Enterococcus faecalis]
>gb|JX017365.1|+|244-1120|ARO:3001989|CTX-M-130 [Escherichia coli]
```

`card.json` is the machine-readable bundle RGI uses (models + sequences + SNPs). Ontology twins (`aro.json` / `.obo` / `.owl`, plus `mo`, `ro`, `ncbi_taxonomy`, `viro`) are the same concepts in other formats, not extra isolates.

ALPAR also ships a small CARD subset under `sr_amr/card_data/` (`snps.txt`, `shortname_antibiotics.tsv`, `shortname_pathogens.tsv`). The files in `data/card/` are the full download.
