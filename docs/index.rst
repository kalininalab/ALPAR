########
ALPAR
########

ALPAR, short for Automated Learning Pipeline for Antimicrobial Resistance, serves as a tool crafted to generate antimicrobial resistance predicting machine learning models. Licensed under the MIT license, it is open source and conveniently accessible on
`GitHub <https://github.com/kalininalab/ALPAR>`_. Installation is made simple through
`conda <https://anaconda.org/kalininalab/ALPAR>`_.

Install
#######

ALPAR is installable from conda, currently only support Linux:.

.. note::
    It is recommended to use `mamba <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_
    for the installation because conda might not be able to resolve the dependencies of ALPAR successfully.

.. raw:: html
    :file: install.html


Quick Start
#######

Example files
--------------------------
Example files can be downloaded from: `Example files <https://github.com/kalininalab/ALPAR-Example/blob/main/example.zip>`_

Automatic pipeline
--------------------------
From genomic files, creates binary mutation and phenotype tables, applies thresholds, creates phylogenetic tree, conducts GWAS analysis, calculates PRPS score and trains machine learning models with conducting feature importance analysis and splitting data aginst information leakage with `DataSAIL <https://github.com/kalininalab/DataSAIL>`_ against all the given antibiotics.

- Input, `-i`: Path of folder that have structure: `input_folder -> antibiotic -> [Resistant, Susceptible]`_

.. code-block:: shell
    input_folder
    ┣ antibiotic1
    ┃ ┣ Resistant
    ┃ ┃ ┣ fasta1.fna
    ┃ ┃ ┗ fasta2.fna
    ┃ ┃ ┗ ...
    ┃ ┗ Susceptible
    ┃ ┃ ┣ fasta3.fna
    ┃ ┃ ┗ fasta4.fna
    ┃ ┃ ┗ ...
    ┗ antibiotic2
    ┃ ┣ Resistant
    ┃ ┃ ┣ fasta2.fna
    ┃ ┃ ┗ fasta5.fna
    ┃ ┃ ┗ ...
    ┃ ┗ Susceptible
    ┃ ┃ ┣ fasta2.fna
    ┃ ┃ ┗ fasta3.fna
    ┃ ┃ ┗ ...
    ┗ ...

- Output, `-o`: Output folder path, where the output will be stored. If path exist, `--overwrite` option can be used to overwrite existing output.

- Reference, `--reference`: Reference file path, accepted file formats are: `.gbk .gbff`

- Custom database (Optional), `--custom_database`: Fasta file path for protein database creation, can be downloaded from `UniProt <https://www.uniprot.org/>`_ accepted file formats are: `.fasta`

Basic usage:
.. code-block:: shell
    alpar automatix -i example/example_files/ -o example/example_output/ `--reference` example/reference.gbff

For more information about the parameters:
.. code-block:: shell
    alpar automatix -h

.. toctree::
    :maxdepth: 1
    :caption: Workflow

    workflow/subcommands
    posters