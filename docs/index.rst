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

Flowchart
#########

.. image:: /imgs/ALPAR_Updated_Flowchart.png

Tool Environments
#################

The ``alpar`` environment only contains the ALPAR package itself. Each external tool runs in its own conda environment, defined in ``sr_amr/envs/``. These are created automatically on first use (``conda env create``) and ALPAR re-runs the relevant step inside them with ``conda run``, so you do not need to install the tools yourself.

.. list-table::
    :header-rows: 1
    :widths: 20 25 55

    * - Environment
      - Used by
      - Key pinned versions
    * - ``alpar-snippy``
      - variant calling
      - snippy 4.6.0
    * - ``alpar-prokka``
      - annotation
      - prokka 1.14.6, perl-xml-simple
    * - ``alpar-bakta``
      - annotation
      - bakta 1.12.0, diamond 2.1.25
    * - ``alpar-panaroo``
      - gene presence-absence
      - panaroo 1.6.0
    * - ``alpar-cd-hit``
      - gene presence-absence
      - cd-hit 4.8.1
    * - ``alpar-mashtree``
      - ``phylogenetic_tree``
      - mashtree 1.4.6, mash 2.3
    * - ``alpar-panacota``
      - ``panacota``
      - PanACoTA 1.4.0, prokka 1.14.6, Python 3.7
    * - ``alpar-pyseer``
      - ``gwas``
      - pyseer 1.4.0
    * - ``alpar-prps``
      - ``prps``
      - ete3 3.1.3
    * - ``alpar-ml``
      - ``ml``
      - scikit-learn 1.8.0, xgboost 3.2.0
    * - ``alpar-datasail``
      - train/test splitting
      - datasail 1.3.0

.. note::
    - Versions are pinned so results do not change between rebuilds. To rebuild an environment, remove it (``conda env remove -n alpar-prokka``) and run ALPAR again.
    - Use conda 24.11; conda 26.x cannot create environments containing pip packages.
    - The main ``alpar`` environment must use Python 3.12 or 3.13. ``alpar-panacota`` runs Python 3.7 and cannot read pickle protocol 5 files written by Python 3.14.

Quick Start
#######

Example files
--------------------------
Example files can be downloaded from: `Example files <https://github.com/kalininalab/ALPAR-Example/blob/main/example.zip>`_

Automatic pipeline
--------------------------
From genomic files, creates binary mutation and phenotype tables, applies thresholds, creates phylogenetic tree, conducts GWAS analysis, calculates PRPS score and trains machine learning models with conducting feature importance analysis and splitting data aginst information leakage with `DataSAIL <https://github.com/kalininalab/DataSAIL>`_ against all the given antibiotics.

- Input, ``-i``: Path of folder that have structure: input_folder -> antibiotic -> [Resistant, Susceptible]

.. code-block:: shell

    input_folder
    ├── antibiotic1
    │   ├── Resistant
    │   │   ├── fasta1.fna
    │   │   ├── fasta2.fna
    │   │   └── ...
    │   └── Susceptible
    │       ├── fasta3.fna
    │       ├── fasta4.fna
    │       └── ...
    ├── antibiotic2
    │   ├── Resistant
    │   │   ├── fasta2.fna
    │   │   ├── fasta5.fna
    │   │   └── ...
    │   └── Susceptible
    │       ├── fasta2.fna
    │       ├── fasta3.fna
    │       └── ...

- Output, ``-o``: Output folder path, where the output will be stored. If path exist, ``--overwrite`` option can be used to overwrite existing output.

- Reference, ``--reference``: Reference file path, accepted file formats are: ``.gbk .gbff``

- Custom database (Optional) [Recommended], ``--custom_database``: Fasta file path for protein database creation, can be downloaded from `UniProt <https://www.uniprot.org/>`_ accepted file formats are: ``.fasta``

Basic usage:
--------------------------

.. code-block:: shell

    alpar automatix -i example/example_files/ -o example/example_output/ --reference example/reference.gbff

For more information about the automatic pipeline parameters:
--------------------------

.. code-block:: shell

    alpar automatix -h

Citation
#######

For citation details, see :doc:`/citation`.  

.. toctree::
    :maxdepth: 1
    :caption: Workflow

    workflow/subcommands
    camda
    posters
    paper
    citation