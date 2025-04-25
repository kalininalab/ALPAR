###################
ALPAR Subcommands
###################

Create Binary Tables
#######

From genomic files, creates binary mutation and phenotype tables.

- **Input**, `-i`: Path of file that contains the path of genomic fasta files per line or path of folder that has the structure: `input_folder -> antibiotic -> [Resistant, Susceptible]`.

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
        └── ...

- **Output**, `-o`: Output folder path, where the output will be stored. If the path exists, the ``--overwrite`` option can be used to overwrite the existing output.

- **Reference**, ``--reference``: Reference file path, accepted file formats are: `.gbk`, `.gbff`.

- **Custom database (Optional)[Highly recommended]**, ``--custom_database``: Fasta file path for protein database creation for Prokka. These protein fasta can be downloaded from `UniProt <https://www.uniprot.org/>`_. Accepted file format is: `.fasta`.

- **Creation of phenotype table (Optional)**:
    - ``--create_phenotype_from_folder`` should be used.
    - `Genomes_folder_path` should have a structure: `input_folder -> antibiotic -> [Resistant, Susceptible] -> genomic fasta files`.

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
        └── ...

- **Threads (Optional)**, ``--threads``: Number of threads to be used. Default value: 1

- **Memory (Optional)**, ``--ram``: Memory to be used. Default value: 4

- **Keep temporary files (Optional)**, ``--keep_temp_files``: Keep temporary files. Default value: False

**Basic usage**:

.. code-block:: shell

    alpar create_binary_tables -i example/example_files/ -o example/example_output/ --reference example/reference.gbff

Binary Table Threshold
#######

Applies a threshold to the binary mutation table and drops columns that have less than the threshold percentage. This option is useful to reduce sequencing errors in the data.

- **Input**, `-i`: Binary mutation table path.

- **Output**, `-o`: Output folder path, where the output will be stored. If the path exists, the ``--overwrite`` option can be used to overwrite the existing output.

- **Threshold percentage (Optional)**, ``--threshold_percentage``: Threshold percentage value to be used to drop columns. If the column sum is less than this value, columns will be deleted from the table. Default value: 0.2

- **Keep temporary files (Optional)**, ``--keep_temp_files``: Keep temporary files. Default value: False

**Basic usage**:

.. code-block:: shell

    alpar binary_tables_threshold -i example/example_output/binary_mutation_table.tsv -o example/example_output/

Phylogenetic Tree
#######

Runs the phylogeny pipeline to create a phylogenetic tree (alignment-free) with MashTree.

- **Input**, `-i`: Text file that contains the path of each strain per line. It can be found in the `create_binary_tables` output path as `strains.txt`.

- **Output**, `-o`: Output folder path, where the output will be stored. If the path exists, the ``--overwrite`` option can be used to overwrite the existing output.

- **Random names dictionary path (Optional)**, ``--random_names_dict``: Random names text file path. If not provided, the strain's original names will be used for the phylogenetic tree.

- **Keep temporary files (Optional)**, ``--keep_temp_files``: Keep temporary files. Default value: False

**Basic usage**:

.. code-block:: shell

    alpar phylogenetic_tree -i example/example_output/strains.txt -o example/example_output/ --random_names_dict example/example_output/random_names.txt

PanACoTA
#######

Runs the PanACoTA pipeline to create a phylogenetic tree (alignment-based). Requires more time and resources than the `phylogenetic_tree` command.

- **Input**, `-i`: Text file that contains the path of each strain per line. It can be found in the `create_binary_tables` output path as `strains.txt`.

- **Output**, `-o`: Output folder path, where the output will be stored. If the path exists, the ``--overwrite`` option can be used to overwrite the existing output.

- **Random names dictionary path (Optional)**, ``--random_names_dict``: Random names text file path. If not provided, the strain's original names will be used for the phylogenetic tree.

- **Keep temporary files (Optional)**, ``--keep_temp_files``: Keep temporary files. Default value: False

**Basic usage**:

.. code-block:: shell

    alpar panacota -i example/example_output/strains.txt -o example/example_output/

GWAS
#######

Runs GWAS analysis to detect important mutations in the data.

- **Input**, `-i`: Binary mutation table path that is created via the `create_binary_tables` command. It can be found in the `create_binary_tables` output path as `binary_mutation_table_with_gene_presence_absence.tsv` or `binary_mutation_table.tsv`. If a threshold is applied, it can be found in the `binary_table_threshold` output path as `binary_mutation_table_threshold_*_percent.tsv`.

- **Phenotype**, `-p`: Binary phenotype table path. It can be found in the `create_binary_tables` output path as `phenotype_table.tsv` if ``--create_phenotype_from_folder`` is used. It can also be created manually and used.

- **Tree**, `-t`: Phylogenetic tree path. It can be found in the `panacota` output path as `phylogenetic_tree.newick` or the `phylogeny` output path as `phylogenetic_tree.tree`.

- **Output**, `-o`: Output folder path, where the output will be stored. If the path exists, the ``--overwrite`` option can be used to overwrite the existing output.

**Basic usage**:

.. code-block:: shell

    alpar gwas -i example/example_output/binary_mutation_table_with_gene_presence_absence.tsv -p example/example_output/phenotype_table.tsv -t example/example_output/phylogeny/phylogenetic_tree.tree -o example_output/

Phylogeny Related Parallelism Score
#######

Runs PRPS (Phylogeny-Related Parallelism Score) to detect mutations that are more likely associated with phylogeny rather than antimicrobial resistance. Introduced in `Yurtseven et al. 2023 <https://doi.org/10.1186/s12866-023-03147-7>`_.

- **Input**, `-i`: Binary mutation table path that is created via the `create_binary_tables` command. It can be found in the `create_binary_tables` output path as `binary_mutation_table_with_gene_presence_absence.tsv` or `binary_mutation_table.tsv`. If a threshold is applied, it can be found in the `binary_table_threshold` output path as `binary_mutation_table_threshold_*_percent.tsv`.

- **Tree**, ``--tree``: Phylogenetic tree path. It can be found in the `panacota` output path as `phylogenetic_tree.newick` or the `phylogeny` output path as `phylogenetic_tree.tree`.

- **Output**, `-o`: Output folder path, where the output will be stored. If the path exists, the ``--overwrite`` option can be used to overwrite the existing output.

- **Threads (Optional)**, ``--threads``: Number of threads to be used. Default value: 1

- **Keep temporary files (Optional)**, ``--keep_temp_files``: Keep temporary files. Default value: False

- **Temporary directory (Optional)**, ``--temp``: Directory where temporary files will be stored. Default value: `temp` forlder in output directory.

**Basic usage**:

.. code-block:: shell

    alpar prps -i example/example_output/binary_mutation_table.tsv -t example/example_output/phylogeny/phylogenetic_tree.tree -o example_output/

Machine Learning
#######

Trains machine learning models with classification algorithms on the data and optimizes them.

Available classification algorithms: Random Forest, Support Vector Machine, and Gradient Boosting.

- **Input**, `-i`: Binary mutation table path that is created via the `create_binary_tables` command. It can be found in the `create_binary_tables` output path as `binary_mutation_table_with_gene_presence_absence.tsv` or `binary_mutation_table.tsv`.

- **Phenotype**, `-p`: Binary phenotype table path. It can be found in the `create_binary_tables` output path as `phenotype_table.tsv` if ``--create_phenotype_from_folder`` is used. It can also be created manually and used.

- **Output**, `-o`: Output folder path, where the output will be stored. If the path exists, the ``--overwrite`` option can be used to overwrite the existing output.

- **Antibiotic**, `-a`: Antibiotic name that the model will be trained on. This should match the name of the column that represents the phenotype in the binary phenotype table. If none is provided, all the columns will be used.

- **Optional arguments**:
    - **Machine learning algorithm**, ``--ml_algorithm``: Classification algorithm to be used. Available options: `[rf, svm, gb]`.
    - **Resampling strategy**, ``--resampling_strategy``: Resampling strategy to be used. Available options: `[holdout, cv]`.
    - **Parameter optimization**, ``--parameter_optimization``: Parameter optimization for the model with `autosklearn <https://automl.github.io/auto-sklearn/master/index.html>`_.
    - **Save model**, ``--save_model``: Save the trained model.
    - **Feature importance analysis**, ``--feature_importance_analysis``: Analyze important features in the model with Gini importance (for RF & GB) or permutation importance (for SVM, RF, and GB).
    - **Datasail**, ``--sail``: Splits data into training and test sets against information leakage to train better models. Requires a text file that contains the path of each strain per line. It can be found in the `create_binary_tables` output path as `strains.txt`.

    More optional arguments can be found in the help page:

    .. code-block:: shell

        alpar ml -h

**Basic usage**:

.. code-block:: shell

    alpar ml -i example/example_output/binary_mutation_table.tsv -p example/example_output/phenotype_table.tsv -o example_output/ -a amikacin