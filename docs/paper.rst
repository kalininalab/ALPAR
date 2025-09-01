###################
ALPAR Paper
###################

Preprint
========
The ALPAR preprint is available on bioRxiv:

ALPAR: Automated Learning Pipeline for Antimicrobial Resistance (2025).  
DOI: https://doi.org/10.1101/2025.07.08.663126  
Direct link: https://www.biorxiv.org/content/10.1101/2025.07.08.663126v1

Suggested Citation
------------------
Yurtseven et al. 2025. ALPAR: Automated Learning Pipeline for Antimicrobial Resistance. bioRxiv. doi:10.1101/2025.07.08.663126

BibTeX
------
.. code-block:: bibtex

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

Abstract
=====================
The field of machine learning in antimicrobial resistance (AMR) research has experienced rapid growth, fueled by advancements in high-throughput genome sequencing and growing capacity of computational resources. However, the complexity and lack of standardized data preparation and bioinformatic analyses present significant challenges, especially to newcomers to the domain.

In response to these challenges, we introduce ALPAR (Automated Learning Pipeline for Antimicrobial Resistance), a comprehensive AMR data analysis tool covering the entire process from processing of raw genomic data to training machine learning models to interpretation of results. Our method relies on a reproducible pipeline that integrates widely used bioinformatics tools, presenting a simplified, automatic workflow specifically tailored for single-reference AMR analysis. Accepting genomic data in the form of FASTA files as input, ALPAR facilitates generation of machine learning-ready data tables and both training of machine learning and execution of genome-wide association studies (GWAS) experiments. Additionally, our tool offers supplementary functionalities such as phylogeny-based analysis of distribution of mutations, enhancing its utility for researchers. Our tool is accessible through the ALPAR GitHub page (https://github.com/kalininalab/ALPAR) and installable via conda (https://anaconda.org/kalininalab/ALPAR).
