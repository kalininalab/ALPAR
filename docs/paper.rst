###################
ALPAR Paper
###################

Publication
===========
The ALPAR paper is published in Bioinformatics:

ALPAR: automated learning pipeline for antimicrobial resistance (2026).  
DOI: https://doi.org/10.1093/bioinformatics/btag545  
Direct link: https://academic.oup.com/bioinformatics/article/42/8/btag545/8739913

The earlier bioRxiv preprint (doi:10.1101/2025.07.08.663126) is superseded by the published article; please cite the article.

Suggested Citation
------------------
Yurtseven A, Joeres R, Kalinina OV. ALPAR: automated learning pipeline for antimicrobial resistance. Bioinformatics 2026; 42(8):btag545. doi:10.1093/bioinformatics/btag545

BibTeX
------
.. code-block:: bibtex

    @article{ALPAR_2026,
      title   = {ALPAR: automated learning pipeline for antimicrobial resistance},
      author  = {Yurtseven, Alper and Joeres, Roman and Kalinina, Olga V.},
      year    = {2026},
      journal = {Bioinformatics},
      volume  = {42},
      number  = {8},
      pages   = {btag545},
      doi     = {10.1093/bioinformatics/btag545},
      url     = {https://doi.org/10.1093/bioinformatics/btag545}
    }

Abstract
=====================
The field of machine learning in antimicrobial resistance (AMR) research has experienced rapid growth, fueled by advancements in high-throughput genome sequencing and growing capacity of computational resources. However, the complexity and lack of standardized data preparation and bioinformatic analyses present significant challenges, especially to newcomers to the domain.

In response to these challenges, we introduce ALPAR (Automated Learning Pipeline for Antimicrobial Resistance), a comprehensive AMR data analysis tool covering the entire process from processing of raw genomic data to training machine learning models to interpretation of results. Our method relies on a reproducible pipeline that integrates widely used bioinformatics tools, presenting a simplified, automatic workflow specifically tailored for single-reference AMR analysis. Accepting genomic data in the form of FASTA files as input, ALPAR facilitates generation of machine learning-ready data tables and both training of machine learning and execution of genome-wide association studies (GWAS) experiments. Additionally, our tool offers supplementary functionalities such as phylogeny-based analysis of distribution of mutations, enhancing its utility for researchers. Our tool is accessible through the ALPAR GitHub page (https://github.com/kalininalab/ALPAR) and installable via conda (https://anaconda.org/kalininalab/ALPAR).
