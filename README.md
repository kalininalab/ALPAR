# SR-AMR
Pipeline of Single Reference Antimicrobial Resistance Project

Single Reference Pipeline
* Rename files with random names (Otherwise pandas creates a problem) -> Save them in text file
* Run snippy on genomic fasta files -> Collect vcf files
* Run prokka on genomic fasta files -> Collect gff files
* Run panaroo on gff files that are created with prokka -> Collect gene presence absence file
* Create binary table from mutations
* Add gene presence absence information to binary table
* From phenotype information -> Create phenotype table
* Outputs : Genotype -> tsv file, Phenotype -> tsv file
- GWAS (Optional) 
* Preapare GWAS input files from Genotype & Phenotype tables
* Run GWAS with pyseer
* Sort GWAS results with column lrt-pvalue
* (Optional) Apply bonferoni correction (# of features / 0.05)
