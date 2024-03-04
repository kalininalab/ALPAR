# SANS ambages

**Symmetric Alignment-free phylogeNomic Splits**  
***--- phylogenomics with Abundance-filter, Multi-threading and Bootstrapping on Amino-acid or GEnomic Sequences***

* Reference-free
* Alignment-free
* Input: assembled genomes / reads, or coding sequences / amino acid sequences
* Output: phylogenetic splits or tree
* **NEW:** Abundance-filter
* **NEW:** Bootstrapping
* **NEW:** Multi-threading

### Dos and Don'ts

* The genomes should not be too diverged. SANS works well on species level.
* Be careful with outliers and outgroups (for the reason above).
* The sequences should not be too short. Provide whole-genome data or as many coding sequences as possible.
* Be careful with viruses (for the reasons above).
* Have a look at the network (weakly compatible or 2-tree). It does not make much sense to extract a tree, if the split network is a hairball.
* Reconstructed phylogenies are unrooted, even though a Newick file (-N) suggests a root.
* In case of problems, contact us (see below).


### Publications

Rempel, A., Wittler, R.: [SANS serif: alignment-free, whole-genome based phylogenetic reconstruction](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab444/6300510). Bioinformatics. (2021).

Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs](https://pub.uni-bielefeld.de/download/2942421/2942423/s13015-020-00164-3.wittler.pdf).
Algorithms for Molecular Biology. 15: 4 (2020).

Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs](http://drops.dagstuhl.de/opus/volltexte/2019/11032/pdf/LIPIcs-WABI-2019-2.pdf).
In: Huber, K. and Gusfield, D. (eds.) Proceedings of WABI 2019. LIPIcs. 143, Schloss Dagstuhl--Leibniz-Zentrum fuer Informatik, Dagstuhl, Germany (2019).

## Table of Contents

* [Requirements](#requirements)
* [Compilation](#compilation)
* [Usage](#usage)
* [Examples](#examples)
* [Performance evaluation on predicted open reading frames](#performance-evaluation-on-predicted-open-reading-frames)
* [Clustering / dereplication of metagenome assembled genomes (MAGs)](#clustering-dereplication-of-metagenome-assembled-genomes-mags)
* [Contact](#contact)
* [License](#license)

## Requirements

For the main program, there are no strict dependencies other than C++ version 14.

However, there are some **optional** features:
* To read in a **colored de Bruijn graph**, SANS uses the API of [Bifrost](https://github.com/pmelsted/bifrost).
* To convert the output into NEXUS format, the provided script requires Python 3.
* To visualize the splits, we recommend the tool [SplitsTree](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree).

**Windows** is currently supported with full basic functionaly but limited features only. Confer to branch "windows".

## Compilation

```
git clone https://gitlab.ub.uni-bielefeld.de/gi/sans.git
cd sans
make
```

By default, the installation creates:
* a binary (*SANS*)

You may want to make the binary (*SANS*) accessible via your *PATH* variable.

**Optional:** If Bifrost should be used, change the SANS makefile accordingly (easy to see how). Please note the installation instructions regarding the default maximum *k*-mer size of Bifrost in its README. If during the compilation, the Bifrost library files are not found, make sure that the corresponding folder is found as include path by the C++ compiler. You may have to add `-I/usr/local/include` (with the corresponding folder) to the compiler flags in the makefile. We also recommend to have a look at the [FAQs of Bifrost](https://github.com/pmelsted/bifrost#faq).



## Usage

Use `SANS --help` to obtain a detailed list of options.


**Input files**

Specify your input by `-i <list>` where `<list>` is either a file-of-files or in kmtricks format. Each file can be in fasta, multiple fasta or fastq format.
- **File-of-files:**
  ```
  genome_a.fa
  genome_b.fa
  ...
  ```
  Files can be in subfolders and/or compressed:
  ```
  dataset_1/genome_a.fa.gz
  dataset_1/genome_b.fa.gz
  ...
  ```
  One genome can also be composed of several files (the first one will be used as identifier in the output):
  ```
  reads_a_forward.fa reads_a_reverse.fa
  genome_b_chr_1.fa genome_b_chr_2.fa
  ...
  ```
- **kmtricks format:**
  In this format, you can specify individual identifiers and, optionally, abundance thresholds (see "read data as input"):
  ```
  genome_A : reads_a_forward.fa ; reads_a_reverse.fa ! 2
  genome_B : genome_b_chr_1.fa ; genome_b_chr_2.fa ! 1
  ...
  ```

**Input paramters**

- genomes/assemblies as input: just use `-i <list>`
- read data as input: to filter out *k*-mers of low abundance, either use `-q 2` (or higher thresholds) to specify a global threshold for all input files, or use the kmtricks file-of-files format to specify (individual) thresholds.
- mix of assemblies and read data as input: use the kmtricks file-of-files format to specify individual thresholds.
- coding sequences as input: add `-a` if input is provided as translated sequences, or add `-c` if translation is required. See usage information (`SANS --help`) for further details.


**Output**
- The output file, specified by `-o <split-file>`, lists all splits line by line, sorted by their weight in a tab-separated format where the first column is the split weight and the remaining entries are the identifiers of genomes split from the others.
- For large data sets, the list of splits can become very long. We recommend to restrict the output for *n* genomes as input to the *10n* strongest splits in the output using `-t 10n`.
- We recommend to filter the splits using `-f <filter>`. Then, the sorted list of splits is greedily filtered, i.e., splits are iterated from strongest to weakest and a split is kept if and only if the filter criterion is met.

If you want a **tree**, use `-f strict`. In this case, `-N <newick-file>` can be used to write the resulting tree into a newick file; instead or additionally to `-o <split-file>`.
If you want a **network**, use one of the following filters:
* `weakly`: a split is kept if it is weakly compatible to all previously filtered splits (see publication for definition of "weak compatibility").
* `2-tree`: two sets of compatible splits (=trees) are maintained. A split is added to the first if possible (compatible); if not to the second if possible.
  `3-tree`: three sets of compatible splits (=trees) are maintained. A split is added to the first if possible (compatible); if not to the second if possible; if not to the third if possible.

To visualize the splits, we recommend the tool [SplitsTree](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree). To open the SANS output in SplitsTree, use `scripts/sans2nexus.py <split-file> <list> > <nexus-file>`. In SplitsTree, after opening the `<nexus-file>`, select "Draw" > "EqualAngle" > "Apply". To produce a PDF using SplitsTree on the command line instead, use `scripts/sans2pdf.py <split-file> <list>`.


**Further parameters**
- To observe the progress of SANS during computation, use `-v` to switch to verbose mode.
- You may want to try different values for the *k*-mer length using `-k <integer>`. On shorter sequences, e.g. virus data, use a smaller *k*, e.g., `-k 11`.
- If your input contains 'N's or other ambiguous IUPAC characters, affected *k*-mers are skipped by default. Option `-x <small_integer>` can be used to replace these with the corresponding DNA or AA bases, considering all possibilities.
- By default, all available threads are used for parallel processing. The number of threads can be limited by `-T <integer>`.


**Bootstrapping**
To assess the robustness of reconstructed splits with respect to noise in the input data, bootstrap replicates can be constructed by randomly varying the observed *k*-mer content. To compare the originally determined splits to, e.g., 1000 bootstrap replicates, use `-b 1000`. An additional output file `<split-file>.bootstrap` containing the bootstrap support values will be created. To include them in the nexus file for visualization, use `scripts/sans2conf_nexus.py <split-file> <split-file>.bootstrap <list> > <nexus-file>`.

To generate a consensus tree from bootstrapped trees, use `-f tree -b 1000 -C`. 
To generate a consensus network from bootstrapped trees, use `-f tree -b 1000 -C weakly`. 
See usage information (`SANS --help`) for further options.



## Examples


1. **Split network from assemblies**
   ```
   SANS -i list.txt -o sans.splits -t 10n -f weakly
   scripts/sans2nexus.py sans.splits list.txt > sans.nexus
   ```

   **Tree in newick format from assemblies**
   ```
   SANS -i list.txt -N sans.new -f strict
   ```
   
   **Split network (and tree) from read data**
   ```
   SANS -i list.txt -o sans.splits -t 10n -f weakly -q 2
   scripts/sans2nexus.py sans.splits list.txt > sans.nexus
   (SANS -i list.txt -s sans.splits -f strict -N sans.new)
   ```


2. **Drosophila example data**
   ```
   # go to example directory
   cd example_data/drosophila

   # download data: whole genome (or coding sequences)
   ./download_WG.sh
   (./download_CDS.sh)

   # compute splits
   ../../SANS -i WG_list.txt -o WG_weakly.splits -f weakly -v
   (../../SANS -i CDS_list.txt -o CDS_weakly.splits -f weakly -v -c)

   # generate PDF (if SplitsTree installed)
   ../../scripts/sans2pdf.py WG_weakly.splits WG_list.txt
   (../../scripts/sans2pdf.py CDS_weakly.splits CDS_list.txt)
      
   # filter for tree
   ../../SANS -i WG_list.txt -s WG_weakly.splits -N WG.new -f strict
   (../../SANS -i CDS_list.txt -s CDS_weakly.splits -N CDS.new -f strict)
   
   # generate consensus network from bootstrapped trees
   ../../SANS -i WG_list.txt -f strict -b 1000 -C weakly -o WG_weakly.consensus.splits -v
   ../../scripts/sans2pdf.py WG_consensus.splits WG/list.txt
   (../../SANS -i CDS_list.txt -f strict -b 1000 -C weakly -o CDS_weakly.consensus.splits -v)
   (../../scripts/sans2pdf.py CDS_consensus.splits CDS_list.txt)

   ```
   
   <img src="example_data/drosophila/WG_consensus.splits.png" style="border:0;" alt="Example network"/>

   
3. **Virus example data**
   ```
   # go to example directory
   cd example_data/prasinoviruses

   # download data
   ./download.sh

   # compute splits
   ../../SANS -i list.txt -o weakly.splits -f weakly -k 11 -v 
   
   # generate PDF (if SplitsTree installed)
   ../../scripts/sans2pdf.py weakly.splits list.txt 
   ```

   <img src="example_data/prasinoviruses/weakly.splits.nexus.png" style="border:0;" alt="Example network"/>

   
## Performance evaluation on predicted open reading frames
SANS-serif can predict phylogenies based on amino acid sequences as input. Processing coding sequences is faster than processing whole genome data. Experiments show that the reconstruction quality is comparable. 

If you want to use selected marker genes, the number of genes should be as high as possible to provide sufficient sequence information for extracting phylogenetic signals.

If the genomes at hand are not annotated, you can use a tool to predict open reading frames. The following experiment shows that the reconstruction quality does not suffer from such a simple pre-processing or even improves, while saving total running time, especially because genomes can easily be pre-processed in parallel. 

The following tools have been used with the parameters shown in the table.

| Tool | Parameters| Reference |
|:--|:--|:--|
| SANS | -k (see below) -m geom2 -t 10n -filter strict [-a (if preprocessed)]| |
| Getorf (EMBOSS)| -find (see below) -t 11 | Gary Williams, 2000 |
| ORFfinder | -n true -g 11 | NCBI |
| Prodigal | -q | V2.6.3,Doug Hyatt, 2016 |

For estimating the reconstruction accuracy, the (weighted) F1 score has been determined as follows.

| Measure ||
|--|--|
|F1-score  | harmonic mean of precision and recall |
|precision | (number of called splits that are also in reference) <br> / (total number of called splits) |
|recall    | (number of reference splits that are also called) <br> / (total number of reference splits) |
|weighted precision | (total weight of called splits that are also in reference) <br> / (total weight of all called splits) |
|weighted recall | (total weight of reference splits that are also called) <br> / (total weight of all reference splits) |

Further information on the datasets can be found in the initial publication of SANS (Wittler, 2019), see above.

### *Salmonella enterica* Para C
220 genomes, k=31

| Preprocessing | Running time <br> pre-processing | Running time <br>SANS | Running time both <br>(parall. pre-proc.,<br> factor 16) | F1-Score <br> (weighted) |
|:--|--:|--:|--:|--:|
| none (whole genome) | -- | **610s** | -- | **0.878 <br> (0.999)** |
| Getorf (-find 0) | 247s | 459s | 706s <br> (474s) | 0.881 <br> (0.999) |
| Getorf (-find 1) | 199s | 339s | 538s <br> (351s) | 0.868 <br> (0.999) |
| ORFfinder | 5195s | 166s | 5361s <br> (491s) | 0.858 <br> (0.997) |
| Prodigal | 6292s | 134s | 6326s <br> (521s) |0.853 <br> (0.998)  |

### *Salmonella enterica* subspecies enterica
2964 genomes, k=21

| Preprocessing | Running time <br> pre-processing | Running time <br>SANS | Running time both <br>(parall. pre-proc.,<br> factor 16) | F1-Score <br> (weighted) |
|:--|--:|--:|--:|--:|
| none (whole genome) | -- | **190min** | -- | **0.587 <br> (0.792)** |
| Getorf (-find 0) | 55min | 220min | 274min <br> (223min) | 0.624 <br> (0.807) |
| Getorf (-find 1) | 48min | 165min | 213min <br> (168min) | 0.620 <br> (0.799) |
| ORFfinder | 1621min | 78min | 1698min <br> (179min) | 0.594 <br> (0.766) |
| Prodigal | 1322min | 50min | 1372min <br> (133min) |0.587 <br> (0.792)  |






## Clustering / dereplication of metagenome assembled genomes (MAGs)

For clustering of highly similar sequences, a tree can be constructed which is then chopped into many small subtrees such that the taxa in each subtree correspond to one cluster.
This procedure has been successfully applied for dereplication of metagenome assembled genomes (MAGs). Here, the input are MAGs, and the goal is to cluster these such that the MAGs in each cluster belong to the same strain. 


The general procedure is:

```
# reconstruct a tree
SANS --input <list_of_files> --newick <tree_to_cluster.new> --filter strict --kmer 15 (--verbose) --window W --top T (see below)

# determine clusters from tree
scripts/newick2clusters.py <tree_to_cluster.new> > <clusters.tsv>
```


Due to the usually very high number of input sequences, we recommend the usage of parameters `--window` (`-w`) and `--top` (`-t`) in order to save time and memory. (The experimental parameter `--window` is not mentioned in the usage, because it can lower the accuracy of reconstructed phylogenies considerably. But in this case, the reconstructed tree does not need to be an accurate phylogeny and the parameter has only reasonable effect on the clustering.)

| Setting | Parameters |
|:--|:--|
| quick | --window 25 --top 50n |
| thoroughly | --window 10 --top 100n |


The tree is chopped into clusters as follows:
- Re-root tree to maximum degree node
- In post order traversal:
  - ignore non-branching node (merge edges)
  - get clusters from sub-trees (recursively)
  - if edge longer than parent edge:
    - remove found clusters from current leaf set
    - remaining leaf set =: new cluster
    

### Dereplication efficiency and accuracy

This dereplication approach has been evaluated on a data set from the CAMI challenge [Meyer et al. Critical Assessment of Metagenome Interpretation - the second round of challenges, bioRxiv, 2021, doi: https://doi.org/10.1101/2021.07.12.451567], a simulated mouse gut metagenome representing 64 metagenome samples. Alexander Sczyrba and Peter Belmann filtered the MAGs, provided them for our clustering and compared the clustering to the gold standard.

| Filter | \# MAGs | \# Strains |
|:--|--:|--:|
| MIMAG medium |  5,786 | 686 |
| MIMAG high   |  2,510 | 349 |
| no filter    | 11,602 | 791 |

For the comparison, each called cluster is mapped to a gold standard cluster with maximum intersection, i.e., maximum agreement of contained MAGs. Then, the purity and completeness of the clusters are determined by investigating the number of correct (TP), false (FP) and missing (FN) MAGs in each cluster:

purity := TP / (TP + FP)

completeness := TP / (TP +  FN)

These per-cluster measures were then averaged (weighted and unweighted). The following table shows the results for the different input and clustering settings.

| Input | Setting | Running time | Memory | average purity<br>(weighted) | average completeness<br>(weighted) |
|:--|:--|--:|--:|--:|--:|
| MIMAG medium | quick        | 11h  |  56G | 0.973 <br> (0.956) | 0.884 <br> (0.998)  |
|              | thoroughly   | 50h  | 130G | 0.972 <br> (0.952) | 0.881 <br> (0.998)  |
| MIMAG high   | quick        | 2h   |  16G | 0.983 <br> (0.978) | 0.890 <br> (0.991)  |
|              | thoroughly   | 6h   |  36G | 0.983 <br> (0.979) | 0.912 <br> (0.993)  |
| no filter    | quick        | 59h  | 127G | 0.996 <br> (0.983) | 0.173 <br> (0.668)  |
|              | thoroughly   | 185h | 290G | 0.995 <br> (0.979) | 0.190 <br> (0.700)  |





## Contact

For any question, feedback, or problem, please feel free to file an issue on this Git repository or write an email and we will get back to you as soon as possible.

[sans-service@cebitec.uni-bielefeld.de](mailto:sans-service@cebitec.uni-bielefeld.de)

SANS is provided as a service of the [German Network for Bioinformatics Infrastructure (de.NBI)](https://www.denbi.de/). We would appriciate if you would participate in the evaluation of SANS by completing this [very short survey](https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=sans).


## License

* The sparse-map library is licensed under the [MIT license](https://github.com/Tessil/sparse-map/blob/master/LICENSE).
* The Bifrost library is licensed under the [BSD-2 license](https://github.com/pmelsted/bifrost/blob/master/LICENSE).
* SANS uses gzstream, licensed under the [LGPL license](https://gitlab.ub.uni-bielefeld.de/gi/sans/blob/master/src/gz/COPYING.LIB)
* SANS is licensed under the [GNU general public license](https://gitlab.ub.uni-bielefeld.de/gi/sans/blob/master/LICENSE).

<img src="https://piwik.cebitec.uni-bielefeld.de/matomo.php?idsite=12&rec=1&action_name=VisitGitLab&url=https://gitlab.ub.uni-bielefeld.de/gi/sans" style="border:0;" alt="" />



