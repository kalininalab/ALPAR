CARD:Epi Download README

Use or reproduction of these materials, in whole or in part, by any commercial 
organization whether or not for non-commercial (including research) or commercial purposes
is prohibited, except with written permission of McMaster University. Commercial uses are
offered only pursuant to a written license and user fee. To obtain permission and begin 
the licensing process, see http://card.mcmaster.ca/about.

CITATION:

Edalatmand A. & T.E. Ta et al. CARD:Epi – Contextualizing antimicrobial resistance 
determinants using deep learning language models. bioRxiv, 2026.08.14.744850.


MODELS (download separately from results files):

The AMR-Epi classifier model (Logistic Regression, TF-IDF Word) can be found in the 
amr_epi_classifier directory.

Model files for Named Entity Recognition (NER) and Relationship Extraction (RE) are
organized by lexicon:

config.json.gz
pytorch_model.bin.gz
tokenizer_config.json.gz

Included models reflect best in class:

ARO ner epoch 4 rate 3e-05 batch_size 32 
ENVO ner epoch 4 rate 5e-05 batch_size 64 
ENVO re epoch 4 rate 5e-05 batch_size 8 
FOODON ner epoch 4 rate 5e-05 batch_size 16 
FOODON re epoch 4 rate 5e-05 batch_size 64 
SO ner epoch 4 rate 1e-05 batch_size 16 
SO re epoch 4 rate 5e-05 batch_size 32 
IDO ner epoch 4 rate 1e-05 batch_size 32 
IDO re epoch 4 rate 3e-05 batch_size 16 
UBERON ner epoch 4 rate 1e-05 batch_size 16 
UBERON re epoch 4 rate 1e-05 batch_size 16 
GAZ ner epoch 4 rate 1e-05 batch_size 32 
GAZ re epoch 4 rate 5e-05 batch_size 32 
NCBI_TAXON ner epoch 4 rate 3e-05 batch_size 8 
NCBI_TAXON re epoch 4 rate 5e-05 batch_size 8 


TRAINING AND TESTING DATA are organized by lexicon:

test.dev - testing data
devel.tsv (NER) or dev.tsv (RE) - validation data
train.tsv - training data


RESULTS (download separately from model files):

card_epi_relationship.tsv:

re_id, unique identifier for each relationship between an ARO term and a epidemiological term
cvterm_id, internal CARD identifier
aro_accession, ARO accession
cvterm_name, name of ARO term
epi_lexicon, lexicon of associated term
epi_normalized_accession, accession of associated term
epi_normalized_term_name, associated term
paper_count, number of papers supporting the association

card_epi_relationship_pmid.tsv:

re_pmid_id, unique identifier for relationship : PMID source
re_id, relationship ID in above table
pmid, PubMed identifier
sentence_index, sentence from the Abstract used by BioBert
pub_year, year of the publication

