###################
ALPAR on CAMDA (Critical Assessment of Massive Data Analysis)
###################

CAMDA 2024 - Winner
====================

Using the ALPAR pipeline, we performed analysis of 5615 bacterial strains from six species given as the training set by the organizers of the `CAMDA anti-microbial resistance prediction challenge 2024 <https://bipress.boku.ac.at/camda-play/the-camda-contest-challenges>`_, to predict AMR status of 1820 bacterial strains from seven different species (*Campylobacter jejuni*, *Campylobacter coli*, *Escherichia coli*, *Klebsiella pneumoniae*, *Neisseria gonorrhoeae*, *Pseudomonas aeruginosa*, *Salmonella enterica*) held as the private test set in that challenge.

.. table:: MCC values of trained models for CAMDA
   :name: MCC_Results_CAMDA

   +-----------------------------+-------------+
   | **Bacterium**              | **MCC Value** |
   +=============================+===============+
   | *Campylobacter jejuni*     | 0.968         |
   +-----------------------------+---------------+
   | *Escherichia coli*         | 0.349         |
   +-----------------------------+---------------+
   | *Klebsiella pneumoniae*    | 0.886         |
   +-----------------------------+---------------+
   | *Neisseria gonorrhoeae*    | 0.935         |
   +-----------------------------+---------------+
   | *Pseudomonas aeruginosa*   | 0.538         |
   +-----------------------------+---------------+
   | *Salmonella enterica*      | 0.706         |
   +-----------------------------+---------------+

We predicted the test set using trained models with the random forest algorithm. All models were trained using the ALPAR Automatix pipeline with the options mentioned above, utilizing species-specific references and protein databases. (*Campylobacter jejuni* model used for both *Campylobacter jejuni* and *Campylobacter coli*). Our predictions achieved a F1-score of 83/100, which was the best performance in the leaderboard.

CAMDA 2025
====================

Work in progress
