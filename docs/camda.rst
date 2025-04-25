###################
ALPAR on CAMDA (Critical Assessment of Massive Data Analysis)
###################

CAMDA 2024 - Winner
====================

Using the ALPAR pipeline, we performed analysis of 5615 bacterial strains from six species given as the training set by the organizers of the `CAMDA anti-microbial resistance prediction challenge 2024 <https://bipress.boku.ac.at/camda-play/the-camda-contest-challenges>`_, to predict AMR status of 1820 bacterial strains from seven different species (*Campylobacter jejuni*, *Campylobacter coli*, *Escherichia coli*, *Klebsiella pneumoniae*, *Neisseria gonorrhoeae*, *Pseudomonas aeruginosa*, *Salmonella enterica*) held as the private test set in that challenge.

.. raw:: html

   <table>
     <thead>
       <tr>
         <th><strong>Bacterium</strong></th>
         <th><strong>MCC Value</strong></th>
       </tr>
     </thead>
     <tbody>
       <tr><td><em>Campylobacter jejuni</em></td><td>0.968</td></tr>
       <tr><td><em>Escherichia coli</em></td><td>0.349</td></tr>
       <tr><td><em>Klebsiella pneumoniae</em></td><td>0.886</td></tr>
       <tr><td><em>Neisseria gonorrhoeae</em></td><td>0.935</td></tr>
       <tr><td><em>Pseudomonas aeruginosa</em></td><td>0.538</td></tr>
       <tr><td><em>Salmonella enterica</em></td><td>0.706</td></tr>
     </tbody>
   </table>

We predicted the test set using trained models with the random forest algorithm. All models were trained using the ALPAR Automatix pipeline with the options mentioned above, utilizing species-specific references and protein databases. (*Campylobacter jejuni* model used for both *Campylobacter jejuni* and *Campylobacter coli*). Our predictions achieved a F1-score of 83/100, which was the best performance in the leaderboard.

CAMDA 2025
====================

Work in progress
