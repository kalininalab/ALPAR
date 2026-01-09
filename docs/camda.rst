###################
ALPAR on CAMDA (Critical Assessment of Massive Data Analysis)
###################

CAMDA 2025 - Third best model
====================

Using the ALPAR pipeline, we performed analysis of 5346 bacterial strains from nine species given as the training set by the organizers of the `CAMDA anti-microbial resistance prediction challenge 2025 <https://bipress.boku.ac.at/camda2025/>`_, to predict AMR status of 5353 bacterial strains from nine different species (*Acinetobacter baumannii*, *Campylobacter jejuni*, *Escherichia coli*, *Klebsiella pneumoniae*, *Neisseria gonorrhoeae*, *Pseudomonas aeruginosa*, *Salmonella enterica*, *Staphylococcus aureus*, *Streptococcus pneumoniae*) held as the private test set in that challenge.

.. raw:: html

   <style>
     .styled-table {
       border-collapse: collapse;
       margin: 1em 0;
       font-size: 16px;
       font-family: sans-serif;
       width: 100%;
       box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
     }

     .styled-table thead tr {
       background-color: #005aa2;
       color: white;
       text-align: left;
     }

     .styled-table th, .styled-table td {
       padding: 12px 15px;
       border: 1px solid #ddd;
     }

     .styled-table tbody tr:nth-child(even) {
       background-color: #f3f3f3;
     }

     .styled-table em {
       font-style: italic;
     }
   </style>

   <table class="styled-table">
     <thead>
       <tr>
         <th>Bacterium</th>
         <th>MCC Value</th>
       </tr>
     </thead>
     <tbody>
       <tr><td><em>Acinetobacter baumannii</em></td><td>0.772</td></tr>
       <tr><td><em>Campylobacter jejuni</em></td><td>0.931</td></tr>
       <tr><td><em>Escherichia coli</em></td><td>0.463</td></tr>
       <tr><td><em>Klebsiella pneumoniae</em></td><td>0.561</td></tr>
       <tr><td><em>Neisseria gonorrhoeae</em></td><td>0.141</td></tr>
       <tr><td><em>Pseudomonas aeruginosa</em></td><td>0.130</td></tr>
       <tr><td><em>Salmonella enterica</em></td><td>0.817</td></tr>
       <tr><td><em>Staphylococcus aureus</em></td><td>0.948</td></tr>
       <tr><td><em>Streptococcus pneumoniae</em></td><td>0.737</td></tr>
     </tbody>
   </table>

CAMDA 2024 - Winner
====================

Using the ALPAR pipeline, we performed analysis of 5615 bacterial strains from six species given as the training set by the organizers of the `CAMDA anti-microbial resistance prediction challenge 2024 <https://bipress.boku.ac.at/camda-play/the-camda-contest-challenges>`_, to predict AMR status of 1820 bacterial strains from seven different species (*Campylobacter jejuni*, *Campylobacter coli*, *Escherichia coli*, *Klebsiella pneumoniae*, *Neisseria gonorrhoeae*, *Pseudomonas aeruginosa*, *Salmonella enterica*) held as the private test set in that challenge.

.. raw:: html

   <style>
     .styled-table {
       border-collapse: collapse;
       margin: 1em 0;
       font-size: 16px;
       font-family: sans-serif;
       width: 100%;
       box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
     }

     .styled-table thead tr {
       background-color: #005aa2;
       color: white;
       text-align: left;
     }

     .styled-table th, .styled-table td {
       padding: 12px 15px;
       border: 1px solid #ddd;
     }

     .styled-table tbody tr:nth-child(even) {
       background-color: #f3f3f3;
     }

     .styled-table em {
       font-style: italic;
     }
   </style>

   <table class="styled-table">
     <thead>
       <tr>
         <th>Bacterium</th>
         <th>MCC Value</th>
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
