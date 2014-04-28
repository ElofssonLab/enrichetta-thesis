enrichetta-thesis
=================
<h3>Conservation folder contains:</h3><br>
<ul>
<li><strong>Conservation_entropy.py</strong> : compute gap-normalized Shannon entropy starting from multiple sequence alignment in Stockholm format</li><br>

<li><strong>properties.py</strong> : analyse the physicochemical properties for each region in the disordered proteins</li><br>
<li><strong>cons_division.py</strong>  : produce two types of division of the results: one divide the data according to the precomputed conservation score (obtained from Conservation_entropy.py); the other divide the data according to the level of conservation of the ordered/disordered regions in the protein (produce 4 groups)</li><br>
</ul>

<h3>Extract_info folder contains:</h3><br>
<ul>
<li><strong>extract_info.py</strong>: for each protein in the dataset extract information (domain,species,ensembl ID, string ID,GO terms) from uniprot txt file</li><br>
<li><strong>features.py</strong> : for each protein extract features like Modified Residues and return a file that combine the localization of the feature and the region where it appear (disorder/order) and the relative entropy of the region</li><br>
</ul>
