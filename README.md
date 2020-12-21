# IGLOSS algorithm combined with Maximal Clique Heuristic

A research project in bioinformatics made in Python 3 as a part of my master's thesis.

## Folder explanation

### Exact

Works with data received from the IGLOSS algorithm found [here](http://compbioserv.math.hr/igloss/) on the Arabidopsis Thaliana TAIR9 proteome (scale 5).
Contains the Exact.py which uses data from files output_none.tsv and output_ID.tsv to run the Exact Bron-Kerbosch Maximal Clique algorithm and evaluate its effectiveness comparing it with the actual positives found in ATBioPositives.txt.
This code uses clique.txt and temp.txt as auxiliary files.
Final metrics can be found in Ex_metrics.png

### Heuristic

Works with the same data as the Exact algorithm.
Contains the Heuristic.py which uses data from files output_none.tsv and output_ID.tsv to run a Maximal Clique Heuristic algorithm and evaluate its effectiveness comparing it with the actual positives found in ATBioPositives.txt.
Final metrics can be found in He_metrics.png.
