# NEASE
NEASE  (Network Enrichment method for Alternative Splicing Events) a tool for functional enrichment of alternative splicing exons/events. 


## General info
The python package NEASE (Network-based Enrichment method for Alternative Splicing Events) first detects protein domains affected by AS and uses an integrated protein-protein and domain-domain interaction network to identify protein interaction partners likely affected by AS. 

Next, NEASE performs gene set overrepresentation analysis and identifies enriched pathways based on affected edges. Furthermore, since the statistical approach is network-based, it also prioritizes (differentially) spliced genes and finds new disease biomarkers candidates in case of aberrant splicing.


![alt text](https://user-images.githubusercontent.com/22538496/122232299-6a25cb00-cebb-11eb-8230-b16b6db81b01.png)



## Installation

To install the package from PyPI please run:

`pip install nease` 

To install the package from git:

`git clone https://github.com/louadi/NEASE.git  && cd NEASE`

`pip install .`


Enjoy your instance of NEASE




## Data input

The standard input of the package is a DataFrame object with the exon coordinates and Ensembl IDs of the genes (also recommended).
- first column  - genes IDs (Only Ensembl gene IDs can be used).
- Second column - start of the exon coordinate.
- Third column  - end of the exon coordinate.
- Fourth column - dPSI (optionally)









| Gene              |   Start   |   End     |dpsi  | 
|-------------------|-----------|-----------|------|
| ENSG00000154263   | 69314431  | 69315425  |-0.10 | 
| ENSG00000154265   | 87411893  | 87412033  | 0.13 | 



The package also supports the output of multiple AS differential detection tools such as rMATs, Whippet and MAJIQ.


If you need help with your tool or need to add the support of more tools, please contact us.



## Main functions and examples

Please note, that all functions are annotated with dockstrings for more details.

Import NEASE package and pandas:

```python
import nease
import pandas as pd
```



### Run NEASE 

table: Data input as DataFrame object as explained in "Data input".


input_type: Either "Standard",'Whippet','rmats','Dexeq'or "MAJIQ".

```python
events=nease.run(table, organism='Human',input_type='MAJIQ')
```


###  General functions
Get statstics of your data.

```python
events.get_stats()
```


Get a list of all affected domains.
```python
events.get_domains()
```



Get a list of all affected linear motifs.

```python
events.get_elm()
```



Get a list of all affected residues and their interactions from the PDB.

```python
events.get_pdb()
```


List of affected interactions from domains and linear motif binding.

```python
events.get_edges()
```

###  NEASE enrichment 

Main function of NEASE

database: a list of pathway databases to run enrichment on it. 


```python
# Supported databases:
database=  ['PharmGKB', 'HumanCyc', 'Wikipathways', 'Reactome','KEGG', 'SMPDB', 'Signalink','NetPath', 'EHMN', 'INOH','BioCarta','PID']

# Run enrichment on Reactome only
events.enrich(database=['Reactome'])
```


###  Pathway specific analysis


#### Get list of genes affecting pathways and their statistics
path_id: Pathway ID. You can find pathways id in the enrichment table results.


```python
events.path_analysis('R-HSA-112314')
```


#### Visualize a pathway in the PPI:

Generate an HTML file with the network [example](https://tender-elion-977996.netlify.app/).

path_id: a string representing the Pathway ID.

file: A string representing a local file path for the html file.

k: float  -  A Float for the algorithm  Fruchterman-Reingold force-directed for nodes positions to be tuned by the user. You might need to run the functions multiple times for the optimal visualizations.
        [more details](https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html).
       

```python
events.Vis_path("R-HSA-5674135",file='AS data/enrichment/',k=0.8)
```





## Tutorials


A step-by-step guide to use NEASE is available [here](https://github.com/louadi/NEASE-tutorials).


A simple example for running NEASE on a standard input:
([Notebook](https://github.com/louadi/NEASE-tutorials/blob/main/DCM_analysis.ipynb)/[Google Colab](https://colab.research.google.com/github/louadi/NEASE-tutorials/blob/main/DCM_analysis.ipynb))




## Cite

If you use NEASE, please cite:

Louadi Z, Elkjaer ML, Klug M, et al. Functional enrichment of alternative splicing events with NEASE reveals insights into tissue identity and diseases. bioRxiv; 2021. DOI: 10.1101/2021.07.14.452376.



## Contact 
Zakaria Louadi: louadi@wzw.tum.de
