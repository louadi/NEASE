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



## Tutorials

A step-by-step guide to use nease ia available [here](https://github.com/louadi/NEASE-tutorials).


## Contact 
Zakaria Louadi: louadi@wzw.tum.de
