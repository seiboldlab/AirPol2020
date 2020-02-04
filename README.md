# AirPol2020

Code repository for the manuscript: "Genome-wide analysis reveals mucociliary remodeling of the airway epithelium induced by urban PM2.5". All analyses were performed with R version 3.5.1. There are two major R scripts in this directory:
1. DE_analysis.R: This script contains all the codes necessary to run differential expression analysis between different groups. DESeq2 R package was used to perform differential expression tests.
2. enrichment_analysis.R: This script takes a list of differentially genes as an input and submits those list into Enrichr website through Enrichr API. It will output enriched pathways and GOs along with their significances. This script requires two additional files to work:
  - enrichrAPI.py: a python script that provide the direct call to Enrichr website.
  - EnrichrFunctions.R: a compilation of wrapper functions to call the python function through R
