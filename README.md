# Opposing transcription factors MYCL and HEY1 mediate the Notch-dependent airway stem cell fate decision

This repository contains analysis scripts for the datasets presented in "Opposing transcription factors MYCL and HEY1 mediate the Notch-dependent airway stem cell fate decision".

Script and data uploads are currently ongoing during the submission process.

Individual scRNA-seq datasets
  - Each dataset was individually analyzed. R scripts for each dataset can be run using output from Cellranger. Includes SoupX and removal of Scrublet identified doublets.

Integrated scRNA-seq datasets
  - main integration script is called seurat_integration.R. 
  - Tdsort, NS, WtAd3, and Ad3 datasets were integrated to create the dataset depicted in Figure 1A.
  - Gmnc, Mcidas, and NT1 datasets were integrated to create the dataset depicted in Figure 3.
 
CUT&RUN dataset
  - CUT&RUN was performed using the scripts in the autoCutandRun repository at https://github.com/lb15/autoCutandRun. Scripts are derived from CUT&RUNtools at https://bitbucket.org/qzhudfci/cutruntools/ and the CUT&TAG analysis tutorial at https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1.
