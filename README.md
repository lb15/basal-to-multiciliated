# Opposing transcription factors MYCL and HEY1 mediate the Notch-dependent airway stem cell fate decision

This repository contains analysis scripts for the datasets presented in "Opposing transcription factors MYCL and HEY1 mediate the Notch-dependent airway stem cell fate decision".

Script and data uploads are currently ongoing during the submission process.

Individual scRNA-seq datasets
  - Each dataset was individually analyzed. R scripts for each dataset can be run using output from Cellranger.

Integration scRNA-seq datasets
  - Integration was performed using the scripts in the seurat_integration repository at github.com/lb15/seurat_integration
  - Scripts that begin with "sub" are the submission scripts to the UCSF HPC Wynton which uses the SGE scheduler. 

CUT&RUN dataset
  - CUT&RUN was performed using the scripts in the autoCutandRun repository at https://github.com/lb15/autoCutandRun. Scripts are derived mostly from CUT&RUNtools at https://bitbucket.org/qzhudfci/cutruntools/ and the CUT&TAG analysis tutorial at https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1.
