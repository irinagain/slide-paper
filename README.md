# slide-paper

R code with functions that were used for simulations and data analysis in the paper "Structural Learning and Integrative Decomposition of Multi-View Data" by Irina Gaynanova and Gen Li


### Core functions used to perform all simulations

**Models_paper.R** - data generation for all models used in simulations
	
**AuxillaryFunctions.R** - standardization functions and rank determination

**SLIDEfunctions.R** - all the functions for SLIDE method that were used in simulations; a stand alone R package is under development.

**JIVEfunctions.R** - wrapper for running JIVE in R via r.jive (both via permutation and with given ranks), also contains wrapper for onestep PCA (first on concatenated dataset, then on residuals) 

**GFAfunctions.R** - wrapper for running GFA in R via CCAGFA


### Functions to run simulations for all methods except AJIVE 

**Sim_d2.R** - two datasets, Section 4.1

**Sim_d2_AC.R** - two datasets, supplement

**Sim_d3.R** - three datasets, Section 4.2

**Sim_d3_noise.R** - three datasets, one dataset is pure noise, supplement

**Sim_d3_perturb.R** - three datasets, low-rank perturbation of globally-shared structure, supplement

**Sim_d3_snr.R** - SLIDE performance as signal to noise ratio changes, supplement


### Functions to run simulations for AJIVE 

For AJIVE, the data from simulations above is saved and stored, and the following matlab code is run:
**Sim_AJIVE.m**. The code relies on using AJIVE as [implemented in MATLAB here](https://github.com/MeileiJiang/AJIVE_Project).


### Functions for the analysis of TCGA BRCA data

We used the TCGA-BRCA data as pre-processed in [Lock &Dunson, 2013](https://academic.oup.com/bioinformatics/article/29/20/2610/276860) with the processing code available from [Eric Lock's website](http://ericfrazerlock.com/software/BCC.zip), the prior link for data access is no longer active, but the original data can still be loaded from [NCI's Genomic Data Commons](https://gdc.cancer.gov/node/877).

Folder **BRCA**:

**BRCA_analysis_d4_JIVE.R** - analysis of BRCA data using JIVE

**BRCA_analysis_d4_SLIDE.R** - analysis of BRCA data using SLIDE

**BRCA_analysis_d4_GFA.R** - analysis of BRCA data using GFA

**BRCA_d4_results.R** - save the output in a more comparable form for later use
