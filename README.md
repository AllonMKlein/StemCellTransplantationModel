# StemCellTransplantationModel

This repository presents the python code for analyzing the clonal data for stem cell transplantation, in which we try to understand how the clonal output from the first stem cell transplantation predicts the outcome for the second transplantation.  The outcome includes the clonal ID for individual cells, and their cell type (i.e., HSC or LSK). See our paper for more details. 

The repository contains the following notebooks
 - *Combining T1 T2.ipynb*, to combine clonal data from the first transplantation and second transplantation. 
 - *clonal_annotation_T1T2_191101.ipynb*, for generating the clonal ID for a given cell. In this step, there is a parameter *dropout*. If it is set to be zero, then there is no dropout, and the code outputs data into a folder called *NoDropoutCorrection*, otherwise, the data is generated in a folder called *DropoutCorrection*.  To make this notebook to work, we assume that 
 
 Depending which dataset to use for downstream analysis,  we can one of the following notebook for analyzing the clonal data, mostly for generating clonal correlations.
 - *Clonal_data_statistics_summary_NoDropoutCorrection.ipynb* 
 - *Clonal_data_statistics_summary_WithDropoutCorrection.ipynb* 
 These two notebooks are the same, except different input dataset.  
 
 Finally, to generate statistical confidence about these clonal correlations, we use a generative model that is described in the theory supplement in our paper. Depending on the input data, this is implemented in 
 - *Simulate_StemCellDynamics_for_pValue_NoDropoutCorrection.ipynb* 
 - *Simulate_StemCellDynamics_for_pValue_WithDropoutCorrection.ipynb* 
