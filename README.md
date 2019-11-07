# Clonal Data Analysis for Stem Cell Transplantation

This repository presents the python code for analyzing the clonal data for stem cell transplantation, in which we try to understand how the clonal output from the first stem cell transplantation predicts the outcome for the second transplantation.  The outcome includes the clonal ID for individual cells, and their cell type (i.e., HSC or Kit+ cell). See our paper for more details. 

There are three components of our code, and they need to be run in the following order to make sure that each notebook has proper input from the upstream computation. Each notebook has its own explanation within the notebook. 

## Generating the clonal annotation 
 - *Combining T1 T2.ipynb*, to combine clonal data from the first transplantation and second transplantation. This can be skipped if a combined dataset already exists, which is the case here. 
 - *clonal_annotation_T1T2_191101.ipynb*, for generating the clonal ID for a given cell. In this step, there is a parameter *dropout*. If it is set to be zero, then there is no clonal barcode dropout, and the notebook outputs data into a folder called *NoDropoutCorrection*, otherwise, the data is generated in a folder called *DropoutCorrection*. 
 
 The clonal annotation code is adapted from the LARRY preprocessing notebook by Caleb Weinreb, but with an improved criterion for clonal clustering.  To run this notebook, there should be a folder *Combined_T1T2* in the directory of this notebook, and the following 3 files that combine the first and second transplantation dataset should be in this folder
- T1T2_cell_bcs_flat.txt: A list of cell barcodes, one barcode name for each cell
- T1T2_samp_id_flat.txt: A list of sample id's, say T1_HSC or T2_Kit, one id for each cell
- T1T2_LARRY_sorted_and_filtered_barcodes.fastq.gz: A fastq file with raw reads, obtained from target sequencing at the clonal barcode regime
In this repository, we have put in *T1T2_cell_bcs_flat.txt* and *T1T2_samp_id_flat.txt*, but not *T1T2_LARRY_sorted_and_filtered_barcodes.fastq.gz*, which is too large. This file can be accessed here: https://www.dropbox.com/s/h31y0fytmsj6zku/T1T2_LARRY_sorted_and_filtered_barcodes.fastq.gz?dl=0

 ## Performing clonal analysis
 Depending on which dataset to use for downstream analysis,  we can use one of the following notebooks for analyzing the clonal data, mostly for generating clonal correlations.
 - *Clonal_data_statistics_summary_NoDropoutCorrection.ipynb* 
 - *Clonal_data_statistics_summary_WithDropoutCorrection.ipynb* 
 These two notebooks are the same, except for different input datasets.  The data and figures will be stored.
 
 ## Computing statistical significance
 Finally, to generate statistical confidence about these clonal correlations, we use a generative graphic model inferred from the data,  as described in the theory supplement in our paper. Depending on the input data, this is implemented in 
 - *Simulate_StemCellDynamics_for_pValue_NoDropoutCorrection.ipynb* 
 - *Simulate_StemCellDynamics_for_pValue_WithDropoutCorrection.ipynb* 
