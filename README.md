# StemCellTransplantationModel

This repository presents the python code for analyzing the clonal data for stem cell transplantation, in which we try to understand how the clonal output from the first stem cell transplantation predicts the outcome for the second transplantation.  The outcome includes the clonal ID for individual cells, and their cell type (i.e., HSC or LSK). See our paper for more details. 

There are three components of our code, and they need to be run in order to make sure that each notebook has proper input from upstream. Each notebook has its own explanation within the notebook. 

## Generating the clonal annotation 
 - *Combining T1 T2.ipynb*, to combine clonal data from the first transplantation and second transplantation. This can be skipped if a combined dataset alreay exists. 
 - *clonal_annotation_T1T2_191101.ipynb*, for generating the clonal ID for a given cell. In this step, there is a parameter *dropout*. If it is set to be zero, then there is no dropout, and the code outputs data into a folder called *NoDropoutCorrection*, otherwise, the data is generated in a folder called *DropoutCorrection*. 
 
 ...This code is adapted from the LARRY preprocessing notebook by Caleb Weinreb. To run this notebook, there should be a folder ...Combined_T1T2 in the directory of this notebook, and the following 3 files should be in this folder
 ... -T1T2_cell_bcs_flat.txt: A list of cell barcodes, one barcode name for each cell
... -T1T2_samp_id_flat.txt: A list of sample id's, say T1_HSC or T2_Kit, one id for each cell
...- T1T2_LARRY_sorted_and_filtered_barcodes.fastq.gz: A fastq file with raw reads, obtained from target sequencing at the clonal barcode regime
 
 ## Performing clonal analysis
 Depending which dataset to use for downstream analysis,  we can one of the following notebook for analyzing the clonal data, mostly for generating clonal correlations.
 - *Clonal_data_statistics_summary_NoDropoutCorrection.ipynb* 
 - *Clonal_data_statistics_summary_WithDropoutCorrection.ipynb* 
 These two notebooks are the same, except different input dataset.  
 
 ## Computing statistical significance
 Finally, to generate statistical confidence about these clonal correlations, we use a generative model that is described in the theory supplement in our paper. Depending on the input data, this is implemented in 
 - *Simulate_StemCellDynamics_for_pValue_NoDropoutCorrection.ipynb* 
 - *Simulate_StemCellDynamics_for_pValue_WithDropoutCorrection.ipynb* 

