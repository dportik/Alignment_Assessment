**Assess multiple sequence alignments**
---------------------------------

**Purpose and Overview:**

To perform assessments of any number of multiple sequence alignments contained in a directory.
The alignments must be in phylip format, and should have unique names with a '.phy' or 
'.phylip' extension to be recognized. This tool will quantify several metrics for sequences and
alignments and produce tab-delimited output summary files.


**Overview and Workflow:**

The python script ***Alignment_Assessment_v2.py*** should be used to assess alignments and generate
output files which will be written to a new directory called ***Alignment_Assessment***. One of the 
main output files can then be used to visualize the data and produce plots
with the accompanying R script ***Assessment_Figures_Basic_v2.R***. The python script will perform
assessments on two levels, described below. To use the python script, the python module numpy
must be installed. The usage of the script is very simple, and only requires including the full 
path to a directory of phylip files:

    python Alignment_Assessment_v2.py Volumes/Data/Alignment_Files/


***Sequences within Alignments:***

For all sequences included in a given alignment, the following information will be summarized:

+ Sequence length [Seq_length]
+ Number of gaps [Gap_Count]
+ Percentage of missing data [Percent_Missing_Data]

This information will be written to an output file called ***[Alignment Name]_Assessment.txt***.

***Summaries of Alignments:***

For each alignment the following metrics will be quantified:

+ Number of taxa [Taxa_No]
+ Sequence length [Seq_Length]
+ Percentage of the alignment composed of gaps [Percent_Gaps_Aln]
+ Number of parsimony informative sites in the alignment [Number_Inform_Sites]
+ Percentage of parsimony informative sites in the alignment[Percent_Inform_sites_Aln]
+ Number of alignment columns containing no gaps [Columns_no_gaps]
+ Number of alignment columns containing gaps [Columns_with_gap]
+ Number of alignment columns composed of more than 20% gaps [Percent_Columns_with_more_20_perc_gaps]
+ Number of alignment columns composed of more than 40% gaps [Percent_Columns_with_more_40_perc_gaps]
+ Number of alignment columns composed of more than 60% gaps [Percent_Columns_with_more_60_perc_gaps]
+ Percentage of N or ? characters in the alignment[Percent_Ns_Qs_Aln]
+ Percentage of missing data (includes -, ?, and N) [Percent_Missing_Data]
+ Number of sequences within alignment having more than 90% missing data [Seqs_with_above90%_missing_data]
+ Number of sequences within alignment having more than 70% missing data [Seqs_with_above70p_missing_data]
+ Number of sequences within alignment having more than 50% missing data [Seqs_with_above50p_missing_data]
+ Number of sequences within alignment having below 50% missing data [Seqs_with_less50p_missing_data]
+ Number of sequences within alignment with no missing data [Seqs_with_no_missing_data]

These calculations will be performed for every alignment included and written to an output 
file called ***Master_Alignment_Assessment.txt***.

***Additional Output Files:***

The ***Log_All_Alignments_file.txt*** will simply write the output on the screen to a file
for quick reference. This will contain brief summaries of the alignment characteristics.

An additional file called ***Log_Summary_file.txt*** summarizes the final statistics of the 
combined alignments, which appears on screen at the end of the python script.

***Visualizing Results:***

The R script ***Assessment_Figures_Basic_v2.R*** can be used to create histograms and scatterplots
to visualize the alignment metrics. The script will need to be edited to indicate the path to the
***Master_Alignment_Assessment.txt*** file. The script will also need to be tailored to create histograms
that match the empirical distributions for your data set. This can be done by editing the 'seq' function
to create an appropriate range of bin sizes and overall number of bins for the histograms. The 'seq' function 
takes three arguments: a minimum range value, maximum range value, and increment value. So, seq(0,20,1) will 
created bins ranging from zero to twenty by increments of one. For all non-percentage data, these regular 
sequences will have to be edited to match your empirical values. Luckily, there are only three requiring editing: 
taxa_seq, seq_seq, and inf_sites. Once the plots are checked for any errors, the commands at the bottom of the
R script can be used to write the plots as PDF files to the working directory.



**Citation Information:**

***Using the pipeline.***

The scripts involved with this pipeline were originally published as part of the following work:

+ *Portik, D.M., Smith, L.L., and K. Bi. 2016. An evaluation of transcriptome-based exon capture for frog phylogenomics across multiple scales of divergence (Class: Amphibia, Order: Anura). Molecular Ecology Resources 16: 1069–1083. https://doi.org/10.1111/1755-0998.12541*

If you use or modify this script for your own purposes, please cite this publication.


**Contact:**

Daniel Portik, PhD

Postdoctoral Researcher

University of Arizona

daniel.portik@gmail.com


