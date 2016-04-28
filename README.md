# Alignment_Assessment

usage: python alignment_assessment_v1.py [full path to directory with phylip alignment files]

NOTE:
The alignments are required to be in standard phylip format
The alignments must have the '.phylip' or '.phy' extension to be recognized!


This script will sort through all 'NAME.phylip' or 'NAME.phy' files in a directory to evaluate
the number of taxa, sequence length, informative sites, proportion of missing data, etc.
This will be performed for each alignment, resulting in an output file
called 'Out_NAME_Assessment.txt' for each alignment, which is tab-delimited. 
You can open these files in Excel or R for additional exploration.

Additional file outputs:
'Master_Alignment_Assessment.txt', 'Log_Summary_file.txt', 'Log_All_Alignments_file.txt'

The 'Master_Alignment_Assessment.txt' summarizes information for ALL alignments. 
The R-script included can be opened and edited with corrected information 
(as noted in the script) to create plots of useful information such as 
average number of taxa, alignment lengths, proportion of missing data,
number of informative sites, and comparisons of alignment lengths and 
informative sites.

The 'Log_All_Alignments_file.txt' will simply write the output on the screen to a file
for quick reference. This will contain info about alignment characteristics, informative
sites, and missing data that appears in your terminal window.

One last file is created, called 'Log_Summary_file.txt', which summarizes the statistics of
all the combined alignments (it is written on the terminal screen at the end of the script).

###############
Dependencies:
Numpy (Numerical Python)
###############
------------------------
written for Python 2.7.3
Dan Portik
daniel.portik@berkeley.edu --> daniel.portik@uta.edu
August 2015
------------------------
