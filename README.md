# Alignment_Assessment


usage: python alignment_assessment_v1.py [directory with phylip files]

Will sort through all 'NAME.phylip' or 'NAME.phy' files in directory to evaluate
the number of taxa, sequence length, and proportion of missing data.
This will be performed for each alignment, resulting in an output file
called 'Out_NAME_Assessment.txt', which is tab-delimited. 
**NOTE: the alignments must have the '.phylip' extension to be recognized!

Additional outputs:
'Master_Alignment_Assessment.txt', 'Log_Summary_file.txt', 'Log_All_Alignments_file.txt'

The 'Master_Alignment_Assessment.txt' summarizes information for ALL alignments. 
This has an R-script that will need should be edited with correct information to 
create nice plots of very useful information, such as average number of taxa, 
alignment lengths, and proportion of missing data across the alignments.

The 'Log_All_Alignments_file.txt' will simply write the output on the screen to a file
for quick reference. This will contain info about alignment characteristics, informative
sites, and missing data.

One last file is created, called 'Log_Summary_file.txt', which summarizes the statistics of
all the combined alignments (written on screen at the end of the script).

###############
Dependencies:
Numpy (Numerical Python)
###############


------------------------
written for Python 2.7.3
Dan Portik
daniel.portik@berkeley.edu
August 2015
------------------------
