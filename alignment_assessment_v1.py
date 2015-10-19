import sys
import os
import subprocess as sp
import shutil
import numpy as np
'''
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
'''

fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

out_dir = "Alignment_Assessment"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

#making empty lists
alignment_file_list = []
seq_length_list = []
taxon_number_list = []
gap_proportion_list = []
abs_gap_list = []
informative_sites_list = []
ungapped_informative_sites_list = []
total_missing_list = []

#create output files
out_file = "Master_Alignment_Assessment.txt"
fh_out = open(out_file, 'a')
fh_out.write("Alignment"+'\t'+"Taxa_No"+'\t'+"Seq_Length"+'\t'+"Percent_Gaps_Aln"+'\t'+"Number_Inform_Sites"+'\t'+"Percent_Inform_sites_Aln"+'\t'+"Columns_no_gaps"+'\t'+"Columns_with_gap"+'\t'+"Percent_Columns_with_more_20_perc_gaps"+'\t'+"Percent_Columns_with_more_40_perc_gaps"+'\t'+"Percent_Columns_with_more_60_perc_gaps"+'\t'+"Percent_Ns_Aln"+'\t'+"Percent_Missing_Data"+'\t'+"Seqs_with_above90p_missing_data"+'\t'+"Seqs_with_above70p_missing_data"+'\t'+"Seqs_with_above50p_missing_data"+'\t'+"Seqs_with_less50p_missing_data"+'\t'+"Seqs_with_no_missing_data"+'\n')

all_log = "Log_All_Alignments_file.txt"
fh_all_log = open(all_log,'a')

def percent_calc(x,y):
    percentage = float( ( float(x) / float(y)) * float(100))
    percentage = np.around(percentage, decimals = 1)
    return percentage


#begin loop for all phylip files
for filetype in os.listdir('.'):
    if filetype.endswith('.phylip') or if filetype.endswith('.phy'):
        file_split = filetype.split('.')
        file_name = file_split[0]
        alignment_file_list.append(file_name)
        with open(filetype, 'r') as temp_f:
            print '\n', "-----------------------------------------------------------------------------------------------------"
            print "Beginning analysis of alignment '{}'.".format(filetype)
            print 
            temp_out = "Out_"+file_name+"_Assessment.txt"
            temp_fh = open(temp_out, 'a')
            temp_fh.write("Sample"+'\t'+"Seq_length"+'\t'+"Gap_Count"+'\t'+"Percent_Missing_Data"+'\n')
            #begin running through file by reading lines into memory
            lines = temp_f.readlines()
            #need to distinguish line numbers, starting with 1
            line_count = int(1)
            #start total counts of bp and gaps
            total_bp = int(0)
            total_gaps = int(0)
            total_Ns = int(0)
            missing = int(0)
            seq_list=[]
            ind_above_missing_90 = int(0)
            ind_above_missing_80 = int(0)
            ind_above_missing_70 = int(0)
            ind_above_missing_60 = int(0)
            ind_above_missing_50 = int(0)
            ind_below_missing_50 = int(0)
            ind_no_missing_data = int(0)
            
            #loop through all file lines
            for line in lines:
                clean = line.strip()
                clean = clean.split()
                #find header line which contains taxa number and bp number
                if line_count < 2:
                    taxa_no = int(clean[0])
                    taxon_number_list.append(taxa_no)
                    bp = int(clean[1])
                    seq_length_list.append(bp)
                #all other lines will be name and sequence
                else:
                    gap_count = int(0)
                    taxon = clean[0]
                    seq = clean[1]
                    seq_list.append(seq)
                    seq_bp = int(len(seq))
                    missing_ind = int(0)
                    #go through every bp in this sequence:
                    for ind_bp in seq:
                        if ind_bp == '-':
                            gap_count += 1
                            total_gaps += 1
                            missing += 1
                            missing_ind += 1
                            total_bp += 1
                        elif ind_bp == 'N':
                            total_Ns += 1
                            total_bp += 1
                            missing += 1
                            missing_ind += 1
                        else:
                            total_bp += 1
                    #individual missing data perc
                    percent_missing = float(( float(missing_ind) / float(seq_bp)) * float(100))
                    if percent_missing >= float(90):
                        ind_above_missing_90 += 1
                    if percent_missing >= float(80):
                        ind_above_missing_80 += 1
                    if percent_missing >= float(70):
                        ind_above_missing_70 += 1
                    if percent_missing >= float(60):
                        ind_above_missing_60 += 1
                    if percent_missing >= float(50):
                        ind_above_missing_50 += 1
                    if percent_missing < float(50):
                        ind_below_missing_50 += 1
                    if percent_missing == float(0):
                        ind_no_missing_data += 1

                    #convert variable to strings for alignment specific output file writing
                    seq_bp = str(seq_bp)
                    gap_count = str(gap_count)
                    percent_missing = str(percent_missing)
                    temp_fh.write(taxon+'\t'+seq_bp+'\t'+gap_count+'\t'+percent_missing+'\n')
                #add one to the line count
                line_count+=1
                
            #Start all counters
            informative_sites = int(0)
            ungapped_informative_sites = int(0)
            gap_free_columns = int(0)
            gapped_columns = int(0)
            p20p_gaps = int(0)
            p40p_gaps = int(0)
            p60p_gaps = int(0)
            
            #transpose sequence list into matrix style (functionally turn into alignment columns)
            seq_matrix = zip(*seq_list)
            
            #iterate through each alignment column
            for aln_colmn in seq_matrix:
                #start all base pair counts at zero
                A_count = int(0)
                T_count = int(0)
                G_count = int(0)
                C_count = int(0)
                t_bases = int(0)
                N_count = int(0)
                gap_count = int(0)
                #count bases
                for base in aln_colmn:
                    base = base.upper()
                    t_bases += 1
                    if base == 'A':
                        A_count+=1
                    elif base == 'T':
                        T_count+=1
                    elif base == 'G':
                        G_count+=1
                    elif base == 'C':
                        C_count+=1
                    elif base == '-':
                        gap_count+=1
                    elif base == 'N':
                        N_count+=1
                    #Consider how to treat ambiguities (split into haplos)
                    elif base =='R':
                        A_count+=1
                        G_count+=1
                    elif base =='Y':
                        C_count+=1
                        T_count+=1
                    elif base =='S':
                        G_count+=1
                        C_count+=1
                    elif base =='W':
                        A_count+=1
                        T_count+=1
                    elif base =='K':
                        G_count+=1
                        T_count+=1
                    elif base =='M':
                        C_count+=1
                        A_count+=1
                        
                #series of statements for assessing alignment column
                #if more than 2 base pair types in a column, count as an informative site (ignoring gaps)
                if ( (A_count >= int(2) and T_count >= int(2))
                     or (A_count >= int(2) and G_count >= int(2))
                     or (A_count >= int(2) and C_count >= int(2))
                     or (T_count >= int(2) and G_count >= int(2))
                     or (T_count >= int(2) and C_count >= int(2))
                     or (G_count >= int(2) and C_count >= int(2)) ):
                    informative_sites += 1
                #do the same, but only if the column has zero gaps
                if ( gap_count == int(0) and ((A_count >= int(2) and T_count >= int(2))
                     or (A_count >= int(2) and G_count >= int(2))
                     or (A_count >= int(2) and C_count >= int(2))
                     or (T_count >= int(2) and G_count >= int(2))
                     or (T_count >= int(2) and C_count >= int(2))
                     or (G_count >= int(2) and C_count >= int(2))) ):
                    ungapped_informative_sites += 1
                #more assessments about number of gapped columns
                if gap_count == int(0):
                    gap_free_columns+=1
                if gap_count > int(0):
                    gapped_columns+=1
                #assess the gaps of columns with percentage thresholds
                perc_gapped_columns = float( ( float(gap_count) / float(bp)) * float(100))
                if perc_gapped_columns >= float(20):
                    p20p_gaps += 1
                if perc_gapped_columns >= float(40):
                    p40p_gaps += 1
                if perc_gapped_columns >= float(60):
                    p60p_gaps += 1
                    
            #calculate percentages using function percent_calc
            perc_inf_sites = percent_calc(informative_sites, bp)
            informative_sites_list.append(perc_inf_sites)
            perc_ungap_inf_sites = percent_calc(ungapped_informative_sites, bp)
            ungapped_informative_sites_list.append(perc_ungap_inf_sites)
            proportion_20 = percent_calc(p20p_gaps, bp)
            proportion_40 = percent_calc(p40p_gaps, bp)
            proportion_60 = percent_calc(p60p_gaps, bp)
            abs_missing_data = percent_calc(total_gaps, total_bp)
            abs_gap_list.append(abs_missing_data)
            abs_Ns = float( ( float(total_Ns) / float(total_bp)) * float(100))
            abs_Ns = np.around(abs_Ns, decimals = 2)
            abs_missing = percent_calc(total_gaps, total_bp)
            total_missing = percent_calc(missing, total_bp)
            total_missing_list.append(total_missing)
            
            #convert variables to strings
            ungapped_informative_sites_count = str(ungapped_informative_sites)
            perc_ungap_inf_sites = str(perc_ungap_inf_sites)
            informative_sites_count = str(informative_sites)
            perc_inf_sites = str(perc_inf_sites)
            total_Ns = str(total_Ns)
            abs_Ns = str(abs_Ns)
            bp = str(bp)
            taxa_no = str(taxa_no)
            total_bp = str(total_bp)
            total_gaps = str(total_gaps)
            abs_missing_data = str(abs_missing_data)
            gap_free_columns = str(gap_free_columns)
            gapped_columns = str(gapped_columns)
            proportion_20 = str(proportion_20)
            proportion_40 = str(proportion_40)
            proportion_60 = str(proportion_60)
            total_missing = str(total_missing)
            ind_above_missing_90 = str(ind_above_missing_90)
            ind_above_missing_70 = str(ind_above_missing_70)
            ind_above_missing_50 = str(ind_above_missing_50)
            ind_below_missing_50 = str(ind_below_missing_50)
            ind_no_missing_data = str(ind_no_missing_data)
            
            tracker_string = str("Alignment"+'\t'+"Taxa_No"+'\t'+"Seq_Length"+'\t'+"Percent_Gaps_Aln"+'\t'+"Number_Inform_Sites"+'\t'+"Percent_Inform_sites_Aln"+'\t'+"Columns_no_gaps"+'\t'+"Columns_with_gap"+'\t'+"Percent_Columns_with_more_20_perc_gaps"+'\t'+"Percent_Columns_with_more_40_perc_gaps"+'\t'+"Percent_Columns_with_more_60_perc_gaps"+'\t'+"Percent_Ns_Aln"+'\t'+"Percent_Missing_Data"+'\t'+"Seqs_with_above90p_missing_data"+'\t'+"Seqs_with_above70p_missing_data"+'\t'+"Seqs_with_above50p_missing_data"+'\t'+"Seqs_with_less50p_missing_data"+'\t'+"Seqs_with_no_missing_data"+'\n')
            
            fh_out.write(file_name+'\t'+taxa_no+'\t'+bp+'\t'+abs_missing_data+'\t'+informative_sites_count+'\t'+perc_inf_sites+'\t'+gap_free_columns+'\t'+gapped_columns+'\t'+proportion_20+'\t'+proportion_40+'\t'+proportion_60+'\t'+abs_Ns+'\t'+total_missing+'\t'+ind_above_missing_90+'\t'+ind_above_missing_70+'\t'+ind_above_missing_50+'\t'+ind_below_missing_50+'\t'+ind_no_missing_data+'\n')

            print "1. Alignment Characteristics:", '\n'
            print "Number of taxa = {}.".format(taxa_no)
            print "Sequence length = {}.".format(bp), '\n'
            print "2. Informative Sites metrics:", '\n'
            print "Number of informative sites = {}.".format(informative_sites_count)
            print "Percentage of informative sites (gaps ignored) = {}%.".format(perc_inf_sites), '\n'
            print "3. Missing Data metrics:", '\n'
            print "Total number of alignment bp = {}.".format(total_bp), '\n'
            print "Total number of N's in alignment = {}.".format(total_Ns)
            print "Absolute proportion of Ns in alignment = {}%.".format(abs_Ns), '\n'
            print "Total number of alignment gaps = {}.".format(total_gaps)
            print "Absolute proportion of gaps across alignment = {}%.".format(abs_missing_data), '\n'
            print "Number of gap-free alignment columns = {}".format(gap_free_columns)
            print "Number of alignment columns containing gaps = {}".format(gapped_columns)
            print "Percentage of alignment columns containing > 20% gaps = {}%".format(proportion_20)
            print "Percentage of alignment columns containing > 40% gaps = {}%".format(proportion_40)
            print "Percentage of alignment columns containing > 60% gaps = {}%".format(proportion_60), '\n'
            print "Absolute proportion of missing data in alignment (N's and gaps) = {}%.".format(total_missing)
            print '\n', "Sequences with greater than 90% missing data = {}.".format(ind_above_missing_90)
            print "Sequences with greater than 70% missing data = {}.".format(ind_above_missing_50)
            print "Sequences with greater than 50% missing data = {}.".format(ind_above_missing_70)
            print "Sequences with less than 50% missing data = {}.".format(ind_below_missing_50)
            print "Sequences with no missing data = {}.".format(ind_no_missing_data)
            print "-----------------------------------------------------------------------------------------------------", '\n'

            #write messy on screen log to an output file
            fh_all_log.write("-----------------------------------------------------------------------------------------------------"+'\n'+"Analysis of alignment: '{}'.".format(filetype)+'\n'+'\n'+"1. Alignment Characteristics:"+'\n'+"Number of taxa = {}.".format(taxa_no)+'\n'+"Sequence length = {}.".format(bp)+'\n'+'\n'+"2. Informative Sites metrics:"+'\n'+'\n'+"Number of informative sites = {}.".format(informative_sites_count)+'\n'+"Percentage of informative sites (gaps ignored) = {}%.".format(perc_inf_sites)+'\n'+'\n'+"3. Missing Data metrics:"+'\n'+'\n'+"Total number of alignment bp = {}.".format(total_bp)+'\n'+'\n'+"Total number of N's in alignment = {}.".format(total_Ns)+'\n'+"Absolute proportion of Ns in alignment = {}%.".format(abs_Ns)+'\n'+'\n'+"Total number of alignment gaps = {}.".format(total_gaps)+'\n'+"Absolute proportion of gaps across alignment = {}%.".format(abs_missing_data)+'\n'+'\n'+"Number of gap-free alignment columns = {}".format(gap_free_columns)+'\n'+"Number of alignment columns containing gaps = {}".format(gapped_columns)+'\n'+"Percentage of alignment columns containing > 20% gaps = {}%".format(proportion_20)+'\n'+"Percentage of alignment columns containing > 40% gaps = {}%".format(proportion_40)+'\n'+"Percentage of alignment columns containing > 60% gaps = {}%".format(proportion_60)+'\n'+'\n'+"Absolute proportion of missing data in alignment (N's and gaps) = {}%.".format(total_missing)+'\t'+"Sequences with greater than 90% missing data = .".format(ind_above_missing_90)+'\t'+"Sequences with greater than 50% missing data = .".format(ind_above_missing_50)+'\t'+"Sequences with less than 50% missing data = {}.".format(ind_below_missing_50)+'\t'+"Sequences with no missing data = .".format(ind_no_missing_data)+'\n'+"-----------------------------------------------------------------------------------------------------"+'\n'+'\n'+'\n')

            #close temp file, move to out_dir
            temp_fh.close()
            
            for filetype2 in os.listdir('.'):
                if filetype2.startswith('Out_'):
                    shutil.move(filetype2, out_dir)
            #continues loop

fh_out.close()          
shutil.move(out_file, out_dir)

#create function for making averages across our bigger lists
def averaging_function(x):
    x_array = np.asarray(x,dtype=np.float64)
    x_avg = np.average(x_array)
    x_avg = np.around(x_avg, decimals = 1)
    avg_string = str(x_avg)
    return avg_string

#give a name to the output of the function (string of the average) and execute function for each list 
seq_length_avg = averaging_function(seq_length_list)
taxon_number_avg = averaging_function(taxon_number_list)
abs_gap_avg = averaging_function(abs_gap_list)
inf_sites_avg = averaging_function(informative_sites_list)
ung_inf_sites_avg = averaging_function(ungapped_informative_sites_list)
total_missing_avg = averaging_function(total_missing_list)
#count files read
number_files = str(len(alignment_file_list))
#create final summary output file
out_file2 = "Log_Summary_file.txt"
fh_out2 = open(out_file2, 'a')
fh_out2.write("Analyzed "+number_files+" alignment files."+'\n'+"Average number of taxa in alignments = {}.".format(taxon_number_avg)+'\n'+"Average number of base pairs across alignments = {}.".format(seq_length_avg)+'\n'+"Average percent of gaps per alignment = {}%.".format(abs_gap_avg)+'\n'+"Average percent of informative sites per alignment (gaps ignored) = {}%.".format(inf_sites_avg)+'\n'+"Average percent of missing data (gaps and N's) across alignments = {}%.".format(total_missing_avg))

#close files, move them
fh_all_log.close()
fh_out2.close()
shutil.move(out_file2, out_dir)
shutil.move(all_log, out_dir)

#print out info in the final summary output file
print '\n', "-----------------------------------------------------------------------------------------------------"                   
print "Finished analysis of", len(alignment_file_list), "alignment files."
print "Average number of taxa in alignments = {}.".format(taxon_number_avg)
print "Average number of base pairs across alignments = {}.".format(seq_length_avg)
print "Average percent of gaps per alignment = {}%.".format(abs_gap_avg)
print "Average percent of informative sites per alignment (gaps ignored) = {}%.".format(inf_sites_avg)
print "Average percent of missing data (gaps and N's) across alignments = {}%.".format(total_missing_avg)
print "-----------------------------------------------------------------------------------------------------", '\n'
print "#####################################################################################################"
print "Check directory: 'Alignment_Assessment/' which contains the following useful files: 'Master_Alignment_Assessment.txt', 'Log_Summary_file.txt', 'Log_All_Alignments_file.txt', and output files for each alignment.", '\n'
print "The included R script can be used to plot information from 'Master_Alignment_Assessment.txt'."
print "#####################################################################################################", '\n'
