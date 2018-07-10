rm(list=ls())
#############################################
#Alignment summaries

#CHANGE THE PATH TO YOUR FILE ex. "User/Phylip_files/Alignment_Assessment/Master_Alignment_Assessment.txt"
data <- read.delim("FULL PATH TO FILE", header = TRUE, sep = "\t")

#check headers
ls(data)

#quick summary of file:
summary(data)

#show first few lines of table
head(data)

#Header format is:
#Alignment	Taxa_No	Seq_Length	Percent_Gaps_Aln	Number_Inform_Sites	Percent_Inform_sites_Aln	Columns_no_gaps	Columns_with_gap	
#Percent_Columns_with_more_20_perc_gaps	Percent_Columns_with_more_40_perc_gaps	Percent_Columns_with_more_60_perc_gaps	
#Percent_Ns_Qs_Aln	Percent_Missing_Data	Seqs_with_above90%_missing_data	Seqs_with_above70%_missing_data	Seqs_with_above50%_missing_data	Seqs_with_less50%_missing_data	
#Seqs_with_no_missing_data


#summary of variables:
summary(data$Taxa_No)
summary(data$Seq_Length)
summary(data$Percent_Gaps_Aln)
summary(data$Number_Inform_Sites)
summary(data$Percent_Inform_sites)
summary(data$Percent_Ns_Qs_Aln)
summary(data$Percent_Missing_Data)
summary(data$Percent_Columns_with_more_20_perc_gaps)
summary(data$Percent_Columns_with_more_40_perc_gaps)
summary(data$Percent_Columns_with_more_60_perc_gaps)
summary(data$Seqs_with_above90p_missing_data)
summary(data$Seqs_with_above70p_missing_data)
summary(data$Seqs_with_above50p_missing_data)
summary(data$Seqs_with_less50p_missing_data)
summary(data$Seqs_with_no_missing_data)

#####################################
#Begin plotting


#FREQUENCY DISTRIBUTIONS
#Edit below numbers to change the frequency distribution set up for your data!!!!
#syntax is: seq(low number, high number, increment)
#-----------------------------------------------------------------------------------------

taxa_seq<- seq(180,265,1)
hist(data$Taxa_No, breaks=taxa_seq, xlab="Number of Taxa", main="Number of Taxa Across Alignments", col="steelblue1")

seq_seq <- seq(10,1500,20)
hist(data$Seq_Length, breaks=seq_seq, main="Alignment Length Distribution", col="steelblue1")

gaps_seq <- seq(0,100,1)
hist(data$Percent_Gaps_Aln, breaks=gaps_seq, main="Percentage of Gaps Per Alignment", col="steelblue1")

inf_sites <- seq(0,1000,5)
hist(data$Number_Inform_Sites, breaks=inf_sites, main="Number of Informative Sites", col="steelblue1")

perc_inf_sites <- seq(0,100,1)
hist(data$Percent_Inform_sites, breaks=perc_inf_sites, main="Percentage of Informative Sites", col="steelblue1")

Ns_seq <- seq(0,100,1)
hist(data$Percent_Ns_Qs_Aln, breaks=Ns_seq, main="Percentage of N's and ?'s Per Alignment", col="steelblue1")

missing_seq <- seq(0,100,1)
hist(data$Percent_Missing_Data, breaks=missing_seq, main="Percentage of Total Missing Data (Gaps, ?'s, & N's) Per Alignment", col="steelblue1")

cols_seq <- seq(0,100,1)
hist(data$Percent_Columns_with_more_20_perc_gaps, breaks=cols_seq, main="Percentage of columns in alignments with >20% gaps", col="steelblue1")
hist(data$Percent_Columns_with_more_40_perc_gaps, breaks=cols_seq, main="Percentage of columns in alignments with >40% gaps", col="steelblue1")
hist(data$Percent_Columns_with_more_60_perc_gaps, breaks=cols_seq, main="Percentage of columns in alignments with >60% gaps", col="steelblue1")

missing_ind <- seq(0,100,1)
hist(data$Seqs_with_above90p_missing_data, breaks=missing_ind, main="Sequences with > 90% missing data", col="steelblue1")
hist(data$Seqs_with_above70p_missing_data, breaks=missing_ind, main="Sequences with > 70% missing data", col="steelblue1")
hist(data$Seqs_with_above50p_missing_data, breaks=missing_ind, main="Sequences with > 50% missing data", col="steelblue1")
hist(data$Seqs_with_less50p_missing_data, breaks=missing_ind, main="Sequences with < 50% missing data", col="steelblue1")
plot(data$Seqs_with_no_missing_data, breaks=missing_ind, main="Sequences with no missing data", col="steelblue1")


#SCATTER PLOTS
#EDIT THE ylim=c(low number, high number) or xlim=c(low number, high number) TO FIT YOUR DATA SET!!!!
#-----------------------------------------------------------------------------------------

#Alignment length vs. Percentage of Gaps Per Alignment
plot(data$Seq_Length, data$Percent_Gaps_Aln, xlab="Alignment Length", ylab="Percentage of Gaps Per Alignment", ylim=c(0, 50), pch=16)

#Alignment length vs. Percentage of N's Per Alignment
plot(data$Seq_Length, data$Percent_Ns_Aln, xlab="Alignment Length", ylab="Percentage of N's and ?'s Per Alignment", ylim=c(0, 50), pch=16)

#Alignment length vs. Percentage of Total Missing Data (Gaps & N's) Per Alignment
plot(data$Seq_Length, data$Percent_Missing_Data, xlab="Alignment Length", ylab="Percentage of Total Missing Data (Gaps, ?'s, & N's) Per Alignment", ylim=c(0, 50), pch=16)

#Alignment length vs. Number of taxa
plot(data$Seq_Length, data$Taxa_No, xlab="Alignment Length", ylab="Number of Taxa", ylim=c(150, 256), pch=16)

#Number of taxa vs. Alignment length
plot(data$Taxa_No, data$Seq_Length, ylab="Alignment Length", xlab="Number of Taxa", xlim=c(150, 256), pch=16)

#Alignment length vs. informative sites
plot(data$Seq_Length, data$Number_Inform_Sites, xlab="Alignment Length", ylab="Informative Sites", xlim=c(0,1300), ylim=c(0, 700), pch=16)
evolrate_regression <- lm(data$Number_Inform_Sites ~ data$Seq_Length)
evolrate_regression
summary(evolrate_regression)
abline(evolrate_regression, lwd = 3, col='red')

#informative sites vs. Alignment length
plot(data$Number_Inform_Sites, data$Seq_Length, ylab="Alignment Length", xlab="Informative Sites", ylim=c(0,1300), xlim=c(0, 700), pch=16)

#alignment files vs. sequences with missing data
plot(data$Alignment, data$Seqs_with_above90p_missing_data, ylab = "Sequences with > 90% missing data", xlab = "Alignment", pch = 16)
plot(data$Seq_Length, data$Seqs_with_above90p_missing_data, ylab = "Sequences with > 90% missing data", xlab = "Sequence Length", pch = 16)

plot(data$Alignment, data$Seqs_with_above70p_missing_data, ylab = "Sequences with > 70% missing data", xlab = "Alignment", pch = 16)
plot(data$Seq_Length, data$Seqs_with_above70p_missing_data, ylab = "Sequences with > 70% missing data", xlab = "Sequence Length", pch = 16)

plot(data$Alignment, data$Seqs_with_above50p_missing_data, ylab = "Sequences with > 50% missing data", xlab = "Alignment", pch = 16)
plot(data$Seq_Length, data$Seqs_with_above50p_missing_data, ylab = "Sequences with > 50% missing data", xlab = "Sequence Length", pch = 16)

plot(data$Alignment, data$Seqs_with_less50p_missing_data, ylab = "Sequences with < 50% missing data", xlab = "Alignment", pch = 16)
plot(data$Seq_Length, data$Seqs_with_less50p_missing_data, ylab = "Sequences with < 50% missing data", xlab = "Sequence Length", pch = 16)

plot(data$Alignment, data$Seqs_with_no_missing_data, ylab = "Sequences with no missing data", xlab = "Alignment", pch = 16)
plot(data$Seq_Length, data$Seqs_with_no_missing_data, ylab = "Sequences with no missing data", xlab = "Sequence Length", pch = 16)

#############################################