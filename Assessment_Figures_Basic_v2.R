#############################################
#Import Alignment Summary File

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
#Percent_Ns_Qs_Aln	Percent_Missing_Data	Seqs_with_above90p_missing_data	Seqs_with_above70p_missing_data	Seqs_with_above50p_missing_data	Seqs_with_less50p_missing_data	
#Seqs_with_no_missing_data


#####################################
#Begin plotting


#FREQUENCY DISTRIBUTIONS
#Edit below numbers to change the frequency distribution set up for your data
#syntax is: seq(low number, high number, increment)
#-----------------------------------------------------------------------------------------

summary(data$Taxa_No)
#edit this seq based on above values!
taxa_seq<- seq(180,265,1)
hist(data$Taxa_No, breaks=taxa_seq, xlab="Number of Taxa", main="Number of Taxa Across Alignments", col="steelblue1")

summary(data$Seq_Length)
#edit this seq based on above values!
seq_seq <- seq(10,1500,20)
hist(data$Seq_Length, breaks=seq_seq, xlab="Sequence Length (bp)", main="Alignment Length Distribution", col="steelblue1")

summary(data$Percent_Gaps_Aln)
gaps_seq <- seq(0,100,1)
hist(data$Percent_Gaps_Aln, breaks=gaps_seq, xlab="Percent", main="Percentage of Gaps Per Alignment", col="steelblue1")

summary(data$Number_Inform_Sites)
#edit this seq based on above values!
inf_sites <- seq(0,1000,5)
hist(data$Number_Inform_Sites, breaks=inf_sites, xlab="Informative Sites", main="Number of Informative Sites", col="steelblue1")

summary(data$Percent_Inform_sites)
perc_inf_sites <- seq(0,100,1)
hist(data$Percent_Inform_sites, breaks=perc_inf_sites, xlab="Percent", main="Percentage of Informative Sites", col="steelblue1")

summary(data$Percent_Ns_Qs_Aln)
Ns_seq <- seq(0,100,1)
hist(data$Percent_Ns_Qs_Aln, breaks=Ns_seq, xlab="Percent", main="Percentage of N's and ?'s Per Alignment", col="steelblue1")

summary(data$Percent_Missing_Data)
missing_seq <- seq(0,100,1)
hist(data$Percent_Missing_Data, breaks=missing_seq, xlab="Percent", main="Percentage of Total Missing Data (-, ?, N) Per Alignment", col="steelblue1")

cols_seq <- seq(0,100,1)
hist(data$Percent_Columns_with_more_20_perc_gaps, breaks=cols_seq, xlab="Percent", main="Percentage of columns in alignments with >20% gaps", col="steelblue1")
hist(data$Percent_Columns_with_more_40_perc_gaps, breaks=cols_seq, xlab="Percent", main="Percentage of columns in alignments with >40% gaps", col="steelblue1")
hist(data$Percent_Columns_with_more_60_perc_gaps, breaks=cols_seq, xlab="Percent", main="Percentage of columns in alignments with >60% gaps", col="steelblue1")

missing_ind <- seq(0,100,1)
hist(data$Seqs_with_above90p_missing_data, breaks=missing_ind, xlab="Number of sequences", main="Sequences with > 90% missing data", col="steelblue1")
hist(data$Seqs_with_above70p_missing_data, breaks=missing_ind, xlab="Number of sequences", main="Sequences with > 70% missing data", col="steelblue1")
hist(data$Seqs_with_above50p_missing_data, breaks=missing_ind, xlab="Number of sequences", main="Sequences with > 50% missing data", col="steelblue1")
hist(data$Seqs_with_less50p_missing_data, breaks=missing_ind, xlab="Number of sequences", main="Sequences with < 50% missing data", col="steelblue1")
plot(data$Seqs_with_no_missing_data, breaks=missing_ind, xlab="Number of sequences", main="Sequences with no missing data", col="steelblue1")


#SCATTER PLOTS
#-----------------------------------------------------------------------------------------

#Alignment length vs. Percentage of Gaps Per Alignment
plot(data$Seq_Length, data$Percent_Gaps_Aln, xlab="Alignment Length", ylab="Percentage of Gaps Per Alignment", pch=16)

#Alignment length vs. Percentage of N's and Q's Per Alignment
plot(data$Seq_Length, data$Percent_Ns_Aln, xlab="Alignment Length", ylab="Percentage of N's and ?'s Per Alignment", pch=16)

#Alignment length vs. Percentage of Total Missing Data (-, ?, N) Per Alignment
plot(data$Seq_Length, data$Percent_Missing_Data, xlab="Alignment Length", ylab="Percentage of Total Missing Data (-, ?, N) Per Alignment", pch=16)

#Alignment length vs. Number of taxa
plot(data$Seq_Length, data$Taxa_No, xlab="Alignment Length", ylab="Number of Taxa", pch=16)

#Number of taxa vs. Alignment length
plot(data$Taxa_No, data$Seq_Length, ylab="Alignment Length", xlab="Number of Taxa",  pch=16)

#Alignment length vs. informative sites
plot(data$Seq_Length, data$Number_Inform_Sites, xlab="Alignment Length", ylab="Informative Sites", pch=16)
evolrate_regression <- lm(data$Number_Inform_Sites ~ data$Seq_Length)
evolrate_regression
summary(evolrate_regression)
abline(evolrate_regression, lwd = 3, col='red')


#############################################
#Automatically save some of the plots using above settings 

#Save individual plots

pdf("Taxon_Numbers.pdf")
hist(data$Taxa_No, breaks=taxa_seq, xlab="Number of Taxa", main="Number of Taxa Across Alignments", col="steelblue1")
dev.off()

pdf("Sequence_Lengths.pdf")
hist(data$Seq_Length, breaks=seq_seq, xlab="Sequence Length (bp)", main="Alignment Length Distribution", col="steelblue1")
dev.off()

pdf("Gaps_Percents.pdf")
hist(data$Percent_Gaps_Aln, breaks=gaps_seq, xlab="Percent", main="Percentage of Gaps Per Alignment", col="steelblue1")
dev.off()

pdf("Informative_Sites_Numbers.pdf")
hist(data$Number_Inform_Sites, breaks=inf_sites, xlab="Informative Sites", main="Number of Informative Sites", col="steelblue1")
dev.off()

pdf("Informative_Sites_Percents.pdf")
hist(data$Percent_Inform_sites, breaks=perc_inf_sites, xlab="Percent", main="Percentage of Informative Sites", col="steelblue1")
dev.off()

pdf("Missing_Data_Percents.pdf")
hist(data$Percent_Missing_Data, breaks=missing_seq, xlab="Percent", main="Percentage of Total Missing Data (-, ?, N) Per Alignment", col="steelblue1")
dev.off()

pdf("Sequence_Lengths_vs_Number_Informative_Sites.pdf")
plot(data$Seq_Length, data$Number_Inform_Sites, xlab="Alignment Length", ylab="Informative Sites", pch=16, abline(evolrate_regression, lwd = 3, col='red'))
dev.off()


#Now save the important plots all in one file
#This file will have 5 plots on one row

pdf("Taxa_SeqLength_MissData_PercInfSites_LvI.pdf", width=35, height=7)
layout(mat = matrix(c(1,2,3,4,5), nrow = 1, ncol = 5), heights = c(2), widths = c(2,2,2,2,2), respect=TRUE)
hist(data$Taxa_No, breaks=taxa_seq, xlab="Number of Taxa", main="Number of Taxa Across Alignments", col="steelblue1", cex.axis=1.6, cex.lab=1.6, cex.main=2)
hist(data$Seq_Length, breaks=seq_seq, xlab="Sequence Length (bp)", main="Alignment Length Distribution", col="steelblue1", cex.axis=1.6, cex.lab=1.6, cex.main=2)
hist(data$Percent_Missing_Data, breaks=missing_seq, xlab="Percent", main="Percentage of Total Missing Data (-, ?, N) Per Alignment", col="steelblue1", cex.axis=1.6, cex.lab=1.6, cex.main=2)
hist(data$Percent_Inform_sites, breaks=perc_inf_sites, xlab="Percent", main="Percentage of Informative Sites", col="steelblue1", cex.axis=1.6, cex.lab=1.6, cex.main=2)
plot(data$Seq_Length, data$Number_Inform_Sites, xlab="Alignment Length", ylab="Informative Sites", pch=16, abline(evolrate_regression, lwd = 3, col='red'), cex.axis=1.6, cex.lab=1.6, cex.main=2)
dev.off()

