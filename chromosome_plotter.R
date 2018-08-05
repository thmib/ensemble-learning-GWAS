#
# Script for building the plot of chromosomes containg SNPs and the coding genes
# Each SNP belong to a chromosome.
# The SNPs could be in coding or non-coding areas.
# 
library(plyr)
setwd('~/Dropbox/crc_new/')

#
# read data
#
data_assoc <- read.table('ML_final_snps_assoc_with_gene.txt', header = F, col.names = c("CHR", "SNP", "A1", "MAF", "P", "Gene"))
data_assoc[which(data_assoc$Gene == "-"),"Gene"] <- NA
head(data_assoc)
summary <- ddply(data_assoc, ~CHR, summarise, Count=length(SNP), Coding=length(which(!is.na(Gene))) )

all_chrs <- c(1:22)
rem <- subset(all_chrs, !(all_chrs %in% unique(data_assoc$CHR)))
rem_df <- data.frame(CHR = rem, Count = rep(0,length(rem)), Coding = rep(0, length(rem)))
final_summary <- rbind(summary, rem_df)
final_summary <- final_summary[order(final_summary$CHR, decreasing = F),]
final_summary$NonCoding <- final_summary$Count - final_summary$Coding
rownames(final_summary) <- final_summary$CHR

#
# Convert the data frame to matrix for the 'barplot' usage
#
cnts <- as.data.frame(t(final_summary))
bar <- cnts[-c(1,2),]
hght <- as.matrix(bar)


pdf('Rplot_snps_chromosome_dist.pdf', width = 9)
mp<- barplot(hght, xlim = range(1,25), las = 1, main="",
             xlab="Chromosome", ylab = "Number of SNPs", 
             col=c("darkblue","yellow"), legend.text = c("Coding", "Non-Coding"), 
             panel.first = grid())
axis(1,at=mp,labels=c(1:22))
dev.off()



#
# Heterogeneous and homogeneous plotter for SNPs
#
final_snps <- read.table("ML_final_snps_assoc_with_gene.txt", header = F, 
                         col.names = c("CHR", "SNP", "A1", "MAF", "Pvalue", "Gene"))
head(final_snps)

data_gbm <- read.table("dataset-2k-turf-gbm.txt", header = T)

data_gbm_case <- subset(data_gbm, Class == 1)

data_gbm_case_sub <- data_gbm_case[, which(colnames(data_gbm_case) %in% paste(final_snps$SNP, final_snps$A1, sep = "_"))]
data_gbm_case_sub$Class <- data_gbm_case$Class

hst1 <- hist(data_gbm_case_sub[,1], plot = T)
hst1$counts
hst1 <- hist(data_gbm_case_sub[,2], plot = F)
hst1$counts

#
# Count the frequency of SNP values among cases
#
counts <- lapply(data_gbm_case_sub[,-ncol(data_gbm_case_sub)], function(x) hist(x, plot = F)$counts[c(1,5,10)])
counts_df <- data.frame(counts)
head(counts_df)


##%%%%%%%%%%%%%%%%%
## Plot 1
## SNPs heterogeneity and homogeneity among cases
##%%%%%%%%%%%%%%%%%
hght <- as.matrix(counts_df[-1,])
xlabs <- names(counts_df)

pdf('Rplot_snps_hetero_homo.pdf', width = 24, height = 10) 
par(mar=c(8,4,4,1)+0.1, mgp=c(3,1,0)) # adjusting margings 
mp <- barplot(hght, xaxt = "n", xlim = range(1,51), las = 1, main="",
             xlab="", ylab = "Number of cases", 
             col=c("darkblue","yellow"), legend.text = c("1", "2"), 
             panel.first = grid())
axis(1, at=mp, labels=F)
text(x=mp, y=-15, labels=xlabs, srt=45, adj=1, xpd=TRUE)
dev.off()


##%%%%%%%%%%%%%%%%%
## Plot 2
## SNPs heterogeneity and homogeneity among cases
##%%%%%%%%%%%%%%%%%
names(counts_df) <- c(1:44)
hght <- as.matrix(counts_df[-1,])
xlabs <- names(counts_df)

pdf('Rplot_snps_hetero_homo_3.pdf', width = 24, height = 10)
mp <- barplot(hght, xaxt = "n", xlim = range(1,50), las = 1, main="", 
              xlab="", ylab = "Number of cases", 
              #space=0.05,
              col=c("darkblue","yellow"), legend.text = c("1", "2"), 
              panel.first = grid())
axis(1, at=mp, labels=F)
text(x=mp, y=-10, 
     labels=xlabs, srt=45, adj=1, xpd=TRUE)
dev.off()


