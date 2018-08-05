#
# The script for imputing PLINK produced crc dataset
# and then applying TuRF feature selection method 
# to reduce the size of the dataset and preparing 
# it for applying ML methods
# 

library(data.table)

setwd('~/Dropbox/crc_new')
getwd()

# -------------------------------------------
#
# Read the raw data
#
data.raw <- read.table('clean-GWA-data.raw', header = TRUE)
dim(data.raw)
head(names(data.raw), 10)


# -------------------------------------------
#
# find the number of missing values for each sample
#
Sample_Miss <- apply(data.raw, 1, function(x) length(which(is.na(x))))
head(Sample_Miss)
SS <- data.frame(Samples = row.names(data.raw), IID = data.raw[,"IID"], Miss = Sample_Miss, Class = data.raw[,"PHENOTYPE"])
head(SS)
SSS <- SS[order(SS[,3], decreasing = T),]
head(SSS)
hist(SSS[,3])
length(which(SSS[,4] == 1)) # controls
length(which(SSS[,4] == 2)) # cases


# -------------------------------------------
#
# selecting the cases
# sort the samples based on missing values count ascending order
# choose the top cases as the number of controls
#
ss_cases <- SSS[which(SSS[,4] == 2),]
ss_cases <- ss_cases[order(ss_cases[,3], decreasing = F),]
head(ss_cases)
hist(ss_cases[,3])
ss_cases_good <- ss_cases[1:nrow(ss_controls),]
max(ss_cases_good[,3])


# -------------------------------------------
#
# selecting the controls
# final samples to be used for further works
# create the dataset with best samples
#
ss_controls <- SSS[which(SSS[,4] == 1),]
final_samples <- rbind(ss_controls, ss_cases_good)
good_data <- subset(data.raw, IID%in%final_samples[,"IID"])
dim(good_data)

# -------------------------------------------
#
# going to imputing SNPs, count the missing values in SNPs
# remove snps with above threshold missing values
#
SNP_Miss <- apply(good_data, 2, function(x) length(which(is.na(x))) )
MS <- data.frame(SNP_Miss)
MS["SNP"] <- row.names(MS)
MS <- MS[order(MS[,1], decreasing = T),]
head(MS, n = 10)
dim(MS)
hist(MS[,1], xlab = "Number of missing values", main = "", panel.first=grid()) # Histogram of missing values
abline(v = 10, col = 2, lty = 2) # vertical threshold
write.table(MS, file = "snps_missing_frequency.txt", sep = "\t", quote = F, row.names = F)


# -------------------------------------------
#
# remove SNPs with more than 1% of population missing data
#
SNP_subset <- MS[which(MS[,1] < 10),] # 1%  threshold, 0.01*nrow(data.raw)
head(SNP_subset)
(ncol(good_data) - nrow(SNP_subset))
final <- good_data[, names(good_data)%in%SNP_subset[,2]]
dim(final)
fwrite(final, file = "clean-GWA-data-notimputed.txt", quote = F, sep = "\t", row.names = F)
unique(final$PHENOTYPE)
write.table(SNP_subset[,2], file = "snp_names.txt", sep = "\t", quote = F, row.names = F)


# -------------------------------------------
#
# impute missing data with the most frequent values in each column
#
data.notimpute <- fread("clean-GWA-data-notimputed.txt", data.table = FALSE)
nas <- which(is.na(data.notimpute), arr.ind = T)
head(nas)
dim(nas)
final_imputed <- data.notimpute
final_imputed[] <- lapply(final_imputed, function(x) ifelse(is.na(x), as.numeric(names(which.max(table(x)))), x))
nas <- which(is.na(final_imputed), arr.ind = T)
fwrite(final_imputed, file = "clean-GWA-data-imputed-withpheno.txt", quote = F, sep = "\t", row.names = F )


# -------------------------------------------
#
# remove phenotype columns 
#
unique(final_imputed$PHENOTYPE)
final_imputed_nopheno <- cbind(final_imputed[,7:ncol(final_imputed)], Class=final_imputed[,6])
tail(names(final_imputed_nopheno))
unique(final_imputed_nopheno$Class)
fwrite(final_imputed_nopheno, file = "clean-GWA-data-final.txt", quote = F, sep = "\t", row.names = F )


# -------------------------------------------
#
# Creating dataset for TuRF
# based on the subset of the best features
#
turf <- read.table('turf-log.txt', header = T)
head(turf)
pdf('Rplot_turf_threshold_2k_2.pdf')
hst <- hist(turf[,2], yaxt = "n", ylab="Frequency", xlab = "Score", main = "", panel.first = grid())
axis(2,at=c(1,20000,40000,60000), labels=c(0,expression(paste(2, "x", 10^4)),expression(paste(4, "x", 10^4)),expression(paste(6, "x", 10^4))),las=1)
abline(v = mean(turf[,2]), col = "red", lty = 2, lwd = 2)
abline(v = median(turf[,2]), col = "blue", lty = 1, lwd = 2)
abline(v = mean(turf[,2]) + 3*sd(turf[,2]), col = "green", lty = 3, lwd = 2) 
legend('topright', legend = c('Mean', 'Median', "Threshold"), col = c("red", "blue", "green"), bty = 'n', lty = c(2,1,3), lwd = 2, border = NA)
dev.off()


# -------------------------------------------
#
# Choose subset of features that have importance score
# greater than mean + 3SD
# then select those features from the cleaned and imputed dataset
# and create the final reduced dataset for ML methods purposes
#
turf_sub <- subset(turf, Score > mean(turf[,2]) + 3*sd(turf[,2]))
head(turf_sub)
turf.dset <- final_imputed_nopheno[, names(final_imputed_nopheno)%in%turf_sub[,1]]
turf.dset[, ncol(turf.dset) + 1] <- final_imputed_nopheno$Class - 1
names(turf.dset)[ncol(turf.dset)] <- 'Class'
nas <- which(is.na(turf.dset), arr.ind = T)
head(nas)
write.table(turf.dset, file = 'dataset-2k-turf-gbm.txt', sep = '\t', row.names = F, quote = F)

caret <- turf.dset
dim(caret)
nas <- which(is.na(caret), arr.ind = T)
head(nas)
(len_controls <- length(which(caret[,ncol(caret)] == 0)))
(len_cases <- length(which(caret[,ncol(caret)] == 1)))
caret[which(caret[,ncol(caret)] == 0), ncol(caret)] <- "C1"
caret[which(caret[,ncol(caret)] == 1), ncol(caret)] <- "C2"
write.table(caret, file = 'dataset-2k-turf-caret.txt', sep = '\t', row.names = F, quote = F)

