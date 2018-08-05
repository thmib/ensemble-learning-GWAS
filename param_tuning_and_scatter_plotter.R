#
# This script reads the scores produced by ML methods for the SNPs 
# to get the most significant ones, we draw the scatter plot of GBM and RF scores
# the x-axis is GBM score, the y-axis is RF score. 
# the top-rightmost SNPs are the most significant ones by both methods
#

setwd('~/Dropbox/crc_new')
getwd()

# -------------------------------------------
#
# Reading the RF scores
#
rf_scores <- read.table('RF_snps_importance_frequency_N2000_M100.txt', header = T, col.names = c("Name", "RF_importance", "RF_count"))
head(rf_scores)


# -------------------------------------------
#
# Histogram of SNPs cumulative score by RF
#
h = hist(x,main = "")
pdf('Rplot_RF_average_score_hist.pdf')
plot(h,freq=FALSE,las=1, xaxt = "n", main = "", xlab = "Average score", ylab = "Density", panel.first=grid())
axis(1, at = c(10,20,30,40,50,60,70), labels = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
dev.off()


# -------------------------------------------
#
# Reading the GBM scores
#
gbm_scores <- read.table('GBM_snps_importance_frequency_N2000_I10_S1.txt', header = T, col.names = c("Name", "GBM_importance", "GBM_count"))
head(gbm_scores)

# -------------------------------------------
#
# Histogram of SNPs cumulative score by GBM
#
h = hist(gbm_scores$GBM_importance, main = "")
pdf('Rplot_GBM_average_score_hist.pdf')
plot(h,freq=FALSE,las=1, main = "", xaxt = "n", xlab = "Average score", ylab = "Density", panel.first=grid())
axis(1, at = c(0,200,400,600), labels = c(0,2,4,6))
dev.off()


# -------------------------------------------
#
# Merging the RF and GBM score files based on SNP name
#
merged <- merge(rf_scores, gbm_scores, "Name")
merged <- merged[order(merged[,2], decreasing = TRUE),]
head(merged)
write.table(merged, file = "ML_Scores_merged_all.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
#
# Drawing the scatter plot of ML scores
#
RF_mean <- max(merged[,2]) / 3
GBM_mean <- max(merged[,4]) / 3
pdf('Rplot_ml_scores_scatter.pdf')
plot(x = merged[,4]/100, y = merged[,2]/100, las = 1, xlab = "GBM importance score", ylab = "RF importance score", main = "", panel.first=grid()) #Scatter plot of SNPs scores by ML methods
abline(a = 80/100, b = -0.24, col = 2, lty = 2)
dev.off()


# -------------------------------------------
#
# Creating the final SNPs subset 
#
final_subset <- subset(merged, RF_importance > -0.24*GBM_importance + 80, select = -c(3,5))
dim(final_subset)
final_subset <- final_subset[order(final_subset[,2], decreasing = TRUE),]
final_subset[,2] <- round(final_subset[,2], 2)
final_subset[,3] <- round(final_subset[,3], 2)
head(final_subset)
write.table(final_subset, file = "ML_Final_snps_by_MLs.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
# 
# Reading PLINK association files
#
assoc <- read.table('7assoc.assoc', header = T)
head(assoc)
assoc_sub <- subset(assoc, paste(assoc$SNP, assoc$A1, sep = "_")%in%final_subset$Name)
head(assoc_sub)
write.table(subset(assoc_sub, select = c(1,2,4,5,9)), file = "ML_final_snps_assoc.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
#
# Creating graphs for the results of parameter tuning tables
# GBM and RF parameter tuning plots are produced for every combiantion
# and put in the same plot for each method
#
Col1 <- rgb(219, 30, 77, maxColorValue = 255) 
Col2 <- rgb(12, 249, 238, maxColorValue = 255)
Col3 <- rgb(30, 90, 219, maxColorValue = 255) 
Col4 <- rgb(219, 121, 30, maxColorValue = 255)
Col5 <- rgb(77, 219, 30, maxColorValue = 255) 
COL <- c(Col1, Col2, Col3, Col4, Col5)
PCH <- c(2,15,19,7,8,3)


# -------------------------------------------
# Tuning plot for RF
#
rf_params <- read.table('RF_params_table.txt', header = T)
head(rf_params)

pdf('Rplot_RF_params_tuning_3.pdf', width = 12)
par(mfrow=c(1,2))
mtry1 <- subset(rf_params, mtry == 100)
plot(cex.lab=1.2, x = mtry1$ntree, y = mtry1$mean, xlim = range(rf_params$ntree), ylim = range(rf_params$mean), pch = PCH[1], col = COL[1], las = 1, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "a)", xlab = "ntree", ylab = "Average accuracy")

mtry2 <- subset(rf_params, mtry == 200)
points(x = mtry2$ntree, y = mtry2$mean, xlim = range(rf_params$ntree), ylim = range(rf_params$mean), pch = PCH[2], col = COL[2], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "", xlab = "ntree", ylab = "Average accuracy")

mtry2 <- subset(rf_params, mtry == 300)
i <- 3
points(x = mtry2$ntree, y = mtry2$mean, xlim = range(rf_params$ntree), ylim = range(rf_params$mean), pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "", xlab = "ntree", ylab = "Average accuracy")

mtry2 <- subset(rf_params, mtry == 500)
i <- 4
points(x = mtry2$ntree, y = mtry2$mean, xlim = range(rf_params$ntree), ylim = range(rf_params$mean), pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "", xlab = "ntree", ylab = "Average accuracy")

mtry2 <- subset(rf_params, mtry == 1000)
i <- 5
points(x = mtry2$ntree, y = mtry2$mean, xlim = range(rf_params$ntree), ylim = range(rf_params$mean), pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "", xlab = "ntree", ylab = "Average accuracy")
P <- c('100', '200', '300', '500', '1000')
legend("topleft", legend = P, col = COL, bg = "white", horiz=F, pch=PCH, pt.lwd = 1, title = "mtry")

# Plot for AUC
mtry1 <- subset(rf_params, mtry == 100)
plot(cex.lab=1.2, x = mtry1$ntree, y = mtry1$auc, xlim = range(rf_params$ntree), ylim = range(rf_params$auc), pch = PCH[1], col = COL[1], las = 1, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "b)", xlab = "ntree", ylab = "Average AUC")

mtry2 <- subset(rf_params, mtry == 200)
points(x = mtry2$ntree, y = mtry2$auc, xlim = range(rf_params$ntree), ylim = range(rf_params$auc), pch = PCH[2], col = COL[2], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "", xlab = "ntree", ylab = "Average AUC")

mtry2 <- subset(rf_params, mtry == 300)
i <- 3
points(x = mtry2$ntree, y = mtry2$auc, xlim = range(rf_params$ntree), ylim = range(rf_params$auc), pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "", xlab = "ntree", ylab = "Average AUC")

mtry2 <- subset(rf_params, mtry == 500)
i <- 4
points(x = mtry2$ntree, y = mtry2$auc, xlim = range(rf_params$ntree), ylim = range(rf_params$auc), pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "", xlab = "ntree", ylab = "Average AUC")

mtry2 <- subset(rf_params, mtry == 1000)
i <- 5
points(x = mtry2$ntree, y = mtry2$auc, xlim = range(rf_params$ntree), ylim = range(rf_params$auc), pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "", xlab = "ntree", ylab = "Average AUC")

P <- c('100', '200', '300', '500', '1000')
legend("topleft", legend = P, col = COL, bg = "white", horiz=F, pch=PCH, pt.lwd = 1, title = "mtry")

dev.off()


# -------------------------------------------
#
# Doing the same plot for GBM
#
gbm_params <- read.table('GBM_params_table.txt', header = T)
head(gbm_params)

pdf('Rplot_GBM_params_tuning.pdf', width = 12, height = 15) 
par(mfrow=c(3,2))
yrange <- range(0.45,0.85)
CEX <- 1.5

inter <- subset(gbm_params, shrinkage == 0.001 & interaction == 1)
i <- 1
plot(cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, cex.main=CEX, x = inter$ntree, y = inter$mean, xlim = range(gbm_params$ntree), ylim = yrange, pch = PCH[i], col = COL[i], las = 1, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "a) shrinkage = 0.001", xlab = "n.trees", ylab = "Average accuracy")

inter <- subset(gbm_params, shrinkage == 0.001 & interaction == 2)
i <- 2
points(x = inter$ntree, y = inter$mean, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid()) 

inter <- subset(gbm_params, shrinkage == 0.001 & interaction == 10)
i <- 3
points(x = inter$ntree, y = inter$mean, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2, panel.first = grid()) 
P <- c('1', '2', '10')
legend("topleft", legend = P, col = COL, bg = "white", horiz=F, pch=PCH, pt.lwd = 1, title = "int.depth")

# The mean AUC plot
inter <- subset(gbm_params, shrinkage == 0.001 & interaction == 1)
i <- 1
plot(cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, cex.main=CEX, x = inter$ntree, y = inter$auc, xlim = range(gbm_params$ntree), ylim = yrange, pch = PCH[i], col = COL[i], las = 1, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "b) shrinkage = 0.001", xlab = "n.trees", ylab = "Average AUC")

inter <- subset(gbm_params, shrinkage == 0.001 & interaction == 2)
i <- 2
points(x = inter$ntree, y = inter$auc, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 

inter <- subset(gbm_params, shrinkage == 0.001 & interaction == 10)
i <- 3
points(x = inter$ntree, y = inter$auc, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 
P <- c('1', '2', '10')
legend("topleft", legend = P, col = COL, bg = "white", horiz=F, pch=PCH, pt.lwd = 1, title = "int.depth")

# Shrinkage 0.01
# The mean accuracy plot
inter <- subset(gbm_params, shrinkage == 0.01 & interaction == 1)
i <- 1
plot(cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, cex.main=CEX, x = inter$ntree, y = inter$mean, xlim = range(gbm_params$ntree), ylim = yrange, pch = PCH[i], col = COL[i], las = 1, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "c) shrinkage = 0.01", xlab = "n.trees", ylab = "Average accuracy")

inter <- subset(gbm_params, shrinkage == 0.01 & interaction == 2)
i <- 2
points(x = inter$ntree, y = inter$mean, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 

inter <- subset(gbm_params, shrinkage == 0.01 & interaction == 10)
i <- 3
points(x = inter$ntree, y = inter$mean, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 
P <- c('1', '2', '10')
legend("topleft", legend = P, col = COL, bg = "white", horiz=F, pch=PCH, pt.lwd = 1, title = "int.depth")

# The mean AUC plot
inter <- subset(gbm_params, shrinkage == 0.01 & interaction == 1)
i <- 1
plot(cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, cex.main=CEX, x = inter$ntree, y = inter$auc, xlim = range(gbm_params$ntree), ylim = yrange, pch = PCH[i], col = COL[i], las = 1, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "d) shrinkage = 0.01", xlab = "n.trees", ylab = "Average AUC")

inter <- subset(gbm_params, shrinkage == 0.01 & interaction == 2)
i <- 2
points(x = inter$ntree, y = inter$auc, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 

inter <- subset(gbm_params, shrinkage == 0.01 & interaction == 10)
i <- 3
points(x = inter$ntree, y = inter$auc, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 
legend("topleft", legend = P, col = COL, bg = "white", horiz=F, pch=PCH, pt.lwd = 1, title = "int.depth")

# Shrinkage 0.1
# The mean accuracy plot
inter <- subset(gbm_params, shrinkage == 0.1 & interaction == 1)
i <- 1
plot(cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, cex.main=CEX, x = inter$ntree, y = inter$mean, xlim = range(gbm_params$ntree), ylim = yrange, pch = PCH[i], col = COL[i], las = 1, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "e) shrinkage = 0.1", xlab = "n.trees", ylab = "Average accuracy")

inter <- subset(gbm_params, shrinkage == 0.1 & interaction == 2)
i <- 2
points(x = inter$ntree, y = inter$mean, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 

inter <- subset(gbm_params, shrinkage == 0.1 & interaction == 10)
i <- 3
points(x = inter$ntree, y = inter$mean, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 
P <- c('1', '2', '10')
legend("topleft", legend = P, col = COL, bg = "white", horiz=F, pch=PCH, pt.lwd = 1, title = "int.depth")

# The mean AUC plot
inter <- subset(gbm_params, shrinkage == 0.1 & interaction == 1)
i <- 1
plot(cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, cex.main=CEX, x = inter$ntree, y = inter$auc, xlim = range(gbm_params$ntree), ylim = yrange, pch = PCH[i], col = COL[i], las = 1, type = "b", lty = 1, lwd = 2, panel.first = grid(), main = "f) shrinkage = 0.1", xlab = "n.trees", ylab = "Average AUC")

inter <- subset(gbm_params, shrinkage == 0.1 & interaction == 2)
i <- 2
points(x = inter$ntree, y = inter$auc, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 

inter <- subset(gbm_params, shrinkage == 0.1 & interaction == 10)
i <- 3
points(x = inter$ntree, y = inter$auc, pch = PCH[i], col = COL[i], las = 2, type = "b", lty = 1, lwd = 2) 
P <- c('1', '2', '10')
legend("topleft", legend = P, col = COL, bg = "white", horiz=F, pch=PCH, pt.lwd = 1, title = "int.depth")

dev.off()


# -------------------------------------------
#
# Saving SNPs final importance score for both methods
#
scores_ml <- read.table('ML_Final_snps_by_MLs.txt', header = T)
head(scores_ml)
scores_ml$RF_importance <- scores_ml$RF_importance/100
scores_ml$GBM_importance <- scores_ml$GBM_importance/100
write.table(scores_ml, file = "ML_Final_snps_by_MLs_Average.txt", row.names = F, quote = F, sep = "\t")
