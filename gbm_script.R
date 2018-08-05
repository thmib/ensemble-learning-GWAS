#
# Script for applying gbm on the CRC data
# 
# This script applies GBM to the data using 10-fold CV
# for different combinations of paramters
# In the end, the combination with the best accuracy is selected
# for further analysis
#
# GBM Tuning parameters:
#     n.trees (# Boosting Iterations)
#     interaction.depth (Max Tree Depth)
#     shrinkage (Shrinkage)
#     n.minobsinnode (Min. Terminal Node Size)
# Required packages: gbm, plyr, caret, pROC
#
library(caret)
library(pROC)
library(gbm)

setwd('~/Dropbox/crc_new')


data_gbm <- read.table('dataset-2k-turf-gbm.txt', header = T)

ntree_values <- c(100,200,300,500,1000,2000)
depth_values <- c(1,2,5,10)
shrink_values <- c(0.001,0.01,0.1)
minobs_values <- c(10) #c(1,5,10)

# -------------------------------------------
#
# Run the GBM for 10-fold CV for 10 iterations
#
GBM_Frame <- data.frame(Model = "GBM", NTREE = 0, INTDEPTH = 0, SHRINK = 0, Iteration = 0, 
                Mean = 0, Median = 0, SD = 0, Specificity = 0, Sensitivity = 0, AUC = 0)
GBM_Frame <- GBM_Frame[-1,]
N <- 100; I <- 1; S <- 0.001; M <- 10;
for (N in ntree_values) {
  for (I in depth_values) {
    for (S in shrink_values) {
      for (M in minobs_values) {
        print(paste("Running for NTREE=", N, ", INT.DEPTH=", I, ", SHRINK=", S, ", MINOBS=", M))
        for (k in 1:10) {
          D <- Normal_GBM(data_gbm, nfold = 10, NTREE = N, INT.DEPTH = I, SHRK = S, MIN.OBS = M, Iter = k)
          GBM_Frame <- rbind(GBM_Frame, D)
          print(tail(GBM_Frame))
        }
      }
    }
  }
}
dim(GBM_Frame)
write.table(GBM_Frame, file = "GBM_runs_all.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
#
# Average over every 10 iteration of different configurations
#
Sindex <- seq(1,nrow(GBM_Frame), 10)
Eindex <- seq(10,nrow(GBM_Frame), 10)
AVG <- sapply(1:length(Sindex), function (x) sapply(GBM_Frame[Sindex[x]:Eindex[x],-1], mean) )
GBM_summary <- data.frame(NTREE = 0, INTDEPTH = 0, SHRINK = 0, Iteration = 0, Mean = 0, Median = 0, SD = 0, Specificity = 0, Sensitivity = 0, AUC = 0)
for (i in 1:ncol(AVG)) {
  GBM_summary <- rbind(GBM_summary, AVG[,i])
}
GBM_summary <- GBM_summary[-1,]
write.table(round(GBM_summary,3), file = "GBM_runs_summary_all.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
#
# based on the best configuration, detect the most significant SNPs from all runs
#
(best <- GBM_summary[which.max(GBM_summary$Mean),])

# -------------------------------------------
#
# use the best model, with 10 times repitition to detect the most significant SNPs
#
rnames <- names(data_gbm)[-ncol(data_gbm)]
stats <- data.frame(name=rnames, importance=numeric(length(rnames)),count=numeric(length(rnames)))
for (j in 1:10) {
  SS <- Top_Features_GBM(data_gbm, nfold = 10, NTREE = best$NTREE, INT.DEPTH = best$INTDEPTH, SHRK = best$SHRINK, MIN.OBS = 10, Iter = 1, stats = NULL)
  for(row in 1:nrow(SS)){
    featureName <- as.character(SS[row,1]);
    importance <- SS[row,2];
    stats[stats$name== featureName,"importance"] <- stats[stats$name==featureName,"importance"] + importance;
    stats[stats$name== featureName,"count"] <- stats[stats$name==featureName,"count"] + SS[SS$name==featureName,"count"];
  }
}
stats <- stats[order(stats[,2], decreasing = TRUE),]
write.table(stats, file = "GBM_snps_importance_frequency_N2000_I10_S1.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
#
# Drawing the histogram of SNPs cumulative score
#
pdf('Rplot_GBM_cumulative_score_hist.pdf')
hist(stats[,2], xlab = "Score", las = 1, main = "Histogram of SNPs cumulative score", panel.first=grid())
dev.off()


# -------------------------------------------
#
# function for getting most important features with GBM
#
Top_Features_GBM <- function(data_gbm, nfold = 10, NTREE = 100, INT.DEPTH = 2, SHRK = 0.01, MIN.OBS = 10, Iter = 1, stats = NULL) {
  if (is.null(stats)) {
    rnames <- colnames(data)
    stats <- data.frame(name=rnames, importance=numeric(length(rnames)),count=numeric(length(rnames)));
  } else {
    print("Stats is not null")
  }
  
  randcv <- sample(rep(1:nfold, length = nrow(data_gbm)), nrow(data_gbm), replace = FALSE)

  for (i in sort(unique(randcv))) {
    print(paste("-----------", i, "----------------"))
    train <- data_gbm[randcv!=i,]
    test <- data_gbm[randcv==i,]
    print(dim(train))
    print(dim(test))
    
    gbm_model <- gbm.fit(x = subset(train, select = -Class),
                         y = train$Class,
                         distribution = "bernoulli",
                         n.trees = NTREE, 
                         interaction.depth = INT.DEPTH, 
                         shrinkage = SHRK,
                         n.minobsinnode = MIN.OBS,
                         # weights = wts,
                         #cv.folds = 10, # or > 1
                         bag.fraction = 0.5,
                         keep.data = FALSE, #???
                         response.name = "Class",
                         verbose = FALSE
    )
    
    pdf('Rplot_GBM_relative_influence.pdf')
    summary(gbm_model, ylab = "SNP", main = "", panel.first=grid()) #Relative influence of SNPs by GBM
    dev.off()
    
    Imp <- relative.influence(gbm_model, n.trees = NTREE)
    Imp_Names <- names(Imp)
    for(row in 1:length(Imp)){
      featureName <- as.character(Imp_Names[row]);
      importance <- Imp[row];
      stats[stats$name== featureName,"importance"] <- stats[stats$name==featureName,"importance"] + importance;
      stats[stats$name== featureName,"count"] <- stats[stats$name==featureName,"count"] + 1;
    }
  }
  return(stats)
}


# -------------------------------------------
#
# applying simple gbm to get the accuracy and AUC
#
Normal_GBM <- function(data_gbm, nfold = 10, NTREE = 100, INT.DEPTH = 2, SHRK = 0.01, MIN.OBS = 10, Iter = 1) {
  randcv <- sample(rep(1:nfold, length = nrow(data_gbm)), nrow(data_gbm), replace = FALSE)
  accu_list <- c()
  specifity <- c()
  sensivity <- c()
  aucs <- c()
  i <- 1
  for (i in sort(unique(randcv))) {
    print(paste("-----------", i, "----------------"))
    train <- data_gbm[randcv!=i,]
    test <- data_gbm[randcv==i,]
    print(dim(train))
    print(dim(test))

    set.seed(845)
    gbm_model <- gbm.fit(x = subset(train, select = -Class),
                         y = train$Class,
                         distribution = "bernoulli",
                         n.trees = NTREE, 
                         interaction.depth = INT.DEPTH, 
                         shrinkage = SHRK,
                         n.minobsinnode = MIN.OBS,
                         # weights = wts,
                         #cv.folds = 10, # or > 1
                         bag.fraction = 0.5,
                         keep.data = FALSE, #???
                         response.name = "Class",
                         verbose = FALSE
    )
    # best.iter <- gbm.perf(gbm_model,method="OOB")

    pred <- predict.gbm(gbm_model, newdata = test[, -ncol(test)], n.trees = NTREE, type = "response")
    
    ROC <- roc(predictor = pred, response = test$Class)
    aucs[i] <- as.numeric(ROC$auc)
    
    #C <- confusionMatrix(pred$predictions, test$Class)
    lab <- ifelse(pred > 0.5, 1, 0)
    C <- confusionMatrix(lab, test$Class)
    C
    print(C$table)  	
    print(paste("Fold", i, ": Accuracy =", C$overall["Accuracy"])) 
    
    accu_list[i] <- C$overall["Accuracy"]
    specifity[i] <- C$byClass["Specificity"]
    sensivity[i] <- C$byClass["Sensitivity"]
        
  }
  print(accu_list)
  D <- data.frame(Model = "GBM", NTREE = NTREE, INTDEPTH = INT.DEPTH, SHRINK = SHRK, Iteration = Iter, 
                  Mean = mean(accu_list), Median = median(accu_list), SD = sd(accu_list),
                  Specificity = mean(specifity), Sensitivity = mean(sensivity), AUC = mean(aucs))
  print(D)
  return(D)
}


