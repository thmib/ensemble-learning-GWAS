#
# A new script for applying RF to the crc dataset
# This script applies RF to the data using 10-fold CV
# for different combinations of paramters
# In the end, the combination with the best accuracy is selected
# for further analysis
#
# RF Tuning parameters:
#     ntrees (# trees)
#     mtry (number of features to randomly selected to make best split at each not of a tree)
#
# Required packages: ranger, pROC
#

library(ranger)
library(caret)
library(pROC)

setwd("~/Documents/crc_new/")
getwd() 

data <- read.table('~/Dropbox/crc_new/dataset-2k-turf-caret.txt', header = T)
dim(data)


# -------------------------------------------
#
# converting SNPs to factors
#
for(i in 1:(ncol(data)-1)) { 
  data[,i] <- factor(data[,i])
}


# -------------------------------------------
#
# defining different ranges of values for parameters
# we used expert knowledge for choosing these values
# run this algorithm with different runs, repeat each configuration 10 times
#
RF_Frame <- data.frame(Model = "RF", MTRY = 0, NTREE = 0, Iteration = 0, Mean = 0, Median = 0, SD = 0, Specificity = 0, Sensitivity = 0, AUC = 0)
RF_Frame <- RF_Frame[-1,]
trees <- c(500,1000,2000)
mtrys <- c(100,200,300,500,1000,2000)
for (N in trees) {
  for (M in mtrys) {
    for (k in 1:10) {
      D <- Normal_RF(data, MTRY = M, NTREE = N, Iter = k)
      RF_Frame <- rbind(RF_Frame, D)
      print(RF_Frame)
    }
  }
}
write.table(RF_Frame, file = "RF_runs.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
#
# calculate the average of accuracies over 10 runs of each configuration
# specify which combination of parameter values has the highest accuracy
#
RF_summary <- data.frame(Model = "RF", NTREE = 0, MTRY = 0, Mean = 0, Median = 0, SD = 0, Specificity = 0, Sensitivity = 0, AUC = 0)
RF_summary <- RF_summary[-1,]
for (N in trees) {
  for (M in mtrys) {
    C1 <- subset(RF_Frame, NTREE == N & MTRY == M)
    D <- data.frame(Model = "RF", NTREE = N, MTRY = M, Mean = mean(C1$Mean), Median = mean(C1$Median), 
                    SD = mean(C1$SD), Specificity = mean(C1$Specificity), Sensitivity = mean(C1$Sensitivity), 
                    AUC = mean(C1$AUC))
    RF_summary <- rbind(RF_summary, D)
  }
}
write.table(RF_summary, file = "RF_summary.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
#
# based on the best configuration, detect the most significant SNPs from all runs
# for each run, determine the top 100 significnat SNPs and record their appearance
#
(best <- RF_summary[which.max(RF_summary$Mean),])
rnames <- colnames(data)
stats <- data.frame(name=rnames, importance=numeric(length(rnames)),count=numeric(length(rnames)))
for (j in 1:10) {
  SS <- Top_Features_RF(data, nfold = 10, MTRY = best$MTRY, NTREE = best$NTREE)
  row <- 1
  for(row in 1:nrow(SS)){
    featureName <- as.character(SS[row,1]);
    importance <- SS[row,2];
    stats[stats$name== featureName,"importance"] <- stats[stats$name==featureName,"importance"] + importance;
    stats[stats$name== featureName,"count"] <- stats[stats$name==featureName,"count"] + SS[SS$name==featureName,"count"];
  }
}
stats <- stats[order(stats[,3], decreasing = TRUE),]
write.table(stats, file = "RF_snps_importance_frequency_N2000_M100.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
#
# the plot of freuqncy of significance of SNPs
#
DEN <- density(stats[,3])
plot(DEN, las = 1, xlab = "Count", ylab = "Density", main = "Density of SNPs significance")
#hist(stats[,3], las = 1, xlab = "Count", ylab = "Frequency", main = "Histogram of SNPs significance")
#plot(stats[,3], las = 1, xlab = "SNP index", ylab = "Frequency", main = "Frequency of significance of SNPs")
abline(v = 75, col = "red", lty = 2)
S <- subset(stats, count > 75)
S <- S[order(S[,3], decreasing = T),]
S[,2] <- S[,2]/100
head(S)
write.table(S, file = "RF_snps_top.txt", quote = F, row.names = F, sep = "\t")


# -------------------------------------------
#
# functions for applying RF
#
Top_Features_RF <- function(data, MTRY = 500, NTREE = 500, nfold = 10, Iter = 10, stats = NULL) {
  if (is.null(stats)) {
    rnames <- colnames(data)
    stats <- data.frame(name=rnames, importance=numeric(length(rnames)),count=numeric(length(rnames)));
  } else {
    print("Stats is not null")
  }
  MTRY = 100; NTREE = 2000; nfold = 10;
  randcv <- sample(rep(1:nfold, length = nrow(data)), nrow(data), replace = FALSE)
  for (i in sort(unique(randcv))) {
    print(paste("-----------", i, "----------------"))
    train <- data[randcv!=i,]
    test <- data[randcv==i,]

    model <- ranger(Class ~ ., data = train, probability = T, mtry = MTRY, num.trees = NTREE, importance = "impurity")
    
    Imp <- importance(model)
    Imp <- Imp[order(Imp, decreasing = TRUE)]
    pdf('Rplot_RF_var_importance.pdf')
    plot(Imp, las = 1, xlab = "SNP index", ylab = "Importance", main = "", panel.first = grid()) #Variables importance
    dev.off()

    Imp <- importance(model)
    X <- data.frame(name = names(Imp), score = Imp)
    for(row in 1:nrow(X)){
      featureName <- as.character(X[row,1]);
      importance <- X[row,2];
      stats[stats$name== featureName,"importance"] <- stats[stats$name==featureName,"importance"] + importance;
      stats[stats$name== featureName,"count"] <- stats[stats$name==featureName,"count"] + 1;
    }
  }
  return(stats)
}


# -------------------------------------------
# applying ranger method to get the accuracy and AUC
# for a run of RF
#
Normal_RF <- function(data, MTRY = 500, NTREE = 500, nfold = 10, Iter = 10) {
  randcv <- sample(rep(1:nfold, length = nrow(data)), nrow(data), replace = FALSE)
  accu_list <- c()
  specifity <- c()
  sensivity <- c()
  aucs <- c()
  for (i in sort(unique(randcv))) {
    print(paste("-----------", i, "----------------"))
    train <- data[randcv!=i,]
    test <- data[randcv==i,]
    print(dim(train))
    print(dim(test))
    
    model <- ranger(Class ~ ., data = train, probability = T, mtry = MTRY, num.trees = NTREE, importance = "impurity")
    pred <- predict(model, data = test[, -ncol(test)])
    head(pred$predictions)
    
    ROC <- roc(predictor = pred$predictions[,2], response = test$Class)
    aucs[i] <- as.numeric(ROC$auc)
    
    #C <- confusionMatrix(pred$predictions, test$Class)
    lab <- ifelse(pred$predictions[,1] > pred$predictions[,2], "C1", "C2")
    C <- confusionMatrix(lab, test$Class)
    C
    print(C$table)  	
    print(paste("Fold", i, ": Accuracy =", C$overall["Accuracy"])) 
    
    accu_list[i] <- C$overall["Accuracy"]
    specifity[i] <- C$byClass["Specificity"]
    sensivity[i] <- C$byClass["Sensitivity"]

  }
  print(accu_list)
  D <- data.frame(Model = "RF", MTRY = MTRY, NTREE = NTREE, Iteration = Iter, 
                  Mean = mean(accu_list), Median = median(accu_list), SD = sd(accu_list),
                  Specificity = mean(specifity), Sensitivity = mean(sensivity), AUC = mean(aucs))
  print(D)
  return(D)
}

