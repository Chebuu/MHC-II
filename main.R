### ###
### Without chosing random non-binding peptides based on a probability distribution of amino acid frequencencies from the known binding peptides, 
### the dnn achieves accuracies upwards of 90%. However, when creating random peptides based on a probability distribution, the accuracy falls to around 75%.

library(Rcpi)
library(protr)
library(RSNNS)
library(deepnet)
library(caret)

source('lib.R')

# data <- read.csv('epitope_table_export copy.csv', header=T, skip=1)
load('MCHIIPeptides.RData')

peptides.positive <- as.character(data[1:100,'Description']) # The amino acid sequences of known MHCII binding peptides
peptides.positive <- peptides.positive[sapply(peptides.positive, protcheck)] # Check validity of peptides as strings


peptides.negative <- randomPeptides(length(peptides.positive), ## Adding randomly generated peptides to the dataset as controls based on the distribution of amino acid frequencies found in peptides.positive. These peptides are assumed to be non-binders (negative samples).
                                    min.aas=8, max.aas=20, # Randomly generated peptides ranging from 8 to 20 amino acids in length
                                    prob=computePriorDist(peptides.positive)) 
peptides.negative <- peptides.negative[sapply(peptides.negative, protcheck)] # Check validity of peptides
peptides.duplicate <- peptides.negative %in% peptides.positive # Check if any randomly generated negative peptides exist in the positive set
if(any(peptides.duplicate)) peptides.negative <- peptides.negative[-(peptides.duplicate)] # Peptides duplicated across the two sets are removed from the negative set.

peptides <- c(peptides.positive, peptides.negative)

## Generate labels for the peptides 
Y.positive <- rep(1, length(peptides.positive))
Y.negative <- rep(0, length(peptides.negative))
Y_ <- c(Y.positive, Y.negative)

## Create descriptors for the peptides and assemble the dataset
if(length(Y_) != length(peptides)){
  stop('Y_ and peptides are not the same length')
}

descriptors <- c("extractProtCTriad", 
                 "extractProtDC", 
                 "extractProtTC", 
                 "extractProtCTDC", 
                 "extractProtCTDT", 
                 "extractProtCTDD")
X <- makeDescriptors(peptides, descriptors) # A dataframe of descriptors, one row for each peptide
dataset <- as.data.frame(cbind(Y_,X)) # A dataframe of labelled descriptors, unshuffled
dataset <- dataset[sample(1:nrow(dataset)),] # Shuffle rows of the dataset
dataset.split <- splitForTrainingAndTest(x=dataset[,-(1)], y=dataset[,1], 0.25)
# dataset.split <- normTrainingAndTestSet(dataset.split)

## Train a model
model.dnn <- dbn.dnn.train(dataset.split$inputsTrain, dataset.split$targetsTrain, 
                           hidden= c(500, 100), activationfun='tanh', learningrate=0.01, momentum=0.5, learningrate_scale=1.1,
                           batchsize=10, hidden_dropout=0, visible_dropout=0, cd=0.9)

## Get the error
model.dnn.errorRate <- nn.test(model.dnn, dataset.split$inputsTest, dataset.split$targetsTest)
sprintf('Model accuracy = %f', 1 - model.dnn.errorRate)

## Predictions for the test set 
predictions.dnn <- round(nn.predict(model.dnn, dataset.split$inputsTest))
confmat <- confusionMatrix(dataset.split$targetsTest, predictions.dnn)

## Diagnostic accuracy
TP <- confmat[1,1]
TN <- confmat[2,2]
FP <- confmat[2,1]
FN <- confmat[1,2]
total <- sum(rowsums(confmat))

ACC <- (TP + TN) / total; ACC
TPR <- TP / (TP + FN); TPR
SPC <- TN / (TN + FP); SPC
PPV <- TP / (TP + FP); PPV
NPV <- TN / (TN + FN); NPV
FPR <- 1 - SPC; FPR
FNR <- 1 - TPR; FNR
FDR <- 1 - PPV; FDR
F1 <- 2 * TP / (2 * TP + FP + FN); F1
MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) + (TN + FP) + (TN + FN)); MCC
informedness <- TPR + SPC - 1; informedness
markedness <- PPV + NPV - 1; markedness

## Train a model with caret
set.seed(1)
seeds <- vector(mode = "list", length = 26)
for(i in 1:length(seeds)) seeds[[i]] <- sample.int(1000, 12)
seeds[[length(seeds)]] <- sample.int(1000, 12)

ctrl <- trainControl(method = "adaptive_boot",
                     repeats = 5,
                     seeds = seeds)

mod <- caret::train(dataset.split$inputsTrain, 
                    dataset.split$targetsTrain,
                    method = "knn",
                    tuneLength = 12,
                    trControl = ctrl)

