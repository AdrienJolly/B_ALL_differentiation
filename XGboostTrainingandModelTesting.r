

#sce is a SingleCellExperiment, a subset of the dataset containing the presumed cluster of origin
#containing exclusively cells with an assigned differentiation class as described in the manuscript


library(scran)
#selection of highly variable genes to be used for training
dec.sce = modelGeneVar(sce)
HVGs = getTopHVGs(dec.sce, n=500)



# separation in training and test set
TestSample = sample(1:ncol(sce),round(0.10*ncol(sce)))
sceTest= sce[,TestSample]
sceTrain = sce[,-TestSample]

# identify cells of the differentiation classes in training set (here 3 classes)
sceTrain_0 = which(sceTrain$class%in%"0")
sceTrain_1 = which(sceTrain$class%in%"1")
sceTrain_2 = which(sceTrain$class%in%"2")



# we use the Shannon entropy (denotes clonal diversity) within each train dataset to decide of sample size per class in the training set

library(entropy)
entropy(as.vector(table(sceTrain[,sceTrain_0]$LinBarcode)))
#[1] 1.355643
entropy(as.vector(table(sceTrain[,sceTrain_1]$LinBarcode)))
#[1] 1.688756
entropy(as.vector(table(sceTrain[,sceTrain_2]$LinBarcode)))
#[1] 1.36805


#For each training iteration, we set the number of cells used for training the smallest class to 100 
#and we adjust the numbers of cells of the other classes based on their expected relative complexity.
#here for instance 1.36805/1.688756 = 0.81 so we set the number of cells for class 2 to 80 cells


library(xgboost)

#we train our model iteratively on the subset of the training set, testing on the remaing part of the training set 
# and we keep the best performing model as final model
nrounds = 50
dis = vector()
dis2 = vector()
samples = matrix(nrow=nrounds,ncol=80+100+80)
for(i in 1:nrounds)
{
  sample0 = sample(sceTrain_0,80)
  sample1 = sample(sceTrain_1,100)
  sample2 = sample(sceTrain_2,80)
  samples[i,] = c(sample0,sample1,sample2)

  sampled = samples[i,]
  dtrain = xgb.DMatrix(data = t(logcounts(sceTrain[HVGs,sampled])), label = sceTrain$class[sampled])
  bstDMatrix = xgboost(data = dtrain, max.depth = 100, eta = 1, nthread = 50, nrounds = 100, objective = "multi:softmax",num_class=3,verbose = 0)
  pred = predict(bstDMatrix, t(logcounts(sceTrain[HVGs,-sampled])))
  
  test=sceTrain[,-sampled]$class
  
  test_0 = which(test%in%0)
  test_1 = which(test%in%1)
  
  test_2 = which(test%in%2)
  
  test_samp_0 = sample(test_0,length(test_1))
  test_samp_2 = sample(test_2,length(test_1))
  
  test_toBetested = c(test_samp_0,test_1,test_samp_2)
  
  err_0 = length(which(abs(pred[test_0]-test[test_0])!=0))/length(test[test_0])
  err_1 = length(which(abs(pred[test_1]-test[test_1])!=0))/length(test[test_1])
  err_2 = length(which(abs(pred[test_2]-test[test_2])!=0))/length(test[test_2])
  
  dis[i] = mean(c(err_1,err_2))
  
}
dis[dis=="NaN"] = 1

minimum = which.min(dis)

sampled = samples[minimum,]

dtrain = xgb.DMatrix(data = t(logcounts(sceTrain[HVGs,sampled])), label = sceTrain$class[sampled])
bstDMatrix = xgboost(data = dtrain, max.depth = 100, eta = 1, nthread = 50, nrounds = 100, objective = "multi:softmax",num_class=3,verbose = 0)

#here we get the balanced accuracy on the real test set (sceTest)

pred = predict(bstDMatrix, t(logcounts(sceTest[HVGs,])))

comparison = rbind(pred,sceTest$class)


error = which(abs(comparison[1,]-comparison[2,])!=0)

recall = 1-table(comparison[2,error])/table(comparison[2,])

names(recall) = c("Class 0","Class 1", "Class 2")

balanced_accuracy = mean(recall)











############################
# here we test significance of prediction with permutations





sceTrainRessamp = sceTrain


ressampledPred = matrix(ncol=3,nrow=1000)


#now we permute the class labels and perform the same procedure as above 
#starting from the entropy based cell number selection

for(j in 1:1000)
{
  sceTrainRessamp$class = sceTrainRessamp$class[sample(1:ncol(sceTrainRessamp),ncol(sceTrainRessamp),replace=FALSE)]
  

  sceTrainRessamp_0 = which(sceTrainRessamp$class%in%"0")
  sceTrainRessamp_1 = which(sceTrainRessamp$class%in%"1")
  sceTrainRessamp_2 = which(sceTrainRessamp$class%in%"2")
  
  entropy0 = entropy(as.vector(table(sceTrainRessamp[,sceTrainRessamp_0]$LinBarcode)))
  
  entropy1 = entropy(as.vector(table(sceTrainRessamp[,sceTrainRessamp_1]$LinBarcode)))
  
  entropy2 = entropy(as.vector(table(sceTrainRessamp[,sceTrainRessamp_2]$LinBarcode)))
  
  n0= round(100*(entropy0/entropy1))
  
  n2 = round(100*(entropy2/entropy1))
  
  
  
  library(xgboost)
  nrounds = 50
  dis = vector()
  dis2 = vector()
  samples = matrix(nrow=nrounds,ncol=n0+100+n2)
  for(i in 1:nrounds)
  {
    sample0 = sample(sceTrainRessamp_0,n0)
    sample1 = sample(sceTrainRessamp_1,100)
    sample2 = sample(sceTrainRessamp_2,n2)
    samples[i,] = c(sample0,sample1,sample2)
    sampled = samples[i,]
    dtrain <- xgb.DMatrix(data = t(logcounts(sceTrainRessamp[HVGs,sampled])), label = sceTrainRessamp$class[sampled])
    bstDMatrix <- xgboost(data = dtrain, max.depth = 100, eta = 1, nthread = 50, nrounds = 100, objective = "multi:softmax",num_class=3,verbose = 0)
    pred <- predict(bstDMatrix, t(logcounts(sceTrainRessamp[HVGs,-sampled])))
    
    test = sceTrainRessamp[,-sampled]$class
    
    test_0 = which(test%in%0)
    test_1 = which(test%in%1)
    
    test_2 = which(test%in%2)
    
    test_samp_0 = sample(test_0,length(test_1))
    test_samp_2 = sample(test_2,length(test_1))
    
    test_toBetested = c(test_samp_0,test_1,test_samp_2)
    
    err_0 = length(which(abs(pred[test_0]-test[test_0])!=0))/length(test[test_0])
    err_1 = length(which(abs(pred[test_1]-test[test_1])!=0))/length(test[test_1])
    err_2 = length(which(abs(pred[test_2]-test[test_2])!=0))/length(test[test_2])
    
    dis[i] = mean(c(err_1,err_2))
  }
  dis[dis=="NaN"]=1
  
  minimum = which.min(dis)
  
  sampled = samples[minimum,]
  
  dtrain = xgb.DMatrix(data = t(logcounts(sceTrainRessamp[HVGs,sampled])), label = sceTrainRessamp$class[sampled])
  bstDMatrix = xgboost(data = dtrain, max.depth = 100, eta = 1, nthread = 30, nrounds = 100, objective = "multi:softmax",num_class=3,verbose = 0)
  
  pred = predict(bstDMatrix, t(logcounts(sceTest[HVGs,])))
  comparison = rbind(pred,sceTest$class)
  error = which(abs(comparison[1,]-comparison[2,])!=0)
  
  ressampledPred[j,] = table(comparison[2,error])/table(comparison[2,])
}

ressampledPred[,1] = NULL

#balanced accuracies of permutations
ressampledPredBalancedAc = rowMeans(1-ressampledPred)

mean(ressampledPredBalancedAc)


avrand = mean(ressampledPredBalancedAc)
sdrand = sd(ressampledPredBalancedAc)

# here we plot the true balanced accuracy versus mean +/- sd of permutations balanced ac.
pdf("BalancedAccuracy.pdf",width=3)
accplot = barplot(c(balanced_accuracy,mean(ressampledPredBalancedAc)),ylim=c(0,1), main= "balanced accuracy",col=c("darkgrey","lightgrey"),names=c("prediction","random"))
arrows(accplot[2],avrand-sdrand,accplot[2],avrand+sdrand,lwd = 1,angle = 90,code = 3)
segments(x0 = accplot[1], y0 = 0.8, x1 = accplot[2], y1 = 0.8, col = "black", lwd = 1.5)
text(x=(accplot[1]+0.5*(accplot[2]-accplot[1])),y=0.85,"**")
dev.off()

