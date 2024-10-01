#====================================================
# Title: Binary Hiking Optimization with Genetic Operators for Gene Selection
# Author: Elnaz Pashaei
# Date: 1/10/2024
# Description: This code implements a binary version of the Hiking Optimization Algorithm
#              with genetic operators (crossover and mutation) for gene selection.
#====================================================

# Clean memory and load libraries

rm(list = ls(all = TRUE)) 
gc()

# Load necessary libraries
library(caret)
library(plsgenomics)
library(farff)
library(klaR)
library(e1071)
library(naivebayes)
library(gamesGA)

# Path to the dataset (ARFF format)
#x <-"~/DARWIN.arff"
x <-"~/Zscore-prostate-DEGS-297-0.8FC.arff"

# Set random seed for reproducibility
set.seed(2)

# Load the dataset
data <- readARFF(x);data <- as.data.frame(data, check.names = TRUE)
class <- data[, ncol(data)]

# Convert class to numeric if using logit.spls
 class <- as.numeric(class)-1; data$class <- class

# Ensure valid factor levels for classification
levels(data$class) <- make.names(levels(data$class))
str(data)

# Define global variables
col = ncol(data)-1
oout <- matrix(0, nrow = 1, ncol = 3)

#====================================================
# Fitness function to evaluate feature subset quality
#====================================================
fitness <- function(star) {
  index <- which(star == 1)
  if (length(index) == 0 | length(index) == 1) {
    return(0)
  }
  subData <- data[, index]
  n2 = ncol(subData)
  
  dimo <- dim(subData)
  subData$class <- data$class
  #----------------------------------
  # Setting k-fold cross-validation
  type <- "k-fold" 
  kf = 5
  
  score = list()
  for (zz in 1:3) {
    if (type == "k-fold") {
      a <- createFolds(subData$class, k = kf, list = FALSE)
      subData$id <- a
      t = kf
      
    } else{
      t <- nrow(subData)
      subData$id <- 1:t
    }
    #============================================================================
    
    for (i in 1:t) {
      train1 <- subData[subData$id != i,-ncol(subData)]# delete id
      test1 <-  subData[subData$id == i,-ncol(subData)]#
      test_lable <- test1$class
      test2 <- test1[, -ncol(test1)]
      train2=train1[,-ncol(train1)]
      
      trainyy <- train1$class
      testyy <- test1[, ncol(test1)]
      
      
      
      model <-logit.spls(Xtrain = train1[, -ncol(train1)],Ytrain = trainyy,Xtest=test2, lambda.ridge = 0.01,lambda.l1 =0.1,ncomp =4,adapt = TRUE,maxIter = 10,svd.decompose = TRUE);
      results =sum(model$hatYtest==test_lable)/length(test_lable);fv=0.8*results+0.2*((col-ncol(train2))/col);accuracy1<-fv;
      #*************
      #model = svm(trainyy ~ .,data = train1[,-ncol(train1)], type = 'C-classification',kernel = "linear", scale = FALSE);pred <-predict(model, test2,type="class")#gamma =0.001, cost = 1, degree=2)
      #model<- naive_bayes(trainyy~., data= train1[,-ncol(train1)]);pred <-predict(model, test2,type="class");
      #--- Fitness function: accuracy and feature selection trade-off (SVM,RF, NB)
      #results =(sum(test_lable==pred)/length(pred)); fv=0.8*results+0.2*((col-ncol(train2))/col); accuracy1<-fv;
      score[[i]] = round(accuracy1, 3)
      
    }
    oout[zz] <- round(mean(unlist(score)), 3)
  }
  return(mean(oout))
  
}
#=====================================================
dtreecv5 <- function(data,xx){
  mydata=data
  class <-mydata[, ncol(mydata)];
  q = 0;q = which(xx>0)
  if((length(q) == 0)||(length(q) == 1)){
    return(0)
  } else {	
    tmp=0;tmp=mydata[,q]
    tmp= cbind(tmp,class)
  }
  
  n2 = ncol(tmp);n <- nrow(tmp);h=n;allaccuracy=matrix(0,nrow=1,ncol=h);
  a <- createFolds(class,k=h, list=FALSE)
  tmp$id <-a;
  #******************************************************************
  for(i in 1:h){
    train <- tmp[tmp$id != i, -ncol(tmp)]; trainyy <- tmp[tmp$id != i,n2];form <- as.formula(trainyy~.)
    test <-  tmp[tmp$id == i, -ncol(tmp)]; testyy <- tmp[tmp$id == i,n2]
    
    #t<- naive_bayes(trainyy~., data= train[,-n2],usekernel = FALSE);pred <- predict(t, test[, -n2],type="class");accuracy1 =0; accuracy1 <- sum(pred == testyy)/length(pred);
    t <- logit.spls(Xtrain=train[,-n2], Ytrain=trainyy, lambda.ridge = 0.01,lambda.l1 =0.1,ncomp =4,adapt = TRUE,maxIter = 10,svd.decompose = TRUE, Xtest=test[, -n2]);accuracy1 =0; accuracy1 <- sum(t$hatYtest == testyy)/length(testyy);
    allaccuracy[1,i] = accuracy1;
  }  
   return(mean(allaccuracy))
}

#====================================================
# Additional functions: Crossover, Mutation, Helpers
#====================================================
# Mutation: Bit Flip
bit_flip_mutation <- function(position, prob) {
  mutated_position <- position
  for (i in seq_along(position)) {
    if (runif(1) < prob) {
      mutated_position[i] <- 1 - position[i]
    }
  }
  return(mutated_position)
}
#========================================================
# Single Point Crossover
SinglePointCrossover <- function(i1,i2)
{
  nVar=length(i1)
  m= sample(1:(nVar-1),size=1)
  O1= c(i1[1:m],i2[(m+1):nVar])
  O2= c(i2[1:m],i1[(m+1):nVar])
  b=rbind(O1,O2)
  return (data.frame(b))
}

#=================================================
# Crossover wrapper
CrossoverX <- function(i1,i2)
{
  i1=array(as.numeric(unlist(i1)))
  i2=array(as.numeric(unlist(i2)))
  oo=SinglePointCrossover(i1,i2);
  return (oo)  
}  

#===================================================
# Count number of features selected
nfeature <- function(vec) {
  #vec <- my.stars[i,]
  vec <- as.numeric(vec)
  a <- table(vec)
  out <- a[names(a) == 1]
  if (length(out) == 0) {
    out = 0
  }
  return(out)
}
#====================================================
# Main Optimization Loop: Binary Hiking Algorithm (BHOA)
#====================================================

#Parameter Setting
dim(data)
row = nrow(data)
N = 50
SearchAgents = N
Tmax = 50
LB = 0
UB = 1
nVar = col
Dim = nVar

best_Pos = matrix(0, nrow = 1, ncol = Dim)
solo = matrix(0, nrow = 1, ncol = Dim)
best_Score = -Inf
Positions = round((matrix(runif((N) * (nVar)), (N))))
newPositions = matrix(0, nrow = N, ncol = Dim)
PDConv = matrix(0, nrow = 1, ncol = Tmax)
numberGene = matrix(0, nrow = 1, ncol = Tmax)
M = 0
t = 0
fit = matrix(0, nrow = 1, ncol = nrow(Positions))
CBest = matrix(0, nrow = 1, ncol = nrow(Positions))

# Initial fitness evaluation
for (i in 1:nrow(Positions)) {
  fit[1, i] = fitness(Positions[i, ])
  fit[is.na(fit)] <- 0
  if (fit[1, i] > best_Score) {
    best_Score = fit[1, i]
    best_Pos = Positions[i, ]
    
  }
}

ptm  <-  proc.time()
sum = 0

# Start the optimization process
while (t < Tmax + 1) {
  if(max(fit) > best_Score || (max(fit) == best_Score && nfeature(Positions[which.max(fit), ])[[1]] < nfeature(best_Pos)[[1]])){
    best_Score = max(fit)
    best_Pos = Positions[which.max(fit), ]
  }
  for (i in 1:SearchAgents) {
    Xini = Positions[i, ]  # obtain initial position of jth hiker
    theta <-sample(0:50, 1)    # randomize elevation angle of hiker
    s <- tan(theta)             # compute slope
    SF <-sample(1:2, 1)        # sweep factor generate either 1 or 2 randomly
    
    Vel= 6 * exp(-3.5 * abs(s + 0.05))  # Compute walking velocity based on Tobler's Hiking Function
    newVel=Vel + runif(nVar) * (best_Pos - SF * Xini)  # determine new position of hiker
    
    newPositions <-Positions[i,] + newVel  # Update position of hiker
    
    newPositions <-pmin(UB, newPositions)   # bound violating to upper bound
    newPositions <-pmax(LB, newPositions)   # bound violating to lower bound
    for (j in (1:Dim)) {
      if (abs(tanh(newPositions[j])) > runif(1)) {
        newPositions[j] = 1
      } else{
        newPositions[j] = 0
      }
    }
    #****************Activate Genetic Operators: if FALSE original BHOA
    if(TRUE){
      if(runif(1)<0.5) {
        P = CrossoverX(best_Pos, newPositions)
        fp1 = fitness(P[1,])
        fp2 = fitness(P[2,])
        
        if ((fp1 > fp2 && fp1 > best_Score) || (fp1 == best_Score && nfeature(P[1, ])[[1]] < nfeature(best_Pos)[[1]])){
          best_Pos = as.numeric(P[1, ])
          best_Score = fp1
          cat("crossover1", '\n');cat(fp1, '\n')
          newPositions = as.numeric(P[1,]);
        } else if ((fp2 > best_Score) || (fp2 == best_Score && nfeature(P[2, ])[[1]] < nfeature(best_Pos)[[1]])){
          best_Pos =as.numeric(P[2, ])
          best_Score = fp2
          cat("crossover2", '\n');cat(fp2, '\n')
          newPositions = as.numeric(P[2,]);
        }
      }
      if(runif(1)>0.5) {
        GOm1 = bit_flip_mutation(best_Pos,0.9+(-0.9*(1 - t)/Tmax-1));fsolom1 = fitness(GOm1);
        GOm2 = gamesGA::mutation(best_Pos, prob =0.01);fsolom2 = fitness(GOm2);
        
        if (fsolom1 > fsolom2) {GOm <- GOm1;fsolom <- fsolom1} else {GOm <- GOm2;fsolom <- fsolom2}
        
        if(fsolom > best_Score || (fsolom == best_Score && nfeature(GOm)[[1]] < nfeature(best_Pos)[[1]])){
          cat("mutation:", fsolom, '\n')
          best_Score = fsolom
          best_Pos = GOm
          newPositions = GOm;
        }
      }
    }
    fnew <- fitness(newPositions)
    if (fnew > fit[1,i]||((fit[1,i] == fnew) &(nfeature(newPositions)[[1]] < nfeature(Positions[i,])[[1]]))) {
      Positions[i, ] <- newPositions         # store best position
      fit[i] <- fnew             # store new fitness 
    }
    
  }
  
  t = t + 1
  cat("At iteration",(t),"the best solution fitness is",(best_Score),"Number of feature:",(sum(best_Pos)),"\n")
  PDConv[t] = best_Score
  numberGene[t] = (sum(best_Pos))
}
e = proc.time() - ptm
cputime = 0
cputime = e[1] + e[2]
# Final outputs
print(best_Pos)
q1 = which(best_Pos > 0)
print(q1)
print(cputime)
print(length(q1))
xx = 0
x111 = best_Pos
xx = best_Pos
PDConv
numberGene

#LOOCV accuracy
acc=dtreecv5(data,xx); cat("LOOCV accuracy:",acc,"\n")

