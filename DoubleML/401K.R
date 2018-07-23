###########################################################################################################################
#  This is an example of ATE estimation of 401(k) eligibility on accumulated assets using Double Machine Learning Methods 
#  References: "Double/Debiased Machine Learning of Treatment and Causal Parameters",  AER P&P 2017     
#              "Double Machine Learning for Treatment and Causal Parameters",  Arxiv 2016               

# Data source: SIPP 1991 (Abadie, 2003)
# Description of the data: the sample selection and variable contruction follow

# Abadie, Alberto (2003), "Semiparametric instrumental variable estimation of treatment response 
# models," Journal of Econometrics, Elsevier, vol. 113(2), pages 231-263, April.

# The variables in the data set include:

# net_tfa:  net total financial assets
# e401:     = 1 if employer offers 401(k)
# age
# inc:      income
# fsize:    family size
# educ:     years of education
# db:       = 1 if indivuduals has defined benefit pension
# marr:     = 1 if married
# twoearn:  = 1 if two-earner household
# pira:     = 1 if individual participates in IRA plan
# hown      = 1 if home owner
###########################################################################################################################

###################### Loading packages ###########################

library(foreign);
library(quantreg);
library(mnormt);
library(gbm);
library(glmnet);
library(MASS);
library(rpart);
library(doParallel)
library(sandwich);
library(hdm);
library(randomForest);
library(nnet)
library(matrixStats)
library(quadprog)
library(ivmodel)
library(xtable)

###################### Loading functions and Data ##############################

rm(list = ls())  # Clear everything out so we're starting clean
source("ML_Functions.R")  
source("Moment_Functions.R")  
options(warn=-1)
set.seed(1211);
cl   <- makeCluster(2, outfile="")

vec.pac= c("foreign", "quantreg", "gbm", "glmnet",
           "MASS", "rpart", "doParallel", "sandwich", "randomForest",
           "nnet", "matrixStats", "xtable", "readstata13", "car", "lfe", "doParallel",
           "caret", "foreach", "multcomp","cowplot", "stringr", "abind")

data  <- read.dta("sipp1991.dta");

################################ Inputs ########################################

# Outcome Variable
y      <- "net_tfa";

# Treatment Indicator
d      <- "e401";  

# Controls
x      <- "age + inc + educ + fsize + marr + twoearn + db + pira + hown" # use this for tree-based methods like forests and boosted trees
xl     <- "(poly(age, 6, raw=TRUE) + poly(inc, 8, raw=TRUE) + poly(educ, 4, raw=TRUE) + poly(fsize, 2, raw=TRUE) + marr + twoearn + db + pira + hown)^2";  # use this for rlasso etc.

# Model names. For a list of available model names in caret package see: http://topepo.github.io/caret/available-models.html
# some available models: svmLinear2, svmLinear, svmPoly, gbm, rf, gamboost, avNNet, pcaNNet, nnet
methods    <- c("glmnet", "rf", "nnet")   

# A list of arguments for models used in the estimation
args       <- list(svmLinear2=list(type='eps-regression'), svmLinear=list(type='nu-svr'), svmPoly=list(type='nu-svr'), gbm=list(verbose=FALSE), rf=list(ntree=1000), gamboost=list(baselearner='btree'), avNNet=list(verbose = 0, linout = TRUE, trace = FALSE), pcaNNet=list(linout = TRUE, trace = FALSE, MaxNWts=100000, maxit=10000), nnet=list(linout = TRUE, trace = FALSE, MaxNWts=100000, maxit=10000))

methodML   <- c("repeatedcv", "none","repeatedcv")   # resampling method for chosing tuning parameters. available options: boot, boot632, cv, LOOCV, LGOCV, repeatedcv, oob, none
tune       <- c(2, NA,2)                             # number of elements per parameter in the grid. the grid size is tune^{number of tuning parameters}.
proces     <- c("range", "range", "range")           # pre-processing method
select     <- c("best",  NA, "best")                 # optimality criteria for choosing tuning parameter in cross validation. available options: best, oneSE, tolerance
cv         <- c(2,NA,2)                               # the number of folds in cross-validation
rep        <- c(2, NA,2)                             # number of iteration in repeated cross-validations

# put all options together
param      <- list(methodML=methodML, tune=tune, proces=proces, select=select, cv=cv, rep=rep)

# If there is a parameter of the model that user doesn't want to choose with cross validation, it should be set using tune_param variable. Below mtry of random forest is set to 5 
# for glmnet we want to choose both tuning parameters using cross validation so it is set to NULL
tune_param <- list(0) 
tune_param[[1]]  <- 0
tune_param[[2]]  <- data.frame(mtry=5)
tune_param[[3]]  <- 0

############## Arguments for DoubleML function:

# number of times DML estimator is repeated
split <- 2

# data:     : data matrix
# y         : outcome variable
# d         : treatment variable
# z         : instrument  (NULL if not IV estimation)
# xx        : set of raw control variables
# xL        : controls for penalized linear methods
# methods   : machine learning methods
# DML       : DML1 or DML2 estimation (DML1, DML2)
# nfold     : number of folds in cross fitting
# est       : estimation methods (IV, LATE, plinear, interactive)
# arguments : list of arguments for machine learning models
# ensemble  : ML methods used ine ensemble method
# silent    : if FALSE print status messages
# trim      : bounds for propensity score trimming


r <- foreach(k = 1:split, .combine='rbind', .inorder=FALSE, .packages=vec.pac) %dopar% { 
  
  dml <- DoubleML(data=data, y=y, d=d, z=NULL, xx=x, xL=xl, methods=methods, DML="DML1", nfold=2, est="plinear", args=args,  silent=FALSE, trim=c(0.01,0.99), param=param, tune_param=tune_param) 
  
  data.frame(t(dml[1,]), t(dml[2,]))
  
}

################################  Process Output to Create Output Table ########################################

result           <- matrix(0,3, length(methods)+1)
colnames(result) <- cbind(t(methods), "best")
rownames(result) <- cbind("Median ATE", "se(median)",  "se")

result[1,]        <- colQuantiles(r[,1:(length(methods)+1)], probs=0.5)
result[2,]        <- colQuantiles(sqrt(r[,(length(methods)+2):ncol(r)]^2+(t(t(r[,1:(length(methods)+1)]) - colQuantiles(r[,1:(length(methods)+1)], probs=0.5))^2)), probs=0.5)
result[3,]        <- colQuantiles(r[,(length(methods)+2):ncol(r)], probs=0.5)

result_table <- round(result, digits = 0)

for(i in 1:ncol(result_table)){
  for(j in seq(2,nrow(result_table),3)){
    
    result_table[j,i] <- paste("(", result_table[j,i], ")", sep="")
    
  }
  for(j in seq(3,nrow(result_table),3)){
    
    result_table[j,i] <- paste("(", result_table[j,i], ")", sep="")
    
  }
}

print(xtable(result_table, digits=3))

