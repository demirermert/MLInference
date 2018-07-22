#--------------------------------------------------------------------------------------------------------------------
# This program obtains empirical results for the paper "Generic Machine Learning Discovery 
# and Classification Analysis of Heterogenous Treatment Effects in Randomized Experiments"
# by V. CHERNOZHUKOV, M. DEMIRER, E. DUFLO, I. FERNANDEZ-VAL
#--------------------------------------------------------------------------------------------------------------------

# Authors: V. CHERNOZHUKOV, M. DEMIRER, E. DUFLO, I. FERNANDEZ-VAL

# This program returns three outputs
# 1) A plot for each outcome variable reporting ATE by groups
# 2) het_results: Latex table reporting ATE and heterogeneity loading coefficient 
# 3) test_results: A latex table reporting  estimated ATE_group1 - ATE_group5 and its p-value
# 4) best: A latex table reporting best ML methods
# 5) MSE: A latex table of mean squared errors from tuning

rm(list=ls(all=TRUE))
vec.pac= c("foreign", "quantreg", "gbm", "glmnet",
           "MASS", "rpart", "doParallel", "sandwich", "randomForest",
           "nnet", "matrixStats", "xtable", "readstata13", "car", "lfe", "doParallel",
           "caret", "foreach", "multcomp","cowplot")

lapply(vec.pac, require, character.only = TRUE) 
source("ML_Functions.R")
ptm <- proc.time()

set.seed(1211);
cl   <- makeCluster(100, outfile="")
registerDoParallel(cl)


####################################### Load and Process Data  #######################################

data        <- read.dta13("data_rep.dta")
data$paire  <- factor(data$paire)
a           <- as.data.frame(model.matrix(~data$paire-1))
colnames(a) <- (substring( names(a), 6, 12))
data <- cbind(data, a)
colnames(data)[which(colnames(data)=="paire")] <- "vil_pair"

####################################### Inputs  #######################################

sim     <- 100     # number of repetitions
K       <- 2       # number of folds
p       <- 5       # number of groups 
thres   <- 0.2     # quantile for most/least affected group
alpha   <- 0.05    # signifigance level

#  dimension of these three vectors should match. If dimension is greater than 1 the program runs heterogeneity estimation separately for each outcome variable
names <- c("Amount of Loans",  "Output", "Profit", "Consumption")    # vector of labels for outcome variables
Y     <- c("loansamt_total", "output_total", "profit_total", "consumption")     # vector of outcome variables
D     <- rep("treatment", length(Y))  # vector of treatment variables

# specify cluster, fixed effect and partition
cluster      <- "demi_paire"       # if no cluster       use    cluster      <- "0"
fixed_effect <- "vil_pair"         # if no fixed_effect  use    fixed_effect <- "0"
partition    <- "demi_paire"       # if no partition     use    partition    <- "0"

# create a vector of control variables
controls     <- c("members_resid_bl", "nadults_resid_bl", "head_age_bl", "act_livestock_bl", "act_business_bl", "borrowed_total_bl", "members_resid_d_bl", "nadults_resid_d_bl", "head_age_d_bl", "act_livestock_d_bl", "act_business_d_bl", "borrowed_total_d_bl", "ccm_resp_activ", "other_resp_activ", "ccm_resp_activ_d", "other_resp_activ_d")
controls     <- c(controls, names(data)[(substring( names(data), 1, 5)=="paire")])

affected       <- c("head_age_bl", "members_resid_bl", "act_business_bl", "borrowed_total_bl")      # characteristics for most/least affected analysis
names_affected <- c("head_age_bl", "members_resid_bl", "act_business_bl", "borrowed_total_bl")      # characteristics for most/least affected analysis

# generate formula for x, xl is for linear models
X <- ""

for(i in 1:length(controls)){
  X <- paste(X, controls[i], "+", sep = "")
}
X  <- substr(X, 1, nchar(X)-1)
XL <- paste("(", X , ")", sep="")
XL <- X

######################################################################################################

if(fixed_effect=="0" & cluster=="0"){
  data <- data[,c(Y, D,controls)]
}

if(fixed_effect=="0" & cluster!="0"){
  data <- data[,c(Y, D,controls, cluster)]
}

if(fixed_effect!="0" & cluster=="0"){
  data <- data[,c(Y, D,controls, fixed_effect)]
}

if(fixed_effect!="0" & cluster!="0"){
  data <- data[,c(Y, D, controls, cluster, fixed_effect)]
}

####################################### ML Inputs  #######################################

# svmPoly    : Support Vector Machines with Polynomial Kernel , package: kernlab, tuning parameters: degree (Polynomial Degree), scale (Scale), C (Cost)
# svmLinear  : Support Vector Machines with Linear Kernel , package: kernlab , C (Cost)
# svmLinear2 : Support Vector Machines with Linear Kernel , package: e1071 , cost (Cost)
# gbm        : Stochastic Gradient Boosting , package: gbm  , n.trees (# Boosting Iterations), interaction.depth (Max Tree Depth), shrinkage (Shrinkage), n.minobsinnode (Min. Terminal Node Size)
# glmnet     : Regularized Generalized Linear Models, package: glmnet , alpha (Mixing Percentage), lambda (Regularization Parameter)
# blackboost : Boosted Tree , package : mboost , mstop (#Trees), maxdepth (Max Tree Depth)
# nnet       : Neural Network , package : nnet  , size (#Hidden Units) , decay (Weight Decay)
# pcaNNet    : Neural Networks with Feature Extraction, package:nnet , size (#Hidden Units) , decay (Weight Decay)
# rpart      : CART, package:rpart, cp (Complexity Parameter)
# rf         : Random Forest,  package:randomForest, mtry (#Randomly Selected Predictors)


# Model names. For a list of available model names in caret package see: http://topepo.github.io/caret/available-models.html
# some available models given above
methods      <- c("glmnet", "gbm", "pcaNNet", "rf")   
method_names <- c("Elastic Net", "Boosting", "Nnet", "Random Forest")


# A list of arguments for models used in the estimation
args         <- list(svmLinear2=list(type='eps-regression'), svmLinear=list(type='nu-svr'), svmPoly=list(type='nu-svr'), gbm=list(verbose=FALSE), rf=list(ntree=1000), gamboost=list(baselearner='btree'), avNNet=list(verbose = 0, linout = TRUE, trace = FALSE), pcaNNet=list(linout = TRUE, trace = FALSE, MaxNWts=100000, maxit=10000), nnet=list(linout = TRUE, trace = FALSE, MaxNWts=100000, maxit=10000))


methodML   <- c("repeatedcv", "repeatedcv", "repeatedcv", "none")   # resampling method for chosing tuning parameters. available options: boot, boot632, cv, LOOCV, LGOCV, repeatedcv, oob, none
tune       <- c(100, 20, 20, NA)                                    # number of elements per parameter in the grid. the grid size is tune^{number of tuning parameters}. 
proces     <- c("range", "range","range", "range")                  # pre-processing method
select     <- c("best", "best","best", NA)                          # optimality criteria for choosing tuning parameter in cross validation. available options: best, oneSE, tolerance 
cv         <- c(2, 2,2, 2)                                          # the number of folds in cross-validation 
rep        <- c(2, 2,2, NA)                                         # number of iteration in repeated cross-validations 


# If there is a parameter of the model that user doesn't want to choose with cross validation, it should be set using tune_param variable. Below mtry of random forest is set to 5 
# for glmnet we want to choose both tuning parameters using cross validation so it is set to NULL

tune_param       <- list(0)
tune_param[[1]]  <- 0
tune_param[[2]]  <- 0
tune_param[[3]]  <- 0
tune_param[[4]]  <- data.frame(mtry=5)

output_name <- paste("range","-","best", "-", 2, "-" ,2, "-",sim,sep="")
name        <- "EL1"

####################################### Estimation  #######################################

r <- foreach(t = 1:sim, .combine='cbind', .inorder=FALSE, .packages=vec.pac) %dopar% { 
  
  set.seed(t);
  
  results       <- matrix(NA,5*length(Y), length(methods))
  results_het   <- matrix(NA,5*length(Y), length(methods))
  results_test  <- matrix(NA,15*length(Y), length(methods))
  results_group <- matrix(NA,15*length(Y), length(methods))
  table.who     <- matrix(NA, (length(affected)*15)*length(Y), length(methods))
  bestML        <- matrix(NA, 2*length(Y), length(methods))
  
  if(partition!="0"){
    ind <- createDataPartition(data[,partition], p = .5, list = FALSE)
    
    datause_raw <- as.data.frame(data[ ind,])
    dataout_raw <- as.data.frame(data[-ind,])
  }
  
  if(partition=="0"){
    split             <- runif(nrow(data))
    cvgroup           <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = TRUE))  
    
    datause_raw       <- as.data.frame(data[cvgroup == 1,])
    dataout_raw       <- as.data.frame(data[cvgroup != 1,])  
  }
  
  for(i in 1:length(Y)){
    
    y      <- Y[i]
    d      <- D[i]
    
    datause   <- data.frame(datause_raw[complete.cases(datause_raw[, c(controls, y, d, affected)]),])
    dataout   <- data.frame(dataout_raw[complete.cases(dataout_raw[, c(controls, y, d, affected)]),])
    
    ind_u <- which(datause[,d]==1)         # treatment indicator
    
    for(l in 1:length(methods)){
      
      if(methods[l]=="glmnet"){   x         <- XL     }
      if(methods[l]!="glmnet"){   x         <- X      }
      if(tune_param[[l]]==0){ f = NULL}
      if(tune_param[[l]]!=0){ f = tune_param[[l]]}
      
      form           <- as.formula(paste(y,"~",x,sep=""));
      
      ############ Estimate Scores using ML ############
      
      md_x         <- rep((nrow(datause[datause[,d]==1,]) + nrow(dataout[dataout[,d]==1,]))/(nrow(datause) + nrow(dataout)), nrow(dataout))  

      
      fitControl   <- trainControl(method = methodML[l], number = cv[l], repeats = rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=select[l])
      arg          <- c(list(form=form, data = datause[ind_u,],  method = methods[l],  tuneGrid = f, trControl = fitControl, preProcess=proces[l], tuneLength=tune[l]), args[[methods[l]]])
      fit.yz1      <- suppressWarnings(do.call(caret::train, arg))
      my_z1x       <- predict(fit.yz1, newdata=dataout, type="raw")
      
      fitControl   <- trainControl(method = methodML[l], number = cv[l], repeats = rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=select[l])
      arg          <- c(list(form=form, data = datause[-ind_u,],  method = methods[l], tuneGrid = f, trControl = fitControl, preProcess=proces[l], tuneLength=tune[l]), args[[methods[l]]])
      fit.yz0      <- suppressWarnings(do.call(caret::train, arg))
      my_z0x       <- predict(fit.yz0, newdata=dataout, type="raw")
      
      ind           <- (md_x>0.01 & md_x<0.99)
      dataout       <- dataout[ind, ]
      my_z1x        <- my_z1x[ind]
      my_z0x        <- my_z0x[ind]
      md_x          <- md_x[ind]
      
      ############################################# ATE by groups ############################################# 
      
      B   <- my_z0x
      S   <- (my_z1x - my_z0x)
      
      S2        <- S+runif(length(S), 0, 0.00001)
      breaks    <- quantile(S2, seq(0,1, 0.2),  include.lowest =T)
      breaks[1] <- breaks[1] - 0.001
      breaks[6] <- breaks[6] + 0.001
      SG        <- cut(S2, breaks = breaks)
      
      SGX       <- model.matrix(~-1+SG)
      DSG       <- data.frame(as.numeric(I(as.numeric(dataout[,d])-md_x))*SGX)
      
      colnames(DSG) <- c("G1", "G2", "G3", "G4", "G5")
      dataout[,c("B", "S", "G1", "G2", "G3", "G4", "G5", "weight")] <- cbind(B, S, DSG$G1, DSG$G2, DSG$G3, DSG$G4, DSG$G5, as.numeric((1/(md_x*(1-md_x)))))
      
      if(var(dataout$B)==0) {dataout$B <- dataout$B + rnorm(length(dataout$B),  mean=0, sd=0.1) }
      if(var(dataout$S)==0) {dataout$S <- dataout$S + rnorm(length(dataout$S),  mean=0, sd=0.1) }
      
      form1 <- as.formula(paste(y, "~", "B+S+G1+G2+G3+G4+G5 | ", fixed_effect, "| 0 |", cluster, sep=""))
      
      a <- tryCatch({
        a <- felm(form1, data=dataout, weights=dataout$weight)  
      },error=function(e){
        cat("ERROR :",methods[l], t, i, "\n")
        form1  <- as.formula(paste(y, "~", "G1+G2+G3+G4+G5 | ", fixed_effect, "| 0 |", cluster, sep=""))
        reg    <- felm(form1, data=dataout, weights=dataout$weight)  
        return(reg)
      }, warning = function(war) {
        cat("WARNING :",methods[l], t, i, "\n")
        form1  <- as.formula(paste(y, "~", "G1+G2+G3+G4+G5 | ", fixed_effect, "| 0 |", cluster, sep=""))
        reg    <- felm(form1, data=dataout, weights=dataout$weight)  
        return(reg)
      })
      reg   <- a
      
      
      coef <- (summary(reg)$coefficients['G5',1])
      pval <- (summary(reg)$coefficients['G5',4])
      results_test[(1+(i-1)*15):(5+((i-1)*15)),l]  <- c(coef, confint(reg, 'G5', level = 1-alpha)[1:2], (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)) )
      
      coef <- (summary(reg)$coefficients['G1',1])
      pval <- (summary(reg)$coefficients['G1',4])
      results_test[(6+(i-1)*15):(10+((i-1)*15)),l]  <- c(coef, confint(reg, 'G1', level = 1-alpha)[1:2], (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)) )
      
      test <- glht(reg, linfct = c("G5-G1==0"))
      coef <- (summary(reg)$coefficients['G5',1]) - (summary(reg)$coefficients['G1',1])
      pval <- summary(test)$test$pvalues[1]
      results_test[(11+(i-1)*15):(15+((i-1)*15)),l] <- c((confint(test,level = 1-alpha))$confint[1:3],(as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)))
      
      mean <- summary(reg)$coef[c('G1','G2','G3','G4','G5'),1]
      sd   <- summary(reg)$coef[c('G1','G2','G3','G4','G5'),2]
      
      crit <- qnorm(1-alpha/(p))
      
      results_group[((i-1)*15+1):((i-1)*15+5),l]   <- sort(mean)
      results_group[((i-1)*15+6):((i-1)*15+10),l]  <- sort(mean +crit*sd)
      results_group[((i-1)*15+11):((i-1)*15+15),l] <- sort(mean -crit*sd)
      
      bestML[(1+(i-1)*2),l]  <- (sum(mean^2)/5)
      
      ################################### Best Linear Prediction Regression  ################################### 
      
      Sd            <- dataout$S- mean(dataout$S)
      dataout$S_ort <- I((as.numeric(dataout[,d])-md_x)*Sd)
      dataout$d_ort <- I((as.numeric(dataout[,d])-md_x))
      
      form1 <- as.formula(paste(y, "~", "B+S+d_ort+S_ort| ", fixed_effect, "| 0 |", cluster, sep=""))
      
      a  <- tryCatch({
        a  <- felm(form1, data=dataout, weights=dataout$weight)   
      },error=function(e){
        cat("ERROR2 :",methods[l], t, i, "\n")
        form1 <- as.formula(paste(y, "~", "d_ort+S_ort| ", fixed_effect, "| 0 |", cluster, sep=""))
        reg   <- felm(form1, data=dataout, weights=dataout$weight)  
        return(reg)
      }, warning = function(war) {
        cat("WARNING2 :",methods[l], t, i, "\n")
        form1 <- as.formula(paste(y, "~", "d_ort+S_ort| ", fixed_effect, "| 0 |", cluster, sep=""))
        reg   <- felm(form1, data=dataout, weights=dataout$weight)  
        return(reg)
      })
      reg <- a 
      
      coef <- (summary(reg)$coefficients['d_ort',1])
      pval <- (summary(reg)$coefficients['d_ort',4])
      results[(1+(i-1)*5):(i*5),l]      <-c(coef, confint(reg, 'd_ort', level = 1-alpha)[1:2],  (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)))
      
      coef <- (summary(reg)$coefficients['S_ort',1])
      pval <- (summary(reg)$coefficients['S_ort',4])
      results_het[(1+(i-1)*5):(i*5),l] <- c(coef, confint(reg, 'S_ort', level = 1-alpha)[1:2],  (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)))
      bestML[(2+(i-1)*2),l]      <- abs(summary(reg)$coefficients['S_ort',1])*sqrt(var(dataout$S))
      
      
      ################################################ Most/Least Affected  ################################################ 
      
      high.effect     <- quantile(dataout$S, 1-thres);
      low.effect      <- quantile(dataout$S, thres);
      dataout$h       <- as.numeric(dataout$S>high.effect)
      dataout$l       <- as.numeric(dataout$S<low.effect)
      
      if(var(dataout$h)==0){dataout$h <- as.numeric(runif(length(dataout$h))<0.1)}
      if(var(dataout$l)==0){dataout$l <- as.numeric(runif(length(dataout$l))<0.1)}
      
      for(m in 1:length(affected)){
        a  <- tryCatch({
          form  <- paste(affected[m],"~h+l-1", sep="")
          reg   <- lm(form, data=dataout[(dataout$h==1)| (dataout$l==1),])    
          coef  <- reg$coefficients['h'] - reg$coefficients['l']
          test  <- glht(reg, linfct = c("h-l==0"))
          
          coef  <- (summary(reg)$coefficients['h',1])
          pval  <- (summary(reg)$coefficients['h',4])
          res1  <- c(coef, confint(reg, 'h', level = 1-alpha)[1:2], (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)))
          
          coef  <- (summary(reg)$coefficients['l',1])
          pval  <- (summary(reg)$coefficients['l',4])
          res2  <- c(coef, confint(reg, 'l', level = 1-alpha)[1:2], (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)))
          
          coef  <- (summary(reg)$coefficients['h',1]) - (summary(reg)$coefficients['l',1])
          pval  <- summary(test)$test$pvalues[1]
          res3  <- c((confint(test,level = 1-alpha))$confint[1:3], (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)))
          a     <- c(res1, res2, res3)
          
        },error=function(e){
          cat("ERROR3 :",methods[l], t, i, "\n")
          
          res1  <- c(mean(dataout[(dataout$h==1),affected[m]]), mean(dataout[(dataout$h==1),affected[m]]), mean(dataout[(dataout$h==1),affected[m]]), 0.5, 0.5)
          res2  <- c(mean(dataout[(dataout$l==1),affected[m]]), mean(dataout[(dataout$l==1),affected[m]]), mean(dataout[(dataout$l==1),affected[m]]), 0.5, 0.5)
          res3  <- c((res1[1] - res2[1]), (res1[1] - res2[1]), (res1[1] - res2[1]), 0.5 , 0.5)
          a     <- c(res1, res2, res3)
          return(a)
        })
        table.who[((i-1)*length(affected)*15+(m-1)*15+1):((i-1)*length(affected)*15+(m)*15),l]   <- a
      }
    }
  }
  res <- c(as.vector(results_test), as.vector(results), as.vector(results_het), as.vector(results_group), as.vector(bestML), as.vector(table.who))
  print(t)
  r <- data.frame(res)
}

####################################### Prepare Latex Tables for BLP and Most/Least Affected  #######################################

results_test  <- array(c(as.matrix(r[1:(15*length(Y)*length(methods)),])), c(15*length(Y),length(methods), sim))
results       <- array(c(as.matrix(r[((15*length(Y)*length(methods))+1):((15+5)*length(Y)*length(methods)),])), c(5*length(Y),length(methods), sim))
results_het   <- array(c(as.matrix(r[(((20)*length(Y)*length(methods))+1):((20+5)*length(Y)*length(methods)),])), c(5*length(Y),length(methods), sim))
results_group <- array(c(as.matrix(r[(((25)*length(Y)*length(methods))+1):((25+15)*length(Y)*length(methods)),])), c(15*length(Y),length(methods), sim))
bestML        <- array(c(as.matrix(r[(((40)*length(Y)*length(methods))+1):((40+2)*length(Y)*length(methods)),])), c(2*length(Y),length(methods), sim))
table.who     <- array(c(as.matrix(r[(((42)*length(Y)*length(methods))+1):((42+length(affected)*15)*length(Y)*length(methods)),])), c(length(affected)*15*length(Y),length(methods), sim))

results_all       <- t(sapply(seq(1:nrow(results[,,1])), function(x) colMedians(t(results[x,,]))))
results_het_all   <- t(sapply(seq(1:nrow(results_het[,,1])), function(x) colMedians(t(results_het[x,,]))))
test_all          <- t(sapply(seq(1:nrow(results_test[,,1])), function(x) colMedians(t(results_test[x,,]))))
table.who_all     <- t(sapply(seq(1:nrow(table.who[,,1])), function(x) colMedians(t(table.who[x,,]))))
group_all         <- t(sapply(seq(1:nrow(results_group[,,1])), function(x) colMedians(t(results_group[x,,]))))
bestML_all        <- t(sapply(seq(1:nrow(bestML[,,1])), function(x) colMedians(t(bestML[x,,]))))

best1 <- order(-bestML_all[1,])[1:2]
best2 <- order(-bestML_all[2,])[1:2]

if(best1[1]!=best2[1]){ best <- c(best1[1],best2[1])}
if(best1[1]==best2[1]){ best <- c(best2[1],best2[2])}


table.who_all2   <- matrix(0,length(affected)*12*length(Y), length(methods))
results_all2     <- matrix(NA,4*length(Y), length(methods))
results_het_all2 <- matrix(NA,4*length(Y), length(methods))
test_all2        <- matrix(NA,12*length(Y), length(methods))

l <- 1

for(i in seq(1, nrow(table.who_all), 5)){
  
  table.who_all2[l:(l+2),]  <- table.who_all[i:(i+2),]
  table.who_all2[l+3,]      <- sapply(seq(1:length(methods)), function(x) min(1,4*min(table.who_all[i+3,x], table.who_all[i+4,x])))
  
  if(l<nrow(results_all2)){
    results_all2[l:(l+2),]       <- results_all[i:(i+2),]
    results_het_all2[l:(l+2),]   <- results_het_all[i:(i+2),]
    
    results_all2[l+3,]     <- sapply(seq(1:length(methods)), function(x) min(1,4*min(results_all[i+3,x], results_all[i+4,x])))
    results_het_all2[l+3,] <- sapply(seq(1:length(methods)), function(x) min(1,4*min(results_het_all[i+3,x], results_het_all[i+4,x])))
  }
  if(l<nrow(test_all2)){
    
    test_all2[l:(l+2),]          <- test_all[i:(i+2),]
    test_all2[l+3,]              <- sapply(seq(1:length(methods)), function(x) min(1,4*min(test_all[i+3,x], test_all[i+4,x])))
  }
  
  l <- l+4
}


results_all     <- round(results_all2, digits = 3)
results_het_all <- round(results_het_all2, digits = 3)
test_all        <- round(test_all2, digits = 3)
bestML_all      <- format(round(bestML_all, pmax(0,4-nchar(floor(abs(bestML_all))))), nsmall= pmax(0,4-nchar(floor(abs(bestML_all)))))
table.who_all   <- round(table.who_all2, digits = 3)

results_test2   <- matrix(0,9*length(Y), length(methods))
result2         <- matrix(0,3*length(Y), length(methods))
result_het2     <- matrix(0,3*length(Y), length(methods))
table.who_all2  <- matrix(0,9*length(Y)*length(affected), length(methods))

seq3 <- seq(1, nrow(table.who_all), 4)
l    <- 1

for(i in seq(1, nrow(table.who_all2), 3)){
  
  k <- seq3[l]
  
  if(i<nrow(result2)){
    result2[i,]       <- format(round(results_all[k,], pmax(0,4-nchar(floor(abs(results_all[k,]))))), nsmall= pmax(0,4-nchar(floor(abs(results_all[k,])))))
    result2[i+1,]     <- sapply(seq(1:ncol(results_all)), function(x) paste("(", format(round(results_all[k+1,x],pmax(0,4-nchar(floor(abs(results_all[k+1,x]))))), nsmall=pmax(0,4-nchar(floor(abs(results_all[k+1,x]))))), ",", format(round(results_all[k+2,x],pmax(0,4-nchar(floor(abs(results_all[k+2,x]))))) , nsmall=pmax(0,4-nchar(floor(abs(results_all[k+2,x]))))), ")", sep=""))
    result2[i+2,]     <- paste("[", format(results_all[k+3,], nsmall = pmax(0,4-nchar(floor(abs(results_all[k+3,]))))), "]", sep="")
    
    result_het2[i,]   <- format(round(results_het_all[k,],max(0,4-nchar(floor(abs(results_het_all[k,]))))) , nsmall=pmax(0,4-nchar(floor(results_het_all[k,]))))
    result_het2[i+1,] <- sapply(seq(1:ncol(results_het_all)), function(x) paste("(", format(round(results_het_all[k+1,x], pmax(0,4-nchar(floor(abs(results_het_all[k+1,x]))))) , nsmall=pmax(0,4-nchar(floor(abs(results_het_all[k+1,x]))))), ",", format(round(results_het_all[k+2,x],pmax(0,4-nchar(floor(abs(results_het_all[k+2,x]))))) , nsmall=pmax(0,4-nchar(floor(abs(results_het_all[k+2,x]))))), ")", sep=""))
    result_het2[i+2,] <- paste("[", format(results_het_all[k+3,], nsmall=max(0,4-nchar(floor(abs(results_het_all[k+3,]))))), "]", sep="")
  }
  
  if(i<nrow(results_test2)){
    results_test2[i,]    <- format(round(test_all[k,],pmax(0,4-nchar(floor(abs(test_all[k,]))))) , nsmall=pmax(0,4-nchar(floor(test_all[k,]))))
    results_test2[i+1,]  <- sapply(seq(1:ncol(test_all)), function(x) paste("(", format(round(test_all[k+1,x], pmax(0,4-nchar(floor(abs(test_all[k+1,x]))))),nsmall=pmax(0,4-nchar(floor(abs(test_all[k+1,x]))))), ",", format(round(test_all[k+2,x],pmax(0,4-nchar(abs(floor(test_all[k+2,x]))))),  nsmall=pmax(0,4-nchar(floor(abs(test_all[k+2,x]))))), ")", sep=""))
    results_test2[i+2,]  <- paste("[", format(test_all[k+3,], nsmall=pmax(0,4-nchar(floor(abs(test_all[k+3,]))))), "]", sep="")
  }
  
  table.who_all2[i,]       <- format(round(table.who_all[k,],pmax(0,4-nchar(floor(abs(table.who_all[k,]))))) ,nsmall=pmax(0,4-nchar(floor(abs(table.who_all[k,])))))
  table.who_all2[i+1,]     <- sapply(seq(1:ncol(table.who_all)), function(x) paste("(", format(round(table.who_all[k+1,x], pmax(0,4-nchar(floor(abs(table.who_all[k+1,x]))))), nsmall=pmax(0,4-nchar(floor(abs(table.who_all[k+1,x]))))), ",", format(round(table.who_all[k+2,x],pmax(0,4-nchar(floor(abs(table.who_all[k+2,x]))))) , nsmall=pmax(0,4-nchar(floor(abs(table.who_all[k+2,x]))))), ")", sep=""))
  if(i%%9==7){  table.who_all2[i+2,]     <- paste("[", format(table.who_all[k+3,], nsmall=pmax(0,4-nchar(floor(abs(table.who_all[k+3,]))))), "]", sep="") }
  if(i%%9!=7){  table.who_all2[i+2,]     <- "-" }
  
  
  l <- l+1
}


CLAN_final         <- matrix(NA, length(Y)*(length(affected)*3+1), length(best)*3)
GATES_final        <- matrix(NA, length(Y)*3, length(best)*3)
BLP_final          <- matrix(NA, length(Y)*3, length(best)*2)
BEST_final         <- bestML_all

rownames_CLAN      <- matrix(NA, nrow(CLAN_final),1)
rownames_GATES     <- matrix(NA, nrow(GATES_final),1)
rownames_BEST      <- matrix(NA, nrow(bestML_all),1)

a  <- 1
b  <- 1
c  <- 1
c2 <- 1

for(l in 1:length(Y)){
  
  rownames_CLAN[a] <- names[l]
  
  a <- a+1
  
  for(i in 1:length(affected)){
    for(j in 1:length(best)){
      
      k <- best[j]
      
      CLAN_final[(a):(a+2),((j-1)*3+1):(j*3)] <- matrix(table.who_all2[(b):(b+8),k], 3, 3)
      
      if(i==1){
        GATES_final[(c):(c+2),((j-1)*3+1):(j*3)] <- matrix(results_test2[(c2):(c2+8),k], 3, 3)
        rownames_GATES[c]   <- names[l]
        BLP_final[(c):(c+2),((j-1)*2+1):(j*2)] <- cbind(result2[(c):(c+2),j], result_het2[(c):(c+2),k])
      }
      
      rownames_CLAN[a]   <- names_affected[i]
    }
    a <- a+3
    b <- b+9
  }
  c  <- c+3
  c2 <- c2+9
  
  rownames_BEST[((l-1)*2+1):((l-1)*2+2)] <- c(names[l], names[l])
  
}

rownames(CLAN_final)   <- rownames_CLAN
rownames(GATES_final)  <- rownames_GATES
rownames(BLP_final)    <- rownames_GATES
rownames(BEST_final)   <- rownames_BEST

colnames(CLAN_final)   <- rep(c("Most Affected", 	"Least Affected",	"Difference"), length(best))
colnames(GATES_final)  <- rep(c("Most Affected", 	"Least Affected",	"Difference"), length(best))
colnames(BLP_final)    <- rep(c("ATE", 	"HET"), length(best))
colnames(BEST_final)   <- method_names

print(xtable(cbind(rownames(BLP_final),BLP_final)), include.rownames=FALSE,file=paste(name,"_BLP","-",output_name,".txt",sep=""), digits=3)
print(xtable(cbind(rownames(GATES_final),GATES_final)), include.rownames=FALSE,file=paste(name,"_GATES","-",output_name,".txt",sep=""), digits=3)
print(xtable(cbind(rownames(BEST_final),BEST_final)), include.rownames=FALSE,file=paste(name,"_BEST","-",output_name,".txt",sep=""), digits=3)
print(xtable(cbind(rownames(CLAN_final),CLAN_final)), include.rownames=FALSE,file=paste(name,"_CLAN","-",output_name,".txt",sep=""), digits=3)


for(i in 1:length(Y)){
  if(length(methods)>1){
    par(mfrow=c(2,2))
  }
  
  y_range     <- 1*range(group_all[(15*(i-1)+6):(15*(i-1)+10),],group_all[(15*(i-1)+11):(15*(i-1)+15),])
  y_range2    <- y_range
  y_range2[1] <- y_range[1]- (y_range[2] -  y_range[1])*0.1
  y_range2[2] <- y_range[2]+ (y_range[2] -  y_range[1])*0.1
  
  result=list(0)
  
  for(j in 1:length(methods)){
    
    ATE <- data.frame( x = c(-Inf, Inf), y = results_all[(4*(i-1)+1),j] , cutoff = factor(50))
    U   <- data.frame( x = c(-Inf, Inf), y = results_all[(4*(i-1)+3),j] , cutoff = factor(50))
    L   <- data.frame( x = c(-Inf, Inf), y = results_all[(4*(i-1)+2),j] , cutoff = factor(50))
    
    df <- data.frame(x =1:5,
                     F =group_all[(15*(i-1)+1):(15*(i-1)+5),j],
                     L =group_all[(15*(i-1)+6):(15*(i-1)+10),j],
                     U =group_all[(15*(i-1)+11):(15*(i-1)+15),j],
                     group = factor(c(2, 2, 2, 2,2)))
    
    result[[j]] <- ggplot() +
      theme_gray(base_size = 14) +
      geom_point(data=df,aes(y = F, x = x, colour='90% CB(GATES)'), size = 3) +
      geom_errorbar(data=df, aes(ymax = U, ymin = L ,x=x, y=F, height = .2, width=0.7, colour="GATES"), show.legend = TRUE) +
      geom_line(aes( x, y, linetype = cutoff, colour='ATE' ),ATE, linetype = 2) +
      geom_line(aes( x, y, linetype = cutoff, colour='90% CB(ATE)' ), U, linetype = 2) +
      geom_line(aes( x, y, linetype = cutoff ), L, linetype = 2, color="red") +
      scale_colour_manual(values = c("red", "black", "blue", "black"),
                          breaks=c('ATE','90% CB(ATE)',"GATES",'90% CB(GATES)'),
                          guide = guide_legend(override.aes = list(
                            linetype = c("dashed", "dashed"  ,"blank", "solid"),
                            shape = c(NA,NA, 16, NA)), ncol =2,byrow=TRUE)) +
      theme(plot.title = element_text(hjust = 0.5,size = 11, face = "bold"), axis.title=element_text(size=10), legend.text=element_text(size=7), legend.key = element_rect(colour = NA, fill = NA), legend.key.size = unit(1, 'lines'), legend.title=element_blank(),legend.justification=c(0,1), legend.position=c(0,1), legend.background=element_rect(fill=alpha('blue', 0)))  +
      ylim(y_range) +
      labs(title=method_names[j], y = "Treatment Effect", x = "Group by Het Score") 
    
  }
  print(Y[i])
  p      <- plot_grid(result[[1]], result[[2]], result[[3]], result[[4]], ncol=2)
  ggsave(paste(name,"_plot","-",output_name,"-",Y[i],".pdf",sep=""), p)
  p_best <- plot_grid(result[[best[1]]], result[[best[2]]], ncol=2)
  ggsave(paste(name,"_plot_best","-",output_name,"-",Y[i],".pdf",sep=""), p_best, width = 10, height = 5)
}

