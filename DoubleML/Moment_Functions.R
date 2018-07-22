#########################################################################################################
#  Program:     Functions for estimating moments for using Machine Learning Methods.                    #
#  References: "Double/Debiased Machine Learning of Treatment and Causal Parameters",  AER P&P 2017     #
#              "Double Machine Learning for Treatment and Causal Parameters",  Arxiv 2016               #
#  by V.Chernozhukov, D. Chetverikov, M. Demirer, E. Duflo, C. Hansen, W. Newey                         #           
#########################################################################################################

source("ML_Functions.R") 

DoubleML <- function(data, y, d,z, xx, xL, methods, DML, nfold, est, args, silent=FALSE, trim, param, tune_param){
  
  if(nfold==1){  print("nfold should be greater than 1")  
    stop() 
  } 
  
  TE        <- matrix(0,1,(length(methods)+1))
  STE       <- matrix(0,1,(length(methods)+1))
  result    <- matrix(0,2,(length(methods)+1))
  result2   <- matrix(0,2,(length(methods)+1))
  cond.comp <- matrix(list(),length(methods),nfold)
  
  for(i in 1:5){
    assign(paste("MSE", i, sep = ""), matrix(0,length(methods)+1,nfold))    
  }
  
  for(i in c("d","y","zx","z","yz1","yz0","dz1","dz0", "dx", "yd1", "yd0")){
    assign(paste(i, "pool", sep = ""), vector("list", (length(methods)+1)))    
  }
  
  binary    <- as.numeric(checkBinary(data[,d]))
  
  if(est=="LATE"){
    
    binary.z <- as.numeric(checkBinary(data[,z])) 
    if(!(binary.z==1)){
      print("instrument is not binary")
      stop()
    } 
    flag    <- if(sum(!(data[data[,z]==0,d]==0))==0) 1 else 0
  }
  
  split <- runif(nrow(data))
  
  cvgroup   <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/nfold)),include.lowest = TRUE))  
  
  for(k in 1:length(methods)){   
    
    if(silent==FALSE){  cat(methods[k],'\n')}
    x <- if(any(c("glmnet")==methods[k])) xL else xx
    
    for(j in 1:nfold){   

      if(silent==FALSE){ cat('  fold',j,'\n')  }
    
      datause = as.data.frame(data[cvgroup != j,])  
      dataout = as.data.frame(data[cvgroup == j,])
      
      if(est=="LATE" && (length(methods)>0)){
        
        cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, x=x, z=z, method=k, methods=methods, plinear=est, xL=xL, binary=binary, flag=flag, args=args, param=param, tune_param=tune_param)
        
        MSE1[k,j]               <- cond.comp[[k,j]]$err.yz0
        MSE2[k,j]               <- cond.comp[[k,j]]$err.yz1
        MSE3[k,j]               <- if(flag==1) 0 else cond.comp[[k,j]]$err.dz0
        MSE4[k,j]               <- cond.comp[[k,j]]$err.dz1
        MSE5[k,j]               <- cond.comp[[k,j]]$err.z
        
        drop                   <- which(cond.comp[[k,j]]$mz_x>trim[1] & cond.comp[[k,j]]$mz_x<trim[2])      
        mz_x                   <- cond.comp[[k,j]]$mz_x[drop]
        my_z1x                 <- cond.comp[[k,j]]$my_z1x[drop]
        my_z0x                 <- cond.comp[[k,j]]$my_z0x[drop]
        md_z1x                 <- cond.comp[[k,j]]$md_z1x[drop]
        md_z0x                 <- if(flag==1) matrix(0,1,length(my_z0x)) else cond.comp[[k,j]]$md_z0x[drop] 
        yout                   <- dataout[drop,y]
        dout                   <- dataout[drop,d]
        zout                   <- dataout[drop,z]
        
        TE[1,k]                <- LATE(yout, dout, zout, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)/nfold + TE[1,k];
        STE[1,k]               <- (1/(nfold^2))*((SE.LATE(yout, dout, zout, my_z1x, my_z0x, mz_x, md_z1x, md_z0x))^2) + STE[1,k];
        
        ypool[[k]]             <- c(ypool[[k]], yout)
        dpool[[k]]             <- c(dpool[[k]], dout)
        zpool[[k]]             <- c(zpool[[k]], zout)
        zxpool[[k]]            <- c(zxpool[[k]], mz_x)
        yz1pool[[k]]           <- c(yz1pool[[k]], my_z1x)
        yz0pool[[k]]           <- c(yz0pool[[k]], my_z0x)
        dz1pool[[k]]           <- c(dz1pool[[k]], md_z1x)
        dz0pool[[k]]           <- c(dz0pool[[k]], md_z0x)
        
        MSE1[(length(methods)+1),j] <- error(mean(datause[datause[,z]==0,y], na.rm = TRUE), dataout[!is.na(dataout[dataout[,z]==0,y]),y])$err
        MSE2[(length(methods)+1),j] <- error(mean(datause[datause[,z]==1,y], na.rm = TRUE), dataout[!is.na(dataout[dataout[,z]==1,y]),y])$err
        MSE3[(length(methods)+1),j] <- if(flag==1) MSE3[(length(methods)+1),j]=0  else error(mean(datause[datause[,z]==0,d], na.rm = TRUE), dataout[!is.na(dataout[dataout[,z]==0,d]),d])$err
        MSE4[(length(methods)+1),j] <- error(mean(datause[datause[,z]==1,d], na.rm = TRUE), dataout[!is.na(dataout[dataout[,z]==1,d]),d])$err
        MSE5[(length(methods)+1),j] <- error(mean(datause[,z], na.rm = TRUE), dataout[!is.na(dataout[,z]),z])$err
      }
      
      if(est=="interactive" && (length(methods)>0)){
        
        if(methods[k]=="Ensemble") { cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y=y, d=d, x=x, method=k, methods=methods,  plinear=est, xL=xL, binary=binary, arguments=arguments, ensemble=ensemble)}
        else{                       cond.comp[[k,j]]  <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, x=x, method=k, methods=methods,  plinear=est, xL=xL, binary=binary, args=args, param=param, tune_param=tune_param);  }

        MSE1[k,j]               <- cond.comp[[k,j]]$err.yd0
        MSE2[k,j]               <- cond.comp[[k,j]]$err.yd1
        MSE3[k,j]               <- cond.comp[[k,j]]$err.d
      
        drop                   <- which(cond.comp[[k,j]]$md_x>trim[1] & cond.comp[[k,j]]$md_x<trim[2])      
        md_x                   <- cond.comp[[k,j]]$md_x[drop]
        my_d1x                 <- cond.comp[[k,j]]$my_d1x[drop]
        my_d0x                 <- cond.comp[[k,j]]$my_d0x[drop]
        yout                   <- dataout[drop,y]
        dout                   <- dataout[drop,d]

        TE[1,k]                <- ATE(yout, dout, my_d1x, my_d0x, md_x)/nfold + TE[1,k];
        STE[1,k]               <- (1/(nfold^2))*((SE.ATE(yout, dout, my_d1x, my_d0x, md_x))^2) + STE[1,k];
        
        ypool[[k]]             <- c(ypool[[k]], yout)
        dpool[[k]]             <- c(dpool[[k]], dout)
        dxpool[[k]]            <- c(zpool[[k]], md_x)
        yd1pool[[k]]           <- c(yd1pool[[k]], my_d1x)
        yd0pool[[k]]           <- c(yd0pool[[k]], my_d0x)

        MSE1[(length(methods)+1),j] <- error(mean(datause[datause[,d]==0,y], na.rm = TRUE), dataout[!is.na(dataout[dataout[,d]==0,y]),y])$err
        MSE2[(length(methods)+1),j] <- error(mean(datause[datause[,d]==1,y], na.rm = TRUE), dataout[!is.na(dataout[dataout[,d]==1,y]),y])$err
        MSE3[(length(methods)+1),j] <- error(mean(datause[,d], na.rm = TRUE), dataout[!is.na(dataout[,d]),d])$err
        
      }
      
      if(est=="plinear" && (length(methods)>0)){
        
        if(methods[k]=="Ensemble") { cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y=y, d=d, x=x, method=k, methods=methods,  plinear=est, xL=xL, binary=binary, arguments=arguments, ensemble=ensemble)}
        else{                        cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, x=x, z=z, method=k, methods=methods,  plinear=est, xL=xL, binary=binary, flag=flag, args=args, param=param, tune_param=tune_param)}

        MSE1[k,j]              <- cond.comp[[k,j]]$err.y
        MSE2[k,j]              <- cond.comp[[k,j]]$err.d

        lm.fit.ry              <- lm(as.matrix(cond.comp[[k,j]]$ry) ~ as.matrix(cond.comp[[k,j]]$rd)-1);
        ate                    <- lm.fit.ry$coef;
        HCV.coefs              <- vcovHC(lm.fit.ry, type = 'HC');
        STE[1,k]               <- (1/(nfold^2))*(diag(HCV.coefs)) +  STE[1,k] 
        TE[1,k]                <- ate/nfold + TE[1,k] ;
        ypool[[k]]             <- c(ypool[[k]], cond.comp[[k,j]]$ry)
        dpool[[k]]             <- c(dpool[[k]], cond.comp[[k,j]]$rd)
        
        
        MSE1[(length(methods)+1),j] <- error(rep(mean(datause[,y], na.rm = TRUE), length(dataout[!is.na(dataout[,y]),y])), dataout[!is.na(dataout[,y]),y])$err
        MSE2[(length(methods)+1),j] <- error(rep(mean(datause[,d], na.rm = TRUE), length(dataout[!is.na(dataout[,d]),d])), dataout[!is.na(dataout[,d]),d])$err
      }
      
      if(est=="IV" && (length(methods)>0)){
        
        if(methods[k]=="Ensemble") { 
          cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y=y, d=d, z=z, x=x, method=k, methods=methods,  plinear=est, xL=xL, binary=binary, arguments=arguments, ensemble=ensemble)
        }
        else{       
          cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, z=z, x=x, method=k, methods=methods,  plinear=est, xL=xL, binary=binary, args=args, param=param, tune_param=tune_param)
        }
        
        MSE1[k,j]              <- cond.comp[[k,j]]$err.y
        MSE2[k,j]              <- cond.comp[[k,j]]$err.d
        MSE3[k,j]              <- cond.comp[[k,j]]$err.z
        
        lm.fit.ry              <- tsls(y=cond.comp[[k,j]]$ry,d=cond.comp[[k,j]]$rd, x=NULL, z=cond.comp[[k,j]]$rz, intercept = FALSE)
        ate                    <- lm.fit.ry$coef[1];
        HCV.coefs              <- lm.fit.ry$vcov[1]
        
        
        STE[1,k]               <- (1/(nfold^2))*((HCV.coefs)) +  STE[1,k] 
        TE[1,k]                <- ate/nfold + TE[1,k] ;
        
        ypool[[k]]             <- c(ypool[[k]], cond.comp[[k,j]]$ry)
        zpool[[k]]             <- c(zpool[[k]], cond.comp[[k,j]]$rz)
        dpool[[k]]             <- c(dpool[[k]], cond.comp[[k,j]]$rd)
        
        MSE1[(length(methods)+1),j] <- NA
        MSE2[(length(methods)+1),j] <- NA
        MSE3[(length(methods)+1),j] <- NA
        
      }
    }  
  }
  if(est=="LATE"){
    
    min1 <- if(length(methods)>1) which.min(rowMeans(MSE1[1:length(methods),])) else which.min(mean(MSE1[1:length(methods),]))
    min2 <- if(length(methods)>1) which.min(rowMeans(MSE2[1:length(methods),])) else which.min(mean(MSE2[1:length(methods),]))
    min3 <- if(length(methods)>1) which.min(rowMeans(MSE3[1:length(methods),])) else which.min(mean(MSE3[1:length(methods),]))
    min4 <- if(length(methods)>1) which.min(rowMeans(MSE4[1:length(methods),])) else which.min(mean(MSE4[1:length(methods),]))
    min5 <- if(length(methods)>1) which.min(rowMeans(MSE5[1:length(methods),])) else which.min(mean(MSE5[1:length(methods),]))
    
    if(silent==FALSE){
      cat('  best methods for E[Y|X, D=0]:',methods[min1],'\n')
      cat('  best methods for E[Y|X, D=1]:',methods[min2],'\n')
      cat('  best methods for E[D|X, Z=0]:',methods[min3],'\n')
      cat('  best methods for E[D|X, Z=1]:',methods[min4],'\n')
      cat('  best methods for E[Z|X]:',methods[min5],'\n')
    }
  }
  
  if(est=="interactive"){
    
    min1 <- if(length(methods)>1) which.min(rowMeans(MSE1[1:length(methods),])) else which.min(mean(MSE1[1:length(methods),]))
    min2 <- if(length(methods)>1) which.min(rowMeans(MSE2[1:length(methods),])) else which.min(mean(MSE2[1:length(methods),]))
    min3 <- if(length(methods)>1) which.min(rowMeans(MSE3[1:length(methods),])) else which.min(mean(MSE3[1:length(methods),]))
    
    if(silent==FALSE){
      cat('  best methods for E[Y|X, D=0]:',methods[min1],'\n')
      cat('  best methods for E[Y|X, D=1]:',methods[min2],'\n')
      cat('  best methods for E[D|X]:',methods[min3],'\n')
    }
  }
  
  if(est=="plinear"){
    
    min1 <- if(length(methods)>1) which.min(rowMeans(MSE1[1:length(methods),])) else which.min(mean(MSE1[1:length(methods),]))
    min2 <- if(length(methods)>1) which.min(rowMeans(MSE2[1:length(methods),])) else which.min(mean(MSE2[1:length(methods),]))
    
    if(silent==FALSE){   
      cat('  best methods for E[Y|X]:',methods[min1],'\n')
      cat('  best methods for E[D|X]:',methods[min2],'\n')
    }    
  }
  
  if(est=="IV"){
    
    min1 <- if(length(methods)>1) which.min(rowMeans(MSE1[1:length(methods),])) else which.min(mean(MSE1[1:length(methods),]))
    min2 <- if(length(methods)>1) which.min(rowMeans(MSE2[1:length(methods),])) else which.min(mean(MSE2[1:length(methods),]))
    min3 <- if(length(methods)>1) which.min(rowMeans(MSE3[1:length(methods),])) else which.min(mean(MSE3[1:length(methods),]))
    
    if(silent==FALSE){   
      cat('  best methods for E[Y|X]:',methods[min1],'\n')
      cat('  best methods for E[D|X]:',methods[min2],'\n')
      cat('  best methods for E[Z|X]:',methods[min3],'\n')
    }    
  }
  
  for(j in 1:nfold){  
    
    datause = as.data.frame(data[cvgroup != j,])  
    dataout = as.data.frame(data[cvgroup == j,])
    
    if(est=="LATE"){
      
      drop                   <- which(cond.comp[[min5,j]]$mz_x>trim[1] & cond.comp[[min5,j]]$mz_x<trim[2])      
      mz_x                   <- cond.comp[[min1,j]]$mz_x[drop]
      my_z1x                 <- cond.comp[[min2,j]]$my_z1x[drop]
      my_z0x                 <- cond.comp[[min3,j]]$my_z0x[drop]
      md_z1x                 <- cond.comp[[min4,j]]$md_z1x[drop]
      if(flag==1){ md_z0x    <- matrix(0,1,length(my_z0x))}
      else{  md_z0x          <- cond.comp[[min5,j]]$md_z0x[drop] }
      
      yout                   <- dataout[drop,y]
      dout                   <- dataout[drop,d]
      zout                   <- dataout[drop,z]
      
      TE[1,(k+1)]            <- LATE(yout, dout, zout, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)/nfold + TE[1,(k+1)];
      STE[1,(k+1)]           <- (1/(nfold^2))*((SE.LATE(yout, dout, zout, my_z1x, my_z0x, mz_x, md_z1x, md_z0x))^2) + STE[1,(k+1)];
      
      ypool[[k+1]]             <- c(ypool[[k+1]], yout)
      dpool[[k+1]]             <- c(dpool[[k+1]], dout)
      zpool[[k+1]]             <- c(zpool[[k+1]], zout)
      zxpool[[k+1]]            <- c(zxpool[[k+1]], mz_x)
      yz1pool[[k+1]]           <- c(yz1pool[[k+1]], my_z1x)
      yz0pool[[k+1]]           <- c(yz0pool[[k+1]], my_z0x)
      dz1pool[[k+1]]           <- c(dz1pool[[k+1]], md_z1x)
      dz0pool[[k+1]]           <- c(dz0pool[[k+1]], md_z0x)
      
    }
    
    
    if(est=="interactive"){
      
      drop                   <- which(cond.comp[[min3,j]]$md_x>trim[1] & cond.comp[[min3,j]]$md_x<trim[2])      
      md_x                   <- cond.comp[[min1,j]]$md_x[drop]
      my_d1x                 <- cond.comp[[min2,j]]$my_d1x[drop]
      my_d0x                 <- cond.comp[[min3,j]]$my_d0x[drop]
      yout                   <- dataout[drop,y]
      dout                   <- dataout[drop,d]
      
      TE[1,(k+1)]            <- ATE(yout, dout, my_d1x, my_d0x, md_x)/nfold + TE[1,(k+1)];
      STE[1,(k+1)]           <- (1/(nfold^2))*((SE.ATE(yout, dout, my_d1x, my_d0x, md_x))^2) + STE[1,(k+1)];

      ypool[[k+1]]             <- c(ypool[[k+1]], yout)
      dpool[[k+1]]             <- c(dpool[[k+1]], dout)
      dxpool[[k+1]]            <- c(dxpool[[k+1]], md_x)
      yd1pool[[k+1]]           <- c(yd1pool[[k+1]], my_d1x)
      yd0pool[[k+1]]          <- c(yd0pool[[k+1]], my_d0x)
      
    }
    
    if(est=="plinear"){

      lm.fit.ry              <- lm(as.matrix(cond.comp[[min1,j]]$ry) ~ as.matrix(cond.comp[[min2,j]]$rd)-1);
      ate                    <- lm.fit.ry$coef;
      HCV.coefs              <- vcovHC(lm.fit.ry, type = 'HC');
      STE[1,(k+1)]           <- (1/(nfold^2))*(diag(HCV.coefs)) +  STE[1,(k+1)] 
      TE[1,(k+1)]            <- ate/nfold + TE[1,(k+1)] ;
      
      ypool[[k+1]]             <- c(ypool[[k+1]], cond.comp[[min1,j]]$ry)
      dpool[[k+1]]             <- c(dpool[[k+1]], cond.comp[[min2,j]]$rd)
    }
    
    if(est=="IV"){
      
      lm.fit.ry              <- tsls(y=cond.comp[[min1,j]]$ry, d=cond.comp[[min2,j]]$rd, x=NULL, z=cond.comp[[min3,j]]$rz, intercept = FALSE)
      ate                    <- lm.fit.ry$coef[1];
      HCV.coefs              <- lm.fit.ry$vcov[1]
      
      
      STE[1,(k+1)]           <- (1/(nfold^2))*((HCV.coefs)) +  STE[1,(k+1)] 
      TE[1,(k+1)]            <- ate/nfold + TE[1,(k+1)] ;
      
      ypool[[k+1]]             <- c(ypool[[k+1]], cond.comp[[min1,j]]$ry)
      zpool[[k+1]]             <- c(zpool[[k+1]], cond.comp[[min3,j]]$rz)
      dpool[[k+1]]             <- c(dpool[[k+1]], cond.comp[[min2,j]]$rd)
      
    }
  }
  
  TE_pool        <- matrix(0,1,(length(methods)+1))
  STE_pool       <- matrix(0,1,(length(methods)+1))
  
  for(k in 1:(length(methods)+1)){ 
    
    if(est=="LATE"){
      
      TE_pool[1,(k)]         <- LATE(ypool[[k]], dpool[[k]], zpool[[k]],  yz1pool[[k]], yz0pool[[k]], zxpool[[k]], dz1pool[[k]], dz0pool[[k]])
      STE_pool[1,(k)]        <- ((SE.LATE(ypool[[k]], dpool[[k]], zpool[[k]],  yz1pool[[k]], yz0pool[[k]], zxpool[[k]], dz1pool[[k]], dz0pool[[k]]))^2) 
      
    }
    
    if(est=="interactive"){
      
      TE_pool[1,(k)]         <- ATE(ypool[[k]], dpool[[k]], yd1pool[[k]], yd0pool[[k]], dxpool[[k]])
      STE_pool[1,(k)]        <- ((SE.ATE(ypool[[k]], dpool[[k]], yd1pool[[k]], yd0pool[[k]], dxpool[[k]]))^2)
      
    }
    
    if(est=="plinear"){

      lm.fit.ry              <- lm(as.matrix(ypool[[k]]) ~ as.matrix(dpool[[k]])-1);
      ate                    <- lm.fit.ry$coef;
      HCV.coefs              <- vcovHC(lm.fit.ry, type = 'HC');
      STE_pool[1,k]          <- (diag(HCV.coefs))
      TE_pool[1,k]           <- ate

    }
    
    if(est=="IV"){
      
      lm.fit.ry              <- tsls(y=ypool[[k]],d=dpool[[k]], x=NULL, z=zpool[[k]], intercept = FALSE)
      ate                    <- lm.fit.ry$coef[1];
      HCV.coefs              <- lm.fit.ry$vcov[1]
      
      STE_pool[1,k]          <- HCV.coefs
      TE_pool[1,k]           <- ate
      
    }
  }

  if(length(methods)==1){
    
    TE_pool[1,(length(methods)+1)]  <- TE_pool[1,(length(methods))]
    STE_pool[1,(length(methods)+1)] <- STE_pool[1,(length(methods))]
    
  }
  colnames(result)   <- c(methods, "best") 
  colnames(result2)  <- c(methods, "best") 
  rownames(MSE1)     <- c(methods, "best") 
  rownames(MSE2)     <- c(methods, "best") 
  rownames(MSE3)     <- c(methods, "best") 
  rownames(result)   <- c("ATE", "se")
  rownames(result2)  <- c("ATE", "se")

  if(DML=="DML1"){
    result[1,]         <- colMeans(TE)
    result[2,]         <- sqrt((STE))
  }
  
  if(DML=="DML2"){
    result[1,]         <- colMeans(TE_pool)
    result[2,]         <- sqrt((STE_pool))
  }
  
  if(est=="plinear"){   
    table <- rbind(result, rowMeans(MSE1), rowMeans(MSE2)) 
    rownames(table)[3:4]   <- c("MSE[Y|X]", "MSE[D|X]") 
  }
  
  if(est=="IV"){   
    table <- rbind(result, rowMeans(MSE1), rowMeans(MSE2) , rowMeans(MSE3)) 
    rownames(table)[3:5]   <- c("MSE[Y|X]", "MSE[D|X]", "MSE[Z|X]") 
  }

  if(est=="interactive"){    
    table <- rbind(result, rowMeans(MSE1), rowMeans(MSE2) , rowMeans(MSE3))   
    rownames(table)[3:5]   <- c("MSE[Y|X, D=0]", "MSE[Y|X, D=1]", "MSE[D|X]")
  }  
  
  if(est=="LATE"){    
    table <- rbind(result, rowMeans(MSE1), rowMeans(MSE2) , rowMeans(MSE3),rowMeans(MSE4) , rowMeans(MSE5))    
    rownames(table)[3:7]   <- c("MSE[Y|X, Z=0]", "MSE[Y|X, Z=1]", "MSE[D|X,Z=0]", "MSE[D|X,Z=1]" ,"MSE[Z|X]")
  }  
  
  colnames(table)[length(methods)+1] = "best"
  return(table)
}  


cond_comp <- function(datause, dataout, y, d, z=NULL, x, method, methods, plinear,xL, binary,flag=0, args, param, tune_param){
  
  ind_u <- if(plinear=="LATE") which(datause[,z]==1) else which(datause[,d]==1)
  ind_o <- if(plinear=="LATE")  which(dataout[,z]==1) else which(dataout[,d]==1)
  
  res   <- list(0)

  l        <- method
  
  if(tune_param[[l]]==0){ tune = NULL}
  else { tune=tune_param[[l]]}
  
  form1  <- as.formula(paste(y, "~", x));
  form2  <- as.formula(paste(d, "~", x));
  form3  <- as.formula(paste(z, "~", x));
  
  if(plinear=="LATE"){
    
    fitControl   <- trainControl(method = param$methodML[l], number = param$cv[l], repeats = param$rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=param$select[l])
    arg          <- c(list(form=form1, data = datause[ind_u,],  method = methods[l],  tuneGrid = tune, trControl = fitControl, preProcess=param$proces[l], tuneLength=param$tune[l]), args[[methods[l]]])
    fit.yz1      <- suppressWarnings(do.call(caret::train, arg))
    res$my_z1x   <- predict(fit.yz1, newdata=dataout, type="raw")
    res$err.yz1  <- min(fit.yz1$results[,'RMSE'])
    
    fitControl   <- trainControl(method = param$methodML[l], number = param$cv[l], repeats = param$rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=param$select[l])
    arg          <- c(list(form=form1, data = datause[-ind_u,],  method = methods[l],  tuneGrid = tune, trControl = fitControl, preProcess=param$proces[l], tuneLength=param$tune[l]), args[[methods[l]]])
    fit.yz0      <- suppressWarnings(do.call(caret::train, arg))
    res$my_z0x       <- predict(fit.yz1, newdata=dataout, type="raw")
    res$err.yz0      <- min(fit.yz1$results[,'RMSE'])
    
    fitControl   <- trainControl(method = param$methodML[l], number = param$cv[l], repeats = param$rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=param$select[l])
    arg          <- c(list(form=form2, data = datause[ind_u,],  method = methods[l],  tuneGrid = tune, trControl = fitControl, preProcess=param$proces[l], tuneLength=param$tune[l]), args[[methods[l]]])
    fit.dz1      <- suppressWarnings(do.call(caret::train, arg))
    res$md_z1x   <- predict(fit.dz1, newdata=dataout, type="raw")
    res$err.dz1  <- min(fit.dz1$results[,'RMSE'])

    
    if(flag==0){    
      fitControl   <- trainControl(method = param$methodML[l], number = param$cv[l], repeats = param$rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=param$select[l])
      arg          <- c(list(form=form2, data = datause[-ind_u,],  method = methods[l],  tuneGrid = tune_param[[l]], trControl = fitControl, preProcess=param$proces[l], tuneLength=param$tune[l]), args[[methods[l]]])
      fit.dz0      <- suppressWarnings(do.call(caret::train, arg))
      res$md_z0x       <- predict(fit.dz0, newdata=dataout, type="raw")
      res$err.dz0      <- min(fit.dz0$results[,'RMSE'])
    } 
  }
  
  
  if(plinear=="interactive"){

    fitControl   <- trainControl(method = param$methodML[l], number = param$cv[l], repeats = param$rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=param$select[l])
    arg          <- c(list(form=form1, data = datause[ind_u,],  method = methods[l],  tuneGrid = tune, trControl = fitControl, preProcess=param$proces[l], tuneLength=param$tune[l]), args[[methods[l]]])
    fit.yd1      <- suppressWarnings(do.call(caret::train, arg))
    res$my_d1x   <- predict(fit.yd1, newdata=dataout, type="raw")
    res$err.yd1  <- min(fit.yd1$results[,'RMSE'])

    fitControl   <- trainControl(method = param$methodML[l], number = param$cv[l], repeats = param$rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=param$select[l])
    arg          <- c(list(form=form1, data = datause[-ind_u,],  method = methods[l],  tuneGrid = tune, trControl = fitControl, preProcess=param$proces[l], tuneLength=param$tune[l]), args[[methods[l]]])
    fit.yd0      <- suppressWarnings(do.call(caret::train, arg))
    res$my_d0x   <- predict(fit.yd0, newdata=dataout, type="raw")
    res$err.yd0  <- min(fit.yd0$results[,'RMSE'])
    
  }

  if(plinear=="plinear" | plinear=="interactive" | plinear=="IV"){
    fitControl   <- trainControl(method = param$methodML[l], number = param$cv[l], repeats = param$rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=param$select[l])
    arg          <- c(list(form=form2, data = datause,  method = methods[l],  tuneGrid = tune, trControl = fitControl, preProcess=param$proces[l], tuneLength=param$tune[l]), args[[methods[l]]])
    fit.d        <- suppressWarnings(do.call(caret::train, arg))
    res$md_x     <- predict(fit.d, newdata=dataout, type="raw")
    res$rd       <- dataout[,d] -  res$md_x
    res$err.d    <- min(fit.d$results[,'RMSE'])
    
  }
    
  if(plinear=="plinear" | plinear=="IV"){  
    
    fitControl   <- trainControl(method = param$methodML[l], number = param$cv[l], repeats = param$rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=param$select[l])
    arg          <- c(list(form=form1, data = datause,  method = methods[l],  tuneGrid = tune, trControl = fitControl, preProcess=param$proces[l], tuneLength=param$tune[l]), args[[methods[l]]])
    fit.y        <- suppressWarnings(do.call(caret::train, arg))
    res$my_x     <- predict(fit.y, newdata=dataout, type="raw")
    res$ry       <- dataout[,y] -  res$my_x
    res$err.y    <- min(fit.y$results[,'RMSE'])
  
  }
  
  if(plinear=="IV" | plinear=="LATE"){
    
    fitControl   <- trainControl(method = param$methodML[l], number = param$cv[l], repeats = param$rep[l], allowParallel = FALSE, verboseIter=FALSE, search="random", selectionFunction=param$select[l])
    arg          <- c(list(form=form3, data = datause,  method = methods[l],  tuneGrid = tune, trControl = fitControl, preProcess=param$proces[l], tuneLength=param$tune[l]), args[[methods[l]]])
    fit.z        <- suppressWarnings(do.call(caret::train, arg))
    res$mz_x     <- predict(fit.z, newdata=dataout, type="raw")
    res$rz       <- dataout[,z] -  res$mz_x
    res$err.z    <- min(fit.z$results[,'RMSE'])
  }

  return(res);
}  








