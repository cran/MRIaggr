#### 1- calc functions ####

# calcCriteriaGR <- function(contrast,groups,W=NULL,sigma=NULL,breaks,
#                            criterion_d1=FALSE,criterion_entropy=TRUE,criterion_Kalinsky=TRUE,criterion_Laboure=TRUE){
#   
#   #### preparation ####
#   n <- length(contrast) 
#   
#   if(criterion_d1==FALSE && criterion_entropy==FALSE && criterion_Kalinsky == FALSE && criterion_Laboure == FALSE){
#     stop("calcCriteriaGR : all criteria must not be simultaneously FALSE \n")
#   }
#   
#   criteria <- c("d1","entropy","Kalinsky","Laboure")[c(criterion_d1,criterion_entropy,criterion_Kalinsky,criterion_Laboure)]
#   
#   df.criterion <- data.frame(matrix(NA,ncol=1+length(criteria),nrow=1))
#   names(df.criterion) <- c("sigma",criteria)
#   if(!is.null(sigma)){df.criterion$sigma <- sigma}
#   
#   #### tests ####
#   if(length(groups)!=n){
#     stop("calcCriteriaGR : \'groups\' must have same length as \'contrast\' \n",       
#          "length(contrast) : ",length(contrast),"\n",
#          "length(groups) : ",length(groups),"\n",)
#   }
#   
#   if(!is.logical(groups)){
#     stop("calcCriteriaGR :  wrong specification of \'groups\' \n",
#          "groups must be a logical vector \n",
#          "is(groups) : ",paste(is(groups),collapse=" "),"\n")
#   }
#   indexT <- which(groups==TRUE)
#   n.T <- length(indexT)
#   indexF <- which(groups==FALSE)
#   n.F <- length(indexF)
#   
#   if(length(indexT)<2 || length(indexF)<2){
#     stop("calcCriteriaGR :  too small groups \n",
#          "sum(groups==TRUE) : ",sum(groups==TRUE),"\n",
#          "sum(groups==FALSE) : ",sum(groups==FALSE),"\n")
#     
#     return(NULL)
#   }
#   
#   #### criterion sur la frontiere ####
#   if(criterion_d1==TRUE){
#     
#     if(is.null(W)){
#       stop("calcCriteriaGR : \'W\' must not be NULL \n")
#     }
#     
#     if(any(dim(W)!=c(n,n))){
#       stop("calcCriteriaGR : wrong dimension of \'W\' \n",
#            "required dimension : ",n," ; ",n,"\n",
#            "dim(W) :",paste(dim(W),collapse=" "),"\n")
#     }
#     
#     
#     # identification de la frontiere
#     nb.voisT <- Matrix::rowSums(W[indexT,indexF]>0)
#     index.frontiereT <- indexT[which(nb.voisT>0)]
#     n.frontiere <- length(index.frontiereT)
#     
#     # iteration sur les points de  la frontiere
#     df.criterion$d1 <- 0 
#     
#     for(iter_frontiere in 1:n.frontiere){
#       indexiter_frontiere <- index.frontiereT[iter_frontiere]
#       indexiter_voisF <- indexF[which(W[indexiter_frontiere,indexF]>0) ]
#       indexiter_voisT <- indexT[which(W[indexiter_frontiere,indexT]>0) ]
#       
#       gradF <- contrast[indexiter_frontiere]-mean(contrast[indexiter_voisF])
#       gradT <- contrast[indexiter_frontiere]-mean(contrast[indexiter_voisT])
#       df.criterion$d1 <- df.criterion$d1 + (abs(gradF) - abs(gradT))    
#     }
#     df.criterion$d1 <- df.criterion$d1/n.frontiere
#     
#   }
#   
#   #### criterion sur l'entropy ####
#   if(criterion_entropy==TRUE){
#     RF <- hist(contrast[indexF],breaks=breaks,plot=FALSE)$counts
#     RT <- hist(contrast[indexT],breaks=breaks,plot=FALSE)$counts
#     
#     RF <- RF[RF>0]
#     RT <- RT[RT>0]
#     
#     df.criterion$entropy <- -mean(RF*log(RF)) - mean(RT*log(RT)) + log(length(indexF)*length(indexT))
#     #df.criterion$entropy <- sum(contrast_entropy[indexT]*log(contrast_entropy[indexT])) + sum( (1-contrast_entropy)[indexF]*log( (1-contrast_entropy)[indexF]))
#   }
#   
#   #### criterion sur la region ####
#   if(criterion_Kalinsky==TRUE || criterion_Laboure==TRUE){
#     moyT <- mean(contrast[indexT])
#     moyF <- mean(contrast[indexF])
#     WT <- sum((contrast[indexT]-moyT)^2)
#     WF <- sum((contrast[indexF]-moyF)^2)
#   }
#   
#   if(criterion_Kalinsky==TRUE){
#     moy <- mean(contrast)
#     
#     B <- (moy-moyT)^2+(moy-moyF)^2 # ddf(B) = 2-1 = 1
#     W <- (WT + WF)/(n.T+n.F-1) # ddf(W) = n-1
#     df.criterion$Kalinsky <- B / W
#   }
#   
#   if(criterion_Laboure==TRUE){
#     df.criterion$Laboure <- 1/sqrt(WT + WF)
#   }
#   
#   return(df.criterion)
# }



calcGR <- function(contrast,W,seed,sigma_max,range=c(-Inf,+Inf),range.seed=c(-Inf,+Inf),                 
                   breaks=100,scale=FALSE,iter_max=100,sd.robust=FALSE,trace=TRUE,
                   history_sigma=FALSE,history_step=FALSE,history_front=FALSE){
  
  res <- initGR(contrast=contrast,W=W,seed=seed,range=range,range.seed=range.seed,
                breaks=breaks,scale=scale,trace=trace,method="calcGR")
  
  contrast <- res$contrast
  breaks <- res$breaks 
  step <- res$step 
  seed <- res$seed 
  
  #### algorithm ####
  if(trace==TRUE){cat("GR (sigma_max=",sigma_max,") : ",sep="")}
   
  resGR <- GRalgo(contrast=contrast,W=W,seed=seed,sigma_max=sigma_max,range=range,
                  breaks=breaks,step=step,operator=if(sd.robust){"mad"}else{"sd"},iter_max=iter_max,                    
                  history_sigma=history_sigma,history_step=history_step,history_front=history_front)
    
  if(trace==TRUE){cat(length(resGR$GR)," observations in ",resGR$iter," iterations \n",sep="")}
  
  if(resGR$iter==0 || (resGR$iter==iter_max && resGR$test.id==FALSE) || resGR$test.break==TRUE){
    cv <- FALSE
    if(trace==TRUE){
      cat("WARNING :  the algorithme has not converged \n")   
      if(resGR$iter==0 || (resGR$iter==iter_max &&  resGR$test.id==FALSE)){
        cat("maximum number of iterations reached (iter_max=",iter_max,")\n",sep="")
      }else{
        cat("algorithm breaks due to insufficient number of observations in GR (\'sigma\' may be insufficient)\n")
      }
    }
    
  }
  
  #### export ####
  resGR$breaks <- breaks

  if(history_front && resGR$iter>1){
    front.split <- lapply(strsplit(resGR$history_front,split=".",fixed=TRUE),
                          function(x){c(as.numeric(x),rep(NA,resGR$iter-length(x)))}
    )
    resGR$Mfront <- matrix(unlist(front.split),ncol=resGR$iter,nrow=length(resGR$GR),byrow=TRUE)
  }else{
    resGR$Mfront <- NULL
  }
  
  resGR$breaks <- breaks
  return(resGR)
  
}

# calcSigmaGR <- function(contrast,W,seed,sigma,range=c(-Inf,+Inf),range.seed=c(-Inf,+Inf),
#                         breaks=100,scale=FALSE,iter_max=100,sd.robust=FALSE,
#                         criterion_d1=FALSE,criterion_entropy=TRUE,criterion_Kalinsky=TRUE,criterion_Laboure=TRUE,
#                         trace=TRUE,mar=rep(2,4),mgp=c(2,0.75,0),main="",
#                         window=FALSE,filename="calcSigmaGR",width=1000,height=700,path=NULL,unit="px",res=NA){
#   
#   
#   #### initialization ####
#   sigma <- sort(sigma)
#   
#   ## GR
#   n <- length(contrast)
#   
#   res.initGR <- initGR(contrast=contrast,W=W,seed=seed,range=range,range.seed=range.seed,
#                        breaks=breaks,scale=scale,trace=trace,method="calcSigmaGR")
#   
#   contrast <- res.initGR$contrast
#   breaks <- res.initGR$breaks 
#   step <- res.initGR$step 
#   seed.sigma <- res.initGR$seed 
#   
#   # window
#   res.initW <- initWindow(window,filename,path,width,height,unit,res,
#                           n.plot=1,mfrow=c(1,1),xlim=NULL,ylim=NULL,method="calcSigmaGR")
#   scale <- res.initW$scale
#   
#   ## specific criteria
#   n.sigma <- length(sigma)
#   
#   criteria <- c("d1","entropy","Kalinsky","Laboure")[c(criterion_d1,criterion_entropy,criterion_Kalinsky,criterion_Laboure)]
#   n.criteria <- length(criteria)
#   if(n.criteria==0){
#     stop("calcSigmaGR : all criteria must not be simultaneously FALSE \n")
#   }
#   
#   df.criterion <- data.frame(matrix(NA,ncol=1+n.criteria,nrow=n.sigma))
#   names(df.criterion) <- c("sigma",criteria)
#   df.criterion$sigma <- sigma
#   
#   
#   #### loop ####
#   list.GR <- sapply(1:n.criteria,function(x){c(list(),NA)})
#   names(list.GR) <- criteria
#   
#   best <- data.frame(matrix(-Inf,nrow=1,ncol=4))
#   names(best) <- criteria
#   
#   if(trace==TRUE){cat("loop over ",n.sigma," sigma : ",sep="")}
#   
#   for(iter_sigma in 1:n.sigma){
#     if(trace==TRUE){cat("*")}
#     
#     ## GR
#     resGR <- GRalgo(contrast=contrast,W=W,seed=seed,sigma_max=sigma[iter_sigma],range=range,
#                     breaks=breaks,step=step,operator=if(sd.robust){"mad"}else{"sd"},iter_max=iter_max,                    
#                     history_sigma=TRUE,history_step=FALSE,history_front=FALSE)
#     
#     ## Compute criteria
#     df.criterion[iter_sigma,] <- calcCriteriaGR(contrast=contrast,groups=(1:n) %in% resGR$GR,W=W,
#                                                 sigma=sigma[iter_sigma],breaks=breaks,
#                                                 criterion_d1=criterion_d1,criterion_entropy=criterion_entropy,
#                                                 criterion_Kalinsky=criterion_Kalinsky,criterion_Laboure=criterion_Laboure)
#     
#     ## save optimum
#     for(iter_criteria in 1:n.criteria){
#       if(!is.na(df.criterion[iter_sigma,iter_criteria+1]) && df.criterion[iter_sigma,iter_criteria+1]>best[iter_criteria]){
#         list.GR[[iter_criteria]] <- resGR$GR
#         best[iter_criteria] <- df.criterion[iter_sigma,iter_criteria+1]
#       }            
#     }
#     
#     ## update
#     if(length(resGR$GR)>0){seed.sigma <- resGR$GR}    
#     if(length(resGR$GR)==n){break}
#   }
#   if(trace==TRUE){cat("\n")}
#   
#   
#   if(!is.null(plot)){
#     
#     df.criterion_scale <- scale(df.criterion[,-1,drop=F])
#     diff_sigma <- min(diff(sigma))
#     diff_sigmaOpt <- seq(-diff_sigma/10,+diff_sigma/10,length.out=n.criteria)
#     
#     initDisplayWindow(window=window,filename=filename,path=path,width=width,height=height,scale=scale,res=res,
#                       mfrow=c(1,1),bg=par()$bg,pty=par()$pty,mar=mar,mgp=mgp)
#     
#     plot(x=sigma,y=rep(NA,n.sigma),ylim=c(0,0.5)+range(df.criterion_scale,na.rm=T),xlab="sigma",ylab="scaled value",main=main)
#     for(iter_criterion in 1:n.criteria){
#       points(sigma,df.criterion_scale[,iter_criterion],col=rainbow(10)[iter_criterion],
#              type="b",pch=20)
#       
#       diff_opt <- abs(df.criterion[,iter_criterion+1]-as.numeric(best[iter_criterion]))
#       abline(v=sigma[which.min(diff_opt)]+diff_sigmaOpt[iter_criterion],col=rainbow(10)[iter_criterion])
#     }    
#     
#     legend("topright",legend=criteria,col=rainbow(10)[1:n.criteria],pch=20,bty="n")
#     
#     switch(window,
#            eps=dev.off(),
#            svg=dev.off(),
#            png=dev.off()
#     )
#   }
#   
#   return(list(df.criterion=df.criterion,
#               list.GR=list.GR,
#               best=best)
#   )
# }

GRalgo <- function(contrast,W,seed,sigma_max,range,breaks,step,operator,iter_max,
                   history_sigma,history_step,history_front){  
  
  #### initialisation ####
  if(history_sigma){
    hist_sigma <- matrix(NA,nrow=iter_max,ncol=2)
  }
  if(history_step){
    hist_stepGR <- rep(0,length(seed))
  }
  if(history_front){
    hist_frontGR <- rep(0,length(seed))
  }
  
  iter <- 0
  test.break <- FALSE  
  
  GR <- seed
  GRsave <- NULL
  
  breaks.GR <- breaks[hist(contrast[GR],breaks=breaks,plot=FALSE)$density>0]
  n.breaks.GR <- length(breaks.GR)
  
  ##### iterations
  while(identical(GR,GRsave)==FALSE && iter < iter_max){
    iter <- iter+1   
    GRsave <- GR  

    # dilatation
    Dnew <- unique(c(GR,which(spam::colSums(W[GR,,drop=FALSE])>0)))
    
    Dnew <- Dnew[ (contrast[Dnew]>=range[1])*(contrast[Dnew]<=range[2]) == 1 ] # supplement perso pour ne pas ajouter des px en dehors du seuil possible

    
    # sigma
    if(length(Dnew)<=1){GR <- Dnew ; test.break = TRUE ; break} # check sigma validity
    sigma <- eval(parse(text=paste(operator,"(contrast[Dnew])",sep="")))
    
    if(history_sigma==TRUE){hist_sigma[iter,1] <- sigma}
    
    if(sigma <= sigma_max){ # homogeneous      
      
      GR <- Dnew      
      
    }else{ # reduce dilatation (les breaks correspondent a la valeur inf)
      
      Dnew_red <- Dnew[ (contrast[Dnew] >= breaks.GR[1])*(contrast[Dnew] <= (breaks.GR[length(breaks.GR)]+step))==1 ] 
      
      # sigma
      if(length(Dnew_red)<=1){ test.break = TRUE ; break}
      sigma <- eval(parse(text=paste(operator,"(contrast[Dnew_red])",sep="")))
      
      if(sigma <= sigma_max){ # homogenous : extensions
        
        Neighborhood <- setdiff(which(spam::colSums(W[Dnew_red,,drop=FALSE])>0),
                                Dnew_red)
        
        test.leftExtension <- (contrast[Neighborhood] >= (breaks.GR[1]-step))*(contrast[Neighborhood] < (breaks.GR[n.breaks.GR]+step))
        leftExtension <- Neighborhood[which(test.leftExtension==1)]
        test.rightExtension  <- (contrast[Neighborhood] >= (breaks.GR[1]))*(contrast[Neighborhood] < (breaks.GR[n.breaks.GR]+2*step))
        rightExtension <- Neighborhood[which(test.rightExtension==1)]
        
        # sigma 
        sigma <- eval(parse(text=paste(operator,"(contrast[c(Dnew_red,leftExtension,rightExtension)])",sep="")))
        
        if(sigma <= sigma_max){ # homogeneous full extension
          
          GR <- c(Dnew_red,unique(c(leftExtension,rightExtension)))
          
        }else{ # homogeneous left, right or no extension
          
          if(length(leftExtension)>0 && breaks.GR[1]-step>=range[1] ){
            sigmaL <- eval(parse(text=paste(operator,"(contrast[c(Dnew_red,leftExtension)])",sep="")))
          }else{sigmaL <- sigma_max+1}
          
          if(length(rightExtension)>0 && breaks.GR[n.breaks.GR]+2*step<=range[2] ){
            sigmaR <- eval(parse(text=paste(operator,"(contrast[c(Dnew_red,rightExtension)])",sep="")))
          }else{sigmaR <- sigma_max+1}
          
          if(sigmaL <= sigma_max && sigmaL <= sigmaR){
            
            GR <- c(Dnew_red,leftExtension)
            
          }else if(sigmaR <= sigma_max){
            
            GR <- c(Dnew_red,rightExtension)     
            
          }else{
            
            GR <- Dnew_red
            
          }
          
        }
      }else{  # contraction
        ordreC <- 1      
        max <- 0
        
        while(max==0){
          if(ordreC>n.breaks.GR){GR <- GRsave ;  test.break = TRUE ; break}
          
          for(iter_C in 0:ordreC){
            
            breaks.min <- breaks.GR[1]+iter_C*step
            breaks.max <- breaks.GR[n.breaks.GR]-(ordreC-iter_C)*step
            
            C <- Dnew_red[ (contrast[Dnew_red] >= breaks.min)*(contrast[Dnew_red] <= breaks.max) == 1 ]
            if(length(C)<=1){GR <- GRsave ;  test.break = TRUE ; break }
            
            sigma.test <- eval(parse(text=paste(operator,"(contrast[C])",sep="")))
            
            if(sigma.test <= sigma_max && length(C)>max){
              max <- length(C)
              GR <- C  
            }
          }
          
          if(max==0){ordreC <- ordreC+1}
          
        }
        
      } # end non homogeneous dilatation      
    } # end non homogeneous
    
    GR <- sort(GR)
    breaks.GR <- breaks[hist(contrast[GR],breaks=breaks,plot=FALSE)$density>0]
    n.breaks.GR <- length(breaks.GR)
    
    if(history_sigma==TRUE){hist_sigma[iter,2] <- sd(contrast[GR])}
    
    if(history_step==TRUE || history_front==TRUE){
      GR_in_GRsave <- which(GR %in% GRsave)
      GRsave_in_GR <- which(GRsave %in% GR)
      
      if(history_step==TRUE){
        hist_stepGR.sauve <- hist_stepGR
        hist_stepGR <- rep(iter,length(GR))
        hist_stepGR[GR_in_GRsave] <- hist_stepGR.sauve[GRsave_in_GR] 
      }
      
      if(history_front==TRUE){
        hist_frontGR.sauve <- hist_frontGR
        hist_frontGR <- rep(-1,length(GR))
        
        # old GR points
        hist_frontGR[GR_in_GRsave] <- hist_frontGR.sauve[GRsave_in_GR]      
        
        # new GR points
        newGR <- setdiff(GR,GRsave)
        if(length(newGR)>0){
          front <- as.factor(calcGroupsW_cpp(W,newGR-1,10000)$group_subset)
          
          # order front by 
          levels(front) <- levels(front)[rank(-table(front)+seq(0,0.1,length.out=length(levels(front))))]
          front <- as.numeric(as.character(front))
          
          if(iter==1){
            hist_frontGR[GR %in% GRsave == 0] <- front
          }else{
            
            origin <- sapply(newGR,function(x){  # privilegie l appartenance aux grands fronts i.e. max(front)
              index_Wn0 <- which(W[x,GRsave]>0)
              valW <- W[x,GRsave[index_Wn0]]
              if(length(valW)>0){max(hist_frontGR.sauve[index_Wn0[which(valW==min(valW))]])}else{""}
            })
            
            newGR2 <- newGR[which(nchar(origin)==0)]
            
            if(length(newGR2)>0){ # cas ou il y a eu une extension il faut recuperer le voisinage precedent
              newGR_noExp <- which(newGR %in% c(leftExtension,rightExtension) == FALSE)
              restrictedGR <- newGR[newGR_noExp]           
              
              origin2 <- sapply(newGR2,function(x){
                index_Wn0 <- which(W[x,restrictedGR]>0)
                W_red <- W[x,restrictedGR[index_Wn0]]
                max(origin[newGR_noExp[index_Wn0[which(W_red == min(W_red))]]])
              })
              origin[newGR %in% newGR2] <- origin2
            }
            
            hist_frontGR[GR %in% newGR] <- paste(origin,front,sep=".")              
          }
        }
        hist_frontGR[GR_in_GRsave] <- paste(hist_frontGR[GR_in_GRsave],".",sep="")         
      } # end if history_front
    } # end if history_step || history_front
    
  } # end while loop 
  
  #### export ####
  return(list(GR=GR,
              test.break=test.break,
              iter=iter,
              test.id=identical(GR,GRsave),
              sigma=if(history_sigma && iter>0){hist_sigma[1:iter,]}else{NULL},
              history_step=if(history_step && iter>0){hist_stepGR}else{NULL},
              history_front=if(history_front && iter>0){hist_frontGR}else{NULL}
  ))
  
}

#### 2- init functions ####

initGR <- function(contrast,W,seed,range,range.seed,
                   breaks,scale,trace,method="initGR"){
 
  #### tests ####
  if(length(seed)==0){
    stop(method," : wrong specification of \'seed\' \n",
         "\'seed\' is empty \n")
  }
  
  if(any(is.na(contrast))){
    stop(method," : wrong specification of \'contrast\' \n",
         "\'contrast\' contains NA \n",
         "number of NA : ",sum(is.na(contrast)),"\n")
  }
  
  if(!is.numeric(contrast)){
    stop(method," : wrong specification of \'contrast\' \n",
         "\'contrast must be a numeric vector \n",
         "is(contrast) : ",paste(is(contrast),collapse=" "),"\n")
  }
  
  if(nrow(W) != length(contrast) || ncol(W) != length(contrast)){
    stop(method," : \'W\' dimension is incoherent with \'contrast\' length \n",
         "length(contrast) : ",length(contrast),"\n",
         "dim(W) : ",paste(dim(W),collapse=" "),"\n")
  }
  
  if(length(range)!=2 || !is.numeric(range) || sum(is.na(range))>0){
    stop(method," : wrong specification of \'range\' \n",
         "\'range\' must be a numeric vector of length 2 without NA \n",
         "proposed \'range\' : ",paste(range,collapse=" ")," \n")
  }
  
  if( range[1] > max(contrast) || range[2] < min(contrast) ){
    stop(method," : wrong specification of \'range\' \n",
         "requested contrast range : ",paste(range,collapse=" "),"\n",
         "range of \'contrast\' : ",paste(range(contrast),collapse=" "),"\n")
  }
  
  if(length(range.seed)!=2 || !is.numeric(range.seed) || sum(is.na(range.seed))>0){
    stop(method," : wrong specification of \'range.seed\' \n",
         "\'range.seed\' must be a numeric vector of length 2 without NA \n",
         "proposed \'range.seed\' : ",paste(range.seed,collapse=" ")," \n")
  }
  
  #### settings ####
  if(scale==TRUE){
    contrast <- as.numeric(scale(contrast))
  }
  
  if(length(breaks)==1){
    breaks <- seq(min(contrast),max(contrast),length.out=breaks)
  }
  step <- breaks[2]-breaks[1]+10^{-12}
  
  # seeds
  if(is.logical(seed)){seed <- which(seed==TRUE)}  
  
  if(any(seed %% 1 != 0) || any(seed<=0) || any(seed>length(contrast))){
    stop(method," : incorrect specification of \'seed\' \n",
         "seed must contains integers between 1 and ",length(contrast),"\n",
         "range(seed) : ",paste(range(seed),collapse=" ")," \n",
         "sum(seed %% 1 !=0) : ",sum(seed %% 1 !=0),"\n")
  }
  
  test.seed <- (contrast[seed] >= range.seed[1])*(contrast[seed] <= range.seed[2])
  if(sum(test.seed)==0){
    stop(method," : incorrect specification of \'range.seed\' \n",
         "requested contrast range : ",paste(range.seed,collapse=" "),"\n",
         "contrast range of the \'seeds\' ",paste(range(contrast[seed]),collapse=" "),"\n")
  }
  if(trace==TRUE){cat("number of valid seeds : ",sum(test.seed)," over ",length(seed)," seeds \n",sep="")}
  seed <- seed[test.seed==1]
  
  if(sum((contrast[seed] >= range[1])*(contrast[seed] <= range[2]))==0){
    stop(method," : incorrect specification of \'seed\' \n",
         "requested contrast range : ",paste(seed,collapse=" "),"\n",
         "contrast range of \'seeds\' ",paste(range(contrast[seed]),collapse=" "),"\n")
  }
  
  #### export ####  
  
  res <- list()
  res$contrast <- contrast
  res$breaks <- breaks
  res$step <- step
  res$seed <- seed
  
  return(res)
}