#
#***************** 0 Objet methodes generiques ********************
#
# A) Selecteurs 
# B) Affectants
# C) Methode 

# if (!isGeneric("summary")){
#  setGeneric("summary", function(object, ...) standardGeneric("summary"))
# }

##### A) Selecteurs ############################################

# Carto3D MRIaggr
setGeneric(name="selectContrast",
           def=function(object,...){standardGeneric("selectContrast")}
)

# MRIaggr
setGeneric(name="selectClinic",
           def=function(object,...){standardGeneric("selectClinic")}
)

# MRIaggr
setGeneric(name="selectCoords",
           def=function(object,...){standardGeneric("selectCoords")}
)

# Carto3D MRIaggr  
setGeneric(name="selectDefault_value",
           def=function(object,...){standardGeneric("selectDefault_value")}
)

# MRIaggr 
setGeneric(name="selectDescStats",
           def=function(object,...){standardGeneric("selectDescStats")}
)

# MRIaggr 
setGeneric(name="selectHistory",
           def=function(object,...){standardGeneric("selectHistory")}
)

# MRIaggr 
setGeneric(name="selectHemispheres",
           def=function(object,...){standardGeneric("selectHemispheres")}
)

# Carto3D MRIaggr 
setGeneric(name="selectIdentifier",
           def=function(object,...){standardGeneric("selectIdentifier")}
)

# MRIaggr 
setGeneric(name="selectMidplane",
           def=function(object,...){standardGeneric("selectMidplane")}
)

# MRIaggr
setGeneric(name="selectN",
           def=function(object,...){standardGeneric("selectN")}
)

# MRIaggr
setGeneric(name="selectNormalization",
           def=function(object,...){standardGeneric("selectNormalization")}
)

# Carto3D MRIaggr 
setGeneric(name="selectParameter",
           def=function(object,...){standardGeneric("selectParameter")}
)

# MRIaggr 
setGeneric(name="selectTable",
           def=function(object,...){standardGeneric("selectTable")}
)

# Carto3D MRIaggr 
setGeneric(name="selectVoxelDim",
           def=function(object,...){standardGeneric("selectVoxelDim")}
)

# Carto3D MRIaggr 
setGeneric(name="selectVoxelSize",
           def=function(object,...){standardGeneric("selectVoxelSize")}
)

##### B) Affectants ############################################


# MRIaggr 
setGeneric(name="affectContrast<-",
           def=function(object,param=NULL,default_value=NULL,overwrite=FALSE,trace=TRUE,value){standardGeneric("affectContrast<-")}
)

# MRIaggr
setGeneric(name="affectClinic<-",
           def=function(object,add=FALSE,overwrite=FALSE,trace=TRUE,value){standardGeneric("affectClinic<-")}
)

# MRIaggr 
setGeneric(name="affectDescStats<-",
           def=function(object,name,overwrite=FALSE,trace=TRUE,value){standardGeneric("affectDescStats<-")}
)

# MRIaggr
setGeneric(name="affectHemisphere<-",
           def=function(object,overwrite=FALSE,trace=TRUE,value){standardGeneric("affectHemisphere<-")}
)

# MRIaggr
setGeneric(name="affectNormalization<-",
           def=function(object,overwrite=FALSE,trace=TRUE,value){standardGeneric("affectNormalization<-")}
)

# MRIaggr 
setGeneric(name="affectTable<-",
           def=function(object,type,overwrite=FALSE,trace=TRUE,value){standardGeneric("affectTable<-")}
)

# MRIaggr
setGeneric(name="supprContrast<-",
           def=function(object,trace=TRUE,value){standardGeneric("supprContrast<-")}
)

# MRIaggr
setGeneric(name="supprDescStats<-",
           def=function(object,trace=TRUE,value){standardGeneric("supprDescStats<-")}
)

##### C) Methodes ############################################

#### calc. ####

# MRIaggr
setGeneric(name ="calcBrainMask",
           def=function(object,...){
             res <-standardGeneric("calcBrainMask")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               affectContrast(object,param="mask",overwrite=res$overwrite,trace=res$trace) <- res$res$best_group
               
               # update history
               object@history <- c(object@history,
                                   list(calcBrainMask=list(call=match.call(),date=date(),mask_name=res$res$mask_name))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                envir=.GlobalEnv)",
                 sep="")))             
             }
             
             return(invisible(res$res))
             
           }            
)

# MRIaggr
setGeneric(name="calcControlateral",
           def=function(object,...){
             
             res <- standardGeneric("calcControlateral")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               newdata <- data.frame(matrix(NA,nrow=selectN(object),ncol=ncol(res$data)))
               names(newdata) <- names(res$data)   
               
               newdata[res$data[,"index"],] <- res$data
               newdata <- newdata[,names(newdata) %in% c("index","i_hemisphere","j_hemisphere","hemisphere") == FALSE]
               
               affectContrast(object,overwrite=res$overwrite,trace=res$trace) <- newdata
               
               # update history
               object@history <- c(object@history,
                                   list(calcControlateral=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                envir=.GlobalEnv)",
                 sep="")))         
             }
             
             res$update.object <- NULL
             res$overwrite <- NULL
             res$trace <- NULL
             
             return(invisible(res))
             
           }
)

# MRIaggr
setGeneric(name ="calcDistMask",
           def=function(object,...){
             res <- standardGeneric("calcDistMask")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               affectContrast(object,overwrite=res$overwrite,trace=res$trace) <- res$res
               
               # update history
               object@history <- c(object@history,
                                   list(calcDistMask=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                envir=.GlobalEnv)",
                 sep="")))             
             }
             
             return(invisible(res$res))
             
           }            
)


# MRIaggr
setGeneric(name="calcDistTissues",
           def=function(object,...){standardGeneric("calcDistTissues")}
)

# MRIaggr
setGeneric(name="calcFilter",
           def=function(object,...){
             res <- standardGeneric("calcFilter")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               param <- setdiff(names(res$res),c("i","j","k"))              
               
               affectContrast(object,overwrite=res$overwrite,trace=res$trace) <- res$res[,param,drop=F]
               
               # update history
               object@history <- c(object@history,
                                   list(calcFilter=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                    envir=.GlobalEnv)",
                 sep="")))                 
             }
             
             res$update.object <- NULL
             res$overwrite <- NULL
             res$trace <- NULL
             
             return(invisible(res))
             
           }
)

# MRIaggr
setGeneric(name="calcGroupsMask",
           def=function(object,...)
           {  res <- standardGeneric("calcGroupsMask")
              
              if(res$update.object==TRUE){
                
                # affect
                nom_object <- as.character(substitute(object))
                affectDescStats(object,name="GroupsLesion",overwrite=res$overwrite,trace=res$trace) <- lapply(res$res,function(x){x$group_size})
                
                # update history
                object@history <- c(object@history,
                                    list(calcGroupsMask=list(call=match.call(),date=date()))
                )
                
                # update object
                eval(parse(text=paste(
                  "assign(\"",nom_object,"\",value=object,
                envir=.GlobalEnv)",
                  sep="")))  
              }
              
              return(invisible(res$res))
              
           }
)

# MRIaggr
setGeneric(name="calcHemisphere",
           def=function(object,...){
             res <- standardGeneric("calcHemisphere")
             
             if(res$update.object==TRUE){
            
               # affect
               nom_object <- as.character(substitute(object))
               affectHemisphere(object,overwrite=res$overwrite,trace=res$trace) <- list(midplane=res$res$midplane,                                                                                            
                                                                                        data=res$res$data)
               
               if(!is.null(res$res$hemispheres)){
                 affectHemisphere(object,overwrite=res$overwrite,trace=res$trace) <- list(hemispheres=res$res$hemispheres)
               }
               
               # update history
               object@history <- c(object@history,
                                   list(calcHemisphere=list(call=match.call(),date=date(),optimum=res$res$optimum[,c("position_i","position_j","angle_rad")]))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                envir=.GlobalEnv)",
                 sep="")))  
             }  
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name="calcROCthreshold",
           def=function(object,...){ 
             res <- standardGeneric("calcROCthreshold")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               affectDescStats(object,name="Mask_threshold",overwrite=res$overwrite,trace=res$trace) <- res$res
               
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                 envir=.GlobalEnv)",
                 sep="")))
               
               # update history
               object@history <- c(object@history,
                                   list(calcROCthreshold=list(call=match.call(),date=date()))
               )
               
               # update object
               res$update.object <- NULL
               res$overwrite <- NULL
             }
             return(invisible(res$res))
           }
)

# MRIaggr
setGeneric(name="calcNormalization",
           def=function(object,...){
             
             res <- standardGeneric("calcNormalization")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               affectNormalization(object,overwrite=res$overwrite,trace=res$trace) <- res$res
               
               # update history
               object@history <- c(object@history,
                                   list(calcNormalization=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                 envir=.GlobalEnv)",
                 sep="")))             
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name="calcRegionalContrast",
           def=function(object,...){
             res <- standardGeneric("calcRegionalContrast")
             
             if(!is.list(res)){
               return(res)
             }else{
               if(res$update.object==TRUE){
                 
                 # affect
                 nom_object <- as.character(substitute(object))
                 affectContrast(object,overwrite=res$overwrite,trace=res$trace) <- res$res
                 
                 # update history
                 object@history <- c(object@history,
                                     list(calcRegionalContrast=list(call=match.call(),date=date()))
                 )
                 
                 # update object
                 eval(parse(text=paste(
                   "assign(\"",nom_object,"\",value=object,
                   envir=.GlobalEnv)",
                   sep="")))
               }
               
               return(invisible(res$res))
             }
             
           }
)

# MRIaggr
setGeneric(name="calcSmoothMask",
           def=function(object,...){
             res <- standardGeneric("calcSmoothMask")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               affectContrast(object,param="mask",overwrite=res$overwrite,trace=res$trace) <- res$res$mask
               
               # update history
               object@history <- c(object@history,
                                   list(calcSmoothMask=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                 envir=.GlobalEnv)",
                 sep="")))   
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name="calcTableHypoReperf",
           def=function(object,...){          
             res <- standardGeneric("calcTableHypoReperf")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               
               if("volume_hypo" %in% names(res$res)){
                 affectTable(object,type="hypoperfusion",overwrite=res$overwrite,trace=res$trace) <- res$res$volume_hypo
               }
               
               if("volume_reperf" %in% names(res$res)){
                 affectTable(object,type="reperfusion",overwrite=res$overwrite,trace=res$trace) <- res$res$volume_reperf
               }
               
               if("pixel" %in% names(res$res)){
                 nom_param <- names(res$res$pixel)
                 param.reperf_pc <- grep(pattern="reperf_pc",nom_param,value=TRUE)
                 param.reperf <- grep(pattern="reperf",nom_param[nom_param %in% param.reperf_pc == FALSE],value=TRUE)
                 param.deperf_pc <- grep(pattern="deperf_pc",nom_param,value=TRUE)
                 param.deperf <- grep(pattern="deperf",nom_param[nom_param %in% param.deperf_pc == FALSE],value=TRUE)
                 param.shift <- grep(pattern="shift",nom_param,value=TRUE)
                 
                 eval(parse(text=paste(
                   "nom_param <- c(\"i\",\"j\",\"k\",",paste(paste("param.",res$param.update,sep=""),collapse=","),")",
                   sep="")))
                 
                 affectContrast(object,overwrite=res$overwrite,trace=res$trace) <- res$res$pixel[,nom_param]
               }
               
               # update history
               object@history <- c(object@history,
                                   list(calcTableHypoReperf=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                 envir=.GlobalEnv)",
                 sep="")))            
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name="calcTableLesion",
           def=function(object,...){ 
             res <- standardGeneric("calcTableLesion")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               affectTable(object,type="lesion",overwrite=res$overwrite,trace=res$trace) <- res$res
               
               # update history
               object@history <- c(object@history,
                                   list(calcTableLesion=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                 envir=.GlobalEnv)",
                 sep="")))             
             }
             
             return(invisible(res$res))
             
           }
)

setGeneric(name ="calcThresholdMRIaggr",
           def=function(object,...){
             res <- standardGeneric("calcThresholdMRIaggr")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               default_value <- data.frame(matrix(TRUE,ncol=length(res$name_newparam)))
               names(default_value) <- res$name_newparam
             
               affectContrast(object,param=res$name_newparam,default_value=default_value,overwrite=res$overwrite,trace=res$trace) <- res$res[,res$name_newparam]
               
               # update history
               object@history <- c(object@history,
                                   list(calcThresholdMRIaggr=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                 envir=.GlobalEnv)",
                 sep="")))             
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name="calcTissueType",
           def=function(object,...){
             res <- standardGeneric("calcTissueType")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               
               affectContrast(object,param=res$name_newparam,overwrite=res$overwrite,trace=res$trace) <- res$res$prob
               
               # update history
               object@history <- c(object@history,
                                   list(calcTissueType=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                 envir=.GlobalEnv)",
                 sep="")))
             }
                          
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name="calcW",
           def=function(object,...){
             res <- standardGeneric("calcW")
             
             if(res$update.object==TRUE){
               
               # affect
               nom_object <- as.character(substitute(object))
               affectDescStats(object,name="W_euclidean",overwrite=res$overwrite,trace=res$trace) <- res$res
               
               # update history
               object@history <- c(object@history,
                                   list(calcW=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                envir=.GlobalEnv)",
                 sep="")))      
               
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name ="outlineMRIaggr",
           def=function(object,...){
             res <- standardGeneric("outlineMRIaggr")
             
             if(res$update.object==TRUE){
            
               # affect
               nom_object <- as.character(substitute(object))
               default_value <- data.frame(TRUE)
               names(default_value) <- res$name_newparam
               
               affectContrast(object,param=res$name_newparam,default_value=default_value,overwrite=res$overwrite,trace=res$trace) <- res$res[,c("i","j","k",res$name_newparam)]
               
               # update history
               object@history <- c(object@history,
                                   list(outlineMRIaggr=list(call=match.call(),date=date()))
               )
               
               # update object
               eval(parse(text=paste(
                 "assign(\"",nom_object,"\",value=object,
                envir=.GlobalEnv)",
                 sep="")))             
             }
             
             return(invisible(res$res))
             
           }
)

#### plot ####

# MRIaggr
setGeneric(name="boxplotMask",
           def=function(object,...){
             standardGeneric("boxplotMask")
           }
)

# MRIaggr
setGeneric(name="heatmapMRIaggr",
           def=function(object,...){
             standardGeneric("heatmapMRIaggr")
           }
)

# Carto3D MRIaggr  
setGeneric(name="multiplot",
           def=function(object,...){
             standardGeneric("multiplot")
           }
)

# MRIaggr
setGeneric(name="pointsHemisphere",
           def=function(object,...){
             standardGeneric("pointsHemisphere")
           }
)

# MRIaggr
setGeneric(name="plotLesion3D",
           def=function(object,...){
             standardGeneric("plotLesion3D")
           }
)

# MRIaggr
setGeneric(name="plotTableLesion",
           def=function(object,...){
             standardGeneric("plotTableLesion")
           }
)

# MRIaggr
setGeneric(name="plotDistClass",
           def=function(object,...){
             standardGeneric("plotDistClass")
           }
)

#### const. ####

# MRIaggr 
setGeneric(name="constCompressMRIaggr",
           def=function(object,...){
             standardGeneric("constCompressMRIaggr")
           }
)

# MRIaggr
setGeneric(name="constReduceMRIaggr",
           def=function(object,...){
             standardGeneric("constReduceMRIaggr")
           }
)

#### init.  ####
# Carto3D MRIaggr
setGeneric(name="initNum",
           def=function(object,...){
             standardGeneric("initNum")
           }
)

# MRIaggr
setGeneric(name="initParameter",
           def=function(object,...){
             standardGeneric("initParameter")
           }
)
