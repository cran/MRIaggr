#**********************************************************************
#**********************************************************************
#*************         0A Class MRIaggr             *******************
#**********************************************************************
#**********************************************************************
#
###### Sommaire #################
# A) Class definition
# B) Selecters 
# C) Allocators
# D) Methods 


##### A) Class definition #############################################################

setClass(
  
  Class="MRIaggr",
  
  representation(
    identifier = "character",   # id patient de la carto
    contrast = "data.frame",      # data.frame integrant toutes les cartos d interet de dim (L*nb,l)
    clinic = "data.frame",
    fieldDim = "data.frame", 
    voxelDim = "data.frame",
    default_value = "data.frame",
    
    history="list",
    normalization = "list",
    hemispheres = "data.frame",
    midplane="data.frame",
    
    table_lesion = "data.frame", 
    table_reperfusion = "data.frame",             # both elements : reperfusion par threshold et par intervalle
    table_hypoperfusion = "data.frame",
    
    ls_descStats = "list"
  ),
  
  validity = function(object){
    # cat("--- MRIaggr : checking --- ")
    
    if(any(c("i","j","k","index") %in% names(object@contrast)==FALSE))
    {stop("validity[MRIaggr] : names of \'@contrast\' are not correct\n",
          "it must contains : \"i\" \"j\" \"k\" \"index\" \n",
          "proposed names : ",paste(names(object@contrast),collapse=" "),"\n")
    }
    
    param_data <- names(object@contrast)[names(object@contrast) %in% c("i","j","k","index","mask")==FALSE]
    param_default_values <- names(object@default_value)[names(object@default_value) %in% c("i","j","k","index")==FALSE]
    
    if(any(param_data %in% param_default_values == FALSE) || any(param_default_values %in% param_data == FALSE))
    {stop("validity[MRIaggr] : names of \'@default_value\' do not match those of @contrast \n",
          "missing \'@default_value\'  : ",paste(param_data[param_data %in% param_default_values == FALSE],collapse=" "),"\n",
          "missing \'@contrast\' : ",paste(param_default_values[param_default_values %in% param_data == FALSE],collapse=" "),"\n")
    }
    
    if(ncol(object@fieldDim)!=3 || any(names(object@fieldDim)!=c("i","j","k")))
    { stop("validity[MRIaggr] : wrong specification of \'fieldDim\' \n",
           "required names : \"i\" \"j\" \"k\" \n",
           "proposed names : ",names(object@fieldDim),"\n")   
    }
    
    if(ncol(object@voxelDim)!=4 || any(names(object@voxelDim)!=c("i","j","k","unit")))
    { stop("validity[MRIaggr] : wrong specification of \'voxelDim\' \n",
           "required names : \"i\" \"j\" \"k\" \"unit\" \n",
           "proposed names : ",names(object@voxelDim),"\n")   
    }
    
    if(!identical(names(object@midplane),c("i","j")))
    {stop("validity[MRIaggr] : names of \'midplane\' are incorrects \n",
          "valid names : \"i\" \"j\" \n",
          "proposed names : ",paste(names(object@midplane),collapse=" "),"\n")
    }
    
    if(identical(names(object@hemispheres),c("left","right"))==FALSE)
    {stop("validity[MRIaggr] : names of \'hemispheres\' are incorrects \n",
          "valid names : \"left\" \"right\" \n",
          "proposed names : ",paste(names(object@hemispheres),collapse=" "),"\n")
    }
    
    if(any(object@hemispheres %in% c("lesion","contralateral","defined","undefined")==FALSE))
    {stop("validity[MRIaggr] : \'hemispheres\' is incorrect \n",
          "valid components : \"lesion\" \"contralateral\" \"defined\" \"undefined\" (must be characters and not factors) \n",
          "proposed components : ",paste(object@hemispheres,collapse=" "),"\n")
    }
    
    # cat(" : valid MRIaggr  \n")
    return(TRUE)} 
  
  
)

# Initiateur

setMethod(
  f="initialize",
  signature="MRIaggr",
  definition=function(.Object,identifier,contrast,clinic,
                      fieldDim,voxelDim,default_value,history,
                      normalization,midplane,hemispheres,
                      table_lesion,table_reperfusion,table_hypoperfusion,table_mask,
                      ls_descStats)
  {    
    if(any(c("i","j","k") %in% names(contrast)==FALSE))
    { stop("initialize[MRIaggr] : wrong specification of \'contrast\' \n",
           "it must contains columns named \"i\" \"j\" \"k\" \n",
           "proposed names : ",paste(names(contrast),collapse=" "),"\n")   
    }
    
    if("index" %in% names(contrast)==FALSE){
      nom_data <- names(contrast)
      contrast <- data.frame(1:nrow(contrast),contrast)  
      names(contrast) <- c("index",nom_data)
    }else{
      contrast$index <- 1:nrow(contrast)
    }
    
    
    if(missing(identifier)){
      stop("initialize[MRIaggr] : \'@identifier\' must be specified \n")
    }
    
    .Object@identifier <- identifier
    
    if(missing(contrast)){
      stop("initialize[MRIaggr] : \'@contrast\' must be specified \n")
    }
    .Object@contrast <- contrast
    
    if(!missing(clinic))
    {.Object@clinic <- clinic}
    
    if(!missing(table_lesion))
    {.Object@table_lesion <- table_lesion}
    
    if(!missing(table_reperfusion))
    {.Object@table_reperfusion <- table_reperfusion}
    
    if(!missing(table_hypoperfusion))
    {.Object@table_hypoperfusion <- table_hypoperfusion}
    
    if(!missing(table_mask))
    {.Object@table_mask <- table_mask}
    
    if(!missing(ls_descStats))
    {.Object@ls_descStats <- ls_descStats}else
    {.Object@ls_descStats <- list()}
    
    if(!missing(history))
    {.Object@history <- history}else
    {.Object@history <- list()}
    
    if(!missing(midplane))
    {.Object@midplane <- midplane}else
    {.Object@midplane <- data.frame(i=NA,j=NA)}
    
    if(!missing(hemispheres))
    {.Object@hemispheres <- hemispheres}else
    {.Object@hemispheres <- data.frame(left="undefined",right="undefined",stringsAsFactors = FALSE)}    
    
    if(!missing(normalization))
    {.Object@normalization <- normalization}
    
    if(!missing(fieldDim))
    {.Object@fieldDim <- fieldDim}
    
    if(!missing(voxelDim))
    {.Object@voxelDim <- voxelDim}
    
    if(!missing(default_value))
    {.Object@default_value <- default_value}
    
    validObject(.Object)
    return(.Object)
  }
)

##### B) Selecters ############################################

# selectContrast
setMethod(f ="selectContrast",
          signature ="MRIaggr",
          definition = function(object,param=NULL,num=NULL,format="data.frame",slice_var="k",
                                coords=FALSE,hemisphere="both",norm_mu=FALSE,norm_sigma=FALSE,
                                na.rm=FALSE,subset=NULL)
          { 
            #### initialisation ####
            # coords
            if(identical(coords,TRUE)){coords <- c("i","j","k")}
            if(identical(coords,FALSE)){coords <- NULL}
            
            # param
            param <- c(coords,param)
            param <- initParameter(object=object,param=param,test=TRUE,init=TRUE,accept.coords=TRUE,
                                   method="selectContrast")  
            paramMRI <- param[param %in% c("index","i","j","k")==FALSE]
            
            
            # subset
            if(!is.null(subset) && length(subset)==1 && is.character(subset)){
              
              initParameter(object=object,param=subset,test=TRUE,init=FALSE,accept.coords=FALSE,
                            arg_name="subset",long_name="subset",method="selectContrast")              
              subset <- object@contrast[,subset]
              
              if(is.logical(subset)==FALSE){
                stop("selectContrast[MRIaggr] : \'subset\'is not of type logical \n",
                     "type of \'subset\' : ",paste(is(subset),collapse=" "),"\n")
              }
              
              subset <- which(subset)
              
            }
            
            # hemisphere
            if(hemisphere %in% c("left","right","both") == FALSE){
              hemisphere <- selectHemispheres(object,hemisphere=hemisphere)
            }
            
            # norm
            if(is.null(norm_mu) || norm_mu %in% c(FALSE,"global","contralateral","global_1slice","contralateral_1slice","global_3slices","contralateral_3slices","default_value") ==FALSE)
            {stop("selectContrast[MRIaggr] : wrong specification of \'norm_mu\' \n",
                  "valid normalisations  : FALSE \"global\" \"contralateral\" \"c_global\" \"c_contralateral\" \"c3_global\" \"c3_contralateral\" \"default_value\" \n",
                  "requested normalisation : ",norm_mu,"\n")}
            
            if(is.null(norm_sigma) || norm_sigma %in% c(FALSE,"global","contralateral","global_1slice","contralateral_1slice","global_3slices","contralateral_3slices") ==FALSE)
            {stop("selectContrast[MRIaggr] : wrong specification of \'norm_sigma\' \n",
                  "valid normalisations  : FALSE \"global\" \"contralateral\" \"c_global\" \"c_contralateral\" \"c3_global\" \"c3_contralateral\" \n",
                  "requested normalisation : ",norm_sigma,"\n")}
            
            # num
            num <- initNum(object,num,test=TRUE,init=TRUE,slice_var=slice_var,method="selectContrast")
            
            #### tests ####
            if( (hemisphere!="both") && ("hemisphere" %in% names(object@contrast) == FALSE) )
            {stop("selectContrast[MRIaggr] : \'hemisphere\' is missing in \'object\' \n",
                  "use the allocHemisphere function to allocate it to \'object\' (and calcHemisphere function if you need to compute it) \n",
                  "the function must be called with the default \'hemisphere\' argument (\"both\") \n")}
            
            if( ( format %in% c("data.frame","matrix","vector")  ) == FALSE)
            {stop("selectContrast[MRIaggr] : wrong specification of \'format\' \n",
                  "valid formats : \"data.frame\" \"matrix\" \"vector\" \n",
                  "requested format : ",format,"\n")}
            
            if( (format == "vector") && (length(param)>1) ){
              stop("selectContrast[MRIaggr] : wrong specification of \'format\' \n",
                   "vector format is not available when several parameters are requested \n",
                   "requested parameters : ",paste(param,collpase=" "),"\n")}
            
            if( hemisphere %in% c("both","right","left") == FALSE)
            {stop("selectContrast[MRIaggr] : wrong specification of \'hemisphere\'  \n",
                  "valid values : \"lesion\" \"contralateral\" \"both\"  \right\" \"left\" \n",
                  "requested value : ",hemisphere,"\n")}
            
            ## data
            res <- object@contrast[,c("index","i","j","k",paramMRI),drop=FALSE]
            
            ## Selection des sites
            test <- rep(TRUE,nrow(object@contrast))
            
            # selection par index
            if(!is.null(subset))
            { 
              if(any( subset %in% object@contrast$index == FALSE)){
                warning("selectContrast[MRIaggr] : proposed \'index\' do not match object@contrast$index \n")
              }
              test_index <- object@contrast$index %in% subset
              test <- test*test_index 
            }
            
            # selection de la cartographie
            if(!is.null(num))
            {test_k <- object@contrast[,slice_var] %in% num
             test <- test*test_k }
            
            # elimination des NA
            if(na.rm)
            { test_Nna <- rowSums(is.na(res))==0
              index_Nna <- which(test_Nna)
              test <- test*test_Nna
            }
            
            # selection de l hemisphere
            if(hemisphere=="right")
            {test_hemi <- object@contrast$hemisphere=="right" 
             test <- test*test_hemi}
            if(hemisphere=="left")
            {test_hemi <- object@contrast$hemisphere=="left"
             test <- test*test_hemi}
            
            # application de la selection
            res <- res[as.logical(test),]
            
            if(norm_mu=="default_value"){
              default_value <- selectDefault_value(object,param=paramMRI,as.numeric=TRUE)              
              res[,paramMRI] <- sweep(res[,paramMRI,drop=FALSE],2,
                                      default_value,FUN="-")
            }
            
            # normalisation
            if( norm_mu %in% c(FALSE,"default_value") == FALSE || norm_sigma!=FALSE){
              
              if(length(selectNormalization(object))==0)
              {stop("selectContrast[MRIaggr] : normalization values are missing in \'object\' \n",
                    "use the calcNormalization  function to compute it and the allocNormalization function to allocate it to 'object' \n")}
              
              if(norm_mu=="global"){
                res[,paramMRI] <-  sweep(res[,paramMRI,drop=FALSE],2,
                                         as.matrix(selectNormalization(object,type="global",mu=TRUE,sigma=FALSE,param=paramMRI,hemisphere="both")),FUN="-")
              } 
              if(norm_mu=="global_1slice"){
                for(iter_k in unique(res$k)){
                  index_k <- which(res$k==iter_k)
                  res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                  as.matrix(selectNormalization(object,type="slice",mu=TRUE,sigma=FALSE,param=paramMRI,hemisphere="both",num=iter_k)),FUN="-")
                }
              }
              if(norm_mu=="global_3slices"){
                for(iter_k in unique(res$k)){
                  index_k <- which(res$k==iter_k)
                  res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                  as.matrix(selectNormalization(object,type="3slices",mu=TRUE,sigma=FALSE,param=paramMRI,hemisphere="both",num=iter_k)),FUN="-")                  
                }
              }              
              
              if(norm_mu %in% c("contralateral","contralateral_1slice","contralateral_3slices") || norm_sigma %in% c("contralateral","contralateral_1slice","contralateral_3slices")){
                
                hemi <- merge(x=res[,names(res)!=c("hemisphere")],
                              y=object@contrast[,c("i","j","k","hemisphere")],
                              by=c("i","j","k"),all.y=FALSE,sort=FALSE)$hemisphere
                
                index_hemiL <- which(hemi=="left")
                index_hemiR <- which(hemi=="right")
                
                if(any("undefined" %in% hemi)){
                  warning("selectContrast[MRIaggr] : the hemisphere is undefined for some observations : by default they are assigned to the left hemisphere \n")
                  index_hemiL <- c(index_hemiL,which(hemi=="undefined"))
                }
                
                if(length(index_hemiL)>0){
                  if(norm_mu=="contralateral"){
                    res[index_hemiL,paramMRI] <-  sweep(res[index_hemiL,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="global",mu=TRUE,sigma=FALSE,param=paramMRI,hemisphere="right")),FUN="-")                                                      
                  }
                  if(norm_mu=="contralateral_1slice"){
                    for(iter_k in unique(res$k)){
                      index_k <- intersect(which(res$k==iter_k),index_hemiL)
                      if(length(index_k)>0)
                        res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="slice",mu=TRUE,sigma=FALSE,param=paramMRI,hemisphere="right",num=iter_k)),FUN="-")
                    }
                  }
                  if(norm_mu=="contralateral_3slices"){
                    for(iter_k in unique(res$k)){
                      index_k <- intersect(which(res$k==iter_k),index_hemiL)
                      if(length(index_k)>0)
                        res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="3slices",mu=TRUE,sigma=FALSE,param=paramMRI,hemisphere="right",num=iter_k)),FUN="-")                      
                    }
                  }
                  
                  if(norm_sigma=="contralateral"){
                    res[index_hemiL,paramMRI] <-  sweep(res[index_hemiL,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="global",mu=FALSE,sigma=TRUE,param=paramMRI,hemisphere="right")),FUN="/")
                  }
                  
                  if(norm_sigma=="contralateral_1slice"){
                    for(iter_k in unique(res$k)){
                      index_k <- intersect(which(res$k==iter_k),index_hemiL)
                      if(length(index_k)>0)
                        res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="slice",mu=FALSE,sigma=TRUE,param=paramMRI,hemisphere="right",num=iter_k)),FUN="/")
                    }
                  }
                  if(norm_sigma=="contralateral_3slices"){
                    for(iter_k in unique(res$k)){
                      index_k <- intersect(which(res$k==iter_k),index_hemiL)
                      if(length(index_k)>0)
                        res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="3slices",mu=FALSE,sigma=TRUE,param=paramMRI,hemisphere="right",num=iter_k)),FUN="/")
                    }
                  }
                }
                
                
                if(length(index_hemiR)>0){
                  if(norm_mu=="contralateral"){
                    res[index_hemiR,paramMRI] <-  sweep(res[index_hemiR,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="global",mu=TRUE,sigma=FALSE,param=paramMRI,hemisphere="left")),FUN="-")
                  }
                  if(norm_mu=="contralateral_1slice"){
                    for(iter_k in unique(res$k)){
                      index_k <- intersect(which(res$k==iter_k),index_hemiR)
                      if(length(index_k)>0)
                        res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="slice",mu=TRUE,sigma=FALSE,param=paramMRI,hemisphere="left",num=iter_k)),FUN="-")
                    }
                  }
                  if(norm_mu=="contralateral_3slices"){
                    for(iter_k in unique(res$k)){
                      index_k <- intersect(which(res$k==iter_k),index_hemiR)
                      if(length(index_k)>0)
                        res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="3slices",mu=TRUE,sigma=FALSE,param=paramMRI,hemisphere="left",num=iter_k)),FUN="-")
                    }
                  }
                  
                  if(norm_sigma=="contralateral"){
                    res[index_hemiR,paramMRI] <-  sweep(res[index_hemiR,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="global",mu=FALSE,sigma=TRUE,param=paramMRI,hemisphere="left")),FUN="/")
                  }
                  if(norm_sigma=="contralateral_1slice"){
                    for(iter_k in unique(res$k)){
                      index_k <- intersect(which(res$k==iter_k),index_hemiR)
                      if(length(index_k)>0)
                        res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="slice",mu=FALSE,sigma=TRUE,param=paramMRI,hemisphere="left",num=iter_k)),FUN="/")
                    }
                  }
                  if(norm_sigma=="contralateral_3slices"){
                    for(iter_k in unique(res$k)){
                      index_k <- intersect(which(res$k==iter_k),index_hemiR)
                      if(length(index_k)>0)
                        res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                        as.matrix(selectNormalization(object,type="3slices",mu=FALSE,sigma=TRUE,param=paramMRI,hemisphere="left",num=iter_k)),FUN="/")
                    }
                  }
                }
              }
            }
            
            if(norm_sigma=="global"){
              res[,paramMRI] <-  sweep(res[,paramMRI,drop=FALSE],2,
                                       as.matrix(selectNormalization(object,type="global",mu=FALSE,sigma=TRUE,param=paramMRI,hemisphere="both")),FUN="/")
            }
            if(norm_sigma=="global_1slice"){
              for(iter_k in unique(res$k)){
                index_k <- which(res$k==iter_k)
                res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                as.matrix(selectNormalization(object,type="slice",mu=FALSE,sigma=TRUE,param=paramMRI,hemisphere="both",num=iter_k)),FUN="/")
              }
            }
            if(norm_sigma=="global_3slices"){
              for(iter_k in unique(res$k)){
                index_k <- which(res$k==iter_k)
                res[index_k,paramMRI] <-  sweep(res[index_k,paramMRI,drop=FALSE],2,
                                                as.matrix(selectNormalization(object,type="3slices",mu=FALSE,sigma=TRUE,param=paramMRI,hemisphere="both",num=iter_k)),FUN="/")
              }
            }
            
            ## Export
            # gestion des coordonnees
            res <- res[,param]
            
            if(format=="vector")
            {return(res)}
            
            if(format=="data.frame")
            { res <- as.data.frame(res)
              names(res) <- param
              return(res)}
            
            if(format=="matrix")
            {res <- as.matrix(res)
             colnames(res) <- NULL
             return(res)}
            
          }
          
)

# selectCoords
setMethod(f ="selectCoords",
          signature ="MRIaggr",
          definition = function(object,coords=c("i","j","k"),spatial_res=c(1,1,1),num=NULL,hemisphere="both",
                                subset=NULL,slice_var="k",format="data.frame")
          { 
            if(any(coords %in% c("i","j","k","index")==FALSE)){
              stop("selectCoords[MRIaggr] : wrong specification of \'coords\' \n",
                   "valid values : \"i\" \"j\" \"k\" or \"index\" \n",
                   "requested values : ",paste(coords,sep=" "),"\n")
            }
            
            if(length(spatial_res)!=3 || is.numeric(spatial_res)==FALSE){
              stop("selectCoords[MRIaggr] : wrong specification of \'spatial_res\' \n",
                   "\'spatial_res\' must be a numeric vector of length 3 \n",
                   "length(spatial_res) : ",length(spatial_res),"\n",
                   "is(spatial_res) : ",paste(is(spatial_res),collapse=""),"\n")
            }
            
            res <- selectContrast(object,param=coords,num=num,hemisphere=hemisphere,
                                  format=format,subset=subset,slice_var=slice_var)
            
            if(spatial_res[1]!=1 && "i" %in% coords){res$i <- res$i * spatial_res[1]}
            if(spatial_res[2]!=1 && "j" %in% coords){res$j <- res$j * spatial_res[2]}
            if(spatial_res[2]!=1 && "k" %in% coords){res$k <- res$k * spatial_res[3]}
            
            return(res)
            
          }
)

setMethod(f ="selectClinic",
          signature ="MRIaggr",
          definition = function(object,param=NULL)
          {
            if(is.null(param)){
              return(object@clinic)
            }
            
            if(any(param %in% names(object@clinic)==FALSE)){
              dist_string <- adist(param[param %in% names(object@clinic)==FALSE],names(object@clinic))
              propositions <- as.vector(unlist(apply(dist_string,1,function(y){names(object@clinic)[which(rank(y)<5)]})))
              stop("selectClinic[MRIaggr] : wrong specification of \'param\' \n",
                   "available parameters  : ",paste(propositions,collapse=" ")," ... \n",
                   "requested parameter : ",paste(param[param %in% names(object@clinic)==FALSE],collapse=" "),"\n")
            }else{
              return(object@clinic[,param,drop=FALSE])
            }
            
            
          }
)

# selectDefault_value
setMethod(f ="selectDefault_value",
          signature ="MRIaggr",
          definition = function(object,param=NULL,as.numeric=FALSE)
          { 
            
            param <- initParameter(object=object,param=param,test=TRUE,init=TRUE,
                                   accept.coords=FALSE,accept.mask=FALSE,accept.index=FALSE,method="selectDefault_value")
            
            default_value <- object@default_value[param]
            if(as.numeric){             
              default_value <- as.numeric(default_value)
              names(default_value) <- param
            }
            
            return(default_value) }
)

# selectDescStats
setMethod(f ="selectDescStats",
          signature ="MRIaggr",
          definition = function(object,name=NULL,
                                subset_W=NULL,hemisphere="both",num=NULL,slice_var="k")
          {  
            if(is.null(name))
            {return(object@ls_descStats) }
            
            if(length(name)>1)
            {stop("selectDescStats[MRIaggr] : \'name\' must be of length 1 \n",
                  "length of the requested names : ",length(name),"\n")}
            
            if(name %in% names(object@ls_descStats) == FALSE)
            {stop("selectDescStats[MRIaggr] : \'name\' not available \n",
                  "available elements in @ls_descStats : ",paste(names(object@ls_descStats),collapse=" "),"\n",
                  "requested element : ",name,"\n")}
            
            if(name %in% c("W_euclidean") == FALSE){
              return(object@ls_descStats[[name]])  
            }
            
            if(name %in% c("W_euclidean")){
              #### test package W
              test.package <- requireNamespace("spam",quietly=TRUE)
              if(test.package==FALSE){
                stop("selectDescStats[MRIaggr] : this function with argument name=\"W_euclidean\" requires to have installed the spam package to work \n")
              }
              
              test.package <- requireNamespace("Matrix",quietly=TRUE)
              if(test.package==FALSE){
                stop("selectDescStats[MRIaggr] : this function with argument name=\"W_euclidean\" requires to have installed the Matrix package to work \n")
              }
              
              if(is.null(subset_W)){
                subset_W <- selectContrast(object,param="index",hemisphere=hemisphere,num=num,slice_var="k",format="vector")
              }
              
              if(min(subset_W) < 1 || max(subset_W) > selectN(object))
              {stop("selectDescStats[MRIaggr] : wrong specification of \'subset_W\' \n",
                    "dim W : ",selectN(object),"\n",
                    "range(subset_W) : ",paste(range(subset_W),collapse=" "),"\n")}
              
              res <- object@ls_descStats[["W_euclidean"]][subset_W,subset_W]
              
              return(res)
            }
            
          }
)

# selectHemispheres
setMethod(f ="selectHemispheres",
          signature ="MRIaggr",
          definition = function(object,hemisphere="both"){
            
            #### tests
            if(length(hemisphere)!=1){
              stop("selectHemispheres[MRIaggr] :  wrong specification of \'hemisphere\' \n",
                   "must be a single character \n",
                   "length(hemisphere) : ",length(hemisphere),"\n")
            }
            
            if(any(hemisphere %in% c("both","left","right","lesion","contralateral") == FALSE)){
              stop("selectHemispheres[MRIaggr] :  wrong specification of \'hemisphere\' \n",
                   "valid hemisphere values : \"both\" \"left\" \"right\" \"lesion\" \"contralateral\" \n",
                   "proposed \'hemisphere\' : ",hemisphere,"\n")
            }
            
            #### extraction
            
            if(hemisphere=="both"){
              return(object@hemispheres)
            }
            
            if(hemisphere=="left"){
              return(object@hemispheres[,"left"])
            }
            
            if(hemisphere=="right"){
              return(object@hemispheres[,"right"])
            }
            
            if(hemisphere=="lesion"){ 
              hemiLesion <- object@hemispheres=="lesion"
              
              if(sum(hemiLesion)==2){hemisphere <- "both"}else
                if(sum(hemiLesion)==0){hemisphere <- NULL}else
                  if(hemiLesion[,"left"]==TRUE){hemisphere <- "left"}else{hemisphere <- "right"}
              
              return(hemisphere)
            }
            
            if(hemisphere=="contralateral"){ 
              hemiContro <- object@hemispheres=="contralateral"
              
              if(sum(hemiContro)==2){hemisphere <- "both"}else
                if(sum(hemiContro)==0){hemisphere <- NULL}else
                  if(hemiContro[,"left"]==TRUE){hemisphere <- "left"}else{hemisphere <- "right"}
              
              return(hemisphere)
            }
            
          }
)

# selectHistory
setMethod(f ="selectHistory",
          signature ="MRIaggr",
          definition = function(object)
          { return(object@history)     }
)

# selectIdentifier
setMethod(f ="selectIdentifier",
          signature ="MRIaggr",
          definition = function(object)
          { return(object@identifier)     }
)

# selectMidplane
setMethod(f ="selectMidplane",
          signature ="MRIaggr",
          definition = function(object)
          {return(object@midplane) }
)

# selectN
setMethod(f ="selectN",
          signature ="MRIaggr",
          definition = function(object,num=NULL,hemisphere="both",subset=NULL)
          { return(nrow(selectCoords(object,num=num,subset=subset,hemisphere=hemisphere))) 
          }
) 

# selectNormalization
setMethod(f ="selectNormalization",
          signature ="MRIaggr",
          definition = function(object,type=NULL,mu=TRUE,sigma=TRUE,hemisphere="both",num=NULL,param=NULL)
          {
            if(!is.null(type) && type %in% c("global","slice","3slices")==FALSE){
              stop("selectNormalization[MRIaggr] : wrong specification of \'type\' \n",
                   "valid types : NULL \"global\" \"slice\" \"3slices\" \n",
                   "proposed type : ",type,"\n")
            }
            
            if(is.null(type))
            {return(object@normalization)}
            
            #             if(hemisphere == "both" && type != "global"){
            #               stop("selectNormalization[MRIaggr] : wrong specification of \'hemisphere\' \n",
            #                    "if type is not \"global\" hemisphere cannot be \"both\" \n",
            #                    "requested hemisphere ; type : ",hemisphere," ; ",type,"\n")
            #             }
            
            if( all( c(mu,sigma) == FALSE) ){
              stop("selectNormalization[MRIaggr] : \'mu\' and \'sigma\' cannot both be FALSE \n",
                   "requested mu ; sigma : ",mu," ; ",sigma,"\n")
            }
            
            if(hemisphere %in% c("both","left","right")==FALSE){
              stop("selectNormalization[MRIaggr] : wrong specification of \'hemisphere\' \n",
                   "valid hemisphere : \"both\" \"left\" \"right\" \n",
                   "requested hemisphere : ",hemisphere,"\n")
            }
            
            
            ## type global
            if(type=="global"){
              
              res_norm <- object@normalization$norm_global
              
              if(is.null(hemisphere)){hemisphere <- c("both","left","right")}
              if(hemisphere=="left"){hemisphere <- "left"}
              if(hemisphere=="right"){hemisphere <- "right"}
              
              stat <- NULL
              if(mu){stat <- c(stat,"mu")}
              if(sigma){stat <- c(stat,"sigma")}
              
              nom_lignes <- as.vector(sapply(stat,function(object){paste(object,"_",hemisphere,sep="")}))
              res_norm <- res_norm[nom_lignes,,drop=FALSE]
              
              if(!is.null(param)){
                if(any(param %in% names(res_norm)==FALSE)){
                  stop("selectNormalization[MRIaggr] : \'param\' not available for normalization \n",
                       "available paramaters : ",paste(names(res_norm),collapse=" ")," \n",
                       "requested but not available parameter",if(sum(param %in% names(res_norm)==FALSE)>1){"s"}," : ",
                       paste(param[param %in% names(res_norm)==FALSE],collapse=" "),"\n")                  
                }                
                res_norm <- res_norm[,param,drop=FALSE]
              }
              
              return(res_norm)
              
            }else{
              
              if( mu==TRUE && sigma==TRUE ){
                stop("selectNormalization[MRIaggr] :  arguments \'mu\' and \'sigma\' cannot be simultaneously TRUE \n",
                     "requested mu ; sigma : ",mu," ; ",sigma,"\n")
              }
              
              if(mu==TRUE){stat <- "Mu"}
              if(sigma==TRUE){stat <- "Sigma"}
              nom <- paste("norm",stat,"_",type,"_",hemisphere,sep="")
              
              if( nom %in% names(object@normalization) == FALSE ){
                stop("selectNormalization[MRIaggr] : normalisation not available \n",
                     "available normalisations : ",paste(names(object@normalization),collapse=" "),"\n",
                     "not available normalisation : ",nom[nom %in% names(object@normalization) == FALSE],"\n")
              }
              
              res_norm <- object@normalization[[nom]]
              
              if(!is.null(param)){
                if(any(param %in% names(res_norm)==FALSE)){
                  stop("selectNormalization[MRIaggr] : wrong specification of \'param\' \n",
                       "available parameters : ",paste(names(res_norm),collapse=" ")," \n",
                       "requested but not available parameter",if(sum(param %in% names(res_norm)==FALSE)>1){"s"}," : ",
                       paste(param[param %in% names(res_norm)==FALSE],collapse=" "),"\n")                  
                }                
                res_norm <- res_norm[,param,drop=FALSE]
              }
              
              if(!is.null(num)){
                if(any(num %in% (1:nrow(res_norm))==FALSE)){
                  stop("selectNormalization[MRIaggr] :  wrong specification of \'num\' \n",
                       "available slices : ",paste(1:nrow(res_norm),collapse=" ")," \n",
                       "requested but not available slice",if(sum(param %in% names(res_norm)==FALSE)>1){"s"}," : ",
                       paste(num[num %in% (1:nrow(res_norm))==FALSE],collapse=" "),"\n")                  
                }                
                res_norm <- res_norm[num,,drop=FALSE]
              }
              
              return(res_norm)
            }
            
          }
)

# selectParameter
setMethod(f ="selectParameter",
          signature ="MRIaggr",
          definition = function(object,type="contrast",mask=TRUE)
          { 
            if(type %in% c("clinic","ls_descStats", "contrast") ==FALSE)
            {stop("selectParameter[MRIaggr] : wrong specification of \'type\' \n",
                  "valid type : \"clinic\" \"ls_descStats\" \"contrast\" \n",
                  "proposed type : ",type ,"\n")}
            
            if(type=="clinic")
            {return(names(object@clinic))}
            
            if(type=="ls_descStats")
            {return(names(object@ls_descStats))}
            
            if(type=="contrast"){
              if(mask){
                index_coords <- which(names(object@contrast) %in% c("i","j","k","index"))
              }else{
                index_coords <- which(names(object@contrast) %in% c("i","j","k","index","mask"))
              }
              return(names(object@contrast)[-index_coords])  
            }
          }
)

# selectTable
setMethod(f ="selectTable",
          signature ="MRIaggr",
          definition = function(object,type,size=FALSE)
          {            
            if( type %in% c("lesion","reperfusion","hypoperfusion") == FALSE )
            {  stop("selectTable[MRIaggr] :  wrong specification of \'type\' \n",
                    "valid types  : \'lesion\' \'reperfusion\' \'hypoperfusion\' \n",
                    "proposed type : ",type,"\n")
            }    
            
            res <- switch(type,
                          "lesion"=object@table_lesion,
                          "reperfusion"=object@table_reperfusion,
                          "hypoperfusion"=object@table_hypoperfusion)
            if(size==TRUE){
              Volume <- prod(selectVoxelDim(object)[1:3])
              if(type=="lesion"){
                res <- res*Volume
              }else{                
                test.V <- sapply(names(res),function(x){substr(x,1,1)})=="V"
                res[,test.V] <- res[,test.V]*Volume                
              }
            }
            
            return(res) 
          }
)

# selectFieldDim
setMethod(f ="selectFieldDim",
          signature ="MRIaggr",
          definition = function(object)
          {return(object@fieldDim) }
)

# selectVoxelDim
setMethod(f ="selectVoxelDim",
          signature ="MRIaggr",
          definition = function(object,unit=TRUE)
          {
            if(unit==TRUE){
              return(object@voxelDim) 
            }else{
              return(object@voxelDim[,c("i","j","k"),drop=FALSE]) 
            }
          }
)

##### C) Allocators ############################################

setReplaceMethod(f ="allocContrast",
                 signature ="MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                 definition = function(object,param=NULL,default_value=NULL,overwrite=FALSE,trace=TRUE,value)
                 { 
                   #### verification optionnelle des coordonnees ####
                   value <- as.data.frame(value)
                   index_coords <- which(names(value) %in% c("i","j","k"))
                   
                   if(length(index_coords)==3)
                   {if(any(value$i != object@contrast$i) || any(value$j != object@contrast$j) || any(value$k != object@contrast$k))
                   {stop("allocContrast[MRIaggr] : coordinates do not match between \'value\' and \'object@contrast\' \n",
                         "number of mismatch by coordinate (i ; j ; k) : ",
                         sum(value$i != object@contrast$i),";",
                         sum(value$j != object@contrast$j),";",
                         sum(value$k != object@contrast$k),"\n")
                   }}
                   
                   if(length(index_coords)>0){value <- value[,-index_coords,drop=FALSE]}
                   
                   #####                    
                   if(is.null(param)){
                     param <- names(value)
                   }else{
                     if(length(param) != ncol(value))
                     {stop("allocContrast[MRIaggr] : mismatch between \'param\' and the number of parameters in \'value\' \n",
                           "length(param) : ",length(param),"(",paste(param,collapse=" "),") \n",
                           "ncol(value) : ",ncol(value),"(",paste(names(value),collapse=" "),") \n")
                     }
                     
                     names(value) <- param
                   }
                   
                   if("index" %in% param){
                     stop("allocContrast[MRIaggr] : cannot allocate a parameter with name \"index\" \n",
                          "\"index\" is a reserved name \n",
                          "column of concern : ",which(param == "index"),"\n")   
                   }
                   
                   if("mask" %in% param){
                     if(dim(value)[2]>1)
                     {stop("allocContrast[MRIaggr] : when assigning the mask only one parameter can be considered \n",
                           "ncol(value) : ",ncol(value),"\n")
                     }
                     
                     if(is.logical(value[,1])==FALSE && is.integer(value[,1])==FALSE)
                     {stop("allocContrast[MRIaggr] : the \'mask\' must be of type logical or integer \n",
                           "is(value)  : ",paste(is(value[,1]),collapse=" "),"\n")
                     }
                   }
                   
                   if("hemisphere" %in% param){                     
                     table_tempo <- table(value[,"hemisphere"])
                     if(length(names(table_tempo))>3 || any(names(table_tempo) %in% c("left","right","undefined") == FALSE)){
                       stop("allocContrast[MRIaggr] : wrong specification of the value of the \"hemisphere\" parameter \n",
                            "only three values are allowed : \"left\" \"right\" \"undefined\" \n",
                            "proposed values  : ",paste(names(table_tempo),collapse=" "),"\n")
                     }
                   }
                   
                   if(nrow(value) != nrow(object@contrast))
                   {stop("allocContrast[MRIaggr] : mismatch between the dimension of \'value\' and \'object@contrast\'  \n",
                         "nrow(object@contrast) : ",nrow(object@contrast),"\n",
                         "nrow(value) : ",nrow(value),"\n")
                   }
                   
                   if(sum(param %in% names(object@contrast))>0 && overwrite==FALSE){ 
                     stop("allocContrast[MRIaggr] : -  some names of \'param\' are already present in \'object@contrast\' \n",
                          "redundant names : ",paste(param[which(param %in% names(object@contrast))],collapse=" "),"\n",
                          "set \'overwrite\' to TRUE to perform this operation \n")
                   }
                   
                   if("mask" %in% param ==FALSE && is.null(default_value)){
                     default_value <- data.frame(matrix("undefined",ncol=length(param),nrow=1))
                     names(default_value) <- param
                   }
                   
                   if("mask" %in% param ==FALSE && (is.data.frame(default_value) == FALSE || sum(names(default_value)!=param)>0 ) ){
                     stop("allocContrast[MRIaggr] : wrong specification of \'default_value\' \n",
                          "\'default_value\' must be a data.frame with names : \"",paste(param,collapse="\" \""),"\" \n",
                          "is(default_value): ",paste(is(default_value),collpase=" "),"\n",
                          "names(default_value) : \"",paste(names(default_value),collapse="\" \""),"\"\n")
                   }
                   
                   param_commun <- selectParameter(object)[selectParameter(object) %in% param]
                   param_nouveaux <- colnames(value)[param %in% selectParameter(object) == FALSE]
                   
                   if(length(param_commun)>0){
                     for(iter_param in param_commun){
                       object@contrast[,iter_param] <- value[,iter_param,drop=FALSE]
                       if(iter_param != "mask"){
                         object@default_value[,iter_param] <- default_value[,iter_param,drop=FALSE]
                       }
                     }
                   }
                   
                   if(length(param_nouveaux)>0){
                     for(iter_param in param_nouveaux){
                       object@contrast <- data.frame(object@contrast,value[,iter_param,drop=FALSE],stringsAsFactors=FALSE)
                       if("mask" %in% param ==FALSE){
                         object@default_value <- data.frame(object@default_value,default_value[,iter_param,drop=FALSE],stringsAsFactors=FALSE)
                       }
                     }
                   }
                   
                   validObject(object)
                   
                   if(trace){
                     
                     if(length(param_commun)==1){                    
                       cat("allocContrast[MRIaggr] : Cartography \"",param_commun,"\" \n",
                           "                         has been updated \n",sep="")                      
                     }
                     if(length(param_commun)>1){                    
                       cat("allocContrast[MRIaggr] : Cartographies \"",paste(param_commun,collapse="\" \""),"\" \n",
                           "                         have been updated \n",sep="")                      
                     }
                     
                     if(length(param_nouveaux)==1){                    
                       cat("allocContrast[MRIaggr] : Cartography \"",param_nouveaux,"\" \n",
                           "                         has been allocated \n",sep="")                      
                     }
                     if(length(param_nouveaux)>1){                    
                       cat("allocContrast[MRIaggr] : Cartographies \"",paste(param_nouveaux,collapse="\" \""),"\" \n",
                           "                         have been allocated \n",sep="")                      
                     }                 
                     
                   }
                   
                   return(object)
                 }
)

setReplaceMethod(f = "allocClinic",
                 signature = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                 definition = function(object,add=FALSE,overwrite=FALSE,trace=TRUE,value)
                 {    
                   if(!is.logical(add)){
                     stop("allocClinic[MRIaggr] : wrong specification of \'add\' \n",
                          "\'add\' must be logical \n",
                          "proposed \'add\' : ",add," \n")
                   }
                   
                   if(!is.logical(overwrite)){
                     stop("allocClinic[MRIaggr] : wrong specification of \'overwrite\' \n",
                          "\'overwrite\' must be logical \n",
                          "proposed \'overwrite\' : ",overwrite," \n")
                   }
                   
                   if(!is.data.frame(value)){
                     stop("allocClinic[MRIaggr] : wrong specification of \'value\' \n",
                          "\'value\' must be a data.frame \n",
                          "is(value) : ",paste(is(value),collapse=" ")," \n")
                   }
                   
                   if(ncol(object@clinic)==0){add <- FALSE}
                   
                   if(add==FALSE){
                     if(trace){sauveClinic <- object@clinic}
                     
                     object@clinic <- value
                     validObject(object)
                     
                     if(trace){
                       cat("allocClinic[MRIaggr] : @clinic has been ",
                           if(sum(!is.na(sauveClinic))==0){"allocated"}else{"updated"},"\n",sep="")                                       
                     }
                     
                   }else{
                     names_alloc <- names(value)[names(value) %in% names(object@clinic)==FALSE]
                     names_replace <- names(value)[names(value) %in% names(object@clinic)]  
                     
                     if(length(names_alloc)>0){                       
                       object@clinic <- data.frame(object@clinic,value)                       
                     }
                     if(length(names_replace)>0){
                       
                       if(overwrite==FALSE){
                         stop("allocClinic[MRIaggr] : \'value\' contains elements already present in \'@clinic\' \n",
                              "clinical parameter",if(length(names_replace)==1){"s"}," already present in @clinic : ",paste(names_replace,collapse=" "),"\n",
                              "set \'overwrite\' to TRUE to overwrite them \n")
                       }
                     }else{
                       object@clinic[,names_replace] <- value[,names_replace]
                     }
                     
                     validObject(object)
                     if(length(names_alloc)>0){
                       cat("allocClinic[MRIaggr] : parameter",if(length(names_alloc)>1){"s"},
                           " \"",paste(names_alloc,collapse="\" \""),"\" \n",
                           "                       ",if(length(names_alloc)==1){"has"}else{"have"}," been added to @clinic \n",sep="")
                     }
                     if(length(names_replace)>0){
                       cat("allocClinic[MRIaggr] : parameter",if(length(names_replace)>1){"s"},
                           " \"",paste(names_replace,collapse="\" \""),"\" \n",
                           "                       ",if(length(names_replace)==1){"has"}else{"have"}," been updated in @clinic \n",sep="")
                     }
                     
                   }
                   
                   return(object)
                 }
)

setReplaceMethod(f ="allocDescStats",
                 signature ="MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                 definition = function(object,name,overwrite=FALSE,trace=TRUE,value)
                 {    
                   if(length(name)>1)
                   {stop("allocDescStats[MRIaggr] : Only one element can be allocated at a time \n",            
                         "length(name) : ",name,"\n")
                   }                   
                   
                   if(name %in% "W_euclidean"){
                     
                     if("dgCMatrix" %in% class(value) == FALSE){
                       stop("allocDescStats[MRIaggr] : wrong specification of \'value\' \n",
                            "\'value\' must be of class dgCMatrix when \'name\'= \"W_euclidean\" \n",
                            "type of value : ",paste(is(value),collapse=" "),"\n")
                     }
                     
                     if(nrow(value)!=selectN(object) || ncol(value)!=selectN(object)){
                       stop("allocDescStats[MRIaggr] : wrong specification of \'value\' \n",
                            "\'value\' have the same dimension than \'object\' when \'name\'= \"W_euclidean\" : ",selectN(object)," ",selectN(object),"\n",
                            " dimension of value : ",paste(dim(value),collapse=" "),"\n")
                     }                     
                   }
                   
                   test.overwrite <- name %in% selectParameter(object,type="ls_descStats")
                   
                   if(sum(test.overwrite)>0 && overwrite==FALSE)
                   {stop("allocDescStats[MRIaggr] : The requested field(s) already exist in object@ls_descStats \n",
                         "Set \'overwrite\' to TRUE to replace the field \n",
                         "already existing fields : ",name[test.overwrite],"\n")
                   }
                   
                   param_commun <- names(object@ls_descStats)[(names(object@ls_descStats) %in% name)]
                   param_nouveaux <- name[(name %in% names(object@ls_descStats) == FALSE )]
                   
                   if(length(param_commun)>0){
                     object@ls_descStats[[param_commun]] <- value 
                   }
                   
                   if(length(param_nouveaux)>0){                    
                     object@ls_descStats <- c(object@ls_descStats,list(value))
                     names(object@ls_descStats)[length(object@ls_descStats)] <- param_nouveaux
                   }
                   
                   validObject(object)
                   
                   if(trace){
                     
                     if(length(param_commun)==1){
                       cat("allocDescStats[MRIaggr] : Element \"",param_commun,"\" \n",
                           "                          has been updated \n",sep="")
                     }
                     if(length(param_commun)>1){
                       cat("allocDescStats[MRIaggr] : Elements \"",paste(param_commun,collapse="\" \""),"\" \n",
                           "                          have been \n",sep="")
                     }
                     
                     if(length(param_nouveaux)==1){
                       cat("allocDescStats[MRIaggr] : Element \"",param_nouveaux,"\" \n",
                           "                          has been allocated \n",sep="")
                     }
                     if(length(param_nouveaux)>1){
                       cat("allocDescStats[MRIaggr] : Elements \"",paste(param_nouveaux,collapse="\" \""),"\" \n",
                           "                          have been allocated \n",sep="")
                     }                     
                   }
                   
                   return(object)
                 }
)

setReplaceMethod(f = "allocHemisphere",
                 signature = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                 definition = function(object,overwrite=FALSE,trace=TRUE,value)
                 {
                   if(trace){
                     sauveMidplane <- object@midplane
                     sauveHemispheres <- object@hemispheres
                   }
                   
                   
                   if(!is.list(value) || is.null(names(value)) || any(names(value) %in% c("midplane","hemispheres","data")==FALSE)){
                     stop("allocHemisphere[MRIaggr] : wrong specification of \'value\' \n",
                          "\'value\' must be a list containing some of the following elements \"midplane\" \"hemispheres\" \"data\" \n",
                          "incorrect names : ",paste(names(value)[names(value) %in% c("midplane","hemispheres","data")==FALSE],collapse=" "),"\n",
                          "is(value)  : ",paste(is(value),collapse=" "),"\n"
                     )
                   }
                   
                   if("midplane" %in% names(value)){
                     if(!is.null(object@midplane) && any(!is.na(object@midplane)) && overwrite==FALSE){
                       stop("allocHemisphere[MRIaggr] : midplane already existing in \'object\' \n",                          
                            "set \'overwrite\' to TRUE to replace it \n")
                     }
                     object@midplane <- value$midplane
                   }
                   
                   if("hemispheres" %in% names(value)){
                     if( (object@hemispheres$right!="undefined" || object@hemispheres$left!="undefined") && overwrite==FALSE){
                       stop("allocHemisphere[MRIaggr] : hemispheres already existing in \'object\' \n",                          
                            "set \'overwrite\' to TRUE to replace it \n")
                     }
                     object@hemispheres <- value$hemispheres
                   }
                   
                   if("data" %in% names(value)){
                     if(!is.data.frame(value$data) || nrow(value$data) != selectN(object)){
                       stop("allocHemisphere[MRIaggr] : wrong specification of the data element of \'value\' \n",                          
                            "data must be a column data.frame with ",selectN(object)," rows \n",
                            "type of data : ",paste(is(value$data),collapse=" "),"\n",
                            "number of rows : ",nrow(value$data),"\n")
                     }
                     
                     if(any(names(value$data) %in% c("hemisphere","i_hemisphere","j_hemisphere") == FALSE)){
                       stop("allocHemisphere[MRIaggr] : wrong specification of the data element of \'value\' \n",                          
                            "data must contains columns named \"hemisphere\" \"i_hemisphere\" \"j_hemisphere\" \n",
                            "proposed column names : ",paste(names(value$data),collapse=" "),"\n")
                     }
                     
                     allocContrast(object,overwrite=overwrite) <- value$data
                   }                 
                   
                   validObject(object)
                   
                   if(trace){
                     if("midplane" %in% names(value)){
                       cat("allocHemisphere[MRIaggr] : @midplane has been ")
                       if(sum(!is.na(sauveMidplane))==0){cat("alllocated \n")}else{cat("updated\n ")}
                     }
                     
                     if("hemispheres" %in% names(value)){
                       cat("allocHemisphere[MRIaggr] : @hemispheres has been ")
                       if(sum(sauveHemispheres!="undefined")==0){cat("allocated \n")}else{cat("updated\n ")}
                     }                       
                   }
                   
                   return(object)
                 }
)

setReplaceMethod(f = "allocNormalization",
                 signature = "MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                 definition = function(object,overwrite=FALSE,trace=TRUE,value)
                 {
                   if(trace){sauveNormalization <- object@normalization}
                   
                   if(overwrite==FALSE && !is.null(object@normalization) && length(object@normalization)>0 && any(!is.na(object@normalization))){
                     stop("allocNormalization[MRIaggr] : normalization already existing in \'object\' \n",                          
                          "set \'overwrite\' to TRUE to replace it \n")
                     
                   }
                   
                   valid_names <- c("norm_global",
                                    "normMu_slice_both","normSigma_slice_both",
                                    "normMu_slice_left","normSigma_slice_left",
                                    "normMu_slice_right","normSigma_slice_right",
                                    "normMu_3slices_both","normSigma_3slices_both",
                                    "normMu_3slices_left","normSigma_3slices_left",
                                    "normMu_3slices_right","normSigma_3slices_right") 
                   
                   if(is.null(names(value)) || any(names(value) %in% valid_names == FALSE)){
                     stop("allocNormalization[MRIaggr] : wrong specification of \'value\' \n",                          
                          "it must be a list containing normalization values in named data.frame \n",
                          "see the calcNormalization function for an example of valid value for the normalization slot \n")
                   }
                   object@normalization <- value
                   
                   validObject(object)
                   if(trace){
                     cat("allocNormalization[MRIaggr] : @normalization has been ")
                     if(length(sauveNormalization)==0){cat("allocated \n")}else{cat("updated\n ")}                    
                   }
                   return(object)
                 }
)

setReplaceMethod(f = "allocTable",
                 signature ="MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                 definition = function(object,type,overwrite=FALSE,trace=TRUE,value)
                 { 
                   if(length(type)!=1){
                     stop("selectTable[MRIaggr] : wrong specification of \'type\' \n",
                          "only one table can be allocated at the same time \n",
                          "length(type) : ",length(type),"\n")                     
                   }
                   
                   if(type %in% c("lesion","reperfusion","hypoperfusion")==FALSE){
                     stop("selectTable[MRIaggr] : wrong specification of \'type\' \n",
                          "valid type : \"lesion\" \"reperfusion\" \"hypoperfusion\" \"mask\" \n",
                          "proposed type : ",type,"\n")                     
                   }
                   
                   if(!is.null(selectTable(object,type=type)) && any(!is.na(selectTable(object,type=type))) && overwrite==FALSE){
                     stop("allocTable[MRIaggr] : Table already existing in \'object\' \n",                          
                          "set \'overwrite\' to TRUE to replace it \n")
                   }
                   
                   if(type=="lesion"){
                     if(trace){saveTable<- object@table_lesion}
                     object@table_lesion <- value
                   }
                   
                   if(type=="reperfusion"){
                     if(trace){saveTable<- object@table_reperfusion}
                     object@table_reperfusion <- value
                   }                 
                   
                   if(type=="hypoperfusion"){
                     if(trace){saveTable<- object@table_hypoperfusion}
                     object@table_hypoperfusion <- value
                   }
                   
                   validObject(object)
                   if(trace){
                     cat(paste("allocTable[MRIaggr] : @table_",type," has been ",sep=""))
                     if(sum(!is.na(saveTable))==0){cat(" allocated \n")}else{cat(" updated \n ")}
                   }
                   return(object)
                 }              
)

setReplaceMethod(f ="supprContrast",
                 signature ="MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                 definition = function(object,trace=TRUE,value)
                 {
                   nom <- initParameter(object=object,param=value,test=TRUE,init=TRUE,accept.coords=FALSE,accept.index=FALSE,
                                        arg_name="value",method="supprContrast")  
                   
                   value.default <- which( (names(selectDefault_value(object)) %in% nom)*(selectDefault_value(object)!="mask") == 1)
                   value <- which( names(object@contrast) %in% nom )
                   
                   object@contrast <- object@contrast[,-value,drop=FALSE]
                   if(length(value.default>0)){
                     object@default_value <- object@default_value[,-value.default,drop=FALSE]
                   }
                   validObject(object)
                   
                   if(trace){
                     cat("supprContrast[MRIaggr] : ")
                     if(length(nom)==1){
                       cat("Cartography \"",nom,"\" \n",
                         "                         has been removed \n",sep="")
                     }else{
                       cat("Cartographies \"",paste(nom,collapse="\" \""),"\" \n",
                         "                         have been removed \n",sep="")                      
                     }                  
                   }
                   
                   return(object)
                 }
)

setReplaceMethod(f ="supprDescStats",
                 signature ="MRIaggr", # penser a modifier les generic functions car c est la ou sont les vrais arguments par defaut
                 definition = function(object,trace=TRUE,value)
                 {
                   nom_descStats <- selectParameter(object,"ls_descStats")
                   
                   if(length(nom_descStats)==0 || any(value %in% nom_descStats==FALSE)){
                     stop("supprDescStats[MRIaggr] : wrong specification of \'value\' \n",
                          "available element in ls_descStats : ",paste(nom_descStats,collapse=" ")," \n",
                          "unavailable requested elements : ",paste(value[value %in% nom_descStats==FALSE],collapse=" "),"\n")       
                   }
                   
                   for(iter_l in value){
                     object@ls_descStats[[iter_l]] <- NULL
                   }
                   
                   validObject(object)
                   
                   if(trace){
                     cat("supprDescStats[MRIaggr] : ")
                     if(length(value)==1){
                       cat("Element \"",value,"\" \n", 
                         "                          has been removed \n",sep="")
                     }else{
                       cat("Elements \"",paste(value,collapse="\" \""),"\" \n",
                         "                          have been \n",sep="")                      
                     }                     
                   }
                   
                   return(object)
                 }
)

#####  D) Methods #############################################################


#### calc. ####

setMethod(f ="calcBrainMask",
          signature = "MRIaggr",
          definition = function(object,param,type="kmeans",
                                th.breaks=100,th.smoothing=TRUE,th.select_optima=1,th.upper=TRUE,plot=TRUE,
                                kmeans.n_groups=2:4,kmeans.Neighborhood=3,
                                skull.param=NULL,skull.n_groups=3,
                                window=FALSE,filename=paste("calcBrainMask",type,object@identifier,sep="_"),
                                width=1000,height=700,path=NULL,unit="px",res=NA,
                                trace=TRUE,update.object=FALSE,overwrite=FALSE){
            
            initParameter(object=object,param=param,test=TRUE,init=FALSE,
                          accept.coords=FALSE,accept.index=FALSE,accept.mask=FALSE,method="calcBrainMask")
            carto <- selectContrast(object,param=param,format="data.frame")
            coords <- selectCoords(object)
            
            p <- length(param)
            if(type %in% c("kmeans","threshold")==FALSE){
              stop("calcBrainMask[MRIaggr] : wrong specification of  \'type\' \n",
                   "valid types : \"kmeans\" \"threshold\" \n",
                   "requested type : ",type,"\n")
            }
            
            if(plot==TRUE){
              scale <- initWindow(window=window,filename=filename,path=path,width=width,height=height,unit=unit,res=res,
                                  n.plot=1,mfrow=c(1,1),xlim=NULL,ylim=NULL,method="calcBrainMask[MRIaggr]")$scale
            }
            
            if(type == "threshold"){
              
              #### test ####
              if(length(param)!=1){
                stop("calcBrainMask[MRIaggr] : wrong specification of  \'param\' \n",
                     "must have length 1 \n",
                     "length(param) : ",length(param),"\n")
              }
              
              #### initialization ####
              res <- list()
              
              if(length(th.breaks)==1){
                breaks <- seq(min(carto[,param]),max(carto[,param]),length.out=th.breaks)  
                
                #                 test.breaks <- sapply(breaks,function(x){sum(carto[,param]>x)})
                #                 seqY <- seq(test.breaks[1],0,length.out=th.breaks)
                #                 breaks <- unique(sapply(seqY,function(x){breaks[which.min(abs(test.breaks-x))]}))
                
                breaks <- c(seq(breaks[1],breaks[2],length.out=12),tail(breaks,length(breaks)-2))
                breaks <- c(head(breaks,length(breaks)-2),seq(breaks[length(breaks)-1],breaks[length(breaks)],length.out=12))
                
                unvalid.breaks <- c(1:8, (length(breaks)-7):length(breaks))
              }else{
                breaks <- th.breaks              
                unvalid.breaks <- NULL
              }
              n.breaks <-  length(breaks)
              
              res$analysis <- matrix(NA,nrow=n.breaks,ncol=5)
              colnames(res$analysis) <- c("threshold","Nb","dNb","dNb.filtered","optima")
              
              ####
              
              res$analysis[,"threshold"] <- breaks
              res$analysis[,"Nb"] <- sapply(res$analysis[,"threshold"],
                                            function(x){sum(carto[,param]>x)}
              )
              
              res$analysis[,"dNb"] <- c(NA,res$analysis[-n.breaks,"Nb"]-res$analysis[-1,"Nb"])
              
              if(th.smoothing>0 && identical(th.smoothing,TRUE)){
                res$analysis[,"dNb.filtered"] <- res$analysis[,"dNb"]                
              }
              if(th.smoothing>0 && !identical(th.smoothing,TRUE)){
                res$analysis[,"dNb.filtered"] <- as.numeric(filter(res$analysis[,"dNb"], dbinom(0:th.smoothing,th.smoothing,0.5), sides=2))                
              }
              
              #### smoothing ####
              cumul <- 0
              bandwidth <- 2
              
              res$analysis[,"dNb"]
              
              
              while(length(cumul)==1 && bandwidth<=15){
                
                test.post <- c(NA,res$analysis[-n.breaks,"dNb.filtered"]>res$analysis[-1,"dNb.filtered"])
                unvalid.post <- c(unvalid.breaks,which(is.na(test.post)),which(is.na(test.post))-1,length(test.post))
                res$analysis[,"optima"] <- 0
                
                # optima
                for(iter_th in setdiff(1:length(test.post),unvalid.post)){
                  if(test.post[iter_th]==test.post[iter_th+1]){
                    cumul[1] <- cumul[1] + 1
                  }else{
                    res$analysis[iter_th,"optima"] <- 1
                    cumul <- c(0,cumul)
                  }
                }
                
                if(cumul[1]<=1){cumul <- cumul[-1]} # first may be truncated
                if(cumul[length(cumul)]<=1){cumul <- cumul[-length(cumul)]} # last may be truncated
                
                # smoothing
                if(identical(th.smoothing,TRUE) && any(cumul<=1)){     
                  bandwidth <- bandwidth+1
                  res$analysis[,"dNb.filtered"] <- as.numeric(filter(res$analysis[,"dNb"], dbinom(0:bandwidth,bandwidth,0.5), sides=2))                               
                  cumul <- 0
                }               
                
              }
              
              #### optima ####
              optima <- which(res$analysis[,"optima"]==1)
              if(length(optima) < th.select_optima){
                stop("calcBrainMask[MRIaggr] : wrong specification of \'th.select_optima\' \n",
                     "only ",length(optima)," optima are available \n",
                     "(",paste(round(res$analysis[optima,"threshold"],2),collapse=" "),") \n",
                     "requested \'th.select_optima\' : ",th.select_optima,"\n")
              } 
              
              res$th_opt <- data.frame(matrix(NA,ncol=length(optima),nrow=2))
              names(res$th_opt) <- 1:length(optima)
              rownames(res$th_opt) <- c("Th","derivative")
              
              res$th_opt["Th",] <- res$analysis[optima,"threshold"]
              res$th_opt["derivative",] <- res$analysis[optima,"dNb.filtered"]
              
              if(trace==TRUE){
                
                if(th.smoothing>0){
                  cat("Derivative has been smoothed with a Gaussian kernel of width ",bandwidth," breaks \n",sep="")
                }
                
                traceSeuil <- as.numeric(res$th_opt["Th",])
                traceDerivative <- as.numeric(res$th_opt["derivative",])
                traceUnit <- nchar(round(traceSeuil))
                
                if(trace==TRUE){cat("Threshold : derivative (selected)\n")}
                for(iter_optima in 1:length(optima)){   
                  
                  diff <- 6-nchar(round(traceSeuil[iter_optima],digit=5-traceUnit[iter_optima]))
                  
                  cat(rep(" ",max(traceUnit)-traceUnit[iter_optima]),
                      round(traceSeuil[iter_optima],digit=5-traceUnit[iter_optima]),rep(" ",diff+traceUnit[iter_optima])," : ",traceDerivative[iter_optima],sep="")
                  if(iter_optima==th.select_optima){cat(rep(" ",max(0,10-nchar(as.numeric(traceDerivative[iter_optima]))) )," (*)",sep="")}
                  cat("\n")
                }
                
              }
              
              #### select observations ####
              
              if(th.upper){
                res$best_group <- selectContrast(object,param=param,format="vector") > as.numeric(res$th_opt["Th",th.select_optima])
              }else{
                res$best_group <-  selectContrast(object,param=param,format="vector") < as.numeric(res$th_opt["Th",th.select_optima])
              }
              
              res$mask_name <- paste(param,"threshold",th.select_optima,sep=".")
              
              #### plot ####
              
              if(plot==TRUE){
                
                initDisplayWindow(window=window,filename=filename,path=path,width=width,height=height,scale=scale,res=res,
                                  mfrow=c(p,2),bg=NULL,pty=NULL,mar=rep(3,4),mgp=c(2,0.5,0))
                
                plot(res$analysis[,"threshold"],res$analysis[,"Nb"],xlab=param,ylab="Nb",main="Number of voxels",type="o")
                points(breaks[optima[-th.select_optima]],res$analysis[optima[-th.select_optima],"Nb"],col=rainbow(length(optima))[-th.select_optima],pch=15)
                points(breaks[optima[th.select_optima]],res$analysis[optima[th.select_optima],"Nb"],col=rainbow(length(optima))[th.select_optima],pch=8,cex=2)
                abline(v=breaks[optima[th.select_optima]],col=rainbow(length(optima))[th.select_optima])
                legend("topright",legend=1:length(optima),col=rainbow(length(optima)),bty="n",pch=20)
                
                plot(res$analysis[,"threshold"],res$analysis[,"dNb.filtered"],xlab=param,ylab=if(bandwidth>0){paste("dNb.filtered - width=",bandwidth,sep="")}else{"dNb"},main="Derivative",type="o")
                points(breaks[optima[-th.select_optima]],res$analysis[optima[-th.select_optima],"dNb.filtered"],col=rainbow(length(optima))[-th.select_optima],pch=15)
                points(breaks[optima[th.select_optima]],res$analysis[optima[th.select_optima],"dNb.filtered"],col=rainbow(length(optima))[th.select_optima],pch=8,cex=2)
                abline(v=breaks[optima[th.select_optima]],col=rainbow(length(optima))[th.select_optima])
                
                if(!is.null(window) && window %in% c("eps","svg","png","pdf")){
                  dev.off()
                }
                
              }
              
            }
            
            
            if(type == "kmeans"){
              
              #### combinaisons possibles
              names_combin <- NULL
              for(iter_group in kmeans.n_groups){
                names_combin <- c(names_combin,
                                  paste("G",iter_group,".=.",2:iter_group,sep=""),
                                  if(iter_group==3){"G3.>=.2"},
                                  if(iter_group>3){paste("G",iter_group,".>=.",seq(3,iter_group-1),sep="")},
                                  if(iter_group>3){paste("G",iter_group,".<=.",seq(3,iter_group-1),sep="")}
                )
              }
              
              #### stockage des resultats 
              res <- list()
              res$kmeans <- list()
              res$potential <- data.frame(matrix(NA,ncol=2,
                                                 nrow=length(names_combin)))
              names(res$potential) <- c("nb_groups","V") 
              rownames(res$potential) <- names_combin
              res$best_V <- -1
              
              #### potentiel
              if(is.null(kmeans.Neighborhood)){
                stop("calcBrainMask[MRIaggr] : \'kmeans.Neighborhood\' must be defined \n",
                     "in order to choose the right kmeans partition \n")                
              }
              
              for(iter_group in 1:length(kmeans.n_groups)){
                
                if(trace){cat(iter_group," brain groups : ",sep="")}
                res$kmeans[[iter_group]] <- kmeans(carto[,param],centers=kmeans.n_groups[iter_group])
                
                conversion <- rank(res$kmeans[[iter_group]]$center)
                res$kmeans[[iter_group]]$cluster <- conversion[res$kmeans[[iter_group]]$cluster]
                res$kmeans[[iter_group]]$size <- table(res$kmeans[[iter_group]]$cluster)
                res$kmeans[[iter_group]]$center <- tapply(carto[,param],res$kmeans[[iter_group]]$cluster,mean)
                
                test.combin <- lapply(strsplit(names_combin,split=".",fixed=TRUE),function(x){x[1]==paste("G",kmeans.n_groups[iter_group],sep="")})
                combin_group <- strsplit(names_combin,split=".",fixed=TRUE)[which(unlist(test.combin))]
                
                for(iter_combin in 1:length(combin_group)){
                  combin_tempo <- combin_group[[iter_combin]]
                  
                  if(combin_tempo[2]=="="){
                    kmeans_logical <- res$kmeans[[iter_group]]$cluster==as.numeric(combin_tempo[3])
                  }
                  if(combin_tempo[2]==">="){
                    kmeans_logical <- res$kmeans[[iter_group]]$cluster>=as.numeric(combin_tempo[3])
                  }
                  if(combin_tempo[2]=="<="){
                    kmeans_logical <- as.logical((res$kmeans[[iter_group]]$cluster!=1)*as.numeric(res$kmeans[[iter_group]]$cluster<=as.numeric(combin_tempo[3])))
                  }
                  
                  if(trace){cat(paste(combin_tempo,collapse=".")," ",sep="")}
                  
                  res$potential[paste(combin_tempo,collapse="."),"nb_groups"] <- kmeans.n_groups[iter_group]
                  
                  V <- calcFilter(df2array(contrast=as.logical(kmeans_logical),coords=coords)$contrast[[1]],
                                  filter=paste("3D_I",kmeans.Neighborhood,sep=""),norm=FALSE)$res
                  
                  res$potential[paste(combin_tempo,collapse="."),"V"] <- mean(V[kmeans_logical==1])
                  
                  
                  #### retain the best
                  if(res$potential[paste(combin_tempo,collapse="."),"V"]>res$best_V){
                    res$best_V <- res$potential[paste(combin_tempo,collapse="."),"V"]
                    res$best_group <- kmeans_logical
                    res$mask_name <- paste(param,"kmeans",paste(combin_tempo,collapse="."),sep=".")
                  }
                  
                }
                
                if(trace){cat("\n")}
              }
              
            }
            
            if(!is.null(skull.param)){
              
              if(trace){cat("skull groups : ")}
              
              initParameter(object=object,param=skull.param,test=TRUE,init=FALSE,
                            accept.coords=FALSE,accept.index=FALSE,accept.mask=FALSE,method="calcBrainMask",
                            arg_name="skull.param",long_name="skull parameter")
              
              carto <- selectContrast(object,param=skull.param,format="data.frame")
              res$potential_skull <- data.frame(matrix(NA,ncol=2,
                                                       nrow=length(skull.n_groups)+1))
              names(res$potential_skull) <- c("nb_groups","V") 
              rownames(res$potential_skull) <- c("row",skull.n_groups)
              res$potential_skull["row",] <- c(NA,res$best_V)
              index_Skull <- NULL
              
              # identification de la boite cranienne
              for(iter_kmeans in 1:length(skull.n_groups)){
                if(trace){cat(skull.n_groups[iter_kmeans]," ",sep="")}  
                # kmeans
                res_kmeans <- kmeans(carto[,skull.param],centers=skull.n_groups[iter_kmeans])
                
                group <- (1:skull.n_groups[iter_kmeans])[-which.max(res_kmeans$size)]
                recouvrement <- sapply(group,
                                       function(x){mean(res$best_group[res_kmeans$cluster==x])}
                )              
                index_Skull_test <- which(res_kmeans$cluster==group[which.min(recouvrement)])
                kmeans_logical <- res$best_group
                kmeans_logical[index_Skull_test] <- FALSE
                
                # evaluation
                V <- calcFilter(df2array(contrast=kmeans_logical,coords=coords)$contrast[[1]],
                                filter=paste("3D_I",kmeans.Neighborhood,sep=""),norm=FALSE)$res
                
                res$potential_skull[iter_kmeans+1,] <- c(skull.n_groups[iter_kmeans],mean(V[kmeans_logical==1]))
                
                if(which.max(res$potential_skull$V)==(iter_kmeans+1)){
                  index_Skull <- index_Skull_test
                }
              }
              if(trace){cat("\n")}
              
              if(trace){cat("intersect mask : ",sum(res$best_group[index_Skull]),"\n",sep="")}
              
              # boite cranienne retiree du masque
              if(length(index_Skull>0)){
                res$best_group[index_Skull] <- FALSE
              }
            }
            
            return(list(res=res,
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite))
          }
)


setMethod(f ="calcContralateral",
          signature ="MRIaggr",
          definition = function(object,param,num=NULL,type="mean",param.ref=NULL,distband=1,lambda=1,
                                trace=TRUE,update.object=FALSE,overwrite=FALSE)
          { if(is.null(num))
          {num <- 1:object@fieldDim$k}
          
          if(any(c("i_hemisphere","j_hemisphere","hemisphere") %in% selectParameter(object) == FALSE))
          {stop("calcContralateral[MRIaggr] : missing elements in \'object@contrast\' \n",
                "missing parameters : \"i_hemisphere\" \"j_hemisphere\" \"hemisphere\" \n",
                "use calcHemisphere method with update.objet = TRUE to compute and allocate these elements \n")}
          
          if(type %in% c("mean","median","1NN_penalised") ==FALSE)
          {stop("calcContralateral[MRIaggr] : wrong specification of \'type\' \n",
                "valid types : \"mean\" \"median\" \"1NN_penalised\" \n",
                "proposed types : ",type,"\n")}
          
          if(type=="1NN_penalised" && is.null(param.ref))
          {stop("calcContralateral[MRIaggr] : argument \'param.ref\' must be specified if \'type\' is \"1NN_penalised\" \n")}
          
          if(type=="1NN_penalised"){
            sd <- sd(selectContrast(object,param=param.ref,num=num,format="vector"),na.rm=TRUE) 
          }else{
            sd <- 0
            param.ref <- param[1]
          }
          
          if(is.null(param))
          {param <- selectParameter(object)}
          
          data <- selectContrast(object,param=unique(c("index",param.ref,param)),num=num,coords=TRUE)
          n.px <- nrow(data)
          p <- length(param)
          
          data_miroir <- data.frame(matrix(NA,nrow=n.px,ncol=length(param)+5))
          names(data_miroir) <- c("index","i_hemisphere","j_hemisphere","k","hemisphere",paste(param,"_contro",sep=""))
          data_miroir$index <- data$index
          data_miroir$k <- data$k
          
          # changement de repere
          data_miroir[,c("i_hemisphere","j_hemisphere","hemisphere")] <- selectContrast(object,param=c("i_hemisphere","j_hemisphere","hemisphere"))
          
          ####### reperage des slices - hemi-left
          index.plot_lesionL <- numeric(0)
          index.plot_controL <- numeric(0)
          
          if(type=="mean"){type_moy <- TRUE}else{type_moy <- FALSE}
          if(type=="median"){type_med <- TRUE}else{type_med <- FALSE}
          if(type=="1NN_penalised"){type_NN <- TRUE}else{type_NN <- FALSE}
          
          if(trace){cat(paste("Left slice : ",sep=""))}
          for(iter_slice in num)
          { 
            if(trace){cat(paste(iter_slice," ",sep=""))}
            # reperage des px de la slice
            index_k <- which(data_miroir$k==iter_slice)
            index_k_lesion <- na.omit(which(data_miroir[index_k,"hemisphere"]=="left"))
            index_k_contro <- na.omit(which(data_miroir[index_k,"hemisphere"]=="right"))
            
            if(length(index_k_lesion)==0 || length(index_k_contro)==0){next}                
            param_C <- unique(c(param,param.ref))
            
            res <- calcContro_cpp(contrast=as.matrix(data[index_k,param_C,drop=FALSE]), 
                                  coords_px=as.matrix(data_miroir[index_k,c("i_hemisphere","j_hemisphere")]), 
                                  index_k=index_k_lesion-1, index_k_contro=index_k_contro-1,
                                  d_lim=distband, lambda=lambda, param_ref=which(names(data[,param_C,drop=FALSE]) == param.ref)-1, var_ref=sd,
                                  type_moy=type_moy, type_med=type_med, type_NN=type_NN, trace=trace)
            data_miroir[index_k[index_k_lesion],paste(param,"_contro",sep="")] <- res$valeur_NormContro[,param_C %in% param]
            index.plot_lesionL <- c(index.plot_lesionL,index_k[index_k_lesion[which(res$index_plot_k)]])
            index.plot_controL <- c(index.plot_controL,index_k[index_k_contro[which(res$index_plot_k_contro)]])                
          }
          
          if(trace){cat("\n")}
          
          
          ######## reperage des slices - lesion right
          index.plot_lesionR <- numeric(0)
          index.plot_controR <- numeric(0)
          
          if(type=="mean"){type_moy <- TRUE}else{type_moy <- FALSE}
          if(type=="median"){type_med <- TRUE}else{type_med <- FALSE}
          if(type=="1NN_penalised"){type_NN <- TRUE}else{type_NN <- FALSE}
          
          if(trace){cat(paste("Right slice : ",sep=""))}
          for(iter_slice in num)
          { 
            if(trace){cat(paste(iter_slice," ",sep=""))}
            # reperage des px de la slice
            index_k <- which(data_miroir$k==iter_slice)
            index_k_lesion <- na.omit(which(data_miroir[index_k,"hemisphere"]=="right"))
            index_k_contro <- na.omit(which(data_miroir[index_k,"hemisphere"]=="left"))
            
            if(length(index_k_lesion)==0 || length(index_k_contro)==0){next}
            
            ##### C
            param_C <- unique(c(param,param.ref))
            
            res <- calcContro_cpp(contrast=as.matrix(data[index_k,param_C,drop=FALSE]), 
                                  coords_px=as.matrix(data_miroir[index_k,c("i_hemisphere","j_hemisphere")]), 
                                  index_k=index_k_lesion-1, index_k_contro=index_k_contro-1,
                                  d_lim=distband, lambda=lambda, param_ref=which(names(data[,param_C,drop=FALSE]) == param.ref)-1, var_ref=sd,
                                  type_moy=type_moy, type_med=type_med, type_NN=type_NN, trace=trace)
            
            data_miroir[index_k[index_k_lesion],paste(param,"_contro",sep="")] <- res$valeur_NormContro[,param_C %in% param]
            index.plot_lesionR <- c(index.plot_lesionR,index_k[index_k_lesion[which(res$index_plot_k)]])
            index.plot_controR <- c(index.plot_controR,index_k[index_k_contro[which(res$index_plot_k_contro)]])                
          }
          if(trace){cat("\n")}
          
          
          return(list(data=cbind(data[,c("i","j")],data_miroir),
                      index_plot=list(index.plot_lesionR=index.plot_lesionR,
                                      index.plot_lesionL=index.plot_lesionL,
                                      index.plot_controR=index.plot_controR,
                                      index.plot_controL=index.plot_controL),
                      trace=trace,
                      update.object=update.object,
                      overwrite=overwrite
          )
          )
          }
)

setMethod(f ="calcDistMask",
          signature ="MRIaggr",
          definition = function(object,mask,name_newparam=paste("dist",mask,sep="_"),
                                spatial_res=c(1,1,1),as.logical=FALSE,Neighborhood="3D_N10",
                                trace=TRUE,update.object=FALSE,overwrite=FALSE){
            
            test.package <- requireNamespace("RANN",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcDistMask[MRIaggr] : this function requires to have installed the RANN package to work \n")
            }
            
            initParameter(object=object,param=mask,test=TRUE,init=FALSE,
                          accept.coords=FALSE,accept.index=FALSE,method="calcDistMask",
                          arg_name="mask",long_name="mask")
            
            p <- length(mask)
            
            if(length(name_newparam)!=p){
              stop("calcDistMask[MRIaggr] : wrong specification of \'name_newparam\' \n",
                   "it must be a numeric vector of length ",p," \n",
                   "length(name_newparam) : ",length(name_newparam),"\n")
            }
            
            if(length(spatial_res)!=3 || is.numeric(spatial_res)==FALSE){
              stop("calcDistMask[MRIaggr] : wrong specification of \'spatial_res\' \n",
                   "it must be a numeric vector of length 3 \n",
                   "proposed 'spatial_res' : ",paste(spatial_res,collapse=" "),"\n")
            }
            
            data <- selectContrast(object,param=mask,coords=TRUE)
            data$i_scaled <- data$i*spatial_res[1]
            data$j_scaled <- data$j*spatial_res[2]
            data$k_scaled<- data$k*spatial_res[3]
            
            res <- data.frame(matrix(NA,nrow=selectN(object),ncol=p))
            names(res) <- name_newparam
            
            for(iter_mask in 1:p){
              if(trace){cat(iter_mask," ",sep="")}
              
              # test logical 
              data[,mask[iter_mask]] <- initMask(object,mask[iter_mask],test=TRUE,init=as.logical,
                                                 arg_name="mask",long_name="mask",method="calcDistMask",format="vector")
              
              index_mask <- which(data[,mask[iter_mask]]==TRUE)
              
              if(length(index_mask)>0){
                
                mask_outline <- pointsOutline(data[index_mask,c("i","j","k")],filter=Neighborhood)
                mask_outline$i <- mask_outline$i*spatial_res[1]
                mask_outline$j <- mask_outline$j*spatial_res[2]
                mask_outline$k <- mask_outline$k*spatial_res[3]
                
                res[index_mask,iter_mask] <- 0
                res[-index_mask,iter_mask] <- RANN::nn2(data=mask_outline,
                                                        query=data[-index_mask,c("i_scaled","j_scaled","k_scaled")],
                                                        k=1)$nn.dists              
              }
              
            }
            if(trace){cat("\n")}
            
            return(list(res=res,
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite))
            
          }
)

setMethod(f = "calcDistTissues",
          signature = "MRIaggr",
          definition = function(object,param,class,num=NULL,hemisphere="both")
          { 
            ### import
            data_class <- selectContrast(object,num=num,param=class,hemisphere=hemisphere,format="data.frame")
            data_class <- apply(data_class,2,as.numeric)
            p <- length(param)
            n.class <- length(class)
            
            test_range <- apply(data_class,2,range)
            if(any(test_range<0) || any(test_range>1)){
              stop("calcDistTissues[MRIaggr] : wrong specification of \'class\' \n",
                   "class parameters must be in [0;1] \n",
                   "proposed parameters range in : ",paste(range(test_range),collapse=" "),"\n",
                   "problematic parameters : ",paste(class[colSums( (test_range<0)+(test_range>1) )>0],collapse=" "),"\n")
            }
            
            data_param <-  selectContrast(object,num=num,param=param,hemisphere=hemisphere,format="data.frame")
            
            ### calcul des moments
            moments <- data.frame(matrix(NA,nrow=n.class,ncol=4*p+1))
            names(moments) <- c(paste(param,"mu",sep="_"),paste(param,"sigma",sep="_"),paste(param,"skewness",sep="_"),paste(param,"kurtosis",sep="_"),"nb.vx")
            rownames(moments) <- class
            
            for(iter_class in 1:n.class){
              moments[class[iter_class],"nb.vx"] <- sum(data_class[,class[iter_class]])
              
              for(iter_p in 1:p){
                
                moments[class[iter_class],paste(param[iter_p],"mu",sep="_")] <- weighted.mean(data_param[,param[iter_p]],w=data_class[,class[iter_class]],na.rm=TRUE)  
                diff_mu <- data_param[,param[iter_p]]-moments[class[iter_class],paste(param[iter_p],"mu",sep="_")]
                
                moments[class[iter_class],paste(param[iter_p],"sigma",sep="_")] <- sqrt(weighted.mean(diff_mu^2,w=data_class[,class[iter_class]],na.rm=TRUE))  
                diff_sigma <- diff_mu/moments[class[iter_class],paste(param[iter_p],"sigma",sep="_")]
                
                moments[class[iter_class],paste(param[iter_p],"skewness",sep="_")] <- weighted.mean(diff_sigma^3,w=data_class[,class[iter_class]],na.rm=TRUE)
                
                moments[class[iter_class],paste(param[iter_p],"kurtosis",sep="_")] <- weighted.mean(diff_sigma^4,w=data_class[,class[iter_class]],na.rm=TRUE)
                
              }
              
            }
            
            return(moments)   
          }
)

setMethod(f ="calcFilter",
          signature ="MRIaggr",
          definition = function(object,param,filter,norm=TRUE,w_contrast=FALSE,na.rm=FALSE,name_newparam=NULL,
                                trace=TRUE,update.object=FALSE,overwrite=FALSE)
          { 
            # preparation
            initParameter(object=object,param=param,test=TRUE,init=FALSE,
                          accept.coords=FALSE,accept.index=FALSE,method="calcFilter")
            p <- length(param)
            
            if(is.null(name_newparam)){
              if(is.character(filter)){
                name_newparam <- paste(param,filter,sep="_")
              }else{
                name_newparam <- paste(param,"filtered",sep="_")
              }
            }
            
            if(length(param)!=length(name_newparam)){
              stop("calcFilter[MRIaggr] : \'param\' and \'name_newparam\' must have the same length \n",
                   "length(param) : ",length(param),"\n ",
                   "length(name_newparam) : ",length(name_newparam),"\n")
            }
            
            carto <- selectContrast(object,param=param,coords=TRUE,format="data.frame")
            carto <- df2array(contrast=carto[,param],
                              range.coords=object@fieldDim,
                              coords=carto[,c("i","j","k")])
            
            res <- data.frame(matrix(NA,ncol=p+3,nrow=selectN(object)))
            names(res) <- c("i","j","k",name_newparam)
            res[,c("i","j","k")] <- selectCoords(object)
            indexData <- res$i+object@fieldDim$i*(res$j-1)+object@fieldDim$i*object@fieldDim$j*(res$k-1)
            
            for(iter_param in 1:p){
              if(trace){cat(param[iter_param]," ",sep="")}
              
              tempo <- calcFilter(carto$contrast[[iter_param]],filter=filter,norm=norm,
                                  w_contrast=w_contrast,na.rm=na.rm)
              Mfilter <- tempo$filter
              index_NNA <- which(!is.na(tempo$res))
              tempo <- array2df(array=tempo$res,
                                name_newparam=name_newparam[iter_param],names_coords=c("i","j","k"))
              
              res[which(indexData %in% index_NNA),name_newparam[iter_param]] <- tempo[,name_newparam[iter_param]] 
            }
            if(trace){cat("\n")}
            
            return(list(res=res,
                        trace=trace,
                        filter=Mfilter,
                        update.object=update.object,
                        overwrite=overwrite))
            
          }
)

setMethod(f ="calcGroupsMask",
          signature ="MRIaggr",
          definition = function(object,mask,as.logical=FALSE,
                                W="ifany",W.range,W.spatial_res=c(1,1,1),
                                trace=TRUE,update.object=FALSE,overwrite=TRUE)
          { 
            #### test package W
            test.package <- requireNamespace("spam",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcGroupsMask[MRIaggr] : this function requires to have installed the spam package to work \n")
            }
            
            test.package <- requireNamespace("Matrix",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcGroupsMask[MRIaggr] : this function requires to have installed the Matrix package to work \n")
            }
            
            #### preliminaries 
            # test parameter
            initParameter(object=object,param=mask,test=TRUE,init=FALSE,
                          accept.coords=FALSE,accept.index=FALSE,method="calcGroupsMask",
                          arg_name="mask",long_name="parameters")
            p <- length(mask)
            
            # coords
            coords <- selectCoords(object)
            if(length(W.spatial_res)!=3 || is.numeric(W.spatial_res)==FALSE){
              stop("calcGroupsMask[MRIaggr] : wrong specification of \'W.spatial_res\' \n",
                   "it must be a numeric vector of length 3 \n",
                   "proposed 'W.spatial_res' : ",paste(W.spatial_res,collapse=" "),"\n")
            }
            
            # data lesion
            carto <- selectContrast(object,param=mask,format="data.frame")
            
            
            carto[,mask] <- initMask(object,mask,test=TRUE,init=as.logical,
                                     arg_name="mask",long_name="mask",method="calcGroupsMask",format="matrix")
            
            if(identical(W,"ifany") && "W_euclidean" %in% selectParameter(object,"ls_descStats")){
              W <- selectDescStats(object,"W_euclidean")
            }else{
              W <- NULL
            }
            
            if(!is.null(W) && "dgCMatrix" %in% class(W) == FALSE){
              stop("calcGroupsMask[MRIaggr] : wrong specification of \'W\' \n",
                   "\'W\' must be of class dgCMatrix \n",
                   "proposed class of \'W\' : ",paste(class(W),collapse=" "),"\n")
            }
            
            if(!is.null(W) && (nrow(W)!=selectN(object) || ncol(W)!=selectN(object))){
              stop("calcGroupsMask[MRIaggr] : wrong specification of \'W\' \n",
                   "required dim of \'W\' : ",selectN(object)," ",selectN(object),"\n",
                   "proposed dim of \'W\' : ",paste(dim(ls_W[[iter_param]]),collapse=" "),"\n")
            }
            
            #### identification des groupes
            if(trace){
              ncharMax <- max(nchar(mask),15)
              cat("mask parameter",rep(" ",ncharMax-14)," : number of observations per spatial group (total) \n",sep="")
            }
            
            res <- list()
            
            for(iter_param in 1:p){
              if(trace){cat(mask[iter_param],rep(" ",ncharMax-nchar(mask[iter_param]))," : ",sep="")}
              
              param_tempo <- mask[iter_param]
              index_N <- which(carto[,param_tempo]==TRUE)                         
              
              if(length(index_N)>1){
                if(is.null(W)){
                  W_lesion <- calcW(coords[index_N,],method="euclidean",range=W.range,
                                    upper=NULL,format="dgCMatrix",row.norm=FALSE,spatial_res=W.spatial_res)
                }else{              
                  W_lesion <- W[index_N,index_N]
                }
                
                res[[iter_param]] <-  calcGroupsW(W_lesion, max_groups=10000)               
                carto[index_N,iter_param] <- res[[iter_param]]$group
                res[[iter_param]]$group <- carto[,iter_param]
                
              }else{
                res[[iter_param]] <- list()
                res[[iter_param]]$group <- rep(0,selectN(object))
                res[[iter_param]]$group_size <- if(length(iter_param)==1){1}else{NA}
                res[[iter_param]]$group_number <- length(iter_param)
                res[[iter_param]]$group_max <- if(length(iter_param)==1){1}else{NA}
                
                if(length(index_N)==1){
                  res[[iter_param]]$group[index_N] <- 1  
                }
              }
              
              if(trace){
                text_group <- paste(res[[iter_param]]$group_size,collapse=" ")
                cat(text_group,rep(" ",max(0,40-nchar(text_group)))," (",length(index_N),") \n",sep="")
              } 
              
            }
            names(res) <- mask
            
            return(list(res=res,
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite)
            )
          }
)

setMethod(f ="calcHemisphere",
          signature ="MRIaggr",
          definition = function(object,param,num=NULL,p=1,subset=NULL,penalty="symmetry",mask=NULL,as.logical=FALSE,n.points=100,
                                gridSearch=TRUE,i_test=seq(-20,20,by=5),angle_test=seq(-30,30,by=5),unit_angle="degree",
                                NelderMead=TRUE,maxit=100,reltol=0.001,
                                plot=TRUE,window=FALSE,filename=paste(object@identifier,"_calcHemisphere",sep=""),
                                width=1000,height=700,path=NULL,unit="px",res=NA,
                                trace=TRUE,update.object=FALSE,overwrite=FALSE)
          { 
            #             gridSearch <- TRUE
            #             NelderMead <- FALSE            
            
            #### preliminary tests ####
            if(gridSearch==FALSE && NelderMead==FALSE){
              stop("calcHemisphere[MRIaggr] : arguments gridSearch and NelderMead should not be simultaneously FALSE \n")
            }
            
            if(!is.numeric(p) || p<=0)
            {stop("calcHemisphere[MRIaggr] : \'p\' must be a positive number \n",
                  "proposed \'p\' : ",p,"\n")
            }
            
            if(length(param)!=1)
            {stop("calcHemisphere[MRIaggr] : \'param\' must have lenght 1 \n",
                  "proposed \'param\' : ",paste(param,collapse=" "),"\n")
            }
            
            if(penalty %in% c("symmetry","asymmetry") == FALSE)
            {stop("calcHemisphere[MRIaggr] : incorrect specification of \'penalty\' \n",
                  "valid values : \"symmetry\" \"asymmetry\" \n",
                  "proposed value : ",penalty,"\n")
            }
            
            if(unit_angle %in% c("radian","degree")==FALSE)
            {stop("calcHemisphere[MRIaggr] : incorrect specification of \'unit_angle\' \n",
                  "valid units : \"radian\" \"degree\" \n",
                  "proposed unit : ",unit_angle,"\n")
            }
            
            if(!is.null(mask)){
              initParameter(object=object,param=mask,test=TRUE,init=FALSE,
                            accept.coords=FALSE,accept.index=FALSE,method="calcHemisphere",
                            arg_name="mask",long_name="mask")
              
              mask <- selectContrast(object,param=mask,format="matrix")
              if(as.logical==TRUE){mask <- apply(mask,2,as.logical)}
              
              if(!is.logical(mask)){
                stop("calcHemisphere[MRIaggr] : type of \'mask\' is not logical \n",
                     "proposed type : ",paste(is(mask),collapse=" "),"\n",
                     "to force the conversion to logical set \'as.logical\'= TRUE \n")
              }
              mask <- rowSums(mask)>0
            }
            
            if(plot==TRUE){
              scale <- initWindow(window=window,filename=filename,path=path,width=width,height=height,unit=unit,res=res,
                                  n.plot=1,mfrow=c(1,1),xlim=NULL,ylim=NULL,method="calcHemisphere[MRIaggr]")$scale
            }

            #### initialisation ####
            
            #### data
            num <- initNum(object,num=num,method="calcHemisphere")
            
            data <- selectContrast(object,param=param,num=num,subset=subset,coord=TRUE)
            data <- data[is.na(data[,param])==FALSE,]
            n <- nrow(data)
            n.num <- length(num)
            ls.indexK <- lapply(num,function(x){which(data$k==x)-1})
            
            sd_data <- sd(data[,param],na.rm=TRUE)
            sdp_data <- sd_data^p
            
            if(n<2)
            {stop("calcHemisphere[MRIaggr] : there is not enougth data to distinguish 2 hemispheres \n",
                  "nrow(data) : ",n,"\n")
            }
            
            #### seed
            deg2rad <- 2*pi/360
            rad2deg <- 360/2*pi
            if(penalty=="symmetry"){penaltyNA <- 3}else{penaltyNA <- 1} # penalised data with no controlateral correspondant
              
            i_median <- median(unique(data$i))
            j_median <- median(unique(data$j))
            angle_median <- 0
            
            df.optimum <- data.frame(position_i=i_median,
                                     position_j=j_median,
                                     angle_rad=angle_median,
                                     asymetrie=NA)  
                                    
            #### grid search ####
            if(gridSearch==TRUE){
              
              # initialisation
              if(unit_angle=="degree"){
                angle_test <- deg2rad*angle_test
              }
                                          
              grid <- expand.grid(i=i_median+i_test,j=df.optimum$position_j,angle=angle_median+angle_test)
              grid$rank_i <- as.numeric(as.factor(rank(grid$i-i_median)))
              grid$rank_angle <- as.numeric(as.factor(rank(grid$angle)))
              grid$rank <- grid$rank_i+grid$rank_angle-2
              
              n.i_test <- length(i_test)
              n.angle_test <- length(angle_test)
              n.grid <- nrow(grid)
              
              grid$penalty <- 0
              grid$nb <- 0
              grid$moy <- 0
              
              if(penalty=="symmetry"){optimum <- -Inf}else{optimum <- +Inf}  
              order_optimum <- +Inf
              
              if(trace==TRUE){
                index_trace <- round(c(1,seq(1,n.grid,length.out=10)[c(-1,-10)],n.grid))
              }
              
              #### computation
              if(trace==T){
                cat("Grid Search : ",n.grid," parametrisations \n",
                    "i     : ",paste(i_test,collapse=" ")," (in voxels) \n",
                   # "j     : ",paste(df.optimum$position_j,collapse=" "),"\n",
                    "angle : ",paste(round(360*angle_test/(2*pi),1),collapse=" ")," (in degrees) \n",sep="")
              }
              
              for(iter_grid in 1:n.grid){ 
                
                if(trace==T && iter_grid %in% index_trace){cat("*")}
          
                res_cpp <- calcHemi_cpp(coordsI=data$i, coordsJ=data$j, ls_indexK=ls.indexK,  n_num=n.num, value=data[,param], n=n,
                                        i_pos=grid[iter_grid,"i"], j_pos=grid[iter_grid,"j"], angle_pos=grid[iter_grid,"angle"],
                                        penaltyNA=penaltyNA, sd_data=sd_data, p=p, symetrie=(penalty=="symmetry"))
#                 cat(res_cpp$numberAssociated," ",res_cpp$pcNA," ",res_cpp$criteria/res_cpp$numberAssociated," ",res_cpp$compromise,"\n")
                grid[iter_grid,"penalty"] <- res_cpp$criteria
                grid[iter_grid,"nb"] <- res_cpp$numberAssociated
                grid[iter_grid,"moy"] <- res_cpp$compromise
                
                if(penalty=="symmetry"){
                  if(res_cpp$compromise > optimum || (res_cpp$compromise == optimum)&&(grid$rank[iter_grid]<order_optimum) ){
                    optimum <- res_cpp$compromise
                    order_optimum <- grid$rank[iter_grid]
                  }
                }else{
                  if(res_cpp$compromise < optimum || (res_cpp$compromise == optimum)&&(grid$rank[iter_grid]<order_optimum) ){
                    optimum <- res_cpp$compromise
                    order_optimum <- grid$rank[iter_grid]
                  }
                }
              }    
              
              index_optimum <- which.min(abs(grid$moy-optimum))
              
              df.optimum$position_i <- grid[index_optimum,"i"]
              df.optimum$angle_rad <- grid[index_optimum,"angle"]
              df.optimum$asymetrie <- optimum    
              
              test.bordI <- df.optimum$position_i %in% c(i_test[1],tail(i_test,1))
              test.bordAngle <- df.optimum$angle_rad %in% c(angle_test[1],tail(angle_test,1))
              
              if(test.bordI || test.bordAngle){
                gridSearch.cv <- FALSE
                if(trace){cat("\n")}
              }else{
                gridSearch.cv <- TRUE
                if(trace){cat(" cv \n")}
              }
              
            }else{
              gridSearch.cv <- NA
            }
            
            #### optimization ####

            optim.control <- list(fnscale=-1,maxit=maxit,reltol=reltol,trace=if(trace>1){trace}else{0})
            
            if(NelderMead==TRUE){
              
              if(trace==TRUE){cat("Nelder-Mead with optim \n")}
              dim_carto <- selectFieldDim(object)
        
              res <- optim(par=list(df.optimum$position_i,df.optimum$angle_rad),
                           fn=function(value){
                             
                             if(abs(value[2])>pi){if(penalty=="symmetry"){return(-Inf)}else{return(Inf)}}
                             
                             calcHemi_cpp(coordsI=data$i, coordsJ=data$j, ls_indexK=ls.indexK,  n_num=n.num, value=data[,param], n=n,
                                          i_pos=value[1], j_pos=df.optimum$position_j, angle_pos=value[2],
                                          penaltyNA=penaltyNA, sd_data=sd_data, p=p, symetrie=(penalty=="symmetry"))$compromise                             
                           },
                           method="Nelder-Mead",control=optim.control
              )
        
              df.optimum$position_i <- res$par[1]
              df.optimum$angle_rad <- res$par[2]
              df.optimum$asymetrie <- res$value         
            }
       

            #### midplane ####
            j <- seq(min(data$j),max(data$j),length.out=n.points)
            midplane <-  data.frame(i=df.optimum$position_i + sin(df.optimum$angle_rad)*(j-df.optimum$position_j),j=j )
            
            #### hemisphere caracterization ####      
            data_miroir <- data.frame(matrix(NA,nrow=selectN(object),ncol=8))
            names(data_miroir) <- c("index","i","j","k","i_hemisphere","j_hemisphere","hemisphere")
            data_miroir[,c("i","j","k")] <- selectCoords(object)
            data_miroir$index <- selectContrast(object,"index",format="vector")
            if(!is.null(mask)){data_miroir$mask <- mask}
            
            # changement de repere
            data_miroir$i_hemisphere <- cos(df.optimum$angle_rad)*(data_miroir$i-df.optimum$position_i) - sin(df.optimum$angle_rad)*(data_miroir$j-df.optimum$position_j)
            data_miroir$j_hemisphere <- sin(df.optimum$angle_rad)*(data_miroir$i-df.optimum$position_i) + cos(df.optimum$angle_rad)*(data_miroir$j-df.optimum$position_j)
            
            data_miroir$hemisphere <- "undefined"
            data_miroir$hemisphere[data_miroir$i_hemisphere>0] <- "left"
            data_miroir$hemisphere[data_miroir$i_hemisphere<0] <- "right"
            
            hemispheres <- data.frame(left="defined",right="defined",stringsAsFactors=FALSE)
            if(!is.null(mask)){              
              
              index_lesion <- which(data_miroir$mask==TRUE)
              table_lesion <- table(data_miroir$hemisphere[index_lesion])
              
              if("left" %in% names(table_lesion)){
                hemispheres$left <- "lesion"
              }else{
                hemispheres$left <- "contralateral"
              }
              
              if("right" %in% names(table_lesion)){
                hemispheres$right <- "lesion"
              }else{
                hemispheres$right <- "contralateral"
              }
            }
           
            #### display ####
            if(plot==TRUE && gridSearch==TRUE){
              
              initDisplayWindow(window=window,filename=filename,path=path,width=width,height=height,scale=scale,res=res,
                                mfrow=c(1,2),bg=NULL,pty=NULL,mar=rep(3,4),mgp=c(2,0.5,0))
              plot.seq_moy <- seq(min(grid$moy),max(c(grid$moy,df.optimum$asymetrie)),length.out=5)
              plot.seq_i <- unique(grid$i)
              plot.seq_angle <- unique(grid$angle)
              
              # optimum selon i
              palette <- rainbow(n.angle_test)
              plot(NA,
                   xlim=range(c(plot.seq_i,df.optimum$position_i)),xlab="i",
                   ylim=c(plot.seq_moy[1],plot.seq_moy[5]),ylab=penalty,
                   axes=FALSE,main="color=angles")
              box()
              
              axis(1,at=plot.seq_i,label=round(plot.seq_i,4))
              axis(2,at=plot.seq_moy,label=round(plot.seq_moy,4))
              
              for(iter_angle in unique(grid$rank_angle)){
                index_angle <- which(grid$rank_angle==iter_angle)
                col <- palette[iter_angle]
                points(grid[index_angle,c("i","moy")],lty=1,col=col,type="l")
              }
              points(df.optimum[,c("position_i","asymetrie")],col="black",cex=2,pch=15)
              
              # optimum selon angle
              palette <- rainbow(n.i_test)
              plot(NA,
                   xlim=range(c(plot.seq_angle,df.optimum$angle_rad)),xlab="angle",
                   ylim=c(plot.seq_moy[1],plot.seq_moy[5]),ylab=penalty,
                   axes=FALSE,main="color=i")
              box()
              axis(1,at=plot.seq_angle,label=round(plot.seq_angle,4))
              axis(2,at=plot.seq_moy,label=round(plot.seq_moy,4))
              
              for(iter_i in unique(grid$rank_i)){
                index_i <- which(grid$rank_i==iter_i)
                col <- palette[iter_i]                                
                points(grid[index_i,c("angle","moy")],lty=1,col=col,type="l")               
              }
              points(df.optimum[,c("angle_rad","asymetrie")],col="black",cex=2,pch=15)
              
              if(!is.null(window) && window %in% c("eps","svg","png","pdf")){dev.off()}
              
            }
            
            #### export ####
            res <- list()
            
            res$data <- data_miroir[,c("i_hemisphere","j_hemisphere","hemisphere")]
            res$hemispheres <- hemispheres
            res$midplane <- midplane
            res$optimum <- df.optimum
            res$grid <- grid
            res$cv <- gridSearch.cv            
            
            return(list(res=res,
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite))
            
          }  
)

setMethod(f ="calcNormalization",
          signature ="MRIaggr",
          definition = function(object,param,mu_type="mean",sigma_type="sd",rm.CSF=FALSE,rm.WM=FALSE,rm.GM=FALSE,
                                trace=TRUE,update.object=FALSE,overwrite=FALSE)
          {
            # controle preliminaire
            if(rm.GM==TRUE && rm.WM==TRUE){
              stop("calcNormalization[MRIaggr] : \'rm.GM\' and \'rm.WM\' cannot be simultaneously TRUE \n",
                   "set at least one of them to FALSE \n")
            }
            
            if(is.character(mu_type)){
              if(mu_type %in% c("min","1Q","median","3Q","max","mean") == FALSE)
              {stop("calcNormalization[MRIaggr] : wrong specification of \'mu_type\' \n",
                    "if \'mu_type\' is a character it must be one of the following :  \"min\", \"1Q\", \"median\", \"3Q\" , \"max\" or \"mean\" \n",
                    "proposed mu_type : ",mu_type,"\n")}
              
              mu_type <- switch(mu_type,
                                "min"=0,
                                "1Q"=0.25,
                                "median"=0.5,
                                "3Q"=0.75,
                                "max"=1)
            }
            
            if(is.numeric(mu_type)){
              
              if(mu_type > 1 || mu_type < 0){
                stop("calcNormalization[MRIaggr] : wrong specification of \'mu_type\' \n",
                     "\'mu_type\' must be between  0 and 1 \n",
                     "proposed \'mu_type\' : ",mu_type,"\n")
              }
              
              mu_norm <- function(x,w){return(quantile(x[w>0.5],probs=mu_type,na.rm=TRUE))}
            }else{
              mu_norm <- function(x,w){weighted.mean(x,w=w,na.rm=TRUE)}
            }
            
            if(sigma_type %in% c("sd","mad") == FALSE){
              stop("calcNormalization[MRIaggr] : wrong specification of \'sd_type\' \n",
                   "valid types :  \"sd\", \"mad\" \n",
                   "proposed type : ",sd_type,"\n")
            }
            
            if(sigma_type=="sd"){
              sigma_norm <- function(x,w,center){sd(x[w>0.5],na.rm=TRUE)}
            }else{
              sigma_norm <- function(x,w,center){mad(x,w=w,center=center,na.rm=TRUE)}
            }
            
            if("hemisphere" %in% selectParameter(object) == FALSE){
              warning("calcNormalization[MRIaggr] : missing \"hemisphere\" parameter in \'object\' \n",
                      "the computations will be incomplete \n",
                      "use the calcHemisphere and calcContralateral function to compute and allocate it \n")
              test.hemi <- FALSE
            }else{
              test.hemi <- TRUE
            }
            
            # extraction des coordonnees
            coords_both <- selectCoords(object,hemisphere="both")
            
            # extraction des parametres
            carto_both <- selectContrast(object,param=param,hemisphere="both",format="matrix")
            p <- length(param)
            
            # extraction des hemispheres
            if(test.hemi){
              index_hemiL <- which(object@contrast$hemisphere=="left")
              index_hemiR <- which(object@contrast$hemisphere=="right")
            }
            
            # suppression des parties majoritairement CSF
            if(rm.CSF==TRUE || rm.WM==TRUE || rm.GM==TRUE){
              
              if(any(c("WM","GM","CSF") %in% selectParameter(object)==FALSE)){
                stop("calcNormalization[MRIaggr] : impossible to remove the CSF \n",
                     c("WM","GM","CSF")[c("WM","GM","CSF") %in% selectParameter(object)==FALSE]," is not available in \'x\' \n",
                     "use the calcTissueType function to obtain compute and allocate these parameters \n")}            
              
              param_tissue <- c("CSF","WM","GM")[c(rm.CSF==FALSE,rm.WM==FALSE,rm.GM==FALSE)]
              w.tissue <- rowSums(selectContrast(object,param=param_tissue,hemisphere="both",format="matrix"))
            }else{
              w.tissue <- rep(1,selectN(object))
            }
            
            #### global
            if(trace){cat("global \n")}
            
            norm_global <- data.frame(matrix(NA,ncol=length(param),nrow=6))
            names(norm_global) <- param
            row.names(norm_global) <- c("mu_both","mu_left","mu_right","sigma_both","sigma_left","sigma_right")
            
            norm_global["mu_both",] <- apply(carto_both,2,function(x){mu_norm(x,w.tissue)})
            if(test.hemi){
              norm_global["mu_left",] <- apply(carto_both[index_hemiL,,drop=FALSE],2,function(x){mu_norm(x,w.tissue[index_hemiL])})
              norm_global["mu_right",] <- apply(carto_both[index_hemiR,,drop=FALSE],2,function(x){mu_norm(x,w.tissue[index_hemiR])})
            }
            norm_global["sigma_both",] <- sapply(1:p,function(x){sigma_norm(x=carto_both[,x],w=w.tissue,center=norm_global["mu_both",x])})
            if(test.hemi){
              norm_global["sigma_left",] <- sapply(1:p,function(x){sigma_norm(x=carto_both[index_hemiL,x],w=w.tissue[index_hemiL],center=norm_global["mu_left",x])}) 
              norm_global["sigma_right",] <- sapply(1:p,function(x){sigma_norm(x=carto_both[index_hemiR,x],w=w.tissue[index_hemiR],center=norm_global["mu_right",x])}) 
            }
            
            #### par coupe
            if(trace){cat("slice : ")}
            
            normMu_slice_both <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normMu_slice_both) <- param
            normMu_slice_left <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normMu_slice_left) <- param
            normMu_slice_right <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normMu_slice_right) <- param
            normSigma_slice_both <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normSigma_slice_both) <- param
            normSigma_slice_left <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normSigma_slice_left) <- param
            normSigma_slice_right <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normSigma_slice_right) <- param
            
            normMu_3slices_both <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normMu_3slices_both) <- param
            normMu_3slices_left <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normMu_3slices_left) <- param
            normMu_3slices_right <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normMu_3slices_right) <- param
            normSigma_3slices_both <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normSigma_3slices_both) <- param
            normSigma_3slices_left <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normSigma_3slices_left) <- param
            normSigma_3slices_right <- data.frame(matrix(NA,ncol=length(param),nrow=object@fieldDim$k))
            names(normSigma_3slices_right) <- param
            
            
            for(iter_k in 1:object@fieldDim$k){
              if(trace){cat(iter_k," ",sep="")}
              
              index_k <- which(coords_both$k==iter_k)
              
              normMu_slice_both[iter_k,] <- apply(carto_both[index_k,,drop=FALSE],2,function(x){mu_norm(x,w.tissue[index_k])})
              normSigma_slice_both[iter_k,] <- sapply(1:p,function(x){sigma_norm(x=carto_both[index_k,x],w=w.tissue[index_k],center=normMu_slice_both[iter_k,x])})               
              
              if(test.hemi){
                normMu_slice_left[iter_k,] <- apply(carto_both[intersect(index_k,index_hemiL),,drop=FALSE],2,function(x){mu_norm(x,w.tissue[intersect(index_k,index_hemiL)])})
                normSigma_slice_left[iter_k,] <- sapply(1:p,function(x){sigma_norm(x=carto_both[intersect(index_k,index_hemiL),x],w=w.tissue[intersect(index_k,index_hemiL)],center=normMu_slice_left[iter_k,x])}) 
                
                normMu_slice_right[iter_k,] <- apply(carto_both[intersect(index_k,index_hemiR),,drop=FALSE],2,function(x){mu_norm(x,w.tissue[intersect(index_k,index_hemiR)])})
                normSigma_slice_right[iter_k,] <- sapply(1:p,function(x){sigma_norm(x=carto_both[intersect(index_k,index_hemiR),x],w=w.tissue[intersect(index_k,index_hemiR)],center=normMu_slice_right[iter_k,x])}) 
              }
              
              index_k <- which(coords_both$k %in% seq(iter_k-1,iter_k+1))
              
              normMu_3slices_both[iter_k,] <- apply(carto_both[index_k,,drop=FALSE],2,function(x){mu_norm(x,w.tissue[index_k])})
              normSigma_3slices_both[iter_k,] <- sapply(1:p,function(x){sigma_norm(x=carto_both[index_k,x],w=w.tissue[index_k],center=normMu_3slices_both[iter_k,x])})               
              
              if(test.hemi){
                normMu_3slices_left[iter_k,] <- apply(carto_both[intersect(index_k,index_hemiL),,drop=FALSE],2,function(x){mu_norm(x,w.tissue[intersect(index_k,index_hemiL)])})
                normSigma_3slices_left[iter_k,] <- sapply(1:p,function(x){sigma_norm(x=carto_both[intersect(index_k,index_hemiL),x],w=w.tissue[intersect(index_k,index_hemiL)],center=normMu_3slices_left[iter_k,x])})
                
                normMu_3slices_right[iter_k,] <- apply(carto_both[intersect(index_k,index_hemiR),,drop=FALSE],2,function(x){mu_norm(x,w.tissue[intersect(index_k,index_hemiR)])})
                normSigma_3slices_right[iter_k,] <- sapply(1:p,function(x){sigma_norm(x=carto_both[intersect(index_k,index_hemiR),x],w=w.tissue[intersect(index_k,index_hemiR)],center=normMu_3slices_right[iter_k,x])}) 
              }
            }
            cat("\n")
            
            res <- list(norm_global=norm_global,
                        normMu_slice_both=normMu_slice_both,
                        normSigma_slice_both=normSigma_slice_both,
                        normMu_slice_left=normMu_slice_left,
                        normSigma_slice_left=normSigma_slice_left,
                        normMu_slice_right=normMu_slice_right,
                        normSigma_slice_right=normSigma_slice_right,
                        normMu_3slices_both=normMu_3slices_both,
                        normSigma_3slices_both=normSigma_3slices_both,
                        normMu_3slices_left=normMu_3slices_left,
                        normSigma_3slices_left=normSigma_3slices_left,
                        normMu_3slices_right=normMu_3slices_right,
                        normSigma_3slices_right=normSigma_3slices_right)
            
            return(list(res=res,
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite)
            )
          }
)

setMethod(f ="calcRegionalContrast",
          signature ="MRIaggr",
          definition = function(object,param,bandwidth,power=2,diagonal=FALSE,                                
                                W="ifany",W.range,W.spatial_res=c(1,1,1),
                                num=NULL,hemisphere="both",
                                trace=TRUE,name_newparam=paste(param,"regional",sep="_"),update.object=FALSE,overwrite=FALSE){
            #### test package W
            test.package <- requireNamespace("spam",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcRegionalContrast[MRIaggr] : this function requires to have installed the spam package to work \n")
            }
            
            test.package <- requireNamespace("Matrix",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcRegionalContrast[MRIaggr] : this function requires to have installed the Matrix package to work \n")
            }
            
            # definition de la matrice de voisinage
            if(identical(W,"ifany") && "W_euclidean" %in% selectParameter(object,"ls_descStats")){
              if(trace){cat(" loading W ... ")}
              
              W <- selectDescStats(object,"W_euclidean",num=num,hemisphere=hemisphere)
              
              # test symetrie
              #               indexP_tempo <- which(W@p>0)[1]-1 # ou which(diff(W@p)>0)
              #               indexI_tempo <- W@i[1]+1
              #               if(abs(W[indexP_tempo,indexI_tempo]-W[indexI_tempo,indexP_tempo])>0){
              if(any(spam::rowSums(W)-spam::colSums(W)>0.001)){
                if(trace==TRUE){
                  cat("\'W\' is not symetrical \n",
                      "it might be because \'W\' is stored as a upper or a lower matrix \n",
                      "\'W\' is replace by W+t(W) \n")
                }
                W <- W+spam::t(W)
              }
              W <- Matrix::drop0(W,is.Csparse=TRUE)      
              
            }else if(is.null(W) || identical(W,"ifany")){
              if(trace){cat(" computing W ... ")}
              
              W <- calcW(selectCoords(object,num=num,hemisphere=hemisphere),spatial_res=W.spatial_res,
                         range=W.range,upper=NULL,format="dgCMatrix",row.norm=FALSE)               
            }
            
            if("dgCMatrix" %in% class(W) == FALSE){
              stop("calcRegionalContrast[MRIaggr] : wrong specification of \'W\' \n",
                   "\'W\' must be of class dgCMatrix \n",
                   "proposed class of \'W\' : ",paste(class(W),collapse=" "),"\n")
            }
            
            # convert distance to weights
            if(diagonal==TRUE){
              diag(W) <- -1
              W@x[W@x>0] <- EDK(W@x[W@x>=0],bandwidth=bandwidth,power=power)  
              W@x[W@x==-1] <- EDK(0,bandwidth=bandwidth,power=power)    
            }else{
              W@x <- EDK(W@x,bandwidth=bandwidth,power=power)                    
            }
            
            # normalization
            Rsum <- spam::rowSums(W)   
            Rsum[is.na(Rsum)] <- -1                
            W <- W/Rsum
            
            if(trace){cat("W ready \n ")}
            
            carto <- selectContrast(object,param=param,num=num,hemisphere=hemisphere,format="matrix")
            p <- length(param)
            
            if( (nrow(W)!=nrow(carto)) || (ncol(W)!=nrow(carto)) ){
              stop("calcRegionalContrast[MRIaggr] : wrong specification of \'W\' \n",
                   "required dim of \'W\' : ",nrow(carto)," ",nrow(carto),"\n",
                   "proposed dim of \'W\' : ",paste(dim(W),collapse=" "),"\n")
            }
            
            if(length(name_newparam)!=p){
              stop("calcRegionalContrast[MRIaggr] : wrong specification of \'name_newparam\' \n",
                   "it must be a numeric vector of length ",p," \n",
                   "length(name_newparam) : ",length(name_newparam),"\n")
            }
            
            # calcul des valeurs regionales
            carto_regional <- data.frame(matrix(NA,nrow=nrow(carto),ncol=p))
            names(carto_regional) <- name_newparam
            
            if(trace){cat("param : ")}
            for(iter_param in 1:p){
              if(trace){cat(param[iter_param]," ",sep="")}
              carto_regional[,iter_param] <- as.matrix(W %*% carto[,iter_param,drop=FALSE])
            }
            if(trace){cat("\n")}
            
            return(list(res=carto_regional,
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite)
            )
            
          }
)

setMethod(f ="calcROCthreshold",
          signature ="MRIaggr",
          definition = function(object,param,mask,as.logical=FALSE,
                                digit=2,plot="ROC_Youden",digit.plot=3,
                                window=FALSE,filename=paste(object@identifier,"calcROCthreshold",plot,sep="_"),
                                width=1000,height=700,path=NULL,unit="px",res=NA,
                                trace=TRUE,update.object=FALSE,overwrite=FALSE){
            
            test.package <- requireNamespace("ROCR",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcROCthreshold[MRIaggr] : this function requires to have installed the ROCR package to work \n")
            }
            
            initParameter(object=object,param=c(mask,param),test=TRUE,init=FALSE,
                          accept.coords=FALSE,accept.index=FALSE,
                          arg_name="mask/param",method="calcROCthreshold")     
            
            if(length(mask)!=length(param)){
              stop("calcROCthreshold[MRIaggr] : \'mask\' and \'param\' must have the same length \n",
                   "length(mask) : ",length(mask),"\n",
                   "length(param) : ",length(param),"\n")
            }
            
            if(plot %in% c(FALSE,"ROC_Youden","ROC_prev","boxplot_Youden","boxplot_prev") == FALSE){
              stop("calcROCthreshold[MRIaggr] : wrong specification of \'plot\'\n",
                   "valid plot : FALSE \"ROC_Youden\" \"ROC_prev\" \"boxplot_Youden\" \"boxplot_prev\" \n",
                   "requested plot : ",plot,"\n")
            }
            
            p <- length(mask)
            
            if(plot!=FALSE){
              
              res.init <- initWindow(window=window,filename=filename,path=path,width=width,height=height,unit=unit,res=res,
                                     n.plot=p,mfrow=NULL,xlim=NULL,ylim=NULL,method="calcROCthreshold[MRIaggr]")
              scale <- res.init$scale
              mfrow <- res.init$mfrow
            }
            
            #### mise en place du JDD ####
            data <- selectContrast(object,param=c(mask,param))
            
            data[,mask] <- initMask(object,mask,test=TRUE,init=as.logical,
                                    arg_name="mask",long_name="mask",method="calcROCthreshold",format="matrix")
            
            if(!is.null(digit)){
              for(iter_param in param){data[,iter_param] <- round(data[,iter_param],digit=digit)}
            }
            
            
            #### etude des thresholds ####
            res.ROC <- data.frame(matrix(NA,ncol=10,nrow=p))
            names(res.ROC) <- c("mask","param","AUC","AUPRC","OptTh_Youden","Se_Youden","Sp_Youden","OptTh_prev","Se_prev","Sp_prev")
            res.ROC[,"mask"] <- mask
            res.ROC[,"param"] <- param
            if(!is.null(plot) && plot %in% c("ROC_Youden","ROC_prev")){
              data_plot <- list()
            }
            
            for(iter_param in 1:p){
              
              if(sum(data[,mask[iter_param]]==TRUE)==0){
                warning("calcROCthreshold[MRIaggr] : only FALSE were found for ",mask[iter_param],"\n")
                next
              }
              
              if(sum(data[,mask[iter_param]]==FALSE)==0){
                warning("calcROCthreshold[MRIaggr] : only TRUE were found for ",mask[iter_param],"\n")
                next
              }
              
              #### calcul du ROC ####
              roc_tempo <- list()
              prediction_tempo <- ROCR::prediction(data[,param[iter_param]], data[,mask[iter_param]])
              performance_tempo <- ROCR::performance(prediction_tempo,x.measure="spec",measure="sens")
              
              roc_tempo$Specificity <- performance_tempo@x.values[[1]]
              roc_tempo$Sensitivity <- performance_tempo@y.values[[1]]
              roc_tempo$Threshold <- performance_tempo@alpha.values[[1]]
              if(!is.null(plot) && plot %in% c("ROC_Youden","ROC_prev")){
                data_plot[[iter_param]] <- data.frame(Specificity=roc_tempo$Specificity,
                                                      Sensitivity=roc_tempo$Sensitivity,
                                                      Threshold=roc_tempo$Threshold)                
              }
              
              performance_tempo <- ROCR::performance(prediction_tempo,x.measure="rec",measure="prec")
              
              roc_tempo$Recall <- performance_tempo@x.values[[1]]
              roc_tempo$Precision <- performance_tempo@y.values[[1]]
              
              #### summary statistics ####
              res.ROC[iter_param,"AUC"] <- ROCR::performance(prediction_tempo,measure="auc")@y.values[[1]]
              res.ROC[iter_param,"AUPRC"] <- calcAUPRC(x=NULL,y=NULL,performance=performance_tempo)
              
              #### seuil optimal : Youden ####
              OptTh <- which.max(roc_tempo$Sensitivity + roc_tempo$Specificity)
              res.ROC$OptTh_Youden[iter_param] <- roc_tempo$Threshold[OptTh]
              res.ROC$Se_Youden[iter_param] <- roc_tempo$Sensitivity[OptTh]
              res.ROC$Sp_Youden[iter_param] <- roc_tempo$Specificity[OptTh]
              
              #### seuil optimal : prevalence ####
              prevalence <- mean(data[,mask[iter_param]])
              OptTh <- which.max(prevalence*roc_tempo$Sensitivity + (1-prevalence)*roc_tempo$Specificity)
              res.ROC$OptTh_prev[iter_param] <- roc_tempo$Threshold[OptTh]
              res.ROC$Se_prev[iter_param] <- roc_tempo$Sensitivity[OptTh]
              res.ROC$Sp_prev[iter_param] <- roc_tempo$Specificity[OptTh]
              
            }
            
            #### display ####     
            if(plot!=FALSE){
              
              if(p>4){
                stop("calcROCthreshold[MRIaggr] : Too many parameters to be ploted \n",
                     "maximum number of parameters : 4 \n ",
                     "length(mask) : ",length(mask),"\n")
              }
              
              initDisplayWindow(window=window,filename=filename,path=path,width=width,height=height,scale=scale,res=res,
                                mfrow=mfrow,bg=NULL,pty=NULL,mar=rep(3,4),mgp=c(1.5,0.5,0))  
              
              for(iter_param in 1:p){             
                
                if(plot %in% c("ROC_Youden","ROC_prev")){
                  prevalence <- mean(data[,mask[iter_param]])
                  
                  plot(1-data_plot[[iter_param]]$Specificity,data_plot[[iter_param]]$Sensitivity,
                       main=paste(mask[iter_param]," ~ ",param[iter_param]," - patient ",object@identifier),
                       ylab="sensitivity",xlab="1-specificity",type="o")
                  points(c(0,1),c(0,1),type="l",col="grey")
                  
                  if(plot=="ROC_prev"){
                    points(1-res.ROC$Sp_prev[iter_param],res.ROC$Se_prev[iter_param],col="red",pch=15,cex=1.5)
                    legend(x=0.2,y=0.6,xjust=0,yjust=1,col="red",pch=15,
                           legend=paste("Th_Youden : ",round(res.ROC$OptTh_prev[iter_param],digit),"\n",
                                        "Se : ",round(res.ROC$Se_prev[iter_param],digit.plot),"\n",
                                        "Sp : ",round(res.ROC$Sp_prev[iter_param],digit.plot),"\n"),
                           bty="n")
                  }
                  
                  if(plot=="ROC_Youden"){
                    points(1-res.ROC$Sp_Youden[iter_param],res.ROC$Se_Youden[iter_param],col="blue",pch=15,cex=1.5)
                    legend(x=0.2,y=0.6,xjust=0,yjust=1,col="blue",pch=15,
                           legend=paste("Th[Y/prev] : ",round(res.ROC$OptTh_Youden[iter_param],digit),"\n",
                                        "Se : ",round(res.ROC$Se_Youden[iter_param],digit.plot),"\n",
                                        "Sp : ",round(res.ROC$Sp_Youden[iter_param],digit.plot),"\n",
                                        "prevalence : ",round(prevalence,digit.plot)),
                           bty="n")
                  }
                }
                
                if(plot %in% c("boxplot_Youden","boxplot_prev")){
                  boxplot(data[,param[iter_param]] ~ data[,mask[iter_param]],
                          ylab=param[iter_param],main=paste(mask[iter_param]," - ",object@identifier))  
                  
                  if(plot=="boxplot_Youden"){
                    abline(h=res.ROC$OptTh_Youden[iter_param],col="red")
                    text(x=1.5,y=1.1*res.ROC$OptTh_Youden[iter_param],col="red",labels="Youden",bty="n")
                  }
                  if(plot=="boxplot_prev"){
                    abline(h=res.ROC$OptTh_prev[iter_param],col="blue")
                    text(x=1.5,y=1.1*res.ROC$OptTh_prev[iter_param],col="blue",labels="prev",bty="n")
                  }
                }
                
              }
              
              if(!is.null(window) && window %in% c("eps","svg","png","pdf")){
                dev.off()
              }
              
            }
            
            return(list(res=res.ROC,
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite))
          }
)

setMethod(f ="calcSmoothMask",
          signature ="MRIaggr",
          definition = function(object,mask="mask",as.logical=FALSE,
                                size_2Dgroup=50,Neighborhood_2D="3D_N8",rm.2Dhole=FALSE,
                                size_3Dgroup="unique",Neighborhood_3D="3D_N10",rm.3Dhole=TRUE,erosion_th=0.75,  
                                Vmask_min=0.25,Vbackground_max=0.75,Neighborhood_V="3D_N10",
                                trace=TRUE,update.object=FALSE,overwrite=FALSE)
          { #### preparation
            if(is.null(size_2Dgroup) || (size_2Dgroup!=FALSE && size_2Dgroup != "unique" && is.numeric(size_2Dgroup)==FALSE))
            {stop("calcSmoothMask[MRIaggr] : wrong specification of \'size_2Dgroup\' \n",
                  "\'size_2Dgroup\' must be FALSE, numeric or be \"unique\" \n",
                  "proposed value  : ",size_2Dgroup,"\n")}
            
            if(is.null(size_3Dgroup) || (size_3Dgroup!=FALSE && size_3Dgroup != "unique" && is.numeric(size_3Dgroup)==FALSE))
            {stop("calcSmoothMask[MRIaggr] : wrong specification of \'size_3Dgroup\' \n",
                  "\'size_3Dgroup\' must be FALSE, numeric or be \"unique\" \n",
                  "proposed value  : ",size_3Dgroup,"\n")}
            
            if(is.null(erosion_th) || (erosion_th!=FALSE &&  is.numeric(erosion_th)==FALSE))
            {stop("calcSmoothMask[MRIaggr] : wrong specification of \'erosion_th\' \n",
                  "\'erosion_th\' must be FALSE or numeric  \n",
                  "proposed value  : ",erosion_th,"\n")}
            
            coords <- selectCoords(object)            
            n <- selectN(object)            
            num <- unique(coords$k)
            n.slices <- length(num)
            
            index_tous <- 1:n            
            
            #### Application du mask            
            index_fond <- NULL
            index_mask <- index_tous
            
            if(is.character(mask)){
              mask <- selectContrast(object,param=mask,format="matrix",coord=FALSE)
            }
            if(as.logical==TRUE){
              mask <- as.logical(mask)
            }
            if(!is.logical(mask)){
              stop("calcSmoothMask[MRIaggr] : \'mask\' is not of type logical \n",
                   "proposed type : ",paste(is(mask),collapse=" "),"\n",
                   "to force the conversion to logical set \'as.logical\'= TRUE \n")
            }    
            
            index_fond <- index_tous[which(mask==0)]
            index_mask.ref <- index_tous[which(mask==1)] 
            index_mask <- index_mask.ref
            
            #### Exclusion des petits groupes 2D
            if( identical(size_2Dgroup,FALSE)==FALSE && length(index_mask)>0)
            { if(trace){cat("rm small2D : ")}       
              group2D <- calcGroupsCoords(coords=coords[index_mask,],
                                          Neighborhood=Neighborhood_2D,trace=trace)
              
              valid_group2D <- numeric()
              for(iter_num in 1:n.slices){
                index_num <- which(group2D$df.group[,"k"]==iter_num)
                group2D_num <- unique(group2D$df.group[index_num,"group"])
                
                if(size_2Dgroup=="unique"){
                  valid_group2D <- c(valid_group2D,
                                     group2D_num[which.max(group2D$group_size[group2D_num])])
                }else{
                  valid_group2D <- c(valid_group2D,
                                     group2D_num[group2D$group_size[group2D_num]>size_2Dgroup])
                }
              }
              index_unvalid <- group2D$df.group$index[group2D$df.group$group %in%  valid_group2D == FALSE]
              
              if(length(index_unvalid>0)){
                index_fond <- sort(union(index_fond,index_mask[index_unvalid]))
                index_mask <- sort(index_mask[-index_unvalid])
              }
              
              if(trace){cat(" ",sum(group2D$group_size[-valid_group2D]), 
                            " vx ; ",length(group2D$group_size[-valid_group2D])," groups \n",sep="")}              
            }
            
            # Rebouchage des trous des trous 2D
            if( identical(rm.2Dhole,FALSE)==FALSE  && length(index_fond)>0)
            {  if(trace){cat("add hole2D : ")}
               
               group2D <- calcGroupsCoords(coords=coords[index_fond,],
                                           Neighborhood=Neighborhood_2D,trace=trace)
               
               valid_group2D <- numeric()
               for(iter_num in 1:n.slices){
                 index_num <- which(group2D$df.group[,"k"]==iter_num)
                 group2D_num <- unique(group2D$df.group[index_num,"group"])
                 
                 valid_group2D <- c(valid_group2D,
                                    group2D_num[which.max(group2D$group_size[group2D_num])])                 
               }
               index_unvalid <- group2D$df.group$index[group2D$df.group$group %in%  valid_group2D == FALSE]
               
               if(length(index_unvalid>0)){
                 index_mask <- sort(union(index_mask,index_fond[index_unvalid]))
                 index_fond <- sort(index_fond[-index_unvalid])
               }
               
               if(trace){cat(" ",sum(group2D$group_size[-valid_group2D]), 
                             " vx ; ",length(group2D$group_size[-valid_group2D])," groups \n",sep="")}
            }  
            
            #### Exclusion des petits groupes 3D
            if( identical(size_3Dgroup,FALSE)==FALSE && length(index_mask)>0)
            { if(trace){cat("rm small3D : ")}
              
              group3D <- calcGroupsCoords(coords=coords[index_mask,],
                                          Neighborhood=Neighborhood_3D,trace=trace)
              
              if(size_3Dgroup=="unique"){
                valid_group3D <- which.max(group3D$group_size)
              }else{
                valid_group3D <- which(group3D$group_size>size_3Dgroup)                          
              }         
              index_unvalid <- group3D$df.group$index[group3D$df.group$group %in% valid_group3D == FALSE]
              
              if(length(index_unvalid>0)){
                index_fond <- sort(union(index_fond,index_mask[index_unvalid]))              
                index_mask <- sort(index_mask[-index_unvalid])              
              }
              
              if(trace){cat(" ",sum(group3D$group_size[-valid_group3D]), 
                            " vx ; ",length(group3D$group_size[-valid_group3D])," groups \n",sep="")}              
            }
            
            #### erosion
            if( identical(erosion_th,FALSE)==FALSE && size_3Dgroup=="unique" && !is.null(size_3Dgroup)  && length(index_mask)>0){
              
              if(trace){cat("erosion : ")}
              
              Amask <- df2array(rep(0,n),
                                coords=coords)$contrast[[1]]              
              
              # identification des observations a eroder
              Amask[index_fond] <- FALSE
              Amask[index_mask] <- TRUE
              Vlocal <- calcFilter(Amask,filter=Neighborhood_V,norm=TRUE)$res
              
              index_erosion <- intersect(which(Vlocal<erosion_th),index_mask)
              
              index_fond <- sort(union(index_fond,index_erosion))
              index_mask <- sort(setdiff(index_mask,index_erosion))
              
              # identification des nouveaux groupes spatiaux
              Amask[index_fond] <- NA
              Amask[index_mask] <- TRUE
              
              group3D <- calcGroupsCoords(array=Amask,
                                          Neighborhood=Neighborhood_3D,trace=trace)
              
              valid_group3D <- which.max(group3D$group_size)
              
              index_unvalid <- group3D$df.group$index[group3D$df.group$group %in%  valid_group3D == FALSE]
              
              if(length(index_unvalid>0)){
                index_fond <- sort(union(index_fond,index_mask[index_unvalid]))
                index_mask <- sort(index_mask[-index_unvalid])              
              }
              
              # gestion des observations erodes
              Amask[index_fond] <- FALSE
              Amask[index_mask] <- TRUE
              Vlocal <- calcFilter(Amask,filter="2D_N8",norm=TRUE)$res
              
              index_fond <- setdiff(index_fond,index_erosion[which(Vlocal[index_erosion]>0)])
              index_fond <- sort(index_fond)              
              index_mask <- union(index_mask,index_erosion[which(Vlocal[index_erosion]>0)])
              index_mask <- sort(index_mask)
              
              if(trace){cat(" ",sum(group3D$group_size[-valid_group3D])-length(intersect(index_fond,index_erosion[which(Vlocal[index_erosion]>0)])), 
                            " vx ; ",length(group3D$group_size[-valid_group3D])," groups \n",sep="")} 
            }
            
            #### Rebouchage des trous des trous 3D
            if( identical(rm.3Dhole,FALSE)==FALSE  && length(index_fond)>0)
            {  if(trace){cat("add hole3D : ")}
               
               group3D <- calcGroupsCoords(coords=coords[index_fond,],
                                           Neighborhood=Neighborhood_3D,trace=trace)
               
               valid_group3D <- which.max(group3D$group_size)                         
               index_unvalid <- group3D$df.group$index[group3D$df.group$group %in%  valid_group3D == FALSE]
               
               if(length(index_unvalid>0)){
                 index_mask <- sort(union(index_mask,index_fond[index_unvalid]))
                 index_fond <- sort(index_fond[-index_unvalid])          
               }
               
               if(trace){cat(" ",sum(group3D$group_size[-valid_group3D]), 
                             " vx ; ",length(group3D$group_size[-valid_group3D])," groups \n",sep="")}
            }  
            
            # lissage du masque
            if( identical(Vmask_min,FALSE)==FALSE  || identical(Vbackground_max,FALSE)==FALSE){
              
              if(trace){cat("smoothing (add/rm) : ")}
              
              Amask <- df2array(rep(0,n),
                                coords=coords)$contrast[[1]]
              iter_max <- 20
              iter <- 1
              n.modif <- 1
              
              while(n.modif>0 && iter < iter_max){
                Amask[index_fond] <- FALSE
                Amask[index_mask] <- TRUE
                n.modif <- 0
                
                Vlocal <- calcFilter(Amask,filter=Neighborhood_V,norm=TRUE)$res
                
                if(!is.null(Vmask_min)){
                  index_rm <-  intersect(index_mask,which(Vlocal<Vmask_min))
                  n.modif <- n.modif + length(index_rm)
                  if(length(index_rm)>0){
                    index_fond <- union(index_fond,index_rm)
                    index_mask <- setdiff(index_mask,index_rm)                
                  }
                }
                
                if(!is.null(Vbackground_max)){
                  index_add <-  intersect(index_fond,which(Vlocal>Vbackground_max))
                  n.modif <- n.modif + length(index_add)
                  if(length(index_add)>0){
                    index_fond <- setdiff(index_fond,index_add)
                    index_mask <- union(index_mask,index_add)                
                  } 
                }
                if(trace){cat("(",iter,") ",length(index_add),"/",length(index_rm),"  ",sep="")}
                
                iter <- iter + 1
              }
              if(trace){cat("\n")}
            }
            
            
            
            # mise en forme pour l export            
            mask <- rep(FALSE,n)
            mask[index_mask] <- TRUE
            
            return(list(res=data.frame(mask=mask,coords),
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite))
          }
)

setMethod(f ="calcTableHypoReperf",
          signature ="MRIaggr",
          definition = function(object,param,time,threshold=1:10,sep="_",norm_mu=FALSE,norm_sigma=FALSE,
                                mask=NULL,as.logical=FALSE,
                                trace=TRUE,param.update="reperf",update.object=FALSE,overwrite=FALSE)
          {
            #### preparation ####
            if(length(time) %in% c(1,2) == FALSE){
              stop("calcTableHypoReperf[MRIaggr] : wrong specification of \'time\'\n",
                   "time must have length 1 (hypoperfusion only) or 2 (hypoperfusion and reperfusion) \n",
                   "length(time) : ",length(time),"\n")
            }
            
            if(any(param.update %in% c("shift","reperf","reperf_pc","deperf","deperf_pc")) == FALSE){
              stop("calcTableHypoReperf[MRIaggr] : wrong specification of \'param.update\'\n",
                   "\'param.update\' must be in : \"shift\" \"reperf\" \"reperf_pc\" \"deperf\" \"deperf_pc\" \n",
                   "proposed \'param.update\' : ",paste(param.update,collapse=" "),"\n")
            }
            
            #### hypoperfusion ####
            
            ## mise en place 
            param_time1 <- paste(param,time[1],sep=sep)
            initParameter(object=object,param=param_time1,test=TRUE,init=FALSE,accept.coords=FALSE,
                          arg_name="param_time1",method="calcTableHypoReperf")
            if(length(time)==2){
              param_time2 <- paste(param,time[2],sep=sep)
              initParameter(object=object,param=param_time2,test=TRUE,init=FALSE,accept.coords=FALSE,
                            arg_name="param_time2",method="calcTableHypoReperf")
            }
            if(!is.null(mask)){
              initParameter(object=object,param=mask,test=TRUE,init=FALSE,accept.coords=FALSE,
                            arg_name="mask",method="calcTableHypoReperf")
            }
            
            n.param <- length(param)
            n.threshold <- length(threshold)
            
            index_mask <- selectContrast(object,param="index",format="vector")
            n.mask <- length(index_mask)
            
            res <- list()
            if(length(time)==1){
              res$volume_hypo <- data.frame(matrix(NA,nrow=n.threshold,ncol=3*n.param+1))
              names(res$volume_hypo ) <- c("threshold",paste("Vhypo.",param_time1,sep=""),
                                           paste("Vmismatch.",param,sep=""),paste("PCmismatch.",param,sep=""))
            }else{
              res$volume_hypo <- data.frame(matrix(NA,nrow=n.threshold,ncol=4*n.param+1))
              names(res$volume_hypo ) <- c("threshold",paste("Vhypo.",param_time1,sep=""),paste("Vhypo.",param_time2,sep=""),
                                           paste("Vmismatch.",param,sep=""),paste("PCmismatch.",param,sep=""))
            }
            row.names(res$volume_hypo) <- c(as.character(threshold))
            res$volume_hypo$threshold <- threshold
            
            ## data hypoperfusion
            dataRaw.time1 <- selectContrast(object,param=param_time1,norm_mu=norm_mu,norm_sigma=norm_sigma)
            
            data.time1 <- apply(dataRaw.time1,2,function(x){
              tempo <- cut(x,breaks=c(-Inf,threshold-1*10^{-12},Inf))
              levels(tempo) <- c(threshold[1]-1,threshold,threshold[n.threshold]+1)
              return(as.numeric(as.character(tempo)))
            })
            dataRaw.time1[dataRaw.time1<0] <- 0
            dataRaw.time1[dataRaw.time1>threshold[n.threshold]] <- threshold[n.threshold]
            
            if(!is.null(mask)){
              data.mask <- initMask(object,mask,test=TRUE,init=as.logical,
                                    arg_name="mask",long_name="mask",method="calcTableHypoReperf")
              
              test.mask <- (data.mask==FALSE)
            }
            
            #### reperfusion ####
            if(length(time)==2){
              ## mise en place
              res$pixel <- data.frame(matrix(NA,nrow=selectN(object),ncol=5*n.param+3))
              res$pixel[,1:3] <- selectCoords(object)
              names(res$pixel) <- c("i","j","k",paste(param,"_shift",sep=""),
                                    paste(param,"_reperf",sep=""),paste(param,"_reperf_pc",sep=""),
                                    paste(param,"_deperf",sep=""),paste(param,"_deperf_pc",sep=""))
              
              res$volume_reperf <- data.frame(matrix(NA,nrow=n.threshold,ncol=12*n.param+1))
              names(res$volume_reperf) <- c("threshold",
                                            paste("Vreperf.",param,sep=""),paste("PCreperf.",param,sep=""),paste("VreperfW.",param,sep=""),paste("PCreperfW.",param,sep=""),
                                            paste("Vdeperf.",param,sep=""),paste("PCdeperf.",param,sep=""),paste("VdeperfW.",param,sep=""),paste("PCdeperfW.",param,sep=""),
                                            paste("Vshift_reperf.",param,sep=""),paste("PCshift_reperf.",param,sep=""),
                                            paste("Vshift_deperf.",param,sep=""),paste("PCshift_deperf.",param,sep=""))
              row.names(res$volume_reperf) <- c(as.character(threshold))
              res$volume_reperf$threshold <- threshold
              
              ## data              
              dataRaw.time2 <- selectContrast(object,param=param_time2,norm_mu=norm_mu,norm_sigma=norm_sigma)
              
              data.time2 <- apply(dataRaw.time2,2,function(x){
                tempo <- cut(x,breaks=c(-Inf,threshold-1*10^{-12},Inf))
                levels(tempo) <- c(threshold[1]-1,threshold,threshold[n.threshold]+1)
                return(as.numeric(as.character(tempo)))
              })
              dataRaw.time2[dataRaw.time2<0] <- 0
              dataRaw.time2[dataRaw.time2>threshold[n.threshold]] <- threshold[n.threshold]
              
              res$pixel[index_mask,paste(param,"_shift",sep="")] <- dataRaw.time2 - dataRaw.time1
            }
            
            #### calcul ####
            
            for(iter_param in param){
              if(trace){cat(iter_param," ",sep="")}
              
              iter_param1 <- paste(iter_param,sep,time[1],sep="")
              index_H0pos <- which(dataRaw.time1[,iter_param1]>0)
              index_H0lim <- which(dataRaw.time1[,iter_param1]<threshold[n.threshold])
              
              if(length(time)==2){
                iter_param2 <- paste(iter_param,sep,time[2],sep="")
                iter_paramShift <- paste(iter_param,"_shift",sep="")
                iter_paramReperf <- paste(iter_param,"_reperf",sep="")
                iter_paramReperfPC <- paste(iter_param,"_reperf_pc",sep="")
                iter_paramDeperf <- paste(iter_param,"_deperf",sep="")
                iter_paramDeperfPC <- paste(iter_param,"_deperf_pc",sep="")
                
                # reperf px              
                index_reperf <- which(res$pixel[index_mask,iter_paramShift]<0)
                res$pixel[index_mask[index_reperf],iter_paramReperf] <- -data.time1[index_reperf,iter_param1]
                res$pixel[index_mask[index_reperf],iter_paramReperfPC] <- res$pixel[index_mask[index_reperf],iter_paramShift]/(data.time1[index_reperf,iter_param1])
                
                # deperf px
                index_deperf <- which(res$pixel[index_mask,iter_paramShift]>0)
                res$pixel[index_mask[index_deperf],iter_paramDeperf] <- data.time1[index_deperf,iter_param1]
                res$pixel[index_mask[index_deperf],iter_paramDeperfPC] <- res$pixel[index_mask[index_deperf],iter_paramShift]/(threshold[n.threshold]-data.time1[index_deperf,iter_param1])
              }
              
              for(iter_threshold in 1:n.threshold){
                test.hypoH0 <- data.time1[,iter_param1]>=threshold[iter_threshold]
                res$volume_hypo[iter_threshold,paste("Vhypo.",iter_param1,sep="")] <- sum(test.hypoH0,na.rm=T)
                
                if(!is.null(mask)){
                  res$volume_hypo[iter_threshold,paste("Vmismatch.",iter_param,sep="")] <- sum( test.hypoH0*test.mask,na.rm=T)
                }
                
                if(length(time)==2){
                  test.hypoH2 <- data.time2[,iter_param2]>=threshold[iter_threshold]
                  res$volume_hypo[iter_threshold,paste("Vhypo.",iter_param2,sep="")] <- sum(test.hypoH2,na.rm=T)
                  
                  test.hyperH0 <- data.time1[,iter_param1]<threshold[iter_threshold]
                  test.hyperH2 <- data.time2[,iter_param2]<threshold[iter_threshold]
                  index_n0Reperf <- intersect(index_H0pos, index_reperf)
                  index_n0Deperf <- intersect(index_H0lim, index_deperf)
                  w.reperf <- rep(0,n.mask)
                  w.reperf[index_n0Reperf] <- -res$pixel[index_mask,iter_paramShift][index_n0Reperf]/dataRaw.time1[index_n0Reperf,iter_param1]
                  w.deperf <- rep(0,n.mask)
                  w.deperf[index_n0Deperf] <- res$pixel[index_mask,iter_paramShift][index_n0Deperf]/(threshold[n.threshold]-dataRaw.time1[index_n0Deperf,iter_param1])
                  
                  # Vreperf et Vdeperf
                  res$volume_reperf[iter_threshold,paste("Vreperf.",iter_param,sep="")] <- sum( test.hypoH0*test.hyperH2, na.rm=T)
                  res$volume_reperf[iter_threshold,paste("VreperfW.",iter_param,sep="")] <- sum( test.hypoH0*test.hyperH2*w.reperf, na.rm=T)
                  
                  #                   if(!is.null(mask)){
                  res$volume_reperf[iter_threshold,paste("Vdeperf.",iter_param,sep="")] <- sum( test.hyperH0*test.hypoH2, na.rm=T)
                  res$volume_reperf[iter_threshold,paste("VdeperfW.",iter_param,sep="")] <- sum( test.hyperH0*test.hypoH2*w.deperf, na.rm=T)
                  #                   }
                  
                  # Vshift
                  res$volume_reperf[iter_threshold,paste("Vshift_reperf.",iter_param,sep="")] <- sum(res$pixel[index_mask,iter_paramShift] <= -threshold[iter_threshold],na.rm=T)
                  res$volume_reperf[iter_threshold,paste("Vshift_deperf.",iter_param,sep="")] <- sum(res$pixel[index_mask,iter_paramShift] >= threshold[iter_threshold],na.rm=T)
                }    
              }
              
              if(!is.null(mask)){
                res$volume_hypo[,paste("PCmismatch.",iter_param,sep="")] <- res$volume_hypo[,paste("Vmismatch.",iter_param,sep="")]/sum(data.mask==TRUE)
              }
              if(length(time)==2){
                res$volume_reperf[,paste("PCreperf.",iter_param,sep="")] <- res$volume_reperf[,paste("Vreperf.",iter_param,sep="")]/res$volume_hypo[,paste("Vhypo.",iter_param1,sep="")]
                res$volume_reperf[,paste("PCdeperf.",iter_param,sep="")] <- res$volume_reperf[,paste("Vdeperf.",iter_param,sep="")]/res$volume_hypo[,paste("Vhypo.",iter_param1,sep="")]
                res$volume_reperf[,paste("PCreperfW.",iter_param,sep="")] <- res$volume_reperf[,paste("VreperfW.",iter_param,sep="")]/res$volume_hypo[,paste("Vhypo.",iter_param1,sep="")]
                res$volume_reperf[,paste("PCdeperfW.",iter_param,sep="")] <- res$volume_reperf[,paste("VdeperfW.",iter_param,sep="")]/res$volume_hypo[,paste("Vhypo.",iter_param1,sep="")]
                res$volume_reperf[,paste("PCshift_reperf.",iter_param,sep="")] <- res$volume_reperf[,paste("Vshift_reperf.",iter_param,sep="")]/res$volume_hypo[,paste("Vhypo.",iter_param1,sep="")]
                res$volume_reperf[,paste("PCshift_deperf.",iter_param,sep="")] <- res$volume_reperf[,paste("Vshift_deperf.",iter_param,sep="")]/res$volume_hypo[,paste("Vhypo.",iter_param1,sep="")]              
              }
            }
            if(trace){cat("\n")}
            
            
            
            #### export ####
            return(list(res=res,
                        trace=trace,
                        update.object=update.object,
                        param.update=param.update,
                        overwrite=overwrite)
            )
          }
)

setMethod(f ="calcTableLesion",
          signature ="MRIaggr",
          definition = function(object,maskN,mask=NULL,as.logical=FALSE,
                                trace=TRUE,update.object=FALSE,overwrite=FALSE)
          {  
            dim <- object@fieldDim
            params <- names(object@contrast)
            
            #### tests preliminaires : maskN
            initParameter(object=object,param=maskN,test=TRUE,init=FALSE,accept.coords=TRUE,
                          arg_name="maskN",long_name="parameters",method="calcTableLesion")          
            n.maskN <- length(maskN)
            
            data <- selectContrast(object,param=maskN,coords="k",format="data.frame")
            
            data[,maskN] <- initMask(object,maskN,test=TRUE,init=as.logical,
                                     arg_name="mask",long_name="mask",method="calcTableLesion",format="matrix")
            
            #### tests preliminaires : mask
            if(length(mask)>0){
              if(is.character(mask)){
                mask <- selectContrast(object,param=mask,format="vector")
              }else{
                mask <- as.vector(mask)
                if(!is.null(mask) && length(mask)!=selectN(object)){
                  stop("calcTableLesion[MRIaggr] : length of \'mask\' incompatible with \'object\' \n",
                       "number of observations in \'object\' : ",selectN(object),"\n",
                       "length(mask) : ",length(mask),"\n")
                }
              }
              if(as.logical==TRUE){mask <- as.logical(mask)}
              if(!is.logical(mask)){
                stop("calcTableLesion[MRIaggr] : type of \'mask\' is not logical \n",
                     "proposed type : ",paste(is(mask),collapse=" "),"\n",
                     "to force the conversion to logical set \'as.logical\'= TRUE \n")
              }
              data <- data.frame(data,mask=mask)
            }
            
            matlevels <- matrix(levels(interaction(maskN, maskN,sep="_outside_")),
                                nrow=n.maskN)
            
            names.res <- c(maskN,
                           matlevels[lower.tri(matlevels)],
                           matlevels[upper.tri(matlevels)]
            )
            
            if(!is.null(mask)){
              names.res <- c("mask",names.res,paste(maskN,"_outside_mask",sep=""))
            }
            
            res <- data.frame(matrix(NA,ncol=length(names.res),nrow=dim$k+1))
            names(res) <- names.res
            rownames(res)[dim$k+1] <- c("total")
            
            
            #### calcul des tables : sans interaction
            for(iter_param in c(if(!is.null(mask)){"mask"},maskN)){
              table_tempo <- tapply(data[,iter_param],data[,"k"],
                                    function(x){table_tempo <- table(x) ; if("TRUE" %in% names(table_tempo)){table_tempo["TRUE"]}else{0}})  
              res[,iter_param] <- c(table_tempo,sum(table_tempo))              
            }
            
            #### calcul des tables : avec interaction
            for(iter_param in c(matlevels[lower.tri(matlevels)],matlevels[upper.tri(matlevels)])){
              param_tempo <- unlist(strsplit(iter_param,split="_outside_"))
              
              table_tempo <- tapply(data[,param_tempo[1]]>data[,param_tempo[2]],data[,"k"],function(x){table(x)["TRUE"]})  
              table_tempo[is.na(table_tempo)] <- 0
              res[,iter_param] <- c(table_tempo,sum(table_tempo))              
            }
            
            #### calcul des tables : avec le mask
            if(!is.null(mask)){
              for(iter_param in maskN){
                param_tempo <- c(iter_param,"mask")                
                table_tempo <- tapply(data[,param_tempo[1]]-data[,param_tempo[2]],data[,"k"],function(x){table(x)["1"]})  
                table_tempo[is.na(table_tempo)] <- 0
                res[,paste(iter_param,"mask",sep="_outside_")] <- c(table_tempo,sum(table_tempo))              
              }
            }
            
            
            #### export
            return(list(res=res,
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite))
            
          }
)

setMethod(f ="calcThresholdMRIaggr",
          signature ="MRIaggr",
          definition = function(object,param,hemisphere="both",rm.CSF=FALSE,threshold=1:10,decreasing=FALSE,
                                GRalgo=FALSE,W="ifany",seed=NULL,as.logical=FALSE,W.range,W.spatial_res=rep(1,3),
                                name_newparam=paste(param,"Th",sep="_"),trace=TRUE,update.object=FALSE,overwrite=FALSE){
            
            #### test ####
            if(length(param)!=length(name_newparam)){
              stop("calcThresholdMRIaggr[MRIaggr] : wrong specification of \'name_newparam\' \n",
                   "must be a character vector of length : ",length(param),"\n",
                   "proposed \'name_newparam\' : ",paste(name_newparam,collapse=" "),"\n")
            }
            
            #### ####
            param_select <- param
            if(hemisphere!="both"){
              param_select <- c(param_select)
            }
            
            if(is.character(seed) && seed %in% selectParameter(object)){
              param_select <- c(param_select,seed)
            }
            
            data <- selectContrast(object,param=c("index",param_select),hemisphere=hemisphere)
            
            if(rm.CSF==TRUE){
              data$CSF <- as.numeric(apply(selectContrast(object,param=c("CSF","GM","WM"),hemisphere=hemisphere),1,which.max)==1)             
              data <- data[data$CSF!=1,]
            }
            
            if(GRalgo==TRUE){
              
              #### test package W
              test.package <- requireNamespace("spam",quietly=TRUE)
              if(test.package==FALSE){
                stop("calcThresholdMRIaggr[MRIaggr] : this function with argument GRalgo=TRUE requires to have installed the spam package to work \n")
              }
              
              test.package <- requireNamespace("Matrix",quietly=TRUE)
              if(test.package==FALSE){
                stop("calcThresholdMRIaggr[MRIaggr] : this function with argument GRalgo=TRUE requires to have installed the Matrix package to work \n")
              }
              
              if(identical(W,"ifany") && "W_euclidean" %in% selectParameter(object,"ls_descStats")){
                W <- selectDescStats(object,"W_euclidean",subset_W=data$index)
                
                if(any(spam::rowSums(W)-spam::colSums(W)>0.001)){
                  if(trace==TRUE){
                    cat("\'W\' is not symetrical \n",
                        "it might be because \'W\' is stored as a upper or a lower matrix \n",
                        "\'W\' is replace by W+t(W) \n")
                  }
                  W <- W+spam::t(W)
                }
              }else{
                if(trace==TRUE){cat("computing W ... \n")}
                coords <- selectCoords(object)[data$index,]                
                W <- calcW(coords,range=W.range,method="euclidean",upper=NULL,format="dgCMatrix",row.norm=TRUE,
                           spatial_res=W.spatial_res)                
              }            
            }
            
            
            if(GRalgo==TRUE && !is.null(seed) && is.character(seed)){              
              data[,seed] <- selectContrast(object,param=seed,format="matrix")[data$index,]                 
            }
            
            res.th <- calcThreshold(contrast=cbind(data,CSF=0,hemisphere=hemisphere),param=param,hemisphere=if(hemisphere=="both"){NULL}else{hemisphere},rm.CSF=rm.CSF,threshold=threshold,decreasing=decreasing,
                                    W=W,GRalgo=GRalgo,as.logical=as.logical,seed=seed,trace=trace)
            #             res.th <- calcThreshold(contrast=data,param=param,hemisphere=NULL,rm.CSF=FALSE,threshold=threshold,decreasing=decreasing,
            #                                     W=W,GRalgo=GRalgo,as.logical=as.logical,seed=seed,trace=trace)
            
            names(res.th) <- name_newparam
            
            #### export
            res <- data.frame(matrix(0,nrow=selectN(object),ncol=ncol(res.th)))
            names(res) <- name_newparam
            res[data$index,] <- res.th
            
            return(list(res=res,
                        trace=trace,
                        name_newparam=name_newparam,
                        update.object=update.object,
                        overwrite=overwrite))
            
          }
)

setMethod(f ="calcTissueType",
          signature ="MRIaggr",
          definition = function(object,param,niter=100,nnei=6,
                                beta=if(sub==TRUE){0.3}else{0.7},sub=TRUE,digit=0,trace=TRUE,
                                name_newparam=c("CSF","GM","WM"),update.object=FALSE,overwrite=FALSE){
            
            test.package <- requireNamespace("mritc",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcTissueType[MRIaggr] : this function requires to have installed the mritc package to work \n")
            }
            
            carto <- selectContrast(object,param=param,format="vector")
            coords <- selectCoords(object)
            
            init  <- mritc::initOtsu(round(carto,digit=digit), m=length(name_newparam)-1)
            
            mask <- df2array(rep(1,selectN(object)),
                             coords,default_value=0)$contrast[[1]]
            
            W <-  mritc::makeMRIspatial(mask, nnei=nnei,sub=sub)
            
            res <- mritc::mritc.bayes(y=carto,
                                      neighbors=W$neighbors,
                                      blocks=W$blocks,
                                      mu=init$mu,
                                      sigma=init$sigma,
                                      sub=sub,
                                      niter=niter,
                                      subvox=W$subvox,
                                      verbose=trace)
            
            order <- order(res$mu,decreasing=FALSE)
            res$prob <- res$prob[,order]
            res$sigma <- res$sigma[order]
            res$mu <- res$mu[order]
            
            return(list(res=res,
                        trace=trace,
                        name_newparam=name_newparam,
                        update.object=update.object,
                        overwrite=overwrite)
            )
          }
)

setMethod(f ="calcW",
          signature ="MRIaggr",
          definition = function(object,range,spatial_res=c(1,1,1),num=NULL,hemisphere="both",subset=NULL,
                                upper=TRUE,format="dgCMatrix",row.norm=FALSE,
                                trace=TRUE,update.object=FALSE,overwrite=FALSE){
            
            #### test package W
            test.package <- requireNamespace("spam",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcW[MRIaggr] : this function requires to have installed the spam package to work \n")
            }
            
            test.package <- requireNamespace("Matrix",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcW[MRIaggr] : this function requires to have installed the Matrix package to work \n")
            }
            
            W <- calcW(object=selectCoords(object,num=num,hemisphere=hemisphere,subset=subset),
                       spatial_res=spatial_res,range=range,upper=upper,format=format,row.norm=row.norm)
            
            W <- Matrix::drop0(W,is.Csparse=TRUE)       
            
            return(list(res=W,
                        trace=trace,
                        update.object=update.object,
                        overwrite=overwrite))
            
          }
)

#### plot. ####
setMethod(f ="boxplotMask",
          signature ="MRIaggr",
          definition = function(object,param,mask,num=NULL,hemisphere="both",norm_mu=FALSE,norm_sigma=FALSE,scale=TRUE,as.logical=FALSE,
                                window=FALSE,ylim=NULL,col=c("white","purple"),main=NULL,mgp=c(2,0.5,0),x.legend="topright",y.legend=NULL,cex.legend=0.8,
                                filename=paste(object@identifier,"boxplotMask",sep="_"),
                                width=1000,height=700,path=NULL,unit="px",res=NA){
            
            ### preparation
            p <- length(param)
            
            if(length(mask)!=1){
              stop("boxplotMask[MRIaggr] : wrong specification of \'mask\'\n",
                   "only one mask can be specified \n",
                   "length(mask) : ",length(mask),"\n")                  
            }
            
            carto <- selectContrast(object,param=param,num=num,hemisphere=hemisphere,norm_mu=norm_mu,norm_sigma=norm_sigma,format="data.frame")
            if(scale==TRUE){carto <- scale(carto)}
            
            carto_mask <- selectContrast(object,param=mask,num=num,hemisphere=hemisphere,format="data.frame")
            
            if(as.logical==TRUE){carto_mask[,1] <- as.logical(carto_mask[,1])}
            if(!is.logical(carto_mask[,1])){
              stop("boxplotMask[MRIaggr] : type of \'mask\' is not logical \n",
                   "proposed type : ",paste(is(carto_mask),collapse=" "),"\n",
                   "to force the conversion to logical set \'as.logical\'= TRUE \n")
            }
            
            carto_ls <- list()
            for(iter_param in 1:p){
              carto_ls[[paste(param[iter_param])]] <- carto[carto_mask[,1]==FALSE,iter_param]
              carto_ls[[paste(param[iter_param],1)]] <- carto[carto_mask[,1]==TRUE,iter_param]
            }
            
            ### test export
            res.init <- initWindow(window=window,filename=filename,path=path,width=width,height=height,unit=unit,res=res,
                                   n.plot=1,mfrow=1,xlim=NULL,ylim=NULL,
                                   method="boxplotMask[MRIaggr]")
            scale <- res.init$scale
            
            #### display                       
            initDisplayWindow(window=window,filename=filename,path=path,width=width,height=height,scale=scale,res=res,
                              mfrow=NULL,bg=NULL,pty=NULL,mar=NULL,mgp=mgp)              
            
            if(is.null(ylim)){
              ylim <- range(carto)
            }
            
            boxplot(carto_ls,ylim=ylim,col=col,main=main,names=rep("",p*2))  
            mtext(text=param,side = 1, line = 1, at = seq(from=1.5,by=2,length.out=p))
            abline(v=seq(from=2.5,by=2,length.out=p-1),lty=2)
            
            legend(x=x.legend,y=y.legend,bty="n",legend=c("other",mask),pch=c(21,21),col=rep("black",2),
                   cex=cex.legend)
            legend(x=x.legend,y=y.legend,bty="n",legend=c("other",mask),pch=c(20,20),col=col,
                   cex=cex.legend)
            
            if(!is.null(window) && window %in% c("eps","svg","png","pdf")){
              dev.off()
            }
            
          }
)

setMethod(f ="heatmapMRIaggr",
          signature ="MRIaggr",
          definition = function(object,param,num=NULL,hemisphere="both",
                                scale=TRUE,method="pearson",points.values=TRUE,type="image",digit=3,breaks=NULL,
                                window=FALSE,col=cm.colors(256),main=NULL,mgp=c(2,0.5,0),mar=c(4,4,1,6),las=1,cex.axis=1,
                                filename=paste(object@identifier,"heatmapMRIaggr",sep="_"),
                                width=1000,height=700,path=NULL,unit="px",res=NA){
            
            ### preparation
            p <- length(param)
            
            carto <- selectContrast(object,param=param,num=num,hemisphere=hemisphere,norm_mu=F,norm_sigma=F,format="matrix")
            if(scale==TRUE){carto <- scale(carto)}
            
            attr(carto,"scaled:center")<-NULL
            attr(carto,"scaled:scale")<-NULL 
            
            carto.corr <- cor(carto,method=method)
            rownames(carto.corr) <- param
            colnames(carto.corr) <- param
            grid <- expand.grid(1:p,1:p)
            
            if(is.null(type) || type %in% c(FALSE,"image","image.plot")==FALSE)
            {stop("heatmapMRIaggr[MRIaggr] : wrong specification of  \'type\' \n",
                  "valid types :  FALSE \"image\" \"image.plot\" \n",
                  "requested type : ",type,"\n")
            }
            
            if(is.null(breaks)){breaks <- seq(min(carto.corr),max(carto.corr),length.out=length(col)+1)}
            
            ### test export
            res.init <- initWindow(window=window,filename=filename,path=path,width=width,height=height,unit=unit,res=res,
                                   n.plot=1,mfrow=1,xlim=NULL,ylim=NULL,
                                   method="boxplotMask[MRIaggr]")
            scale <- res.init$scale
            
            
            #### display
            initDisplayWindow(window=window,filename=filename,path=path,width=width,height=height,scale=scale,res=res,
                              mfrow=NULL,bg=NULL,pty=NULL,mar=NULL,mgp=mgp)
            
            if(type == "image"){
              image(1:p,1:p,carto.corr,col=col,breaks=breaks,
                    main=main,axes=FALSE,xlab="",ylab="")
              axis(1,at=0:(p+1),label=c("",param,""),las=las,cex.axis=cex.axis)
              axis(2,at=0:(p+1),label=c("",param,""),las=las,cex.axis=cex.axis)
              
              if(points.values==TRUE){
                text(x=grid[,1],y=grid[,2],
                     round(as.vector(carto.corr),digit=digit)
                )
              }   
            }
            
            if(type == "image.plot"){
              
              par(mar=mar)
              image(1:p,1:p,carto.corr,col=col,breaks=breaks,
                    main=main,axes=FALSE,xlab="",ylab="")
              axis(1,at=0:(p+1),label=c("",param,""),las=las,cex.axis=cex.axis)
              axis(2,at=0:(p+1),label=c("",param,""),las=las,cex.axis=cex.axis)
              
              if(points.values==TRUE){
                text(x=grid[,1],y=grid[,2],
                     round(as.vector(carto.corr),digit=digit)
                )
              }
              
              test.package <- requireNamespace("fields",quietly=TRUE)
              if(test.package==FALSE){
                stop("heatmapMRIaggr[MRIaggr] : this function with argument type=\"image.plot\" requires to have installed the fields package to work \n")
              }
              
              fields::image.plot(1:p,1:p,carto.corr,col=col,legend.only=TRUE)
              
            }
            
            
            if(!is.null(window) && window %in% c("eps","svg","png","pdf")){
              dev.off()
            }
            
            return(invisible(carto.corr))
            
            
          }
)

setMethod(f ="multiplot",
          signature ="MRIaggr",
          definition = function(object,param,num=NULL,index1=NULL,index2=NULL,index3=NULL,midplane=FALSE,
                                slice_var="k",hemisphere="both",norm_mu=FALSE,norm_sigma=FALSE,as.logical=FALSE,
                                breaks=50,type.breaks="range",palette="terrain.colors",col=NULL,pch=NULL,cex=1,
                                col.NA="lightyellow",pch.NA=8,col.midplane="red",xlim=NULL,ylim=NULL,axes=TRUE,
                                window=FALSE,legend=TRUE,mfrow=NULL,mar=rep(1.5,4),mgp=c(2,0.5,0),pty=NULL,asp=1,bg="lightblue",
                                xlab="",ylab="",main=NULL,num.main=TRUE,cex.main=1.5,
                                quantiles.legend=TRUE,digit.legend=3,cex.legend=1.5,mar.legend=c(2,7,2,2),main.legend=param,
                                filename="multiplot",width=1000,height=500,path=NULL,unit="px",res=NA){
            
            if(is.null(main)){
              main <- paste(paste(param,collapse="-")," : ", object@identifier," - slice ",slice_var,sep="")
            }
            
            #### data ####
            num <- initNum(object=object,num=num,slice_var=slice_var,method="multiplot")
            n.plot <- length(num)
            if(!is.null(legend) && legend==TRUE){n.plot <- n.plot + 1}
            
            coords <- selectCoords(object,hemisphere=hemisphere,num=num,format="data.frame")
            contrast <- selectContrast(object,param=param,hemisphere=hemisphere,num=num,norm_mu=norm_mu,norm_sigma=norm_sigma,format="matrix",slice_var=slice_var)
            contrast <- data.frame(apply(contrast,2,as.numeric) )
            n.px <- nrow(contrast)
            
            if(n.px==0){
              stop("multiplot[data.frame] : \'data\' contains no observation when restricted to \'num\' and \'hemisphere\' \n",
                   "requested slices : ",paste(num,collapse=" ")," \n",
                   "requested hemisphere : ",hemisphere," \n")
            }
            
            #### initialization and tests ####
            
            ## windows
            mar.init <- par()$mar
            
            res.init <- initWindow(window=window,filename=filename,path=path,width=width,height=height,unit=unit,res=res,
                                   n.plot=n.plot,mfrow=mfrow,xlim=xlim,ylim=ylim,
                                   method="multiplot[MRIaggr]")
            scale <- res.init$scale
            mfrow <- res.init$mfrow
            n.graph_par_window <- res.init$n.graph_par_window
            xlim.plot <- res.init$xlim.plot
            ylim.plot <- res.init$ylim.plot
            
            ## color and breaks   
            res.init <- initCol(contrast=contrast,coords=coords,param=param,pch=pch,col=col,palette=palette,breaks=breaks,legend=legend,type.breaks=type.breaks,
                                method="multiplot[MRIaggr]")
            
            contrast <- res.init$contrast 
            palette <- res.init$palette 
            breaks <- res.init$breaks
            col <- res.init$col
            pch <- res.init$pch
            palette_sauve <- res.init$palette_sauve
            breaks_sauve <- res.init$breaks_sauve
            index_duplicated <- res.init$index_duplicated
            index_order <- res.init$index_order
            param <- res.init$param
            
            ## index  
            if(!is.null(index1)){
              res.init <- initIndex(object=object,index=index1,num=num,hemisphere=hemisphere,as.logical=as.logical,
                                    method="multiplot[MRIaggr]",indexNum=1,
                                    cex.default=1,pch.default=20,col.default="red",filter_default="2D_N4")
              index1 <- res.init$coords
              indexindex1 <- res.init$index
              pch_index1 <- res.init$pch
              cex_index1 <- res.init$cex
              col_index1 <- res.init$col
            }
            
            if(!is.null(index2)){
              res.init <- initIndex(object=object,index=index2,num=num,hemisphere=hemisphere,as.logical=as.logical,
                                    method="multiplot[MRIaggr]",indexNum=2,
                                    cex.default=1,pch.default=21,col.default="purple",filter_default="2D_N4")
              index2 <- res.init$coords
              indexindex2 <- res.init$index
              pch_index2 <- res.init$pch
              cex_index2 <- res.init$cex
              col_index2 <- res.init$col
            }
            
            if(!is.null(index3)){
              
              test.intersection1 <- length(index3)==1 && index3=="intersection"
              test.intersection2 <- is.list(index3) && "coords" %in% names(index3) && length(index3$coords)==1 && index3$coords=="intersection"
              
              if(test.intersection1 || test.intersection2){
                
                if(test.intersection1){index3 <- list(coords="intersection")}
                if(is.null(index1)){
                  stop("multiplot[MRIaggr] : \'index3\' cannot request \"interaction\" parameter if \'index1\' is missing  \n")
                }
                if(is.null(index2)){
                  stop("multiplot[MRIaggr] : \'index3\' cannot request \"interaction\" parameter if \'index2\' is missing  \n")
                }
                
                test1.intersect <- indexindex1 %in% indexindex2
                test2.intersect <- indexindex2 %in% indexindex1
                
                index1_sauve <- index1
                index3$coords <- index1_sauve$coords[test1.intersect==TRUE,c("i","j","k")]
                index1 <- index2$coords[test1.intersect==FALSE,c("i","j","k")]
                index2 <- index1_sauve$coords[test2.intersect==FALSE,c("i","j","k")]                      
              }
              
              res.init <- initIndex(object=object,index=index3,num=num,hemisphere=hemisphere,as.logical=as.logical,
                                    method="multiplot[MRIaggr]",indexNum=3,
                                    cex.default=1,pch.default=22,col.default="green",filter_default="2D_N4")
              index3 <- res.init$coords 
              pch_index3 <- res.init$pch
              cex_index3 <- res.init$cex
              col_index3 <- res.init$col     
            }
            
            #### display plot ####
            compteur <- 1
            plot_var <- setdiff(c("i","j","k"),slice_var)
            
            for(iter_num in 1:length(num)){
              
              # device
              if(!is.null(window) && compteur == 1){
                
                if(window %in% c("eps","svg","png","pdf") && iter_num>1){dev.off()}
                filename_all <- paste(object@identifier,"_",filename,"_",param,"(",slice_var,"slice",num[iter_num],"-",min(max(num),num[iter_num+n.graph_par_window-1],na.rm=TRUE),"_",hemisphere,")",sep="")                  
                initDisplayWindow(window=window,filename=filename_all,path=path,width=width,height=height,scale=scale,res=res,
                                  mfrow=mfrow,bg=bg,pty=pty,mar=mar,mgp=mgp)
                
              }
              
              # data
              index_k <- which(coords[,slice_var]==num[iter_num])
              contrastK <- contrast[index_k,,drop=FALSE]
              coordsK <- coords[index_k,plot_var,drop=FALSE]
              if(is.null(col)){colK <- NULL}else{colK <- col[index_k]}
              if(num.main==TRUE){mainK <- paste(main,num[iter_num],sep="")}else{mainK <- main}
              
              # xlim-ylim
              if(is.null(xlim)){                
                if(is.null(coordsK) || nrow(coordsK)==0){
                  xlim.plot <- NULL
                }else{
                  xlim.plot <- c(min(coordsK[,plot_var[1]])-0.5,max(coordsK[,plot_var[1]])+0.5)
                }                 
              }
              
              if(is.null(ylim)){
                if(is.null(coordsK) || nrow(coordsK)==0){
                  ylim.plot <- NULL
                }else{
                  ylim.plot <- c(min(coordsK[,plot_var[2]])-0.5,max(coordsK[,plot_var[2]])+0.5)
                }
              }
              
              plot.test <- plotMRI(contrast=contrastK,coords=coordsK[,plot_var],breaks=breaks,col=colK,palette=palette,
                                   asp=asp,
                                   xlim=xlim.plot,ylim=ylim.plot,pch=pch,cex=cex,axes=axes,col.NA=col.NA,pch.NA=pch.NA,
                                   xlab=xlab,ylab=ylab,main=mainK,cex.main=cex.main)
              
              if(midplane==TRUE){
                pointsHemisphere(object,col=col.midplane)
              }
              
              if(!is.null(index1) && any(index1[,slice_var]==num[iter_num])){
                points(index1[index1[,slice_var]==num[iter_num],plot_var,drop=FALSE],
                       pch=pch_index1,cex=cex_index1,col=col_index1)
              }
              
              if(!is.null(index2) && any(index2[,slice_var]==num[iter_num])){
                points(index2[index2[,slice_var]==num[iter_num],plot_var,drop=FALSE],
                       pch=pch_index2,cex=cex_index2,col=col_index2)
              }
              
              if(!is.null(index3) && any(index3[,slice_var]==num[iter_num])){
                points(index3[index3[,slice_var]==num[iter_num],plot_var,drop=FALSE],
                       pch=pch_index3,cex=cex_index3,col=col_index3)
              }
              
              compteur <- compteur + 1
              if(compteur > n.graph_par_window)
              {compteur <- 1}
              
            }
            
            #### legend ####
            if(is.null(legend)){ # affiche sur une fenetre a part
              compteur <- 1
              mfrow <- c(1,1)
              legend <- TRUE
            }
            
            if(legend==TRUE){
              
              if(!is.null(window) && compteur == 1){
                
                if(window %in% c("eps","svg","png","pdf") && iter_num>1){dev.off()}
                filename_all <- paste(object@identifier,"_",filename,"_",param,"(",slice_var,"slice",min(num),"-",max(num),"_",hemisphere,")-legend",sep="")
                
                initDisplayWindow(window=window,filename=filename_all,path=path,width=width,height=height,scale=scale,res=res,
                                  mfrow=mfrow,bg=bg,pty=pty,mar=NULL,mgp=mgp)
                
              }
              
              if(!is.null(index_duplicated) && !is.null(index_order)){ # palette rgb ou hsv
                col <- col[index_duplicated][index_order]                  
                contrast <- contrast[index_duplicated,1,drop=FALSE][index_order,1,drop=FALSE]
                coords <- coords[index_duplicated,,drop=FALSE][index_order,,drop=FALSE]
                quantiles.legend <- NULL
                
                index_min <- which.min(abs(contrast[,1]-1))
                col[seq(index_min,index_min+round(index_min/25))] <- NA
                index_min <- which.min(abs(contrast[,1]-2))
                col[seq(index_min,index_min+round(index_min/50))] <- NA
              }else{
                if(quantiles.legend==TRUE){
                  quantiles.legend <- quantile(contrast[,1],na.rm=TRUE)
                }else{
                  quantiles.legend <- quantiles.legend
                }
              }
              
              plot.test <- legendMRI(breaks=breaks_sauve,palette=palette_sauve,mar=mar.legend,
                                     cex=cex.legend,main=paste(main.legend,collapse=" "),cex.main=cex.main,quantiles=quantiles.legend,digit=digit.legend)
              
            }
            
            
            
            if(!is.null(window) && window %in% c("eps","svg","png","pdf")){
              dev.off()
            }
            par(mar=mar.init)
            
          }
)

setMethod(f ="plotDistClass",
          signature ="MRIaggr",
          definition = function(object,param,class,num=NULL,hemisphere="both",norm_mu=FALSE,norm_sigma=FALSE,
                                bw.adjust=1,kernel="gaussian",from=NULL,to=NULL,ylim=NULL,
                                window=FALSE,col=1:6,main=NULL,mgp=c(2,0.5,0),type="l",pch=20,lwd=1,x.legend="topright",y.legend=NULL,cex.legend=0.8,
                                filename=paste(object@identifier,"plotDistClass",sep="_"),
                                width=1000,height=700,path=NULL,unit="px",res=NA){
            
            ### preparation
            n.class <- length(class)
            
            if(length(param)!=1){
              stop("plotDistClass[MRIaggr] : wrong specification of \'param\'\n",
                   "only one parameter can be specified \n",
                   "length(param) : ",length(param),"\n")                  
            }
            
            carto <- cbind(selectContrast(object,param=param,num=num,hemisphere=hemisphere,norm_mu=norm_mu,norm_sigma=norm_sigma,format="data.frame"),
                           selectContrast(object,param=class,num=num,hemisphere=hemisphere,format="data.frame"))
            if(is.null(from)){
              from <- min(carto[,param])
            }
            if(is.null(to)){
              to <- max(carto[,param])
            }
            
            density <- list()
            for(iter_class in 1:n.class){
              
              density[[iter_class]] <- density(carto[,param],adjust=bw.adjust,
                                               kernel=kernel,
                                               weights=carto[,class[iter_class]]/sum(carto[,class[iter_class]]),
                                               from=from,to=to)
            }
            max_density <- max(unlist(lapply(density,function(x){max(x$y)})))
            
            ### test export
            res.init <- initWindow(window=window,filename=filename,path=path,width=width,height=height,unit=unit,res=res,
                                   n.plot=1,mfrow=1,xlim=NULL,ylim=NULL,
                                   method="plotDistClass[MRIaggr]")
            scale <- res.init$scale
            
            #### display
            if(!is.null(window)){
              
              initDisplayWindow(window=window,filename=filename,path=path,width=width,height=height,scale=scale,res=res,
                                mfrow=NULL,bg=NULL,pty=NULL,mar=NULL,mgp=mgp)      
              
              if(is.null(ylim)){
                ylim <- c(0,max_density)
              }
              
              plot(NA,NA,xlim=c(from,to),xlab=param,ylab="density",main=main,
                   ylim=ylim)  
              
              for(iter_class in 1:n.class){
                points(density[[iter_class]]$x,density[[iter_class]]$y,
                       type=type,pch=pch,col=col[iter_class],lwd=lwd)  
              }
              
              legend(x=x.legend,y=y.legend,bty="n",legend=class,pch=pch,col=col,
                     cex=cex.legend)
            }
          }
)

setMethod(f ="pointsHemisphere",
          signature ="MRIaggr",
          definition = function(object,col="red",lwd=2,lty=1)
          { 
            midplane <- object@midplane[,c("i","j")]
            points(midplane,type="l",col=col,lwd=lwd,lty=lty)
          }
)

# plotLesion3D
setMethod(f ="plotLesion3D",
          signature ="MRIaggr",
          definition = function(object,mask,edge=FALSE,Neighborhood="3D_N6",as.logical=FALSE,spatial_res=c(1,1,1),
                                xlim=NULL,ylim=NULL,zlim=NULL,type.plot="shapelist3d",px_max=10000,
                                radius=1,type="s",col="red",col.edge="black")
          { test.package <- requireNamespace("rgl",quietly=TRUE)
            if(test.package==FALSE){
              stop("plotLesion3D[MRIaggr] : this function requires to have installed the rgl package to work \n")
            }
            
            # import
            carto <- selectContrast(object,param=mask,coords=TRUE)
            
            if(length(spatial_res)!=3 || is.numeric(spatial_res)==FALSE){
              stop("plotLesion3D[MRIaggr] : wrong specification of \'spatial_res\' \n",
                   "it must be a numeric vector of length 3 \n",
                   "proposed 'spatial_res' : ",paste(spatial_res,collapse=" "),"\n")
            }
            
            carto[,"i"] <- carto[,"i"]*spatial_res[1]
            carto[,"j"] <- carto[,"j"]*spatial_res[2]
            carto[,"k"] <- carto[,"k"]*spatial_res[3]
            
            if(as.logical==TRUE){carto[,mask] <- as.logical(carto[,mask])}
            if(!is.logical(carto[,mask])){
              stop("plotLesion3D[MRIaggr] : type of \'mask\' is not logical \n",
                   "proposed type : ",paste(is(carto[,mask]),collapse=" "),"\n",
                   "to force the conversion to logical set \'as.logical\'= TRUE \n")
            }
            
            index_mask <- which(carto[,mask]==1)
            n.mask <- length(index_mask)
            
            if(edge){
              carto_Wmask <- calcFilter(object,param=mask,filter=Neighborhood,norm=FALSE,trace=FALSE)
              nb_neighbor <- max(carto_Wmask[,paste(mask,Neighborhood,sep="_")])
              index_core <- which(carto_Wmask[,paste(mask,Neighborhood,sep="_")]==nb_neighbor)
              index_edge <- setdiff(index_mask,index_core)
            }
            
            # test
            
            if(type.plot %in% c("plot3d","shapelist3d") ==FALSE)
            {stop("plotLesion3D[MRIaggr] : wrong specification of \'type.plot\' \n",
                  "valid types : \"plot3d\" \"shapelist3d\" \n",
                  "proposed type.plot : ",type.plot,"\n")}
            
            if(n.mask>px_max)
            {stop("plotLesion3D[MRIaggr] : the number of points to be displayed exceed the limit \'pt_max\' \n",
                  "number of points to be diplayed : ",test.nb_pt,"\n",
                  "limit  : ",px_max,"\n")}
            
            if(length(mask)!=1)
            {stop("plotLesion3D[MRIaggr] : argument \'mask\' must not contain several values \n",
                  "number of mask requested : ",length(mask),"\n")}
            
            
            # preparation 
            
            if(is.null(xlim))
            { xlim <- range(carto[index_mask,"i"])}
            
            if(is.null(ylim))
            { ylim <- range(carto[index_mask,"j"])}
            
            if(is.null(zlim))
            { zlim <- range(carto[index_mask,"k"])}
            
            # trace
            if(type.plot=="plot3d"){
              
              if(edge==FALSE){
                rgl::plot3d(carto[index_mask,c("i","j","k")],xlim=xlim,ylim=ylim,zlim=zlim,
                            col=col,type=type,radius=radius)
              }else{
                rgl::plot3d(carto[index_core,c("i","j","k")],xlim=xlim,ylim=ylim,zlim=zlim,
                            col=col,type=type,radius=radius)
                rgl::points3d(carto[index_edge,c("i","j","k")],
                              col=col.edge)
              }
            }
            
            if(type.plot=="shapelist3d"){
              M <- diag(spatial_res)
              if(edge==FALSE){
                rgl::shapelist3d(rgl::cube3d(trans=M), x=carto[index_mask,"i"], y=carto[index_mask,"j"], z=carto[index_mask,"k"],
                                 size=radius,xlim=xlim,ylim=ylim,zlim=zlim,col=col)                
                rgl::axes3d(c('x','y','z'))
                rgl::title3d(main=mask,xlab="i",ylab="j",zlab="k")
                
              }else{
                rgl::shapelist3d(rgl::cube3d(trans=M), x=carto[index_edge,"i"], y=carto[index_edge,"j"], z=carto[index_edge,"k"],
                                 size=radius,col=col.edge)
                rgl::axes3d(c('x','y','z'))
                rgl::title3d(main=mask,xlab="i",ylab="j",zlab="k")
                
              }
            }
            
          }
)

# plotTableLesion
setMethod(f ="plotTableLesion",
          signature ="MRIaggr",
          definition = function(object,mask,num=NULL,type="matplot",
                                window=FALSE,col=1:5,lty=1:5,lwd=1,mgp=c(2,0.5,0),mar=rep(3,4),
                                main=paste("lesion - ",object@identifier,sep=""), cex.legend=1,cex.main=1,cex.axis=1,cex.lab=1,
                                filename=paste(object@identifier,"plotTableLesion",sep="_"),
                                width=1000,height=700,path=NULL,unit="px",res=NA)
          { 
            
            num <- initNum(object,num=num,method="plotTableLesion")
            
            table_lesion <- selectTable(object,type="lesion")
            n.num <- length(num)
            p <- length(mask)
            
            if(any(mask %in% names(table_lesion)==FALSE)){
              stop("plotTableLesion[MRIaggr] : wrong specification of \'mask\' \n",
                   "valid masks : ",paste(names(table_lesion),collapse=" "),"\n",
                   "invalid requested masks : ",mask[mask %in% names(table_lesion)==FALSE],"\n")
            }
            
            if(type %in% c("matplot","evolution") == FALSE){
              stop("plotTableLesion[MRIaggr] : wrong specification of \'type\' \n",
                   "valid types : \"matplot\" \"evolution\" \n",
                   "requested type : ",type,"\n")
            }
            
            #### test export
            res.init <- initWindow(window=window,filename=filename,path=path,width=width,height=height,unit=unit,res=res,
                                   n.plot=1,mfrow=1,xlim=NULL,ylim=NULL,
                                   method="boxplotMask[MRIaggr]")
            scale <- res.init$scale
            
            #### display
            initDisplayWindow(window=window,filename=filename,path=path,width=width,height=height,scale=scale,res=res,
                              mfrow=NULL,bg=NULL,pty=NULL,mar=NULL,mgp=mgp)      
            
            if(type=="matplot"){
              matplot(table_lesion[num,mask],
                      type="l",xlab="k",ylab="observation number",lty=lty,lwd=lwd,col=col,main=main,axes=FALSE,
                      cex.main=cex.main,cex.lab=cex.lab)
              axis(1,at=1:n.num,label=num,cex.axis=cex.axis)
              axis(2,cex.axis=cex.axis)
              
              legend("topleft",legend=mask,lty=lty,col=col,bty="n",cex=cex.legend)
            }else{
              if(!is.null(window)){
                par(mfrow=c(1,p))
                if(!is.null(mar)){par(mar=mar)}
              }
              
              lwd.init <- lwd
              col.init <- col
              for(iter_p in 1:p){
                mask_max <- max(table_lesion[num,mask[iter_p]])
                col <- rep(col.init[1],n.num)
                
                if(iter_p >1){
                  croissance <-  table_lesion[num,mask[iter_p]] - table_sauve
                  lwd <- lwd.init + 10*abs(croissance)/max(abs(croissance))
                  col[croissance>0] <- col.init[2]
                  col[croissance<0] <- col.init[3]
                }
                
                plot(NA,NA,xlim=c(0,mask_max),ylim=range(num),
                     xlab=mask[iter_p],ylab="k",main=main,
                     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab)
                segments(y0=num,
                         y1=num,
                         x0=0,
                         x1=table_lesion[num,mask[iter_p]],
                         lwd=lwd,
                         col=col)         
                
                table_sauve <- table_lesion[num,mask[iter_p]]                
              }
            }
            
            if(!is.null(window) && window %in% c("eps","svg","png","pdf")){
              dev.off()
            }
            
            
          }
)

setMethod(f ="outlineMRIaggr",
          signature ="MRIaggr",
          definition = function(object,param,index1=NULL,num=NULL,hemisphere="both",as.logical=FALSE,
                                xlim=NULL,ylim=NULL,legend=FALSE,palette="terrain.colors",col=NULL,breaks=25,
                                fill=TRUE,n=50,sequential=TRUE,min_dist=1,operator_index1="none",
                                col.outline=c("blue","red","grey"),pch=20,cex=c(0.75,1,0.75),
                                trace=TRUE,name_newparam="userMask",update.object=FALSE,overwrite=FALSE){
            
            
            num <- initNum(object=object,num=num,method="outlineMRIaggr")
            param <- initParameter(object=object,param=param,test=TRUE,init=TRUE,method="outlineMRIaggr")      
            
            newparam <- selectCoords(object,coords=c("i","j","k","index"))
            edge <- matrix(NA,nrow=0,ncol=3)
            surface <- matrix(NA,nrow=0,ncol=3)
            
            if(length(name_newparam)!=1){
              stop("outlineMRIaggr[MRIaggr] : wrong specification of \'name_newparam\' \n",
                   "must be a character (length 1)  \n",
                   "length(name_newparam) : ",length(name_newparam),"\n")
            }
            
            if(is.null(index1)){operator_index1 <- "none"}
            if(operator_index1 %in% c("difference","union","intersection","none") == FALSE){
              stop("outlineMRIaggr[MRIaggr] : wrong specification of \'operator_index1\' \n",
                   "valid operators : \"difference\" \"union\" \"intersection\" \"none\" \n",
                   "proposed operator : ",operator_index1,"\n")
            }
            
            if(!is.null(index1)){
              index1 <- initIndex(object=object,index=index1,num=num,hemisphere=hemisphere,as.logical=as.logical,indexNum=1,
                                  cex.default=1,pch.default=20,col.default="black",filter_default="2D_N4",method="outlineMRIaggr")
            }
            
            graphics.off()
            for(iter_num in 1:length(num)){
              
              # display 
              repeat.slice <- 1
              while(!is.na(repeat.slice) && repeat.slice==1){
                
                multiplot(object,param=param,num=num[iter_num],hemisphere=hemisphere,xlim=xlim,ylim=ylim,index1=index1,
                          col=col,breaks=breaks,palette=palette,legend=legend)
                
                res_outline <- outline(n=n,sequential=sequential,min_dist=min_dist,
                                       col=col.outline,pch=pch,cex=cex)
                
                if(is.null(res_outline$edge)){
                  repeat.slice <- readline("Do you want to start again ? (answer 1)\n Otherwise the slice will be skipped  ")                
                }else{
                  repeat.slice <- FALSE
                }
                
              }
              
              if(!is.null(res_outline$edge)){
                edge <- rbind(edge,
                              cbind(res_outline$edge[,c("i","j")],k=num[iter_num]))
              }
              
              if(!is.null(res_outline$surface)){
                surface <- rbind(surface,
                                 cbind(res_outline$surface[,c("i","j")],k=num[iter_num]))
              }
            }
            
            newparam <- merge(x=newparam,y=data.frame(edge,edge=TRUE),by=c("i","j","k"),all.x=TRUE,all.y=FALSE)
            newparam <- merge(x=newparam,y=data.frame(surface,surface=TRUE),by=c("i","j","k"),all.x=TRUE,all.y=FALSE)
            
            newparam$surface[is.na(newparam$surface)] <- FALSE
            newparam$edge[is.na(newparam$edge)] <- FALSE
            
            if(operator_index1 != "none"){
              newparam <- merge(x=newparam,y=data.frame(index1$coords,index1=TRUE),by=c("i","j","k"),all.x=TRUE,all.y=FALSE)
              newparam$index1[is.na(newparam$index1)] <- FALSE
              
              if(operator_index1 == "difference"){
                newparam$edge <- as.logical(newparam$edge+newparam$index1-2*newparam$edge*newparam$index1)
                newparam$surface <- as.logical(newparam$surface+newparam$index1-2*newparam$surface*newparam$index1)
              }
              if(operator_index1 == "intersection"){
                newparam$edge <- newparam$edge*newparam$index1 
                newparam$surface <- newparam$surface*newparam$index1 
              }
              if(operator_index1 == "union"){
                newparam$edge <- as.logical(newparam$edge+newparam$index1)
                newparam$surface <- as.logical(newparam$surface+newparam$index1)
              }
            }
            
            newparam <- newparam[order(newparam$index),]
            
            if(fill==TRUE){  
              newparam <- eval(parse(text=paste("data.frame(newparam,",name_newparam,"=newparam[,\"surface\"])",sep="")))
            }else{
              newparam <- eval(parse(text=paste("data.frame(newparam,",name_newparam,"=newparam[,\"edge\"])",sep="")))
            }
            
            return(list(res=newparam,
                        trace=trace,
                        name_newparam=name_newparam,
                        update.object=update.object,
                        overwrite=overwrite)
            )
          }           
)

setMethod(f ="show",
          signature ="MRIaggr",
          definition = function(object){
            
            summary(object,trace=0)
            
          }
)


setMethod(f ="summary",
          signature ="MRIaggr",
          definition = function(object,param=FALSE,clinic=FALSE,descStats=FALSE,history=FALSE,trace=1){
            if(!is.logical(param) || length(param)!=1){
              stop("summary[MRIaggr] : incorrect specification of \'param\' \n",
                   "must be a logical value \n",
                   "length(param) : ",length(param)," | class(param) : ",class(param),"\n")
            }
            
            if(!is.logical(clinic) || length(clinic)!=1){
              stop("summary[MRIaggr] : incorrect specification of \'clinic\' \n",
                   "must be a logical value \n",
                   "length(clinic) : ",length(clinic)," | class(clinic) : ",class(clinic),"\n")
            }
            
            if(!is.logical(descStats) || length(descStats)!=1){
              stop("summary[MRIaggr] : incorrect specification of \'descStats\' \n",
                   "must be a logical value \n",
                   "length(descStats) : ",length(descStats)," | class(descStats) : ",class(descStats),"\n")
            }
            
            # load 
            p <- length(selectParameter(object))
            n <- selectN(object)
            Id <- selectIdentifier(object)            
            D <- selectFieldDim(object)
            Size <- as.matrix(selectVoxelDim(object))
            hemispheres <- selectHemispheres(object)
            midplane <- selectMidplane(object)
            
            ls_descStats <- selectDescStats(object)
            length.descStats <- length(ls_descStats)
            
            df.clinic <- selectClinic(object)
            ncol.clinic <- ncol(df.clinic)
            
            test.history <- selectHistory(object)
            
            default_value <- selectDefault_value(object)
            ncol.TableLesion <- ncol(selectTable(object,type="lesion"))
            ncol.TableReperfusion <- ncol(selectTable(object,type="reperfusion"))
            ncol.TableHypoperfusion <- ncol(selectTable(object,type="hypoperfusion"))
            
            test.normalization <- length(selectNormalization(object))>0
            
            # print
            cat("Object of class \'MRIaggr\' with identifier : ",Id,"\n",sep="")
            
            cat("  # image dimensions (i,j,k) : ",paste(D,collapse="x")," voxels \n",sep="")
            
            cat("  # voxel dimensions (i,j,k unit) : ",paste(Size[1:3],collapse="x")," ",Size[4]," \n",sep="")
            
            if(trace==1 || history==TRUE){
              cat("  # history : ",length(test.history)," calculations performed on the object \n",sep="")
              if(history==TRUE){
                cat(paste(names(test.history),collapse=" ")," \n",sep=" ")
              }
            }
            
            cat("  # number of observations : ",n," \n",sep="")
            
            cat("  # number of contrast parameters : ",p," \n",sep="")
            if(param==TRUE){
              str(selectContrast(object),max.level=1,give.attr=FALSE)
            }
            
            if(trace==1 || ncol.TableLesion>0 || ncol.TableHypoperfusion>0 || ncol.TableReperfusion>0){
              cat("  # tables : lesion ",if(ncol.TableLesion>0){paste(ncol.TableLesion," columns",sep="")}else{"empty"},
                  ", hypoperfusion ",if(ncol.TableHypoperfusion>0){paste(ncol.TableHypoperfusion," columns",sep="")}else{"empty"},
                  ", reperfusion ",if(ncol.TableReperfusion>0){paste(ncol.TableReperfusion," columns",sep="")}else{"empty"},"\n",sep="")
            }
            
            if(trace==1 || test.normalization>0){
              cat("  # normalization values : ",if(test.normalization>0){"present"}else{"empty"},"\n",sep="")
            }
            
            if(trace==1 || hemispheres$left!="undefined" || hemispheres$right!="undefined"){
              cat("  # hemisphere left : ",hemispheres$left," | hemisphere right : ",hemispheres$right,"\n",sep="")
            }
            
            if(trace==1 || sum(!is.na(midplane))>0 ){
              cat("  # midplane coordinates : ",if(sum(!is.na(midplane))>0){"present"}else{"empty"},"\n",sep="")
            }
            
            if(trace==1 || length.descStats>0 ){
              cat("  # ls_descStats : ",if(length.descStats==0){"empty"}else{if(descStats==FALSE){paste(length.descStats," elements",sep="")}},"\n",sep="")
              if(length.descStats>0 && descStats==TRUE){
                str(ls_descStats,max.level=1,give.attr=FALSE)              
              }
            }
            
            if(trace==1 || ncol.clinic>0 ){
              cat("  # clinic : ",if(ncol.clinic==0){"empty"}else{if(clinic==FALSE){paste(ncol.clinic," parameters",sep="")}},"\n",sep="")
              if(ncol.clinic>0 && clinic==TRUE){
                str(df.clinic,max.level=1,give.attr=FALSE)              
              }
            }
            
          }
)

#### const. ####

setMethod(f ="constCompressMRIaggr",
          signature ="MRIaggr",
          definition = function(object,factor,param=NULL,
                                mask=NULL,threshold=0.49,trace=FALSE)
          { 
            data.initial <- object@contrast
            n_i.initial <- object@fieldDim$i
            n_j.initial <- object@fieldDim$j
            n_slices <- object@fieldDim$k
            param <- initParameter(object,param=param,init=TRUE,test=TRUE,
                                   accept.coords=FALSE,accept.index=FALSE,method="constCompressMRIaggr")
            n.param <- length(param)
            
            #### tests ####
            
            if(factor<=1){
              stop("constCompressMRIaggr[MRIaggr] : wrong specification of \'factor\' \n",
                   "\'factor\' must be over 1 \n",
                   "proposed \'factor\' : ",factor,"\n")
            }
            
            if(object@fieldDim$i %% factor>0 || object@fieldDim$j %% factor>0){
              stop("constCompressMRIaggr[MRIaggr] : wrong specification of \'factor\' \n",
                   "@fieldDim$i or @fieldDim$j is not a multiple of \'factor\' \n",
                   "@fieldDim$i/factor : ",object@fieldDim$i/factor,"\n",
                   "@fieldDim$j/factor : ",object@fieldDim$j/factor,"\n")
            }
            
            test.character <- sapply(1:n.param,function(x){is.character(data.initial[,param[x]])})            
            if(any(test.character)){              
              stop("constCompressMRIaggr[MRIaggr] : wrong specification of \'param\' \n",
                   "there ",if(sum(test.character)==1){"is one parameter which contains"}else{"are parameter which contain"}," character values \n",
                   "(\"",paste(param[test.character],collapse="\" \""),"\") \n")
            }
            
            #### initialization ####
            
            n_i.final <- object@fieldDim$i/factor
            n_j.final <- object@fieldDim$j/factor
            n_px.px <- factor^2
            
            data.final <- data.frame(matrix(NA,nrow=n_i.final*n_j.final*n_slices,ncol=ncol(object@contrast)))
            names(data.final) <- names(data.initial)
            
            #### mise en place ####
            seq_i <- matrix(NA,ncol=n_i.final,nrow=factor)
            for(iter_k in 1:factor){
              seq_i[iter_k,] <- seq(iter_k,by=factor,length.out=n_i.final) 
            }
            
            seq_j <- matrix(NA,ncol=n_j.final,nrow=factor)
            for(iter_k in 1:factor){
              seq_j[iter_k,] <- seq(iter_k,by=factor,length.out=n_j.final) 
            }
            
            index_carto <- 0
            
            #### loop ####
            if(trace){cat("slice : ")}
            
            for(iter_slice in 1:n_slices){
              if(trace){cat(iter_slice," ",sep="")}
              
              data.slice <- df2array(contrast=selectContrast(object,num=iter_slice,param=param,format="data.frame"),
                                     coords=selectCoords(object,num=iter_slice,c("i","j")),
                                     range.coords=c(n_i.initial,n_j.initial))$contrast
              
              for(iter_param in 1:n.param){              
                nom_param <- param[iter_param]
                matrix.init_tempo <- data.slice[[iter_param]]
                matrixNA.init_tempo <- is.na(matrix.init_tempo)
                matrix.col_tempo <- matrix(NA,nrow=n_i.final,ncol=n_j.initial)
                matrixNA.col_tempo <- matrix(NA,nrow=n_i.final,ncol=n_j.initial)
                matrix.final_tempo <- matrix(NA,nrow=n_i.final,ncol=n_j.final)
                matrixNA.final_tempo <- matrix(NA,nrow=n_i.final,ncol=n_j.final)
                
                for(iter_i in 1:n_i.final){
                  matrix.col_tempo[iter_i,] <- colSums(matrix.init_tempo[seq_i[,iter_i],])
                  matrixNA.col_tempo[iter_i,] <- colSums(matrixNA.init_tempo[seq_i[,iter_i],])
                }
                
                for(iter_j in 1:n_j.final){
                  matrix.final_tempo[,iter_j] <- rowSums(matrix.col_tempo[,seq_j[,iter_j]])
                  matrixNA.final_tempo[,iter_j] <- rowSums(matrixNA.col_tempo[,seq_j[,iter_j]])
                }
                
                index_NA <- matrixNA.final_tempo > n_px.px/2
                
                matrix.final_tempo <- matrix.final_tempo/(n_px.px - matrixNA.final_tempo)
                
                if(sum(index_NA)>0){matrix.final_tempo[index_NA] <- NA}
                new_data <- array2df(matrix.final_tempo,name_newparam=nom_param,names_coords=c("i","j"), na.rm = FALSE)
                
                
                if(iter_param==1){index_carto <- max(index_carto) + 1:nrow(new_data)}                
                data.final[index_carto,nom_param] <- new_data[,nom_param]
                
                if(iter_param==1){
                  data.final[index_carto,c("i","j")] <- new_data[,c("i","j")]  
                  data.final[index_carto,"k"] <- iter_slice
                }
                
              }
            }
            if(trace){cat("\n")}
            
            # Reglement des parametres binaires :
            if(!is.null(mask)){
              for(iter_param in mask){
                data.final[,iter_param] <- as.numeric(data.final[,iter_param] > threshold)
              }
            }
            
            #### export ####
            y <- new(Class= "MRIaggr",
                     identifier = object@identifier,
                     contrast=data.final,
                     clinic=object@clinic,                     
                     fieldDim=data.frame(i=n_i.final,j=n_j.final,k=object@fieldDim$k),
                     voxelDim=data.frame(i=object@voxelDim$i*factor,j=object@voxelDim$j*factor,k=object@voxelDim$k,unit=object@voxelDim$unit, stringsAsFactors = FALSE),
                     default_value=object@default_value,
                     history=c(object@history,
                               list(constCompressMRIaggr=list(call=match.call(call = sys.call(sys.parent())),date=date()))
                     ),
                     normalization=object@normalization,
                     midplane=object@midplane/factor,
                     hemispheres=object@hemispheres,
                     table_lesion=object@table_lesion,
                     table_reperfusion=object@table_reperfusion,
                     table_hypoperfusion=object@table_hypoperfusion,
                     ls_descStats=object@ls_descStats
            )
            
            if(trace){
              cat("constCompressMRIaggr[MRIaggr] : MRIaggr has been compressed from (x,y) = ",object@fieldDim$i," ",object@fieldDim$j,"\n",
                  "                                                              to (x,y) = ",n_i.final," ",n_j.final,"\n",sep="")
            }
            
            return(y)           
            
          }
)

setMethod(f ="constReduceMRIaggr",
          signature ="MRIaggr",
          definition = function(object,mask,as.logical=FALSE,keep.index=TRUE)
          { 
            if(length(mask)==1 && is.character(mask)){                            
              mask <- selectContrast(object,param=mask,format="vector")              
              supprContrast(object) <- "mask"              
            }else{
              if(selectN(object)!=length(mask))
              {stop("constReduceMRIaggr[MRIaggr] : \'mask\' is incompatible with \'object@contrast\' \n",
                    "selectN(object) : ",selectN(object),"\n",
                    "length(mask) : ",length(mask),"\n" )}
            }
            if(as.logical==TRUE){mask <- as.logical(mask)}
            if(!is.logical(mask)){
              stop("constReduceMRIaggr[MRIaggr] : type of \'mask\' is not logical \n",
                   "proposed type : ",paste(is(mask),collapse=" "),"\n",
                   "to force the conversion to logical set \'as.logical\'= TRUE \n")
            }
            
            index_mask <- which(mask)
            
            if(keep.index==TRUE){
              allocDescStats(object,name="index_sauve") <- selectContrast(object,param="index",format="vector")
            }
            
            
            y <- new(Class= "MRIaggr",
                     identifier = object@identifier,
                     contrast=object@contrast[index_mask,,drop=FALSE],
                     clinic=object@clinic,
                     fieldDim=object@fieldDim,
                     voxelDim=object@voxelDim,
                     default_value=object@default_value,
                     history=c(object@history,
                               list(constReduceMRIaggr=list(call=match.call(call = sys.call(sys.parent())),date=date()))
                     ),
                     normalization=object@normalization,
                     midplane=object@midplane,
                     hemispheres=object@hemispheres,
                     table_lesion=object@table_lesion,
                     table_reperfusion=object@table_reperfusion,
                     table_hypoperfusion=object@table_hypoperfusion,
                     ls_descStats=object@ls_descStats)
            
            cat("constReduceMRIaggr[MRIaggr] : MRIaggr_red has been created \n")
            
            return(y)
          }
)

setMethod(f ="writeMRIaggr",
          signature ="MRIaggr",
          definition = function(object,param,num=NULL,norm_mu=FALSE,norm_sigma=FALSE,range.coords=NULL,default_value=NA,
                                file, format, gzipped=TRUE, verbose=FALSE, size = "NA_integer_")
          { 
            data <- selectContrast(object,param=param,num=num,norm_mu=norm_mu,norm_sigma=norm_sigma)
            coords <- selectCoords(object,num=num,format="matrix")
            
            array <- df2array(data,coords=coords,format="any",
                              default_value=default_value,range.coords=range.coords)$contrast[[1]]
            
            writeMRI(data=array, file=file, format=format, gzipped=gzipped, verbose=verbose, size = size)              
            
          }
)


#### init. ####

setMethod(f ="initNum",
          signature ="MRIaggr",
          definition = function(object,num,test=TRUE,init=TRUE,slice_var="k",method)
          { 
            if(init==TRUE){
              if(is.null(num))
              {num <- seq(1,object@fieldDim[,slice_var])}
            }
            
            if(test==TRUE){
              
              if(length(slice_var)!=1 || slice_var %in% c("i","j","k") == FALSE){
                stop(method,"[MRIaggr] :  wrong specification of \'slice_var\' \n",
                     "valid values : \"i\" \"j\" \"k\" \n",
                     "proposed value : ",paste(slice_var,collapse=""),"\n")
              }
              
              if(any(num %in% seq(1,object@fieldDim[,slice_var]) == FALSE) || length(unique(num))!=length(num))
              {stop(method,"[MRIaggr] :  wrong specification of \'num\' (slice ",slice_var,") \n",
                    "valid values : ",paste(seq(1,object@fieldDim[,slice_var]),collapse=" "),"\n",
                    "wrong proposed values : ",paste(num[num %in% seq(1,object@fieldDim[,slice_var]) == FALSE],collpase=" "),"\n")}
            }
            
            return(num=num)
          }
)

setMethod(f ="initParameter",
          signature ="MRIaggr",
          definition = function(object,param,test=TRUE,init=FALSE,accept.coords=TRUE,accept.mask=TRUE,accept.index=TRUE,
                                arg_name="param",long_name="parameters",method)
          { 
            if(init==TRUE){
              if(is.null(param)){
                param <- selectParameter(object,mask=accept.mask)
              }
              if(is.numeric(param)){
                param <- selectParameter(object,mask=accept.mask)[param]
              }
            }
            
            coords <- NULL
            if(accept.coords==TRUE){
              coords <- c(coords,c("i","j","k"))
            }
            if(accept.mask==TRUE && "mask" %in% names(object@contrast)){
              coords <- c(coords,c("mask"))
            }
            if(accept.index==TRUE){
              coords <- c(coords,c("index"))
            }
            
            param_data <- c(coords,selectParameter(object,mask=FALSE))
            if(test==TRUE){
              if(any(param %in% param_data == FALSE) || length(unique(param))!=length(param))
              {stop(method,"[MRIaggr] : wrong specification of \'",arg_name,"\' \n",
                    "unknown ",long_name," : ",paste(param[param %in% param_data == FALSE],collapse=" ")," \n",
                    "duplicated ",long_name," : ",paste(unique(param[duplicated(param)]),collapse=" "),"\n",
                    "available ",long_name," : ",paste(param_data,collapse=" "),"\n")
              }
            }
            
            return(invisible(param))
          }
)

setMethod(f ="initMask",
          signature ="MRIaggr",
          definition = function(object,mask,test=TRUE,init=TRUE,format="data.frame",
                                arg_name="mask",long_name="mask",method)
          { 
            
            n.mask <- length(mask)
            M.mask <- selectContrast(object,param=mask,format="matrix")
            
            # test logical 
            if(init==TRUE){
              M.mask <- apply(M.mask,2,function(x){as.logical(x)})
            }
            
            test.logical <- apply(M.mask,2,is.logical)
            
            if(test==TRUE){
              if(any(test.logical==FALSE)){
                stop(method,"[MRIaggr] : wrong specification of \'",arg_name,"\' \n",
                     "non logical values in : \'",paste(mask[test.logical==FALSE],collapse="\' \'"),"\' \n",
                     "proposed values : ",paste( round(head(unique(as.vector( M.mask[,which(test.logical==FALSE)] ))),2), collapse=" " ),"\n",
                     "to force the conversion to logical set \'as.logical\'= TRUE \n")
              }
            }
            
            if(format=="data.frame"){
              M.mask <- as.data.frame(M.mask)
              names(M.mask) <- mask
            }
            
            if(format=="vector"){
              M.mask <- as.vector(M.mask)              
            }
            
            return(invisible(M.mask))
          }
)
