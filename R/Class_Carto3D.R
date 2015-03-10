#**********************************************************************
#**********************************************************************
#*************         0A Classe Carto3D           *******************
#**********************************************************************
#**********************************************************************
#
###### Sommaire #################
# A) Class Definition
# B) Selecters
# C) Affecters
# D) Methods

##### A) Class Definition #############################################################

setClass(
  
  Class="Carto3D",
  
  representation(
    identifier  = "character",   # id patient 
    parameter = "character",     # parametre d intensite de la carto
    contrast = "data.frame",      # matrice integrant toutes les cartos d interet de dim (L,nb,l)
    voxelDim = "data.frame",    # vecteur contenant les dimensions des cartos 2D (l,L,nb)
    voxelSize = "data.frame",    # vecteur contenant les dimensions des cartos 2D (l,L,nb)
    default_value ="character"
  ),
  
  validity = function(object){
    #cat("--- Carto3D : checking --- ")
    
    
    if(ncol(object@voxelDim)!=3 || any(names(object@voxelDim)!=c("i","j","k")))
    { stop("validity[Carto3D] : wrong specification of \'voxelDim\' \n",
           "required names : \"i\" \"j\" \"k\" \n",
           "proposed names : ",names(object@voxelDim),"\n")   
    }
    
    if(ncol(object@voxelSize)!=4 || any(names(object@voxelSize)!=c("i","j","k","unit")))
    { stop("validity[Carto3D] : wrong specification of \'voxelSize\' \n",
           "required names : \"i\" \"j\" \"k\" \"unit\" \n",
           "proposed names : ",names(object@voxelSize),"\n")   
    }
    
    if(ncol(object@contrast) != 4 )
    {  stop("validity[Carto3D] : wrong specification of \'contrast\' \n",
            "required number of columns : 4 \n",
            "proposed number of columns : ",object@voxelDim$i,"\n")
    }
    
    if(nrow(object@contrast) != prod(object@voxelDim) )
    {  stop("validity[Carto3D] : wrong specification of \'contrast\' \n",
            "required number of lines (prod(voxelDim))  : ",prod(object@voxelDim),"\n",
            "proposed number of lines : ",nrow(object@contrast),"\n")
    }    
    
    if( any(names(object@contrast) != c("i","j","k",object@parameter)) )
    {  stop("validity[Carto3D] : wrong specification of \'contrast\' \n",
            "required names  : \"i\" \"j\" \"k\" \"",object@parameter,"\" \n",
            "proposed names : ",paste(names(object@contrast),collapse=" "),"\n")
    }
    
    
    #cat(" : Valid Carto3D \n")
    return(TRUE)} 
  
  
)

# Initiateur

setMethod(
  f="initialize",
  signature="Carto3D",
  definition=function(.Object,identifier ,parameter, contrast, voxelDim, voxelSize, default_value){
    
    if(missing(identifier ))
    {stop("initialize[Carto3D] : \'identifier\' is missing \n")}
    .Object@identifier  <- identifier 
    
    if(missing(parameter)){
      if(!missing(contrast) && is.data.frame(contrast) && ncol(contrast)==4){
        parameter <- names(contrast)[4]
      }else{
        stop("initialize[Carto3D] : \'parameter\' is missing \n")
      }
    }
    .Object@parameter <- parameter
       
    if(missing(contrast))
    {stop("initialize[Carto3D] : \'contrast\' is missing \n")}
    
    
    if(!is.data.frame(contrast)){
      
      if(length(dim(contrast)) %in% 3:4 == FALSE)
      {stop("initialize[Carto3D] : wrong specification of \'contrast\'\n",
            "required dimension of \'contrast\' : 3 or 4 \n",
            "proposed dimension of \'contrast\' : ",length(dim(contrast)),"\n")}
      
      
      if(length(dim(contrast))==4 && dim(contrast)[4]==1)
      {contrast <- contrast[,,,1,drop=TRUE]}
      
      if(missing(voxelDim))
      {voxelDim <- data.frame(i=dim(contrast)[1],j=dim(contrast)[2],k=dim(contrast)[3]) }
      
      if(missing(voxelSize))
      {voxelSize <- data.frame(i=NA,j=NA,k=NA,unit=NA, stringsAsFactors = FALSE) }
      
      contrast <- data.frame(cbind(which(array(TRUE,dim=dim(contrast)),arr.ind=TRUE),
                                  as.vector(contrast)))
      names(contrast) <- c("i","j","k",parameter)
    }
    
    .Object@contrast <- contrast
    
    if(missing(voxelDim))
    {stop("initialize[Carto3D] : \'voxelDim\' is missing \n")}
    .Object@voxelDim <- voxelDim
    
    .Object@voxelSize <- voxelSize
    
    if(!missing(default_value))
    {.Object@default_value <- default_value}else{
      .Object@default_value <- "NA"
    }
    
    validObject(.Object)
    return(.Object)
  }
)

#####  B) Selecters #############################################################

# selectContrast
setMethod(f ="selectContrast", 
          signature = "Carto3D",
          definition = function(object,num=NULL,na.rm=FALSE,coords=TRUE,
                                format="any")
          {     
            
            num <- initNum(object=object,num=num)
            if( format %in% c("any","matrix","data.frame") == FALSE )
            {  stop("selectContrast[Carto3D] : wrong specification of \'format\' \n",
                    "valid formats : \"any\" \"matrix\" \"data.frame\" \n",
                    "requested format : ",format,"\n")
            }           
            
            res <- object@contrast[object@contrast$k %in% num,]
            
            if(na.rm==TRUE){
              res <- res[is.na(res[,selectParameter(object)])==FALSE,]
            }
            
            if(coords==FALSE){
              res <- res[,selectParameter(object)]
            }
            
            if(format=="matrix")
            { res <- as.matrix(res)  }
            
            if(format=="data.frame")
            { res <- as.data.frame(res)  }
            
            return(res)
            
          }
)

# selectCoords
setMethod(f ="selectCoords",
          signature ="Carto3D",
          definition = function(object,coords=c("i","j","k"),num=NULL,format="any")
          { 
            if(any(coords %in% c("i","j","k")==FALSE)){
              stop("selectCoords[Carto3D] : wrong specification of \'coords\' \n",
                   "must be one of the following : \"i\" \"j\" ou \"k\" \n",
                   "proposed \'coords\' : ",paste(coords,sep=" "),"\n")
            }
            
            if( format %in% c("any","matrix","data.frame") == FALSE )
            {  stop("selectCoords[Carto3D] : wrong specification of \'format\' \n",
                    "valid formats : \"any\" \"matrix\" \"data.frame\" \n",
                    "requested format : ",format,"\n")
            } 
            
            res <- selectContrast(object,coords=TRUE,num=num,
                                format=format)[,coords]
            
            
            if(format=="matrix")
            { res <- as.matrix(res)  }
            
            if(format=="data.frame")
            { res <- as.data.frame(res)  }
            
            return(res)
            
          }
)

# selectDefault_value
setMethod(f ="selectDefault_value",
          signature ="Carto3D",
          definition = function(object)
          { return(object@default_value)  }
)

# selectIdentifier 
setMethod(f ="selectIdentifier",
          signature ="Carto3D",
          definition = function(object)
          { return(object@identifier )   }
)

# selectParameter
setMethod(f ="selectParameter",
          signature ="Carto3D",
          definition = function(object)
          {return(object@parameter)    }
)

# selectVoxelDim
setMethod(f ="selectVoxelDim",
          signature ="Carto3D",
          definition = function(object)
          { return(object@voxelDim)  }
)

#####  C) Affecters #############################################################




#####  D) Methods #############################################################


#### plot ####

setMethod(f ="multiplot",
          signature ="Carto3D",
          definition = function(object,num=NULL,
                                breaks=50,type.breaks="range",palette="terrain.colors",col=NULL,pch=NULL,cex=1,
                                col.NA="lightyellow",pch.NA=8,xlim=NULL,ylim=NULL,axes=TRUE,
                                window=FALSE,legend=TRUE,mfrow=NULL,mar=rep(1.5,4),mgp=c(2,0.5,0),pty=NULL,asp=1,bg="lightblue",
                                xlab="",ylab="",main=NULL,num.main=TRUE,cex.main=1.5,
                                quantiles.legend=TRUE,digit.legend=3,cex.legend=1.5,mar.legend=c(2,7,2,2),
                                filename="multiplot",width=1000,height=700,path=NULL,unit="px",res=NA)
          { 
            if(is.null(main)){
              main <- paste(selectParameter(object)," : ",selectIdentifier(object)," - slice",sep="")
            }
            if(filename=="auto"){
              filename <- paste("multiplot",selectIdentifier(object),"_",selectParameter(object),sep="")
            }
            
            multiplot(object=selectCoords(object=object,num=num,format="data.frame"),
                         contrast=selectContrast(object=object,num=num,coords=FALSE),
                         breaks=breaks,type.breaks=type.breaks,palette=palette,col=col,pch=pch,cex=cex,
                         col.NA=col.NA,pch.NA=pch.NA,xlim=xlim,ylim=ylim,axes=axes,
                         window=window,legend=legend,mfrow=mfrow,mar=mar,mgp=mgp,pty=pty,asp=asp,bg=bg,
                         xlab=xlab,ylab=ylab,main=main,num.main=num.main,cex.main=cex.main,
                         quantiles.legend=quantiles.legend,digit.legend=digit.legend,cex.legend=cex.legend,mar.legend=mar.legend,main.legend=selectParameter(object),    
                         filename=filename,width=width,height=height,path=path,unit=unit,res=res)
          }
)

#### init. ####

setMethod(f ="initNum",
          signature ="Carto3D",
          definition = function(object,num,test=TRUE,init=TRUE)
          {
            if(init==TRUE){
              if(is.null(num))
              {num <- seq(1,object@voxelDim$k)}
            }
            
            if(test==TRUE){
              if(any(num %in% seq(1,object@voxelDim$k) == FALSE) || length(unique(num))!=length(num))
              {stop("initNum[Carto3D] : wrong specification of \'num\' \n",
                    "valid values : ",paste(seq(1,object@voxelDim$k),collapse=" "),"\n",
                    "requested values : ",paste(num,collpase=" "),"\n")}
            }
            
            return(num=num)
          }
)

