#### 1- Construction functions for the objects Carto3D and MRIaggr ####

constCarto3D <- function(array,identifier,param,default_value=NULL,
                         pos_default_value=c(1,1,1),voxelDim=NULL,rm.array=FALSE){
  
  nom.array <- as.character(substitute(array))
  
  
  #### conversion du format ####
  if(class(array) == "anlz"){ # analyse object
    if(is.null(voxelDim)){
      voxelDim <- data.frame(i=NA,j=NA,k=NA,unit=array@vox_units, stringsAsFactors = FALSE)
      voxelDim[,c("i","j","k")] <- array@pixdim[2:4] 
    }
    array <- array@.Data
  }else if(class(array) %in% c("nifti","niftiExtension","niftiAuditTrail")){ # nifti object
    if(is.null(voxelDim)){
      voxelDim <- data.frame(i=NA,j=NA,k=NA,unit=oro.nifti::convert.units(oro.nifti::xyzt2space(array@xyzt_units)), stringsAsFactors = FALSE)
      voxelDim[,c("i","j","k")] <- array@pixdim[2:4]    
    }
    array <- array@.Data
  }else if(is.list(array) && length(array)==2 && "img" %in% names(array)){ # dicom object
    array <- array$img
    voxelDim <- data.frame(i=NA,j=NA,k=NA,unit=NA, stringsAsFactors = FALSE)
  }else{ # array object
    if(class(array) %in% c("matrix","array") == FALSE){
      stop("constCarto3D : wrong specification of \'array\' \n",
           "\'array\' must be either an \"anlz\", \"nifti\" or an array object \n",
           "class(array) : ",class(array),"\n")
      
    }
    voxelDim <- data.frame(i=NA,j=NA,k=NA,unit=NA, stringsAsFactors = FALSE)
  }
  
  if(length(dim(array))==4 && dim(array)[4]==1){
    array <- array[,,,1,drop=TRUE]
  }
  if(length(dim(array))==2){
    array <- array(array,dim=c(dim(array),1))
  }
  
  if(length(dim(array))!=3){
    stop("constCarto3D : wrong specification of \'array\' \n",
         "\'array\' must have 3 dimensions \n",
         "proposed \'array\' has ",length(dim(array))," dimensions : ",paste(dim(array),collapse=" ")," (parameter : ",param,") \n")
  }
  
  ##### valeurs par defaut ####
  if(is.null(default_value)){
    
    if(is.vector(pos_default_value)==TRUE){pos_default_value <- matrix(pos_default_value,nrow=1)}
    
    if(ncol(pos_default_value)!=length(dim(array))){
      stop("constCarto3D : \'pos_default_value\' and \'array\' dimensions are inconsistent \n",
           "\'array\' has ",length(dim(array))," dimensions \n",
           "\'pos_default_value\' has : ",ncol(pos_default_value)," dimensions \n")          
    }
    
    default_value <- names(which.max(table(array[pos_default_value]))[1,drop=FALSE])
  }
  
  ##### construction de la carto3D ####  
  res <- new(Class="Carto3D",
             identifier = identifier,
             parameter = param,
             voxelDim=voxelDim,
             default_value = as.character(default_value),
             contrast = array)
  
  #### nettoyage ####
  if(rm.array){
    remove(list=nom.array,envir=globalenv())
  }
  
  #### export ####
  return(res)
  
}

constMRIaggr <- function(ls.array,identifier,param,default_value=NULL,
                         pos_default_value=c(1,1,1),
                         tol=10^{-10},voxelDim=NULL,
                         trace=TRUE,rm.ls.array=FALSE){
  
  nom.ls.array <- as.character(substitute(ls.array))
  
  if(class(ls.array) %in% c("anlz","nifti","niftiExtension","niftiAuditTrail") || (is.list(ls.array) && length(ls.array)==2 && "img" %in% names(ls.array)) ){
    ls.array <- list(ls.array)
  }
  
  if(is.list(ls.array)==FALSE){
    ls.array <- list(ls.array)
  }
  M <- length(ls.array)
  
  if(is.null(param)){
    param <- names(ls.array)
  }
  
  #### test ####
  test.class <- unlist(lapply(ls.array,function(x){class(x) %in% c("anlz","nifti","niftiExtension","niftiAuditTrail","list","array","matrix")}))
  if(any(test.class==FALSE)){
    stop("constMRIaggr : wrong specification of \'ls.array\' \n",
         "\'ls.array\' must be a list of \"anlz\" or \"nifti\" or \"list\" or \"array\" \n",
         "elements of \'ls.array\' that are not of type array : ",paste(which(test.class==FALSE),collapse=" ")," \n")
  }
  
  if(length(identifier)>1){
    stop("constMRIaggr :  wrong specification of \'identifier\' \n",
         "\'identifier\' should have length 1 \n",
         "length of \'identifier\' : ",length(identifier)," \n")
  }
  
  if(length(param)!=M){
    stop("constMRIaggr : \'param\' must have the same length as \'ls.array\' \n",
         "length of \'ls.array\' : ",M," \n",
         "lenght of \'param\' : ",length(param)," \n")
  }
  
  if(!is.null(default_value) && length(default_value)!=M){
    stop("constMRIaggr :  \'default_value\' must have the same length as \'ls.array\' \n",
         "length of \'ls.array\' : ",M," \n",
         "lenght of \'default_value\' : ",length(default_value)," \n")
  }
  
  #### transformation en carto3D
  ls.Carto3D <- list()
  for(iter_m in 1:M){
    
    if(is.null(default_value)){
      ls.Carto3D[[iter_m]] <- constCarto3D(ls.array[[iter_m]],
                                           identifier=identifier,param=param[iter_m],default_value=NULL,
                                           pos_default_value=pos_default_value,voxelDim=voxelDim)
      
    }else{
      ls.Carto3D[[iter_m]] <- constCarto3D(ls.array[[iter_m]],
                                           identifier=identifier,param=param[iter_m],default_value=default_value[iter_m],
                                           pos_default_value=NULL,voxelDim=voxelDim)      
    }
    
  }
  
  #### transformation en MRIaggr
  MRIaggr <- Carto3D2MRIaggr(ls.Carto3D=ls.Carto3D,rm.Carto3D=FALSE,tol=tol,
                                      trace=trace)
  
  #### nettoyage 
  if(rm.ls.array){
    rm(list=nom.ls.array,envir=globalenv())
  }
  
  #### export
  return(MRIaggr)
}  

constLatex <- function(dir,filename=NULL,identifier=NULL,param=NULL,table=NULL,extra_text=NULL,
                        subsection=NULL,label=NULL,
                        width=0.9,trim=c(0,0,0,0),plotPerPage=3,
                        width.legend=0.35,trim.legend=c(0,0,0,0),
                        title="",date="",author="",trace=TRUE){
  
  ls.text <- list()
  
  #### preparation ####
  
  ## dossiers ou sont present les graphiques a afficher
  res_init <- initDirPat_constLatex(dir=dir,param=param,identifier=identifier,trace=trace)
  
  dirs_plot <- res_init$dirs_plot
  n.dirs_plot <- res_init$n.dirs_plot
  names_dirs <- res_init$names_dirs
  param <- res_init$param
  identifier <- res_init$identifier
  n.identifier <- res_init$n.identifier

  ## gestion des sections, subsections, subsubsection, label
  res_init <- initSection_constLatex(dirs_plot=dirs_plot,names_dirs=names_dirs,param=param,
                                     subsection=subsection,label=label)
 
  subsection <- res_init$subsection
  index_subsection <- res_init$index_subsection
  label <- res_init$label
  param <- res_init$param
   
  ## initialisation des parametres latex
  if(!is.list(trim)){
    trim <- lapply(1:100,function(x){trim})
  }
  if(length(width)==1){
    width <- rep(width,100)
  }
  
  ## initialisation de donnes cliniques
  if(!is.null(table)){
    
    if(class(table) != "list"){table <- list(table)}
    n.table <- length(table)

    ## conversion de vector en data.frame (si necessaire)
    table <- lapply(table,function(x){  
      if(is.null(dim(x))){ 
        res <- data.frame(matrix(NA,ncol=length(x),nrow=1),stringsAsFactors=FALSE)
        res[1,] <- x ; names(res) <- names(x) ; res
      }else{x}}
    )
    
    ## conversion des factors en character    
    table <- lapply(table,function(x){ 
      df.tempo <- data.frame(rapply(x, as.character, classes="factor", how="replace"),stringsAsFactors=FALSE)
      names(df.tempo) <- names(x)
      return(df.tempo)
    })
    
  }
  
  ## initialisation du texte a afficher  
  if(!is.null(extra_text) && length(extra_text) != length(identifier)){
    stop("constLatex : \'extra_text\' and \'identifier\' must have same length \n",
         " length(identifier)  : ",length(identifier)," \n",
         " length(extra_text) : ",length(extra_text)," \n") 
  }
  
  ## display 
  if(trace==TRUE){
    space.index_subsection <- initSpace(c("num",index_subsection))
    space.subsection <- initSpace(c("subsection",subsection))
    
    cat("%%%% num",space.index_subsection[1]," : subsection",space.subsection[1]," | directory \n",sep="")
    cat("% (label) \n\n",sep="")
    
    for(iter_subsection in 1:length(subsection)){
      cat("%%%% ",iter_subsection,space.index_subsection[iter_subsection+1]," : ",subsection[iter_subsection],sep="")
    
      cat(paste(space.subsection[iter_subsection+1]," | /",names_dirs[index_subsection==iter_subsection],"\n",
                "% label : \"",gsub("\n","",label[[iter_subsection]],fixed=TRUE),"\" \n\n",
                sep=""),sep="")      
    }
    
    cat("% graph of patient ")       
  }
  
  #### preambule ####
  text.preamble <- c("%endTrace \n \n \\documentclass[a4paper]{article} \n \n \n",
                     "%%%%%%%%%%%%%%%%% Preamble %%%%%%%%%%%%%%%%% \n",
                     "\\usepackage[utf8]{inputenc} \n",
                     "\\usepackage{amssymb} \n",
                     "\\usepackage{amsmath} \n",
                     "\\usepackage{titlesec} \n",
                     "\\usepackage{geometry} \n",
                     "\\usepackage{enumitem} \n",
                     "\\usepackage{graphicx} \n",
                     "\\usepackage{color} \n",
                     "\\usepackage{xspace} \n",
                     "\\usepackage{hyperref} \n",
                     "\\usepackage{caption} \n \n \n",
                     "%%%% Margin %%%% \n",
                     "\\geometry{\n",
                     "a4paper,\n",
                     "total={210mm,297mm},\n",
                     "left=10mm,\n",
                     "right=10mm,\n",
                     "top=5mm,\n",
                     "bottom=12.5mm,\n",
                     "}\n",
                     "\\setlength{\\textfloatsep}{5pt} % space between last top float or first bottom float and the text \n",
                     "\\setlength{\\intextsep}{5pt} % space left on top and bottom of an in-text float \n",
                     "\\setlength{\\abovecaptionskip}{10.0pt} % space above caption \n",
                     "\\setlength{\\belowcaptionskip}{0pt} % space below caption \n",
                     "\\titlespacing\\section{0pt}{0ex}{5ex} % {left space}{above space}{below space} \n",
                     "\\titlespacing\\subsection{25pt}{3ex}{1ex} % {left space}{above space}{below space} \n \n \n",                     
                     "%%%% Margin %%%% \n",
                     paste("\\graphicspath{{./",tail(strsplit(dir,split="/",fixed=TRUE)[[1]],1),"/}} \n \n \n",sep=""),
                     "%%%% Title %%%% \n",
                     paste("\\title{",title,"} \n",sep=""),
                     paste("\\date{",date,"} \n",sep=""),
                     paste("\\author{",author,"} \n \n \n",sep=""),
                     "%%%%%%%%%%%%%%%%% End Preamble %%%%%%%%%%%%%%%%% \n \n \n"
  )
  
  #### begin doc ####
  text.begin <- c("%","\n",
                  "\\begin{document} \n \n ",
                  "\\maketitle \n \n ",
                  "\\setcounter{tocdepth}{1} \n",
                  "\\tableofcontents \n",
                  "\\setcounter{tocdepth}{3} \n",
                  "\\clearpage \n"
  ) 
 
  #### boucle sur les patients ####
  for(iter_pat in 1:n.identifier){
    
    if(trace==TRUE){cat(iter_pat,"(",identifier[iter_pat],") ",sep="")}
    
    #### definition des graphiques a afficher pour chaque patient 
    res_init <- initPlot_constLatex(dir=dir,names_dirs=names_dirs,dirs_plot=dirs_plot,table=table,
                                     identifier=identifier[iter_pat],plotPerPage=plotPerPage)
     
    ls.plot <- res_init$ls.plot
    ls.newplot <- res_init$ls.newplot
    ls.endplot <- res_init$ls.endplot
    ls.legend <- res_init$ls.legend
    ls.slices <- res_init$ls.slices
    ls.table <- res_init$ls.table
   
     #### affichage clinique ####
    ls.text[[iter_pat]] <- paste("%","\n \n \\section",paste("{patient ",identifier[iter_pat],"} \n \n ",sep=""),sep="")
   
    if(length(ls.table)>0){ # donnees cliniques disponibles pour le patient
      
      for(iter_table in 1:length(ls.table)){
        
        if(ncol(ls.table[[iter_table]])==0){next}
        
        ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],
                                 "\\begin{table}[!h] \n",
                                 "\\centering \n",
                                 paste("\\begin{tabular}{|",paste(rep("c",ncol(ls.table[[iter_table]])),collapse=""),"|} \n",sep=""),
                                 paste(paste(names(ls.table[[iter_table]]),collapse=" & ")," \\\\ \\hline \n",sep=""))
        
        for(iter_rowtable in 1:nrow(ls.table[[iter_table]])){
          ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],
                                   paste(paste(ls.table[[iter_table]][iter_rowtable,],collapse=" & ")," \\\\  \n",sep=""))
        }
        
        ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],
                                 "\\end{tabular} \n",
                                 "\\end{table} \n \n",
                                 "\\smallskip  \n \n")
                                 
      }
    }
    
    #### extra text ####
    if(!is.null(extra_text)){
      ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],extra_text[[iter_pat]]) 
    }
    
    if(length(ls.table)>0 || !is.null(extra_text)){ # si donnees cliniques 
      ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],'\n \n \\clearpage \n \n')    
    }else{ # sinon si la page est pleine
      if(unlist(lapply(ls.plot,function(x){ if(all(is.na(x))){NULL}else{length(x)}} ))[1]>=plotPerPage){
      ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],'\n \n \\vspace{-5ex} \n \n')
      }
    }
    
    #### affichage des images ####
    count.subsection <- 0
    
    for(iter_subsection in 1:length(subsection)){
      
      ## pas de graphique on saute la subsection
      if(all(is.na(ls.plot[[iter_subsection]]))){next}else{count.subsection <- count.subsection + 1} 
      
      ## debut de la subsection
      ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],
                               paste("\\subsection",paste("{",subsection[iter_subsection],"} \n  ",sep=""),sep="")
      )
            
      ## initialisation
      iter_subsubsection <- 0
      n.plot_tempo <- length(ls.plot[[iter_subsection]])
      count_plot <- 1
      count_figure <- 0
      
#       if(length(ls.table)==0 && is.null(extra_text) && count.subsection==1 && n.plot_tempo>=plotPerPage){ # si section sur la meme page que subsection et premiere subsection et la figure arrive en bas
#         ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],'\n \n \\vspace{-3ex} \n \n')
#       }
      
      for(iter_plot in 1:n.plot_tempo){ #### gestion des graphiques
        
        ## debut d une nouvelle page
        if(ls.newplot[[iter_subsection]][iter_plot]==1){
                   
          ## debut d une nouvelle figure
          count_figure <- count_figure + 1 
          ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],
                                   "\\begin{figure}[!h] \n",
                                   "\\centering \n"
          ) 
          
        } 
        
        ## parametres latex entre legend et image
        if(ls.legend[[iter_subsection]][iter_plot]==TRUE){
          trim_tempo <- trim.legend
          width_tempo <- width.legend
        }else{
          trim_tempo <- trim[[count_plot]]
          width_tempo <- width[[count_plot]]
        }
        
        ## affichage de l image 
        ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],
                                 paste("\\includegraphics[trim= ",trim_tempo[1],"mm ",trim_tempo[2],"mm ",trim_tempo[1],"mm ",trim_tempo[1],"mm,clip,width=",width_tempo,"\\textwidth]",sep=""),
                                 paste("{",ls.plot[[iter_subsection]][iter_plot],"} \n",sep="")
        )
        
        ## affichage de la legende
        if(ls.endplot[[iter_subsection]][iter_plot]==1){
         
          ls.text[[iter_pat]] <- c(ls.text[[iter_pat]],
                                   paste("\\caption[foo bar]{Patient ",identifier[iter_pat]," ",ls.slices[[iter_subsection]][count_figure]," - ",label[[iter_subsection]],"} \n",sep=""),
                                   #paste("\\label{fig:Sweave_",identifier[iter_pat],"_",param[iter_subsection],"} \n",sep=""),
                                   "\\vspace{-1cm} \n \n",
                                   "\\end{figure} \n \n",                                   
                                   '\\clearpage \n \n'
          )          
        }
        count_plot <- count_plot + 1         
        
        
      } # iteration plot
    } # iteration subsection
  } # iteration patient
  if(trace==TRUE){cat("\n")}
  
  #### end doc ####
  text.end <- c("%","\n \n \\end{document}  ")    
  
  #### create Latex file ####
  if(!is.null(filename)){
     
    file.create(paste(filename,"tex",sep="."))
    
    sink(paste(filename,"tex",sep="."))
    cat(text.preamble,sep="")
    cat(text.begin,sep="") 
    for(iter_list in 1:length(ls.text)){
       cat(ls.text[[iter_list]],sep="")
    }
    cat(text.end,sep="")
    sink()
    
  }
  
  #### export ####
  
  return(list(text.preamble=text.preamble,
              text.begin=text.begin,
              ls.text=ls.text,
              text.end=text.end)
  )
}

#### 2- Conversion functions ####
array2df <- function(array,coords=NULL,
                     name_newparam="res",names_coords=letters[9:(8+ncol(coords))],na.rm=TRUE){
  
  ### preparation 
  p <- length(dim(array))  
  n <- length(array)  
  index_array <- arrayInd(1:n,.dim=dim(array))
  if(na.rm){
    index_array <- index_array[is.na(array)==FALSE,]
  }
  
  if(is.null(coords)){
    coords <- index_array
  }else{
    if(any(dim(coords) != dim(index_array))){
      stop("array2df : incorrect dimension for \'coords\' : \n",
           "lenght(array[!is.na(array)]),length(dim(array)) =",paste(dim(index_array),collapse=" "),"\n",
           "dim(coords) = ",paste(dim(coords),collapse=" "),"\n")
    }
  }
  
  if(length(names_coords) != p){
    stop("array2df : wrong specification of \'names_coords\' : \n",
         "required number of names : ",p,"\n",
         "proposed number of names : ",length(names_coords),"\n")
  }
  
  ### integration des donnees 
  data <- data.frame(coords,array[index_array])
  names(data) <- c(names_coords,name_newparam)
  
  ### export
  return(data)
}

Carto3D2MRIaggr <- function(ls.Carto3D,rm.Carto3D=FALSE,tol=10^{-10},
                            num=NULL,trace=TRUE){
  
  nom.ls.Carto3D <-  as.character(substitute(ls.Carto3D))[-1]
  
  if(!is.list(ls.Carto3D)){
    ls.Carto3D <- as.list(ls.Carto3D)
  }
  
  #### test preliminaire ####
  if(any(unlist(lapply(ls.Carto3D,function(x){class(x)=="Carto3D"}))==FALSE)){
    stop("Carto3D2MRIaggr : wrong specificition of \'ls.array\' \n",
         "\'ls.array\' must be a list of \"Carto3D\" objects \n",
         "elements of \'ls.array\' that are not \"Carto3D\" objects : ",paste(which(unlist(lapply(ls.Carto3D,function(x){class(x)=="Carto3D"}))==FALSE),collapse=" ")," \n")    
  }
  
  #### preparation ####
  
  if(is.null(num)){
    num <- seq(1,selectFieldDim(ls.Carto3D[[1]])$k,by=1)
  }
  fieldDim <- data.frame(i=selectFieldDim(ls.Carto3D[[1]])$i,
                         j=selectFieldDim(ls.Carto3D[[1]])$j,
                         k=length(num), stringsAsFactors = FALSE)
  
  default_value <- data.frame(matrix("NA",ncol=length(ls.Carto3D),nrow=1),stringsAsFactors=FALSE)
  identifiant <- selectIdentifier(ls.Carto3D[[1]])
  size <- ls.Carto3D[[1]]@voxelDim
  
  data_global <- data.frame(matrix(NA,
                                   nrow=selectN(ls.Carto3D[[1]],num=num),
                                   ncol=length(ls.Carto3D)))
  data_global[,1:3] <- selectCoords(ls.Carto3D[[1]],num=num,format="data.frame")
  names(data_global)[1:3] <- c("i","j","k")
  
  #### assemblage ####
  if(trace){cat("Merging : ")}
  for(iter in 1:length(ls.Carto3D)){
    
    if(trace){cat("(",iter,")",sep="")}
    
    data_tempo <- selectContrast(ls.Carto3D[[iter]],num=num,coords=TRUE)
    nom_param <- selectParameter(ls.Carto3D[[iter]])
    if(trace){cat(" ",nom_param," ",sep="")}
    
    default_value[iter] <- selectDefault_value(ls.Carto3D[[iter]])
    names(default_value)[iter] <- nom_param
    
    if(sum(abs(data_global[,1:3]-data_tempo[,1:3]))>tol)
    {stop("Carto3D2MRIaggr - some \"Carto3D\" objects of \'ls.Carto3D\' have different coordinates \n",
          "Difference between \'ls.Carto3D\'[[1]] and \'ls.Carto3D\'[[",iter,"]]  \n")}
    
    if(selectIdentifier(ls.Carto3D[[iter]])!=identifiant)
    {stop("Carto3D2MRIaggr -  some \"Carto3D\" objects of \'ls.Carto3D objects correspond to different patients \n",
          "patient of \'ls.Carto3D\'[[1]] : ",selectIdentifier(ls.Carto3D[[1]]),"\n",
          "patient of \'ls.Carto3D\'[[",iter,"]] : ",selectIdentifier(ls.Carto3D[[iter]]),"\n")}
    
    if(nom_param %in% names(data_global))
    {warning("Carto3D2MRIaggr - some \"Carto3D\" objects of \'ls.Carto3D\' correspond to the same parameter \n",
             "common parameter : ",nom_param,"\n",
             "index of the elements in \'ls.Carto3D\' : ",paste(which(names(data_global) == nom_param)-2,collapse=" "),") \n")}
    
    data_global[,iter+3] <- data_tempo[,4]
    names(data_global)[iter+3] <- nom_param
  } 
  if(trace){cat("\n")}
  
  #### nettoyage
  if(rm.Carto3D){
    rm(list=nom.ls.Carto3D,envir=globalenv())
  }
  
  #### export ####
  res <- new(Class="MRIaggr",
             identifier=identifiant,
             contrast=data_global,
             default_value=default_value,
             fieldDim=fieldDim,
             voxelDim=size)
  
  return(res)
}

df2array <- function(contrast,coords,format="any",default_value=NA,
                     range.coords=NULL){
  
  ### preparation 
  if(class(coords) %in% c("matrix","data.frame")==FALSE)
  {stop("df2array : \'coords\' must be of class \"matrix\" or \"data.frame\"  \n",
        "type of \'coords\' : ",class(coords),"\n")}
  
  contrast <- as.data.frame(contrast)
  M <- ncol(contrast)
  p <- ncol(coords)
  n <- nrow(coords)
  
  if(n != nrow(contrast)){
    stop("df2array : \'coords\' and \'contrast\' have inconsistent dimensions \n",
         "number of rows of \'coords\' : ",n,"\n",
         "number of rows of \'contrast\' : ",nrow(contrast),"\n")}
  
  if(p>3){
    stop("df2array : wrong specification of \'coords\' \n",
         "\'coords\' can have up to 3 columns \n",
         "number of columns of \'coords\' : ",p,"\n")
  }
  
  if( format %in% c("any","matrix","data.frame","list") == FALSE )
  {  stop("df2array : wrong specification of \'format\' \n",
          "valid formats  :  \"any\" \"matrix\" \"data.frame\" \"list\" \n",
          "requested format : ",format,"\n")
  }  
  
  if(!is.null(range.coords) && p != length(range.coords)){
    stop("df2array : wrong specification of \'range.coords\' \n",
         "if not NULL, \'range.coords\' must be a vector of length ",p," \n",
         "is(range.coords) : ",paste(is(range.coords),collapse=" ")," \n",
         "length(range.coords) : ",length(range.coords),"\n")
  }  
  
  #### definition des coordonnees des points dans le repere de la matrice
  scale <- rep(0,p)
  if(is.null(range.coords)){ 
    Mcoords <- apply(coords,2,
                     function(x){scale_tempo = min(x)-1;
                                 return(c(scale_tempo,x-scale_tempo))}
    )
    scale <- Mcoords[1,]
    Mcoords <- Mcoords[-1,,drop=FALSE]    
  }else{
    Mcoords <- coords
  }
  
  if(any(Mcoords %% 1 !=0)){    
    Mcoords <- apply(Mcoords,2,function(x){as.numeric(as.factor(x))})
  }
  
  dimnames <- names(Mcoords)
  Mcoords <- as.matrix(Mcoords)
  if(is.null(range.coords)){
    range.coords <- unlist(apply(Mcoords,2,function(x){max(x)}))
  }
  if(any(range.coords < apply(Mcoords,2,max))){
    stop("df2array : wrong specification of \'range.coords\' \n",
         "\'range.coords\' must be at least ",paste(apply(Mcoords,2,max),collapse=" ")," \n",
         "requested \'range.coords\' : ",paste(range.coords,collapse=" "),"\n")    
  }
 
  #### index correspondant aux points dans la matrice
  
  Mindex <- eval(parse(text=paste("1+",paste( c(1,if(ncol(Mcoords)>1){range.coords[1]},if(ncol(Mcoords)>2){range.coords[1]*range.coords[2]}),
                                              "*(Mcoords[,",1:ncol(Mcoords),"]-1)",sep="",collapse="+"),sep="")))
  
  ### integration des donnees
  dataM <- list()
  
  for(iter_m in 1:M){
    
    dataM[[iter_m]] <- array(default_value,dim=range.coords,dimnames=dimnames)
    dataM[[iter_m]][Mindex] <- contrast[,iter_m]
    
    if(format=="matrix" && p==2){
      dataM[[iter_m]] <- as.matrix(dataM[[iter_m]])
    }
    if(format=="data.frame" && p==2){
      dataM[[iter_m]] <- as.data.frame(dataM[[iter_m]])
    }
  }
  names(dataM) <- names(contrast)
  
  ### mise en forme
  if(format %in% c("matrix","data.frame")  && M==1){
    dataM <- dataM[[1]]
  }
  
  unique_coords <- lapply(1:ncol(Mcoords),function(x){scale[x]+seq(min(Mcoords[,x]),max(Mcoords[,x]),by=1)})
  names(unique_coords) <- names(coords)
  
  return(list(contrast=dataM,
              coords=coords,
              unique_coords=unique_coords))
  
  
  
}

readMRI <- function (file, format, na.value=0,
                     what = "numeric", size = "NA_integer_", dim = NULL,
                     SPM=FALSE, reorient = FALSE, flipud=FALSE){
  
  if(format %in% c("rawb.gz","analyze","nifti","dicom") == FALSE){
    stop("readMRI : incorrect \'format\' : \n",
         "available format : \"rawb.gz\" \"analyze\" \"nifti\" \"dicom\" \n",
         "proposed format : ",format,"\n")
  }
  
  if (format == "rawb.gz") {
    if (is.null(dim) || !is.vector(dim) || length(dim) != 3){
        stop("readMRI : incorrect \'dim\' \n",
             "dim must be a vector of length 3 \n",
             "proposed dim : ",paste(dim,collapse=" "),"\n")
    }
    f <- gzfile(file, open = "rb")
    on.exit(close(f))
    data <- readBin(con=f, what=what, n=prod(dim), size = size)
    res <- array(data, dim)
  }
  else if (format == "analyze") {
    test.package <- requireNamespace("oro.nifti",quietly=TRUE)
    if(test.package==FALSE){
      stop("readMRI : this function with argument format=\"oro.nifti\" requires to have installed the oro.nifti package to work \n")
    }
    
    res <- oro.nifti::readANALYZE(file,SPM=SPM)
  }
  else if (format == "nifti") {
    test.package <- requireNamespace("oro.nifti",quietly=TRUE)
    if(test.package==FALSE){
      stop("readMRI : this function with argument format=\"oro.nifti\" requires to have installed the oro.nifti package to work \n")
    }
    
    res <- oro.nifti::readNIfTI(file,reorient=reorient)
  }
  else if (format == "dicom")  {
    test.package <- requireNamespace("oro.dicom",quietly=TRUE)
    if(test.package==FALSE){
      stop("readMRI : this function with argument format=\"oro.dicom\" requires to have installed the oro.dicom package to work \n")
    }
    
    res <- oro.dicom::readDICOMFile(file,flipud=flipud)
  }
  
  if(!is.na(na.value)){
    
    if(format=="dicom"){
      test.na <- is.na(res$img)
      if(any(test.na)){res$img[test.na] <- na.value}        
    }else{
      test.na <- is.na(res)
      if(any(test.na)){res[test.na] <- na.value}
    }
    
  }
  
  return(res)
  
}

writeMRI <- function (data, file, format, gzipped=TRUE, verbose=FALSE, size = "NA_integer_"){
  
  if(format %in% c("rawb.gz","analyze","nifti","dicom") == FALSE){
    stop("writeMRI : incorrect \'format\' : \n",
         "available format : \"rawb.gz\" \"analyze\" \"nifti\" \"dicom\" \n",
         "proposed format : ",format,"\n")
  }
  
  objClass <- class(data)
  if (objClass != "array" || length(dim(data)) != 3){
    stop("writeMRI : incorrect \'data\' : \n",
         "data has to be a 3 dimensional array  \n",
         "proposed data (dimension) : ",paste(is(data),collapse=" ")," (",paste(dim(data),collapse=" "),") \n")
  }
    
  if (format == "rawb.gz") {
    f <- gzfile(file, open = "wb")
    writeBin(con=f, object=as.vector(data), size = size)
    on.exit(close(f))
  }
  else if (format == "analyze") {
    test.package <- requireNamespace("oro.nifti",quietly=TRUE)
    if(test.package==FALSE){
      stop("writeMRI : this function with argument format=\"oro.nifti\" requires to have installed the oro.nifti package to work \n")
    }
    
    oro.nifti::writeANALYZE(as(data, "anlz"),filename=file,gzipped=gzipped,verbose=verbose)
  }
  else if (format == "nifti") {
    test.package <- requireNamespace("oro.nifti",quietly=TRUE)
    if(test.package==FALSE){
      stop("writeMRI : this function with argument format=\"oro.nifti\" requires to have installed the oro.nifti package to work \n")
    }
    
    oro.nifti::writeNIfTI(as(data, "nifti"),filename=file,gzipped=gzipped,verbose=verbose)
  }
  else if (format == "dicom")  {
    test.package <- requireNamespace("oro.dicom",quietly=TRUE)
    if(test.package==FALSE){
      stop("writeMRI : this function with argument format=\"oro.dicom\" requires to have installed the oro.dicom package to work \n")
    }
    
    cat("writeMRI for dicom files is not implemented \n")
    invisible(return(FALSE))    
  }

}

#### 3- Calculation functions ####

calcAUPRC <- function(x,y,subdivisions=10000,performance=NULL){
 
  test.package <- requireNamespace("ROCR",quietly=TRUE)
  if(test.package==FALSE){
    stop("calcAUPRC : this function requires to have installed the ROCR package to work \n")
  }
  
  if(is.null(performance)){    
    performance <- ROCR::performance(ROCR::prediction(x,y),x.measure="rec",measure="prec")
  }else{
    if(class(performance)!="performance"){
      stop("calcAUPRC : wrong specification of \'performance\' \n",
           "\'performance\' must an object \"prediction\" \n",
           "is(performance) : ",paste(is(performance),collapse=" "),"\n")   
    }
    
    if( any(c("Precision","Recall") %in% c(performance@x.name,performance@y.name) == FALSE) ){
      stop("calcAUPRC : wrong specification of \'performance\' \n",
           "\'performance\' must contains \"Precision\" and \"Recall\" performance measures \n",
           "measures in \'performance\' : \"",performance@x.name,"\" \"",performance@y.name,"\" \n")   
    }
  }
  
  M.PR <- cbind(performance@x.values[[1]], 
                performance@y.values[[1]])
  colnames(M.PR) <- c(performance@x.name,performance@y.name)
  M.PR <- M.PR[c(-1,-nrow(M.PR)),,drop=FALSE]
  M.PR <- cbind(Recall=unique(M.PR[,"Recall"]),
                Precision=tapply(M.PR[,"Precision"],M.PR[,"Recall"],function(x){max(x)}))
  if(nrow(M.PR)==1){return(M.PR[,"Precision"])} # cas avec un seuil 
  
  f  <- approxfun(x=M.PR[,"Recall"],y=M.PR[,"Precision"])
  f.int <- integrate(f, min(M.PR[,"Recall"]), max(M.PR[,"Recall"]),subdivisions=subdivisions)$value
  return(f.int/(max(M.PR[,"Recall"])-min(M.PR[,"Recall"])))
}

calcGroupsCoords <- function(coords,array=NULL,Neighborhood,max_groups=10000,trace=TRUE){
  
  if(is.null(array)){
    if(any(coords<0)){
      stop("calcGroupsCoords : wrong specification of  \'coords\' \n",
           "\'coords must\' be positve \n",
           "min(coords) : ",min(coords),"\n")
    }
    
    if(nrow(coords)==0){
      return(list(ls.group=NULL,
                  df.group=NULL,
                  group_size=0)
      )
    }
    
    if(nrow(coords)==1){
      return(list(ls.group=list(1),
                  df.group=data.frame(coords,index=1,group=1),
                  group_size=1)
      )
    }
    
    coords_NNA <- apply(coords,2,function(x){x-min(x)})
    dim_coordsNNA <- apply(coords_NNA,2,max)+1
    p <- length(dim_coordsNNA)
    index_NNA <- rowSums(sweep(coords_NNA,2,c(1,cumprod(dim_coordsNNA)[-p]),FUN="*"))
  }else{
    coords_NNA <- which(!is.na(array),arr.ind=TRUE)-1
    
    if(nrow(coords_NNA)==0){
      return(list(ls.group=NULL,
                  df.group=NULL,
                  group_size=0)
      )
    }
    
    if(nrow(coords_NNA)==1){
      return(list(ls.group=list(1),
                  df.group=data.frame(coords_NNA,index=1,group=1),
                  group_size=1)
      )
    }
    
    dim_coordsNNA <- dim(array)
    p <- length(dim_coordsNNA)
    index_NNA <- which(!is.na(array))-1
  }
  
  #### definition de Neighborhood
  if(length(Neighborhood)==1 && is.character(Neighborhood)){
    Neighborhood <- initNeighborhood(Neighborhood,method="calcGroupsCoords")
  }
  
  if(ncol(Neighborhood)!=p){
    stop("calcGroupsCoords : wrong specification of  \'Neighborhood\' \n",
         "required number of columns : ",p,"\n",
         "proposed number of columns : ",ncol(Neighborhood),"\n")
  }
  
  #### Rcpp
  res_cpp <-  calcGroupsCoords_cpp(coords_NNA=coords_NNA,
                                            index_NNA=index_NNA,
                                            Neighborhood=Neighborhood,
                                            coords_max=dim_coordsNNA,
                                            max_groups=max_groups,
                                            trace=trace)
  
  if(length(res_cpp$group_size)==max_groups && sum(res_cpp$group_size)<nrow(res_cpp$group)){
    warning("calcGroupsCoords : maximum number of groups reached \n",
            "the identification of the spatial groups may not be complet \n",
            "set \'max_groups\' higher to allow more groups \n")
  }
  
  #### export  
  ls.group <- list()
  for(iter_group in 1:length(res_cpp$group_size)){
    ls.group[[iter_group]] <- res_cpp$group[res_cpp$group[,2]==iter_group,1]+1
  }
  
  if(is.null(array)){
    df.group <- as.data.frame(cbind(coords[res_cpp$group[,1]+1,],res_cpp$group))
    names(df.group) <- c(letters[8+1:p],"index","group")
    df.group$index <- df.group$index+1
  }else{
    df.group <- NULL
  }
  
  return(list(ls.group=ls.group,
              df.group=df.group,
              group_size=res_cpp$group_size)
  )
}

calcThreshold <- function(contrast,param,hemisphere=NULL,rm.CSF=FALSE,threshold=1:10,decreasing=FALSE,
                          GRalgo=FALSE,W=NULL,seed=NULL,as.logical=FALSE,trace=TRUE){
  
  p <- length(param)
  n <- nrow(contrast)
  if(decreasing==TRUE){
    threshold <- - threshold
  }
  threshold <- sort(threshold)
  n.threshold <- length(threshold)
  
  if(is.matrix(contrast) == FALSE && is.data.frame(contrast) == FALSE){
    stop("calcThreshold : wrong specification of \'contrast\' \n",
         "\'contrast\' must be a matrix or a data.frame \n",
         "is \'contrast\' : ",paste(is(contrast),collapse=" ")," \n")
  }
  
  if(any(param %in% names(contrast) == FALSE)){
    stop("calcThreshold : wrong specification of \'param\' \n",
         "\'contrast\' does not contains all param \n",
         "missing param : ",paste(param[param %in% names(contrast) == FALSE],collapse=" ")," \n")
  }
  
  
  if(is.character(rm.CSF)){
    param.CSF <- rm.CSF
    rm.CSF <- TRUE
  }else{
    param.CSF <- "CSF"
  }
  
  if(rm.CSF==TRUE && (param.CSF %in% names(contrast) == FALSE)){
    stop("calcThreshold : wrong specification of \'contrast\' \n",
         "contrast must be contains a column \"",param.CSF,"\" if rm.CSF=TRUE \n",
         "column names of \'contrast\' : ",paste(names(contrast),collapse=" ")," \n")
  }
  
  if(!is.null(hemisphere)){
    
    if("hemisphere" %in% names(contrast) == FALSE){
      stop("calcThreshold : wrong specification of \'contrast\' \n",
           "contrast must be contains a columns \"hemisphere\" if \'hemisphere\' is not NULL \n",
           "column names of \'contrast\' : ",paste(names(contrast),collapse=" ")," \n")
    }
    
    if(hemisphere %in% unique(contrast$hemisphere) == FALSE){
      stop("calcThreshold : wrong specification of \'hemisphere\' \n",
           "\'hemisphere\' does not match element of the hemisphere column in contrast \n",
           "unique(contrast$hemisphere) : ",unique(contrast$hemisphere),"\n",
           "proposed \'hemisphere\' : ",hemisphere," \n")
    }
    
  } 
  
  if(GRalgo==TRUE){
    #### test package W
	test.package <- requireNamespace("spam",quietly=TRUE)
    if(test.package==FALSE){
              stop("calcThreshold : this function requires to have installed the spam package to work \n")
    }
			
	test.package <- requireNamespace("Matrix",quietly=TRUE)
    if(test.package==FALSE){
       stop("calcThreshold : this function requires to have installed the Matrix package to work \n")
    }
			
    if(is.character(seed) && all(seed %in% names(contrast))){
      seed <- contrast[,seed,drop=FALSE] 
      
      if(as.logical==TRUE){seed <- apply(seed,2,as.logical)}
      if(any(apply(seed,2,is.logical)==FALSE)){
        stop("calcThreshold : type of \'seed\' is not logical \n",
             "proposed type : ",paste(apply(seed,2,class),collapse=" "),"\n",
             "to force the conversion to logical set \'as.logical\'= TRUE \n")
      }    
      seed <-  which(rowSums(seed)>0)
    }
    
    if(any(seed<0) || any(seed>n) || any(seed %% 1 !=0)){
      stop("calcThreshold : wrong specification of \'seed\' \n",
           "it must be intergers between 1 and ",n," indicating the seed sites \n",
           "proposed \'seed\' (range | non integer) : ",paste(range(seed),collapse=" ")," | ",sum(seed %% 1 !=0)," \n")
    }
    
    if(is.null(W) || nrow(W)!=n || ncol(W)!=n){
      stop("calcThreshold : \'W\' dimensions does not match \'contrast\' number of rows \n",
           "nrow(contrast)  : ",nrow(contrast)," \n",
           "dim(W)  : ",paste(dim(W),collapse=" ")," \n")
    }
    
    #     if(!is.null(coords)){
    #       if(nrow(coords) != n){
    #         stop("calcThreshold : \'coords\' number of rows does not match \'contrast\' number of rows \n",
    #              "nrow(contrast)  : ",nrow(contrast)," \n",
    #              "nrow(coords) : ",nrow(coords)," \n")
    #       }
    #       
    #       if(!is.numeric(sd.iter) || sd.iter<=0){
    #         stop("calcThreshold : wrong specification of \'sd.iter\' \n",
    #              "\'sd.iter\' must be a positive number \n",
    #              "proposed \'sd.iter\' : ",sd.iter," \n")
    #       }
    #     }
    
  }
  
  if(decreasing == TRUE){    
    contrast[,param] <- - contrast[,param]
  }
  
  #### step 1 - CSF et hemi
  if(trace==TRUE){cat("Step 1 : ")}
  index.perf <- 1:n
  
  if(!is.null(hemisphere)){
    if(trace==TRUE){cat("keep only the \"",hemisphere,"\" hemisphere",sep="")}
    index.perf <- intersect(index.perf,which(contrast$hemisphere %in% hemisphere))    
  } else{
    if(trace==TRUE){cat("keep both hemipheres")}
  }
  
  if(rm.CSF==TRUE){
    if(trace==TRUE){cat(", remove CSF ")}
    index.perf <- intersect(index.perf,which(contrast[,param.CSF]<0.5))
  } else{
    if(trace==TRUE){cat(", keep CSF ")}
  }
  
  contrast <- contrast[,param,drop=FALSE]
  min_perf <- threshold[1]-1
  contrast[-index.perf,] <- min_perf
  if(trace==TRUE){cat("\n")}
  
  #### etape 2 - seuillage
  if(trace==TRUE){cat("Step 2 : thresholding ")}
  tempo <- matrix(min_perf,nrow=n,ncol=p)
  colnames(tempo) <- param
  
  
  for(iter_threshold in threshold){
    if(trace==TRUE){cat("*")}
    for(iter_param in 1:p){
      tempo[contrast[,param[iter_param]]>=iter_threshold,iter_param] <- iter_threshold      
    }
  }
  
  contrast <- tempo
  rm(tempo)
  gc()
  if(trace==TRUE){cat("\n")}
  
  #### etape 3 - GR
  if(GRalgo==TRUE){
    if(trace==TRUE){cat("Step 3 : Growing Region \n",sep="")}
    
    for(iter_param in param){
      if(trace==TRUE){cat(iter_param," ",sep="")}
      
      for(iter_threshold in 1:n.threshold){
        if(trace==TRUE){cat("*")}
        
        index_threshold <- which(contrast[,iter_param]>=threshold[iter_threshold])
        
        contrastBis <- rep(0,n)
        contrastBis[index_threshold] <- 100
        
        resGR <- GRalgo(contrastBis,W=W,seed=seed[seed>=threshold[iter_threshold]],
                        sigma_max=0.0001,range=c(threshold[iter_threshold],Inf),breaks=seq(-1,109,by=10),step=10,
                        operator="sd",iter_max=1000,
                        history_sigma=FALSE,history_step=FALSE,history_front=FALSE)
        
        #         if(!is.null(coords) && iter_threshold==1){ # post treatment
        # 
        #           # identify all spatial groups
        #           resSpatialGroups <- calcGroupsW(W[resGR$GR,resGR$GR])
        #           
        #           # search all points close to the barycenter
        #           newSeed <- sapply(1:resSpatialGroups$group_number,
        #                             function(iter_group){
        #                               index_group <- which(resSpatialGroups$group==iter_group)
        #                               indexG_group <- resGR$GR[index_group]
        #                               
        #                               # each hypoperfusion area has an influence proportional to its degree
        #                               prevalence <- table(contrast[indexG_group,iter_param])/length(indexG_group)
        #                               prevalence[prevalence<0.001] <- 0.001
        #                               contrast_radius <- as.factor(contrast[indexG_group,iter_param])
        #                               levels(contrast_radius) <- which(as.character(threshold) %in%  levels(contrast_radius))/prevalence
        #                               
        #                               distBary <- calcRadius_cpp(as.matrix(coords[indexG_group,]), 
        #                                                          as.numeric(as.character(contrast_radius)), 10^{-12}, rep(T,length(index_group)), F)$distance
        #                               return(index_group[which.min(distBary)])
        #                             })
        #           
        #           
        #           # new GR from the central seeds
        #           n.GR <- length(resGR$GR)
        #           resGR.iter <- GRalgo(contrast=rep(1,n.GR),W=W[resGR$GR,resGR$GR,drop=FALSE],
        #                                seed=as.vector(unlist(newSeed)),sigma_max=1,range=c(0,2),
        #                                breaks=c(-1,0,1,2),step=1,operator="sd",iter_max=1000,
        #                                history_sigma=FALSE,history_step=TRUE,history_front=TRUE) 
        #           
        #           if(length(resGR.iter$GR)!=n.GR){
        #             stop("calcThreshold : spatial post-treatement failed \n",
        #                  "some sites were not linked with the barycenter sites of the GR \n",
        #                  "possible mispecification of \'W\' \n")
        #           }
        #           # normalize step
        #           sd_step <- sd(resGR.iter$history_step,na.rm=T)
        #           
        #           # summarize front
        #           front.split <- lapply(strsplit(resGR.iter$history_front,split=".",fixed=TRUE),
        #                                 function(x){c(as.numeric(x),rep(NA,resGR.iter$iter-length(x)))}
        #           )
        #           Mfront <- matrix(unlist(front.split),ncol=resGR.iter$iter,nrow=length(resGR.iter$GR),byrow=TRUE)
        #           
        #           index.rm  <- which(resGR.iter$history_step > sd_step*sd.iter)
        #           
        #           if(length(index.rm)>0){           
        #             resGR$GR <- resGR$GR[-index.rm]
        #           }
        #           
        #         }  
        
        
        if(length(setdiff(index_threshold,resGR$GR))==0){next}
        
        if(iter_threshold == 1){
          contrast[setdiff(index_threshold,resGR$GR),iter_param] <- min_perf
        }else{
          contrast[setdiff(index_threshold,resGR$GR),iter_param] <- threshold[iter_threshold-1]
        }
        
        
        
      }
      if(trace==TRUE){cat("\n")}
    }
    if(trace==TRUE){cat("\n")}
  }
  
  if(decreasing == TRUE){
    contrast <- - contrast  
  }
  names(contrast) <- param
  
  return(as.data.frame(contrast))
}

calcGroupsW <- function(W,subset=NULL,max_groups=10000){ 
  
  #### test package W
  test.package <- requireNamespace("spam",quietly=TRUE)
  if(test.package==FALSE){
           stop("calcGroupsW : this function requires to have installed the spam package to work \n")
  }
			
	test.package <- requireNamespace("Matrix",quietly=TRUE)
    if(test.package==FALSE){
       stop("calcGroupsW : this function requires to have installed the Matrix package to work \n")
    }
	
  if("dgCMatrix" %in% is(W)==F){
    stop("calcGroupsW : wrong specification of \'W\' \n",
         "\'W\' must be of class dgCMatrix \n",
         "is(W) : ",paste(is(W),collapse=" ")," \n")
  }
  
  if(is.null(subset)){
    test.subset <- FALSE
    subset <- 0:(nrow(W)-1)
  }else{
    test.subset <- TRUE
    subset <- subset-1
  }
  
  resCpp <- calcGroupsW_cpp(W=W,subset=subset,max_groups=max_groups)
  
  if(resCpp$nbObsRes>0){
    warning("calcGroupsW : maximum number of groups reached \n",
            "the identification of the spatial groups may not be complet \n",
            "set \'max_groups\' higher to allow more groups \n")
  }
  
  # export
  if(test.subset){
    group <- rep(NA,nrow(W))
    group[subset+1] <- resCpp$group_subset
  }
  
  res <- list(group=if(test.subset){group}else{resCpp$group_subset},
              group_subset=if(test.subset){resCpp$group_subset}else{NULL},
              group_size=resCpp$group_size,
              group_number = resCpp$group_length,
              group_max = which.max(resCpp$group_size)              
  )
  
  return(res)
}

setMethod(f ="calcW",
          signature ="data.frame",
          definition = function(object,range,method="euclidean",upper=NULL,format="dgCMatrix",row.norm=FALSE,
                                spatial_res=rep(1,ncol(object)))
          { 
			#### test package W
            test.package <- requireNamespace("spam",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcW : this function requires to have installed the spam package to work \n")
            }
			
			test.package <- requireNamespace("Matrix",quietly=TRUE)
            if(test.package==FALSE){
              stop("calcW : this function requires to have installed the Matrix package to work \n")
            }
			
            if(format %in% c("spam","dgCMatrix")==FALSE){
              stop("calcW[data.frame] : wrong specification of \'format\' \n",
                   "valid format : \"spam\" \"dgCMatrix\" \n",
                   "requested format : ",format," \n")
            }
            
            if(is.null(range) || is.na(range) || range<=0 ){
              stop("calcW[data.frame] : wrong specification of \'range\' \n",
                   "range must be positive \n",
                   "requested range : ",range," \n")
            }
            
            p <- ncol(object)
            
            if(length(spatial_res)!=p || is.numeric(spatial_res)==FALSE){
              stop("calcW[data.frame] : wrong specification of \'spatial_res\' \n",
                   "it must be a numeric vector of length ",p," \n",
                   "proposed 'spatial_res' : ",paste(spatial_res,collapse=" "),"\n")
            }
            
            for(iter_c in 1:p){
              object[,iter_c] <-  object[,iter_c]*spatial_res[iter_c]
            }
            
            
            W <- spam::nearest.dist(object,method="euclidean",delta=range,upper=upper)
            
            if(format=="dgCMatrix"){
              W <- spam::as.dgCMatrix.spam(W)
              if(row.norm){
                pSum <- spam::rowSums(W)
                pSum[pSum==0] <- -1
                W <- W/pSum
              }
              W <- Matrix::drop0(W,is.Csparse=TRUE)
            }      
            
            return(list(res=W,
                        update.object=FALSE,
                        overwrite=FALSE)
            )
          }
)


EDK <-function(x, bandwidth,power=2){ 
  1/(2*pi*bandwidth^2)^(1/power)*exp(-(x/bandwidth)^power) 
}


#### 3- Display functions ####

legendMRI <- function(breaks,palette,mar,cex,cex.main,main,quantiles,digit){
  
  ### definition de scale 
  if(length(breaks)>8){
    at_legend <- seq(min(breaks),max(breaks),length.out=8)
  }else{
    at_legend <- breaks
  }
  at_legend <- signif(at_legend,digit)
  n.legend <- length(at_legend)
  n.breaks <- length(breaks)
  seq_scale <- 1:(n.breaks-1)
  
  if(!is.null(mar)){par(mar=mar)}
  
  image(1,(breaks[-n.breaks]+breaks[-1])/2,rbind(seq_scale),col=palette,
        axes=FALSE,ylab="",xlab="")
  
  title(main,cex.main=cex.main)
  axis(2,cex.axis=cex,at=c(at_legend[-n.legend],0.99*at_legend[n.legend]),labels=at_legend,las=2)            
  
  if(!is.null(quantiles)){              
    abline(h=quantiles,col=c("red","blue","green","blue","red"),lwd=3,lty=c(1,2,3,2,1))
  }
  
  # export
  return(invisible(TRUE))
}

setMethod(f ="multiplot",
          signature ="data.frame",
          definition = function(object,contrast=NULL,num=NULL,index1=NULL,index2=NULL,index3=NULL,
                                breaks=50,type.breaks="range",palette="terrain.colors",col=NULL,pch=NULL,cex=1,
                                col.NA="lightyellow",pch.NA=8,xlim=NULL,ylim=NULL,axes=TRUE,
                                window=FALSE,legend=TRUE,mfrow=NULL,mar=rep(1.5,4),mgp=c(2,0.5,0),pty=NULL,asp=1,bg="lightblue",
                                xlab="",ylab="", main="slice ",num.main=TRUE,cex.main=1.5,
                                quantiles.legend=TRUE, digit.legend=3,cex.legend=1.5,mar.legend=c(2,7,2,2),main.legend=names(contrast),    
                                filename="multiplot",width=1000,height=500,path=NULL,unit="px",res=NA){
            
            #### test ####
            if(is.data.frame(contrast)){
              param <- names(contrast)
              contrast <- as.matrix(contrast)
            }else if(is.matrix(contrast)){
              param <- NULL
            }else if(is.vector(contrast)){
              contrast <- as.matrix(contrast)
              param <- NULL
            }else if(is.null(contrast)){
              contrast <- matrix(1,nrow=nrow(object),ncol=1)
              param <- NULL
            }else{
              stop("multiplot[data.frame] : wrong specification of \'contrast\' \n",
                   "must be either NULL or data.frame or matrix or vector \n",
                   "is(contrast) : ",paste(is(contrast),collapse=" "),"\n")
            }
            
            if(ncol(object)!=3){
              stop("multiplot[data.frame] : \'object\' must be a 3 column matrix with the observation coordinates \n",
                   "ncol(object) : ",ncol(object),"\n")
            }
            
            if(is.null(num)){num <- unique(object[,3])}
            
            if(nrow(contrast)!=nrow(object)){
              stop("multiplot[data.frame] : \'object\' and  \'contrast\' must have the same number of rows \n",
                   "nrow contrast : ",nrow(contrast),"\n",
                   "nrow object : ",nrow(object),"\n")
            }
            
            if(any(num %in% object[,3])==FALSE){
              stop("multiplot[data.frame] : \'object\' inconsistent with \'num\' \n",
                   "requested slices : ",paste(num,collapse=" ")," \n",
                   "slices in \'object\' : ",paste(unique(object[,3]),collapse=" ")," \n")
            }
            
            n.plot <- length(num)
            if(!is.null(legend) && legend==TRUE){n.plot <- n.plot + 1}
            
            contrast <- contrast[num %in% object[,3],,drop=FALSE]
            object <- object[num %in% object[,3],]
            n.px <- nrow(contrast)
            
            #### initialization and tests ####
            
            ## windows
            mar.init <- par()$mar
            
            res.init <- initWindow(window=window,filename=filename,path=path,width=width,height=height,unit=unit,res=res,
                                            n.plot=n.plot,mfrow=mfrow,xlim=xlim,ylim=ylim,
                                            method="multiplot[data.frame]")
            scale <- res.init$scale
            mfrow <- res.init$mfrow
            n.graph_par_window <- res.init$n.graph_par_window
            xlim.plot <- res.init$xlim.plot
            ylim.plot <- res.init$ylim.plot
            
            ## color and breaks   
            res.init <- initCol(contrast=contrast,coords=object,param=param,pch=pch,col=col,palette=palette,breaks=breaks,legend=legend,type.breaks=type.breaks,
                                         method="multiplot[data.frame]")
            
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
              res.init <- initIndex(object=NULL,index=index1,num=num,hemisphere=hemisphere,as.logical=as.logical,
                                             method="multiplot[data.frame]",indexNum=1,
                                             cex.default=1,pch.default=20,col.default="red",filter_default="2D_N4")
              index1 <- res.init$coords
              indexindex1 <- res.init$index
              pch_index1 <- res.init$pch
              cex_index1 <- res.init$cex
              col_index1 <- res.init$col
            }
            
            if(!is.null(index2)){
              res.init <- initIndex(object=NULL,index=index2,num=num,hemisphere=hemisphere,as.logical=as.logical,
                                             method="multiplot[data.frame]",indexNum=2,
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
                  stop("multiplot[data.frame] : \'index3\' cannot request \"interaction\" parameter if \'index1\' is missing  \n")
                }
                if(is.null(index2)){
                  stop("multiplot[data.frame] : \'index3\' cannot request \"interaction\" parameter if \'index2\' is missing  \n")
                }
                
                test1.intersect <- indexindex1 %in% indexindex2
                test2.intersect <- indexindex2 %in% indexindex1
                
                index1_sauve <- index1
                index3$coords <- index1_sauve$coords[test1.intersect==TRUE,c("i","j","k")]
                index1 <- index2$coords[test1.intersect==FALSE,c("i","j","k")]
                index2 <- index1_sauve$coords[test2.intersect==FALSE,c("i","j","k")]
                
              }
              
              res.init <- initIndex(object=NULL,index=index3,num=num,hemisphere=hemisphere,as.logical=as.logical,
                                             method="multiplot[data.frame]",indexNum=3,
                                             cex.default=1,pch.default=22,col.default="green",filter_default="2D_N4")
              index3 <- res.init$coords   
              pch_index3 <- res.init$pch
              cex_index3 <- res.init$cex
              col_index3 <- res.init$col     
              
            }
            
            
            #### display plot ####
            compteur <- 1
            
            for(iter_num in 1:length(num)){
              
              # device
              if(!is.null(window) && compteur == 1){
                
                if(window %in% c("png","eps","svg","pdf") && iter_num>1){dev.off()}
                filename_all <- paste(filename,"_",param,"(slice",num[iter_num],"-",min(max(num),num[iter_num+n.graph_par_window-1],na.rm=TRUE),")",sep="")                  
                
                initDisplayWindow(window=window,filename=filename_all,path=path,width=width,height=height,scale=scale,res=res,
                                           mfrow=mfrow,bg=bg,pty=pty,mar=mar,mgp=mgp)
                
              }
              
              # data
              index_k <- which(object[,3]==num[iter_num])
              contrastK <- contrast[index_k,,drop=FALSE]
              coordsK <- object[index_k,,drop=FALSE]
              if(is.null(col)){colK <- NULL}else{colK <- col[index_k]}
              if(num.main==TRUE){mainK <- paste(main," ",num[iter_num],if(!is.null(param)){paste(" (",param,")",sep="")},sep="")}else{mainK <- main}
              
              # xlim-ylim
              if(is.null(xlim)){                
                if(is.null(coordsK) || nrow(coordsK)==0){
                  xlim.plot <- NULL
                }else{
                  xlim.plot <- c(min(coordsK[,1])-0.5,max(coordsK[,1])+0.5)
                }                 
              }
              
              if(is.null(ylim)){
                if(is.null(coordsK) || nrow(coordsK)==0){
                  ylim.plot <- NULL
                }else{
                  ylim.plot <- c(min(coordsK[,2])-0.5,max(coordsK[,2])+0.5)
                }
              }
              
              plot.test <- plotMRI(contrast=contrastK,coords=coordsK[,1:2],breaks=breaks,col=colK,palette=palette,
                                   asp=asp,
                                   xlim=xlim.plot,ylim=ylim.plot,pch=pch,cex=cex,axes=axes,col.NA=col.NA,pch.NA=pch.NA,
                                   xlab=xlab,ylab=ylab,main=mainK,cex.main=cex.main)
              
              if(!is.null(index1) && any(index1[,3]==num[iter_num])){
                points(index1[index1[,3]==num[iter_num],1:2,drop=FALSE],
                       pch=pch_index1,cex=cex_index1,col=col_index1)
              }
              
              if(!is.null(index2) && any(index2[,3]==num[iter_num])){
                points(index2[index2[,3]==num[iter_num],1:2,drop=FALSE],
                       pch=pch_index2,cex=cex_index2,col=col_index2)
              }
              
              if(!is.null(index3) && any(index3[,3]==num[iter_num])){
                points(index3[index3[,3]==num[iter_num],1:2,drop=FALSE],
                       pch=pch_index3,cex=cex_index3,col=col_index3)
              }
              
              compteur <- compteur + 1
              if(compteur > n.graph_par_window)
              {compteur <- 1}
              
            }            
            
            #### legend ####
            if(is.null(legend)){
              compteur <- 1
              mfrow <- c(1,1)
              legend <- TRUE
            }
     
            if(legend==TRUE){
              
              if(!is.null(window) && compteur == 1){
                
                if(window %in% c("png","eps","svg","pdf") && iter_num>1){dev.off()}
                filename_all <- paste(filename,"_",param,"(slice",min(num),"-",max(num),")-legend",sep="")
                
                initDisplayWindow(window=window,filename=filename_all,width=width,height=height,scale=scale,res=res,
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

outline <- function(n=50,sequential=TRUE,min_dist=1,
                    col=c("blue","red","grey"),pch=20,cex=c(0.75,1,0.75)){
  
  test.package <- requireNamespace("RANN",quietly=TRUE)
  if(test.package==FALSE){
    stop("outline : this function requires to have installed the RANN package to work \n")
  }
  
  #### one shot outline 
  if(sequential==FALSE){
    res_locator <- locator(n=n)
    res_locator$x <- c(res_locator$x,res_locator$x[1])
    res_locator$y <- c(res_locator$y,res_locator$y[1])
    n <- length(res_locator$x)
    
    # round
    res_locator$x <- round(res_locator$x)
    res_locator$y <- round(res_locator$y)  
    
    
    df_points <- matrix(NA,nrow=0,ncol=4)
    
    for(iter_point in 1:(n-1)){
      n_tempo <- max(1+abs(res_locator$x[iter_point]-res_locator$x[iter_point+1]),
                     1+abs(res_locator$y[iter_point]-res_locator$y[iter_point+1]))
      n_tempo <- ceiling(n_tempo)
      i_tempo <- seq(res_locator$x[iter_point],res_locator$x[iter_point+1],length.out=n_tempo)
      j_tempo <- seq(res_locator$y[iter_point],res_locator$y[iter_point+1],length.out=n_tempo)
      
      i_tempo <- round(i_tempo)
      j_tempo <- round(j_tempo)    
      
      if(n_tempo>1){
        df_points <- rbind(df_points,
                           cbind(i=i_tempo[-n_tempo],j=j_tempo[-n_tempo],edge=iter_point,points=c(1,rep(0,n_tempo-2))))
      }
    }
    
  }
  
  
  #### sequential outline
  if(sequential==TRUE){
    
    iter <- 1
    dist <- min_dist+1
    df_points <- matrix(NA,nrow=0,ncol=4)
    
    while(iter <= n && dist>min_dist){
      res_locator <- locator(n=1)
      
      # if user enter echap directly
      if(is.null(res_locator)){
        return(list(edge=NULL,
                    surface=NULL)
        )
      }
      
      if(iter>1){
        dist <- sqrt( (df_points[1,"i"]-res_locator$x)^2 + (df_points[1,"j"]-res_locator$y)^2)
        
        res_locator$x <- c(df_points[nrow(df_points),"i"],round(res_locator$x))
        res_locator$y <- c(df_points[nrow(df_points),"j"],round(res_locator$y))
        
        n_tempo <- max(1+abs(res_locator$x[1]-res_locator$x[2]),
                       1+abs(res_locator$y[1]-res_locator$y[2]))
        n_tempo <- ceiling(n_tempo)
        if(n_tempo==1){i_tempo <- c() ; next} # cas ou l on selectionne le meme points
        
        i_tempo <- round(seq(res_locator$x[1],res_locator$x[2],length.out=n_tempo)[-1])
        j_tempo <- round(seq(res_locator$y[1],res_locator$y[2],length.out=n_tempo)[-1])
        
        points <- c(rep(0,n_tempo-2),1)
      }else{
        
        # round
        res_locator$x <- round(res_locator$x)
        res_locator$y <- round(res_locator$y)  
        
        i_tempo <- res_locator$x
        j_tempo <- res_locator$y
        points <- 1  
      }
      
      if(dist<1){
        if(length(i_tempo)>1){
          df_points <- rbind(df_points,
                             cbind(i=i_tempo[-(n_tempo-1)],j=j_tempo[-(n_tempo-1)],edge=iter,points=0))
        }
      }else{
        df_points <- rbind(df_points,
                           cbind(i=i_tempo,j=j_tempo,edge=iter,points=points))
      }
      
      points(i_tempo,j_tempo,
             col=col[points+1],pch=pch,cex=cex[points+1])
      
      
      
      iter <- iter + 1
    }
    
    # relier au premier point   
    n_tempo <- max(1+abs(df_points[nrow(df_points),"i"]-df_points[1,"i"]),
                   1+abs(df_points[nrow(df_points),"j"]-df_points[1,"j"]))
    n_tempo <- ceiling(n_tempo)
    
    if(n_tempo>2){
      i_tempo <- round(seq(df_points[nrow(df_points),"i"],df_points[1,"i"],length.out=n_tempo))
      j_tempo <- round(seq(df_points[nrow(df_points),"j"],df_points[1,"j"],length.out=n_tempo))
      
      df_points <- rbind(df_points,
                         cbind(i=i_tempo[c(-1,-n_tempo)],j=j_tempo[c(-1,-n_tempo)],edge=iter,points=0))
      
      points(i_tempo[c(-1,-n_tempo)],j_tempo[c(-1,-n_tempo)],
             col=col[points+1],pch=pch,cex=cex[points+1])
    }
  }
  df_points <- as.data.frame(df_points)
  
  # enlever d eventuels dupliques
  df_points <- df_points[duplicated(df_points[,c("i","j")])==FALSE,]
  
  
  #### filling the outline
  df_points <- cbind(df_points,valid=TRUE)
  j_level <- unique(df_points[,"j"])
  i_tempo <- c()
  j_tempo <- c()  
  
  for(iter_j in 1:length(j_level)){
    
    index_j <- which(df_points[,"j"]==j_level[iter_j])
    
    if(length(index_j)==1){
      i_tempo <-  c(i_tempo,df_points[index_j,"i"])
      j_tempo <-  c(j_tempo,df_points[index_j,"j"])
    }else{
      index_j <- index_j[order(df_points[index_j,"i"])]
      
      # enlever les doublons d un meme trait
      diff_tempo <- c(FALSE,abs(diff(index_j)) %in% c(1,nrow(df_points)-1))
      
      while(TRUE %in% diff_tempo && length(diff_tempo)>2){   
        pos_FALSE <- which(diff_tempo==TRUE)[1]
        if(pos_FALSE %% 2==0){
          df_points$valid[index_j[pos_FALSE]] <- FALSE
          index_j <- index_j[-pos_FALSE]          
        }else{          
          df_points$valid[index_j[pos_FALSE-1]] <- FALSE
          index_j <- index_j[-(pos_FALSE-1)]
        }
        diff_tempo <- diff_tempo[-pos_FALSE]
      }
      
      # enlever les V inverses
      df_Vpoints <- cbind(index=1:nrow(df_points),df_points)[df_points$valid==TRUE,]
      index_Vj <- sapply(index_j,function(x){which(df_Vpoints$index==x)})
      index_m1 <- c(nrow(df_Vpoints),1:nrow(df_Vpoints)-1)
      index_p1 <- c(2:nrow(df_Vpoints),1)
      
      diff_tempo <- rbind(df_Vpoints[,"i"]-df_Vpoints[index_m1,"i"],
                          df_Vpoints[,"j"]-df_Vpoints[index_m1,"j"],
                          df_Vpoints[,"i"]-df_Vpoints[index_p1,"i"],
                          df_Vpoints[,"j"]-df_Vpoints[index_p1,"j"])[,index_Vj,drop=FALSE]
      
      test_V <- apply(diff_tempo,2,function(x){
        test.H <- all(x==c(1,-1,-1,-1)) + all(x==c(0,-1,-1,-1)) + all(x==c(0,-1,0,-1)) + all(x==c(1,-1,0,-1))
        test.antiH <- all(x==c(-1,-1,1,-1)) + all(x==c(-1,-1,0,-1)) + all(x==c(0,-1,1,-1))
        return(test.H+test.antiH)       
      }) 
      
      test_Vinv <- apply(diff_tempo,2,function(x){
        test.H <- all(x==c(1,1,-1,1)) + all(x==c(0,1,-1,1)) + all(x==c(0,1,0,1)) + all(x==c(1,1,0,1))
        test.antiH <- all(x==c(-1,1,1,1)) + all(x==c(-1,1,0,1)) + all(x==c(0,1,1,1))
        return(test.H+test.antiH)        
      }) 
      
      i_tempo <- c(i_tempo,df_points[index_j[test_V==1],"i"])
      j_tempo <- c(j_tempo,df_points[index_j[test_V==1],"j"])
      index_j <- index_j[c(test_Vinv+test_V)==0]
      
      n.i <- length(index_j)
      
      # remplissage
      if(n.i>0){
        for(iter_i in 1:floor(n.i/2)){
          seq_i <- seq(df_points[index_j[2*(iter_i-1)+1],"i"],df_points[index_j[2*iter_i],"i"],by=1)
          i_tempo <-  c(i_tempo,seq_i)
          j_tempo <-  c(j_tempo,rep(j_level[iter_j],length(seq_i)))
        }
      }
      
    }
    
  }
  
  df_fill <- data.frame(i=i_tempo,j=j_tempo)
  
  if(sequential==TRUE){
    test <- RANN::nn2(data=df_points[,c("i","j")],query=df_fill[,c("i","j")],k=1)
    points(i_tempo[test$nn.dists>0],j_tempo[test$nn.dists>0],
           col=col[3],pch=pch,cex=cex[3])
  }
  
  
  return(list(edge=df_points,
              surface=df_fill)
  )
}

plotMRI <- function(contrast,coords,breaks,palette,col,asp,
                    xlim,ylim,pch,cex,axes,col.NA,pch.NA,xlab,ylab,main,cex.main){
  
  #### gestion des coupes ####
  if(is.null(main) && !is.null(names(contrast)))
  {main <- names(contrast)}
  
  if(nrow(contrast)==0 || sum(is.na(contrast)==FALSE) == 0 || is.null(contrast)){
    
    if(is.logical(xlim) || is.logical(ylim))
    {plot(1,1,col=col,xlab="",ylab="",axes=axes)}else
    {plot(1,1,col=col,xlim=xlim,ylim=ylim,xlab="",ylab="",axes=axes)}
    legend("center",c("Missing","data"),cex=2,pch=62,bty="n")
    
    title(main=main)
    axis(1)
    axis(2)
    
    return(invisible(FALSE))
  }
  
  #### mis en forme des donnees ####
  if(is.null(col)){
    array <- df2array(contrast=contrast,coords=coords,format="matrix",default_value=NA)
    spdf <- array$contrast
    coords_x <- array$unique_coords[[1]]
    coords_y <- array$unique_coords[[2]]
  }
  
  #### taille du pixel ####
  if(!is.null(col)){
    
    unique_x <- sort(unique(coords[,1]))
    n_x <- ( unique_x[length(unique_x)]-unique_x[1] ) # /  min(diff(unique_x))
    
    unique_y <- sort(unique(coords[,2]))
    n_y <- ( unique_y[length(unique_y)]-unique_y[1] ) # /  min(diff(unique_y))
    
    cex <- max(3*par("pin")/(par("cin")*max(n_y, n_x)))*cex
    #     cex=100*min(par()$plt[2]-par()$plt[1],par()$plt[4]-par()$plt[3])/max(xlim[2]-xlim[1],ylim[2]-ylim[1]) 
  }
  
  #### affichage #### 
  if(is.null(col)){
    eval(parse(text=paste(
      "image(x=coords_x,y=coords_y,z=spdf,
      breaks=breaks,col=palette,
      ",if(!is.null(xlim)){"xlim=xlim,"},if(!is.null(ylim)){"ylim=ylim,"},
      "xlab=xlab,ylab=ylab,axes=axes,asp=asp)",
      sep="")))
    
  }else{
    
    plot(coords,
         pch=pch,cex=cex,col=col,xlim=xlim,ylim=ylim,
         axes=axes,asp=asp,xlab=xlab,ylab=ylab)
    
  }
  
  title(main=main,cex.main=cex.main)
  # affichage des NA
  if(!is.null(col.NA)){
    test_NA <- which(is.na(contrast))
    
    if(sum(test_NA)>0){
      points(coords[test_NA,],col=col.NA,pch=pch.NA,cex=cex[1]/5)
    }
    
  }
  
  # export
  return(invisible(TRUE))
  
}

pointsOutline <- function(coords,array=NULL,filter="2D_N4"){
  
  if(is.null(array)){
    p <- ncol(coords)
    n.neighbors <- as.numeric(paste(strsplit(filter,split="",fixed=TRUE)[[1]][-(1:4)],collapse=""))
    if(nrow(coords)<=n.neighbors){
      return(coords)
    }      
    
    array <- df2array(contrast=rep(1,nrow(coords)),
                               coords=coords)$contrast[[1]]    
  }else{
    coords <- array2df(array)[,1:length(dim(array))]
  }
  
  res <- calcFilter(array,filter=filter,norm=FALSE)
  index_array <- which(!is.na(array))
  index_outline <- which(res$res<max(res$res,na.rm=T))
  
  coords <- coords[index_array %in% index_outline,]
  
  return(coords)
  
}


#### 4- Filtering functions ####

setMethod(f ="calcFilter",
          signature ="array",
          definition = function(object,filter,norm=TRUE,w_contrast=FALSE,na.rm=FALSE){
            
            #### test
            if(!is.logical(norm)){
              stop("calcFilter[array] : wrong specification of \'norm\' \n",
                   "\'norm\' should be logical \n",
                   "is(norm) : ",paste(is(norm),collapse=" "),"\n")
            }
            
            if(!is.logical(w_contrast)){
              stop("calcFilter[array] : wrong specification of \'w_contrast\' \n",
                   "\'w_contrast\' should be logical \n",
                   "is(w_contrast) : ",paste(is(w_contrast),collapse=" "),"\n")
            }
            
            if(!is.logical(na.rm)){
              stop("calcFilter[array] : wrong specification of \'na.rm\' \n",
                   "\'na.rm\' should be logical \n",
                   "is(na.rm) : ",paste(is(na.rm),collapse=" "),"\n")
            }
            
            if(length(dim(object)) %in% c(2,3) == FALSE){
              stop("calcFilter[array] : wrong specification of \'object\' \n",
                   "\'object\' must be an array of dimension 2 or 3 \n",
                   "dim(object) : ",paste(dim(object),collapse=" "),"\n")
            }
            
            #### Filtres disponibles 
            if(length(filter)==1 && is.character(filter)){
              
              filter_split <- strsplit(filter,split="")[[1]]
              
              if(filter_split[[4]] == "N"){
                
                res <- initNeighborhood(filter,method="calcFilter")
                filter_split <- list(ncol(res),"D","_",filter_split[[4]],nrow(res))
                
                if(filter_split[[1]]==2){
                  filter <- matrix(NA,nrow=3,ncol=3)
                  for(iter in 1:filter_split[[5]]){
                    filter[res[iter,1]+2,res[iter,2]+2] <- 1  
                  }
                }else{
                  filter <- array(NA,dim=c(3,3,3))
                  for(iter in 1:filter_split[[5]]){
                    filter[res[iter,1]+2,res[iter,2]+2,res[iter,3]+2] <- 1  
                  }
                }
                
              }else{
                
                res <- initFilter(filter,method="calcFilter")
                filter <- res$filter
                filter_split <- res$filter_split
                
              }
            }else{ # Perso
              filter_split <- c(length(dim(filter)),"D","_","P",dim(filter)[1])
            }
            
            if(filter_split[[4]]=="S" && na.rm==FALSE){
              warning("calcFilter[array] : gradient values may be incorrect at the edges \n",
                      "set \'na.rm\' to TRUE to remove these values \n")
            }
            if(filter_split[[4]]=="M" && w_contrast==TRUE){
              warning("calcFilter[array] : there is no edge preserving median filtering \n",
                      "\'w_contrast\' will be ignored \n")
            }
            
            ### preparation
            p <- dim(filter)
            p.dim <- length(p)
            p_ref <- sapply(1:p.dim,function(x){median(1:p[x])})
            if(p.dim > 3){
              stop("calcFilter[array] : wrong specification of \'filter\' \n",
                   "can only handle 1D 2D and 3D filters \n",
                   "dimension of the proposed filter : ",p.dim,"\n")
            }
            
            M.dim <- dim(object)
            Mres <- array(NA,dim=dim(object))
            
            Ind.operateur_compr <- which(as.vector(filter)!=0)
            Vec.operateur_compr <- as.vector(filter)[Ind.operateur_compr]
            
            # filtrage 2D
            if(filter_split[[4]] %in% c("N","G","I","P","S") && length(M.dim)==2){
              
              resCpp <- filtrage2D_cpp(M_data=object,
                                       M_operateur=filter,
                                       index_data=which(!is.na(object),arr.ind=TRUE)-1,
                                       w_contrast=w_contrast,
                                       na_rm=na.rm)
              
              if(norm==TRUE){
                Mres <- resCpp$Mres/resCpp$Wres
              }else{
                Mres <- resCpp$Mres
              }
              
            }
            
            
            # filtrage 3D
            if(filter_split[[4]] %in% c("N","G","I","P","S") && length(M.dim)==3){
              if(filter_split[[1]]==2){
                for(iter_k in 1:(M.dim[3])){
                 
                  resCpp <- filtrage2D_cpp(M_data=object[,,iter_k],
                                           M_operateur=filter,
                                           index_data=which(!is.na(object[,,iter_k]),arr.ind=TRUE)-1,
                                           w_contrast=w_contrast,
                                           na_rm=na.rm)
                  
                  if(norm==TRUE){
                    Mres[,,iter_k] <- resCpp$Mres/resCpp$Wres
                  }else{
                    Mres[,,iter_k] <- resCpp$Mres
                  }      
                }
              }else{  
      
                resCpp <- filtrage3D_cpp(Vec_data=as.vector(object),p_data=dim(object),
                                         Vec_operateur=as.vector(filter),p_operateur=dim(filter),
                                         index_data=which(!is.na(object),arr.ind=TRUE)-1,
                                         w_contrast=w_contrast,
                                         na_rm=na.rm)                
                
                if(norm==TRUE){
                  Mres <- resCpp$Mres/resCpp$Wres
                }else{
                  Mres <- resCpp$Mres
                }
              }
              
            }           
            
            # filtrage median 2D
            if(filter_split[[4]] == "M" && length(M.dim)==2){
              Mres <- filtrage2Dmed_cpp(M_data=object,
                                        M_operateur=filter,
                                        index_data=which(!is.na(object),arr.ind=TRUE)-1,
                                        na_rm=na.rm)
            }
            
            # filtrage median 3D
            if(filter_split[[4]] == "M" && length(M.dim)==3){
              if(filter_split[[1]]==2){
                for(iter_k in 1:(M.dim[3])){
                  Mres[,,iter_k] <- filtrage2Dmed_cpp(M_data=object[,,iter_k],
                                                      M_operateur=filter,
                                                      index_data=which(!is.na(object[,,iter_k]),arr.ind=TRUE)-1,
                                                      na_rm=na.rm)
                }
              }else{
                
                Mres <- filtrage3Dmed_cpp(Vec_data=as.vector(object),p_data=dim(object),
                                          Vec_operateur=as.vector(filter),p_operateur=dim(filter),
                                          index_data=which(!is.na(object),arr.ind=TRUE)-1,
                                          na_rm=na.rm) 
              }
            }
            
            
            ### export
            return(list(res=Mres,
                        filter=filter,
                        update.object=FALSE,
                        overwrite=FALSE)
            )
          }
)

#### 5- init functions ####

initCol <- function(contrast,coords,param=NULL,pch,col,palette,breaks,legend,type.breaks,method){
  
  # > image is used in all but three cases
  # col has been specifed by the user
  # there are more than one parameter (color have to be defined with a multiparametric palette)
  # user has specified a pch value
  
  p <- ncol(contrast)
  n.px <- nrow(contrast)
  if(is.null(param)){param <- names(contrast)}
  
  #### tests ####
  if(type.breaks %in% c("range","quantile","range_center") ==FALSE){
    stop(method," : wrong specification of \'type.breaks\' \n",
         "valid values : \"range\" \"range_center\" \"quantile\" \n",
         "proposed value : ",type.breaks,"\n")
  }
  
  if(!is.null(legend) && !is.logical(legend)){
    stop(method," : wrong specification of \'legend\' \n",
         "valid values : NULL FALSE TRUE \n",
         "proposed value : ",legend,"\n")
  }
  
  if(!is.numeric(breaks)){
    stop(method," : wrong specificaiton of \'breaks\' \n",
         "\'breaks\' must be either a numeric indicating the number of breaks or a numeric vector containing the break values \n",
         "proposed \'break\' is(break) : ",paste(is(breaks),collapse=" "),"\n")
  }
  
  if((p %in% 1:3) == FALSE){
    if(method=="multiplot[MRIaggr]"){
      stop(method," : wrong specification of \'param\' \n",
           "\'param\' must contains between 1 and 3 parameters \n",
           "proposed parameters (nb) : ",paste(names(contrast),collapse=" ")," (",ncol(contrast),") \n")  
    }else{
      stop(method," : wrong specification of \'contrast\' \n",
           "\'contrast\' must have between 1 and 3 columns \n",
           "columns of \'data\' (ncol) : ",paste(names(contrast),collapse=" ")," (",ncol(contrast),") \n")    
    }
  }
  
  
  #### dealing with multiple parameters ####
  index_duplicated <- NULL
  index_order <- NULL
  
  if(p %in% c(2:3)){
    
    if(length(palette)!=1 || any(palette %in% c("rgb","hsv") == FALSE)){
      stop(method," : wrong specification of \'palette\' \n",
           "available palette for 2 or 3 imaging parameters : \"rgb\" \"hsv\"\n",
           "proposd palette : ",paste(palette,collapse=" "),"\n")
    }
    
    if(any(contrast>1) || any(contrast<0)){
      if(method=="multiplot[MRIaggr]"){
        stop(method," : wrong specification of \'param\' \n",
             "parameters must take value in [0;1] if several parameters are intended to be displayed on the same map \n",                   
             "current range : ",paste(range(contrast),collapse=" "),"\n")
      }else{
        stop(method," : wrong specification of \'contrast\' \n",
             "contrast must contains values be in [0;1] if several parameters are intended to be displayed on the same map \n",                   
             "current range : ",paste(range(contrast),collapse=" "),"\n")
      }
    }
    
    if(ncol(contrast)==2){      
      col <- eval(parse(text=paste("grDevices::",palette,"(contrast[,1],contrast[,2],1)",sep="")))      
      contrast <- contrast[,1,drop=FALSE]
    }else{
      
      col <- eval(parse(text=paste("grDevices::",palette,"(contrast[,1],contrast[,2],contrast[,3])",sep="")))
      
      # ordering dataset by Tier
      order1 <- order(contrast[,1],decreasing=TRUE)[seq(1,n.px,length.out=round(n.px/3))]            
      order2 <- setdiff(order(contrast[,2],decreasing=TRUE),order1)[seq(1,n.px-round(n.px/3),length.out=round(n.px/3))]
      order3 <- setdiff(order(contrast[,3],decreasing=TRUE),c(order1,order2))[seq(1,n.px-2*round(n.px/3),length.out=n.px-2*round(n.px/3))]
      
      contrast[order2,1] <- contrast[order2,2]+1
      contrast[order3,1] <- contrast[order3,3]+2*1
      contrast <- contrast[,1,drop=FALSE]            
      
      index_duplicated <- which(duplicated(contrast[,1])==FALSE)
      index_order <- order(contrast[index_duplicated,1])
    }
    
    param <- paste(param,collapse="-")
  }
  
  #### come back to one parameter
  if(is.null(col)){
    
    ## test palette
    available_palette <- c("grey.colors","gray.colors","rainbow","heat.colors","terrain.colors","topo.colors","cm.colors","tim.colors")
    
    if(length(palette)==1 && palette %in% available_palette == FALSE){
      unique_tempo <-  unique(contrast[,1])
      if(length(unique_tempo)==1){
        palette <- rep("red",2)
        breaks <- c(unique_tempo-1*10^{-12},unique_tempo,unique_tempo+1*10^{-12})
      }else{
        stop(method," : wrong specification of \'palette\' \n",
             "available palette for 1 parameter : \"",paste(available_palette,collapse="\" \""),"\" \n",
             "proposd palette : ",palette,"\n")
      }
    }
    
    ## infinte values
    index.infinite <- which(is.infinite(contrast[,1]))
    
    if(length(index.infinite)>0){
      maxInf <-  max(c(contrast[-index.infinite,1],99999))
      if(method=="multiplot[MRIaggr]"){
        warning(method," : \'param\' values in \'object\' contains Inf values \n",
                "they are set to ",maxInf,"\n")     
      }else{
        warning(method," : \'contrast\' values contains Inf values \n",
                "they are set to ",maxInf,"\n")          
      }
      contrast[index.infinite,1] <- maxInf   
    }
    
    ## breaks 
    if(length(breaks)==1){
      range_data <- range(na.omit(contrast[,1]))
      
      if(length(palette)>1){breaks <- length(palette)+1}
      
      breaks_sauve <- switch(type.breaks,
                             "range"=seq(range_data[1],range_data[2],length.out=min(500,breaks)),
                             "range_center"=seq(-max(abs(range_data)),max(abs(range_data)),length.out=min(500,breaks)),
                             "quantile"=quantile(contrast[,1],probs=seq(0,1,length.out=min(500,breaks)))
      )
      breaks_sauve <- unique(breaks_sauve)
      if(length(breaks_sauve)==1){breaks_sauve <- c(breaks_sauve-1*10^{-12},breaks_sauve,breaks_sauve+1*10^{-12})}
      
      breaks <- breaks_sauve
      breaks[1] <- breaks[1]-1*10^{-12}
      breaks[length(breaks)] <- breaks[length(breaks)]+1*10^{-12}
      
    }else{
      breaks <- sort(breaks)
      breaks_sauve <- breaks
    }
    
    ## palette
    if(length(palette)==1 && palette %in% available_palette){
      
      palette <- eval(parse(text=paste(palette,"(length(breaks)-1)",sep="")))
    }
    palette_sauve <- palette
    
    ## integrate extreme contrast in range
    
    if(max(breaks)<=max(contrast[,1],na.rm=TRUE)){
      breaks[length(breaks)] <-max(contrast[,1],na.rm=TRUE)+1*10^{-12}
    }
    if(min(breaks)>=min(contrast[,1],na.rm=TRUE)){
      breaks[1] <- min(contrast[,1],na.rm=TRUE)-1*10^{-12}           
    }
    
    ## check breaks - palette
    if(length(breaks)!=length(palette)+1){
      stop(method," : \'breaks\' and \'palette\' are incompatible \n",
           "length(palette) must be equal to length(breaks)+1 \n",
           "length(palette) : ",length(palette)," \n",
           "length(breaks) : ",length(breaks)," \n")
    }
    
    if( (length(unique(breaks))!=length(breaks)) || sum( abs(sort(breaks)-breaks)>1e-16 )>0){
      stop(method," : the elements of \'breaks\' must be distinct and in ascending order \n",
           "length(breaks) : ",length(breaks),"\n",
           "length(unique(breaks)) : ",length(unique(breaks)),"\n",
           "rank(breaks) : ",paste(rank(breaks),collapse=" "),"\n")
    }
    
    if(!is.null(pch) || any(coords %% 1 >0)){
      col <- palette[findInterval(contrast[,1], breaks, all.inside=TRUE)]
      col[is.na(contrast[,1])] <- NA
      if(is.null(pch)){pch <- 15}
    }
    
  }else{
    palette_sauve <- unique(col[order(contrast[,1])])
    breaks_sauve <- c(min(contrast,na.rm=TRUE)-1*10^{-12},
                      seq(min(contrast[,1],na.rm=TRUE),max(contrast[,1],na.rm=TRUE),length.out=length(palette_sauve)-1),
                      max(contrast,na.rm=TRUE)+1*10^{-12})
    
    if(any(is.na(contrast[,1]))){
      col[is.na(contrast[,1])] <- NA
    }
    
    if(is.null(pch)){pch <- 15}
  }
  
  #### export ####
  res <- list()
  res$contrast <- contrast
  res$palette <- palette
  res$breaks <- breaks
  res$col <- col
  res$pch <- pch
  res$palette_sauve <- palette_sauve
  res$breaks_sauve <- breaks_sauve
  res$index_duplicated <- index_duplicated
  res$index_order <- index_order
  res$param <- param
  return(res)
  
}

initDisplayWindow <- function(window,filename,path,width,height,scale,res,
                              mfrow,bg,pty,mar,mgp){
  
  if(window %in% c("png","eps","svg","pdf")){
    switch(window,
           "eps"=postscript(file=paste(path,filename,".eps",sep=""),width=width*scale/90,height=height*scale/90,horizontal = FALSE,onefile=FALSE,paper = "special"),
           "svg"=svg(filename=paste(path,filename,".svg",sep=""),width=width*scale/90,height=height*scale/90,onefile=FALSE),
           "png"=png(filename=paste(path,filename,".png",sep=""),width=width*scale,height=height*scale,res=res),
           "pdf"=pdf(file=paste(path,filename,".pdf",sep=""),width=width*scale/90,height=height*scale/90,onefile=FALSE,paper = "special")
    )
  }
  
  if(window==TRUE){dev.new()}
  
  if(!is.null(mfrow)){par(mfrow=mfrow)}
  if(!is.null(bg)){par(bg=bg)}
  if(!is.null(pty)){par(pty=pty)}
  if(!is.null(mar)){par(mar=mar)}   
  if(!is.null(mgp)){par(mgp=mgp)}     
  
  
}

initFilter <- function(filter,method){
  
  filter_name <- as.character(filter)
  filter_split <- as.list(strsplit(filter_name,split="")[[1]])
  
  #### tests ####
  # 1
  if(filter_split[[1]] %in% c("2","3")==FALSE){
    stop(method," : wrong specification of \'filter\' \n",
         "valid filter[1] : \"2\" \"3\" \n",
         "proposed filter[1] :  ",filter_split[[1]],"\n")
  }
  filter_split[[1]] <- as.numeric(filter_split[[1]])
  
  # 2-3
  if(filter_split[[2]] != "D" || filter_split[[3]] != "_" ){
    stop(method," : wrong specification of \'filter\' \n",
         "valid filter[2:3] : \"D_\" \n",
         "proposed filter[2:3] :  ",filter_split[[2]]," ",filter_split[[3]],"\n")
  }
  
  # 4
  if(filter_split[[4]] %in% c("G","M","S","I") ==FALSE){
    stop(method," : wrong specification of \'filter\' \n",
         "valid filter[4] : \"G\", \"M\", \"S\" or \"I\" \n",
         "proposed filter[4] :  ",filter_split[[4]],"\n")
  }
  
  # 5-6
  if(filter_split[[4]]=="S"){  
    
    if(filter_split[[1]]==2 && filter_split[[5]] %in% c("x","y") == FALSE){
      stop(method," : wrong specification of \'filter\' \n",
           "if the fourth letter of \'filter\' is \"S\" then the fifth must be \"x\" or \"y\" \n",
           "proposed letter: ",filter_split[[5]],"\n")
    }
    if(filter_split[[1]]==3 && filter_split[[5]] %in% c("x","y","z") == FALSE){
      stop(method," : wrong specification of \'filter\' \n",
           "if the fourth letter of \'filter\' is \"S\" then the fifth must be \"x\", \"y\" or \"z\" \n",
           "proposed letter: ",filter_split[[5]],"\n")
    }
    
    if(length(filter_split)!=5){
      stop(method," : wrong specification of \'filter\' \n",
           "if the fourth letter of \'filter\' is \"S\" then it must contains only five letters \n",
           "proposed nomber of letter: ",length(filter_split),"\n")
    }
    
  }else{
    
    if(length(filter_split)==5){
      if(filter_split[[5]] %in% c("3","5","7","9") == FALSE){
        stop(method," : wrong specification of \'filter\' \n",
             "the fifth letter of \'filter\' must correspond to \"3\", \"5\", \"7\", or \"9\" \n",
             "proposed letter: ",filter_split[[5]],"\n")
      }
      filter_split[[5]] <- as.numeric(filter_split[[5]])
    }else if(length(filter_split)==6){
      if(filter_split[[5]] %in% c("1","3","5","7","9") == FALSE){
        stop(method," : wrong specification of \'filter\' \n",
             "the fifth letter of \'filter\' be odd \n",
             "proposed letter: ",filter_split[[5]],"\n")
      }
      
      if(filter_split[[6]] %in% as.character(1:9) == FALSE){
        stop(method," : wrong specification of \'filter\' \n",            
             "the six letter of \'filter\' must correspond to \"1\", \"2\" ... \"9\" \n",
             "proposed letter: ",filter_split[[6]],"\n")
      }
      
      filter_split[[5]] <- as.numeric(filter_split[[6]])+10*as.numeric(filter_split[[5]])
      filter_split[[6]] <- NULL
    }else{
      stop(method," : wrong specification of \'filter\' \n",
           "\'filter\' must contains 5 or 6 letters \n",
           "nb of letters proposed : ",length(filter_split),"\n")
    }
    
  }
  
  
  #### creation of the filters ####  
  if(filter_split[[1]]==2){
    if(filter_split[[4]] == "I"){ # immediate neighborhood
      filter <- matrix(0,nrow=filter_split[[5]],ncol=filter_split[[5]])
      index_1 <- rowSums((which(filter==0,arr.ind=TRUE)-median(1:filter_split[[5]]))^2)<=filter_split[[5]]
      filter[index_1] <- 1
    }
    if(filter_split[[4]]=="G"){ # gaussian
      val_Filter <- dbinom(0:(filter_split[[5]]-1),filter_split[[5]]-1,0.5)
      filter <- matrix(val_Filter,nrow=filter_split[[5]]) %*% matrix(val_Filter,ncol=filter_split[[5]])
    }
    if(filter_split[[4]]=="S"){ # sobel
      if(filter_split[[5]]=="x"){filter <- c(1,2,1) %*% t(c(1,0,-1))}
      if(filter_split[[5]]=="y"){filter <- c(1,0,-1) %*% t(c(1,2,1))}      
    }
    if(filter_split[[4]] == "M"){ # median
      filter <- matrix(1,nrow=filter_split[[5]],ncol=filter_split[[5]])
    }     
  }
  if(filter_split[[1]]==3){
    if(filter_split[[4]] == "I"){
      filter <- array(0,dim=rep(filter_split[[5]],3))     
      index_1 <- rowSums((which(filter==0,arr.ind=TRUE)-median(1:filter_split[[5]]))^2)<=filter_split[[5]]
      filter[index_1] <- 1
    }
    if(filter_split[[4]]=="G"){
      val_Filter <- dbinom(0:(filter_split[[5]]-1),filter_split[[5]]-1,0.5)
      filter <- array(val_Filter,dim=filter_split[[5]]) %o% array(val_Filter,dim=filter_split[[5]]) %o% array(val_Filter,dim=filter_split[[5]])
    }
    if(filter_split[[4]]=="S"){      
      if(filter_split[[5]]=="x"){filter <- c(1,2,1) %o% c(1,0,-1) %o% c(0,1,0)}
      if(filter_split[[5]]=="y"){filter <- c(1,0,-1) %o% c(1,2,1) %o% c(0,1,0)}
      if(filter_split[[5]]=="z"){filter <- c(0,1,0) %o% c(1,2,1) %o% c(1,0,-1)}
    }  
    if(filter_split[[4]] == "M"){
      filter <- array(1,dim=rep(filter_split[[5]],3))  
    }      
  }
  
  return(list(filter=filter,
              filter_split=filter_split)
  )
}

initIndex <- function(object,index,num,hemisphere="both",as.logical,indexNum=NULL,
                      cex.default,pch.default,col.default,filter_default,method){
  
  #### initialization ####
  if(length(index)==1 && is.character(index)){
    index <- list(coords=index)
  }
  if(is.matrix(index)){
    index <- as.data.frame(index)
  }
  if(is.data.frame(index)){
    index <- list(coords=index)
  }
  
  #### parameter ####
  if(is.list(index) && length(index$coords)==1 && is.character(index$coords)){
    
    param_index <- index$coords
    
    if(class(object)=="MRIaggr"){
      initParameter(object,param=param_index,test=TRUE,init=FALSE,accept.coords=FALSE,accept.mask=TRUE,accept.index=FALSE,
                    arg_name="param",long_name="parameters",method="multiplot")        
      index$coords <- selectContrast(object,param=param_index,num=num,hemisphere=hemisphere,coords=TRUE)
    }else{        
      stop(method," : wrong specification of \'index",indexNum,"\' \n",
           "\'index",indexNum,"\' it must be a list containing an element named \"coords\" containing the coordinates of the points to display \n",
           "index",indexNum,"$coords : ",index$coords,"\n")
    }
    
    if(as.logical==TRUE){index$coords[,param_index] <- as.logical(index$coords[,param_index])}
    
    if(is.logical(index$coords[,param_index])==FALSE){
      stop(method," : \'index",indexNum,"$coords\'is not of type logical \n",
           "requested parameter : ",param_index,"\n",
           "type of ",param_index," values : ",paste(is(index$coords[,param_index]),collapse=" "),"\n",
           "to force the conversion to logical set \'as.logical\'= TRUE \n")
    }
    
    index$index <- selectContrast(object,param="index",num=num,hemisphere=hemisphere,format="vector")[index$coords[,param_index]==TRUE]
    index$coords <- index$coords[index$coords[,param_index]==TRUE,c("i","j","k")]             
  }
  
  if(!is.list(index) || "coords" %in% names(index) == FALSE ){
    stop(method," : wrong specification of \'index",indexNum,"\' \n",
         "\'index",indexNum,"\' must be a list containing an element named \"coords\" containing the coordinates of the points to display \n",
         "names(index",indexNum,") : ",paste(names(index),collapse=" "),"\n")
  }
  
  if(class(object)=="MRIaggr"){
    if(length(index$coords) == 1 || ncol(index$coords)!=3 || any(names(index$coords) %in% c("i","j","k") ==FALSE)){
      stop(method," : wrong specification of \'index",indexNum,"\' \n",
           "\'index",indexNum,"$coords\' must have 3 columns named \"i\" \"j\" \"k\" \n",
           "names(index",indexNum,"$coords) : ",paste(names(index),collapse=" "),"\n")
    }
  }else{  
    if(length(index$coords) == 1 || ncol(index$coords)!=3){
      stop(method," : wrong specification of \'index",indexNum,"\' \n",
           "\'index",indexNum,"$coords\' must have 3 columns \n",
           "names(index",indexNum,"$coords) : ",paste(names(index),collapse=" "),"\n")
    }
  }  
  
  if("cex" %in% names(index) == FALSE){index$cex <- cex.default}
  if("pch" %in% names(index) == FALSE){index$pch <- pch.default}
  if("col" %in% names(index) == FALSE){index$col <- col.default}
  
  if("outline" %in% names(index) == TRUE &&  index$outline==TRUE){
    if("filter" %in% names(index) == TRUE){filter <- index$filter}else{filter <- filter_default}
    index$coords <-  pointsOutline(index$coords,filter=filter)
  }
  
  #### export ####
  return(index)
  
}

initNeighborhood <- function(Neighborhood,method){
  
  filter_name <- as.character(Neighborhood)
  valid_names <- c("2D_N4","2D_N8","3D_N4","3D_N6","3D_N8","3D_N10","3D_N18","3D_N26")
  
  if(Neighborhood %in% valid_names ==FALSE){
    stop(method," : wrong specification of  \'Neighborhood\' \n",
         "valid names : \"2D_N4\" \"2D_N8\" \"3D_N4\" \"3D_N6\" \"3D_N8\" \"3D_N10\" \"3D_N18\" \"3D_N26\" \n",
         "proposed name : ",Neighborhood,"\n")
  }
  
  Neighborhood_split <- unlist(strsplit(Neighborhood,split=""))
  
  p.Neighborhood <- as.numeric(Neighborhood_split[[1]])
  n.Neighborhood <- if(length(Neighborhood_split)==5){
    as.numeric(Neighborhood_split[5])
  }else{
    sum(as.numeric(Neighborhood_split[5:6])*c(10,1))
  }
  
  Neighborhood <- matrix(0,nrow=n.Neighborhood,ncol=p.Neighborhood)
  
  Neighborhood[1:4,1:2] <- rbind(c(-1,0),
                                 c(0,-1),
                                 c(1,0),
                                 c(0,1))
  
  if(n.Neighborhood %in% c(8,10,18,26)){
    Neighborhood[5:8,1:2] <- rbind(c(-1,-1),
                                   c(1,1),
                                   c(-1,1),
                                   c(1,-1)
    )
  }
  
  if(n.Neighborhood %in% c(6,10,18,26)){
    row_tempo <-  min(which(rowSums(abs(Neighborhood))==0))
    Neighborhood[seq(row_tempo,row_tempo+1),] <- rbind(c(0,0,1),
                                                       c(0,0,-1)
    )
  }
  
  if(n.Neighborhood %in% c(18,26)){
    row_tempo <-  min(which(rowSums(abs(Neighborhood))==0))
    Neighborhood[seq(row_tempo,row_tempo+7),] <- rbind(c(1,0,1),
                                                       c(0,1,1),
                                                       c(-1,0,1),
                                                       c(0,-1,1),
                                                       c(1,0,-1),
                                                       c(0,1,-1),
                                                       c(-1,0,-1),
                                                       c(0,-1,-1)
    )
  }
  
  if(n.Neighborhood == 26){
    row_tempo <-  min(which(rowSums(abs(Neighborhood))==0))
    Neighborhood[seq(row_tempo,row_tempo+7),] <- rbind(c(1,1,1),
                                                       c(-1,1,1),
                                                       c(-1,-1,1),
                                                       c(1,-1,1),
                                                       c(1,1,-1),
                                                       c(-1,1,-1),
                                                       c(-1,-1,-1),
                                                       c(1,-1,-1)
    )
  }
  
  return(Neighborhood)
}

initWindow <- function(window,filename,path,width,height,unit,res,
                       n.plot,mfrow,xlim,ylim,method){
  
  #### tests ####
  if(!is.null(window) && window %in% c(TRUE,FALSE,"png","eps","svg","pdf") ==FALSE){
    stop(method," : wrong specification of \'window\' \n",
         "valid values : NULL TRUE FALSE \"png\" \"eps\" \"svg\" \"pdf\" \n",
         "proposed value : ",window,"\n")
  }
  
  if(length(filename)!=1 || is.character(filename)==FALSE){
    stop(method," : wrong specification of \'filename\' \n",
         "\'filename\' must be have length 1 and be a character \n",
         "length(filename) : ",length(filename),"\n",
         "is(filename) : ",paste(is(filename),collapse=" "),"\n")
  }
  
  if(length(width)!=1 || is.numeric(width)==FALSE){
    stop(method," : wrong specification of \'filename\' \n",
         "\'width\' must be have length 1 and be a numeric \n",
         "length(width) : ",length(width),"\n",
         "is(width) : ",paste(is(width),collapse=" "),"\n")
  }
  
  if(!is.null(path) && (length(path)!=1 || is.character(path)==FALSE)){
    stop(method," : wrong specification of \'path\' \n",
         "\'path\' must be have length 1 and be a character \n",
         "length(path) : ",length(path),"\n",
         "is(path) : ",paste(is(path),collapse=" "),"\n")
  }
  
  if(!is.null(path) && substr(path,start=nchar(path),stop=nchar(path))!="/"){
    warning(method," : possible bad specification of \'path\' \n",
            "\'path\' should end with a fsep (e.g. \"/\") \n",
            "proposed path : ",path,"\n")
  }
  
  if(unit %in% c("px","in","cm","mm") ==FALSE){
    stop("multiplot : wrong specification of \'unit\' \n",
         "valid values : \"px\" \"in\" \"cm\" \"mm\" \n",
         "requested value : ",unit,"\n")
  }
  
  if(length(res)!=1 || (is.numeric(res)==FALSE && is.na(res)==FALSE)){
    stop(method," : wrong specification of \'filename\' \n",
         "\'height\' must be have length 1 and be a numeric or NA  \n",
         "length(height) : ",length(height),"\n",
         "is(height) : ",paste(is(height),collapse=" "),"\n")
  }
  
  #### initialization ####
  scale <- switch(unit,
                  "px"=1,
                  "in"=90,
                  "cm"=35.43,
                  "mm"=3.543)
  
  if(is.null(mfrow)){
    if(n.plot <= 9){
      mfrow <- c(ceiling(n.plot/ceiling(sqrt(n.plot))),ceiling(sqrt(n.plot)))
    }else{
      mfrow <- c(3,3)
    }              
  }else{
    mfrow <- mfrow
  }
  
  n.graph_par_window <- prod(mfrow)
  
  if(!is.null(xlim)){xlim.plot <- xlim}else{xlim.plot <- NULL}
  if(!is.null(ylim)){ylim.plot <- ylim}else{ylim.plot <- NULL}
  
  #### export ####
  res <- list()
  res$scale <- scale
  res$mfrow <- mfrow
  res$n.graph_par_window <- n.graph_par_window
  res$xlim.plot <- xlim.plot
  res$ylim.plot <- ylim.plot
  
  return(res)  
  
}


initSpace <- function(character){
  
  nchar.character <- lapply(character,function(x){nchar(x)})
  nchar_max.character <- max(unlist(nchar.character))
  
  space.character <- sapply(character,function(x){do.call(paste,c(as.list(rep(" ",nchar_max.character-nchar(x))),sep=""))})
  space.character <- unlist(lapply(space.character,function(x){if(length(x)==0){""}else{x}}))
  return(space.character)
}

initDirPat_constLatex <- function(dir,param,identifier,trace){ 

  ## definition de dir
  if(length(dir)==1 && is.character(dir)){
    dirs_plot <- list.dirs(dir)[-1]
  }else{
    stop("constLatex : wrong specification of \'dir\' \n",
         "\'dir\' must be a single character \n",
         "is(dir) : ",paste(is(dir),collapse=" "),"\n",
         "length(dir) : ",length(dir),"\n") 
  }  
  n.dirs_plot <- length(dirs_plot)
  if(n.dirs_plot==0){
    stop("constLatex : no directory inside the root directory \n",
         "\'dir\' must contain one directory per parameter \n",
         "proposed directory : ",dir,"\n") 
  }  
  
  names_dirs <- unlist(lapply(strsplit(dirs_plot,split="/"),function(x){x[length(x)]}))
  
  if(trace==TRUE){
    cat("%% ",n.dirs_plot," directories found in the root directory \n",sep="")
  }
  
  ## dir restrected to param
  if(!is.null(param)){
    if(any(param %in% names_dirs == FALSE)){
      stop("constLatex : wrong specification of \'param\' \n",
           "available parameters : \"",paste(names_dirs[names_dirs %in% param == FALSE],collapse="\" \""),"\" \n",
           "not available parameters : \"",paste(param[param %in% names_dirs == FALSE],collapse="\" \""),"\" \n")       
    }
    
    dirs_plot <- dirs_plot[sapply(param,function(x){which(names_dirs ==x)})]
    n.dirs_plot <- length(dirs_plot)
    names_dirs <- unlist(lapply(strsplit(dirs_plot,split="/"),function(x){x[length(x)]}))    
    if(trace==TRUE){
      cat("% ",n.dirs_plot," directories will be considered according to \'param\' argument \n",sep="")
    }
  }
  
  ## param
  if(is.null(param)){
    param <- names_dirs
  }
  
  if(length(param) != length(names_dirs)){
    stop("constLatex : \'param\' and \'names_dirs\' must have same length \n",
         " length(param)  : ",length(param)," \n",
         " length(names_dirs) : ",length(names_dirs)," \n",
         " param : ",paste(param,collapse=" ")," \n",
         " names_dirs : ",paste(names_dirs,collapse=" ")," \n")
  }
  

  ## identifiants
  Allidentifier <- unique(unlist(lapply(dirs_plot,function(x){res <- strsplit(list.files(x),split="_",fixed=TRUE);
                                             unique(unlist(lapply(res,"[",1)))})))
    
  if(is.null(identifier)){
    identifier <- Allidentifier  
  }
  n.identifier <- length(identifier)
  
  if(trace==TRUE){
    cat("%% ",length(Allidentifier)," identifiers found in the subdirectories \n",sep="")
    cat("% ",n.identifier," identifiers will be considered according to \'param\' argument \n \n",sep="")
  }
  
  if(n.identifier==0){
    stop("constLatex : no graphics inside the root directory \n",
         "proposed dir : ",dir,"\n",
         "proposed dirs_plot : ",paste(dirs_plot,collapse="\n                     "),"\n")
  }  
  
  #### export ####
  return(list(dirs_plot=dirs_plot,
              n.dirs_plot=n.dirs_plot,
              names_dirs=names_dirs,
              param=param,
              identifier=identifier,
              n.identifier=n.identifier))
  
}

initSection_constLatex <- function(dirs_plot,names_dirs,param,subsection,label){
  
  #### initialisation de subsection avec les noms des dossiers figures
  if(is.null(subsection)){
    subsection <- gsub(pattern="_",replacement=" ",x=names_dirs)
  }
  n.subsection <- length(subsection)
  
  index_subsection <- 1:n.subsection
     
  #### initialisation des legendes des graphiques avec les noms des dossiers figures
 
  if(is.null(label)){
    label <- gsub(pattern="_",replacement=" ",x=names_dirs)
  }
  
  if(is.list(label) == TRUE){
    stop("constLatex : wrong specification of \'label\'\n",
         "\'label\' must be a vector and not a list \n",
         "is(label) : ",paste(is(label),collapse=" ")," \n") 
  }
  
  if(length(param) != length(label)){
    stop("constLatex : \'param\' and \'label\' must have same length \n",
         " length(param)  : ",length(param)," \n",
         " length(label) : ",length(label)," \n") 
  }
  
  #### export ####
  return(list(subsection=subsection,
              index_subsection=index_subsection,
              label=label,
              param=param))
}


initPlot_constLatex <- function(dir,names_dirs,dirs_plot,identifier,table,
                                 plotPerPage){
  
  n.dirs_plot <- length(dirs_plot)
  
  ls.plot <- list() # list of the plot to display for each patient
  ls.newplot <- list() # list of the indicator for the creation of a new figure in latex
  ls.endplot <- list() # list of the indicator for the end of a figure in latex
  ls.legend <- list() # list of the legend files
  ls.slices <- list() # list of the first-last slices for each figure
 
  #### iteration sur les sous dossiers
  for(iter_dir in 1:n.dirs_plot){ 

    split_tempo <- strsplit(list.files(dirs_plot[iter_dir]),split="_",fixed=TRUE)
    
    if(length(split_tempo)==0){ # pas de fichier dans le repertoire 
      ls.legend[[iter_dir]] <- NA
      ls.plot[[iter_dir]] <- NA
      ls.newplot[[iter_dir]] <- NA 
      ls.endplot[[iter_dir]] <- NA 
      ls.slices[[iter_dir]] <- NA 
      next      
    } 
  
    index_plot <- which(unlist(lapply(split_tempo, 
                                      function(x){x[1]==identifier})))      
    n.index_plot <- length(index_plot) # nb de fichiers images du patient
    
    if(n.index_plot==0){ # pas de fichier dans le repertoire correspondants au patient      
      ls.legend[[iter_dir]] <- NA
      ls.plot[[iter_dir]] <- NA
      ls.newplot[[iter_dir]] <- NA 
      ls.endplot[[iter_dir]] <- NA 
      ls.slices[[iter_dir]] <- NA 
      next      
    } 
      
      legend_tempo <- unlist(lapply(split_tempo[index_plot],
                                    function(x){length(grep(pattern="legend",x=x))>0}))
      ls.legend[[iter_dir]] <- legend_tempo
      
      if(length(index_plot)==0){
        stop("constLatex : an legend file was found but with no corresponding image \n",
             "directory : ",dirs_plot[iter_dir]," \n",
             "legend file : \"",list.files(dirs_plot[iter_dir])[legend_tempo>0],"\" \n") 
      }
    
       # test fichier multiplot
      position_slices <- unlist(lapply(split_tempo[index_plot],
                                       function(x){grep(x,pattern="slice",fixed=TRUE)}))
      
      if(length(position_slices)==0){ # si ce n est pas un fichier multiplot
        ls.legend[[iter_dir]] <- rep(FALSE,n.index_plot)
        ls.plot[[iter_dir]] <- paste(tail(strsplit(dir,split="/",fixed=TRUE)[[1]],1),"/",
                                     names_dirs[iter_dir],"/",
                                     list.files(dirs_plot[iter_dir])[index_plot],
                                     sep="")
        ls.newplot[[iter_dir]] <- seq(0,n.index_plot-1) %% plotPerPage == 0
        ls.endplot[[iter_dir]] <- rep(FALSE,n.index_plot)
        ls.endplot[[iter_dir]][intersect(which(ls.newplot[[iter_dir]])-1,1:n.index_plot)] <- TRUE
        ls.endplot[[iter_dir]][n.index_plot] <- TRUE
        ls.slices[[iter_dir]] <- rep("",n.index_plot)        
        next
      }
      
        # extraction des coupes
      slice_tempo <- lapply(1:n.index_plot,
                            function(x){strsplit(split_tempo[[index_plot[x]]][position_slices[x]],
                                                 split="slice",fixed=TRUE)[[1]][[2]]})
      
      # verifie qu il s agit du meme plot
      if(length(index_plot)>1){
        test.plot <- sapply(1:n.index_plot,function(x){
          x_tempo <- gsub(pattern="-legend",replacement="",fixed=TRUE,list.files(dirs_plot[iter_dir])[index_plot[x]])
          gsub(pattern=slice_tempo[[x]],replacement="",fixed=TRUE,x_tempo)
        })
        
        if(any(test.plot!=test.plot[1])){
          stop("constLatex : different types of plot in the same directory \n",
               "directory : ",dirs_plot[iter_dir]," \n",
               "files : \"",list.files(dirs_plot[iter_dir])[index_plot[1]],"\" \n",
               "        \"",paste(list.files(dirs_plot[iter_dir])[index_plot[which(test.plot!=test.plot[1])]],collapse="\"\n        \""),"\"\n")
        }
      }      
      
      # mettre en ordre les figures suivant les coupes (avec les figures en premier et la legend en dernier)
      order_tempo <- order(as.numeric(unlist(lapply(slice_tempo,
                                                    function(x){strsplit(x,split="-",fixed=TRUE)[[1]][1]})))+10^6*legend_tempo)
      index_plot <- index_plot[order_tempo]
      
            
      # store results
      ls.legend[[iter_dir]] <- c(ls.legend[[iter_dir]][legend_tempo==0],ls.legend[[iter_dir]][legend_tempo>0])          
      ls.plot[[iter_dir]] <- paste(tail(strsplit(dir,split="/",fixed=TRUE)[[1]],1),"/",
                                   names_dirs[iter_dir],"/",
                                   list.files(dirs_plot[iter_dir])[index_plot],
                                   sep="")
      ls.newplot[[iter_dir]] <- seq(0,n.index_plot-1) %% plotPerPage == 0
      ls.endplot[[iter_dir]] <- rep(FALSE,n.index_plot)
      ls.endplot[[iter_dir]][intersect(which(ls.newplot[[iter_dir]])-1,1:n.index_plot)] <- TRUE
      ls.endplot[[iter_dir]][n.index_plot] <- TRUE
      
      slice_tempo.first <- unlist(lapply(slice_tempo[order_tempo],function(x){strsplit(x,split="-",fixed=TRUE)[[1]][[1]]}))
      slice_tempo.last <- unlist(lapply(slice_tempo[order_tempo],function(x){strsplit(x,split="-",fixed=TRUE)[[1]][[2]]}))
      
      ls.slices[[iter_dir]] <- c(paste("( slices ",slice_tempo.first[ls.newplot[[iter_dir]]==TRUE],"-",
                                       slice_tempo.last[ls.endplot[[iter_dir]]==TRUE]," ",
                                 if(any(legend_tempo)){"with legend "},")",sep="")
      )
  }

  names(ls.legend) <- names_dirs
  names(ls.newplot) <- names_dirs    
  names(ls.endplot) <- names_dirs    
  names(ls.plot) <- names_dirs
  names(ls.slices) <- names_dirs

  #### data clinique
  if(!is.null(table)){
    ## identification du patient
    ls.table <- lapply(table,function(x){
      table_tempo <- x[x[,1]==identifier,-1,drop=FALSE];
      names(table_tempo) <- gsub(pattern="_",replacement="\\_",x=names(table_tempo))
      names(table_tempo) <- gsub(pattern="%",replacement="\\%",x=names(table_tempo))
      return(table_tempo)
      }
    )
  
  }else{
    ls.table <- NULL
  }

  return(list(ls.legend=ls.legend,
              ls.newplot=ls.newplot,
              ls.endplot=ls.endplot,
              ls.plot=ls.plot,
              ls.slices=ls.slices,
              ls.table=ls.table))
}