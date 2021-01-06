
#################### MetaAnalysisDE ####################
MetaDE<-function(data, clin.data, data.type, resp.type, 
                 mixed=NULL, mix.type=NULL,
                 response,covariate=NULL,ind.method, meta.method, 
                 select.group=NULL , ref.level=NULL, paired=NULL,
                 rth=NULL, REM.type=NULL,
                 asymptotic=NULL, tail='abs',
                 parametric=TRUE, nperm=NULL, seed=12345,...) {
  
  ## Call the packages required 
  
  #library(impute)
  #library(Biobase)
  #library(combinat)
  
  ##Input  
  ##data: list of data matrices, genes on the rows, samples on the columns;
  ##clin.data: list of clinical data frames, samples on the rows, clinical
  ##covariates on the columns;
  ##data.type = c("continuous","discrete");
  ##resp.type = c("twoclass", "multiclass", "continuous", "survival");
  ##response: a character (column name of clin.data) indicating the phenotype 
  ##of interest (survival (2 column names));
  ##ind.method = c("limma","sam","limmaVoom","edgeR","DESeq2","pearsonr",
  ##"spearmanr","logrank");
  ##meta.method=c("maxP","maxP.OC","minP","minP.OC","Fisher","Fisher.OC","AW",
  ##"roP","roP.OC","Stouffer","Stouffer.OC","SR","PR","minMCC","FEM","REM",
  ##"rankProd");
  ##select.group: specify the two group names for comparison;
  ##ref.level: specify the reference level of the group factor;
  ##rth: rth smallest p-value;
  ##REM.type = c("HS","HO", "DL", "SJ", "EB", "RML" );
  ##paired: logical vector of size K indicating whether paired design;
  ##asymptotic: logical whether asymptotic dist should be used for indi part;
  ##tail = c("low", "high", "abs");
  ##parametric: indicate whether the parametric methods is chosen to 
  ##calculate the p-values in meta-analysis
  ##nperm: number of permutations;
  
  K<-length(data) # number of studies
  G<-nrow(data[[1]]) # number of matched genes
  response.list <- covariate.list <- vector('list',K)
  if(!is.null(covariate)) covariate.list <-NULL
  for (k in 1:K) {
    response.list[[k]] <- clin.data[[k]][,response]
    if(!is.null(covariate)) {
      covariate.list[[k]] <-	clin.data[[k]][,covariate]
    }
  }  
  ##verify the correct specification of individual method
  ind.method<- match.arg(ind.method,c("limma","sam","limmaVoom","edgeR",
                                      "DESeq2","pearsonr","spearmanr","logrank"),several.ok=TRUE)
  ##verify the correct specification of meta-analysis method
  meta.method<-match.arg(meta.method,c("maxP","maxP.OC","minP","minP.OC",
                                       "Fisher","Fisher.OC","AW","roP","roP.OC","Stouffer","Stouffer.OC",
                                       "SR","PR","minMCC","FEM","REM","rankProd"),several.ok = TRUE)
  ##check dimensions and size of argument
  check.dim(data,response.list,ind.method,meta.method,paired)
  
  ##check the meta analysis method for the corresponding resp.type & verbose
  check.metamethod(data,response.list,resp.type=resp.type,ind.method,meta.method,
                   rth,REM.type, paired) 
  
  check.parametric(meta.method, parametric) ##check if parametric is ok
  check.tail(meta.method,tail) ##check the tail 
  data<-check.exp(data)  # check the gene names for expression data
  
  ## by resp.type    
  if(resp.type == "twoclass") {
    full_dat<-vector(mode = "list", length = K)
    N <- n <- c()
    for(i in 1:K){
      groupLabel <- as.factor(response.list[[i]])
      groupName=levels(groupLabel)  
      name <- select.group
      l<- groupLabel[groupLabel %in% name]
      l<- droplevels(l)
      if(!is.null(ref.level)) {
        l <- relevel(l,ref=ref.level)
      }
      y<-data[[i]][,groupLabel %in% name]
      full_dat[[i]][[1]] <- y
      full_dat[[i]][[2]] <- l
      nns<- get.sample.label.number(l)
      N<-c(N,nns$N) #sample size per study 
      n<-append(n,nns$n) #sample size per label per study
    }
  } else {
    K<-length(data)
    full_dat<-vector(mode = "list", length = K)
    N <- n <- c()
    for(i in 1:K){
      groupLabel <- as.factor(response.list[[i]])
      y<-data[[i]]
      full_dat[[i]][[1]] <- y
      full_dat[[i]][[2]] <- groupLabel
      nns<- get.sample.label.number(groupLabel)
      N<-c(N,nns$N) #sample size per study 
      n<-append(n,nns$n) #sample size per label per study
    }
  }
  
  ## by meta.method    
  if ("minMCC"%in%meta.method) { 
    K<-length(data)
    dat<-lbl<-list()
    for(i in 1:K){
      l <- as.factor(response.list[[i]])
      if(!is.null(select.group)) {
        dat[[i]]<-data[[i]][,l %in% select.group]
        lbl[[i]]<-l[l %in% select.group]
      } else{
        dat[[i]]<-data[[i]]
        lbl[[i]]<-l
      }
      
      if(!is.null(ref.level)) {
        l <- relevel(l,ref=ref.level)
      }
      
    }    
    if(data.type=="continuous") {
      res<-get.minMCC(dat=dat,lbl=lbl,nperm=nperm)
    } else if (data.type=='discrete') {	  
      for(k in 1:K){
        temp_dat <- cpm(dat[[k]],log=T,prior.count=0.25)
        dat[[k]] <- data.matrix(temp_dat)
      }
      res<-get.minMCC(dat=dat,lbl=lbl,nperm=nperm)
    } 
    
    raw.data<- full_dat 
    res<-list(meta.analysis=res,raw.data=raw.data)
    attr(res$meta.analysis,"nperstudy")<- N
    attr(res$meta.analysis,"nperlabelperstudy")<- n
    attr(res$meta.analysis,"data.type")<- data.type
    attr(res$meta.analysis,"meta.analysis")<- meta.method
    attr(res$meta.analysis,"response.type")<- resp.type
    attr(res$meta.analysis,"group.name")<- select.group 
    attr(res$meta.analysis,"ref.group")<- ref.level
    if(is.null(names(data))) { 
      attr(res$meta.analysis,"dataset.name") <- 
        paste("dataset",1:K,sep="") } else {
          attr(res$meta.analysis,"dataset.name") <- names(data) 
        }    
  }  # end of minMCC
  
  if ("rankProd"%in%meta.method){ #twoclass rank-based
    K<-length(data)
    dat<-lbl<-list()
    for(i in 1:K){
      groupLabel <- as.factor(response.list[[i]])
      groupName=levels(groupLabel)	
      name <- select.group
      l<- groupLabel[groupLabel %in% name]
      l<- droplevels(l)
      if(!is.null(ref.level)) {
        l <- relevel(l,ref=ref.level)
      }
      y<-data[[i]][,groupLabel %in% name]
      dat[[i]]<-y
      lbl[[i]]<-l
    }    
    if(data.type=="continuous") {
      res<-get.RP(dat=dat,lbl=lbl,nperm=nperm)
    } else if (data.type=="discrete") {      
      for(k in 1:K){
        temp_dat <- cpm(dat[[k]],log=T,prior.count=0.25)
        dat[[k]] <- data.matrix(temp_dat)
      }		     
      res <-get.RP(dat=dat,lbl=lbl,nperm=nperm)
    }
    raw.data<- full_dat 
    res<-list(meta.analysis=res,raw.data=raw.data)
    attr(res$meta.analysis,"nperstudy")<- N
    attr(res$meta.analysis,"nperlabelperstudy")<- n
    attr(res$meta.analysis,"data.type")<- data.type
    attr(res$meta.analysis,"response.type")<- resp.type
    attr(res$meta.analysis,"group.name")<- select.group 
    attr(res$meta.analysis,"ref.group")<- ref.level
    if(is.null(names(data))) { 
      attr(res$meta.analysis,"dataset.name") <- 
        paste("dataset",1:K,sep="") } else {
          attr(res$meta.analysis,"dataset.name") <- names(data) 
        }      
  } # end of rankPR
  
  if ("FEM"%in%meta.method||"REM"%in%meta.method){ 
    ## effect size model, two classes allowed only
    if(data.type=="continuous") {
      for(k in 1:K){
        full_dat[[k]][[1]] <- data.matrix(full_dat[[k]][[1]])
      }
      ind.res<-ind.cal.ES(full_dat,paired=paired,nperm=nperm)
    } else if (data.type=="discrete") {       
      for(k in 1:K){
        temp_dat <- cpm(full_dat[[k]][[1]],log=T,prior.count=0.25)
        full_dat[[k]][[1]] <- data.matrix(temp_dat)
      }
      ind.res<-ind.cal.ES(full_dat,paired=paired,nperm=nperm)
    }
    meta.res<-MetaDE.ES(ind.res,meta.method=meta.method,REM.type=REM.type)
    raw.data<- full_dat 
    res<-list(meta.analysis=meta.res,ind.ES=ind.res$ES,
              ind.Var=ind.res$Var,raw.data=raw.data)
    attr(res$meta.analysis,"nperstudy")<- N
    attr(res$meta.analysis,"nperlabelperstudy")<- n
    attr(res$meta.analysis,"data.type")<- data.type
    attr(res$meta.analysis,"response.type")<- resp.type
    attr(res$meta.analysis,"group.name")<- select.group 
    attr(res$meta.analysis,"ref.group")<- ref.level
    if(is.null(names(data))) { 
      attr(res$meta.analysis,"dataset.name") <- 
        paste("dataset",1:K,sep="") } else {
          attr(res$meta.analysis,"dataset.name") <- names(data) 
        }  
  } # end of effect size model  
  
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))==0)  {
    ## the other p-value combination methods
    if(!is.null(mixed)){
      if(mixed==T) {
        data.type <- mix.type
      }
    } 
    ind.res<-Indi.DE.Analysis(data=data, clin.data= clin.data, 
                              data.type=data.type, resp.type=resp.type, 
                              response=response, covariate = covariate,
                              ind.method= ind.method,
                              select.group = select.group,
                              ref.level=ref.level,paired=paired,
                              asymptotic = asymptotic, nperm=nperm,
                              tail=tail, seed=seed)
    if (parametric) {
      ind.res$bp<-NULL
      cat("Parametric method was used instead of permutation\n")
    } else cat("Permutation was used instead of the parametric method\n")
    nm<-length(meta.method)
    meta.res<-MetaDE.pvalue(ind.res,meta.method=meta.method,rth=rth,
                            parametric=parametric)$meta.analysis
    
    if(all(data.type=="continuous")){
      for(k in 1:K){
        full_dat[[k]][[1]] <- data.matrix(full_dat[[k]][[1]]) 	
      }	  	
      raw.data<- full_dat 
    } else if (all(data.type=="discrete")){
      for(k in 1:K){
        temp_dat <- cpm(full_dat[[k]][[1]],log=T,prior.count=0.25)
        full_dat[[k]][[1]] <- data.matrix(temp_dat)
      }
      raw.data<- full_dat   
    } else {     
      discrete.index <- which(data.type=="discrete")
      for(k in 1:K){
        if(k %in% discrete.index) {
          temp_dat <- cpm(full_dat[[k]][[1]],log=T,prior.count=0.25)
          full_dat[[k]][[1]] <- data.matrix(temp_dat)
        } else{
          full_dat[[k]][[1]] <- data.matrix(full_dat[[k]][[1]])
        }	 
      }
      raw.data<- full_dat 
    }   
    
    if(resp.type %in% c("twoclass") ){
      ind.stat<- ind.res$log2FC
    }
    if(resp.type %in% c("multiclass") ){
      ind.stat<- NULL
    } 
    if(resp.type %in% c("continuous") ){
      ind.stat<- ind.res$stat
    }
    if(resp.type %in% c("survival") ){
      ind.stat<- ind.res$stat
    }
    
    res<-list(meta.analysis=meta.res,ind.stat=ind.stat,ind.p=ind.res$p,
              raw.data=raw.data)
    
    attr(res$meta.analysis,"nperstudy")<- N
    attr(res$meta.analysis,"nperlabelperstudy")<- n
    attr(res$meta.analysis,"data.type")<- data.type
    attr(res$meta.analysis,"response.type")<- resp.type
    attr(res$meta.analysis,"group.name")<- select.group 
    attr(res$meta.analysis,"ref.group")<- ref.level
    if(is.null(names(data))) { 
      attr(res$meta.analysis,"dataset.name") <- 
        paste("dataset",1:K,sep="") } else {
          attr(res$meta.analysis,"dataset.name") <- names(data) 
        }   
  } ## end of p-value combination methods
  
  #  if("minMCC"%in%meta.method){
  #		class(res)<-"MetaDE.minMCC"
  #	}else if("FEM"%in%meta.method|"REM"%in%meta.method){
  #		class(res)<-"MetaDE.ES"
  #	}else{
  #		class(res)<-"MetaDE.pvalue"
  #	}
  return(res)
} # end of MetaDE function

##Output 
##stat, pvalue, FDR and AW.weight (if AW method is used)

MetaDE.minMCC<-function(x,nperm=100)
{
  x<-check.exp(x)
  K<-length(x)
  dat<-lbl<-list()
  N<-n<-NULL
  for(i in 1:K){
    dat[[i]]<-x[[i]][[1]]
    lbl[[i]]<-x[[i]][[2]]
    nns<-get.sample.label.number2(lbl[[i]])
    N<-c(N,nns$N)
    n<-c(n,nns$n)
  }
  meta.res<-get.minMCC(dat,lbl,nperm=nperm)
  colnames(meta.res$stat)<-colnames(meta.res$pval)<-
    colnames(meta.res$FDR)<-"minMCC"
  attr(meta.res,"nperstudy")<-N
  attr(meta.res,"nperlabelperstudy")<-n 
  res<-list(meta.analysis=meta.res,raw.data=x)
  attr(res$meta.analysis,"nperstudy")<-N
  res$ind.stat<-NA
  res$ind.p<-NA
  #class(res)<-"MetaDE.minMCC"
  return(res)
}

MetaDE.ES<-function(x,meta.method,REM.type) {
  #meta.method<-match.arg(meta.method,c("FEM","REM"))
  K<-ncol(x$ES)
  n <- attr(x,"nperstudy")
  if(meta.method=="REM"){
    res<-get.REM(x$ES,x$Var,n=n, REM.type=REM.type, pe=x$perm.ES,pv=x$perm.Var)
    tempFDR<-matrix(res$FDR,ncol=1)
    rownames(tempFDR)<-rownames(x$ES)
    colnames(tempFDR)<-meta.method
    meta.res<-list(mu.hat=res$mu.hat,mu.var=res$mu.var,Qval=res$Qval,
                   Qpval=res$Qpval,tau2=res$tau2,zval=res$zval,pval=res$pval,
                   FDR=tempFDR)  
  } else{
    res<-get.FEM(x$ES,x$Var,n=n, x$perm.ES,x$perm.Var)
    tempFDR<-matrix(res$FDR,ncol=1)
    rownames(tempFDR)<-rownames(x$ES)
    colnames(tempFDR)<-meta.method
    meta.res<-list(mu.hat=res$mu.hat,mu.var=res$mu.var,zval=res$zval,
                   pval=res$pval,FDR=tempFDR)	
  }	
  attr(meta.res,"nstudy")<-K
  attr(meta.res,"meta.method")<-meta.method 
  #class(meta.res)<-"MetaDE.ES"
  return(meta.res)
}

get.fisher<-function(p,bp=NULL,miss.tol=0.3) {
  k<-ncol(p)
  pval<-stat<-rep(NA,nrow(p))
  if(!is.null(bp)){
    rnum<-which(apply(p,1,function(x) !any(is.na(x))))
    Ug<-matrix(NA,nrow(p),1)
    Ug[rnum,1]<-as.matrix(rowSums(-2*log(p[rnum,])))
    Ubg<-matrix(rowSums(-2*log(bp)),nrow(p),nrow(bp)/nrow(p))
    pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="high")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=Ug, pval=pval, FDR=qval)
  } else {
    rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<miss.tol)
    pval[rnum]<-apply(p[rnum,],1,function(x) {
      ## calculate Fisher p-value 
      pchisq(-2*sum(log(x),na.rm=T),lower.tail=FALSE,2*sum(!is.na(x)) ) 
    })
    qval<-p.adjust(pval,method="BH")
    stat[rnum]<-apply(p[rnum,],1,function(x)-2*sum(log(x),na.rm=T))
    res<-list(stat=stat,pval=pval,FDR=qval)
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}

get.fisher.OC<-function(p,bp=NULL) {
  K<-ncol(p)
  rnum<-which(apply(p,1,function(x) !any(is.na(x))))
  pval<-rep(NA,nrow(p))
  Ug<-matrix(NA,nrow(p),1)
  Ug[rnum,1]<-pmax(rowSums(-log(p))[rnum],rowSums(-log(1-p))[rnum])
  if (is.null(bp)) {
    warning("there're no parametric results for Fisher's OC method,
                  we will use simulation to estimate the p values")
    bp<-matrix(runif(500*nrow(p)*K),500*nrow(p),K) 
  }
  Ubg<-matrix(pmax(rowSums(-log(bp)),rowSums(-log(1-bp))),
              nrow(p),nrow(bp)/nrow(p))
  pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="high")
  qval<-p.adjust(pval,method="BH")
  res<-list(stat=Ug,pval=pval,FDR=qval)
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}

get.maxP<-function(p,bp=NULL,miss.tol=0.3) {
  k<-ncol(p)
  pval<-stat<-rep(NA,nrow(p))
  if(!is.null(bp)){
    rnum<-which(apply(p,1,function(x) !any(is.na(x))))
    Ug<-matrix(NA,nrow(p),1)
    Ug[rnum,1]<-as.matrix(apply(p[rnum,],1,max))
    Ubg<-matrix(apply(bp,1,max),nrow(p),nrow(bp)/nrow(p))
    pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=Ug, pval=pval, FDR=qval)
  }else{
    rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<=miss.tol)
    pval[rnum]<-apply(p[rnum,],1,function(x) {
      ## calculate maxP p-value 
      pbeta(max(x,na.rm=T),sum(!is.na(x)),1)
    })
    stat[rnum]<-apply(p[rnum,],1,function(x) max(x,na.rm=T))
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=stat,pval=pval,FDR=qval)
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}


get.maxP.OC<-function(p,bp=NULL,miss.tol=0.3) {
  k<-ncol(p)
  pval<-stat<-rep(NA,nrow(p))
  if(!is.null(bp)){
    rnum<-which(apply(p,1,function(x) !any(is.na(x))))
    Ug<-matrix(NA,nrow(p),1)
    Ug[rnum,1]<-apply(p[rnum,],1,function(x)min(sort(x)[k],sort(1-x)[k]))
    Ubg<-matrix(apply(bp,1,function(x)min(sort(x)[k],sort(1-x)[k])),
                nrow(p),nrow(bp)/nrow(p))
    pval[rnum]<-perm.p(Ug[rnum,],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=c(Ug), pval=pval, FDR=qval)
  }else{
    rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<=miss.tol)
    stat[rnum]<-apply(p[rnum,],1,function(x) {
      min(sort(x)[sum(!is.na(x))],sort(1-x)[sum(!is.na(1-x))])
    })
    k0<-apply(p[rnum,],1,function(x)sum(!is.na(x)))
    pval[rnum]<-ifelse(stat[rnum]>0.5,2*stat[rnum]^k0-(2*stat[rnum]-1)^k0,
                       2*stat[rnum]^k0)	
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=stat,pval=pval,FDR=qval)
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}

get.minP<-function(p,bp=NULL,miss.tol=0.3){
  k<-ncol(p)
  pval<-stat<-rep(NA,nrow(p))
  if (!is.null(bp)){
    rnum<-which(apply(p,1,function(x) !any(is.na(x))))
    Ug<-matrix(NA,nrow(p),1)
    Ug[rnum,1]<-apply(p[rnum,],1,min)
    Ubg<-matrix(apply(bp,1,min),nrow(p),nrow(bp)/nrow(p))
    pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=c(Ug), pval=pval, FDR=qval)
  }else{
    rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<=miss.tol)
    pval[rnum]<-apply(p[rnum,],1,function(x) {
      ## calculate minP p-value 
      pbeta(min(x,na.rm=T),1,sum(!is.na(x)))
    })
    stat[rnum]<-apply(p[rnum,],1,function(x)min(x,na.rm=T))
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=stat,pval=pval,FDR=qval)
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}

get.minP.OC<-function(p,bp=NULL,miss.tol = 0.3) {
  K<-ncol(p)
  pval<-stat<-rep(NA,nrow(p))
  if (!is.null(bp)){
    rnum<-which(apply(p,1,function(x) !any(is.na(x))))
    Ug<-matrix(NA,nrow(p),1)
    Ug[rnum,1]<-as.matrix(apply(cbind(p,1-p)[rnum,],1,min))
    Ubg<-matrix(apply(cbind(bp,1-bp),1,min),nrow(p),nrow(bp)/nrow(p))
    pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=c(Ug),pval=pval,FDR=qval)
  }else{
    rnum<-which(apply(p,1,function(x) sum(is.na(x))/K)<=miss.tol)
    stat[rnum]<-apply(cbind(p,1-p)[rnum,],1,min,na.rm=T)
    K0<-apply(p[rnum,],1,function(x)sum(!is.na(x)))
    pval[rnum]<-ifelse(stat[rnum]>=0.5,1,1-(1-2*stat[rnum])^K0)
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=stat,pval=pval,FDR=qval)
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}

get.roP<-function(p,bp=NULL,rth) {
  k<-ncol(p)
  rnum<-which(apply(p,1,function(x) !any(is.na(x))))
  pval<-stat<-rep(NA,nrow(p))
  if(!is.null(bp)){
    p<-t(apply(p,1, sort,na.last = T))
    bp<-t(apply(bp,1,sort,na.last = T))
    Ug<-matrix(NA,nrow(p),1)
    Ug[rnum,1]<-p[rnum,rth]	
    Ubg<-matrix(bp[,rth],nrow(p),nrow(bp)/nrow(p))
    pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=Ug,pval=pval, FDR=qval)
  }else{
    pval[rnum]<-apply(p[rnum,],1,function(x) {
      ## calculate rOP p-value 
      pbeta(x[order(x)][rth],rth,k-rth+1)
    }) 
    stat[rnum]<-apply(p[rnum,],1,function(x)x[order(x)][rth])			
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=stat,pval=pval,FDR=qval)	
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}

get.roP.OC<-function(p,bp=NULL,rth) {
  k<-ncol(p)
  rnum<-which(apply(p,1,function(x) !any(is.na(x))))
  pval<-stat<-rep(NA,nrow(p))
  if(!is.null(bp)){
    Ug<-matrix(NA,nrow(p),1)
    Ug[rnum,1]<-apply(p[rnum,],1,function(x)min(sort(x)[rth],sort(1-x)[rth]))
    Ubg<-matrix(apply(bp,1,function(x)min(sort(x)[rth],sort(1-x)[rth])),nrow(p),
                nrow(bp)/nrow(p))
    pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=c(Ug), pval=pval, FDR=qval)
  }else{
    stat[rnum]<-apply(p[rnum,],1,function(x)min(sort(x)[rth],sort(1-x)[rth]))
    pval[rnum]<-sapply(stat[rnum],function(x)CDF.rop(x,rth,k))	
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=stat,pval=pval,FDR=qval)
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}


get.Stouff<-function(p,bp=NULL,miss.tol=0.3){
  k<-ncol(p)
  pval<-stat<-rep(NA,ncol(p))
  if(!is.null(bp)){
    rnum<-which(apply(p,1,function(x) !any(is.na(x))))
    Ug<-matrix(NA,nrow(p),1)
    Ug[rnum,1]<-apply(p[rnum,],1,function(x)sum(qnorm(x))/sqrt(k))
    Ubg<-matrix(apply(bp,1,function(x)sum(qnorm(x))/sqrt(k)),nrow(p),
                nrow(bp)/nrow(p))
    pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="abs")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=c(Ug), pval=pval, FDR=qval)
  }else{
    rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<=miss.tol)
    pval[rnum]<-apply(p[rnum,],1,function(x) {
      ## calculate Stouffer p-value 
      2*pnorm(abs(sum(qnorm(x),na.rm=T)/sqrt(sum(!is.na(x)))),lower.tail=FALSE  ) 
    })
    stat[rnum]<-apply(p[rnum,],1,function(x) {
      sum(qnorm(x),na.rm=T)/sqrt(sum(!is.na(x)))
    })
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=stat,pval=pval,FDR=qval)
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}

get.Stouff.OC<-function(p,bp=NULL){
  k<-ncol(p)
  rnum<-which(apply(p,1,function(x) !any(is.na(x))))
  pval<-rep(NA,nrow(p),1)
  Ug<-UL<-UR<-matrix(NA,nrow(p),1)
  UL[rnum,1]<-apply(p[rnum,],1,function(x)sum(qnorm(x))/sqrt(k))
  UR[rnum,1]<-apply((1-p)[rnum,],1,function(x)sum(qnorm(x))/sqrt(k))
  Ug[rnum,1]<-pmax(UL,UR)[rnum]
  if(!is.null(bp)){
    UbL<-as.matrix(apply(bp,1,function(x)sum(qnorm(x))/sqrt(k)))
    UbR<-as.matrix(apply(1-bp,1,function(x)sum(qnorm(x))/sqrt(k)))
    Ubg<-pmax(UbL,UbR)
    pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="high")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=c(Ug), pval=pval, FDR=qval)
  }else{
    pval[rnum]<-sapply(Ug[rnum,1],function(x)ifelse(x>0,2*(1-pnorm(x)),0))
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=c(Ug),pval=pval,FDR=qval)
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}

########## new version AW-fisher ###############

get.AW <- function(p.values, AW.type= "original", bp=NULL, logInput = FALSE, 
                   log=FALSE, weight.matrix=TRUE) {
  #gene.names<-rownames(p.values) 
  if(NCOL(p.values) == 1) p.values= t(p.values)
  n = NCOL(p.values)
  
  if(logInput){
    pCumLog = -p.values
  } else {
    pCumLog= -log(p.values)		
  }
  
  if(weight.matrix) {
    orderMatrix = t(apply(p.values, 1, order))
    for(i in 1:NROW(p.values)) {
      pCumLog[i,] = cumsum(pCumLog[i,orderMatrix[i,]])
    }
  } else {
    for(i in 1:NROW(p.values)) {
      pCumLog[i,] = cumsum(pCumLog[i,order(-pCumLog[i,])])
    }
  }
  
  sumWeight = rep(1, NROW(p.values))
  if( AW.type == "original") {
    bestStat = -pchisq(2*pCumLog[,1], 2, lower.tail=F, log=T)
    for(i in 2:n) {
      statNew = -pchisq(2*pCumLog[,i], 2*i, lower.tail=F, log=T)
      sumWeight[statNew > bestStat] = i
      bestStat = pmax(bestStat, statNew)
    }
  } else if(AW.type == "uncond") {
    bestStat = -pbeta(exp(-pCumLog[,1]), 1, n, log=T)
    statNew =  -pchisq(2*pCumLog[,n], 2*n, lower.tail=F, log=T)
    sumWeight[statNew > bestStat] = n
    bestStat = pmax(bestStat, statNew)
    
    if(n>2) for(i in 2:(n-1)) { 
      statNew = uncondTopDist(2*pCumLog[,i],n,i)
      sumWeight[statNew > bestStat] = i
      bestStat = pmax(bestStat, statNew)
    }
    
  } else if ( AW.type == "cond") {
    bestStat = -pchisq(2*pCumLog[,n], 2*n, lower.tail=F, log=T)
    for(i in 1:(n-1)) {
      statNew = -pchisq(2*(pCumLog[,i]-i*(pCumLog[,i+1]-pCumLog[,i])), 2*i, 
                        lower.tail=F, log=T)
      sumWeight[statNew > bestStat] = i
      bestStat = pmax(bestStat, statNew)
    }
  } else {
    stop("Incorrect method!")
  }
  
  if(weight.matrix) {
    for(i in 1:NROW(orderMatrix)) {
      orderMatrix[i,orderMatrix[i,1:sumWeight[i]]] = 0
    }
    orderMatrix[orderMatrix != 0] = 1
    
    pval=aw.fisher.stat(bestStat, n,  AW.type, log)
    qval= p.adjust(pval,method='BH')
    res <- list(stat=bestStat, pval=pval, FDR=qval, AW.weight=1-orderMatrix, 
                sum.weight = sumWeight)
  } else {  	
    pval=aw.fisher.stat(bestStat, n, method, log)
    qval= p.adjust(pval,method='BH')
    res <- list(stat=bestStat, pval=pval, FDR=qval,  sum.weight = sumWeight)
  }
  return(res)  
}

### Other Internal Functions #############
aw.fisher.stat <- function(pstat, n, method, log=FALSE) {
  index = match(n, awFisherData[["nList"]])
  #    if(!is.na(index)) {
  #        quant = awFisherData[[method]][index,]
  #    } else {
  ### smooth estimation over n
  numN = NCOL(awFisherData[[method]])
  quant = rep(0, numN)
  for(i in 1:numN) {
    f = splinefun(c(1, awFisherData[["nList"]]),c(awFisherData[["logPTarget"]][i], 
                                                  awFisherData[[method]][,i]))
    quant[i] = f(n)
  }
  #    }
  
  ##### Estimating ###########
  f = splinefun(c(0,quant), c(0,awFisherData[["logPTarget"]]), method="monoH.FC")
  if(log)  -f(pstat)  else exp(-f(pstat))
}


preDefList = exp(log(501)*(1:50)/50)-1
uncondTopDist <- function(stat, n, which) {
  estimates = pvalueTop(A=preDefList, n=n, m=which)
  f = splinefun(c(0,preDefList), c(0,estimates), method="monoH.FC")
  f(stat)
}

intFunc <- function(p,A,n,m) {
  pchisq(A+2*m*log(p), 2*m, lower.tail=F)*dbeta(p, m+1, n-m)
}
pvalueTop <- Vectorize(function(A,n,m) {
  pLim = exp(-A/(2*m))
  -log(integrate(intFunc, lower=pLim,upper=1,A = A, n=n,m=m, stop.on.error=F,
                 rel.tol=1e-3)$value 
       + pbeta(pLim, m+1, n-m))}
)

############################

cal.MCC<-function(dt1,dt2,l1,l2) {
  l1<-unclass(factor(l1))
  l2<-unclass(factor(l2))
  K<-nlevels(l1)
  n1<-table(factor(l1,levels=unique(l1)))
  n2<-table(factor(l2,levels=unique(l2)))
  ind1<-diag(rep(1,length(n1)))[rep(1:nlevels(l1),n1),] 
  ind2<-diag(rep(1,length(n2)))[rep(1:nlevels(l2),n2),]
  xk.<-dt1%*%ind1%*%diag(1/n1)
  yk.<-dt2%*%ind2%*%diag(1/n2)
  x..<-rowMeans(xk.)
  y..<-rowMeans(yk.)
  sxk.yk.<-rowSums(xk.*yk.)
  num<-1/K*sxk.yk.-x..*y..
  sumx2<-dt1^2%*%ind1
  sumy2<-dt2^2%*%ind2
  vx<-1/K*rowSums((sumx2-xk.^2)%*%diag(1/(n1-1)))-x..^2
  vy<-1/K*rowSums((sumy2-yk.^2)%*%diag(1/(n2-1)))-y..^2
  den<-sqrt(vx*vy)
  r<-num/den
  return(r)
}

cal.minMCC<-function(dat,lbl) {
  K<-length(dat)  
  if (K==2) min.MCC<-cal.MCC(dat[[1]],dat[[2]],
                             factor(lbl[[1]]),factor(lbl[[2]]))
  else {
    allcomb<-combn(1:K,2) 
    pair.mcc<-NULL 
    for (i in 1:ncol(allcomb)) {
      dt1<-dat[[allcomb[1,i]]] 
      dt2<-dat[[allcomb[2,i]]]
      l1<-lbl[[allcomb[1,i]]]
      l2<-lbl[[allcomb[2,i]]]
      pair.mcc<-cbind(pair.mcc,cal.MCC(dt1,dt2,factor(l1),factor(l2)))
    }
    row.names(pair.mcc)<-row.names(dat[[1]])
    min.MCC<-apply(pair.mcc,1,min) 
  }
  return(min.MCC)
}

get.minMCC<-function(dat,lbl,nperm) {
  gene.names<-row.names(dat[[1]]) 
  rnum<-unlist(lapply(dat,function(x) which(apply(x,1,function(x) 
    any(!is.na(x))))))
  perm.mcc.perm<-function(dat,lbl) {
    perm.d<-list()
    for (i in 1:length(dat)) {
      perm.d[[i]]<-dat[[i]][,sample(1:ncol(dat[[i]]))]
    }
    perm.r<-cal.minMCC(perm.d,lbl)
    return(perm.r)
  }
  Ug<-cal.minMCC(dat,lbl)
  if (!is.null(nperm)) {
    Ubg<-replicate(nperm,perm.mcc.perm(dat=dat,lbl=lbl))
  } else {
    stop("there're no parametric results for MCC statistic,you need to 
          specify a number to nperm")   
  }
  pval<-qval<-matrix(NA,nrow(dat[[1]]),1)
  pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum],tail="high")
  qval[rnum]<-p.adjust(pval[rnum],method="BH")
  names(Ug)<-row.names(pval)<-row.names(qval)<-gene.names
  res<-list(stat=as.matrix(Ug),pval=as.matrix(pval),FDR=as.matrix(qval))
  attr(res,"nstudy")<-length(dat)
  attr(res,"meta.method")<-"minMCC"
  return(res)
}

perm.p<-function(stat,perm,tail) {
  G<-length(stat)
  B<-length(perm)/G
  if(tail=="low"){
    r = rank(c(stat, as.vector(perm)),ties.method="max")[1:G]
    r2 = rank(c(stat),ties.method="max")
    r3 = r - r2
    p = r3/(B*G)
  }
  if(tail=="high"){
    r = rank(c(stat, as.vector(perm)),ties.method="min")[1:G]
    r2 = rank(stat,ties.method="max")
    r3 = r - r2
    p = 1-r3/(B*G)
  }
  if(tail=="abs"){
    r = rank(c(abs(stat), abs(as.vector(perm))),ties.method="min")[1:G]
    r2 = rank(c(abs(stat)),ties.method="max")
    r3 = r - r2
    p = 1-r3/(B*G)
  }
  p[p==0]<-1e-20
  p[p==1]<-1-1e-10
  return(p)
}


emp<-function(mat,tail) {
  B<-ncol(mat)
  G<-nrow(mat)
  if (tail=="low"){
    s<-matrix(rank(mat,ties.method="max"),G,B)
    p<-s/G/B}
  if (tail=="high"){
    s<-matrix(rank(mat,ties.method="min"),G,B)
    p<-1-(s-1)/G/B}
  if (tail=="abs"){
    s<-matrix(rank(abs(mat),ties.method="min"),G,B)
    p<-1-(s-1)/G/B}
  p[p==0]<-1e-20
  p[p==1]<-1-1e-10
  return(p)
}

CDF.rop<-function(z,r,n){
  require(combinat)
  if (r>=.5*(n+1)){
    pval<-ifelse(z>=0.5&z<1,
                 1-sum(sapply((n-r+1):(r-1),function(y) 
                   sum(sapply((n-r+1):(n-y),function(x)
                     dmnom(c(y,x,n-y-x),n,c(1-z,1-z,2*z-1))))))				
                 ,2*(1-pbinom(r-1,n,z)))}
  else{
    pval<-ifelse(z>=0&z<=0.5,
                 1-sum(sapply((0):(r-1),function(y)sum(sapply(0:(r-1),function(x)
                   dmnom(c(y,x,n-y-x),n,c(z,z,1-2*z)))))),				
                 1)}		
  return(pval)
}

get.SR<-function(p,bp=NULL){
  k<-ncol(p)
  rnum<-which(apply(p,1,function(x) !any(is.na(x))))
  pval<-Ug<-rep(NA,nrow(p))
  Ug[rnum]<-rowSums(apply(p,2,rank))[rnum]
  if(!is.null(bp)){
    nperm<-nrow(bp)/nrow(p)
    Ubg<-matrix(NA,nrow(p),nperm)
    for(j in 1:nperm){
      Ubg[rnum,j]<-rowSums(apply(bp[((j-1)*nrow(p)+1):(j*nrow(p)),],
                                 2,rank))[rnum]
    }
    pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")
  } else{
    nperm=500
    bp<-matrix(runif(500*nrow(p)*k),500*nrow(p),k)
    Ubg<-matrix(NA,nrow(p),nperm)
    for(j in 1:nperm){
      Ubg[rnum,j]<-rowSums(apply(bp[((j-1)*nrow(p)+1):(j*nrow(p)),],
                                 2,rank))[rnum]
    }
    pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")}
  res<-list(stat=Ug, pval=pval, FDR=qval)	
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}


get.PR<-function(p,bp=NULL){
  rowProds<-function (x, ...){
    s <- (x == 0)
    s <- rowSums(s,na.rm=T)
    ok <- (s == 0)
    rm(s)
    x <- x[ok, , drop = FALSE]
    y <- vector(mode(x), nrow(x))
    s <- (x < 0)
    s <- rowSums(s,na.rm=T)
    s <- (s%%2)
    s <- c(+1, -1)[s + 1]
    x <- abs(x)
    x <- log(x)
    x <- rowSums(x, ...)
    x <- exp(x)
    x <- s * x
    y[ok] <- x
    rm(ok, s, x)
    y
  }
  k<-ncol(p)
  pval<-Ug<-rep(NA,nrow(p))
  rnum<-which(apply(p,1,function(x) !any(is.na(x))))
  Ug[rnum]<-rowProds(apply(p,2,rank),na.rm=T)[rnum]
  if(!is.null(bp)){
    nperm<-nrow(bp)/nrow(p)
    Ubg<-matrix(NA,nrow(p),nperm)
    for(j in 1:nperm){
      Ubg[rnum,j]<-rowProds(apply(bp[((j-1)*nrow(p)+1):(j*nrow(p)),],2,rank),
                            na.rm=T)[rnum]
    }
    pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")
  }else{
    nperm=500
    bp<-matrix(runif(500*nrow(p)*k),500*nrow(p),k)
    Ubg<-matrix(NA,nrow(p),nperm)
    for(j in 1:nperm){
      Ubg[rnum,j]<-rowProds(apply(bp[((j-1)*nrow(p)+1):(j*nrow(p)),],2,rank),
                            na.rm=T)[rnum]
    }
    pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")}
  res<-list(stat=Ug, pval=pval, FDR=qval)	
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}


cal.ES<-function(y,l,paired=FALSE){
  l<-unclass(factor(l))
  n<-table(factor(l))
  if(paired){
    if (n[1]!=n[2]) {
      stop("The study is not paired design")
    }
    DM<-y[,l==2] ## disease is label = 2
    CM<-y[,l==1] ## control is label = 1, originally ref level
    ydiff<-DM-CM
    den<-sqrt(1/(n[1]-1)*(rowSums(ydiff^2))-1/(n[1]^2-n[1])*(rowSums(ydiff))^2)
    t<-rowMeans(ydiff)/(den/sqrt(n[1]))
    rnum<-n[1]*rowSums(DM*CM)-rowSums(DM)*rowSums(CM)
    rden<-sqrt((n[1]*rowSums(DM^2)-(rowSums(DM))^2)*(n[1]*rowSums(CM^2)-
                                                       (rowSums(CM))^2))
    r<-rnum/rden
    d<-t*sqrt(2*(1-r)/n[1])
    m<-n[1]-1
    cm = gamma(min(m/2, 100))/(sqrt(m/2) * gamma(min((m - 1)/2,100) ))
    dprime=cm*d
    vard=(2*(1-r)/n[1])*((n[1]-1)/(n[1]-3))*
      (1+n[1]*dprime^2/(2*(1-r)))-dprime^2/cm^2
    vardprime=cm^2*vard
  }else{
    ind<-diag(rep(1,length(n)))[l,]
    ym<-y%*%ind%*%diag(1/n) 
    ntilde<-1/sum(1/n)
    m=sum(n)-2
    cm = gamma(min(m/2,100))/(sqrt(m/2) * gamma(min((m - 1)/2,100) ))
    s<-sqrt((1/(sum(n)-2)*((y^2%*%ind)%*%diag(1/(n-1))-
                             ym^2%*%diag(n/(n-1)))%*%(n-1)))
    d<-(ym[,2]-ym[,1])/s
    dprime=d-3*d/(4*(sum(n)-2)-1)
    terme1=1/ntilde
    vard = terme1 + d^2 * (terme1 * ntilde - 1/cm^2)
    vardprime=sum(1/n)+dprime^2/(2*sum(n))			
  }
  result = cbind( dprime, vardprime)
  colnames(result) = c( "dprime", "vardprime")
  rownames(result)<-rownames(y)
  result
}

get.ES<-function(x,paired){   
  K<-length(x) 
  ES.m<-Var.m<- N<-n<-NULL  
  for (k in 1:K){
    y<-x[[k]][[1]]
    l<-x[[k]][[2]]
    temp<-cal.ES(y,l,paired[k])
    ES.m<-cbind(ES.m,temp[,"dprime"])
    Var.m<-cbind(Var.m,temp[,"vardprime"])
    N<-c(N,length(l))
    n<-c(n,table(l))
  }
  rownames(ES.m)<-rownames(y)
  rownames(Var.m)<-rownames(y)
  colnames(ES.m)<-paste("study",1:K,sep="")
  colnames(Var.m)<-paste("study",1:K,sep="")
  res<-list(ES=ES.m,Var=Var.m)
  attr(res,"nperstudy")<-N
  attr(res,"nperlabelperstudy")<-n 
  return(res)	
}

ind.cal.ES<-function(x,paired,nperm=NULL){
  K<-length(x)
  res<-get.ES(x,paired=paired)
  if(!is.null(nperm)){
    perm.ES<-perm.Var<-NULL
    for(i in 1:nperm){
      for(k in 1:K){
        x[[k]][[2]]<-perm.lab(x[[k]][[2]],paired[k])
      }
      tempRes<-get.ES(x,paired=paired)
      perm.ES<-rbind(perm.ES,tempRes$ES)
      perm.Var<-rbind(perm.Var,tempRes$Var)
    }
  }else{
    perm.ES<-perm.Var<-NULL
  }
  if(is.null(names(x))){colnames(res$ES)<-colnames(res$Var)<-paste("dataset",
                                                                   1:K,sep="")
  }else{colnames(res$ES)<-colnames(res$Var)<-names(x)}
  result<-list(ES=res$ES,Var=res$Var,perm.ES=perm.ES,perm.Var=perm.Var)
  attr(result,"nperstudy")<-attr(res,"nperstudy")
  attr(result,"nperlabelperstudy")<-attr(res,"nperlabelperstudy") 
  return(result)
}

get.tau2 <- function(em, vm, k,  n, REM.type, threshold=10^-5, maxiter = 100) {  
  ## em, vm are all matrices of G rows and K columns;
  ## n is a vector of length K;
  ## em: observed effect size, vm: variance of effect size, n: sample size;
  ## k: number of studies;
  if (REM.type == "HS") {
    wt <- matrix(rep(n,nrow(em)),byrow=T,nrow=nrow(em),ncol=length(n))
    n_bar <- mean(n)
    n_sum <- sum(n)
    em_mean <- rowMeans(em)
    num <- rowSums(wt*(em-em_mean)^2)
    denom <- n_sum
    
    tau2 <- num/denom  - ((n_bar-1)/(n_bar-3)) * (4/n_bar) * (1+em_mean^2/8)
    res <- list(wt,tau2)
    names(res) <- c('weight','tau2')
    return(res)
  }
  if (REM.type == "HO"){
    wt <- 1/vm
    em_mean <- rowMeans(em)
    s1 <- rowSums((em - em_mean)^2)
    temp <- s1/(k-1) - rowMeans(vm)
    tau2 <-  pmax(temp,0)
    res <- list(wt,tau2)
    names(res) <- c('weight','tau2')
    return(res)
  }
  if (REM.type == "DL"){
    wt <- 1/vm
    theta <- rowSums(wt*em) / rowSums(wt)
    num <- rowSums(wt*(em - theta)^2) - (k-1)
    s1 <- rowSums(wt) 
    denom <- s1 - rowSums(wt^2)/ s1
    temp <- num  / denom
    tau2 <-  pmax(temp,0)
    res <- list(wt,tau2)
    names(res) <- c('weight','tau2')
    return(res)
  }
  if (REM.type == "SJ"){
    em_mean <- rowMeans(em)
    tau0 <-rowMeans((em - em_mean)^2)
    wt <- vm/tau0 + 1
    num <- rowSums(em / wt)
    denom <- rowSums(1/wt)
    theta <- num  / denom
    tau2 <- (1/(k-1))*rowSums((em-theta)^2/wt)
    res <- list(wt,tau2)
    names(res) <- c('weight','tau2')
    return(res)
  }
  
  if (REM.type == "EB"){
    tau2<- compute.tau2.EB(em , vm, 1/vm ,k) 
    iter <- 1
    change <- threshold + 1
    while (max(change) > threshold) {
      iter <- iter + 1
      print(iter)
      tau2.old <- tau2 
      wt <- compute.wt.EB(vm, tau2)
      if (any(is.infinite(wt))) 
        stop("Division by zero when computing the inverse variance weights.")
      tau2 <- compute.tau2.EB(em, vm, wt, k) 
      change <- abs(tau2.old - tau2)
      if (iter > maxiter) {
        break
      }
    }
    wt <- compute.wt.EB(vm, tau2)
    res <- list(wt,tau2)
    names(res) <- c('weight','tau2')
    return(res)
  }
  
  if (REM.type == "RML"){
    tau2<- compute.tau2.RML(em , vm, 1/vm ,k) 
    iter <- 1
    change <- threshold + 1
    while (max(change) > threshold) {
      iter <- iter + 1
      print(iter)
      tau2.old <- tau2 
      wt <- compute.wt.RML(vm, tau2)
      if (any(is.infinite(wt))) 
        stop("Division by zero when computing the inverse variance weights.")
      tau2 <- compute.tau2.RML(em, vm, wt, k) 
      change <- abs(tau2.old - tau2)
      if (iter > maxiter) {
        break
      }
    }
    wt <- compute.wt.RML(vm, tau2)
    res <- list(wt,tau2)
    names(res) <- c('weight','tau2')
    return(res)
  }  
}

compute.tau2.EB <- function(em, vm, wt, k) {
  num <- rowSums(em / wt)
  denom <- rowSums(1/wt)
  theta <- num  / denom
  temp <- rowSums(wt^2*( (k/(k-1)) * (em - theta)^2 - vm)) / rowSums(wt^2)
  tau2 <- pmax(0,temp)
  return(tau2)
}

compute.wt.EB <- function(vm, tau2) {
  wt <- 1/(vm + tau2)
  return(wt)
}

compute.tau2.RML <- function(em, vm, wt, k) {
  num <- rowSums(wt* em)
  denom <- rowSums(wt)
  theta <- num  / denom
  #temp <- rowSums(wt^2*( ( (em - theta)^2 + 1 )/(rowSums(wt) - vm)) ) /  
  #rowSums(wt^2)
  temp <- rowSums(wt^2*(  (em - theta)^2  - vm) ) /rowSums(wt^2) + 1/rowSums(wt)
  tau2 <- pmax(0,temp)
  return(tau2)
}

compute.wt.RML <- function(vm, tau2) {
  wt <- 1/(vm + tau2)
  return(wt)
}


get.Q<-function(em, wt){
  temp1 <- wt * em
  mu.hat <- rowSums(temp1)/rowSums(wt)
  Q <- rowSums(wt * (em - mu.hat)^2)
  return(Q)
}

get.FEM<-function(em,vm,n, pe=NULL,pv=NULL){
  wt<-1/vm
  mu.hat<-rowSums(wt*em)/rowSums(wt)
  mu.var<-1/rowSums(wt)
  z.score<-mu.hat/sqrt(mu.var)
  if(!is.null(pe)&!is.null(pv)){
    rnum<-which(apply(em,1,function(x) !any(is.na(x))))
    Z0<-matrix(get.REM2(pe,pv, n=n, REM.type="HO")$zval,nrow(em),
               nrow(pe)/nrow(em))
    z.p<-rep(NA,nrow(em))
    z.p[rnum]<-perm.p(z.score[rnum],Z0[rnum,],"abs")
  }else{
    z.p<-2*(1-pnorm(abs(z.score)))
  }
  qval<-p.adjust(z.p,method="BH")
  res<-list(mu.hat=mu.hat,mu.var=mu.var,zval=z.score,pval=z.p,FDR=qval)
  return(res)
}

get.REM2<-function(em,vm, n, REM.type){
  k<-ncol(em)
  tau2.res <-get.tau2(em = em ,vm = vm ,k=k, n=n, REM.type=REM.type)
  tau2 <- tau2.res$tau2
  wt <- tau2.res$weight
  Q.val<-get.Q(em = em , wt = wt)
  temp.wt<-1/(vm+tau2)  
  mu.hat<-rowSums(temp.wt*em)/rowSums(temp.wt)
  mu.var<-1/rowSums(temp.wt)
  Qpval <- pchisq(Q.val, df = k - 1, lower.tail = FALSE)
  z.score<-mu.hat/sqrt(mu.var)
  z.p<-2*(1-pnorm(abs(z.score)))
  qval<-p.adjust(z.p,method="BH")
  res<-list(mu.hat=mu.hat,mu.var=mu.var,Qval=Q.val,Qpval=Qpval,tau2=tau2,
            zval=z.score,pval=z.p,FDR=qval)
  return(res)
}

get.REM<-function(em, vm, n, REM.type, pe=NULL,pv=NULL){
  k<-ncol(em)
  tau2.res <-get.tau2(em = em ,vm = vm ,k=k, n=n, REM.type=REM.type)
  tau2 <- tau2.res$tau2
  wt <- tau2.res$weight
  Q.val<-get.Q(em = em ,wt = wt)
  temp.wt<-1/(vm+tau2)	
  mu.hat<-rowSums(temp.wt*em)/rowSums(temp.wt)
  mu.var<-1/rowSums(temp.wt)
  Qpval <- pchisq(Q.val, df = k - 1, lower.tail = FALSE)
  z.score<-get.REM2(em = em ,vm = vm, n=n, REM.type="HO")$zval
  if(!is.null(pe)&!is.null(pv)){
    rnum<-which(apply(em,1,function(x) !any(is.na(x))))
    Z0<-matrix(get.REM2(pe,pv, n=n, REM.type=REM.type)$zval,nrow(em),
               nrow(pe)/nrow(em))
    z.p<-rep(NA,nrow(em))
    z.p[rnum]<-perm.p(z.score[rnum],Z0[rnum,],"abs")
  }else{
    z.p<-2*(1-pnorm(abs(z.score)))
  }
  qval<-p.adjust(z.p,method="BH")
  res<-list(mu.hat=mu.hat,mu.var=mu.var,Qval=Q.val,Qpval=Qpval,tau2=tau2,
            zval=z.score,pval=z.p,FDR=qval)
  return(res)
}


get.RP<-function(dat,lbl, nperm = 100, logged = TRUE) {
  gene.names<-row.names(dat[[1]]) 
  num.perm<-nperm
  num.ori<-length(lbl) 
  num.gene<-nrow(dat[[1]]) 
  nk<-unlist(lapply(lbl,function(x) length(x))) 
  origin<-rep(1:num.ori,nk)
  
  data<-cl<-NULL
  for (k in 1:num.ori)
  {
    data<-cbind(data,dat[[k]])
    cl<-c(cl,lbl[[k]])
  }
  
  total.sam = length(origin)
  total.sam1<-ncol(data)
  if (total.sam != total.sam1) 
    stop("the lbl number does not match the dat columns")
  
  data.pre = OriginxyCall(data=data, cl=cl, origin=origin)
  y = data.pre$data2
  
  data = as.matrix(data)
  mode(data) = "numeric"
  NA.genes <- NULL
  if (any(is.na(data))) {
    NA.genes <- unique(ceiling(which(is.na(t(data)))/ncol(data)))
    cat("Warning: There are", length(NA.genes), "genes with at least one 
            missing value.", 
        "\n", "\n")
    if (na.rm) 
      data[NA.genes, ] <- NaReplace2(data[NA.genes, ],origin)
    if (!na.rm) 
      cat(" This value is not used to compute rank product.", 
          "\n", "\n")
  }
  
  if (!is.null(y)) {
    num.class = 2
    data1.all = data.pre$data1
    data2.all = data.pre$data2
    fold.change = NULL
    for (l in 1:num.ori) {
      data1 = as.matrix(data1.all[[l]])
      data2 = as.matrix(data2.all[[l]])
      data1.ave = apply(data1, 1, mean)
      data2.ave = apply(data2, 1, mean)
      if (logged) {
        fold.change = cbind(fold.change,(data1.ave-data2.ave))
      }
      else {
        fold.change = cbind(fold.change, (data1.ave/data2.ave))
      }
      rm(data1, data2, data1.ave, data2.ave)
    }
    ave.fold.change = apply(fold.change, 1, mean)
  }
  if (is.null(y)) {
    num.class = 1
    data1.all = data.pre$data1
    data2.all = data.pre$data2
    fold.change = NULL
    for (l in 1:num.ori) {
      data1 = as.matrix(data1.all[[l]])
      fold.change = cbind(fold.change, apply(data1, 1, 
                                             mean))
      rm(data1)
    }
    ave.fold.change = apply(fold.change, 1, mean)
  }
  RP.ori.out.upin2 = RankProd2(data1.all, data2.all, num.ori, 
                               num.gene, logged, num.class, rev.sorting = FALSE)
  RP.ori.upin2 = RP.ori.out.upin2$RP       
  rank.ori.upin2 = rank(RP.ori.upin2)
  RP.ori.out.downin2 = RankProd2(data1.all, data2.all, num.ori, 
                                 num.gene, logged, num.class, rev.sorting = TRUE)
  RP.ori.downin2 = RP.ori.out.downin2$RP   
  rank.ori.downin2 = rank(RP.ori.downin2)
  RP.perm.upin2 <- matrix(NA, num.gene, num.perm)
  RP.perm.downin2 <- matrix(NA, num.gene, num.perm)
  cat("Starting ", num.perm, "permutations...", "\n")
  for (p in 1:num.perm) {
    new.data.temp = NewdataCom(data1.all, data2.all, num.ori,num.class)
    new.data1.all = new.data.temp$new.data1.all
    new.data2.all = new.data.temp$new.data2.all
    temp1 = RankProd2(new.data1.all, new.data2.all, num.ori,num.gene, 
                      logged, num.class, rev.sorting = FALSE)
    RP.perm.upin2[, p] = temp1$RP
    rm(temp1)
    temp2 = RankProd2(new.data1.all, new.data2.all, num.ori,num.gene, 
                      logged, num.class, rev.sorting = TRUE)
    RP.perm.downin2[, p] = temp2$RP
    rm(temp2)
  }
  pval.upin2<-perm.p(RP.ori.upin2,RP.perm.upin2,tail="low") 
  pval.downin2<-perm.p(RP.ori.downin2,RP.perm.downin2,tail="low") 
  qval.upin2<-p.adjust(pval.upin2,method="BH")
  qval.downin2<-p.adjust(pval.downin2,method="BH")
  names(RP.ori.upin2)<-names(pval.upin2)<-names(qval.upin2)<- 
    names(RP.ori.downin2)<-names(pval.downin2)<-names(qval.downin2)<-gene.names
  
  res<-list(meta.stat.up=RP.ori.upin2,pval.up=pval.upin2,FDR.up=qval.upin2,
            meta.stat.down=RP.ori.downin2,pval.down=pval.downin2,
            FDR.down=qval.downin2, AveFC = ave.fold.change)    
  return(res)
}

OriginxyCall<-function (data, cl, origin, sum = FALSE) {
  lev <- unique(cl)
  uni.cl <- length(lev)
  if (uni.cl > 2) 
    stop("There is something wrong with the classlabels")
  ori.lev <- unique(origin)
  uni.ori <- length(ori.lev)
  cat(" The data is from ", uni.ori, "different origins \n \n")
  if (min(ori.lev) != 1 | max(ori.lev) != uni.ori) {
    cat("Warning: origins labels are not numbered from 1 to ", 
        uni.ori, "\n", "\n")
  }
  if (uni.cl == 1) {
    if (sum) {
      cat("Rank Sum analysis for one-class case", "\n", 
          "\n")
    }
    else {
      cat("Rank Product analysis for one-class case", "\n", 
          "\n")
    }
    if (lev != 1) {
      cat("warning: Expected classlabel is 1, cl will hus be set to 1.", 
          "\n", "\n")
      cl = rep(1, length(cl))
    }
    data2 <- NULL
    data1 <- vector("list", uni.ori)
    for (c in 1:uni.ori) {
      data1[[c]] = data[, origin == ori.lev[[c]]]
    }
  }
  if (uni.cl == 2) {
    if (sum) {
      cat("Rank Sum analysis for two-class case", "\n", 
          "\n")
    }
    else {
      cat("Rank Product analysis for two-class case", "\n", 
          "\n")
    }
    if (min(lev) != 0 | max(lev) != 1) {
      cat("Warning: Expected classlabels are 0 and 1. cl will 
                thus be set to 0 and 1.", 
          "\n", "\n")
      cl[which(cl == min(lev))] <- 0
      cl[which(cl == max(lev))] <- 1
    }
    data2 <- vector("list", uni.ori)
    data1 <- vector("list", uni.ori)
    for (c in 1:uni.ori) {
      index1 <- which(origin == ori.lev[[c]] & cl == 0)
      index2 <- which(origin == ori.lev[[c]] & cl == 1)
      if (length(index1) == 0 | length(index1) == 0) 
        stop("Error: data from different origins should contain 
                     data from both classs")
      data1[[c]] <- data[, index1]
      data2[[c]] <- data[, index2]
      rm(index1, index2)
    }
  }
  list(data1 = data1, data2 = data2)
}

RankProd2<-function (data1.all, data2.all, num.ori, num.gene, logged, 
                     num.class, rev.sorting) {
  num.rank.all = 0
  rank.rep.all = t(t(1:num.gene))
  for (l in 1:num.ori) {
    data1 = data1.all[[l]]
    data2 = data2.all[[l]]
    data1 = as.matrix(data1)
    if (num.class == 2) {
      data2 = as.matrix(data2)
    }
    temp = RankComp(data1, data2, logged, num.class, rev.sorting)
    rank.rep.all = cbind(rank.rep.all, temp$rank)
    num.rank.all = num.rank.all + temp$num.rank
    rm(temp)
  }
  rank.all = rank.rep.all[, -1]
  rank.prod.temp = rank.all^(1/num.rank.all)
  rank.prod = apply(rank.prod.temp, 1, prod)
  rank.prod[num.rank.all == 0] = NA
  list(RP = rank.prod, rank.all = rank.all)
}

RankComp<-function (data1, data2, logged, num.class, rev.sorting) 
{
  num.gene = dim(data1)[1]
  if (num.class == 2) {
    if (rev.sorting) {
      data1.wk <- data2
      data2.wk <- data1
    }
    else {
      data1.wk <- data1
      data2.wk <- data2
    }
    k1 = dim(data1.wk)[2]
    k2 = dim(data2.wk)[2]
    num.rep = k1 * k2
    data.rep = matrix(NA, num.gene, num.rep)
    if (logged) {
      for (k in 1:k1) {
        temp = ((k - 1) * k2 + 1):(k * k2)
        data.rep[, temp] = data1.wk[, k] - data2.wk
      }
    }
    else {
      for (k in 1:k1) {
        temp = ((k - 1) * k2 + 1):(k * k2)
        data.rep[, temp] = data1.wk[, k]/data2.wk
      }
    }
    rank.rep = apply(data.rep, 2, rank)
  }
  if (num.class == 1) {
    data.rep = data1
    if (rev.sorting) {
      num.rep = dim(data1)[2]
      rank.temp = matrix(NA, num.gene, num.rep)
      for (r in 1:num.rep) {
        rank.temp[, r] = rank(data1[, r], na.last = FALSE)
      }
      rank.rep = (num.gene + 1) - rank.temp
    }
    else {
      rank.rep = apply(data1, 2, rank)     ### where rank is obtained
    }
  }
  rank.rep[is.na(data.rep)] = 1
  num.rank = apply(is.na(data.rep) == FALSE, 1, sum)
  list(rank = rank.rep, num.rank = num.rank)
}

NewdataCom <-function(data1.all,data2.all,num.ori,num.class)
{
  new.data1.all <- vector("list",num.ori)
  new.data2.all <- vector("list",num.ori)  
  for ( l in 1:num.ori ){
    
    data1 <- as.matrix(data1.all[[l]])
    if (num.class == 2) {data2 <- as.matrix(data2.all[[l]])}
    
    temp.data <- Newdata(data1,data2,num.class)
    new.data1.all[[l]] <- temp.data$new.data1
    new.data2.all[[l]] <- temp.data$new.data2
  }
  if(num.class == 1) {new.data2.all <- NULL} 
  list(new.data1.all = new.data1.all,new.data2.all = new.data2.all)
}

Newdata<-function(data1,data2,num.class)
{
  
  k1 <- dim(data1)[2]
  num.gene <- dim(data1)[1]
  new.data1 <- matrix(NA,num.gene,k1)
  
  for (k_1 in 1:k1)
  {
    temp <- sample(1:num.gene,num.gene,replace = FALSE, prob = NULL)
    new.data1[,k_1] <- data1[temp,k_1];
    rm(temp)
  }
  rm(k_1,k1)
  
  if (num.class == 2) {
    k2 <- dim(data2)[2]
    new.data2 <- matrix(NA,num.gene,k2) 
    for (k_2 in 1:k2) {   
      temp <- sample(1:num.gene,num.gene,replace = FALSE, prob = NULL)
      new.data2[,k_2] <- data2[temp,k_2];
    }
    rm(k_2,k2)
  }
  
  if (num.class == 1) { new.data2=NULL}
  list(new.data1 = new.data1,new.data2 = new.data2)
  
}

######################
####  Data for AW.new
######################

awFisherData = list(logPTarget=-log(c(0.99, 0.95,0.9,0.8,0.7,0.5,0.3,
                                      1e-1,1e-2,1e-3,1e-5,1e-7,1e-10, 
                                      1e-15, 1e-25, 1e-80,1e-200)),
                    original=matrix(c(0.0991299,0.2307506,0.3616958,0.4908891,
                                      0.6012938,0.7007693,0.8010211,0.8852742,0.9674844,
                                      1.113856,1.304995,1.612286,2.274926,3.816328,
                                      6.379592,10.03372,15.8871,28.39492,50.09994,107.3815,
                                      0.2532749,0.4589056,0.6426983,0.7975563,0.9350853,
                                      1.056027,1.165338,1.267054,1.36635,1.545324,
                                      1.810291,2.243907,3.162265,5.079755,8.110184,
                                      12.3064,18.90731,32.46575,55.64708,115.5134,
                                      0.3785702,0.6240129,0.8244336,0.9922922,1.142877,
                                      1.27297,1.397758,1.508944,1.621461,1.83459,
                                      2.140719,2.646127,3.666674,5.782382,9.040087,
                                      13.49952,20.3788,34.47629,58.35517,119.4172,
                                      0.5911318,0.8747565,1.10027,1.286399,1.455918,
                                      1.606457,1.752032,1.885532,2.017578,2.267131,
                                      2.638745,3.232473,4.391949,6.732669,10.29312,
                                      15.05824,22.34432,37.068,61.78113,124.354,
                                      0.7887108,1.100477,1.346518,1.550477,1.736986,
                                      1.905264,2.067322,2.21877,2.36439,2.644688,
                                      3.053511,3.714803,4.978453,7.494532,11.26784,
                                      16.27226,23.84434,38.99281,64.31994,128.0403,
                                      1.21693,1.582346,1.865415,2.105853,2.32561,
                                      2.519655,2.706602,2.884535,3.060112,3.404064,
                                      3.876032,4.631569,6.076903,8.897185,13.01535,
                                      18.40227,26.42127,42.35352,68.65255,134.2008,
                                      1.812172,2.224871,2.54456,2.815906,3.066485,
                                      3.299503,3.521331,3.727825,3.931819,4.326293,
                                      4.869171,5.727013,7.353057,10.47594,14.96737,
                                      20.72703,29.19919,45.86845,73.2044,140.5048,
                                      3.006446,3.498343,3.882972,4.202215,4.504013,
                                      4.763724,5.026574,5.270521,5.514426,5.981952,
                                      6.619728,7.632539,9.526118,13.06491,18.04471,
                                      24.38708,33.53347,51.32331,80.00293,149.9905,
                                      5.412572,5.978331,6.421937,6.818137,7.178373,
                                      7.489069,7.806217,8.101734,8.40125,8.966433,
                                      9.724746,10.9367,13.20626,17.32554,22.98528,
                                      30.12292,40.1709,59.49315,90.07466,163.6995,
                                      7.771121,8.360937,8.862098,9.283587,9.705136,
                                      10.04135,10.41136,10.72468,11.07541,11.6783,
                                      12.54049,13.87575,16.40208,20.96055,27.11554,
                                      34.78738,45.51271,65.95685,97.89573,174.1557,
                                      12.44521,13.09599,13.6526,14.10011,14.54778,
                                      14.986,15.30801,15.77355,16.13772,16.79953,17.82213,
                                      19.38547,22.27196,27.42194,34.3746,42.84719,
                                      54.5786,76.72491,110.749,191.039,17.09525,17.75282,
                                      18.3716,18.88985,19.37994,19.78108,20.14133,20.75287,
                                      21.01495,21.71198,22.83626,24.6096,27.71989,33.2676,
                                      40.9529,49.94472,62.55643,85.9504,121.8311,205.3097,
                                      24.03921,24.68517,25.40527,25.95085,26.45911,
                                      26.87557,27.39865,27.67894,28.19278,29.0289,
                                      30.23832,32.32301,35.47315,41.61328,49.97596,
                                      59.59598,73.60026,98.42317,136.9592,224.2329,
                                      35.56544,36.24341,36.8793,37.71307,38.04985,
                                      39.3861,38.91839,39.42456,39.97605,40.79525,42.19192,
                                      44.34679,47.92825,54.43455,64.27763,74.33595,
                                      90.52811,116.6218,161.774,248.6012,
                                      58.59329,59.29901,59.84986,60.87042,61.81504,
                                      61.85231,61.88179,62.51328,63.15607,64.2352,65.50203,
                                      68.01901,72.74197,79.66,89.91935,101.7426,118.9691,
                                      150.7868,192.5185,292.6912,185.1815,186.014,
                                      186.6892,187.0414,187.5587,187.9156,189.0134,
                                      189.5565,190.0352,192.4073,192.6559,194.5913,
                                      201.7345,209.2415,223.7879,239.5384,258.8478,
                                      297.156,352.8421,474.5017,461.2872,462.3472,462.4968,
                                      462.8381,463.8254,463.8606,464.7668,464.4796,
                                      465.2994,465.9452,467.9054,471.0704,474.8733,
                                      484.9073,498.1799,512.5592,541.8607,585.8783,
                                      651.6489,798.9819), ncol=17),
                    cond=matrix(c(0.03810959,0.0621685,0.08131234,0.09920804,0.111354,
                                  0.122629,0.1349937,0.1432366,0.1520866,
                                  0.1679167,0.184224,0.2087301,0.2434008,0.2871737,0.331107,
                                  0.3677778,0.4102335,0.4527058,0.5024291,0.5630171,
                                  0.127912,0.1818843,0.2207695,0.2516294,0.2775123,0.2983333,
                                  0.3172555,0.337089,0.3547296,0.3786171,
                                  0.4092174,0.4455226,0.4985854,0.5779685,0.6441868,
                                  0.69751,0.7509269,0.827722,0.8875773,0.9727183,
                                  0.2186853,0.2913149,0.3406531,0.3804417,0.4133604,
                                  0.4433051,0.4685546,0.4876554,0.5078706,0.541914,0.579527,
                                  0.6256311,0.6896221,0.7820427,0.8602653,0.9251285,0.9848993,
                                  1.076998,1.144793,1.237292,0.3893207,0.4826316,0.5499885,
                                  0.6051856,0.6505681,0.6844553,0.7160721,0.7424228,0.7674958,
                                  0.8070018,0.8569188,0.9147883,0.9964052,1.10293,1.196161,
                                  1.270607,1.347474,1.446223,1.527901,1.636886,
                                  0.5585807,0.6742425,0.7560492,0.8207702,0.8732926,
                                  0.9129082,0.9493195,0.9797121,1.008026,1.05531,
                                  1.111715,1.180343,1.269962,1.391778,1.492729,1.57213,
                                  1.662787,1.768504,1.850884,1.979766,0.95566,1.110985,
                                  1.216068,1.29076,1.355766,1.405023,1.445523,1.484516,
                                  1.518811,1.580653,1.645448,1.726037,1.839056,1.978376,
                                  2.099252,2.193948,2.290192,2.411119,2.51721,2.651896,
                                  1.52536,1.713372,1.832813,1.918937,1.993254,2.056027,
                                  2.107726,2.152054,2.201163,2.271805,2.348409,
                                  2.440764,2.568491,2.726969,2.870715,2.971539,3.078847,
                                  3.214829,3.323123,3.47035,2.70562,2.929678,3.080541,
                                  3.193354,3.277625,3.338019,3.4043,3.450513,3.511548,
                                  3.599451,3.674784,3.804746,3.947966,4.123648,4.289046,
                                  4.406467,4.52903,4.700287,4.825234,4.960964,5.077873,
                                  5.34332,5.526001,5.657974,5.772755,5.848972,5.931133,
                                  6.003743,6.089314,6.176835,6.254821,6.350939,6.530929,
                                  6.826697,6.97927,7.100861,7.447727,7.758332,7.515395,
                                  7.543793,7.415679,7.711341,7.889183,8.047575,8.198342,
                                  8.252384,8.337884,8.467926,8.543158,8.638408,
                                  8.671341,8.793482,9.047245,9.584596,9.452446,9.45232,
                                  9.892516,10.00523,9.830353,10.04764,12.0838,12.39359,
                                  12.6006,12.75797,12.90931,13.02841,13.0343,13.156,
                                  13.22006,13.26682,13.39121,13.35111,13.73286,14.03898,
                                  14.02277,14.26102,14.27311,14.44403,14.33591,14.14482,
                                  16.70767,17.00769,17.31219,17.4414,17.57655,17.72012,
                                  17.67593,17.99483,17.70582,17.76762,18.03914,18.0793,
                                  18.0617,18.34354,19.13029,20.76382,18.30004,18.7526,
                                  18.83723,18.87949,23.64351,23.91838,24.24045,24.47294,
                                  24.67672,24.56047,24.69437,24.52458,24.53108,24.68174,
                                  24.84996,25.85694,25.28258,25.2019,25.52906,25.05026,
                                  25.52343,26.48087,25.91837,25.62153,35.17247,35.44539,
                                  35.75883,36.02557,35.93328,36.91222,35.91838,35.94448,
                                  36.09843,35.89931,37.04063,36.88587,36.29341,36.53636,
                                  42.25099,36.58584,36.77084,36.50841,41.56533,37.13722,
                                  58.20533,58.44961,58.59536,59.20801,59.48873,59.02293,
                                  58.63219,58.86448,58.81296,58.91426,59.27569,59.34844,
                                  60.65453,59.62144,59.50907,59.61764,59.74752,59.77,
                                  59.80256,66.04942,184.7801,185.2772,185.0311,185.1535,
                                  185.1581,185.0703,185.2399,185.5662,185.5465,186.0258,
                                  185.8041,185.8345,185.7155,186.1403,186.0538,186.0937,
                                  186.0333,186.2792,186.202,186.6681,460.8942,461.4731,
                                  461.2796,461.2338,461.6603,461.3963,461.3641,461.362,
                                  461.3927,462.4369,462.3063,461.5865,461.6614,461.646,
                                  461.8454,462.9842,462.1843,462.2884,462.3873,463.0256), 
                                ncol=17),
                    uncond=matrix(c(0.01167915,0.01380246,0.01517486,0.01701253,0.01780689,
                                    0.01869392,0.02015334,0.02096671,0.02257714,0.02432555,
                                    0.02654229,0.02890041,0.03275772,0.04085687,0.04920061,
                                    0.0568219,0.06522525,0.09740171,0.1110982,0.1383659,
                                    0.06097223,0.06811177,0.07404343,0.07932669,0.08442349,
                                    0.0883944,0.09298591,0.0964361,0.09975088,0.1070759,
                                    0.1121707,0.1199021,0.1349449,0.1580033,0.1787628,
                                    0.1968977,0.2150478,0.2727937,0.2969215,0.342652,
                                    0.1227619,0.1360761,0.1473208,0.1565375,0.1659376,
                                    0.1737538,0.1820478,0.1874003,0.1911773,0.1991744,
                                    0.2086036,0.2226894,0.245687,0.2794811,0.3077283,
                                    0.3285149,0.3587663,0.4260081,0.4548423,0.5025304,
                                    0.2569734,0.2807185,0.3004544,0.3166578,0.3310639,
                                    0.3423149,0.3539376,0.3638423,0.3707911,0.3849804,
                                    0.4027341,0.4267229,0.4548959,0.4978543,0.5393592,
                                    0.5693082,0.6173452,0.6846767,0.7185779,0.7757726,
                                    0.4063705,0.4418796,0.4701133,0.4912607,0.510468,
                                    0.528245,0.5394942,0.5477942,0.5610035,0.5793259,
                                    0.6046962,0.6299256,0.6682264,0.7182808,0.7674903,
                                    0.8059049,0.8593518,0.9323056,0.9668378,1.033236,
                                    0.7807272,0.843123,0.8808574,0.9126672,0.9396796,
                                    0.960439,0.978473,0.9922291,1.008917,1.038971,1.071566,
                                    1.107042,1.151322,1.220174,1.283471,1.330622,1.381731,
                                    1.465399,1.510273,1.570922,1.342772,1.422975,1.47796,
                                    1.514661,1.556413,1.586022,1.608228,1.630456,1.652725,
                                    1.693232,1.725,1.76851,1.824637,1.904191,1.980362,
                                    2.024488,2.089696,2.179498,2.226941,2.267269,
                                    2.529015,2.652573,2.730094,2.785313,2.832947,2.862264,
                                    2.886662,2.909889,2.945636,2.987453,3.027188,
                                    3.086178,3.145357,3.244686,3.322164,3.375402,3.441574,
                                    3.537447,3.594048,3.606768,4.940525,5.100422,5.211391,
                                    5.288603,5.345469,5.379895,5.410683,5.44296,5.506372,
                                    5.538015,5.614381,5.612352,5.736223,5.811657,5.860137,
                                    5.914662,6.229093,6.139462,6.136908,6.083464,
                                    7.308872,7.500962,7.616301,7.699705,7.776036,7.817443,
                                    7.831599,7.91947,7.970223,8.03098,8.046011,
                                    8.037557,8.229807,8.207456,8.27931,8.353512,8.695223,
                                    8.395307,8.444783,8.332947,11.99451,12.21962,12.37417,
                                    12.44409,12.53803,12.5802,12.60703,12.69658,12.79327,
                                    12.80341,12.72086,12.68056,12.90105,12.95976,12.90918,
                                    13.10627,12.74418,12.65242,12.97392,12.54703,
                                    16.64509,16.87218,17.06891,17.19426,17.28926,17.30334,
                                    17.25235,17.70002,17.41296,17.31092,17.42144,
                                    17.36365,17.58575,17.48862,18.22938,17.31067,17.21504,
                                    17.31552,16.94165,17.2159,23.59017,23.81475,24.08688,
                                    24.22897,24.25407,24.15445,24.6323,24.20535,24.21415,
                                    24.25998,24.17379,24.42562,24.20963,24.07561,24.17858,
                                    24.3309,24.21982,23.96657,23.18492,23.19917,35.12209,
                                    35.38765,35.54676,35.93301,35.72465,37.15706,35.62476,
                                    35.65893,35.79171,35.74076,35.81368,35.84007,35.68559,
                                    35.55025,35.65699,35.59003,36.08693,35.01463,34.55093,
                                    34.54058,58.16783,58.45219,58.56312,59.10713,59.24837,
                                    58.97294,58.55741,58.68723,58.78797,58.88322,58.60795,
                                    59.02749,58.99561,58.86774,58.73229,58.48511,57.83821,
                                    57.56131,57.55227,57.56808,184.7882,185.2093,185.3756,
                                    185.2619,185.2936,185.2547,185.5489,185.53,185.5679,
                                    186.3204,185.5461,185.3334,186.3847,184.6643,184.2206,
                                    184.1881,184.2057,184.22,184.2071,184.2049,
                                    460.981,461.4302,461.0973,460.9452,461.0463,460.7887,
                                    460.7704,460.6247,460.6181,460.5617,460.5465,460.5333,
                                    460.538,460.5324,460.5328,460.5398,460.51,460.52,460.531,
                                    460.5183),ncol=17),
                    nList = c(2:10,12,15,20,30,50,80,120,180,300,500,1000),
                    totalN = 100000)



#################### Utility ####################
## Utility functions for plotting and other non stats purposes (Internal)
## Author: Tianzhou Ma
## Institution: University of pittsburgh
## Date: 08/29/2016

#----------------------------------------------# 
# Matrix manipulation methods 
#----------------------------------------------# 
# Flip matrix (upside-down) 
flip.matrix <- function(x) { 
  mirror.matrix(rotate180.matrix(x)) 
}   

# Mirror matrix (left-right) 
mirror.matrix <- function(x) { 
  xx <- as.data.frame(x); 
  xx <- rev(xx); 
  xx <- as.matrix(xx); 
  xx; 
} 

# Rotate matrix 180 clockworks 
rotate180.matrix <- function(x) { 
  xx <- rev(x); 
  dim(xx) <- dim(x); 
  xx; 
} 

#-----------------------------------------------------#
# generate nPr    with repetition                     #
# for generating all possible weights                 #
#-----------------------------------------------------#  
permut<-function (n, r) { 
  v<-1:n
  sub <- function(n, r, v) {
    if (r == 1) matrix(v, n, 1)
    else if (n == 1) matrix(v, 1, r)
    else {
      inner <- Recall(n, r - 1, v)
      cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
                                                ncol = ncol(inner), nrow = nrow(inner) * n, byrow = TRUE))
    }
  }
  sub(n, r, v[1:n])
}
#-----------------------------------------------------#
#        wt=possible weights for each datasets
#        n=# of datasets                              #
#-----------------------------------------------------#
gen.weights<-function(wt,n) {
  comb<-permut(n=length(wt),r=n)
  weight<-matrix(wt[comb],ncol=n)
  return(weight[-1,])
}

#--------------------------------------------------------------#
# plot a matrix of gene expression data with rows are genes
# columns are samples
# n: # of samples in each study
# ni: # of samples in each class of each study
#--------------------------------------------------------------#
plot.matrix<-function(mat,color="GR") {
  n<-attr(mat,"n")
  ni<-attr(mat,"ni")
  label<-attr(mat,"label")
  dataset.name <- attr(mat,"dataset.name")
  group.name <- attr(mat,"group.name")
  ref.group <- attr(mat,"ref.group")
  nc<-ncol(mat)
  nr<-nrow(mat)
  cexCol1<-1/log(nc)
  cexCol2<-1/log10(nc)
  cexRow<-1/log(nr,4)
  K<-length(n)
  #mycol<- c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", 
  #"#238443", "#006837", "#004529")
  colfun<-colorRampPalette(c("green","black","red"))
  if(color=="BY") colfun<-colorRampPalette(c("blue","black","yellow"))
  mycol<-colfun(16)
  mval<-min(max(abs(mat),na.rm=T),3)
  xcut<-seq(-mval,mval,len=length(mycol)-1)
  xcut<-c(-100,xcut,100)
  m <- matrix(1:2, 2, 1)
  nf<-layout(m, heights=c(6, 1))
  par(mar=c(2,3,2,5))
  image(x=1:nc,y=1:nr,z=t(flip.matrix(mat)),col=mycol,axes=FALSE,xlab = "", 
        ylab = "", breaks=xcut)
  #axis(3,1:nc,labels=label,las=2,line=-0.5,tick=0,cex.axis=cexCol1)
  axis(3,cumsum(ni)-ni/2+0.5,labels=rep(group.name,K),las=1,line = -0.5,
       tick=0,cex.axis=cexCol2*2,font=2)
  axis(4,nr:1,labels=(row.names(mat)), las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow*2,font=2)
  axis(1,cumsum(n)-n/2+0.5,labels=dataset.name,las=1,line = -1,
       tick=0,cex.axis=cexCol2*2,font=2)
  #---distinguish studies----#
  abline(v=cumsum(n)+0.5,lwd=2,col="white")
  #---distinguish classes----#
  if (is.null(ncol(label))) abline(v=cumsum(ni)+0.5,lwd=2,col="white",lty=2)
  #---------if AW method we add category information on the plot---------#
  if (!is.null(attr(mat,'category'))) {
    cat<-attr(mat,'category')
    at.line <-cumsum(table(cat))+0.5
    sum.pos <- cumsum(table(cat))
    at <- rep(NA, length(sum.pos))
    for (i in 2:length(sum.pos)){
      at[i] <- sum.pos[i-1] + (sum.pos[i] - sum.pos[i-1])/2 + 0.5
    }
    at[1] <- sum.pos[1] + (sum.pos[1] - 0)/2 + 0.5
    axis(2,at-0.5,labels=rev(unique(cat)),tick = 0, las=1,
         cex.axis = (cexRow+0.2)*2,font=2)
    abline(h=at.line,lwd=2,col="white")
  }
  #----add legend---------------#
  l<-length(xcut)
  image(1:(l-1),0.5,matrix(xcut[-1],nrow=l-1,ncol=1),col=mycol,breaks=xcut,
        axes=F,xlab="",ylab="")
  marcas<-(0:(l-1))+0.5
  axis(1,marcas,round(xcut,1),tick=0.5,cex.axis=cexCol2*2,line=-0.5,
       font=2)
}

#-----------------------------------------------------------#
# order genes for plotting obtained from all methods except AW  #
#-----------------------------------------------------------#
order.genes.simple<-function(dat) {
  r<-hclust(dist(dat))$order
  plot.genes<-dat[r,]
  return(plot.genes)
}

#---------------------------------------------------------------------------#
# output genes according to the order of categories defined from AW.weight  #
#---------------------------------------------------------------------------#
order.genes.AW<-function(dat,AW.weight) {
  K<-ncol(AW.weight)
  #wt.group<-gen.weights(c(0,1),K) # generate all possible weights
  wt.group<-do.call(expand.grid, rep(list(c(0, 1)), K))[-1,]
  wt.group<-wt.group[order(apply(wt.group,1,sum)),]
  row.names(wt.group)<-nrow(wt.group):1 
  
  ng<-nrow(dat) # number of significant genes
  group<-rep(NA,ng)
  
  #---------order the OW categories nicely ------------------#
  for (i in 1:ng) {
    for (j in 1:nrow(wt.group)) {
      if (sum(AW.weight[i,]==wt.group[j,])==K) {
        group[i]<-row.names(wt.group)[j]
        next
      }
    }
  }
  #perform hierarchical clustering in each weight group
  plot.genes.ordered<-NULL
  for (i in sort(as.numeric(names(table(group))))) {
    x<-subset(dat,group==i)
    if (nrow(x)>2) {
      newx<-x[hclust(dist(x))$order,]
    } else newx<-x
    plot.genes.ordered<-rbind(plot.genes.ordered,newx)
  }
  attr(plot.genes.ordered,"category")<-
    apply(wt.group,1,paste,collapse=',')[sort(group)]
  return(plot.genes.ordered)
}

#-----------------------------------------------------------------------------#
# check gene names
#-----------------------------------------------------------------------------#
check.exp<-function(x)
{
  if (is.null(row.names(x[[1]]))) {
    K<-length(x)
    ng<-nrow(x[[1]])
    for (k in 1:K) {
      row.names(x[[k]])<-paste("gene",1:ng)
    }
  }
  return(x)
}

#------------------------------------------------------------------------------#
# check dimensions and size of argument
#------------------------------------------------------------------------------#
check.dim<-function(x,y,ind.method,meta.method,paired){
  K<-length(x)
  nperstudy<-sapply(x,function(y) ncol(y))
  nlabels<-sapply(y,function(z) length(z))						
  if(sum(nperstudy==nlabels)!=K) {
    stop(cat("The number of samples does not match with the dimension of 
             labels in study(s)",paste((1:K)[nperstudy!=nlabels],"",
                                       collapse=","),"!"))
  }
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1) {
    if(length(ind.method)!=K)stop(paste('Argument "ind.method" should be 
                                        a character vector of size',K))
  }
  if(("REM"%in%meta.method|"FEM"%in%meta.method)&(length(paired)!=K)) {
    stop(paste('Argument "paired" should be a logical vecter of size',K))
  }
}

#-----------------------------------------------------------------------------#
#   check if parametric is ok                                                 #
#-----------------------------------------------------------------------------#
check.parametric<-function(meta.method,parametric)
{
  if (parametric==TRUE&sum(meta.method%in%c("SR","PR","rankProd",
                                            "Fisher.OC","minMCC"))>0) {
    stop(paste("There is no parametric result for",meta.method))
  }
}

#-----------------------------------------------------------------------------#
#   check tail                                         #
#-----------------------------------------------------------------------------#
check.tail<-function(meta.method,tail)
{
  if (tail=='abs'&length(grep('.OC',meta.method))>0) {
    stop(paste("If you chose",meta.method,",then you should specify the 'tail' 
               to be either 'high' or 'low'"))
  }
}

#-----------------------------------------------------------------------------#
#   check individual methods and response                                       
#-----------------------------------------------------------------------------#
check.indmethod<-function(y, resp.type, data.type,ind.method, 
                          select.group,tail,paired) {  
  K<-length(y)
  for(k in 1:K) {
    #  if(data.type[k] == "continuous") {
    #   if(!(ind.method[k] %in%c('limma','sam',"pearsonr","spearmanr",'logrank'))) {
    #    stop ("Incorrect method for microarray or RNAseq FPKM")
    #   }
    #  } else if (data.type[k] =='discrete') {
    #     if(!(ind.method[k] %in%c('edgeR','DESeq2',"spearmanr",'limmaVoom') ) ) {
    #       stop ("Incorrect method for RNAseq count")
    #     } 
    #    }   
    if(resp.type %in% c("twoclass", "multiclass")) {
      if(!(ind.method[k] %in%c("limma","sam","limmaVoom","edgeR","DESeq2") ) ) {
        stop ("Incorrect method for the response")
      }  
    }
    if(resp.type =="twoclass") {
      if(nlevels(as.factor(y[[k]]))!=2 && length(select.group)!=2 ) {
        stop (cat(resp.type," requires two levels ") )
      }
      if (is.null(paired)) {
        stop(paste("you need to specify the logical value for 'paired' 
               for study", k))
      }
      
    }     
    if(resp.type=="multiclass") {
      if(nlevels(as.factor(y[[k]])) <=2) {
        stop (cat(resp.type," requires at least three levels") )
      }
      if(tail!="abs") {
        stop (cat(resp.type," cannot perform one-sided test") )
      }
    }  
    if(resp.type %in% c("continuous") ) {
      if(!(ind.method[k] %in%c("pearsonr","spearmanr") ) ) { 
        stop ("Incorrect method for the response") 
      }     
      if(!is.numeric(y[[k]])) {
        stop (cat(resp.type," response has to be numeric") )
      }
    }
    else if (resp.type %in% c("survival") ) {
      if(!(ind.method[k] %in%c("logrank") ) ) {
        stop ("Incorrect method for the response")
      }
      if(is.null(y[[k]][,1])) {
        stop( cat("dataset",k, "is missing the event time")) 
      }
      if(is.null(y[[k]][,2])) {
        stop( cat("dataset",k, "is missing the censoring status"))
      } 
      
      if(!is.numeric(y[[k]][,1])) {
        stop( cat("dataset", k, ": Survival time has to be numeric")) 
      }
      
      if(!(y[[k]][,2] %in% c(0,1)| y[[k]][,2] %in% c('0','1') )) {
        stop( cat("dataset", k, ": Censoring status has to be either 0 or 1")) 
      }
      
    }
  }  
}

#----------------------------------------------------------------------------#
# check meta methods
#----------------------------------------------------------------------------#

check.metamethod<- function(x,y, resp.type,ind.method,meta.method,rth=NULL,
                            REM.type=NULL,paired=NULL) {
  ## x is the data, y is the response
  K<-length(x)
  cat("Please make sure the following is correct:\n")
  cat("*You input",K,"studies\n")
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1) {
    cat("*You selected",ind.method,"for your",K,"studies respectively\n")
  }
  if (sum(paired)>0) cat("*Some of the studies are paired design\n")
  if (is.null(paired)) cat("*They are not paired design\n")
  cat("*",meta.method, "was chosen to combine the",K,"studies,respectively\n")
  if (length(meta.method)>1 & 
      sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))>0) {
    stop("Sorry, we currently do not allow multiple choices of meta.method 
         for 'FEM','REM','rankProd','minMCC'")
  }
  
  if (resp.type %in% c("continuous", "survival")) {
    if (meta.method %in% c("FEM","REM","minMCC","rankProd") ) {
      stop("The meta.methods of 'FEM','REM','rankProd','minMCC' are not 
              applicable to 'continuous','survival' response") 
    }
  } 
  
  #  if(sum(meta.method%in%c("maxP.OC","minP.OC","Fisher.OC","roP.OC",
  #                          "Stouffer.OC"))>0) {
  #  	  if (tail == c("abs") ) {
  #         stop("One-side corrected method ") 
  #  	   }
  #  }	
  
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1) {
    for(i in 1:K){
      check.label.for.meta(y[[i]],method=ind.method[i], 
                           meta.method=meta.method,k=i)
    }
  }
  
  if (length(grep('roP',meta.method))!=0&is.null(rth)) {
    stop("You should specify rth=XXX, when you choose roP method")
  }
  if (length(grep('REM',meta.method))!=0&is.null(REM.type)) {
    stop("You should specify REM.type=XXX, when you choose REM method")
  }
  if (!is.null(rth)){
    if (length(grep('roP',meta.method))!=0&&length(x)<rth) {
      stop("rth shouldn't be larger than the number of datasets")
    }
  } 
  
  if (!is.null(paired)&length(paired)<K) {
    stop(paste("you need to specify a vector of logical value for 'paired' 
               for all",K,"studies"))
  }
  
  ANOVA <- (resp.type == "multiclass")
  if (length(grep('maxP',meta.method))==0&&ANOVA&&meta.method!="minMCC") {
    warning("maxP is suggested for ANOVA model")
  }
  
  if (meta.method == "rankProd") {
    warning("rankProd method is time consuming")
  }
}

check.label.for.meta<-function(L,method,meta.method,k)
{
  if(!is.null(dim(L)))stop(cat("Please check whether the dimension in study",k,
                               " is matched to individual method",method,"?"))
  nL<-nlevels(as.factor(L))
  if ( nL<2) stop("All individual methods requires at least two levels")
  if (length(grep('minMCC',meta.method))!=0)
  {
    if (nL==2) warning("minMCC could test for two classes. But for better 
                      performance, try limma+maxP or sam+maxP")
    if (nL<2) stop(paste(meta.method,"minMCC method requires at least two 
                        levels"))
  }
  if (sum(table(L)<=1)==1) stop("<= one sample in the group, can not do the test 
                               or check labels")
  if (!is.null(method)&length(grep('minMCC',meta.method))!=0) 
    warning(paste("minMCC is a method that can not be combined with",method,". 
                  We'll perform minMCC only."))
  if (!is.null(method)&length(grep('rankProd',meta.method))!=0) 
    warning(paste("rankProd is a method that can not be combined with",method,".
                  We'll perform rankProd only."))
  
}


#------------------------------------------------------------------------------#
#perm.lab: function to permute the labels of disease status
# x: labels
# paired: a logical to specify whether the data is pair-designed or not
#------------------------------------------------------------------------------#
perm.lab<-function(x,paired=FALSE){
  New.x<-x	
  if(paired){
    templab<-rbinom(length(x)/2,1,0.5)
    index.d<-which(x==levels(x)[2])
    index.c<-which(x==levels(x)[1])
    dx<-x[index.d]
    cx<-x[index.c]
    New.dx<-dx
    New.cx<-cx
    New.dx[which(templab==1)]<-cx[which(templab==1)]
    New.cx[which(templab==1)]<-dx[which(templab==1)]
    New.x[index.d]<-New.dx
    New.x[index.c]<-New.cx
  }else{
    New.x<-sample(x)
  }
  return(New.x)
}

#=============================================================================#
# summary the DE number in a table
# pm: the p-value matrix
# p.cut: a numeric vector of p-values at which the DE numbers are counted 
# q.cut: a numeric vector of q-values at which the DE numbers are counted
# method: a vector of character string specifying the method
#-----------------------------------------------------------------------------#
count.DEnumber<-function(result,p.cut,q.cut){
  if(class(result)=="MetaDE.pvalue"){
    pm<-cbind(result$ind.p,result$meta.analysis$pval) 
  }else if(class(result)=="MetaDE.ES"){
    pm<-cbind(result$meta.analysis$pval)
    colnames(pm)<-attr(result$meta.analysis,"meta.method")
  }else if(class(result)=="MetaDE.minMCC"){
    pm<-cbind(result$meta.analysis$pval)
    colnames(pm)<-attr(result$meta.analysis,"meta.method")
  }else{
    pm<-result
  }
  qm<-cbind(apply(pm,2,function(x)p.adjust(x,method="BH")))
  table.p<-matrix(NA,length(p.cut),ncol(pm))
  for(i in 1:length(p.cut)){
    table.p[i,]<-apply(pm,2,function(x)sum(x<=p.cut[i],na.rm=T))
  }
  table.q<-matrix(NA,length(q.cut),ncol(pm))
  for(i in 1:length(q.cut)){
    table.q[i,]<-apply(qm,2,function(x)sum(x<=q.cut[i],na.rm=T))
  }       
  rownames(table.p)<-paste("p=",p.cut,sep="")
  rownames(table.q)<-paste("FDR=",q.cut,sep="")
  colnames(table.p)<-colnames(table.q)<-colnames(pm)
  return(list(pval.table=table.p,FDR.table=table.q))
}

draw.DEnumber<-function(result,maxcut,mlty=NULL,mcol=NULL,mlwd=NULL,mpch=NULL,
                        FDR=TRUE){
  if(class(result)=="MetaDE.pvalue"){
    pm<-cbind(result$ind.p,result$meta.analysis$pval) 
  }else if(class(result)=="MetaDE.ES"){
    pm<-cbind(result$meta.analysis$pval)
    colnames(pm)<-attr(result$meta.analysis,"meta.method")
  }else if(class(result)=="MetaDE.minMCC"){
    pm<-cbind(result$meta.analysis$pval)
    colnames(pm)<-attr(result$meta.analysis,"meta.method")
  }else{
    pm<-result
  }
  method<-colnames(pm)
  if(FDR) pm<-cbind(apply(pm,2,function(x)p.adjust(x,method="BH")))
  maxp<-max(pm,na.rm=T)
  if (maxcut>maxp)
  {
    cat("Your maximum cut point exceeds the maximum of observed 
              p-value/FDR\n",
        "we will use",maxp,"as the maximum cut point\n")
    maxcut<-maxp
  }
  ns<-ncol(pm)
  ymax<-max(apply(pm,2,function(x)sum(x<=maxcut,na.rm=T)))
  if(is.null(mlty))mlty=1:ns
  if(is.null(mcol))mcol=1:ns
  if(is.null(mlwd))mlwd=rep(2,ns)
  if(is.null(mpch))mpch=1:ns
  xlab0<-ifelse(FDR,"FDR cut-off","p-value cut-off")
  #----------get an optimal place to draw the symbols--------#
  get.c<-function(cut,pm)
  {
    s<-apply(pm,2,function(x,y) sum(x<=y,na.rm=T),y=cut)        
    return(sum(dist(cbind(cut,s))))
  }
  mycut<-as.matrix(seq(0,maxcut,length=20))
  dis<-apply(mycut,1,get.c,pm=pm)
  minx.pos<-mycut[which.max(dis)]
  
  plot(c(0,maxcut),c(1,ymax),type='n',xlab=xlab0,ylab="Significant tests")
  for(i in 1:ns){		
    y.pos<-sum(pm[,i]<=minx.pos,na.rm=T)
    if(y.pos==0){x.pos<-minx.pos
    }else{
      x.pos<-sort(pm[,i])[y.pos]}
    points(x.pos,y.pos,pch=mpch[i],col=mcol[i],lwd=3)
    lines(sort(pm[,i]),rank(sort(pm[,i]),ties.method="max"),lty=mlty[i],
          col=mcol[i],lwd=mlwd[i])
    
  }
  legend("topleft",method,lty=mlty,lwd=mlwd,col=mcol,bty='n',pch=mpch)
}



#################### SingleStudyDE ####################
Indi.DE.Analysis <- function(data, clin.data, data.type, resp.type, 
                             response,covariate=NULL,ind.method, 
                             select.group=NULL, ref.level=NULL, 
                             paired=NULL,asymptotic=NULL,nperm=NULL, 
                             tail="abs", seed=12345,
                             ... ) {
  
  ## Call the packages required 
  
  #library(survival)
  #library(limma)
  #library(samr)
  #library(edgeR)
  #library(DESeq2)
  
  set.seed(seed)
  
  ##Input  
  ##data: list of data matrices, genes on the rows, samples on the columns;
  ##clin.data: list of clinical data frames, samples on the rows, clinical
  ##covariates on the columns;
  ##data.type = c("continuous","discrete");
  ##resp.type = c("twoclass", "multiclass", "continuous", "survival");
  ##response: a character (column name of clin.data) indicating the phenotype 
  ##of interest (survival (2 column names));
  ##ind.method = c("limma","sam","limmaVoom","edgeR","DESeq2","pearsonr",
  ##"spearmanr","logrank");
  ##select.group: specify the two group names for comparison;
  ##ref.level: specify the reference level of the group factor;
  ##paired: logical indicating whether paired design;
  ##asymptotic: logical whether asymptotic dist should be used;
  ##nperm: number of permutations;
  ##tail = c("low", "high", "abs");
  
  
  K<-length(data) # number of studies
  G<-nrow(data[[1]]) # number of matched genes
  response.list <- covariate.list <- vector('list',K)
  if(!is.null(covariate)) covariate.list <-NULL
  for (k in 1:K) {
    response.list[[k]] <- clin.data[[k]][,response]
    if(!is.null(covariate)) {
      covariate.list[[k]] <-	clin.data[[k]][,covariate]
    }
  }  
  ##check the individual method for the corresponding resp.type and data.type,
  ##in addition, check the consistency between resp.type and response.
  
  check.indmethod(response.list, resp.type, data.type, ind.method, 
                  select.group, tail, paired) 
  
  ##verify the correct specification of individual method
  ind.method<- match.arg(ind.method,c("limma","sam","limmaVoom","edgeR",
                                      "DESeq2","pearsonr","spearmanr","logrank"),
                         several.ok=TRUE)
  
  data <-check.exp(data)  # check the gene names for expression data
  
  ##Association with groups (DE, ANOVA, etc.)  
  if(resp.type %in% c("twoclass","multiclass") )  {
    log2FC<- lfcSE <-  pvalue <- NULL #prespecify the output matrix
    N<-n<-NULL
    for (k in 1:K) { 
      groupLabel <- as.factor(response.list[[k]])
      if(!is.null(select.group)) {
        l<- groupLabel[groupLabel %in% select.group]
        y<-data[[k]][,groupLabel %in% select.group]
        l<- droplevels(l)
        if(is.null(covariate.list[[k]])) {
          c<- NULL
        }      
        if(!is.null(covariate.list[[k]]) && is.vector(covariate.list[[k]]))  {
          c <- covariate.list[[k]][groupLabel %in% select.group]
        } 
        if(!is.null(covariate.list[[k]]) && !is.vector(covariate.list[[k]])) {
          c <- covariate.list[[k]][groupLabel %in% select.group,]
        }
      } else{
        l<- groupLabel
        y<- data[[k]]
      }
      if(!is.null(ref.level)) {
        l <- relevel(l,ref=ref.level)
      } 
      
      ANOVA<- (resp.type == "multiclass") #indicate whether ANOVA should be performed          
      ind.res<-switch(ind.method[k],
                      limma={get.limma(y=y,l=l,c=c,name=name, 
                                       ANOVA=ANOVA,tail=tail,paired=paired[k])},
                      sam={get.sam(y=y,l=l,name=name,seed=seed,
                                   ANOVA=ANOVA,tail=tail,paired=paired[k])},
                      limmaVoom={get.limmaVoom(y=y,l=l,c=c,name=name,
                                               ANOVA=ANOVA,tail=tail,paired=paired[k])},
                      edgeR={get.edgeR(y=y,l=l,c=c,name=name,
                                       ANOVA=ANOVA,tail=tail,paired=paired[k])},
                      DESeq2={get.DESeq2(y=y,l=l,c=c,name=name,
                                         ANOVA=ANOVA,tail=tail,paired=paired[k])})
      
      pvalue<- cbind(pvalue,ind.res$pvalue) #p value
      
      if(!ANOVA) {
        if(any(grepl('log2FC',names(ind.res))) ) {
          log2FC<- cbind(log2FC,ind.res$log2FC) #log2 fold change 
        } else {
          log2FC<- cbind(log2FC,rep(NA,nrow(y)))
        }  
        
        if(any(grepl('lfcSE',names(ind.res))) )  {
          lfcSE <- cbind(lfcSE,ind.res$lfcSE) #standard error of log2FC 
        } else {
          lfcSE <- cbind(lfcSE,rep(NA,nrow(y)))
        }  
      }
      
      cat("dataset",k,"is done\n")
      nns<-get.sample.label.number(l)
      N<-c(N,nns$N) #sample size per study 
      n<-append(n,nns$n) #sample size per label per study
    } # end of K study for loop
    
    if(!ANOVA) { 	  
      if(is.null(names(data))) {
        colnames(log2FC) <-colnames(lfcSE)<-colnames(pvalue)<-
          paste("dataset",1:K,sep="") # naming studies
      } else {
        colnames(log2FC) <-colnames(lfcSE)<-colnames(pvalue)<-names(data) 
      }
      
      if(is.null(rownames(data[[1]])) ) {
        rownames(log2FC) <- rownames(lfcSE)<-rownames(pvalue)<-
          paste("gene",1:G,sep="") # naming genes
      } else {
        rownames(log2FC) <-rownames(lfcSE)<-rownames(pvalue)<-rownames(data[[1]]) 
      }
      all.res<-list(log2FC=log2FC,lfcSE=lfcSE,p=pvalue) 
      
    } else {
      all.res<-list(p=pvalue)
    }
    attr(all.res$p,"nperstudy")<-N
    attr(all.res$p,"nperlabelperstudy")<-n 
    attr(all.res$p,"data.type")<- data.type
    attr(all.res$p,"individual.analysis")<-ind.method
    attr(all.res$p,"response.type")<- resp.type
    attr(all.res$p,"alter.hypothesis")<- tail
    return(all.res)
  } #end of twoclass, multiclass comparison   
  
  if(resp.type == c("continuous") )  {
    P<-stat<-P.perm<-NULL #prespecify the output matrix
    N<- NULL
    for (k in 1:K) {   
      y<- data[[k]] 
      r<- response.list[[k]]  
      ind.res<-switch(ind.method[k],
                      pearsonr={get.p.pearsonr(y,r,asymptotic=asymptotic, 
                                               nperm=nperm,tail=tail)},
                      spearmanr={get.p.spearmanr(y,r,asymptotic= asymptotic,
                                                 nperm=nperm,tail=tail)})  
      P<- cbind(P,(ind.res$obs)[,2]) #p value
      stat<-cbind(stat,(ind.res$obs[,1])) # observed stats
      if (!is.null(nperm)) {
        P.perm<-cbind(P.perm,c(ind.res$perm)) #p value from permutation 
      }
      cat("dataset",k,"is done\n")
      N<-c(N,ncol(y))	  
    } # end of K study for loop
    
    if(is.null(names(data))) {
      colnames(stat)<-colnames(P)<-paste("dataset",1:K,sep="") # naming studies
    } else {
      colnames(stat)<-colnames(P)<-names(data) 
    }
    if (!is.null(nperm)) {
      colnames(P.perm)<-paste("dataset",1:K,sep="") # naming studies
    }
    
    all.res<-list(stat=stat,p=P,bp=P.perm)
    attr(all.res$p,"nperstudy")<-N
    attr(all.res$p,"data.type")<- data.type
    attr(all.res$p,"individual.analysis")<-ind.method
    attr(all.res$p,"response.type")<- resp.type
    attr(all.res$p,"alter.hypothesis")<- tail
    return(all.res)
  } #end of association with continuous phenotype
  
  if(resp.type == c("survival") )  {
    P<-stat<-P.perm<-NULL #prespecify the output matrix
    N<- NULL
    for (k in 1:K) {  
      y<- data[[k]] 
      r<- response.list[[k]]    
      ind.res<- get.p.logrank(y,r[,1], r[,2],asymptotic= asymptotic,
                              nperm=nperm,tail=tail)     
      P<- cbind(P,(ind.res$obs)[,2]) #p value
      stat<-cbind(stat,(ind.res$obs[,1])) # observed stats
      if (!is.null(nperm)) {
        P.perm<-cbind(P.perm,c(ind.res$perm)) #p value from permutation 
      }
      cat("dataset",k,"is done\n")
      N<-c(N,ncol(y))	    
    } # end of K study for loop
    
    if(is.null(names(data))) {
      colnames(stat)<-colnames(P)<-paste("dataset",1:K,sep="") #naming studies
    } else {
      colnames(stat)<-colnames(P)<-names(data)
    }
    if (!is.null(nperm)) {
      colnames(P.perm)<-paste("dataset",1:K,sep="") # naming studies
    }
    all.res<-list(stat=stat,p=P,bp=P.perm)
    attr(all.res$p,"nperstudy")<-N
    attr(all.res$p,"data.type")<- data.type
    attr(all.res$p,"individual.analysis")<-ind.method
    attr(all.res$p,"response.type")<- resp.type
    attr(all.res$p,"alter.hypothesis")<- tail
    return(all.res)    
  } #end of association with survival 
  
} # end of Indi.DE.Analysis function 

## Output 
## twoclass: pvalue, log2FC and its SE (if available) 
## multiclass: pvalue 
## continuous, survival: stat, p, bp (perm.p)


get.sample.label.number<-function(lbl) {
  N<-length(lbl) #sample size per study 
  n<- table(lbl) #sample size per label per study
  return(list(N=N,n=n))
}

get.limma<-function(y,l,c,name,ANOVA,tail, paired){
  ## y: intensity matrix, l: group label, c: clinical data, name: group name
  if(!ANOVA) {
    if(is.null(c)) {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        design <-model.matrix(~ l + subject) 
      } else{
        design <-model.matrix(~ l)  # design matrix
      } 
    } else {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        lc <- data.frame(l=l,c=c,s=subject)
        design <-model.matrix(~ l + c + s,data=lc) 
      } else{
        lc <- data.frame(l=l,c=c)
        design <-model.matrix(~ l + c,data=lc) # design matrix
      }    
    }  
    #log2FC <- rowMeans(y[,l==name[1]])-rowMeans(y[,l==name[2]]) 
    fit <-lmFit(y, design)
    ebFit<-eBayes(fit)
    out.table <- topTable(ebFit,coef=2, number=Inf, sort.by='none')
    log2FC <- out.table$logFC
    lfcSE <- sqrt(ebFit$s2.post) * fit$stdev.unscaled[,2]
    #stat <- out.table$t
    p <- as.numeric(out.table$P.Value)
    dir <- sign(log2FC)
    if (tail=="high") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(y/2) } else{
          return(1-y/2)	
        }, x=dir, y=p)
    }
    if (tail=="low") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(1-y/2) } else{
          return(y/2)	
        }, x=dir, y=p)
    }    	 
    
    if (tail=="abs")  pvalue <- p  
    
    #pvalue <-out.table$P.Value
    limma.out <- list(log2FC,lfcSE,pvalue)
    names(limma.out) <- c('log2FC','lfcSE','pvalue')
  } else {
    design <-model.matrix(~ -1+l)  # design matrix
    colnames(design)[1:nlevels(l)] = make.names(levels(l))
    combination <- apply(combn(make.names(levels(l)),2),2,function(x) {
      paste(x,collapse='-')
    } )
    cont.matrix<-makeContrasts(contrasts=combination,levels = colnames(design))
    fit <-lmFit(y, design)
    fit2 <-contrasts.fit(fit, contrasts=cont.matrix)
    ebFit<-eBayes(fit2)
    pvalue <- ebFit$F.p.value
    limma.out <- list(pvalue)
    names(limma.out) <- c('pvalue')
  }
  return(limma.out)
}

get.sam<-function(y,l,name, seed, ANOVA, tail, paired) {
  ## y: intensity matrix, l: group label, name: group name 
  options(verbose=F)
  if(!ANOVA) {
    #log2FC <- rowMeans(y[,l==name[1]])-rowMeans(y[,l==name[2]]) 
    ## prepare the sam data
    sam.data<-list(x=y,y=l, geneid= 1:nrow(y), genenames=rownames(y),logged2=T)
    ## run the sam permutation test 
    samr.obj<-samr(sam.data,  resp.type=ifelse(paired,
                                               "Two class paired","Two class unpaired"), nperms=100, 
                   random.seed=seed)
    #samr.obj<-samr(sam.data,  resp.type="Two class unpaired", nperms=100, 
    #               random.seed=seed)
    p<- as.numeric(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
    log2FC <- as.numeric(log2(samr.obj$foldchange))
    dir <- sign(log2FC)
    if (tail=="high") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(y/2) } else{
          return(1-y/2)	
        }, x=dir, y=p)
    }
    if (tail=="low") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(1-y/2) } else{
          return(y/2)	
        }, x=dir, y=p)
    }    	 
    
    if (tail=="abs")  pvalue <- p  
    
    sam.out <- list(log2FC,pvalue)
    names(sam.out) <- c('log2FC','pvalue')
  } else {
    sam.data<-list(x=y,y=l, geneid= 1:nrow(y), genenames=rownames(y),logged2=T)
    samr.obj<-samr(sam.data,resp.type="Multiclass", nperms=100,random.seed=seed)
    pvalue<- as.numeric(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
    sam.out <- list(pvalue)
    names(sam.out) <- c('pvalue')
  }
  options(verbose=T)
  return(sam.out)
}

get.limmaVoom <- function(y,l,c,name, ANOVA, tail, paired) {
  ## y: count matrix, l: group label, c: clinical data, name: group name 
  dge <- DGEList(counts=y)
  dge <- calcNormFactors(dge)  
  if(!ANOVA){
    if(is.null(c)) {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        design <-model.matrix(~ l + subject) 
      } else{
        design <-model.matrix(~ l)  # design matrix
      }  
    } else {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        lc <- data.frame(l=l,c=c,s=subject)
        design <-model.matrix(~ l + c + s,data=lc) 
      }else{
        lc <- data.frame(l=l,c=c)
        design <-model.matrix(~ l + c,data=lc) # design matrix
      }    
    }  
    v <- voom(dge,design,plot=FALSE,normalize="quantile") # voom normalization
    fit <-lmFit(v, design)
    ebFit<-eBayes(fit)
    out.table <- topTable(ebFit,coef=2, number=Inf, sort.by='none')
    log2FC <- out.table$logFC
    lfcSE <- sqrt(ebFit$var.prior)
    #stat <- out.table$t
    p <- as.numeric(out.table$P.Value)
    dir <- sign(log2FC)
    if (tail=="high") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(y/2) } else{
          return(1-y/2)	
        }, x=dir, y=p)
    }
    if (tail=="low") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(1-y/2) } else{
          return(y/2)	
        }, x=dir, y=p)
    }    	 
    
    if (tail=="abs")  pvalue <- p  
    
    limmaVoom.out <- list(log2FC,lfcSE,pvalue)
    names(limmaVoom.out) <- c('log2FC','lfcSE','pvalue')
  } else {
    design <-model.matrix(~ -1+l) 
    colnames(design)[1:nlevels(l)] = make.names(levels(l))
    combination <- apply(combn(make.names(levels(l)),2),2,function(x) {
      paste(x,collapse='-')
    } )
    cont.matrix<-makeContrasts(contrasts=combination,levels = colnames(design) )
    v <- voom(dge,design,plot=FALSE,normalize="quantile") # voom normalization
    fit <-lmFit(v, design)
    fit2 <-contrasts.fit(fit, contrasts=cont.matrix)
    ebFit<-eBayes(fit2)
    pvalue <- ebFit$F.p.value
    limmaVoom.out <- list(pvalue)
    names(limmaVoom.out) <- c('pvalue')
  }
  return(limmaVoom.out)
}

get.edgeR <- function(y,l,c,name, ANOVA, tail, paired) {
  ## y: count matrix, l: group label, c: clinical data, name: group name 
  dat <- DGEList(counts=y)
  dat <- calcNormFactors(dat)
  dat=estimateGLMCommonDisp(dat) #estimate common dispersion
  dat=estimateGLMTrendedDisp(dat) #estiamte trended dispersion
  dat=estimateGLMTagwiseDisp(dat) #estimate tagwise dispersion
  #dispersion <- dat$tagwise.dispersion
  
  if(!ANOVA){
    if(is.null(c)) {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        design <-model.matrix(~ l + subject) 
      } else{
        design <-model.matrix(~ l)  # design matrix
      }  
    } else {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        lc <- data.frame(l=l,c=c,s=subject)
        design <-model.matrix(~ l + c + s,data=lc) 
      }else{
        lc <- data.frame(l=l,c=c)
        design <-model.matrix(~ l + c,data=lc) # design matrix
      }    
    }      
    fit <- glmFit(dat, design)
    lrt <- glmLRT(fit, coef=2)
    out.table <-  lrt$table
    log2FC <- out.table$logFC
    
    p <- as.numeric(out.table$PValue)
    dir <- sign(log2FC)
    if (tail=="high") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(y/2) } else{
          return(1-y/2)	
        }, x=dir, y=p)
    }
    if (tail=="low") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(1-y/2) } else{
          return(y/2)	
        }, x=dir, y=p)
    }    	 
    
    if (tail=="abs")  pvalue <- p  
    
    edgeR.out <- list(log2FC,pvalue)
    names( edgeR.out) <- c('log2FC','pvalue')
  } else {
    if(is.null(c)) {
      design <-model.matrix(~ l)  # design matrix
    } else{
      lc <- data.frame(l=l,c=c)
      design <-model.matrix(~ l + c,data=lc) # design matrix
    }      
    fit <- glmFit(dat, design)
    lrt <- glmLRT(fit, coef=2:(nlevels(l)))
    out.table <-  lrt$table
    #log2FC <- out.table$logFC
    pvalue <- out.table$PValue
    edgeR.out <- list(pvalue)
    names(edgeR.out) <- c('pvalue') 
  }
  return(edgeR.out)
}

get.DESeq2 <- function(y,l,c,name, ANOVA, tail, paired) {
  ## y: count matrix, l: group label, c: clinical data, name: group name 
  #library(SummarizedExperiment)  
  if(!ANOVA){   
    if(is.null(c)) {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        design <-model.matrix(~ l + subject)  # design matrix
        colData <- data.frame(l=l, s=subject)
        colnames(colData) <-colnames(design)[-1]         
      }else{
        design <-model.matrix(~ l)  # design matrix
        colData <- data.frame(l=l)
        colnames(colData) <-colnames(design)[-1] 
      }  
    } else {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        lc <- data.frame(l=l,c=c,s=subject) 
        design <-model.matrix(~ l + c + s,data=lc)  # design matrix
        colData <- lc
        colnames(colData) <-colnames(design)[-1]           
      } else{
        lc <- data.frame(l=l,c=c) 
        design <-model.matrix(~ l + c,data=lc)  # design matrix
        colData <- lc
        colnames(colData) <-colnames(design)[-1] 
      } 
    }
    
    ddsMat <- DESeqDataSetFromMatrix(countData = y,
                                     colData = colData,
                                     design = as.formula(
                                       paste(" ~ ",paste(colnames(colData), collapse=" + ") ) 
                                     )  
    )
    ddsMat <- DESeq(ddsMat)
    res <- results(ddsMat,contrast=c(colnames(colData)[1],levels(l)[2],
                                     levels(l)[1]) )
    log2FC <- as.numeric(res$log2FoldChange)
    lfcSE <- as.numeric(res$lfcSE)
    
    p <- as.numeric(res$pvalue)
    dir <- sign(log2FC)
    if (tail=="high") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(y/2) } else{
          return(1-y/2)	
        }, x=dir, y=p)
    }
    if (tail=="low") { 
      pvalue <- mapply(function(x,y) if(x==1){
        return(1-y/2) } else{
          return(y/2)	
        }, x=dir, y=p)
    }    	 
    
    if (tail=="abs")  pvalue <- p  
    
    DESeq2.out <- list(log2FC,lfcSE,pvalue)
    names(DESeq2.out) <- c('log2FC','lfcSE','pvalue')  
  } else {     
    if(is.null(c)) {
      design <-model.matrix(~ l) 
      colData <- data.frame(group=l)
    } else{
      lc <- data.frame(l=l,c=c)
      design <-model.matrix(~ l + c,data=lc)
      colData <- lc
      colnames(colData) <-colnames(design)[-1] 
    }
    
    ddsMat <- DESeqDataSetFromMatrix(countData = y,
                                     colData = colData,
                                     design = as.formula(
                                       paste(" ~ ",paste(colnames(colData), collapse=" + ") )
                                     )  
    )    
    ddsMat <- DESeq(ddsMat,test='LRT',reduced= as.formula(
      paste(" ~ ",paste(colnames(colData)[-c(1:(nlevels(l)-1))], collapse=" + "))  
    )  
    )
    res <- results(ddsMat)
    pvalue <-  as.numeric(res$pvalue)
    DESeq2.out <- list(pvalue)
    names(DESeq2.out) <- c('pvalue')      
  }  
  return(DESeq2.out)
}

#--------calculate r statistic (pearson's correlation)for all genes-----#
#---note: l is continuous------------#
cal.pearsonr<-function(y,l) {
  stopifnot(length(l)==ncol(y), is.matrix(y))
  n<-length(l)
  num<-n*y%*%l-rowSums(y)*sum(l)
  den<-sqrt(n*rowSums(y^2)-rowSums(y)^2)*sqrt(n*sum(l^2)-sum(l)^2)
  r<-num/den
  r
}


#-------get p value for pearson correlation r using permutation-------#

get.p.pearsonr<-function(y,l,asymptotic, nperm,tail) {
  ## y: expression matrix , l: continuous outcome, tail: specify the Ha
  rnum<-which(apply(y,1,function(x) !any(is.na(x))))
  stat.obs<-p<-q<-rep(NA,nrow(y))
  names(stat.obs)<-names(p)<-rownames(y)
  stat.obs[rnum]<-cal.pearsonr(y[rnum,],l) #observed stat
  if (asymptotic==F) {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)
    stat.perm[rnum,]<-replicate(nperm,cal.pearsonr(y[rnum,],l)) #stat from perm
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail) #p from perm
    #q[rnum]<-p.adjust(p[rnum],method="BH")
    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
  } else {
    n<-length(l)
    t<-stat.obs*sqrt((n-2)/(1-stat.obs^2)) #fisher's transformation
    if (tail=="low") p[rnum]<-pt(t,n-2)
    if (tail=="high") p[rnum]<-1-pt(t,n-2)
    if (tail=="abs") p[rnum]<-2*(pmin(pt(t,n-2),1-pt(t,n-2)))
    res<-list(obs=cbind(stat.obs,p))
  }
  return(res)
}

#--------calculate r statistic (spearman's correlation)for all genes------#
#---note: l is continuous------------#
cal.spearmanr<-function(y,l) {
  stopifnot(length(l)==ncol(y), is.matrix(y))
  n<-length(l)
  y<-t(apply(y,1,rank))
  l<-rank(l)
  num<-n*y%*%l-rowSums(y)*sum(l)
  den<-sqrt(n*rowSums(y^2)-rowSums(y)^2)*sqrt(n*sum(l^2)-sum(l)^2)
  r<-num/den
  r
}

#-------get p value for spearman correlation r using permutation----------#
get.p.spearmanr<-function (y,l,asymptotic, nperm,tail) {
  ## y: expression matrix , l: continuous outcome, tail: specify the Ha
  rnum<-which(apply(y,1,function(x) !any(is.na(x))))
  stat.obs<-p<-q<-rep(NA,nrow(y))
  names(stat.obs)<-names(p)<-rownames(y)
  
  stat.obs[rnum]<-cal.spearmanr(y[rnum,],l)#observed stat
  
  if (asymptotic==F)
  {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)
    
    stat.perm[rnum,]<-replicate(nperm,cal.spearmanr(y[rnum,],l))# stat from perm
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail) #p from perm
    #q[rnum]<-p.adjust(p[rnum],method="BH")
    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
  } else {
    n<-length(l)
    t<-stat.obs*sqrt((n-2)/(1-stat.obs^2)) #fisher's transformation
    if (tail=="low") p[rnum]<-pt(t[rnum],n-2)
    if (tail=="high") p[rnum]<-1-pt(t[rnum],n-2)
    if (tail=="abs") p[rnum]<-2*(pmin(pt(t[rnum],n-2),1-pt(t[rnum],n-2)))
    res<-list(obs=cbind(stat.obs,p))
  }
  return(res)
}

#-------get logrank z statistic-----------------#
cal.logrank<-function(y,time,event)
{
  get.stat<-function(x,time,event) {  
    stat<-summary(coxph(Surv(time,event)~x),method="breslow")$sctest[1]
    stat
  }
  z<-apply(y,1,get.stat,time=time,event=event)
  z
}
cal.p.logrank<-function(y,time,event)
{
  get.p<-function(x,time,event) {  
    p<-summary(coxph(Surv(time,event)~x,method="breslow"))$sctest[3]
    p
  }
  p<-apply(y,1,get.p,time=time,event=event)
  p
}

#-------get p value for logrank z using permutation-----------------#
get.p.logrank<-function(y,time,event,asymptotic,nperm,tail)
{
  rnum<-which(apply(y,1,function(x) !any(is.na(x))))
  stat.obs<-p<-pp<-rep(NA,nrow(y))
  
  stat.obs[rnum]<-cal.logrank(y[rnum,],time=time,event=event)
  
  names(stat.obs)<-names(p)<-rownames(y)
  if(asymptotic== F)
  {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)
    stat.perm[rnum,]<-replicate(nperm,cal.logrank(y[rnum,],time=time,
                                                  event=event))
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail)
    #q[rnum]<-p.adjust(p[rnum],method="BH")    
    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
  } else {
    p[rnum]<-cal.p.logrank(y[rnum,],time=time,event=event)
    if (tail=="low")  pp[rnum]<-p
    if (tail=="high") pp[rnum]<-1-p
    if (tail=="abs")  pp[rnum]<-2*(pmin(p,1-p))
    res<-list(obs=cbind(stat.obs,p=pp))
  }
  res
}

#################### MetaDE.pvalue ####################
MetaDE.pvalue <-function(x,meta.method,rth=NULL,parametric=TRUE) {
  #meta.method<-match.arg(meta.method,several.ok = TRUE)
  check.parametric(meta.method,parametric)
  K<-ncol(x$p)
  if (parametric) x$bp<-NULL     
  nm<-length(meta.method)
  meta.res<-list(stat=NA,pval=NA,FDR=NA,AW.weight=NA)
  meta.res$stat<-meta.res$pval<-meta.res$FDR<-matrix(NA,nrow(x$p),nm)
  for( i in 1:nm){
    temp<-switch(meta.method[i],
                 maxP={get.maxP(x$p,x$bp)},minP={get.minP(x$p,x$bp)},
                 Fisher={get.fisher(x$p,x$bp)},roP={get.roP(x$p,x$bp,rth=rth)},
                 AW={get.AW(x$p)},
                 Fisher.OC={get.fisher.OC(x$p,x$bp)},
                 maxP.OC={get.maxP.OC(x$p,x$bp)},
                 minP.OC={get.minP.OC(x$p,x$bp)},
                 roP.OC={get.roP.OC(x$p,x$bp,rth=rth)},
                 Stouffer={get.Stouff(x$p,x$bp)},
                 Stouffer.OC={get.Stouff.OC(x$p,x$bp)},
                 SR={get.SR(x$p,x$bp)},PR={get.PR(x$p,x$bp)})
    meta.res$stat[,i]<-temp$stat
    meta.res$pval[,i]<-temp$pval
    meta.res$FDR[,i]<-temp$FDR
    if(meta.method[i]=="AW"){
      meta.res$AW.weight<-temp$AW.weight
    }
  }
  colnames(meta.res$stat)<-colnames(meta.res$pval)<-colnames(meta.res$FDR)<-
    meta.method
  rownames(meta.res$stat)<-rownames(meta.res$pval)<-
    rownames(meta.res$FDR)<-rownames(x$p)   
  attr(meta.res,"nstudy")<-K
  attr(meta.res,"meta.method")<-meta.method 
  res<-list(meta.analysis=meta.res,ind.p=x$p)	 
  #class(res)<-"MetaDE.pvalue"
  return(res)
}

