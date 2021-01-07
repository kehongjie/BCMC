##' Calculate the permutation-based p-value and FDR adjusted p-value (also 
##' known as q-value) after BCMC.
##'
##' This function use a permutation-based test to calculate the p-value and 
##' FDR adjusted p-value (also known as q-value) after doing the BCMC.
##' Details can be found in the paper. 
##'
##' @import parallel
##' 
##' @param B: Number of permutations. The default is 100.
##' @param data.exp: The gene expression data. It should be a list in which
##' every element represents a study. Within each element (study), the data is
##' supposed to be matrix where each row is a gene and each column is a sample.
##' @param data.clin: The clinical data about whether a sample is case or control.
##' It should be a list in which every element represents a study. Within each
##' element (study), the data is supposed to be a data frame with only one column
##' and each row is a sample (same order as in \code{data.exp}). In each data
##' frame, code control as 0 and case a 1. 
##' @param res.bcmc: The results of running BCMC. Typically the output of 
##' \code{bcmc} function. 
##' @param parallel: A logical value indicating whether the parallel computing
##' should be done. The default is FALSE. 
##' @param num.cores: If \code{parallel} is TRUE, this is the number of cores 
##' to be used in the paralle computing.
##' 
##' @return A list with following components:
##' \item{pvalue}{A matrix of permutation p-values where each row is a gene. The 
##' first column is for up-regulated statistc, the second column is for down
##' -regulated statistic, and the third colum is for the dominant statistic.}
##' \item{qvalue}{A matrix of permutation q-values where each row is a gene. The 
##' first column is for up-regulated statistc, the second column is for down
##' -regulated statistic, and the third colum is for the dominant statistic.}
##' 
##' @export
##' @examples
##' data("SimulDE")
##' result_bcmc <- bcmc(data.exp=SimulDE$express, data.clin=SimulDE$clin)
##' result_perm <- perm.bcmc(B=5, data.exp=SimulDE$express, data.clin=SimulDE$clin, 
##'                          res.bcmc=result_bcmc, parallel=FALSE)
##' names(result_perm)
##' head(result_perm$pvalue)


perm.bcmc <- function(B=100, data.exp, data.clin, res.bcmc, 
                      parallel=FALSE, num.cores=NULL) {
  
  K <- length(data.exp)
  n <- ncol(data.exp[[1]])
  G <- nrow(data.exp[[1]])
  studyNames <- 1:K
  
  ## permutate the data
  permB_list_data <- vector("list", B)
  permB_list_label <- vector("list", B)
  for (b in 1:B){
    for (s in 1:K) {
      permB_list_data[[b]][[s]] <- data.exp[[s]][,sample(n)]
      # rownames(permB_list_data[[b]][[s]]) <- c(1:G)
      permB_list_label[[b]][[s]] <- data.clin[[s]] 
    }
  }
  
  if (parallel==TRUE) {
    res_perm <- parallel::mcmapply(FUN = bcmc, data.exp=permB_list_data, data.clin=permB_list_label, 
                         mc.cores = numCores, SIMPLIFY = F)
  } else {
    res_perm <- mapply(FUN = bcmc, data.exp=permB_list_data, data.clin=permB_list_label, 
                       SIMPLIFY = F)
  }

  pqval <- calc.pq(res.bcmc=res.bcmc, res.perm=res_perm) ## permutation p&q-values
  return(pqval)
}



calc.pq <- function(res.bcmc, res.perm) {
  B <- length(res.perm)
  G <- nrow(res.perm[[1]]$lfc)
  
  ############### for best Tw+ ###############
  ## p-value
  permB_posTwpvalue<-vector("list",B)
  for (b in 1:B){
    permB_posTwpvalue[[b]]<-matrix(NA,ncol=G,nrow=G)
    for (g in 1:G) {
      # one_g <- as.numeric(res.perm[[b]][["pomax"]] >= res.bcmc[["pomax"]][g])
      one_g <- as.numeric(res.perm[[b]][["Rg"]][,1] >= res.bcmc[["Rg"]][g,1])
      permB_posTwpvalue[[b]][g, ] <- one_g 
    }
  }
  part1_posG<-vector("list",G)
  part2_posGB<-rep(NA,G)
  posTwpvalue_G<-rep(NA,G) ## **p-value**, for each gene
  for (g in 1:G){
    part1_posG[[g]]<-rep(NA,B)
    for (b in 1:B){
      part1_posG[[g]][b]<-sum(permB_posTwpvalue[[b]][g,])
    }
    part2_posGB[g]<-sum(part1_posG[[g]])
    posTwpvalue_G[g]<-part2_posGB[g]/(B*G)
  }
  
  ## pi hat 
  pihat_posTw_index<-rep(NA,G)
  for (i in 1:G){
    if(posTwpvalue_G[i] >= 0.5 & posTwpvalue_G[i] <=1){
      pihat_posTw_index[i]<-1
    }else{
      pihat_posTw_index[i]<-0
    }
  }
  pihat_posTw<-sum(pihat_posTw_index)/(G*0.5) ## **pai0 hat**
  
  ## calculate q-value
  qvalue_Gpos<-rep(NA,G) ## q-value, for each gene
  part2Gpos<-vector("list",G)
  for (i in 1:G){
    # part2Gpos[[i]] <- as.numeric(res.bcmc[,27] >= res.bcmc[i,27])
    part2Gpos[[i]] <- as.numeric(res.bcmc[["Rg"]][,1] >= res.bcmc[["Rg"]][i,1])
    qvalue_Gpos[i]<-min(1,(pihat_posTw*part2_posGB[i])/(B*sum(part2Gpos[[i]])))
  }
  
  
  ############### for best Tw- ############### 
  ## p-value
  permB_negTwpvalue<-vector("list",B)
  for (b in 1:B){
    permB_negTwpvalue[[b]]<-matrix(NA,ncol=G,nrow=G)
    for (g in 1:G) {
      # one_g <- as.numeric(res.perm[[b]][["nemax"]] >= res.bcmc[["nemax"]][g])
      one_g <- as.numeric(res.perm[[b]][["Rg"]][,2] >= res.bcmc[["Rg"]][g,2])
      permB_negTwpvalue[[b]][g, ] <- one_g 
    }
  }
  part1_negG<-vector("list",G)
  part2_negGB<-rep(NA,G)
  negTwpvalue_G<-rep(NA,G)
  for (g in 1:G){
    part1_negG[[g]]<-rep(NA,B)
    for (b in 1:B){
      part1_negG[[g]][b]<-sum(permB_negTwpvalue[[b]][g,])
    }
    part2_negGB[g]<-sum(part1_negG[[g]])
    negTwpvalue_G[g]<-part2_negGB[g]/(B*G)
  }
  
  ## pi hat 
  pihat_negTw_index<-rep(NA,G)
  for (i in 1:G){
    if(negTwpvalue_G[i] >= 0.5 & negTwpvalue_G[i] <= 1){
      pihat_negTw_index[i]<-1
    }else{
      pihat_negTw_index[i]<-0
    }
    pihat_negTw<-sum(pihat_negTw_index)/(G*0.5)
  }
  
  ## calculate q-value
  qvalue_Gneg<-rep(NA,G)
  part2Gneg<-vector("list",G)
  for (i in 1:G){
    # part2Gneg[[i]] <- as.numeric(res.bcmc[,34] >= res.bcmc[i,34])
    part2Gneg[[i]] <- as.numeric(res.bcmc[["Rg"]][,2] >= res.bcmc[["Rg"]][i,2])
    qvalue_Gneg[i]<-min(1,(pihat_negTw*part2_negGB[i])/(B*sum(part2Gneg[[i]])))
  }
  
  
  ########## for best Tw between Tw- and Tw+ ########## 
  ## p-values
  # cal_maxrg <- apply(res.bcmc[,c(27,34)], 1, max)
  cal_maxrg <- as.vector(res.bcmc[["Rg"]][,3])
  permB_pnTwpvalue<-vector("list",B)
  for (b in 1:B){
    # maxrg_vec <- rbind(as.vector(res.perm[[b]][["pomax"]]), as.vector(res.perm[[b]][["nemax"]]))
    # maxrg_vec <- apply(maxrg_vec, 2, max)
    maxrg_vec <- as.vector(res.perm[[b]][["Rg"]][,3])
    permB_pnTwpvalue[[b]]<-matrix(NA,ncol=G,nrow=G)
    for (g in 1:G) {
      one_g <- as.numeric(maxrg_vec >= cal_maxrg[g])
      permB_pnTwpvalue[[b]][g, ] <- one_g 
    }
  }
  part1_pnG<-vector("list",G)
  part2_pnGB<-rep(NA,G)
  pnTwpvalue_G<-rep(NA,G)
  for (g in 1:G){
    part1_pnG[[g]]<-rep(NA,B)
    for (b in 1:B){
      part1_pnG[[g]][b]<-sum(permB_pnTwpvalue[[b]][g,])
    }
    part2_pnGB[g]<-sum(part1_pnG[[g]])
    pnTwpvalue_G[g]<-part2_pnGB[g]/(B*G)
  }
  
  ## pi hat 
  pihat_pnTw_index<-rep(NA,G)
  for (i in 1:G){
    if(pnTwpvalue_G[i] >= 0.5 & pnTwpvalue_G[i] <= 1){
      pihat_pnTw_index[i]<-1
    }else{
      pihat_pnTw_index[i]<-0
    }
    pihat_pnTw<-sum(pihat_pnTw_index)/(G*0.5)
  }
  
  ## calculate q-value
  qvalue_Gpn<-rep(NA,G)
  part2Gpn<-vector("list",G)
  for (i in 1:G){
    # maxrg_vec <- rbind(as.vector(res.bcmc[,27]), as.vector(res.bcmc[,34]))
    # maxrg_vec <- apply(maxrg_vec, 2, max)
    # part2Gpn[[i]] <- as.numeric(maxrg_vec >= max(res.bcmc[i,27],res.bcmc[i,34]))
    part2Gpn[[i]] <- as.numeric(res.bcmc[["Rg"]][,3] >= res.bcmc[["Rg"]][i,3])
    qvalue_Gpn[i]<-min(1,(pihat_pnTw*part2_pnGB[i])/(B*sum(part2Gpn[[i]])))
  }
  
  #################### return #################### 
  pvalue <- cbind(pos=posTwpvalue_G, neg=negTwpvalue_G, max=pnTwpvalue_G)
  qvalue <- cbind(pos=qvalue_Gpos, neg=qvalue_Gneg, max=qvalue_Gpn)
  output <- list(pvalue=pvalue, qvalue=qvalue)
  return(output)
  
}







