##' BCMC (Biomarker Categorization in Meta-analysis by Concordance)
##' for biomarker detection and categorization.
##'
##' This function will first run the \code{MetaDE} to perform DE analysis, and
##' then implement the meta-analysis method BCMC (Biomarker Categorization
##' in Meta-analysis by Concordance) for biomarker detection and categorization.
##'
##' @import Biobase
##' @import impute
##' @import stringr
##' @import gtools
##' @import limma
##' 
##' @param data.exp: The gene expression data. It should be a list in which
##' every element represents a study. Within each element (study), the data is
##' supposed to be matrix where each row is a gene and each column is a sample.
##' @param data.clin: The clinical data about whether a sample is case or control.
##' It should be a list in which every element represents a study. Within each
##' element (study), the data is supposed to be a data frame with only one column
##' and each row is a sample (same order as in \code{data.exp}). In each data
##' frame, code control as 0 and case a 1.
##' @param meta.de: A logical value indicating whether the \code{MetaDE} should
##' be run before the BCMC. The default is TRUE. If FALSE, users need to provide
##' the log fold change and p-value from DE analysis as arguments \code{de.lfc}
##' and \code{de.pval}.
##' @param de.lfc: If \code{meta.de} is FALSE, this should be a matrix of log fold
##' change from DE analysis. Each row should be a gene and each column should be
##' a study. Otherwise, the deafult is TRUE and \code{MetaDE} will be run for
##' DE analysis before the BCMC.
##' @param de.pval: If \code{meta.de} is FALSE, this should be a matrix of p-value
##' from DE analysis. Each row should be a gene and each column should be
##' a study. Otherwise, the deafult is TRUE and \code{MetaDE} will be run for
##' DE analysis before the BCMC.
##'
##' @return A list with following components:
##' \item{lfc}{The matrix of log fold change from DE analysis. Each row is a gene
##' and each column is a study.}
##' \item{pval}{The matrix of p-value from DE analysis. Each row is a gene
##' and each column is a study.}
##' \item{pos.wp}{The matrix of predicted up-regulated weight pattern. Each row
##' is a gene and each column is a study.}
##' \item{neg.wp}{The matrix of predicted down-regulated weight pattern. Each row
##' is a gene and each column is a study.}
##' \item{max.wp}{The matrix of predicted dominant weight pattern. Each row
##' is a gene and each column is a study. This is found by comparing the two BCMC
##' statistics for up-regualted and down-regulated.}
##' \item{Rg}{A matrix of BCMC statistcs, as the \eqn{R_g} in the paper. Each row
##' is a gene. The first column is for up-regulated pattern, the second column
##' is for down-regulated pattern, and the third column is the larger one of the
##' previous two statistics. See the paper for more details.}
##' 
##' @export
##' @examples
##' data("SimulDE")
##' result_bcmc <- bcmc(data.exp=SimulDE$express, data.clin=SimulDE$clin)
##' names(result_bcmc)
##' head(result_bcmc$Rg)
##' head(result_bcmc$pos.wp)


bcmc <- function(data.exp, data.clin, meta.de=TRUE, de.lfc=NULL, de.pval=NULL) {

  K <- length(data.exp)
  n <- ncol(data.exp[[1]])
  G <- nrow(data.exp[[1]])
  studyNames <- 1:K

  ##### run MetaDE and AW #####
  if (meta.de) {
    sim_AWmeta <- MetaDE(data=data.exp, clin.data=data.clin,
                         data.type="continuous", resp.type = "twoclass",
                         response='label', covariate = NULL, ind.method=rep("limma", K), meta.method="AW",
                         select.group = c("0","1"), ref.level="0",
                         paired=rep(FALSE, K), rth=NULL, REM.type=NULL, tail='abs', parametric=TRUE)

    de.lfc <- sim_AWmeta[["ind.stat"]] ## gene*study matrix after DE
    colnames(de.lfc) <- paste("lfc.", 1:K, sep="")
    de.pval <- sim_AWmeta[["ind.p"]] ## gene*study matrix after DE
    colnames(de.pval) <- paste("pval.", 1:K, sep="")
  }

  ##### run BCMC #####
  sim_w <- c(0,1)
  ## 2^K combination of 0/1
  simW_random <- data.matrix(permutations(n=length(sim_w), r=K, v=sim_w, repeats.allowed = T))

  ## all possible weight patterns
  simw <- vector("list", 2^K) ## list of combintation 2^K, each is a matrix of G*K
  for (T_wp in 1:nrow(simW_random)){
    simw[[T_wp]]<-data.matrix(rep.row(simW_random[T_wp,], G))
    colnames(simw[[T_wp]]) <- studyNames
  }


  ## search for best weight in our method
  ## Apply sim_posTw function
  sim.data_posTw <- matrix(0, nrow=G, ncol=2^K) ## now a G*2^K matrix
  for (T_wp in 1:length(simw)){
    for (T_dataset in 1:G) {
      sim.data_posTw[T_dataset,T_wp] <- sim_posTw(de.lfc[T_dataset,], simw[[T_wp]][T_dataset,],
                                                  de.pval[T_dataset,])
    }
  }
  ## Apply sim_negTw function
  sim.data_negTw <- matrix(0, nrow=G, ncol=2^K)
  for (T_wp in 1:length(simw)){
    for (T_dataset in 1:G) {
      sim.data_negTw[T_dataset,T_wp] <- sim_negTw(de.lfc[T_dataset,], simw[[T_wp]][T_dataset,],
                                                  de.pval[T_dataset,])
    }
  }

  ## make all NA to be 0
  sim.data_posTw[is.na(sim.data_posTw)] <- 0
  sim.data_negTw[is.na(sim.data_negTw)] <- 0

  ## find best Tw+
  simTw_posmax<-apply(sim.data_posTw,1,max) ## the largest T+
  simTw_which.pomax<-apply(sim.data_posTw,1,which.max) ## which combination of weight maximize
  # sim.data_posTw<-cbind.data.frame(sim.data_posTw, pomax=simTw_posmax,
  #                                  which.pomax=simTw_which.pomax) ## now a G*(2^K+2) matrix
  # Twpos_res <- data.frame(pomax=simTw_posmax, which.pomax=simTw_which.pomax) ## a G*2 matrix of max Tw+ information
  Twpos_res <- matrix(NA, G, K+1) ## G*(K+1) matrix of max Tw+ and argmax weight
  colnames(Twpos_res) <- c("posmax", paste("study", 1:K, sep=""))
  Twpos_res[,1] <- simTw_posmax
  for (i in 1:G){
    # sim.data_posTw[i,c(35:39)]<-simW_random[sim.data_posTw$which.pomax[i], ]
    # colnames(sim.data_posTw)[35:39]<-c("poswp.1","poswp.2","poswp.3","poswp.4","poswp.5")
    Twpos_res[i,-1] <- simW_random[simTw_which.pomax[i], ] ## save the argmax weight
  }

  ## find best Tw-
  simTw_negmax<-apply(sim.data_negTw,1,max)
  simTw_which.nemax<-apply(sim.data_negTw,1,which.max)
  # sim.data_negTw <- cbind.data.frame(sim.data_negTw, nemax=simTw_negmax,
  #                                    which.nemax=simTw_which.nemax)
  # Twneg_res <- data.frame(nemax=simTw_negmax, which.nemax=simTw_which.nemax) ## a G*2 matrix of max Tw- information
  Twneg_res <- matrix(NA, G, K+1) ## G*(K+1) matrix of max Tw- and argmax weight
  colnames(Twneg_res) <- c("negmax", paste("study", 1:K, sep=""))
  Twneg_res[,1] <- simTw_negmax
  for (i in 1:G){
    # sim.data_negTw[i,c(35:39)]<-simW_random[sim.data_negTw$which.nemax[i],c(1:5)]
    # colnames(sim.data_negTw)[35:39]<-c("negwp.1","negwp.2","negwp.3","negwp.4","negwp.5")
    Twneg_res[i,-1] <- simW_random[simTw_which.nemax[i],]
  }

  ## compare Tw+ and Tw-
  Twmax_res <- matrix(NA, G, K+1) ## G*(K+1) matrix of max(Tw+,Tw-) and argmax weight
  colnames(Twmax_res) <- c("pnmax", paste("study", 1:K, sep=""))
  for (i in 1:G){
    if (Twpos_res[i,1] > Twneg_res[i,1]) {
      Twmax_res[i,] <- Twpos_res[i,]
    } else {
      Twmax_res[i,] <- Twneg_res[i,]
    }
  }

  ##### summarize all the results #####
  Rg <- cbind(Twpos_res[,1], Twneg_res[,1], Twmax_res[,1])
  colnames(Rg) <- c("Rg.pos", "Rg.neg", "Rg.max")
  output <- list()
  output[["lfc"]] <- de.lfc
  output[["pval"]] <- de.pval
  output[["pos.wp"]] <- Twpos_res[,-1]
  output[["neg.wp"]] <- Twneg_res[,-1]
  output[["max.wp"]] <- Twmax_res[,-1]
  output[["Rg"]] <- Rg

  return(output)
}

