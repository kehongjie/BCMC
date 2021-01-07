
## FUNCTION: repeat a row vector n times to a matrix
rep.row <- function(x,n) {
  matrix(rep(x,each=n),nrow=n)
}

## FUNCTION: extract numbers from a string
numextract <- function(string) { 
  stringr::str_extract(string, "\\-*\\d+\\.*\\d*")
} 


## FUNCTION: calculate Tw+ 
#### beta: effect size estimates (log2FC)
#### w: weight of 0 or 1
#### p: p-values
sim_posTw <- function(beta,w,p) {
  ## Tw = sum(w_k*beta_k*w_k'*beta_k'*|log10p_k+log10p_k'|)/sum(w)
  beta_new <- beta[which(beta>0)]
  beta.index <- as.numeric(as.character(numextract(names(beta_new))))
  p_new<-p[beta.index]
  w_new <- w[beta.index]
  
  w_out <- outer(w_new,w_new,FUN="*")
  w_pair <- w_out[upper.tri(w_out)]
  
  beta_out <- outer(beta_new,beta_new,FUN="*")
  beta_pair <- beta_out[upper.tri(beta_out)]
  
  p_out <- outer(log10(p_new),log10(p_new),FUN="+")
  p_pair <- p_out[upper.tri(p_out)]
  
  part1<-sum(w_pair*beta_pair*abs(p_pair))
  part2<-sum(w_new) 
  return(part1/part2)
}

## FUNCTION: calculate Tw-
sim_negTw <- function(beta,w,p) {
  ## Tw = sum(w_k*beta_k*w_k'*beta_k'*|log10p_k+log10p_k'|)/sum(w)
  beta_new <- beta[which(beta<0)]
  beta.index <- as.numeric(as.character(numextract(names(beta_new))))
  p_new<-p[beta.index]
  w_new <- w[beta.index]
  
  w_out <- outer(w_new,w_new,FUN="*")
  w_pair <- w_out[upper.tri(w_out)]
  
  beta_out <- outer(beta_new,beta_new,FUN="*")
  beta_pair <- beta_out[upper.tri(beta_out)]
  
  p_out <- outer(log10(p_new),log10(p_new),FUN="+")
  p_pair <- p_out[upper.tri(p_out)]
  
  part1<-sum(w_pair*beta_pair*abs(p_pair))
  part2<-sum(w_new)
  return(part1/part2)
}

