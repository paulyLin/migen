mmik <-function(cts,disc,ki){
  
  N <- length(disc)
  
  if (N !=length(cts)){
    stop("X and Y must have the same lengths")
  }
  disc <- factor(disc)
  
  disc <-disc[order(cts)]
  cts <-sort(cts)
  
  m <- length(levels(disc))
  Ni <- table(disc)
  
#   for (i in 1:m) {
#     if (any(ki[[i]]>=Ni[i])) {
#       stop("each ki must be smaller than the size of the corresponding category")
#     }
#   }
  
  dkn <- function(S,k) {
    d<-knn.dist(S,k)[,k]
    Y2 <- setdiff(cts,S)
    ls<-rowSums(abs(sapply(Y2, function(x) x-S))<d)+k
    return(ls)
  }
  
  lsum <- rep(NA,m)
  for (i in 1:m){
    li <- which(disc==levels(disc)[i])
    lsum[i] <- sum(digamma(dkn(cts[li],ki[i])))
  }
  lsum <- sum(lsum)
  
  ksum <- sum(Ni*(digamma(ki)-digamma(Ni)))
  
  MI <- digamma(N) + (ksum - lsum)/N
  
  return(MI)
}
