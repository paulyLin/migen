
## bandwidth version - varying bandwidth

mmik <-function(cts, disc, bw=rep(0,length(unique(disc))),skewnessCorrection=TRUE){
  
  N <- length(disc)
  
  if (N !=length(cts)){
    stop("disc and cts must have the same lengths")
  }
  disc<-factor(disc)
  
  disc <-disc[order(cts)]
  cts <-sort(cts)
  
  cts <- cts-min(cts)+0.5

# The fllowing function uses codes from   
#  http://rstudio-pubs-static.s3.amazonaws.com/1563_1ae2544c0e324b9bb7f6e63cf8f9e098.html

  if(skewnessCorrection){
    skew.score <- function(c, x) (skewness(log(x + c)))^2
    best.c <- optimise(skew.score, c(0, 20), x = cts)$minimum
    cts <- log(cts + best.c)
  }
  
  m <- length(levels(disc))
  Ni <- table(disc)

  disclevels <- list()
  for (i in 1:m){
    disclevels[[i]] <- which(disc==levels(disc)[i])
  }
  
  if(bw[1]==0){
    for(i in 1:m){
      bw[i]<-sd(cts[disclevels[[i]]])
    }
  }
  
  k<-rep(NaN,N)
  ka<-rep(NaN,N)
  d<-rep(NaN,N)
  for(i in 1:N){
    wl <- which(disc[i]==levels(disc))
    disci<- disclevels[[wl]]
    index <- which(cts[disci]==cts[i])
    disciPoints <- cts[disci][which(abs(cts[disci]-cts[i]) < bw[wl])]
    k[i] <- max(length(disciPoints) - 1, 1)
    d[i]<- max(abs(cts[i] - disciPoints))
    if(d[i]==0 | is.nan(d[i])){
      d[i]<-max(abs(cts[disci][c(index-1,index+1)]-cts[i]), na.rm = TRUE)
    }
    ka[i]<-sum(abs(cts-cts[i])<d[i])
  }

  P3<-mean(digamma(k)-digamma(ka))
  
  Nisum <- sum(Ni*(-digamma(Ni)))/N
  
  MI <- digamma(N) + Nisum+P3
  
  return(MI)
}
