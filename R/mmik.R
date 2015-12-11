
## bandwidth version - varying bandwidth

mmik <-function(cts, disc, bw=rep(0,length(unique(disc))),skewnessCorrection=TRUE){
  
  N <- length(disc)
  
  if (N !=length(cts)){
    stop("disc and cts must have the same lengths")
  }
  disc<-factor(disc)
  
  sortOrder <- order(cts)
  disc <-disc[sortOrder]
  cts <-cts[sortOrder]
  
  cts <- cts-min(cts)+0.5

# The fllowing function uses codes from   
#  http://rstudio-pubs-static.s3.amazonaws.com/1563_1ae2544c0e324b9bb7f6e63cf8f9e098.html

  if(skewnessCorrection & abs(skewness(cts)) > 0.5){
    skew.score <- function(c, x) (skewness(log(x + c)))^2
    best.c <- optimise(skew.score, c(0, 20), x = cts)$minimum
    cts <- log(cts + best.c)
  }
  
  m <- length(levels(disc))
  Ni <- table(disc)
  
  ctslevels <- split(cts, disc)

  if(bw[1]==0){
    for(i in 1:m){
      bw[i]<-sd(ctslevels[[i]])
    }
  }
  
  k<-rep(NaN,N)
  ka<-rep(NaN,N)
  d<-rep(NaN,N)
  
  kk <- 1
  for (i in 1:m){
    wl <- ctslevels[[i]]
    for (j in 1:length(wl)){
      point <- wl[j]
      disciPoints <- wl[which(abs(wl-point) < bw[i])]
      k[kk] <- max(length(disciPoints) - 1, 1)
      d[kk]<- max(abs(point - disciPoints))
      if(d[kk]==0 | is.nan(d[kk])){
        d[kk]<-max(abs(wl[c(j-1,j+1)]-point), na.rm = TRUE)
      }
      ka[kk]<-sum(abs(cts-point)<d[kk]) # Is this d[i] a mistake?!?!?!?
      kk <- kk + 1
    }
  }

  P3<-mean(digamma(k)-digamma(ka))
  
  Nisum <- sum(Ni*(-digamma(Ni)))/N
  
  MI <- digamma(N) + Nisum+P3
  
  return(MI)
}

library(Rcpp)
sourceCpp("mmik.cpp")

mmikc <-function(cts, disc, bw=rep(0,length(unique(disc))),skewnessCorrection=TRUE){
  
  N <- length(disc)
  
  if (N !=length(cts)){
    stop("disc and cts must have the same length")
  }
  disc<-factor(disc)
  
  sortOrder <- order(cts)
  disc <-disc[sortOrder]
  cts <-cts[sortOrder]
  
  cts <- cts-min(cts)+0.5

# The fllowing function uses codes from   
#  http://rstudio-pubs-static.s3.amazonaws.com/1563_1ae2544c0e324b9bb7f6e63cf8f9e098.html

  if(skewnessCorrection & abs(skewness(cts)) > 0.5){
    skew.score <- function(c, x) (skewness(log(x + c)))^2
    best.c <- optimise(skew.score, c(0, 20), x = cts)$minimum
    cts <- log(cts + best.c)
  }
  
  m <- length(levels(disc))
  Ni <- table(disc)
  
  ctslevels <- split(cts, disc)

  return(mmikCpp(ctslevels, bw, N, cts, Ni))
}
