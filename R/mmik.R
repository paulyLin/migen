mmik <-function(cts, disc, bw=rep(0,length(unique(disc))),skewnessCorrection=TRUE){

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

  return(mmikCpp(ctslevels, bw, N, cts, Ni))
}

