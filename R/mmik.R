mmik <-function(cts, 
                disc, 
                bw = rep(0, length(unique(disc))),
                skewnessCorrection = TRUE,
                na.rm = FALSE)
{

  N <- length(disc)

  if (N != length(cts))
  {
    stop("cts and disc must have the same length")
  }
  disc <- factor(disc)

  if(na.rm)
  {
      okcts <- !is.na(cts)
      okdisc <- !is.na(disc)
      cts <- cts[okcts & okdisc]
      disc <- disc[okcts & okdisc]
  }

  sortOrder <- order(cts)
  disc <- disc[sortOrder]
  cts <- cts[sortOrder]

  cts <- cts - min(cts) + 0.5

# The following function is based on code from
#  http://rstudio-pubs-static.s3.amazonaws.com/1563_1ae2544c0e324b9bb7f6e63cf8f9e098.html

  if(skewnessCorrection && abs(skewness(cts)) > 0.5)
  {
    skew.score <- function(c, x) (skewness(log(x + c)))^2
    best.c <- optimise(skew.score, c(0, 20), x = cts)$minimum
    cts <- log(cts + best.c)
  }

  Ni <- table(disc)
  ctslevels <- split(cts, disc)

  return(mmikCpp(ctslevels, bw, N, cts, Ni))
}

