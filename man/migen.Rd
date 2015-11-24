\name{migen}
\alias{migen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mutual Information - General
}
\description{
This package calcuates mutual information through either
nearest k neighbor algorithm or kernel smoothing.
}

\author{
Paul Lin, Chris Pardy
}

\examples{
## normal indpendent data
  w2 <- c(1,1,1)
  probs <- w2 / sum(w2)
  w <- rowSums(rmultinom(100, 1, prob=probs))
  X<-rnorm(sum(w),0,1)
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))

## Using k neighbor algorithm
  mpmik(X,Y)
## Using kernel smoothing
  mpmis(X,Y)
  
## US state data  
  x<-state.x77[,1]
  y<-state.x77[,2]
  plot(x,y)
  cor(x,y)
  cmik(x,y)
  cmis(x,y)
  
## US arrests data
  x<-USArrests[,1]
  y<-USArrests[,2]
  plot(x,y)
  cor(x,y)
  cmik(x,y)
  cmis(x,y)

## Spot Patterns
  x <-c(rnorm(50,-10,3), rnorm(50,0,3),rnorm(50,2,3))
  y <-c(rnorm(50,-5,3), rnorm(50,5,3),rnorm(50,-10,3)) 
  plot(x,y)
  cor(x,y)
  cmik(x,y)
  cmis(x,y)
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
