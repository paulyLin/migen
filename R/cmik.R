# bandwidth version ##########

cmik <-function(X,Y,bw=c(0.8,0.8),kmax=floor(sqrt(length(X)))){
  
  dX <- duplicated(X)
  dY <- duplicated(Y)
  X[dX]<-X[dX]+rnorm(length(X[dX]),0,0.001)
  Y[dY]<-Y[dY]+rnorm(length(Y[dY]),0,0.001)  
   
  X<-scale(X)
  Y<-scale(Y)

  N<-length(X)

  k<- rep(NaN,N)
  kN <- rep(NaN,N)
  eD <- rep(NaN,N)
  pM <- rep(NaN,N)
  s <- rep(NaN,N)
  sN<-rep(NaN,N)
  t <- rep(NaN,N)
  tN<-rep(NaN,N)
  l <- rep(NaN,N)
  m <- rep(NaN,N)
  
  for (i in 1:N){
    s[i]<-sum(abs(X-X[i])<bw[1])
    sN[i]<-order(abs(X-X[i]))[s[i]+1]
    t[i]<-sum(abs(Y-Y[i])<bw[2])
    tN[i]<-order(abs(Y-Y[i]))[t[i]+1]
    M <- cbind(abs(X-X[i]),abs(Y-Y[i]))
    d <- apply(M,1,max)
    k[i]<-min(sum(d<d[sN[i]]),sum(d<d[tN[i]]),kmax)
    kN[i]<-order(d)[k[i]+1]
    eD[i] <-sort(d)[k[i]+1]
    pM[i] <- which.max(M[kN[i],])
    l[i]<-sum(abs(X-X[i])<eD[i])
    m[i]<-sum(abs(Y-Y[i])<eD[i])
  }
    
 MI <- digamma(N)+mean(digamma(k))-mean(digamma(l)+digamma(m))
   
  return(MI)
}

