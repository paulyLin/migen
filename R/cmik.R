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


# Need to decide best sorting approach, including the possibility
# of partial sorting or just selecting nth largest element.
cmik2 <-function(X, 
                 Y,
                 bw = c(1, 1),
                 kmax = floor(sqrt(length(X))),
                 tiebreak = TRUE,
                 scale.data = TRUE)
{

    if (tiebreak)
    {

        dX <- duplicated(X)
        dY <- duplicated(Y)
        X[dX] <- X[dX] + rnorm(length(X[dX]), 0, 0.001)
        Y[dY] <- Y[dY] + rnorm(length(Y[dY]), 0, 0.001)  
    }

    if (scale.data)
    {
        X<-as.vector(scale(X))
        Y<-as.vector(scale(Y))
    }

    N<-length(X)

    k<- rep(NaN,N)
    # kN <- rep(NaN,N)
    eD <- rep(NaN,N)
    # pM <- rep(NaN,N)
    s <- rep(NaN,N)
    sN<-rep(NaN,N)
    t <- rep(NaN,N)
    tN<-rep(NaN,N)
    l <- rep(NaN,N)
    m <- rep(NaN,N)

    # Matrices of differences from each point
    # (Doesn't help speed in R)
    xdiffs <- abs(outer(X, X, "-"))
    ydiffs <- abs(outer(Y, Y, "-"))
    # These matrices are symmetrical (due to abs()), which can be used
    # to speed up the C++ code (as in mpmi Fortran code).

    # We could do a similar thing with the max norm distances.
    # Could be a big improvement due to symmetry.
    distmat <- pmax(xdiffs, ydiffs)
    # This gives a huge speed up in R.

    for (i in 1:N)
    {
        # dXi <- abs(X-X[i]) # == xdiffs[, i]
        # dYi <- abs(Y-Y[i])

        # N.B., if we order the distances at the beginning then
        # all the counting becomes much faster because we can break
        # the loop as soon as we exceed the distance (i.e., bw or eD).
        #
        # Can maybe include kmax in the loop used to find sN and tN.

        # Count number within bandwidth
        s[i] <- sum(xdiffs[, i] < bw[1])
        # Get index of s[i]th neighbour
        # (Can be outside bandwidth!)
        # I.e., this is index of closest point
        # outside bandwidth.

        # Stop problems latter:
        s[i] <- max(2, s[i])

        sN[i] <- order(xdiffs[, i])[s[i]] # Get rid of +1 here
        # for C++ probably best to just sort before finding s[i]

        # C++ can use nth_element !!!! This could be very fast.

        # Same as for X
        t[i] <- sum(ydiffs[, i] < bw[2])
        t[i] <- max(2, t[i])

        tN[i] <- order(ydiffs[, i])[t[i]]

        # M <- cbind(xdiffs[, i], ydiffs[, i])
        # Max norm distance to every other point (X,Y pair)
        # d <- apply(M, 1, max)
        d <- distmat[, i]
        
        # For each point, k is the number of points closer
        # to it (using the max norm) that are closer than
        # the distance to the closest point outside the bandwidth.

        # I.e., k is going to be either kmax, s[i] or t[i]?
        # (NO)

        xeD <- d[sN[i]]
        yeD <- d[tN[i]]

        k1 <- sum(d < xeD)
        k2 <- sum(d < yeD)

        if (kmax < k1 && kmax < k2)
        {
            k[i] <- kmax
            eD[i] <- sort(d, partial = 1:(k[i] + 1))[k[i] + 1] 
            # Can use C++ std::partial_sort
            # EVEN BETTER: use std::nth_element
            # Maybe there is sugar
            # (R's sorting is probably crap, but that needs to be
            # balanced against any overhead that might exist when 
            # applying a C++ std function to an R vector. Need to test.)
        } else
        {
            which_k <- which.min(c(k1, k2))
            k[i] <- c(k1, k2)[which_k]
            eD[i] <- c(xeD, yeD)[which_k]
        }

        # k[i] <- min(sum(d < d[sN[i]]), sum(d < d[tN[i]]), kmax)

        # In C++, can skip straight from X and Y to next step as soon
        # as we count neighbours > kmax.
       
        # Maybe we don't need s, sN, t and tN.

        # Alternatively we may not need further calculation to 
        # get eD.

        # kN[i] <- order(d)[k[i] + 1]
        # eD[i] <- sort(d)[k[i] + 1] 
        # if k!=kmax then eD[i] will be either d[sN[i]] or d[tN[i]]
        # pM[i] <- which.max(M[kN[i], ])
        l[i] <- sum(xdiffs[, i] < eD[i])
        m[i] <- sum(ydiffs[, i] < eD[i])
        # This can be faster if these are sorted
        # probably a shortcut here too (maybe not)
    }

    MI <- digamma(N) + mean(digamma(k)) - mean(digamma(l) + digamma(m))

    return(MI)
}
    cppout <- cmikCpp(bw, N, X, Y, kmax) 
    cppout$s
    cppout$t
    cppout$sN + 1 # C++ indexes start at 0.
    cppout$tN + 1 # C++ indexes start at 0.
    cppout$eD

library(Rcpp)
sourceCpp("cmik.cpp")

cmikc <-function(X, 
                 Y,
                 bw = c(1, 1),
                 kmax = floor(sqrt(length(X))),
                 tiebreak = TRUE,
                 scale.data = TRUE)
{

    if (tiebreak)
    {

        dX <- duplicated(X)
        dY <- duplicated(Y)
        X[dX] <- X[dX] + rnorm(length(X[dX]), 0, 0.001)
        Y[dY] <- Y[dY] + rnorm(length(Y[dY]), 0, 0.001)  
    }

    if (scale.data)
    {
        X<-as.vector(scale(X))
        Y<-as.vector(scale(Y))
    }

    N<-length(X)

    return(cmikCpp(bw, N, X, Y, kmax)$MI)
}

# Timing
x <- rnorm(1000)
y <- rnorm(1000)

system.time(cmik(x, y))
# Before factoring:
#    user  system elapsed 
#   1.416   0.000   1.414
# After (using dXi and dYi):
#    user  system elapsed 
#   1.387   0.000   1.384 
# 
# Hardly any difference

# Outer product doesn't seem to help either (WORSE!)
#   user  system elapsed 
#  1.453   0.000   1.452
