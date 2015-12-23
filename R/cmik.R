# bandwidth version ##########

# library(Rcpp)
# sourceCpp("cmik.cpp")

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

    return(cmikCpp(bw, N, X, Y, kmax))
}

