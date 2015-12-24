cmik <-function(X, 
                Y,
                bw = c(1, 1),
                kmax = floor(sqrt(length(X))),
                tiebreak = TRUE,
                scale.data = TRUE, 
                na.rm = FALSE)
{
    N <- length(X)

    if (N != length(Y))
    {
        stop("X and Y must have the same length")
    }

    if(na.rm)
    {
        okX <- !is.na(X)
        okY <- !is.na(Y)
        X <- X[okX & okY]
        Y <- Y[okX & okY]
    }

    if (tiebreak)
    {

        dX <- duplicated(X)
        dY <- duplicated(Y)
        X[dX] <- X[dX] + rnorm(length(X[dX]), 0, 0.001)
        Y[dY] <- Y[dY] + rnorm(length(Y[dY]), 0, 0.001)  
    }

    if (scale.data)
    {
        X <- as.vector(scale(X))
        Y <- as.vector(scale(Y))
    }

    return(cmikCpp(bw, N, X, Y, kmax))
}

