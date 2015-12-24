library(migen)
library(MASS)
library(parmigene)

set.seed(11)

## simulation bivariate random normal variables ###########

N <-100
p<-1000

sC <- rep(NaN,p)
kP <- rep(NaN,p)
kL <- rep(NaN,p)

par(mfrow=c(2,2))

## no correlation ############################################

rho <- 0

for (i in 1:p){
  z<-mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),2,2))
  x<-z[,1]
  y<-z[,2]
  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y,k=sqrt(N))
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names = c("cmis", "knnmi","cmik"))
x1<--log(1-rho^2)/2
abline(h=x1,col="RED")
mtext(paste(signif(mean(abs(sC-x1)),digits = 2), "    ",
            signif(mean(abs(kP-x1)),digits=2), "    ",
            signif(mean(abs(kL-x1)),digits=2)),
      line=2,side=1,cex=0.8)

## moderately correlated ###################################

rho <- 0.5

for (i in 1:p){
  z<-mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),2,2))
  x<-z[,1]
  y<-z[,2]
  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y,k=sqrt(N))
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names=c("cmis", "knnmi","cmik"))
x1<--log(1-rho^2)/2
abline(h=x1,col="RED")
mtext(paste(signif(mean(abs(sC-x1)),digits = 2), "    ",
            signif(mean(abs(kP-x1)),digits=2), "    ",
            signif(mean(abs(kL-x1)),digits=2)),
      line=2,side=1,cex=0.8)

## highly correlated ############################################

rho <- 0.8

for (i in 1:p){
  z<-mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),2,2))
  x<-z[,1]
  y<-z[,2]
  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y,k=sqrt(N))
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names=c("cmis", "knnmi","cmik"))
x1<--log(1-rho^2)/2
abline(h=x1,col="RED")
mtext(paste(signif(mean(abs(sC-x1)),digits = 2),"    ",
            signif(mean(abs(kP-x1)),digits=2),"    ",
            signif(mean(abs(kL-x1)),digits=2)),
      line=2,side=1,cex = 0.8)

## very highly correlated ###################################

rho <- 0.95

for (i in 1:p){
  z<-mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),2,2))
  x<-z[,1]
  y<-z[,2]
  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y,k=sqrt(N))
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names=c("cmis", "knnmi","cmik"))
x1<--log(1-rho^2)/2
abline(h=x1,col="RED")
mtext(paste(signif(mean(abs(sC-x1)),digits = 2), "    ",
            signif(mean(abs(kP-x1)),digits=2), "    ",
            signif(mean(abs(kL-x1)),digits=2)),
      line=2,side=1,cex = 0.8)

## Exponential #######################
par(mfrow=c(1,2))

## independent exponential ####################################

sC<-rep(NaN,p)
kP<-rep(NaN,p)
kL<-rep(NaN,p)

for (i in 1:p){
  x<-rexp(N)
  y <- rexp(N)

  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y,k=sqrt(N))
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names = c("cmis", "knnmi","cmik"))
x1<-0
abline(h=x1,col="RED")
mtext(paste(signif(mean(abs(sC-x1)),digits = 2), "     ",
            signif(mean(abs(kP-x1)),digits=2),"     ",
            signif(mean(abs(kL-x1)),digits=2)),
      line=2,side=1)

## exponential dependent #################################

sC<-rep(NaN,p)
kP<-rep(NaN,p)
kL<-rep(NaN,p)

for (i in 1:p){
  x <- rexp(N)
  y <- x+rexp(N)

  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y,k=sqrt(N))
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names = c("cmis", "knnmi","cmik"))
x2= 1-digamma(2)
abline(h=x2,col="RED")
mtext(paste(signif(mean(abs(sC-x2)),digits = 2), "     ",
            signif(mean(abs(kP-x2)),digits=2), "     ",
            signif(mean(abs(kL-x2)),digits=2)),
      line=2,side=1)


###### Time comparisons #########################
X <- rexp(N)
Y <- x+rexp(N)

system.time(for (i in 1:1000) cmis(X, Y))
system.time(for (i in 1:1000) knnmi(X, Y,k=sqrt(N)))
system.time(for (i in 1:1000) cmik(X, Y))




