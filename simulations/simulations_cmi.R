# bandwidth version ###########

set.seed(11)

library(parmigene)

## simulation bivariate random normal variables ###########

N <-100
p<-500

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
  kP[i]<-knnmi(x,y)
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names = c("cmis", "knnmi","cmik"))
x1<--log(1-rho^2)/2
abline(h=x1,col="RED")
mtext(paste("CP:",signif(mean(abs(sC-x1)),digits = 2),
            "P:",signif(mean(abs(kP-x1)),digits=2),
            "PL",signif(mean(abs(kL-x1)),digits=2)),
      line=3,side=1)

## moderately correlated ###################################

rho <- 0.5

for (i in 1:p){
  z<-mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),2,2))
  x<-z[,1]
  y<-z[,2]
  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y)
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names=c("cmis", "knnmi","cmik"))
x1<--log(1-rho^2)/2
abline(h=x1,col="RED")
mtext(paste("CP:",signif(mean(abs(sC-x1)),digits = 2),
            "P:",signif(mean(abs(kP-x1)),digits=2),
            "PL",signif(mean(abs(kL-x1)),digits=2)),
      line=3,side=1)

## highly correlated ############################################

rho <- 0.8

for (i in 1:p){
  z<-mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),2,2))
  x<-z[,1]
  y<-z[,2]
  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y)
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names=c("cmis", "knnmi","cmik"))
x1<--log(1-rho^2)/2
abline(h=x1,col="RED")
mtext(paste("CP:",signif(mean(abs(sC-x1)),digits = 2),
            "P:",signif(mean(abs(kP-x1)),digits=2),
            "PL",signif(mean(abs(kL-x1)),digits=2)),
      line=3,side=1)

## very highly correlated ###################################

rho <- 0.95

for (i in 1:p){
  z<-mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),2,2))
  x<-z[,1]
  y<-z[,2]
  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y)
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names=c("cmis", "knnmi","cmik"))
x1<--log(1-rho^2)/2
abline(h=x1,col="RED")
mtext(paste("CP:",signif(mean(abs(sC-x1)),digits = 2),
            "P:",signif(mean(abs(kP-x1)),digits=2),
            "PL",signif(mean(abs(kL-x1)),digits=2)),
      line=3,side=1)


## independent exponential ####################################

sC<-rep(NaN,p)
kP<-rep(NaN,p)
kL<-rep(NaN,p)

for (i in 1:p){
  x<-rexp(N)
  y <- rexp(N)
  #  x<-rank(x)
  #  y<-rank(y)
  #   x<-log(x)
  #   y<-log(y)
  #   x<-qnorm(ecdf(x)(x))
  #   y<-qnorm(ecdf(y)(y))
  #   x<-x[x!=Inf]
  #   y<-y[y!=Inf]

  #   skew.score <- function(c, x) (skewness(log(x + c)))^2
  #   best.c <- optimise(skew.score, c(0, 20), x = x)$minimum
  #   x <- log(x + best.c)
  #   best.c <- optimise(skew.score, c(0, 20), x = y)$minimum
  #   y <- log(y + best.c)
  #

  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y)
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names = c("cmis", "knnmi","cmik"))
x1<-0
abline(h=x1,col="RED")
mtext(paste("CP:",signif(mean(abs(sC-x1)),digits = 2),
            "P:",signif(mean(abs(kP-x1)),digits=2),
            "PL",signif(mean(abs(kL-x1)),digits=2)),
      line=3,side=1)


## exponential dependent #################################

sC<-rep(NaN,p)
kP<-rep(NaN,p)
kL<-rep(NaN,p)

for (i in 1:p){
  x<-rexp(N)
  y <- x+rexp(N)
#  y<-x
#  x<-rank(x)
#  y<-rank(y)
#   x<-log(x)
#   y<-log(y)
#   x<-qnorm(ecdf(x)(x))
#   y<-qnorm(ecdf(y)(y))
#   x<-x[x!=Inf]
#   y<-y[y!=Inf]
#   skew.score <- function(c, x) (skewness(log(x + c)))^2
#   best.c <- optimise(skew.score, c(0, 20), x = x)$minimum
#   x <- log(x + best.c)
#   best.c <- optimise(skew.score, c(0, 20), x = y)$minimum
#   y <- log(y + best.c)

#   lambda<-BoxCox.lambda(x)
#   x<-BoxCox(x,lambda)
#   lambda<-BoxCox.lambda(y)
#   y<-BoxCox(y,lambda)

  sC[i] <- cmis(x,y)
  kP[i]<-knnmi(x,y)
  kL[i]<-cmik(x,y)
}

boxplot(sC,kP,kL,names = c("cmis", "knnmi","cmik"))
x1= 1-digamma(2)
abline(h=x1,col="RED")
mtext(paste("CP:",signif(mean(abs(sC-x1)),digits = 2),
            "P:",signif(mean(abs(kP-x1)),digits=2),
            "PL",signif(mean(abs(kL-x1)),digits=2)),
      line=3,side=1)

## spotting patterns ####################################
set.seed(17)
par(mfrow=c(4,2))

x<-state.x77[,1]
y<-state.x77[,2]
plot(x,y)
barplot(c(knnmi(x,y),cmik(sqrt(x),sqrt(y)),cmi.pw(x,y)$bcmi,cor(x,y)),
        names=c("knnmi","cmik","cmi","cor"))

x<-USArrests[,1]
y<-USArrests[,2]
plot(x,y)
barplot(c(knnmi(x,y),cmik(x,y),cmi.pw(x,y)$bcmi,cor(x,y)),
        names=c("knnmi","cmik","cmi","cor"))



X <-c(rnorm(50,-10,1), rnorm(50,0,1),rnorm(50,10,1))
Y <-c(rnorm(50,0,1), rnorm(50,0,1),rnorm(50,0,1))

R <- matrix(c(cos(pi/4),-sin(pi/4),sin(pi/4),cos(pi/4)),2,2)
RR <- R%*%rbind(as.numeric(scale(X)),as.numeric(scale(Y)))
x<-RR[1,]
y<-RR[2,]
plot(x,y)
barplot(c(knnmi(x,y),cmik(x,y),cmi.pw(x,y)$bcmi,cor(x,y)),
        names=c("knnmi","cmik","cmi","cor"))


x <-c(rnorm(50,-10,3), rnorm(50,0,3),rnorm(50,2,3))
y <-c(rnorm(50,-5,3), rnorm(50,5,3),rnorm(50,-5,3))
plot(x,y)
barplot(c(knnmi(x,y),cmik(x,y),cmi.pw(x,y)$bcmi,cor(x,y)),
        names=c("knnmi","cmik","cmi","cor"),ylim=c(-0.5,0.7))



abline(lm(y~x))
