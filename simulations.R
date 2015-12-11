library(migen)

## TMI calculates the target MI
#### set integration limit ########
w2 <- c(1,1,1)
probs <- w2 / sum(w2)

il<-20

TMI <- function(x){
  u <- function(y){
    z<-0
    for(i in 1:nrow(x)){
      z<-z+w2[i]/sum(w2)*dnorm(y,x[i,1],x[i,2])
    }
    return(z)
  }
  TMI <- 0
  for(i in 1:nrow(x)){
    ui <- function(y){
      z<-dnorm(y,x[i,1],x[i,2])*log(dnorm(y,x[i,1],x[i,2])/u(y))
      return(z)
    }
    TMI <- TMI+w2[i]/sum(w2)*integrate(ui,lower=-il,upper=il)$value
  }
  return(TMI)
}

################################################################

x1 <- TMI(cbind(c(0,0,0),rep(1,3)))
x2 <- TMI(cbind(c(-1,0,1),rep(1,3)))
x3 <- TMI(cbind(c(-5,0,5),rep(1,3)))

################################################################
set.seed(11)
p<- 10000
s<- 100

mip1 <- rep(NA,p)
mip2 <- rep(NA,p)

start1 <-Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-rnorm(sum(w),0,1)
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip1[i]<-mmis(X,Y)
}
end1 <-Sys.time()

start2 <- Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-rnorm(sum(w),0,1)
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip2[i]<-mmik(X,Y)
}
end2 <- Sys.time()

mip3 <- rep(NA,p)
mip4 <- rep(NA,p)

start3 <-Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-c(rnorm(w[1],-1,1), rnorm(w[2],0,1), rnorm(w[3],1,1))
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip3[i]<-mmis(X,Y)
}
end3 <-Sys.time()

start4 <- Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-c(rnorm(w[1],-1,1), rnorm(w[2],0,1), rnorm(w[3],1,1))
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip4[i]<-mmik(X,Y)
}
end4 <- Sys.time()

mip5 <- rep(NA,p)
mip6 <- rep(NA,p)

start5 <-Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-c(rnorm(w[1],-5,1), rnorm(w[2],0,1), rnorm(w[3],5,1))
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip5[i]<-mmis(X,Y)
}
end5 <-Sys.time()

start6 <- Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-c(rnorm(w[1],-5,1), rnorm(w[2],0,1), rnorm(w[3],5,1))
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip6[i]<-mmik(X,Y)
}
end6 <- Sys.time()

## plot ################################
par(mfrow=c(1,3))

boxplot(mip1,mip2,
        main="No Separation")
# mtext(paste("C:", signif(end1-start1,digits = 2),
#             "P:", signif(end2-start2,digits=2)),
#       line=2,side=1)
mtext(paste("s:", signif(mean(abs(mip1-x1)),digits = 2),
            "k:", signif(mean(abs(mip2-x1)),digits=2)),
      line=3,side=1)
abline(h=x1,col="RED")

boxplot(mip3,mip4,
        main="Difficult to Separate")
# mtext(paste("C:", signif(end3-start3,digits = 2),
#         "P:", signif(end4-start4,digits=2)),
#         line=2,side=1)
mtext(paste("s:", signif(mean(abs(mip3-x2)),digits = 2),
            "k:", signif(mean(abs(mip4-x2)),digits=2)),
      line=3,side=1)
abline(h=x2,col="RED")

boxplot(mip5,mip6,
        main="Easily Separated")
# mtext(paste("C:", signif(end5-start5,digits = 2),
#             "P:", signif(end6-start6,digits=2)),
#       line=2,side=1)
mtext(paste("s:", signif(mean(abs(mip5-x3)),digits = 2),
            "k:", signif(mean(abs(mip6-x3)),digits=2)),
      line=3,side=1)
abline(h=x3,col="RED")

# MSE:
# mean((x-target)^2)
# 


#### Gamma ##################################
w2 <- c(1,1,1)
probs <- w2 / sum(w2)

il<-20

TMI <- function(x){
  u <- function(y){
    z<-0
    for(i in 1:nrow(x)){
      z<-z+w2[i]/sum(w2)*dgamma(y,shape = x[i,1], scale = x[i,2])
    }
    return(z)
  }
  TMI <- 0
  for(i in 1:nrow(x)){
    ui <- function(y){
      z<-dgamma(y,shape=x[i,1],scale=x[i,2])*log(dgamma(y,shape=x[i,1],scale=x[i,2])/u(y))
      return(z)
    }
    TMI <- TMI+w2[i]/sum(w2)*integrate(ui,lower=0,upper=il)$value
  }
  return(TMI)
}

################################################################

means <- c(1,1,1)
vars <- rep(1,3)
scales1<-vars/means
shapes1 <- (means+vars)/(scales1*(1+scales1))
x1 <- TMI(cbind(shapes1,scales1))

means <- c(1,2,3)
vars <- rep(1,3)
scales2<-vars/means
shapes2 <- (means+vars)/(scales2*(1+scales2))
x2 <- TMI(cbind(shapes2,scales2))

means <- c(1,5,10)
vars <- c(3,2,1)
scales3<-vars/means
shapes3 <- (means+vars)/(scales3*(1+scales3))
x3 <- TMI(cbind(shapes3,scales3))

################################################################
set.seed(12)
p<- 10000
s<- 100

#X<-c(rep("A",w[1]*s), rep("B",w[2]*s), rep("C",w[3]*s))

mip1 <- rep(NA,p)
mip2 <- rep(NA,p)

start1 <-Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-rgamma(sum(w),shape = shapes1[1],scale=scales1[1])
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip1[i]<-mmis(X,Y)
}
end1 <-Sys.time()

start2 <- Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-rgamma(sum(w),shape = shapes1[1],scale=scales1[1])
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip2[i]<-mmik(X,Y)
}
end2 <- Sys.time()

mip3 <- rep(NA,p)
mip4 <- rep(NA,p)

start3 <-Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-c(rgamma(w[1],shape=shapes2[1],scale = scales2[1]), 
       rgamma(w[2],shape=shapes2[2],scale = scales2[2]),
       rgamma(w[3],shape=shapes2[3],scale = scales2[3]))
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip3[i]<-mmis(X,Y)
}
end3 <-Sys.time()

start4 <- Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-c(rgamma(w[1],shape=shapes2[1],scale = scales2[1]), 
       rgamma(w[2],shape=shapes2[2],scale = scales2[2]),
       rgamma(w[3],shape=shapes2[3],scale = scales2[3]))
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip4[i]<-mmik(X,Y)
}
end4 <- Sys.time()

mip5 <- rep(NA,p)
mip6 <- rep(NA,p)

start5 <-Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-c(rgamma(w[1],shape=shapes3[1],scale = scales3[1]), 
       rgamma(w[2],shape=shapes3[2],scale = scales3[2]),
       rgamma(w[3],shape=shapes3[3],scale = scales3[3]))
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip5[i]<-mmis(X,Y)
}
end5 <-Sys.time()

start6 <- Sys.time()
for(i in 1:p) {
  w <- rowSums(rmultinom(s, 1, prob=probs))
  X<-c(rgamma(w[1],shape=shapes3[1],scale = scales3[1]), 
       rgamma(w[2],shape=shapes3[2],scale = scales3[2]),
       rgamma(w[3],shape=shapes3[3],scale = scales3[3]))
  Y<-c(rep("A",w[1]), rep("B",w[2]), rep("C",w[3]))
  mip6[i]<-mmik(X,Y)
}
end6 <- Sys.time()

## plot ################################
par(mfrow=c(1,3))

boxplot(mip1,mip2,
        main="No Separation")
# mtext(paste("C:", signif(end1-start1,digits = 2),
#             "P:", signif(end2-start2,digits=2)),
#       line=2,side=1)
mtext(paste("s:", signif(mean(abs(mip1-x1)),digits = 2),
            "k:", signif(mean(abs(mip2-x1)),digits=2)),
      line=3,side=1)
abline(h=x1,col="RED")

boxplot(mip3,mip4,
        main="Difficult to Separate")
# mtext(paste("C:", signif(end3-start3,digits = 2),
#             "P:", signif(end4-start4,digits=2)),
#       line=2,side=1)
mtext(paste("s:", signif(mean(abs(mip3-x2)),digits = 2),
            "k:", signif(mean(abs(mip4-x2)),digits=2)),
      line=3,side=1)
abline(h=x2,col="RED")

boxplot(mip5,mip6,
        main="Easily Separated")
# mtext(paste("C:", signif(end5-start5,digits = 2),
#             "P:", signif(end6-start6,digits=2)),
#       line=2,side=1)
mtext(paste("s:", signif(mean(abs(mip5-x3)),digits = 2),
            "k:", signif(mean(abs(mip6-x3)),digits=2)),
      line=3,side=1)
abline(h=x3,col="RED")
