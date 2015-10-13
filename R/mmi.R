mmi <- function(cts, disc, method="Kernel-Smooth",
                level=3L,na.rm=FALSE, h, bc,
                ki=floor(table(factor(disc))/2)){
  if (method=="Kernel-Smooth"){
    result <- mmis.pw(cts,disc)
    return(result)
  } else if (method=="k-Neighbor"){
    MI <- mmik(cts,disc,ki)
    return(MI)
  }
}