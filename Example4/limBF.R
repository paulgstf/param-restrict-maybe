limBF <- function(q) {

  mc0 <- updt()

  lmbf.ba <- 0
  tmp <- bwd.B(q)
  if (!is.null(tmp)) {
    lmbf.ba <- 1/mean(mc0$cmptB)
  }

  sn.gr <- seq(from=.5,to=1,by=.005)
  tmp <- rep(NA, length(sn.gr))
  for (i in 1:length(sn.gr)) {
    tmp[i] <- !is.null(bwd.D(q, sn.gr[i]))
  }
  lmbf.da <- mean(tmp)/mean(mc0$cmptD)
  
  tmp <- c(lmbf.da,1, lmbf.ba)
  tmp <- tmp/sum(tmp)
  names(tmp) <- c("M00","M10","M01")
  
  list(ba=lmbf.ba, da=lmbf.da, pspr=tmp)
}
