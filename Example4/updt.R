### function to do the MCMC under the cell prob parameterization
### keep what's needed to do all 4 inferences plus weighting

updt <- function(dat=NULL, NSIM=5000, sn.lo=.5) {

  tmp <- rep(0,8)
  if (!is.null(dat)) {
    tmp <- tmp +
           as.vector(table(factor(dat$y,levels=0:1),
                           factor(dat$x.str,levels=0:1),
                           factor(dat$c,levels=0:1)))
  }

  mc.opt <- MCmultinomdirichlet(tmp, alpha0=rep(1,8), mc=NSIM)

  delA.opt <- cmptB.opt <- delB.opt <- snB.opt <- dstC.opt <- cmptD.opt <- delD.opt <-
  rep(NA, NSIM)

  snD.opt <- runif(NSIM, sn.lo, 1)

  for (i in 1:NSIM) {
    tmp <- array(mc.opt[i,], dim=c(2,2,2))
    q <- list(
      c = sum(tmp[,,2]),
      xstr..c = apply(tmp[,2,], 2, sum) / apply(tmp, 3, sum),
      y..xstr.c = tmp[2,,] / apply(tmp, c(2,3), sum)
    )

    ## compute target under model A
    delA.opt[i] <- (1-q$c)*(q$y..xstr.c[2,1]-q$y..xstr.c[1,1]) +
                       q$c*(q$y..xstr.c[2,2]-q$y..xstr.c[1,2])

    ### compatible with model B?
    p <- bwd.B(q, sn.lo=sn.lo)
    cmptB.opt[i] <- !is.null(p)
    if (cmptB.opt[i]) {
      delB.opt[i] <- p$del
       snB.opt[i] <- p$sn
    }

    ## compute model C terms
    dstC.opt[i] <- ( (q$y..xstr.c[2,1]-q$y..xstr.c[1,1]) -
                     (q$y..xstr.c[2,2]-q$y..xstr.c[1,2]) )^2

    ## compute model D terms
    p <- bwd.D(q, snD.opt[i])
    cmptD.opt[i] <- !is.null(p)
    if (cmptD.opt[i]) {
      delD.opt[i] <- (1-p$c)*(p$y..x.c[2,1]-p$y..x.c[1,1]) +
                       p$c*(p$y..x.c[2,2]-p$y..x.c[1,2])
    }
  }

list(mc=mc.opt, delA=delA.opt,
     cmptB=cmptB.opt, snB=snB.opt, delB=delB.opt,
     dstC=dstC.opt,
     cmptD=cmptD.opt, delD=delD.opt)
}
