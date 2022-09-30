bwd.B <- function(q,tol=.005, sn.lo=0.5) {

  p <- list()

  gr.sn <- seq(from=sn.lo, to=1, length=10000)


  ### use the c=0 term to express del as a function of sn
  cvl <- 1
  tmp.a <- q$y..xstr.c[1,cvl]*(1-q$xstr..c[cvl])
  tmp.b <- 1-q$xstr..c[cvl]/gr.sn
  tmp.c <- ((1-gr.sn)/gr.sn)*q$xstr..c[cvl]*q$y..xstr.c[2,cvl]
  gr.del <- q$y..xstr.c[2,cvl] + (tmp.c-tmp.a)/tmp.b

  ### now solve to make the c=1 term correct
  cvl <- 2
  tmp.a <- q$y..xstr.c[1,cvl]*(1-q$xstr..c[cvl])
  tmp.b <- 1-q$xstr..c[cvl]/gr.sn
  tmp.c <- ((1-gr.sn)/gr.sn)*q$xstr..c[cvl]*q$y..xstr.c[2,cvl]

  objctv <- tmp.a - tmp.b*(q$y..xstr.c[2,cvl] - gr.del) - tmp.c
  tmp <- which.min(abs(objctv))

  p$sn <- gr.sn[tmp]
  p$del <- gr.del[tmp]

  ### fill in remaining components
  p$c <- q$c
  p$y..x0.c <- q$y..xstr.c[2,] - rep(p$del,2)
  p$x..c <- q$xstr..c / p$sn

  ### consider legit soln if
  if (!( (tmp>1)&&
         (tmp<length(gr.sn))&&
         ((objctv[tmp-1]*objctv[tmp+1])<0)&&
         (abs(p$del)<1)&&
         (min(c(p$y..x0.c,p$y..x0+p$del))>0)&&
         (max(c(p$y..x0.c,p$y..x0+p$del))<1)&&
         (max(p$x..c)<1) )) {
    p <- NULL
  }
p
}
