bwd.D <- function(q, sn) {

  p <- list(y..x.c=matrix(NA,2,2))

  if (max(q$xstr..c) < sn) {
    p$c <- q$c
    p$x..c <- q$xstr..c/sn
    p$y..x.c[2,] <- q$y..xstr.c[2,]

    tmp <- q$y..xstr.c[1,] - q$y..xstr.c[2,]*((1-sn)/sn)*q$xstr..c/(1-q$xstr..c)
    tmp <- tmp * (1-q$xstr..c)/(1-q$xstr..c/sn)

    if ((min(tmp)>0)&&(max(tmp)<1)) {
      p$y..x.c[1,] <- tmp
    } else {
      p <- NULL
    }
  } else {
    p <- NULL
  }
p
}
