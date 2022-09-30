limPM <- function(q) {

  ## M00
  sn.gr <- seq(from=.5,to=1,by=.005)
  tmp <- rep(NA, length(sn.gr))
  for (i in 1:length(sn.gr)) {
    p <- bwd.D(q, sn.gr[i])
    if (!is.null(p)) {
      tmp[i] <- (1-p$c)*(p$y..x.c[2,1]-p$y..x.c[1,1]) +
                    p$c*(p$y..x.c[2,2]-p$y..x.c[1,2])
    }
  }
  ans.00 <- mean(tmp[!is.na(tmp)])
  
  ## M10 
  ans.10 <- (1-q$c)*(q$y..xstr.c[2,1]-q$y..xstr.c[1,1]) +
                    q$c*(q$y..xstr.c[2,2]-q$y..xstr.c[1,2])
  
  ## M01
  p <- bwd.B(q)
  if (!is.null(p)) { 
    ans.01 <- p$del
  } else {
    ans.01 <- NA
  }  

list(M00=ans.00, M10=ans.10, M01=ans.01)  
}

