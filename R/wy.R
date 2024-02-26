wy <- function(y, x, d, wx = 0.1, wy_seq = seq(0.1, 1, by = .1), wh = 1.5,
               B = 500, xdensity = "normal", method = "FM") {
  space <- "pdf"
  if (wy_seq[1] == 0) {
    stop("Error!: h Sequence should not start zero")
  }
  if (method == "FM" || method == "CM") {
    hy <- wy_seq
    dj <- matrix(0, nrow = B, ncol = 1)
    dist.r <- matrix(0, nrow = length(hy), ncol = 1)
    y <- as.matrix(y)
    # create progress bar
    pb <- txtProgressBar(min = 0, max = length(hy), style = 3)

    p <- ncol(x)
    n <- nrow(x)
    for (j in 1:length(hy)) {
      Hy <- hy[j]
      xy.dr <- itdr(y, x, d, wx = wx, wy = Hy, wh = wh, space = "pdf", xdensity = xdensity, method = method)
      s_d <- xy.dr$eta_hat
      dataframe <- data.frame(y, x)
      boost.df <- dataframe[sample(nrow(x), n, replace = TRUE), ]
      for (jj in 1:B) {
        boost.df <- dataframe[sample(nrow(x), n, replace = TRUE), ]
        y.boostrap <- boost.df[, 1]
        x.boostrap <- boost.df[, -c(1)]
        xy.dr <- itdr(y.boostrap, x.boostrap, d, wx = wx, wy = Hy, wh = wh, space = "pdf", xdensity = xdensity, method = method)
        s_dj <- xy.dr$eta_hat
        dist.dj <- dsp(s_dj, s_d)
        dj[jj] <- dist.dj$r
      }
      dist.r[j, 1] <- mean(dj)
      setTxtProgressBar(pb, j)
    }
    close(pb)
    disttab <- data.frame(hy = wy_seq, dbar = dist.r)
    hy.hat <- wy_seq[which.min(dist.r)]
  } else {
    stop("Error!:  method should be either 'FM' or 'CM'")
  }
  list(dis_wy = disttab, wy.hat = hy.hat)
}
