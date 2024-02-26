wh <- function(y, x, d, wx = 0.1, wy = 1, wh_seq = seq(0.1, 3, by = .1),
               B = 500, space = "mean", method = "FM") {
  if (wh_seq[1] == 0) {
    stop("Error!: h Sequence should not start zero")
  }
  if (method == "FM" || method == "CM") {
    xdensity <- "kernel"
    h <- wh_seq
    dj <- matrix(0, nrow = B, ncol = 1)
    dist.r <- matrix(0, nrow = length(h), ncol = 1)
    y <- as.matrix(y)
    # create progress bar
    pb <- txtProgressBar(min = 0, max = length(h), style = 3)

    p <- ncol(x)
    n <- nrow(x)
    for (j in 1:length(h)) {
      H <- h[j]
      xy.dr <- itdr(y, x, d, wx = wx, wy = wy, wh = H, space = space, xdensity = xdensity, method = method)
      s_d <- xy.dr$eta_hat
      dataframe <- data.frame(y, x)
      boost.df <- dataframe[sample(nrow(x), n, replace = TRUE), ]
      for (jj in 1:B) {
        boost.df <- dataframe[sample(nrow(x), n, replace = TRUE), ]
        y.boostrap <- boost.df[, 1]
        x.boostrap <- boost.df[, -c(1)]
        xy.dr <- itdr(y.boostrap, x.boostrap, d, wx = wx, wy = wy, wh = H, space = space, xdensity = xdensity, method)
        s_dj <- xy.dr$eta_hat
        dist.dj <- dsp(s_dj, s_d)
        dj[jj] <- dist.dj$r
      }
      dist.r[j, 1] <- mean(dj)
      setTxtProgressBar(pb, j)
    }
    close(pb)
    disttab <- data.frame(h = wh_seq, dbar = dist.r)
    h.hat <- wh_seq[which.min(dist.r)]
  } else {
    stop("Error!:  method should be either 'FM' or 'CM'")
  }
  list(dis_h = disttab, h.hat = h.hat)
}
