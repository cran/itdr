wx <- function(
    y, x, d, wx_seq = seq(0.1, 5, by = .1), wy = 1, wh = 1.5, B = 500, space = "mean",
    xdensity = "normal", method = "FM") {
  if (wx_seq[1] == 0) {
    stop("Error!: hx/sw2 Sequence should not start with zero")
  }
  if (method == "FM" || method == "CM") {
    hx <- wx_seq
    dj <- matrix(0, nrow = B, ncol = 1)
    dist.r <- matrix(0, nrow = length(hx), ncol = 1)
    y <- as.matrix(y)
    # create progress bar
    pb <- txtProgressBar(min = 0, max = length(hx), style = 3)

    p <- ncol(x)
    n <- nrow(x)
    for (j in 1:length(hx)) {
      Hx <- hx[j]
      xy.dr <- itdr(
        y, x, d,
        wx = Hx, wy = wy, wh = wh, space = space,
        xdensity = xdensity, method = method
      )

      s_d <- xy.dr$eta_hat
      dataframe <- data.frame(y, x)
      boost.df <- dataframe[sample(nrow(x), n, replace = TRUE), ]
      for (jj in 1:B) {
        boost.df <- dataframe[sample(nrow(x), n, replace = TRUE), ]
        y.boostrap <- boost.df[, 1]
        x.boostrap <- boost.df[, -c(1)]
        xy.dr <- itdr(y.boostrap, x.boostrap, d, wx = Hx, wy = wy, wh = wh, space = space, xdensity = xdensity, method = method)
        s_dj <- xy.dr$eta_hat
        dist.dj <- dsp(s_dj, s_d)
        dj[jj] <- dist.dj$r
      }
      dist.r[j, 1] <- mean(dj)
      setTxtProgressBar(pb, j)
    }
    close(pb)
    disttab <- data.frame(hx = wx_seq, dbar = dist.r)
    hx.hat <- wx_seq[which.min(dist.r)]
  } else {
    stop("Error!:  method should be either 'FM' or 'CM'")
  }
  list(dis_wx = disttab, wx.hat = hx.hat)
}
