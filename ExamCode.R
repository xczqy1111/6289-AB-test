

# study 2
require('carat')
# imbalance
ib <- function(X.1, X.0, TA) {
  Ts <- ifelse(TA == 'A', 1, -1)
  oi <- abs(sum(Ts))
  p <- ncol(X.1)
  mis <- c()
  for (i in 1:p) {
    mis <- rbind(mis, matrix(abs(tapply(Ts, X.1[, i], sum))))
  }
  q <- ncol(X.0)
  wsis <- c()
  for (j in 1:q) {
    wsis <- rbind(wsis, abs(sum(Ts[which(X.0[, j] == 1)])))
  }
  list(oi, mis, wsis)
} 
set.seed(666)
nsim <- 1000
ns <- c(30, 60, 120)
np <- c()
Dn <- c()
Dnk <- c()
DnkI <- c()
for (n in ns) {
  np.n <- c()
  Dn.n <- c()
  Dnk.n <- c()
  DnkI.n <- c()
  for (i in 1:nsim) {
    X.0 <- t(rmultinom(n = n, size = 1, 
                     prob = c(3 / 20, 10 / 20, 
                              5 / 20, 2 / 20)))
    np.n <- rbind(np.n, apply(X.0, 2, sum))
    X1 <- X.0[, 1] + X.0[, 2]
    X2 <- X.0[, 1] + X.0[, 3]
    X.1 <- cbind(X1, X2)
    colnames(X.1) <- c('medical', 'smoke')
    # complete randomization
    res.cr <- rbinom(n, 1, prob = .5)
    res.cr <- ifelse(res.cr == 1, 'A', 'B')
    # Pocock and Simon's procedure
    res.poc <- PocSimMIN(data = X.1, weight = c(.5, .5),
                         p = .8)
    # PBR
    res.pbr <- StrPBR(data = X.1, bsize = 4)
    # Hu and Hu
    res.HH <- HuHuCAR(data = X.1, 
                      omega = c(.2, .2, .3, .3), 
                      p = .8)
    
    # imbalance
    im.cr <- ib(X.1 = X.1, X.0 = X.0, TA = res.cr)
    im.poc <- ib(X.1 = X.1, X.0 = X.0, TA = res.poc$assignments)
    im.pbr <- ib(X.1 = X.1, X.0 = X.0, TA = res.pbr$assignments)
    im.HH <- ib(X.1 = X.1, X.0 = X.0, TA = res.HH$assignments)
    Dn.n <- rbind(Dn.n, c(im.cr[[1]], im.poc[[1]],
                          im.pbr[[1]], im.HH[[1]]))
    Dnk.n <- rbind(Dnk.n, c(im.cr[[2]], im.poc[[2]], 
                            im.pbr[[2]], im.HH[[2]]))
    DnkI.n <- rbind(DnkI.n, c(im.cr[[3]], im.poc[[3]], 
                             im.pbr[[3]], im.HH[[3]]))
  }
  np <- rbind(np, apply(np.n, 2, mean))
  Dn <- rbind(Dn, apply(Dn.n, 2, 
                        function(x) 
                          c(mean(x), 
                            quantile(x, c(.5, .95)))))
  Dnk <- rbind(Dnk, matrix(apply(Dnk.n, 2, mean),
                           ncol = 4, byrow = FALSE))
  tab.ki <- matrix(apply(DnkI.n, 2, mean), 
                   ncol = 4, byrow = FALSE)
  tab.ki <- rbind(tab.ki, apply(tab.ki, 2, mean))
  DnkI <- rbind(DnkI, tab.ki)
}

# type I error
mu1 <- mu2 <- 0
beta1 <- 2; beta2 <- -1
res <- c()
for (n in ns) {
  res.n <- vector('list', length = 4)
  for (i in 1:nsim) {
    # covariates
    X.0 <- t(rmultinom(n = n, size = 1, 
                       prob = c(3 / 20, 10 / 20, 
                                5 / 20, 2 / 20)))
    np.n <- rbind(np.n, apply(X.0, 2, sum))
    X1 <- X.0[, 1] + X.0[, 2]
    X2 <- X.0[, 1] + X.0[, 3]
    X.1 <- cbind(X1, X2)
    colnames(X.1) <- c('medical', 'smoke')
    # response
    Ys <- c()
    for (j in 1:n) {
      Ys[j] <- beta1 * X.1[j, 1] + beta2 * X.1[j, 2] + rnorm(1)
    }
    # complete randomization
    res.cr <- rbinom(n, 1, prob = .5)
    # Pocock and Simon's procedure
    res.poc <- PocSimMIN(data = X.1, weight = c(.5, .5),
                         p = .8)
    res.poc <- ifelse(res.poc$assignments == 'A', 
                      1, 0)
    # PBR
    res.pbr <- StrPBR(data = X.1, bsize = 4)
    res.pbr <- ifelse(res.pbr$assignments == 'A', 
                      1, 0)
    # Hu and Hu
    res.HH <- HuHuCAR(data = X.1, 
                      omega = c(.2, .2, .3, .3), 
                      p = .8)
    res.HH <- ifelse(res.HH$assignments == 'A', 
                     1, 0)
    
    TAs <- list(res.cr, res.poc, res.pbr, res.HH)
    for (z in 1:4) {
      ta <- TAs[[z]]
      lmod1 <- lm(Ys ~ ta + X.1[, 1])
      lmod2 <- lm(Ys ~ ta + X.1[, 2])
      lmod3 <- lm(Ys ~ ta + X.1[, 1] + X.1[, 2])
      res.n[[z]] <- rbind(res.n[[z]],
                          c(ifelse(t.test(Ys[which(ta == 1)], 
                                    Ys[which(ta == 0)])$p.value <= .05,
                             yes = 1, no = 0),
                            ifelse(coef(summary(lmod1))[2, 4] <= .05, 
                             1, 0),
                            ifelse(coef(summary(lmod2))[2, 4] <= .05, 
                             1, 0),
                            ifelse(coef(summary(lmod3))[2, 4] <= .05, 
                             1, 0)))
    }
  }
  for (z in 1:4) {
    res.n[[z]] <- apply(res.n[[z]], 2, mean)
  }
  res <- c(res, res.n)
}
res <- do.call('rbind', res)
res <- res[c(1, 5, 9, 2, 6, 10, 3, 7, 11,
             4, 8, 12), ]
xtable(res, digits = 3)

# power
mu1s <- seq(from = 0, to = 1, by = .1)
beta1 <- 2; beta2 <- -1
res <- vector('list', length = 4)
for (n in ns) {
  res.mu1 <- vector('list', length = 4)
  for (mu1 in mu1s) {
    res.n <- vector('list', length = 4)
    for (i in 1:nsim) {
      # covariates
      X.0 <- t(rmultinom(n = n, size = 1, 
                         prob = c(3 / 20, 10 / 20, 
                                  5 / 20, 2 / 20)))
      X1 <- X.0[, 1] + X.0[, 2]
      X2 <- X.0[, 1] + X.0[, 3]
      X.1 <- cbind(X1, X2)
      colnames(X.1) <- c('medical', 'smoke')
      # complete randomization
      res.cr <- rbinom(n, 1, prob = .5)
      # Pocock and Simon's procedure
      res.poc <- PocSimMIN(data = X.1, weight = c(.5, .5),
                           p = .8)
      res.poc <- ifelse(res.poc$assignments == 'A', 
                        1, 0)
      # PBR
      res.pbr <- StrPBR(data = X.1, bsize = 4)
      res.pbr <- ifelse(res.pbr$assignments == 'A', 
                        1, 0)
      # Hu and Hu
      res.HH <- HuHuCAR(data = X.1, 
                        omega = c(.2, .2, .3, .3), 
                        p = .8)
      res.HH <- ifelse(res.HH$assignments == 'A', 
                       1, 0)
      
      TAs <- list(res.cr, res.poc, res.pbr, res.HH)
      for (z in 1:4) {
        ta <- TAs[[z]]
        # response
        Ys.z <- c()
        for (j in 1:n) {
          Ys.z[j] <- mu1 * ta[j] + 
            beta1 * X.1[j, 1] + beta2 * X.1[j, 2] + rnorm(1)
        }
        lmod1 <- lm(Ys.z ~ ta + X.1[, 1])
        lmod2 <- lm(Ys.z ~ ta + X.1[, 2])
        lmod3 <- lm(Ys.z ~ ta + X.1[, 1] + X.1[, 2])
        res.n[[z]] <- rbind(res.n[[z]],
                            c(ifelse(t.test(Ys.z[which(ta == 1)], 
                                            Ys.z[which(ta == 0)])$p.value <= .05,
                                     yes = 1, no = 0),
                              ifelse(coef(summary(lmod1))[2, 4] <= .05, 
                                     1, 0),
                              ifelse(coef(summary(lmod2))[2, 4] <= .05, 
                                     1, 0),
                              ifelse(coef(summary(lmod3))[2, 4] <= .05, 
                                     1, 0)))
      }
    }
    # compute power for each procedure
    for (z in 1:4) {
      res.n[[z]] <- apply(res.n[[z]], 2, mean)
    }
    # row bind different mu1
    for (z in 1:4) {
      res.mu1[[z]] <- rbind(res.mu1[[z]], res.n[[z]])
    }
  }
  # column bind different n
  for (z in 1:4) {
    res[[z]] <- cbind(res[[z]], res.mu1[[z]])
  }
}
tabs <- c()
for (j in 1:4) {
  tab.i <- c()
  for (i in 1:3) {
    for (z in 1:4) {
      tab.i <- cbind(tab.i, res[[z]][, j + (i - 1) * 4])
    }
  }
  tabs[[j]] <- tab.i
}
# t-test
xtable(tabs[[1]], digits = 3)
# lm(Z1)
xtable(tabs[[2]], digits = 3)
# lm(Z2)
xtable(tabs[[3]], digits = 3)
# lm(Z1, Z2)
xtable(tabs[[4]], digits = 3)

