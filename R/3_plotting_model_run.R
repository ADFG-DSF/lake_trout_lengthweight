library(jagsUI)
library(jagshelper)
library(tidyverse)


load(file="Rdata/laketrout_sampling_formodel.Rdata")

load(file="interim_posts/int_Winf_modelrun.Rdata")
niter/1000
ncores


### a few helper functions and color vectors for plotting
estcols <- ifelse(morphometry$make_estimates, 4, 3)
datacols <- ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_W, 4,
                   ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 2,
                          ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 3,
                                 ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 5, 1))))
datarank <- rank(ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_W, 1,
                        ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 2,
                               ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 3,
                                      ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 4, 5)))), ties.method = "first")
datalegend <- function(loc="topleft") {
  legend(loc, legend=c("Wq, then", "Ages, then", "Lq, then", "Area, then","nothing"),
         lwd=3, col=c(4,2,3,5,1), cex=.5)
}

### should try tabulating datarank stuff more fully
# datacols uses:
# - presence of Wq (whichlakes_W)
# - presence of Age relationship (whichlakes_Lt)
# - presence of Lq (whichlakes_L)
# - presence of Area (whichlakes_LA)
nn <- 1:nrow(morphometry)
with(int_Winf_data,
     table(nn %in% whichlakes_W,
           nn %in% whichlakes_Lt))
## ^^ not convinced this is useful



### ------------------- convergence diagnostics  ------------------ ###

# how many nodes (or dim) per each named parameter
nbyname(int_Winf_jags_out)

# how many named parameters
length(nbyname(int_Winf_jags_out))

# plotting the Rhat values of all parameters
par(mfrow=c(1,1))
plotRhats(int_Winf_jags_out)

# trace plots for the nodes with the worst Rhat value for each named parameter
traceworstRhat(int_Winf_jags_out, parmfrow=c(3,3))

graphics.off()



### ------------------- exploring model appropriateness  ------------------ ###

# comparing distributions of priors to the associated parameter
comparepriors(int_Winf_jags_out, parmfrow=c(2,2))

# plotting the correlation between all parameters
par(mfrow=c(1,1))
plotcor_jags(int_Winf_jags_out, p=c("Linf", "mu_LA", "t0", "k", "b0", "b1", "sig","gam","lam","eta","zeta"))

graphics.off()


# cross plots between parameters ------ FILL THIS IN BETTER

# length 59
# t0, k, sig_L, qLpp, sig_W, qWpp

# length 62
# Linf, Winf, b0, b1, mu_b0, mu_b1, b0_interp

# length 1
# sig_Lt
# mu_t0, sig_t0, mu_k, sig_k
# sig_LA
# gam, lam
# eta_L, zeta_L
# eta_W, zeta_W
# sig_WL
# sig_b0, sig_b1
# b0_int, b1_int

par(mfrow=c(1,1))
crossplot(int_Winf_jags_out, p=c("t0","k"), drawblob=TRUE, col = "random")
crossplot(int_Winf_jags_out, p=c("sig_L","sig_W"), drawblob=TRUE, col = "random")
crossplot(int_Winf_jags_out, p=c("sig_L","qLpp"), drawblob=TRUE, col = "random")
crossplot(int_Winf_jags_out, p=c("sig_W","qWpp"), drawblob=TRUE, col = "random")
crossplot(int_Winf_jags_out, p=c("Linf","Winf"), drawblob=TRUE, col = "random")
crossplot(int_Winf_jags_out, p=c("b0","b1"), drawx=TRUE, col = "random")
crossplot(int_Winf_jags_out, p=c("mu_b0","mu_b1"), drawx=TRUE, col = "random")

par(mfrow=c(2,2))
crossplot(int_Winf_jags_out, p=c("sig_Lt","sig_LA"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("sig_Lt","sig_WL"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("mu_t0","sig_t0"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("mu_k","sig_k"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("mu_t0","mu_k"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("sig_t0","sig_k"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("gam","lam"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("eta_L","zeta_L"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("eta_W","zeta_W"), drawblob=TRUE)
crossplot(int_Winf_jags_out, p=c("sig_b0","sig_b1"), drawblob=TRUE)

## ^^ keep filling these in


### Posterior predictive for qL
par(mfrow=c(1,2))
qq_postpred(int_Winf_jags_out$sims.list$qLpp[,int_Winf_data$whichlakes_L],
            y=int_Winf_data$qL[int_Winf_data$whichlakes_L],
            main="Post predictive for qL")
ts_postpred(int_Winf_jags_out$sims.list$qLpp[,int_Winf_data$whichlakes_L],
            y=int_Winf_data$qL[int_Winf_data$whichlakes_L],
            x=int_Winf_data$nL[int_Winf_data$whichlakes_L], log="x",
            xlab="n Length for quantiles")
par(mfcol=c(3,3))
plot_postpred(int_Winf_jags_out$sims.list$qLpp[,int_Winf_data$whichlakes_L],
              y=int_Winf_data$qL[int_Winf_data$whichlakes_L],
              x=log(int_Winf_data$nL[int_Winf_data$whichlakes_L]))


### plot how sig_L varies by lake
par(mfrow=c(1,2))
caterpillar(int_Winf_jags_out, "sig_L", col=estcols)
envelope(int_Winf_jags_out$sims.list$sig_L[, int_Winf_data$whichlakes_L],
         x=int_Winf_data$nL[int_Winf_data$whichlakes_L], log="x",
         xlab="n Length for quantiles", ylab="sig_L")



### Posterior predictive for qW
par(mfrow=c(1,2))
qq_postpred(int_Winf_jags_out$sims.list$qWpp[,int_Winf_data$whichlakes_W],
            y=int_Winf_data$qW[int_Winf_data$whichlakes_W],
            main="Post predictive for qW")
ts_postpred(int_Winf_jags_out$sims.list$qWpp[,int_Winf_data$whichlakes_W],
            y=int_Winf_data$qW[int_Winf_data$whichlakes_W],
            x=int_Winf_data$nW[int_Winf_data$whichlakes_W], log="x",
            xlab="n Wength for quantiles")
par(mfcol=c(3,3))
plot_postpred(int_Winf_jags_out$sims.list$qWpp[,int_Winf_data$whichlakes_W],
              y=int_Winf_data$qW[int_Winf_data$whichlakes_W],
              x=log(int_Winf_data$nW[int_Winf_data$whichlakes_W]))


### plot how sig_W varies by lake
par(mfrow=c(1,2))
caterpillar(int_Winf_jags_out, "sig_W", col=estcols)
envelope(int_Winf_jags_out$sims.list$sig_W[, int_Winf_data$whichlakes_W],
         x=int_Winf_data$nW[int_Winf_data$whichlakes_W], log="x",
         xlab="n Weight for quantiles", ylab="sig_W")





### ------------------- deep-dive into results  ------------------ ###


### LVB parameters
par(mfrow=c(1,2))
caterpillar(int_Winf_jags_out, "t0",
            col=estcols)
caterpillar(int_Winf_jags_out, "k",
            col=estcols)





## length ~ age (plots for all lakes where there are lengths and ages)
par(mfrow=c(3,3))
for(j in int_Winf_data$whichlakes_Lt) {
  envelope(int_Winf_jags_out$sims.list$Lfit[,,j], main=lakenames[j],
           ylim=range(0,int_Winf_data$L[!is.na(int_Winf_data$Age)]),
           # col=ifelse(j %in% int_Winf_data$whichlakes_L, 2, 4),
           col=datacols[j])
  points(int_Winf_data$Age[int_Winf_data$lake==j], int_Winf_data$L[int_Winf_data$lake==j])
  # legend("topleft", legend=c("L ~ Age only", "qL also"),
  #        fill=adjustcolor(c(4,2), alpha.f=.3), border=c(4,2), cex=.6)
  datalegend()
}

par(mfrow=c(1,1))
cols <- adjustcolor(rcolors(ncol(int_Winf_jags_out$q50$Lfit)), red.f = .5, blue.f=.5, green.f=.5)
plot(NA, xlim=c(0, 70), ylim=c(0,1000), xlab="Age", ylab="Length")
for(j in int_Winf_data$whichlakes_Lt) {
  envelope(int_Winf_jags_out$sims.list$Lfit[,,j],
           add=TRUE, col=adjustcolor(cols[j], alpha.f=.5))
}
for(j in int_Winf_data$whichlakes_Lt) {
  lines(int_Winf_jags_out$q50$Lfit[,j],
        col=cols[j], lwd=2)
}
# text(x=rep(max(int_Winf_data$Agefit), length(int_Winf_data$whichlakes_Lt)),
#      y=int_Winf_jags_out$q50$Lfit[nrow(int_Winf_jags_out$q50$Lfit),],
#      labels=lakenames[int_Winf_data$whichlakes_Lt], pos=4)
text(x=rep(max(int_Winf_data$Agefit), ncol(int_Winf_jags_out$q50$Lfit)),
     y=int_Winf_jags_out$q50$Lfit[nrow(int_Winf_jags_out$q50$Lfit),],
     labels=lakenames[1:ncol(int_Winf_jags_out$q50$Lfit)], pos=4, col=cols)



caterpillar1 <- function (df, p = NULL, x = NA, row = NULL, column = NULL, median = TRUE,
                          mean = FALSE, ci = c(0.5, 0.95), lwd = 1, col = 4, add = FALSE,
                          xlab = "", ylab = "", main = NULL, ylim = NULL, xax = NA, do_xax=TRUE,
                          transform = c("none", "exp", "expit"), medlwd = lwd, medwd = 1,
                          ...)
{
  if (!inherits(df, "jagsUI") & !inherits(df, c("matrix", "data.frame",
                                                "numeric", "integer"))) {
    stop("Input must be a data.frame or output from jagsUI::jags() plus parameter name")
  }
  if (inherits(df, "jagsUI") & length(p) != 1)
    stop("Need single parameter name in p= argument")
  if (class(df)[1] == "jagsUI" & !is.null(p)) {
    simslist <- df$sims.list
    if (all(names(simslist) != p))
      stop("No parameters with matching names")
    theparm <- simslist[names(simslist) == p][[1]]
    if (is.null(row) & is.null(column)) {
      df <- theparm
    }
    else {
      theparm <- simslist[names(simslist) == p][[1]]
      if (!is.null(row))
        df <- theparm[, row, ]
      if (!is.null(column))
        df <- theparm[, , column]
    }
    if (is.null(main) & length(p) == 1) {
      if (is.null(row) & is.null(column)) {
        main <- p
      }
      else {
        if (!is.null(row))
          main <- paste0(p, "[", row, ",]")
        if (!is.null(column))
          main <- paste0(p, "[,", column, "]")
      }
    }
  }
  if (is.null(main))
    main <- ""
  df <- as.matrix(df)
  transform <- match.arg(transform)
  if (transform == "exp")
    df <- exp(df)
  if (transform == "expit")
    df <- expit(df)
  ci <- rev(sort(ci))
  loq <- apply(df, 2, quantile, p = (1 - ci)/2, na.rm = T)
  hiq <- apply(df, 2, quantile, p = 1 - (1 - ci)/2, na.rm = T)
  if (length(ci) == 1) {
    loq <- t(as.matrix(loq))
    hiq <- t(as.matrix(hiq))
  }
  med <- apply(df, 2, median, na.rm = T)
  if (all(is.na(x)))
    x <- 1:ncol(df)
  d <- ifelse(length(x) > 1, diff(range(x, na.rm = TRUE))/length(x),
              1)
  nn <- ncol(df)
  if (all(is.na(xax)))
    xax <- names(df)
  lwds <- (1 + 2 * (1:length(ci) - 1)) * lwd
  if (!add) {
    if (is.null(ylim))
      ylim <- range(loq, hiq, na.rm = T)
    xlims <- range(x - (0.5 * d * (1 + (length(x) == 1))),
                   x + (0.5 * d * (1 + (length(x) == 1))), na.rm=TRUE)
    plot(NA, type = "l", xlim = xlims, xlab = xlab, ylab = ylab,
         main = main, ylim = ylim, xaxt = "n", ... = ...)
    if(do_xax) axis(1, x, labels = xax, las = list(...)$las)
  }
  if (median) {
    segments(x0 = x - 0.2 * d * medwd, x1 = x + 0.2 * d *
               medwd, y0 = med, y1 = med, col = col, lwd = medlwd,
             lend = 1)
  }
  if (mean)
    points(x, colMeans(df, na.rm = T), pch = 16, col = col)
  for (i in 1:length(ci)) segments(x0 = x, x1 = x, y0 = loq[i,
  ], y1 = hiq[i, ], col = col, lwd = lwds[i], lend = 1)
}   ##### ADD THIS PATCH TO JAGSHELPER

caterpillar_plus <- function(df, p, x=NA, col=datacols, ...) {
  nn <- length(df$q50[[p]])
  if(all(is.na(x))) x <- 1:nn
  caterpillar1(df=df, p=p, col=col, x=x, do_xax=FALSE, ...=...)
  axis(side=1, axTicks(1))
  text(x=x[1:nn], y=df$q50[[p]],
       labels=lakenames[1:nn],
       cex=.6, col=adjustcolor(col[1:nn], alpha.f=.6, blue.f=.8, red.f=.8, green.f=.8), srt=90)
  datalegend()
}

par(mfrow=c(1,2))

length10 <- int_Winf_jags_out$sims.list$Lfit[,which(int_Winf_data$Agefit==10),]
length20 <- int_Winf_jags_out$sims.list$Lfit[,which(int_Winf_data$Agefit==20),]

caterpillar(length10, main="Mean length at age 10", ylim=c(0,1000))
text(x=1:ncol(length10), y=apply(length10, 2, median),
     labels=lakenames[1:ncol(length10)],
     cex=.6, col=adjustcolor(4, alpha.f=.6, blue.f=.8, red.f=.8, green.f=.8), srt=90)

caterpillar(length20, main="Mean length at age 20", ylim=c(0,1000))
text(x=1:ncol(length20), y=apply(length20, 2, median),
     labels=lakenames[1:ncol(length20)],
     cex=.6, col=adjustcolor(4, alpha.f=.6, blue.f=.8, red.f=.8, green.f=.8), srt=90)




##### Linf #####

# vs data availability
caterpillar_plus(p="Linf", x=datarank,
                 df=int_Winf_jags_out, col=datacols)

# vs qL value
caterpillar_plus(p="Linf", x=int_Winf_data$qL, xlab="qL",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)

# vs Area
plot(x=morphometry$SurfaceArea_h, y=int_Winf_jags_out$q50$Linf,
     pch=16, col=datacols, log="x", ylim=c(0, max(int_Winf_jags_out$q97.5$Linf)),
     xlab="Surface Area (HA)", ylab="Linf (mm)")
caterpillar_plus(p="Linf", x=morphometry$SurfaceArea_h, #x=int_Winf_data$logareac,
                 df=int_Winf_jags_out, col=datacols,
                 add=TRUE, median=FALSE)
curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))), add=TRUE, lty=2)

# vs Lester Linf
caterpillar_plus(p="Linf", x=morphometry$L_inf_lester, xlab="Lester Linf",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)




##### Winf #####

# vs data availability
caterpillar_plus(p="Winf", x=datarank,
                 df=int_Winf_jags_out, col=datacols)

# vs qW value
caterpillar_plus(p="Winf", x=int_Winf_data$qW, xlab="qW",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)

# vs Area
plot(x=morphometry$SurfaceArea_h, y=int_Winf_jags_out$q50$Winf,
     pch=16, col=datacols, log="x", ylim=c(0, max(int_Winf_jags_out$q97.5$Winf, na.rm=TRUE)),
     xlab="Surface Area (HA)", ylab="Winf (kg)")
caterpillar_plus(p="Winf", x=morphometry$SurfaceArea_h, #x=int_Winf_data$logareac,
                 df=int_Winf_jags_out, col=datacols,
                 add=TRUE, median=FALSE)
curve(exp(-19.56 +
            3.2*log(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))))),
      add=TRUE, lty=2)

# vs Lester Linf
caterpillar_plus(p="Winf", x=morphometry$W_inf_lester, xlab="Lester Winf",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)



## possibly try this with
## - k and t0
## - b0 and b1
# versus Area, Elevation, Latitude, Temp, Linf
par(mfrow=c(1,2))
caterpillar_plus(p="k", x=NA,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="t0", x=NA,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="k", x=log(morphometry$SurfaceArea_h),
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="t0", x=log(morphometry$SurfaceArea_h),,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="k", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="t0", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="k", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="t0", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="b0", x=NA,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="b1", x=NA,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="b0", x=log(morphometry$SurfaceArea_h),
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="b1", x=log(morphometry$SurfaceArea_h),,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="b0", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="b1", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="b0", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="b1", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols)



hasLW <- which(laketrout_Winf$n_Length > 0 & laketrout_Winf$n_Weight > 0)

plot(int_Winf_jags_out$q50$b0[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW])
summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW]))
plot(int_Winf_jags_out$q50$b1[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW])
summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW]))

plot(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$`Temp (C)`[hasLW])
summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$`Temp (C)`[hasLW]))
plot(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$`Temp (C)`[hasLW])
summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$`Temp (C)`[hasLW]))

plot(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$Latitude_WGS84[hasLW])
summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$Latitude_WGS84[hasLW]))
plot(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$Latitude_WGS84[hasLW])
summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$Latitude_WGS84[hasLW]))

summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW] +
             morphometry$`Temp (C)`[hasLW] +
             morphometry$Latitude_WGS84[hasLW]))

summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW] +
             morphometry$`Temp (C)`[hasLW] +
             morphometry$Latitude_WGS84[hasLW]))

### ^^^^^ this stuff is handled much better elsewhere




## ------- figure out better what is happening here


### Linf by data availability
par(mfrow=c(2,2))
caterpillar(int_Winf_jags_out, "Linf", col=3+int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt)
legend("topleft", legend=c("has AGES", "no AGES"), col=4:3, lwd=3, cex=0.5)
caterpillar(int_Winf_jags_out, "Linf", col=3+int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA)
legend("topleft", legend=c("has AREA", "no AREA"), col=4:3, lwd=3, cex=0.5)
caterpillar(int_Winf_jags_out, "Linf", col=3+int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L)
legend("topleft", legend=c("has Lq", "no Lq"), col=4:3, lwd=3, cex=0.5)
caterpillar(int_Winf_jags_out, "Linf", col=3+int_Winf_data$alllakes %in% int_Winf_data$whichlakes_W)
legend("topleft", legend=c("has Wq", "no Wq"), col=4:3, lwd=3, cex=0.5)

# # ages and lq
# caterpillar(int_Winf_jags_out, "Linf", col=3+(int_Winf_data$alllakes %in%
#               intersect(int_Winf_data$whichlakes_L, int_Winf_data$whichlakes_Lt)))
# legend("topleft", legend=c("has Lq and AGES", "no"), col=4:3, lwd=3, cex=0.5)
#
# # ages and lq and area
# caterpillar(int_Winf_jags_out, "Linf", col=3+(int_Winf_data$alllakes %in%
#                                              intersect(int_Winf_data$whichlakes_L,
#                                                        intersect(int_Winf_data$whichlakes_Lt,
#                                                        int_Winf_data$whichlakes_LA))))
# legend("topleft", legend=c("has Lq and AGES and AREA", "no"), col=4:3, lwd=3, cex=0.5)

# ages, then lq, then area
par(mfrow=c(1,1))
theorder <- ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 1,
                   ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 2,
                          ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 3, 4))) %>% rank(ties.method = "first")
caterpillar(int_Winf_jags_out, "Linf", x=theorder,
            col=ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 4,
                       ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 2,
                              ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 3, 1))))
legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing"), col=c(4,2,3,1), lwd=3, cex=0.5)
points(x=theorder[int_Winf_data$whichlakes_W], y=int_Winf_jags_out$q50$Linf[int_Winf_data$whichlakes_W])#, col=int_Winf_data$alllakes %in% int_Winf_data$whichlakes_W)
points(x=theorder[morphometry$make_estimates], y=int_Winf_jags_out$q50$Linf[morphometry$make_estimates], pch=16)

points(x=theorder, y=morphometry$L_inf_lester, pch="x")

caterpillar(int_Winf_jags_out, "Linf", x=datarank, col=datacols)
datalegend()

## ---

plot(NA, xlim=range(int_Winf_data$qL, na.rm=T),
     ylim=range(int_Winf_jags_out$q2.5$Linf, int_Winf_jags_out$q97.5$Linf, na.rm=TRUE),
     main="Linf", xlab="qL", ylab="")
caterpillar(int_Winf_jags_out, "Linf", col=ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 4,
                                                  ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 2,
                                                         ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 3, 1))),
            x=int_Winf_data$qL,
            add=TRUE)
legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing"), col=c(4,2,3,1), lwd=3, cex=0.5)
abline(0,1, lty=2)
# points(x=int_Winf_data$qL[int_Winf_data$whichlakes_LA],
#        y=int_Winf_jags_out$q50$mu_LA[int_Winf_data$whichlakes_LA])


plot(NA, xlim=range(int_Winf_data$qL, na.rm=T),
     ylim=range(int_Winf_jags_out$q2.5$Linf, int_Winf_jags_out$q97.5$Linf, na.rm=TRUE),
     main="Linf", xlab="qL", ylab="")
caterpillar(int_Winf_jags_out, "Linf", x=int_Winf_data$qL, col=datacols, add=T)
abline(0,1, lty=2)
datalegend()

#---

plot(NA, xlim=range(int_Winf_data$qW, na.rm=T),
     ylim=range(int_Winf_jags_out$q2.5$Winf, int_Winf_jags_out$q97.5$Winf, na.rm=TRUE),
     main="Winf", xlab="qW", ylab="")
caterpillar(int_Winf_jags_out, "Winf", x=int_Winf_data$qW, col=datacols, add=T)
abline(0,1, lty=2)
datalegend()


# ----------

par(mfrow=c(1,1))
plot(NA,
     ylim=range(int_Winf_jags_out$q2.5$Linf, int_Winf_jags_out$q97.5$Linf, na.rm=TRUE),
     xlim=range(log(int_Winf_data$Area), -1, na.rm=TRUE),
     main="Linf", xlab="Area (ha)", ylab="", xaxt="n")
lbls <- c(1, 5, 10, 50, 100, 500, 1000, 5000)
axis(side=1, at=log(lbls), labels=lbls)
thecols <- ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_Lt, 4,
                  ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_L, 2,
                         ifelse(int_Winf_data$alllakes %in% int_Winf_data$whichlakes_LA, 3, 1)))
caterpillar(int_Winf_jags_out, "Linf", col=datacols,#col=thecols,
            x=ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)),
            add=TRUE)
points(x=ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)),
       y=int_Winf_jags_out$q50$Linf,
       col=datacols, pch=16)
curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + x))), add=TRUE, lty=2)
# legend("topleft", legend=c("AGES, then","Lq, then","AREA, then","nothing","Lester priors"),
#        col=c(4,2,3,1,1), lwd=c(rep(3,4),1), lty=c(rep(1,4),2), cex=0.5)
datalegend()

# ---

# this bit was commented out
points(x=(ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)))[int_Winf_data$whichlakes_L],
       y=int_Winf_data$qL[int_Winf_data$whichlakes_L])
points(x=(ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)))[morphometry$make_estimates],
       y=int_Winf_jags_out$q50$Linf[morphometry$make_estimates], pch=16)

points(x=ifelse(is.na(int_Winf_data$Area), -1, log(int_Winf_data$Area)), y=morphometry$L_inf_lester, pch="x")

# more of the same, maybe these thoughts can be combined?
caterpillar(int_Winf_jags_out$sims.list$Linf[,int_Winf_data$whichlakes_L],
            x=int_Winf_data$qL[int_Winf_data$whichlakes_L],
            xlab="qL", main="subset with qL")
abline(0,1,lty=3)
caterpillar(int_Winf_jags_out$sims.list$Linf[,int_Winf_data$whichlakes_Lt],
            x=int_Winf_data$qL[int_Winf_data$whichlakes_Lt],
            xlab="qL", main="subset with Ages")
abline(0,1,lty=3)
int <- intersect(int_Winf_data$whichlakes_Lt, int_Winf_data$whichlakes_L)
caterpillar(int_Winf_jags_out$sims.list$Linf[,int],
            x=int_Winf_data$qL[int],
            xlab="qL", main="subset with qL and Ages")
abline(0,1,lty=3)


## logW ~ logL for each lake
## b0 ~ lat & area
## b1 ~ lat & area
## mu_b0 ~ lat & area
## mu_b1 ~ lat & area



### ------------------- compiling and saving output  ------------------ ###


### filling in summary output to save to external file

ff <- function(x, n=nrow(laketrout_Winf)) { # f for fill
  c(x, rep(NA, n - length(x)))
}

laketrout_Winf$Winf_est <- int_Winf_jags_out$q50$Winf %>% ff
laketrout_Winf$Winf_se <- int_Winf_jags_out$sd$Winf %>% ff
laketrout_Winf$Winf_cilo <- int_Winf_jags_out$q2.5$Winf %>% ff
laketrout_Winf$Winf_cihi <- int_Winf_jags_out$q97.5$Winf %>% ff

laketrout_Winf$Linf_est <- int_Winf_jags_out$q50$Linf %>% ff
laketrout_Winf$Linf_se <- int_Winf_jags_out$sd$Linf %>% ff
laketrout_Winf$Linf_cilo <- int_Winf_jags_out$q2.5$Linf %>% ff
laketrout_Winf$Linf_cihi <- int_Winf_jags_out$q97.5$Linf %>% ff

laketrout_Winf$t0_est <- int_Winf_jags_out$q50$t0 %>% ff
laketrout_Winf$t0_se <- int_Winf_jags_out$sd$t0 %>% ff
laketrout_Winf$k_est <- int_Winf_jags_out$q50$k %>% ff
laketrout_Winf$k_se <- int_Winf_jags_out$sd$k %>% ff

laketrout_Winf$b0_est <- int_Winf_jags_out$q50$b0_interp %>% ff
laketrout_Winf$b0_se <- int_Winf_jags_out$sd$b0_interp %>% ff
laketrout_Winf$b1_est <- int_Winf_jags_out$q50$b1 %>% ff
laketrout_Winf$b1_se <- int_Winf_jags_out$sd$b1 %>% ff

save_results
if(save_results) {
  write.csv(laketrout_Winf, file="Winf_estimates.csv")
}

##### add convergence diagnostics etc here!!

## make sure data object is appropriate
## make sure model output is complete (plots etc?)
## - pull out plotting stuff into its own 3_ script
