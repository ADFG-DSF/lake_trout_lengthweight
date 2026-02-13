library(jagsUI)           # for running JAGS models
library(jagshelper)       # for plotting JAGS output
library(tidyverse)        # for data manipulation
library(sf)               # for spatial projections & transformations
library(ggrepel)          # for plotting


### ------- loading data and custom functions -------- ###

## loading input data for the JAGS model
load(file="Rdata/laketrout_sampling_formodel.Rdata")

## loading the JAGS output itself (this is a big file)
load(file="interim_posts/int_Winf_modelrun.Rdata")

## displaying how many MCMC iterations and how many chains
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

### more custom functions...
## defining a tweaked version of jagshelper::caterpillar for use in future plots
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

## a wrapper function that draws caterpillar plus lake names
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


# A good-enough sideways version of caterpillar (data is on x-axis)
side_cat <- function(x, a, xlab="") { # a will be an index to highlight

  ci <- c(0.5, 0.95)

  loq <- apply(x, 2, quantile, p = (1 - ci)/2, na.rm = T)
  hiq <- apply(x, 2, quantile, p = 1 - (1 - ci)/2, na.rm = T)
  meds <- apply(x, 2, median, na.rm=T)

  plot(NA, xlim=range(0, loq, hiq, na.rm=TRUE), ylim=c(ncol(x),1),
       ylab="", xlab=xlab, yaxt="n")

  segments(x0=loq[1,], x1=hiq[1,], y0=1:ncol(x),
           lend=1, lwd=3, col=adjustcolor(1, alpha.f=.2))
  segments(x0=loq[2,], x1=hiq[2,], y0=1:ncol(x),
           lend=1, lwd=1, col=adjustcolor(1, alpha.f=.2))

  segments(x0=loq[1,a], x1=hiq[1,a], y0=a,
           lend=1, lwd=6, col=1)
  segments(x0=loq[2,a], x1=hiq[2,a], y0=a,
           lend=1, lwd=2, col=1)
  segments(x0=meds[a], x1=meds[a], y0=a-1, y1=a+1,
           lend=1, lwd=2, col=1)

  # grid(nx=NULL, ny=NA)
}

# a quick function to add individual sideways caterpillar bars, given a single (MC)MC vector
catbars <- function(x, h, col=1) {
  ci <- c(0.5, 0.95)
  lolo <- quantile(x, p = (1 - ci)/2, na.rm = T)
  hihi <- quantile(x, p = 1 - (1 - ci)/2, na.rm = T)
  segments(x0=lolo[1], x1=hihi[1], y0=h, lwd=3, lend=1, col=col)
  segments(x0=lolo[2], x1=hihi[2], y0=h, lwd=1, lend=1, col=col)
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
crossplot(int_Winf_jags_out, p=c("Linf","Winf"), col = "random")
crossplot(int_Winf_jags_out, p=c("Linf","Winf"), drawx=TRUE, drawcross=FALSE, col = "random")
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
caterpillar1(int_Winf_jags_out, "t0", do_xax=FALSE,
            col=estcols)
caterpillar1(int_Winf_jags_out, "k", do_xax=FALSE,
            col=estcols)




### Mapping parameter estimates spatially - are there spatial stories to tell?

## creating an Alaska state outline
AK <- map_data("world") %>%
  filter(region=="USA") %>%
  filter(subregion=="Alaska") %>%
  cbind(sf::sf_project(pts=.[,1:2], to="+proj=aea +lat_1=55 +lat_2=65
    +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
    +ellps=GRS80")) %>%
  rename(x="1", y="2")

## weird little helper function
ff <- function(x, n=nrow(laketrout_Winf)) { # f for fill
  c(x, rep(NA, n - length(x)))
}

## creating a data frame for plotting, and pre-loading parameter estimates
morph1 <- cbind(morphometry, sf::sf_project(pts=morphometry[,4:5], to="+proj=aea +lat_1=55 +lat_2=65
    +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
    +ellps=GRS80")) %>%
  rename(x="1", y="2") %>%
  mutate(b0 = ff(int_Winf_jags_out$q50$b0),
         b1 = ff(int_Winf_jags_out$q50$b1),
         k = ff(int_Winf_jags_out$q50$k),
         t0 = ff(int_Winf_jags_out$q50$t0),
         Winf = ff(int_Winf_jags_out$q50$Winf),
         Linf = ff(int_Winf_jags_out$q50$Linf))

## creating a base map with the AK state outline
basemap <- ggplot() +
  geom_polygon(data=AK, mapping=aes(x=x, y=y, group=group), fill="white", col="black") +
  # geom_point(data=morph1,
  #            mapping=aes(x=x, y=y)) +
  coord_fixed() +
  theme_bw() +
  # theme(axis.title = element_blank(), axis.ticks = element_blank()) +
  scale_x_continuous(labels = NULL, breaks = NULL) + labs(x = "") +
  scale_y_continuous(labels = NULL, breaks = NULL) + labs(y = "") +
  scale_color_gradientn(colors = c("blue", "grey95", "red"))

## a wrapper function to make it easier to map one parameter at a time
ggAK <- function(theparm, digits=3, NAvec=NULL) {
  toplot <- select(morph1, all_of(theparm))[,1]
  minmax <- rep(NA, length(toplot))
  minmax[which.min(toplot)] <- paste(lakenames[which.min(toplot)],
                                     round(min(toplot, na.rm = TRUE), digits=digits))
  minmax[which.max(toplot)] <- paste(lakenames[which.max(toplot)],
                                     round(max(toplot, na.rm = TRUE), digits=digits))
  morphplot <- mutate(morph1, toplot = toplot,
                      minmax = minmax)
  if(!is.null(NAvec)) {
    morphplot$toplot[NAvec] <- NA
  }
  justmm <- filter(morphplot, !is.na(minmax))

  basemap +
    geom_point(data=justmm,
               mapping=aes(x=x, y=y), fill="white", pch=21, size=2.5) +
    # geom_point(data=morphplot,
    #            mapping=aes(x=x, y=y, col=toplot)) +
    geom_point(data=filter(morphplot, !is.na(toplot)),
               mapping=aes(x=x, y=y, col=toplot)) +
    geom_point(data=filter(morphplot, is.na(toplot)),
               mapping=aes(x=x, y=y), pch="+") +
    geom_text_repel(data=morphplot,
                    mapping=aes(x=x, y=y, label=minmax),
                    size=2) +
    labs(color=theparm)
}

wl_lakes <- subset(laketrout, !is.na(Weight_g) & !is.na(ForkLength_mm))$LakeNum %>%
  unique %>%
  sort %>%
  as.character %>%
  as.numeric
not_wl_lakes <- !((1:nrow(morphometry)) %in% wl_lakes)

## Actually plotting!
par(mfrow=c(1,1))
ggAK("b0")
ggAK("b0", NAvec = not_wl_lakes)
ggAK("b1", digits=2)
ggAK("b1", digits=2, NAvec = not_wl_lakes)
ggAK("k") +
  scale_color_gradientn(colors = c("blue", "grey95", "red"), trans="log10")
ggAK("t0")
ggAK("Linf", digits=0)
ggAK("Winf", digits=2) +
  scale_color_gradientn(colors = c("blue", "grey95", "red"), trans="log10")




## This is retained in case we would like to revisit interactive mapping with Leaflet

# library(leaflet)
# thecol <- ifelse(is.na(ff(int_Winf_jags_out$q50$b1)), "white",
#                  ifelse(ff(int_Winf_jags_out$q50$b1) < 3.15, "blue",
#                         ifelse(ff(int_Winf_jags_out$q50$b1) > 3.3, "red", "grey")))
# leaflet(morphometry) %>%
#   addTiles() %>%
#   addCircles(lng=~Longitude_WGS84,
#              lat=~Latitude_WGS84,
#              # color=colorRampPalette(c("blue","white","red"))(nrow(morphometry))[rank(ff(int_Winf_jags_out$q50$b1))]
#              color=thecol
#              )



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

# graphics.off()


## overlaying length ~ age on a single plot (MAYBE NOT USEFUL LIKE THIS)
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



## caterpillar plots of modeled mean length at 10 YEARS and 20 YEARS age
par(mfrow=c(1,2))

length10 <- int_Winf_jags_out$sims.list$Lfit[,which(int_Winf_data$Agefit==10),]
length20 <- int_Winf_jags_out$sims.list$Lfit[,which(int_Winf_data$Agefit==20),]

caterpillar1(length10, main="Mean length at age 10", ylim=c(0,1000), do_xax = FALSE)
text(x=1:ncol(length10), y=apply(length10, 2, median),
     labels=lakenames[1:ncol(length10)],
     cex=.6, col=adjustcolor(4, alpha.f=.6, blue.f=.8, red.f=.8, green.f=.8), srt=90)

caterpillar1(length20, main="Mean length at age 20", ylim=c(0,1000), do_xax = FALSE)
text(x=1:ncol(length20), y=apply(length20, 2, median),
     labels=lakenames[1:ncol(length20)],
     cex=.6, col=adjustcolor(4, alpha.f=.6, blue.f=.8, red.f=.8, green.f=.8), srt=90)



#######  weight ~ length (plots for all lakes where there are weights and lengths)
par(mfrow=c(2,2)) # maybe make this more appropriate for a bigger screen
for(i in seq_along(lakenames)) {
  # W ~ Length
  plot(NA, xlab="Length (mm)", ylab="Weight (kg)", main=lakenames[i], #log="xy",
       xlim=range(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE),
       ylim=range(exp(int_Winf_data$logW), na.rm=TRUE),
       # xlim=quantile(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE, p=c(.01,.999)),
       # ylim=quantile(exp(int_Winf_data$logW), na.rm=TRUE, p=c(.01,.999)),
  )
  for(j in 1:length(lakenames)) {
    curve(exp(int_Winf_jags_out$q50$b0_interp[j])*x^int_Winf_jags_out$q50$b1[j],
          add=TRUE, col=adjustcolor(1, alpha.f = .1))
  }
  points(x=int_Winf_data$L[int_Winf_data$lake==i],
         y=exp(int_Winf_data$logW)[int_Winf_data$lake==i])
  lvec <- seq(from=min(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE),
              to=max(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE),
              length.out=20)
  envelope(exp(int_Winf_jags_out$sims.list$b0_interp[,i]) *
             t(outer(lvec, int_Winf_jags_out$sims.list$b1[,i], "^")),
           x=lvec, add=TRUE)

  plot(NA, xlab="Length (mm)", ylab="Weight (kg)", main=lakenames[i], log="xy",
       xlim=range(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE),
       ylim=range(exp(int_Winf_data$logW), na.rm=TRUE),
       # xlim=quantile(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE, p=c(.01,.999)),
       # ylim=quantile(exp(int_Winf_data$logW), na.rm=TRUE, p=c(.01,.999)),
  )
  for(j in 1:length(lakenames)) {
    curve(exp(int_Winf_jags_out$q50$b0_interp[j])*x^int_Winf_jags_out$q50$b1[j],
          add=TRUE, col=adjustcolor(1, alpha.f = .1))
  }
  points(x=int_Winf_data$L[int_Winf_data$lake==i],
         y=exp(int_Winf_data$logW)[int_Winf_data$lake==i])
  lvec <- seq(from=min(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE),
              to=max(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE),
              length.out=20)
  envelope(exp(int_Winf_jags_out$sims.list$b0_interp[,i]) *
             t(outer(lvec, int_Winf_jags_out$sims.list$b1[,i], "^")),
           x=lvec, add=TRUE)
}

# graphics.off()



## caterpillar plots of modeled median weight at 300mm and 600mm??? <-- perhaps betterer lengths
par(mfrow=c(1,1))

makelabels <- function(x) {
  text(x=1:ncol(x), y=apply(x, 2, median),
       labels=lakenames[1:ncol(x)],
       cex=.6, col=adjustcolor(4, alpha.f=.6, blue.f=.8, red.f=.8, green.f=.8), srt=90)
}
maxes <- tapply(laketrout$ForkLength_mm, factor(laketrout$LakeName, levels=lakenames), max, na.rm=TRUE)

df300 <- exp(int_Winf_jags_out$sims.list$b0_interp)*
  300^int_Winf_jags_out$sims.list$b1
caterpillar1(df300, main="Median weight (kg) at 300mm", do_xax = FALSE,
            col=ifelse(is.na(maxes), adjustcolor(4, alpha.f=.2),
                       ifelse(maxes < 300, adjustcolor(4, alpha.f=.6), 4)))
makelabels(df300)

df600 <- exp(int_Winf_jags_out$sims.list$b0_interp)*
  600^int_Winf_jags_out$sims.list$b1
caterpillar1(df600, main="Median weight (kg) at 600mm", do_xax = FALSE,
            col=ifelse(is.na(maxes), adjustcolor(4, alpha.f=.2),
                       ifelse(maxes < 600, adjustcolor(4, alpha.f=.6), 4)))
makelabels(df600)

df800 <- exp(int_Winf_jags_out$sims.list$b0_interp)*
  800^int_Winf_jags_out$sims.list$b1
caterpillar1(df800, main="Median weight (kg) at 800mm", do_xax = FALSE,
            col=ifelse(is.na(maxes), adjustcolor(4, alpha.f=.2),
                       ifelse(maxes < 800, adjustcolor(4, alpha.f=.6), 4)))
makelabels(df800)


## mapping modeled median weight
morph1 <- morph1 %>%
  mutate(W300 = ifelse(is.na(maxes), NA,
                       ifelse(maxes < 300, NA,
                              apply(df300, 2, median))),
         W600 = ifelse(is.na(maxes), NA,
                       ifelse(maxes < 600, NA,
                              apply(df600, 2, median))),
         W800 = ifelse(is.na(maxes), NA,
                       ifelse(maxes < 800, NA,
                              apply(df800, 2, median))),
         x = morph1$x,
         y = morph1$y)


ggAK("W300")
ggAK("W600")
ggAK("W800")


## plotting by latitude
caterpillar1(df300[,!is.na(morph1$W300)], main="Median weight (kg) at 300mm",
             x=morphometry$Latitude_WGS84[!is.na(morph1$W300)],
             xlab="Latitude", do_xax = FALSE)
axis(side=1, axTicks(1))
caterpillar1(df600[,!is.na(morph1$W600)], main="Median weight (kg) at 600mm",
             x=morphometry$Latitude_WGS84[!is.na(morph1$W600)],
             xlab="Latitude", do_xax = FALSE)
axis(side=1, axTicks(1))
caterpillar1(df800[,!is.na(morph1$W800)], main="Median weight (kg) at 800mm",
             x=morphometry$Latitude_WGS84[!is.na(morph1$W800)],
             xlab="Latitude", do_xax = FALSE)
axis(side=1, axTicks(1))


# ## plotting by Log Area
# caterpillar1(df300[,!is.na(morph1$W300)], main="Median weight (kg) at 300mm",
#              x=log(morphometry$SurfaceArea_h)[!is.na(morph1$W300)],
#              xlab="Log Area", do_xax = FALSE)
# axis(side=1, axTicks(1))
# caterpillar1(df600[,!is.na(morph1$W600)], main="Median weight (kg) at 600mm",
#              x=log(morphometry$SurfaceArea_h)[!is.na(morph1$W600)],
#              xlab="Log Area", do_xax = FALSE)
# axis(side=1, axTicks(1))
# caterpillar1(df800[,!is.na(morph1$W800)], main="Median weight (kg) at 800mm",
#              x=log(morphometry$SurfaceArea_h)[!is.na(morph1$W800)],
#              xlab="Log Area", do_xax = FALSE)
# axis(side=1, axTicks(1))
#
#
# ## plotting by temperature
# caterpillar1(df300[,!is.na(morph1$W300)], main="Median weight (kg) at 300mm",
#              x=morphometry$`Temp (C)`[!is.na(morph1$W300)],
#              xlab="Temp (C)", do_xax = FALSE)
# axis(side=1, axTicks(1))
# caterpillar1(df600[,!is.na(morph1$W600)], main="Median weight (kg) at 600mm",
#              x=morphometry$`Temp (C)`[!is.na(morph1$W600)],
#              xlab="Temp (C)", do_xax = FALSE)
# axis(side=1, axTicks(1))
# caterpillar1(df800[,!is.na(morph1$W800)], main="Median weight (kg) at 800mm",
#              x=morphometry$`Temp (C)`[!is.na(morph1$W800)],
#              xlab="Temp (C)", do_xax = FALSE)
# axis(side=1, axTicks(1))






##### ------ Plotting asymptotic length  Linf ----- #####

par(mfrow=c(1,1))

# vs data availability
caterpillar_plus(p="Linf", x=datarank,
                 df=int_Winf_jags_out, col=datacols)

caterpillar1(p="Linf", x=datarank, do_xax=FALSE,
            df=int_Winf_jags_out, col=datacols)
axis(side=1, axTicks(1))
points(x=datarank[morphometry$make_estimates],
       y=int_Winf_jags_out$q50$Linf[morphometry$make_estimates],
       pch=16)
legend("topleft", pch=16, legend="Estimates of Winf desired", cex=0.5)


# vs qL value
caterpillar_plus(p="Linf", x=int_Winf_data$qL, xlab="qL",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)

caterpillar1(p="Linf", x=int_Winf_data$qL, do_xax=FALSE,
             df=int_Winf_jags_out, col=datacols, xlab="qL")
axis(side=1, axTicks(1))
points(x=int_Winf_data$qL[morphometry$make_estimates],
       y=int_Winf_jags_out$q50$Linf[morphometry$make_estimates],
       pch=16)
abline(0, 1, lty=3)
legend("topleft", pch=16, legend="Estimates of Winf desired", cex=0.5)


# vs Area
plot(x=morphometry$SurfaceArea_h, y=int_Winf_jags_out$q50$Linf,
     pch=16, col=datacols, log="x", ylim=c(0, max(int_Winf_jags_out$q97.5$Linf)),
     xlab="Surface Area (HA)", ylab="Linf (mm)")
caterpillar_plus(p="Linf", x=morphometry$SurfaceArea_h, #x=int_Winf_data$logareac,
                 df=int_Winf_jags_out, col=datacols,
                 add=TRUE, median=FALSE)
curve(int_Winf_jags_out$q50$gam * (1 - exp(-int_Winf_jags_out$q50$lam * (1 + log(x)))), add=TRUE, lty=2)
curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))), add=TRUE, lty=3)
legend("topright", lty=c(2,3), legend=c("Post median","Lester"), cex=.5)

plot(x=morphometry$SurfaceArea_h, y=int_Winf_jags_out$q50$Linf,
     pch=16, col=datacols, log="x", ylim=c(0, max(int_Winf_jags_out$q97.5$Linf)),
     xlab="Surface Area (HA)", ylab="Linf (mm)")
caterpillar1(p="Linf", x=morphometry$SurfaceArea_h, #x=int_Winf_data$logareac,
                 df=int_Winf_jags_out, col=datacols,
                 add=TRUE, median=FALSE)
points(x=morphometry$SurfaceArea_h[morphometry$make_estimates],
       y=int_Winf_jags_out$q50$Linf[morphometry$make_estimates],
       pch=16)
curve(int_Winf_jags_out$q50$gam * (1 - exp(-int_Winf_jags_out$q50$lam * (1 + log(x)))), add=TRUE, lty=2)
curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))), add=TRUE, lty=3)
legend("topright", lty=c(2,3), legend=c("Post median","Lester"), cex=.5)
legend("topleft", pch=16, legend="Estimates of Winf desired", cex=0.5)

# # vs Lester Linf
# caterpillar_plus(p="Linf", x=morphometry$L_inf_lester, xlab="Lester Linf",
#                  df=int_Winf_jags_out, col=datacols)
# abline(0, 1, lty=3)




##### ------ Plotting asymptotic weight Winf ----- #####

# vs data availability
caterpillar_plus(p="Winf", x=datarank,
                 df=int_Winf_jags_out, col=datacols)

caterpillar1(p="Winf", x=datarank, do_xax=FALSE,
             df=int_Winf_jags_out, col=datacols)
axis(side=1, axTicks(1))
points(x=datarank[morphometry$make_estimates],
       y=int_Winf_jags_out$q50$Winf[morphometry$make_estimates],
       pch=16)
legend("topleft", pch=16, legend="Estimates of Winf desired", cex=0.5)


# vs qW value
caterpillar_plus(p="Winf", x=int_Winf_data$qW, xlab="qW",
                 df=int_Winf_jags_out, col=datacols)
abline(0, 1, lty=3)

caterpillar1(p="Winf", x=int_Winf_data$qW, do_xax=FALSE,
             df=int_Winf_jags_out, col=datacols, xlab="qW")
axis(side=1, axTicks(1))
points(x=int_Winf_data$qW[morphometry$make_estimates],
       y=int_Winf_jags_out$q50$Winf[morphometry$make_estimates],
       pch=16)
abline(0, 1, lty=3)
legend("topleft", pch=16, legend="Estimates of Winf desired", cex=0.5)

# vs Area
plot(x=morphometry$SurfaceArea_h, y=int_Winf_jags_out$q50$Winf,
     pch=16, col=datacols, log="x", ylim=c(0, max(int_Winf_jags_out$q97.5$Winf, na.rm=TRUE)),
     xlab="Surface Area (HA)", ylab="Winf (kg)")
caterpillar_plus(p="Winf", x=morphometry$SurfaceArea_h, #x=int_Winf_data$logareac,
                 df=int_Winf_jags_out, col=datacols,
                 add=TRUE, median=FALSE)
curve(exp(-19.56 +
            3.2*log(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))))),
      add=TRUE, lty=3)

plot(x=morphometry$SurfaceArea_h, y=int_Winf_jags_out$q50$Winf,
     pch=16, col=datacols, log="x", ylim=c(0, max(int_Winf_jags_out$q97.5$Winf, na.rm=TRUE)),
     xlab="Surface Area (HA)", ylab="Winf (kg)")
caterpillar1(p="Winf", x=morphometry$SurfaceArea_h, #x=int_Winf_data$logareac,
                 df=int_Winf_jags_out, col=datacols,
                 add=TRUE, median=FALSE)
points(x=morphometry$SurfaceArea_h[morphometry$make_estimates],
       y=int_Winf_jags_out$q50$Winf[morphometry$make_estimates],
       pch=16)
curve(exp(-19.56 +
            3.2*log(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))))),
      add=TRUE, lty=3)
legend("topleft", pch=c(16,NA), lty=c(NA,3), legend=c("Estimates of Winf desired","Lester"), cex=0.5)

# # vs Lester Winf
# caterpillar_plus(p="Winf", x=morphometry$W_inf_lester, xlab="Lester Winf",
#                  df=int_Winf_jags_out, col=datacols)
# abline(0, 1, lty=3)



## Caterpillar plot of
## - k and t0
## - b0 and b1
## versus Area, Temp, Latitude, Linf
## NOT SURE IF THIS IS USEFUL
par(mfrow=c(1,2))
caterpillar_plus(p="k", x=NA,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="t0", x=NA,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="k", x=log(morphometry$SurfaceArea_h),
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Surface area (ha)")
caterpillar_plus(p="t0", x=log(morphometry$SurfaceArea_h),,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Surface area (ha)")

caterpillar_plus(p="k", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Mean temp (C)")
caterpillar_plus(p="t0", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Mean temp (C)")

caterpillar_plus(p="k", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Latitude")
caterpillar_plus(p="t0", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Latitude")

caterpillar_plus(p="k", x=int_Winf_jags_out$q50$Linf,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Linf (mm)")
caterpillar_plus(p="t0", x=int_Winf_jags_out$q50$Linf,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Linf (mm)")


caterpillar_plus(p="b0", x=NA,
                 df=int_Winf_jags_out, col=datacols)
caterpillar_plus(p="b1", x=NA,
                 df=int_Winf_jags_out, col=datacols)

caterpillar_plus(p="b0", x=log(morphometry$SurfaceArea_h),
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Surface area (ha)")
caterpillar_plus(p="b1", x=log(morphometry$SurfaceArea_h),,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Surface area (ha)")

caterpillar_plus(p="b0", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Mean temp (C)")
caterpillar_plus(p="b1", x=morphometry$`Temp (C)`,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Mean temp (C)")

caterpillar_plus(p="b0", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Latitude")
caterpillar_plus(p="b1", x=morphometry$Latitude_WGS84,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Latitude")

caterpillar_plus(p="b0", x=int_Winf_jags_out$q50$Linf,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Linf (mm)")
caterpillar_plus(p="b1", x=int_Winf_jags_out$q50$Linf,
                 df=int_Winf_jags_out, col=datacols,
                 xlab="Linf (mm)")


# ## Regressions of point estimates of parameters, versus site variables.
# ## This is handled more completely in 2a_LWmodel_exploration.R
# ## but is retained here for completeness.
#
# hasLW <- which(laketrout_Winf$n_Length > 0 & laketrout_Winf$n_Weight > 0)
#
# plot(int_Winf_jags_out$q50$b0[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW])
# summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW]))
# plot(int_Winf_jags_out$q50$b1[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW])
# summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW]))
#
# plot(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$`Temp (C)`[hasLW])
# summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$`Temp (C)`[hasLW]))
# plot(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$`Temp (C)`[hasLW])
# summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$`Temp (C)`[hasLW]))
#
# plot(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$Latitude_WGS84[hasLW])
# summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ morphometry$Latitude_WGS84[hasLW]))
# plot(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$Latitude_WGS84[hasLW])
# summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ morphometry$Latitude_WGS84[hasLW]))
#
# summary(lm(int_Winf_jags_out$q50$b0[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW] +
#              morphometry$`Temp (C)`[hasLW] +
#              morphometry$Latitude_WGS84[hasLW]))
#
# summary(lm(int_Winf_jags_out$q50$b1[hasLW] ~ log(morphometry$SurfaceArea_h)[hasLW] +
#              morphometry$`Temp (C)`[hasLW] +
#              morphometry$Latitude_WGS84[hasLW]))



## --- A COLLECTION OF SUMMARY PLOTS FOR EACH LAKE --- ##


## Actually plotting -- THIS WILL MAKE A LOT OF PLOTS, BE WARNED!
par(mfrow=c(2,2))
for(i in seq_along(lakenames)) {
  par(xpd=FALSE)

  # L ~ Age
  plot(NA, xlab="Age", ylab="Length (mm)", main=lakenames[i],
       xlim=c(0, max(int_Winf_data$Age, na.rm=TRUE)),
       ylim=c(0, max(int_Winf_data$L[!is.na(int_Winf_data$Age)], na.rm=TRUE)))
  for(j in 1:ncol(int_Winf_jags_out$q50$Lfit)) {
    lines(int_Winf_jags_out$q50$Lfit[,j], col=adjustcolor(1, alpha.f = .1))
  }
  if(i <= dim(int_Winf_jags_out$sims.list$Lfit)[3]) {
    envelope(int_Winf_jags_out$sims.list$Lfit[,,i], add=TRUE)
  }
  points(x=int_Winf_data$Age[int_Winf_data$lake==i],
         y=int_Winf_data$L[int_Winf_data$lake==i])


  # W ~ Length
  plot(NA, xlab="Length (mm)", ylab="Weight (kg)", main=lakenames[i], #log="xy",
       xlim=range(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE),
       ylim=range(exp(int_Winf_data$logW), na.rm=TRUE),
       # xlim=quantile(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE, p=c(.01,.999)),
       # ylim=quantile(exp(int_Winf_data$logW), na.rm=TRUE, p=c(.01,.999)),
       )
  for(j in 1:length(lakenames)) {
    curve(exp(int_Winf_jags_out$q50$b0_interp[j])*x^int_Winf_jags_out$q50$b1[j],
          add=TRUE, col=adjustcolor(1, alpha.f = .1))
  }
  points(x=int_Winf_data$L[int_Winf_data$lake==i],
         y=exp(int_Winf_data$logW)[int_Winf_data$lake==i])
  lvec <- seq(from=min(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE),
              to=max(int_Winf_data$L[!is.na(int_Winf_data$logW)], na.rm=TRUE),
              length.out=20)
  envelope(exp(int_Winf_jags_out$sims.list$b0_interp[,i]) *
                t(outer(lvec, int_Winf_jags_out$sims.list$b1[,i], "^")),
              x=lvec, add=TRUE)


  # Linf by data sources
  side_cat(x=int_Winf_jags_out$sims.list$Linf, a=i, xlab="Linf (mm)")

  par(xpd=NA)
  cols <- c(1, 2, 4) # ages, lengths, area

  segments(x0=int_Winf_data$qL[i], y0=-5-1, y1=1+i, lty=2, col=cols[2])
  sig_L <- int_Winf_jags_out$sims.list$sig_L
  if(i <= dim(int_Winf_jags_out$sims.list$Lfit)[3]) {
    qLvec <- rnorm(n=nrow(sig_L), mean=int_Winf_data$qL[i], sd=sig_L[,i])
    catbars(x=qLvec, h=-5, col=cols[2])
  }

  logmu_vec <- log(int_Winf_jags_out$sims.list$gam * (1 - exp(-int_Winf_jags_out$sims.list$lam * (1 + log(int_Winf_data$Area[i])))))
  sigLAvec <- int_Winf_jags_out$sims.list$sig_LA
  areaLvec <- rlnorm(n=length(logmu_vec), meanlog = logmu_vec, sdlog = sigLAvec)
  # abline(v=median(areaLvec), lty=2)
  segments(x0=median(areaLvec), y0=-3-1, y1=1+i, lty=2, col=cols[3])
  catbars(x=areaLvec, h=-3, col=cols[3])

  text(x=c(0, 0), y=c(-7, -5, -3)-1, pos=4, cex=0.5, col=cols,
       labels=paste(c("n ages =", "n lengths =", "area"),
                    c(laketrout_Winf$n_Age[i], laketrout_Winf$n_Length[i], "")))


  # Winf by data sources
  side_cat(x=int_Winf_jags_out$sims.list$Winf, a=i, xlab="Winf (kg)")

  par(xpd=NA)
  cols <- c(1, 2, 3, 4) # ages, lengths, area

  segments(x0=int_Winf_data$qW[i], y0=-7-1, y1=1+i, lty=2, col=cols[2])
  sig_W <- int_Winf_jags_out$sims.list$sig_W
  if(i <= dim(int_Winf_jags_out$sims.list$Lfit)[3]) {
    qWvec <- rnorm(n=nrow(sig_W), mean=int_Winf_data$qW[i], sd=sig_W[,i])
    catbars(x=qWvec, h=-7, col=cols[2])
  }

  qWLvec <- exp(int_Winf_jags_out$sims.list$b0_interp[,i])*qLvec^int_Winf_jags_out$sims.list$b1[,i]
  segments(x0=median(qWLvec), y0=-5-1, y1=1+i, lty=2, col=cols[3])
  catbars(x=qWLvec, h=-5, col=cols[3])

  areaWLvec <- exp(int_Winf_jags_out$sims.list$b0_interp[,i])*areaLvec^int_Winf_jags_out$sims.list$b1[,i]
  segments(x0=median(areaWLvec), y0=-3-1, y1=1+i, lty=2, col=cols[4])
  catbars(x=areaWLvec, h=-3, col=cols[4])

  text(x=c(0, 0), y=c(-9, -7, -5, -3)-1, pos=4, cex=0.5, col=cols,
       labels=paste(c("n ages =", "n weights =", "n lengths =", "area"),
                    c(laketrout_Winf$n_Age[i], laketrout_Winf$n_Weight[i], laketrout_Winf$n_Length[i], "")))
}

# graphics.off()


