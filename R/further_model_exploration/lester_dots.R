library(jpeg)
# Read the JPG file
img <- readJPEG("raw_data/lester_linf_area.jpg", native = TRUE)

pars <- par(no.readonly = TRUE)

# Create an empty plot area
plot(0:1, 0:1, par(mar=c(0,0,0,0)))#, type = "n", ann = FALSE, axes = FALSE)

# Render the image inside the plot window
rasterImage(img, 0, 0, 1, 1)

# locator(n=2, type="p", pch="+")

corners <- locator(n=2, type="p", pch="+")
# should be x=log(c(5,500000)), y=c(1000,0)

dots <- locator(n=129, type="p", pch="+")


origx <- corners$x[1]
origy <- corners$y[2]
points(origx, origy)

rangex <- abs(diff(corners$x))
rangey <- abs(diff(corners$y))

scalex <- abs(diff(log(c(5,500000))))/rangex
scaley <- 1000/rangey

datay <- (dots$y-origy)*scaley
datax <- exp((dots$x-origx)*scalex + log(5))

par(pars)
plot(datax, datay, log="x")
abline(v=c(50,50000, 500000))
save(datax, datay, file="raw_data/dotsfromdots.Rdata")



load(file="Rdata/laketrout_sampling_formodel.Rdata")
load(file="interim_posts/int_Winf_modelrun.Rdata")
# load(file="interim_posts/int_Winf_modelrun_26_05_30.Rdata")
library(magrittr)
library(jagshelper)
# plot(datax, datay, log="x")

par(mar=c(3, 3, 1, 1)+0.1)
par(family="serif")
par(mgp=c(2, 1, 0))

# y_cat2(int_Winf_jags_out$sims.list$Linf,
#        xval=morphometry$SurfaceArea_h,
#        highlight=morphometry$make_estimates,
#        x0=TRUE,
#        log="x",
#        xlab="Surface area (ha)",
#        ylab="Asymptotic length (mm FL)")

datalakes <- rep(FALSE, length(int_Winf_jags_out$q50$k))
datalakes[union(union(int_Winf_data$whichlakes_L,
      int_Winf_data$whichlakes_W),
      int_Winf_data$whichlakes_Lt)] <- TRUE

caterpillar(int_Winf_jags_out, p="Linf",
            x = morphometry$SurfaceArea_h,
            # col = plotcol,
            col = ifelse(!datalakes, "grey", 1),
            log="x", main="", ylim_add=0,
            xlim=range(morphometry$SurfaceArea_h, datax, na.rm=TRUE),
            xlab="Surface area (ha)",
            # ylab="Asymptotic length (mm FL)",
            ylab=expression(paste(L[infinity], " (mm FL)")))
curve(int_Winf_jags_out$q50$gam * (1 - exp(-int_Winf_jags_out$q50$lam * (1 + log(x)))), add=TRUE, lty=2)
curve(int_Winf_data$gam_lester * (1 - exp(-int_Winf_data$lam_lester * (1 + log(x)))), add=TRUE, lty=3)

points(datax, datay)




# let's try re-fitting the Lester dots!!
# specify model, which is written to a temporary file
dots_jags <- tempfile()
cat('model {
  for(i in 1:n) {
    # y[i] ~ dlnorm(logmu[i], tau)
    # logmu[i] <- log(gam * (1 - exp(-lam * (1 + log(x[i])))))
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- gam * (1 - exp(-lam * (1 + log(x[i]))))
  }

  tau <- pow(sig, -2)
  # sig ~ dunif(0, 3)
  sig ~ dunif(0, 300)

  gam ~ dnorm(gam_lester, pow(cv_gam_lester*gam_lester, -2))T(0.01,)
  lam ~ dnorm(lam_lester, pow(cv_lam_lester*lam_lester, -2))T(0.01,)

}', file=dots_jags)

# bundle data to pass into JAGS
dots_data <- list(x=datax,
                  y=datay,
                  n=length(datax),
                  gam_lester=int_Winf_data$gam_lester,
                  lam_lester=int_Winf_data$lam_lester,
                  cv_gam_lester=5*int_Winf_data$cv_gam_lester,
                  cv_lam_lester=5*int_Winf_data$cv_lam_lester)

# JAGS controls
niter <- 100000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  dots_jags_out <- jagsUI::jags(model.file=dots_jags, data=dots_data,
                                parameters.to.save=c("sig","gam","lam"),
                                n.chains=ncores, parallel=T, n.iter=niter,
                                n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# nbyname(dots_jags_out)
plotRhats(dots_jags_out)
traceworstRhat(dots_jags_out, parmfrow = c(3, 3))

dots_jags_out$q50
dots_jags_out$sd
