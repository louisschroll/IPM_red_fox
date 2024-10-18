


sim_DS_data <- function (type = c("line", "point"),
                         nsites = 100,
                         mean.lambda = 200,
                         mean.sigma = 30,
                         B = 100,
                         discard0 = TRUE){
  type <- match.arg(type)
  sigma <- rnorm(n = nsites, mean = mean.sigma, sd = 1)
  N <- rpois(nsites, lambda)
  N.true <- N
  data <- NULL
  for (i in 1:nsites) {
    if (N[i] == 0) {
      data <- rbind(data, c(i, NA, NA, NA, NA))
      next
    }
    if (type == "line") {
      distances <- runif(N[i], 0, B)
      detection_prob <- exp(-distances * distances / (2 * (sigma[i] ^ 2)))
      y <- rbinom(N[i], 1, detection_prob)
      u <- v <- rep(NA, N[i])
      distances <- distances[y == 1]
      u <- u[y == 1]
      v <- v[y == 1]
      y <- y[y == 1]
    }
    if (type == "point") {
      u <- runif(N[i], 0, 2 * B)
      v <- runif(N[i], 0, 2 * B)
      distances <- sqrt((u - B) ^ 2 + (v - B) ^ 2)
      N.true[i] <- sum(distances <= B)
      detection_prob <- exp(-distances * distances / (2 * (sigma[i] ^ 2)))
      pp <- ifelse(distances <= B, 1, 0) * detection_prob
      y <- rbinom(N[i], 1, pp)
      u <- u[y == 1]
      v <- v[y == 1]
      distances <- distances[y == 1]
      y <- y[y == 1]
    }
    if (sum(y) > 0)
      data <- rbind(data, cbind(rep(i, sum(y)), y, u, v, distances))
    else
      data <- rbind(data, c(i, NA, NA, NA, NA))
  }
  colnames(data) <- c("site", "y", "u", "v", "d")
  if (discard0)
    data <- data[!is.na(data[, 2]), ]
  
  list(
    type = type,
    nsites = nsites,
    lambda = lambda,
    sigma = sigma,
    B = B,
    data = data,
    N = N,
    N.true = N.true
  )
}


plot_DS_data <- function(DS_data){
  with (DS_data, {
    if (type == "line") {
      op <- par(mfrow = c(1, 3))
      on.exit(par(op))
      tryPlot <- try({
        hist(
          data[, "d"],
          col = "lightblue",
          breaks = 20,
          main = "Frequency of distances",
          xlab = "Distance"
        )
        ttt <- table(data[, 1])
        n <- rep(0, nsites)
        n[as.numeric(rownames(ttt))] <- ttt
        plot(habitat, n, main = "Observed counts (n) vs. habitat")
        plot(wind, n, main = "Observed counts (n) vs. wind speed")
      }, silent = TRUE)
      if (inherits(tryPlot, "try-error"))
        tryPlotError(tryPlot)
    }
    if (type == "point") {
      op <- par(mfrow = c(2, 2))
      on.exit(par(op))
      tryPlot <- try({
        plot(
          data[, "u"],
          data[, "v"],
          pch = 16,
          main = "Located individuals in point transects",
          xlim = c(0, 2 * B),
          ylim = c(0, 2 * B),
          col = data[, 1],
          asp = 1
        )
        points(B,
               B,
               pch = "+",
               cex = 3,
               col = "black")
        draw.circle(B, B, B)
        hist(
          data[, "d"],
          col = "lightblue",
          breaks = 20,
          main = "Frequency of distances",
          xlab = "Distance"
        )
        ttt <- table(data[, 1])
        n <- rep(0, nsites)
        n[as.numeric(rownames(ttt))] <- ttt
        plot(habitat, n, main = "Observed counts (n) vs. habitat")
        plot(wind, n, main = "Observed counts (n) vs. wind speed")
      }, silent = TRUE)
      if (inherits(tryPlot, "try-error"))
        tryPlotError(tryPlot)
    }
  })
}

simdata <- sim_DS_data()  
plot_DS_data(simdata)

# Recreate line transect data set
set.seed(1234)
tmp <- simdata          # Line transect (default)
attach(tmp)
# Data augmentation: add a bunch of "pseudo-individuals"
nz <- 500                        # Augment by 500
nind <- nrow(data)
y <- c(data[,2], rep(0, nz))     # Augmented detection indicator y
site <- c(data[,1], rep(NA, nz)) # Augmented site indicator,
# unknown (i.e., NA) for augmented inds.
d <- c(data[,5], rep(NA,nz))     # Augmented distance data (with NAs)

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, habitat=habitat, wind=wind, B=B,
                      nind=nind, nz=nz, y=y, d=d, site=site) )
win.data$site                    # unknown site cov. for augmented inds.


# BUGS model for line transect HDS (NOT point transects!)
cat("
model{
  # Prior distributions
  beta0 ~ dunif(-10,10)   # Intercept of lambda-habitat regression
  beta1 ~ dunif(-10,10)   # Slope of log(lambda) on habitat
  alpha0 ~ dunif(-10,10)  # Intercept of log(sigma) (half-normal scale)
  alpha1 ~ dunif(-10,10)  # Slope of log(sigma) on wind

  # psi is a derived parameter under DA for stratified populations
  psi <- sum(lambda[]) / (nind+nz)

  # 'Likelihood' (sort of...)
  for(i in 1:(nind + nz)){                 # i is index for individuals
    z[i] ~ dbern(psi)                    # Data augmentation variables
    d[i] ~ dunif(0, B)                   # distance uniformly distributed
    p[i] <- exp(-d[i]*d[i]/(2*sigma[site[i]]*sigma[site[i]])) # Det. function
    mu[i] <- z[i]* p[i]                  # 'straw man' for WinBUGS
    y[i] ~ dbern(mu[i])                  # basic Bernoulli random variable
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
  }

  # Linear models for abundance and for detection
  for(s in 1:nsites){                    # s is index for sites
    # Model for abundance
    # next line not necessary, but allows to make predictions
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0 + beta1*habitat[s] # Linear model abundance
    site.probs[s] <- lambda[s] / sum(lambda[])

    # Linear model for detection
     log(sigma[s]) <- alpha0 + alpha1*wind[s]
  }
  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[])
  area <- nsites*1*2*B   # Unit length == 1, half-width = B
  D <- Ntotal/area
}
",fill=TRUE , file = "model1.txt")


# Inits
zst <- c(rep(1, sum(y)), rep(0, nz)) # ... and for DA variables
inits <- function(){list(beta0=0, beta1=0, alpha0=0, alpha1=0, z=zst)}

# Parameters to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3


library(jagsUI)       # never forget to load jagsUI
out1 <- jags(win.data, inits, params, "model1.txt",
             # n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)
             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel=TRUE)  # ~~~~ speeds up testing

# Summarize posterior output
print(out1, 2)
sum(tmp$N.true)
