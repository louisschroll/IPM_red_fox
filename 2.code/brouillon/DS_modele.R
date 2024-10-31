


# Recreate line transect data set
set.seed(1234)
tmp <- simdata          # Line transect (default)

data <- simdata$data
nsites <- simdata$nsites
dist_max <- simdata$dist_max
# Data augmentation: add a bunch of "pseudo-individuals"
nz <- 200                        # Augment by 500
nind <- nrow(data)
y <- c(data[, 2], rep(0, nz))     # Augmented detection indicator y
y[is.na(y)] <- 0
site <- c(data[, 1], rep(NA, nz)) # Augmented site indicator,
# unknown (i.e., NA) for augmented inds.
d <- c(data[, 3], rep(NA, nz))     # Augmented distance data (with NAs)

# Bundle and summarize data set
str(
  win.data <- list(
    nsites = nsites,
    dist_max = dist_max,
    nind = nind,
    nz = nz,
    y = y,
    d = d,
    site = site
  )
)
win.data$site                    # unknown site cov. for augmented inds.


# Jags model for line transect HDS (NOT point transects!)
cat(
  "
model{
  # Prior distributions
  for (s in 1:nsites){
    sigma[s] ~ dunif(0, 1000)
    lambda[s] ~ dunif(0, 1000)
  }

  # psi is a derived parameter under DA for stratified populations
  psi <- sum(lambda[1:nsites]) / (nind + nz)

  # 'Likelihood' 
  for(i in 1:(nind + nz)){                                    # i is index for individuals
    z[i] ~ dbern(psi)                                         # Data augmentation variables
    d[i] ~ dunif(0, dist_max)                                 # distance uniformly distributed
    p[i] <- exp(-d[i]*d[i]/(2*sigma[site[i]]*sigma[site[i]])) # Detection function
    y[i] ~ dbern(z[i] * p[i])                                 # basic Bernoulli random variable
    site[i] ~ dcat(site.probs[1:nsites])                      # Population distribution among sites
  }

  # Linear models for abundance and for detection
  for(s in 1:nsites){                                         # s is index for sites
    # Model for abundance
    # next line not necessary, but allows to make predictions
    N[s] ~ dpois(lambda[s])                                   # Realized abundance at site s
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])

    # Linear model for detection
    # log(sigma[s]) <- alpha0 + alpha1*wind[s]
  }
  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[1:(nind + nz)])
  area <- nsites*1*2*dist_max   # Unit length == 1, half-width = dist_max
  D <- Ntotal / area
}
", fill=TRUE, file = "model1.txt")


# Inits
zst <- c(rep(1, sum(y)), rep(0, nz)) # ... and for DA variables
inits <- function() {
  list(
    z = zst,
    sigma = rep(100, nsites),
    lambda = rep(20, nsites)
  )
}

# Parameters to save
params <- c("psi", "Ntotal", "D")

# MCMC settings
ni <- 12000
nb <- 2000
nt <- 2
nc <- 1


library(jagsUI)       # never forget to load jagsUI
out1 <- jags(
  win.data,
  inits,
  params,
  "model1.txt",
  # n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)
  n.thin = nt,
  n.chains = nc,
  n.burnin = nb,
  n.iter = ni,
  parallel = TRUE
)  # ~~~~ speeds up testing

# Summarize posterior output
print(out1, 2)
sum(tmp$N.true)
