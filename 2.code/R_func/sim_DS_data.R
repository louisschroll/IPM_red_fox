#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/IPM_red_fox/2.code/sim_DS_data.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-21
#'
#' Script description:
#'
#'
#' -------------------------------------------------------------------------------


sim_DS_data <- function (type = c("line", "point"),
                         nsites = 50,
                         Ntot = 200,
                         mean.sigma = 0.15,
                         dist_max = 0.6,
                         transect_len = 2,
                         size_study_area = 250,
                         discard0 = TRUE){
  type <- match.arg(type)
  size_area_sampled <- 2 * dist_max * transect_len
  
  lambda <- Ntot * size_area_sampled / size_study_area
  N <- rpois(nsites, lambda)
  N.true <- N
  sigma <- rep(mean.sigma, nsites)
  data <- NULL
  for (i in 1:nsites) {
    if (N[i] == 0) {
      data <- rbind(data, c(i, NA, NA, NA, NA))
      next
    }
    if (type == "line") {
      d <- runif(N[i], 0, dist_max)
      p <- exp(-d * d / (2 * (sigma[i] ^ 2)))
      y <- rbinom(N[i], 1, p)
      u <- v <- rep(NA, N[i])
      d <- d[y == 1]
      u <- u[y == 1]
      v <- v[y == 1]
      y <- y[y == 1]
    }
    # if (type == "point") {
    #   u <- runif(N[i], 0, 2 * dist_max)
    #   v <- runif(N[i], 0, 2 * dist_max)
    #   d <- sqrt((u - dist_max) ^ 2 + (v - dist_max) ^ 2)
    #   N.true[i] <- sum(d <= dist_max)
    #   p <- exp(-d * d / (2 * (sigma[i] ^ 2)))
    #   pp <- ifelse(d <= dist_max, 1, 0) * p
    #   y <- rbinom(N[i], 1, pp)
    #   u <- u[y == 1]
    #   v <- v[y == 1]
    #   d <- d[y == 1]
    #   y <- y[y == 1]
    # }
    if (sum(y) > 0)
      data <- rbind(data, cbind(rep(i, sum(y)), y, u, v, d))
    else
      data <- rbind(data, c(i, NA, NA, NA, NA))
  }
  colnames(data) <- c("site", "y", "u", "v", "d")
  if (discard0)
    data <- data[!is.na(data[, 2]), ]
  
  list(
    type = type,
    nsites = nsites,
    mean.sigma = mean.sigma,
    dist_max = dist_max,
    transect_len = transect_len,
    size_study_area = size_study_area,
    data = data,
    N = N,
    N.true = N.true
  )
}

