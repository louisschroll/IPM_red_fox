
plot_DS_data <- function(DS_data){
  with (DS_data, {
    if (type == "line") {
      op <- par(mfrow = c(1, 3))
      on.exit(par(op))
      tryPlot <- try({
        hist(
          data[, "distance"],
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
          xlim = c(0, 2 * dist_max),
          ylim = c(0, 2 * dist_max),
          col = data[, 1],
          asp = 1
        )
        points(dist_max,
               dist_max,
               pch = "+",
               cex = 3,
               col = "black")
        draw.circle(dist_max, dist_max, dist_max)
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
    }
  })
}
