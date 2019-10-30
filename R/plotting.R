# turn a (continuous, ascending) colorBrewer palette into a color ramp function
pal <- function (name, flip = FALSE) {
  cols <- brewer.pal(9, name)
  if (flip) cols <- rev(cols)
  colorRampPalette(cols)
}

# plot a map for a single covariate
plot_cov <- function (coords, values, pal, title) {

  # build a raster
  data <- data.frame(values)
  coordinates(data) <- coords
  spdf <- SpatialPixelsDataFrame(data, tol = 0.0001, data = data.frame(data))
  ras <- raster(spdf, values = TRUE)

  # plot
  plot(ras, col = pal(100),
       legend = FALSE,
       axes = FALSE,
       box = FALSE)

  # add title
  mtext(
    title,
    side = 3,
    line = -5,
    cex = 1.4,
    adj = 0.75
  )

}

# plot all the covariates
plot_covs <- function (occurrence, filename) {

  on.exit(dev.off())
  png(filename,
      width = 1000 * 2,
      height = 700 * 4,
      pointsize = 30)

  covs <- colnames(occurrence)[-(1:5)]
  coords <- occurrence[, c("lon", "lat")]

  palettes <- list(
    min_temp_july = pal("YlOrRd"),
    summer_smd = pal("YlGnBu"),
    winter_smd = pal("YlGnBu"),
    prop_acidic = pal("PuRd"),
    prop_fertile = pal("YlGn"),
    fire_return = pal("Reds", TRUE),
    mean_annual_precip = pal("Blues")
  )

  titles <- list(
    min_temp_july = "minimum temperature in July",
    summer_smd = "soil moisture days in summer",
    winter_smd = "soil moisture days in winter",
    prop_acidic = "proportion soil acidic",
    prop_fertile = "proportion soil highly fertile",
    fire_return = "mean fire return time",
    mean_annual_precip = "mean annual precipitation"
  )

  # temporarily set the plot panel configuration
  mfrow <- par()$mfrow
  on.exit(par(mfrow = mfrow), add = TRUE)
  par(mfrow = c(4, 2))
  for (cov in covs) {
    plot_cov(
      coords = coords,
      values = occurrence[, cov],
      pal = palettes[[cov]],
      title = titles[[cov]]
    )
  }

}

# plot presence/absence
plot_occ <- function (occurrence, filename) {

  on.exit(dev.off())
  png(filename,
      width = 1000,
      height = 700,
      pointsize = 30)

  coords <- occurrence[, c("lon", "lat")]
  pres <- occurrence$presence

  bg_pal <- function (n) rep(grey(0.9), n)

  # plot a grey background
  plot_cov(
    coords = coords,
    values = rep(0.9, nrow(coords)),
    pal = bg_pal,
    title = ""
  )

  # add on presence/absence points
  points(
    coords,
    pch = 16,
    col = ifelse(pres, "blue", "green"),
    cex = 0.2
  )

  legend(
    21, -31.5,
    c("present", "absent"),
    pch = 16,
    col = c("blue", "green"),
    box.lty = 0,
    horiz = TRUE
  )

}

# background for a forest plot
forest_base <- function (params) {
  locs <- seq(1, 0, length.out = nrow(params))
  plot(locs ~ mean,
       data = params,
       type = "n",
       ylab = "",
       xlab = "",
       axes = FALSE,
       bty = "n",
       ylim = c(0, 1),
       xlim = range(c(params$lower, params$upper)))
  abline(v = 0, lwd = 1, col = grey(0.5), lty = 3)
  axis(side = 1, tcl = -0.3)
  axis(side = 2, lty = 0, at = locs, labels = params$param, las = 2)
  title(xlab = "estimate", col.lab = grey(0.3))
}

# add points and bars to existing forest plot
forest_points <- function(params, col = grey(0.8), scaling = 1) {

  locs <- seq(1, 0, length.out = nrow(params))
  for (i in seq_len(nrow(params))) {
    lines(x = c(params$lower[i], params$upper[i]),
          y = cbind(locs[i], locs[i]),
          lwd = 1.5 * scaling ^ 2,
          col = col)
  }

  points(locs ~ params$mean,
         pch = 16,
         cex  = 1 * sqrt(scaling),
         col = col)

}

# Produce forest plots for parameters in the protea DSDM, with a panerl each for the
# informative and noninformative priors.  Priors in grey in the backgound,
# posteriors in the foreground.
plot_estimates <- function(
  naive_draws_summary,
  informative_draws_summary,
  filename
) {

  on.exit(dev.off())
  png(filename,
      width = 1000,
      height = 700,
      pointsize = 30)

  par(mfrow = c(1, 2))

  # naive prior plot
  forest_base(naive_params)
  title(main = "with naive priors", col.main = grey(0.3))
  forest_points(naive_priors, col = grey(0.9), scaling = 2)
  forest_points(naive_params, col = grey(0.3), scaling = 1)

  # informative prior plot
  forest_base(informative_params)
  title(main = "with informative priors", col.main = grey(0.3))
  forest_points(informative_priors, col = grey(0.9), scaling = 2)
  forest_points(informative_params, col = grey(0.3), scaling = 1)

}
