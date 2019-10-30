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
forest_plot_base <- function (params) {

}

# ad points and bars to a forest plot
forest_plot_add <- function (params) {

  lines()
  points()

}

# Produce forest plots for parameters in the protea DSDM, with a panerl each for the
# informative and noninformative priors.  Priors in grey in the backgound,
# posteriors in the foreground.
plot_estimates <- function (protea_estimates, filename) {

  on.exit(dev.off())
  png(filename,
      width = 1000,
      height = 700,
      pointsize = 30)

  # non-informative prior plot
  par(mfrow = c(1, 2))

  forest_plot_base()
  forest_plot_add(noninformative_params)

  forest_plot_base()
  forest_plot_add(informative_params)

  dev.off()

}
