# generate mini figures for schematic figure
set.seed(2019-09-24)

# covariate data
library(raster)
library(fields)

# download a shapefule for Victoria, Australia and convert it into a raster
states <- raster::getData("GADM", country = "AUS", level = 1)
vic <- states[states$NAME_1 == "Victoria", ]
r <- raster(vic)
res(r) <- c(0.1, 0.1)
template <- rasterize(vic, r)
idx <- which(!is.na(getValues(template)))

# generate covariates (Gaussian processes with exponential covariance, as usual)
coords <- xyFromCell(template, idx)
n_cells <- length(idx)
v <- matrix(rnorm(n_cells * 3), n_cells, 3)
d <- rdist.earth(coords)
Sigma <- exp(-d / 40) + diag(n_cells) * 1e-3
L <- t(chol(Sigma))
f <- L %*% v
covs <- stack(template, template, template)
covs[idx] <- f

# vital rate rasters
fec <- surv <- template
X <- cbind(1, covs[idx])
beta_fec <- c(0, rnorm(3, 0, 0.1))
fec[idx] <- exp(X %*% beta_fec)[, 1]
beta_surv <- c(-1, rnorm(3, 0, 1))
surv[idx] <- plogis(X %*% beta_surv)[, 1]

# distribution data
# analytic solution to thedominant eigenvalue for a simple MPM
lambda <- sqrt(surv * (fec  * 0.7 + 1))
prob <- 1 - exp(-0.5 * lambda)
loc <- sample.int(n_cells, 100)
pa <- rbinom(length(loc), 1, prob[idx][loc])
points <- cbind(coords[loc, ], pa)

# predicted distribution
pred <- template
m <- glm(pa ~ ., data = as.data.frame(X[loc, -1]), family = binomial)
pred[idx] <- predict(m, as.data.frame(X[, -1]), type = "response")

# plots
plot_raster <- function (r, col = cols()) {
  image(
    r,
    maxpixels = Inf,
    col = col,
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    bg = NA
  )
}
cols <- function (pal = viridis::magma, n = 1000, rev = TRUE, fraction = 0.5) {
  cols <- pal(n)
  if (rev) {
    cols <- rev(cols)
  }
  range <- seq_len(floor(n * fraction))
  cols[range]
}

render <- function (raster, cols, name, points = NULL) {
  op <- par(no.readonly = TRUE)
  on.exit({dev.off(); par(op)})
  png(paste0(name, ".png"), width = 500, height = 500, bg = NA)
  par(mar = rep(0.1, 4))
  plot_raster(raster, cols)
  if (!is.null(points)) {
    points(points[, 1:2], pch = 21, lwd = 3, cex = 1.5,
           bg = ifelse(points[, 3], "black", "white"),
           col = grey(0.3))
  }
}

# covariates
render(covs[[1]], cols(), "cov1")
render(covs[[2]], cols(viridis::viridis, rev = FALSE), "cov2")
render(covs[[3]], cols(viridis::plasma, rev = FALSE, fraction = 0.8), "cov3")

# vital rates
fec_pal <- colorRampPalette(c(grey(0.9), "magenta"))
render(fec, fec_pal(1000), "fec")
surv_pal <- colorRampPalette(c(grey(0.9), "seagreen"))
render(surv, surv_pal(1000), "surv")

# occurrence data
render(template, grey(0.95), "occ", points = points)

# prediction
pred_pal <- colorRampPalette(c(grey(0.9), "slateblue"))
render(pred ^ 2, pred_pal(1000), "pred")


