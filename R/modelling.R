# modelling functions for both the BBS and protea analyses

# given a vector greta array, where each element corresponds to a cell in the
# raster, an mcmc.list object of posterior samples, return a raster giving the
# posterior mean of the greta array in each cells
map_variable <- function (greta_array, draws, raster_template) {
  vals <- calculate(greta_array, draws)
  vals_mat <- as.matrix(vals)
  mean <- colMeans(vals_mat)
  map <- raster_template$raster
  map[raster_template$idx] <- mean
  map
}

# use all covariates in both models
surv_covs <- fec_covs <- all_covs <- c("pcMix", "pcDec", "pcCon", "bio1", "bio6", "bio12")

# calculate survival, given covariates
get_survival <- function (X, params) {
  survival_logit <- X %*% params$beta_survival
  survival_adult <- ilogit(survival_logit)
  survival_juvenile <- ilogit(survival_logit - params$gamma)
  list(adult = survival_adult,
       juvenile = survival_juvenile)
}

get_fecundity <- function (X, params) {
  fecundity_log <- X %*% params$beta_fecundity
  exp(fecundity_log)
}

get_lambda <- function(survival, fecundity, analytic = FALSE) {

  # use either the analytic solution (specific to this simple model) or the
  # numeric solution (which is more general to DSDMs)
  if (analytic) {

    lambda <- fecundity * survival$juvenile + survival$adult

  } else {

    # build n x m x m 3D array containing Leslie matrices for different locations
    top_row <- cbind(fecundity * survival$juvenile,
                     fecundity * survival$adult)
    bottom_row <- cbind(survival$juvenile,
                        survival$adult)
    matrices <- abind(top_row, bottom_row, along = 3)

    # loop across observations, iterating matrices to get implied intrinsic growth
    iterated <- greta.dynamics::iterate_matrix(matrices)
    lambda <- iterated$lambda

  }

  lambda

}

# calculate the carrying capacity from lambda and scaling parameters, rectified
# to always be positive (and differentiable)

# probability of observed occupancy from growth rate and ancillary parameters
get_prob <- function (lambda, params) {

  # (relative) Beverton-Holt theoretical carrying capacity (Ricker would be
  # lambda * alpha)
  K <- (lambda - 1) * params$alpha

  # rectify to positive relative abundance for all values
  Lambda <- log1pe(K)

  # probability of detection under Poisson sampling
  1 - exp(-Lambda)

}

get_design_matrix <- function(formula, vals) {
  vals_df <- as.data.frame(vals)
  model.matrix(formula, vals_df)
}

# randomly split data into training and test datasets, and compute other useful data things
prep_data <- function (formula, occ, covs, maps_coords) {

  # template raster with minimal mask and index to cells with values
  raster_template <- list(
    raster = covs[[7]] * 0,
    idx = which(!is.na(getValues(covs[[7]])))
  )

  # scale covs to have mean 0 and sd 1 over the whole study area. Handle the
  # prior via an affine transformation beta prior
  covs <- scale(covs)

  # get covariate values for maps locations
  maps_vals <- extract(covs, maps_coords[, c("lon", "lat")])
  x_maps <- get_design_matrix(formula, maps_vals)

  # for all pixel locations
  all_vals <- extract(covs, raster_template$idx)
  x_predict <- get_design_matrix(formula, all_vals)

  # and for all BBS locations
  coords <- occ[, c("Longitude", "Latitude")]
  bbs_vals <- extract(covs, coords)
  x_bbs <- get_design_matrix(formula, bbs_vals)

  # combine, drop NAs, then split again
  all <- cbind(occ[, -(1:2)], x_bbs)
  all <- na.omit(all)

  pa <- all$PA
  coords <- all[, c("Longitude", "Latitude")]
  x_bbs <- all[, colnames(x_bbs)]

  # subset, remove a text column, and convert to a matrix
  train_idx <- sample.int(nrow(all), 500)

  list(
    x_train = x_bbs[train_idx, ],
    x_test = x_bbs[-train_idx, ],
    pa_train = pa[train_idx],
    pa_test = pa[-train_idx],
    coords_train = coords[train_idx, ],
    coords_test = coords[-train_idx, ],
    covs = covs,
    x_maps = x_maps,
    x_predict = x_predict,
    raster_template = raster_template
  )

}

# given the mean and standard deviation of a lognormally-distributed random variable,
# compute the mean and variance of its log
lognormal_prior <- function(mean, sd) {

  var <- sd ^ 2
  list(
    mean = log((mean ^ 2) / sqrt(var + mean ^ 2)),
    sd = sqrt(log(1 + var / (mean ^ 2)))
  )

}

# given the mean and standard deviation of a logitnormally-distributed random variable,
# compute the mean and variance of its logit (a normal random variable)
logitnormal_prior <- function(mean, sd) {

  # given the parameters of the distribution, estimate the moments and hence the
  # mean and variance via the quasi Monte Carlo method listed on wikipedia
  logitnormal_mean_var <- function (mu, sd, K = 1000) {
    i <- 1:(K - 1)
    p <- plogis(qnorm(i/K, mu, sd))
    m1 <- mean(p)
    m2 <- mean(p ^ 2)
    c(mean = m1, variance = m2 - m1 ^ 2)
  }

  # compute the root mean squared error between the logitnormal mean and
  # variance from this (transformed) set of parameters, and the observed mean
  # and variance
  obj <- function (par, observed) {
    mu <- par[1]
    sd <- exp(par[2])
    expected <- logitnormal_mean_var(mu, sd)
    mean((expected - observed) ^ 2)
  }

  # optimise this to find the parameter set that minimised the mean squared error
  op <- optim(par = c(0, 0), obj, observed = c(mean, sd ^ 2))

  # return the mean and variance on the logit scale
  list(mean = op$par[1],
       sd = exp(op$par[2]))

}

# estimate the prior for gamma, from the difference between two normal
# distributions. this is another normal distribution, with mean mu_A - mu_J, and
# variance var_A + var_J
gamma_prior <- function (logit_S_A,
                         logit_S_J,
                         n = 1e5) {

  list(mean = logit_S_A$mean - logit_S_J$mean,
       sd = sqrt(logit_S_A$sd ^ 2 + logit_S_J$sd ^ 2))
}

# given a vector of prior means and standard deviations for some variables z,
# where z = X \beta, compute the mean and covariance of a prior over beta that
# implies this distribution over z.
beta_prior <- function(X, mean, sd) {

  # use weights to increase variance for each datapoint to spread
  # prior over multiple datapoints, without adding too much extra data

  # convert the means and sds into MVN parameters
  n <- nrow(X)

  if (length(mean) == 1) {
    mean <- rep(mean, n)
  }
  if (length(sd) == 1) {
    sd <- rep(sd, n)
  }

  var <- sd ^ 2
  cov <- diag(var)

  # invert X
  iX <- MASS::ginv(X)

  # solve for beta
  list(mean = iX %*% mean,
       Sigma = iX %*% cov %*% t(iX))

}

# greta variables with prior distributions for the regression coefficients for
# the two submodels, using whitened representation to reduce correlations in the
# posteriors
beta_parameter <- function(beta_prior) {
  L <- t(chol(beta_prior$Sigma))
  beta_raw <- normal(0, 1, dim = length(beta_prior$mean))
  beta_prior$mean + L %*% beta_raw
}


build_bbs_model <- function (data_list, analytic = FALSE) {

  # get thee parameters of normally-distributed vital rate priors on the link
  # scale, based on estimates from Ryu et al.
  log_fecundity_prior <- lognormal_prior(0.58, 0.25)
  logit_adult_survival_prior <- logitnormal_prior(0.71, 0.195)
  logit_juvenile_survival_prior <- logitnormal_prior(0.489, 0.195)

  # compute MVN priors for survival and fecundity betas (including intercept)
  priors <- list(
    beta_fecundity = beta_prior(
      data_list$x_maps,
      log_fecundity_prior$mean,
      log_fecundity_prior$sd
    ),
    beta_survival = beta_prior(
      data_list$x_maps,
      logit_adult_survival_prior$mean,
      logit_adult_survival_prior$sd
    ),
    gamma = gamma_prior(
      logit_adult_survival_prior,
      logit_juvenile_survival_prior
    )
  )

  # parameters
  params <- list(
    # submodel regression coefficients
    beta_fecundity = beta_parameter(priors$beta_fecundity),
    beta_survival = beta_parameter(priors$beta_survival),
    # no prior on the combined detection/competition effect
    alpha = variable(lower = 0),
    # uncertainty in the logit-difference between adult and juvenile survival
    gamma = normal(
      priors$gamma$mean,
      priors$gamma$sd
    )
  )

  # get all the components
  survival <- get_survival(data_list$x_train, params)
  fecundity <- get_fecundity(data_list$x_train, params)
  lambdas <- get_lambda(survival, fecundity, analytic = analytic)
  prob <- get_prob(lambdas, params)

  # likelihood for presence-absence data
  presence <- data_list$pa_train
  distribution(presence) <- bernoulli(prob)

  # fit the model
  m <- model(
    params$beta_fecundity,
    params$beta_survival,
    params$alpha,
    params$gamma
  )

  list(
    model = m,
    params = params,
    priors = priors
  )

}

make_bbs_predictions <- function (model_list, draws_list, data_list) {

  params <- model_list$params
  draws <- draws_list$draws
  raster_template <- data_list$raster_template

  # get prediction greta arrays
  survival_predict <- get_survival(data_list$x_predict, params)
  fecundity_predict <- get_fecundity(data_list$x_predict, params)
  lambdas_predict <- get_lambda(survival_predict, fecundity_predict, analytic = TRUE)
  prob_predict <- get_prob(lambdas_predict, params)

  # map posterior mean estimates of these
  prob_pres_map <- map_variable(prob_predict, draws, raster_template)
  fec_map <- map_variable(fecundity_predict, draws, raster_template)
  surv_map <- map_variable(survival_predict$adult, draws, raster_template)

  # decrease survival by 50% across N America
  survival_predict_low <- lapply(survival_predict, "*", 0.75)
  lambdas_low_surv <- get_lambda(survival_predict_low, fecundity_predict, analytic = TRUE)

    # decrease fecundity by 50% across N Americ
  lambdas_low_fec <- get_lambda(survival_predict, fecundity_predict * 0.5, analytic = TRUE)

  # compute the posterior probability of lambda exceeding 1, then plot
  in_niche <- lambdas_predict >= 1
  in_niche_low_surv <- lambdas_low_surv >= 1
  in_niche_low_fec <- lambdas_low_fec >= 1
  prob_in_niche_map <- map_variable(in_niche, draws, raster_template)
  prob_in_niche_low_surv_map <- map_variable(in_niche_low_surv, draws, raster_template)
  prob_in_niche_low_fec_map <- map_variable(in_niche_low_fec, draws, raster_template)

  list(
    prob_in_niche = prob_in_niche_map,
    prob_present = prob_pres_map,
    fecundity = fec_map,
    survival = surv_map,
    prob_in_niche_low_survival = prob_in_niche_low_surv_map,
    prob_in_niche_low_fecundity = prob_in_niche_low_fec_map
  )

}

plot_bbs_maps <- function (predictions, data_list) {

  # get range limits
  predictions$range <- predictions$prob_in_niche > 0.5
  predictions$low_surv_range <- predictions$prob_in_niche_low_survival > 0.05
  predictions$low_fec_range <- predictions$prob_in_niche_low_fecundity > 0.05

  # set up plotting colours and scales
  Dark2 <- brewer.pal(3, "Dark2")
  YlGnBu <- brewer.pal(9, "YlGnBu")[c(1, 1, 2, 3, 4, 5, 6, 7, 7)]
  YlGn <- brewer.pal(9, "YlGn")
  RdPu <- brewer.pal(9, "RdPu")
  prob_cols <- colorRampPalette(YlGnBu)(1000)
  surv_cols <- colorRampPalette(YlGn)(1000)
  fec_cols <- colorRampPalette(RdPu)(1000)

  legend_width <- 2
  title_cex <- 1.3
  title_col <- grey(0.4)
  title_line <- -0.5

  baseline_col <- "lightblue"
  survival_col <- surv_cols[750]
  fecundity_col <- fec_cols[400]

  png("figures/bbs_fig.png",
      width = 2400,
      height = 1600,
      pointsize = 45)

  par(mfrow = c(2, 2),
      oma = c(0, 1, 1, 1),
      mar = c(1, 1, 1, 1))

  # probability of presence from the DDM
  plot(predictions$prob_present,
       zlim = c(0, 1),
       col = prob_cols,
       axes = FALSE,
       box = FALSE,
       legend.width = legend_width,
       maxpixels = Inf)

  points(data_list$coords_train,
         pch = 21,
         cex = 0.7,
         bg = "white",
         lwd = 4,
         col = grey(0.4))

  points(data_list$coords_train[data_list$pa_train == 1, ],
         pch = 16,
         cex = 0.7,
         col = "black")

  title(main = "probability of presence",
        cex.main = title_cex,
        col.main = title_col,
        line = title_line)
  mtext("A", 3, adj = 0, cex = 1.2)

  # plot fecundity

  # clip extrapolation from south east, for plotting only
  fec_map_plot <- predictions$fecundity
  fec_max <- 8
  fec_map_plot[fec_map_plot > fec_max] <- fec_max

  plot(fec_map_plot,
       zlim = c(0, maxValue(fec_map_plot)),
       axes = FALSE,
       box = FALSE,
       maxpixels = Inf,
       col = fec_cols,
       legend.width = legend_width)

  title(main = "fecundity",
        cex.main = title_cex,
        col.main = title_col,
        line = title_line)
  mtext("B", 3, adj = 0, cex = 1.2)

  # plot adult survival
  plot(predictions$survival,
       zlim = c(0, 1),
       axes = FALSE,
       box = FALSE,
       col = surv_cols,
       maxpixels = Inf,
       legend.width = legend_width)

  title(main = "adult survival",
        cex.main = title_cex,
        col.main = title_col,
        line = title_line)
  mtext("C", 3, adj = 0, cex = 1.2)

  # range contractions due to demographic changes (e.g. introduced competitor or
  # pathogen); a form of spatial sensitivity analysis

  # expected range under current conditions
  plot(predictions$range,
       col = c(grey(0.9), baseline_col),
       legend = FALSE,
       axes = FALSE,
       maxpixels = Inf,
       box = FALSE)

  # add range under decreased fertility
  predictions$low_fec_range[predictions$low_fec_range == 0] <- NA
  plot(predictions$low_fec_range,
       col = fecundity_col,
       add = TRUE,
       maxpixels = Inf,
       legend = FALSE)

  # add range under decreased survival
  predictions$low_surv_range[predictions$low_surv_range == 0] <- NA
  plot(predictions$low_surv_range,
       col = survival_col,
       add = TRUE,
       maxpixels = Inf,
       legend = FALSE)

  title(main = "range limits",
        cex.main = title_cex,
        col.main = title_col,
        line = title_line)
  mtext("D", 3, adj = 0, cex = 1.2)

  text(x = -71, y = 35.5,
       labels = "baseline",
       col = "lightblue3",
       cex = 1.4,
       xpd = NA)

  text(x = -71, y = 32,
       labels = "low fecundity",
       col = fecundity_col,
       cex = 1.4,
       xpd = NA)

  text(x = -71, y = 29,
       labels = "low survival",
       col = survival_col,
       cex = 1.4,
       xpd = NA)

  dev.off()

}

loglik <- function (pred, pa) {
  sum(dbinom(pa, 1, pred, log = TRUE))
}

auc <- function(pred, pa) {
  as.numeric(pROC::auc(pa, pred))
}

# compare predictive ability of DSDM and GLM, and correlation between predictions.
compare_bbs_predictions <- function (predictions, data_list) {

  prob_raster <- predictions$prob_present

  # extract predictions from DSDM
  dsdm_pred <- raster::extract(prob_raster, data_list$coords_test)

  # make predictions from GLM (no intercept as it's in all)
  glm <- glm(data_list$pa_train ~ . -1,
             data = as.data.frame(data_list$x_train),
             family = stats::binomial)
  glm_pred <- predict(glm,
                      newdata = as.data.frame(data_list$x_test),
                      type = "response")

  list(
    glm_dsdm_correl = cor(glm_pred, dsdm_pred),
    glm_pa_correl = cor(glm_pred, data_list$pa_test),
    dsdm_pa_correl = cor(dsdm_pred, data_list$pa_test),
    glm_auc = auc(glm_pred, data_list$pa_test),
    dsdm_auc = auc(dsdm_pred, data_list$pa_test)
  )

}

op <- greta::.internals$nodes$constructors$op

# normalise a matrix rowwise in exponential space. Ie. Given a matrix x, return a
# matrix y where:
# exp(y[i, ]) = exp(x[i, ]) / sum(exp(x[i, ]))
# or:
# y[i, ] <- log(exp(x[i, ]) / sum(exp(x[i, ])))
# or:
# y[i, ] <- x[i, ] -  log(sum(exp(x[i, ])))
log_normalise <- function (x) {

  tf_log_normalise <- function (x) {
    # skipping batch dimension, do log_sum_exp on rows, then broadcast subtract it
    sums <- tensorflow::tf$reduce_logsumexp(x, axis = 2L, keepdims = TRUE)
    x - sums
  }

  op("log_normalise", x, tf_operation = "tf_log_normalise")

}

# load the parameter estimates from Merow et al. and format them as the
# parameters of normal priors
protea_prior <- function (model = c("seedling_survival", "adult_survival", "growth",
                             "flowering", "seedheads", "seeds", "germination", "offspring_size"),
                   type = c("environment", "size"),
                   informative = TRUE) {

  # match the arguments
  model <- match.arg(model)
  type <- match.arg(type)

  # load parameter estimates
  estimates <- read.csv("data/clean/parameter_estimates.csv")

  # subset to model of interest (and drop model column)
  estimates <- estimates[estimates$model == model, -1]

  # if they want the size coefficient, check it's available and use that if so. retun a NULL if not
  if (type == "size") {
    estimates <- estimates[estimates$parameter == "size", ]
    if (nrow(estimates) == 0) {
      return(NULL)
    }
  } else {
    # otherwise remove the size coefficient (if it was there at all)
    estimates <- estimates[estimates$parameter != "size", ]
  }

  # pull out means
  means <- estimates$mean

  # get approximate standard deviations for these parameters
  sds <- with(
    estimates,
    ((upper - mean) + (mean - lower) / 2) / 1.96
  )

  names(means) <- names(sds) <- estimates$parameter

  # if we wanted non-informative ones, replace them now (so that we have the names)
  if (!informative) {
    means[] <- 0
    sds[] <- 10
  }

  list(mean = means,
       sd = sds)

}

# return a quadratic version of a variable, for building design matrices
quad <- function (x) {
  call <- match.call()
  name <- as.character(call[[2]])
  mat <- cbind(x, x ^ 2)
  colnames(mat) <- c("_raw", "_squared")
  mat
}

# given a formula, return a greta array of regression coefficients, and a vector
# of them multiplied by their corresponding design matrix. Retrun both in a list
get_env_effect <- function (formula, prior, data) {

  design <- model.matrix(formula, data)

  # check the prior names line up with the design matrix
  prior_names <- names(prior$mean)
  data_names <- colnames(design)
  if (!identical(prior_names, data_names)) {
    stop ("name mismatch between prior and design matrix",
          call. = FALSE)
  }

  beta <- normal(prior$mean, prior$sd)

  # # for intercept-only models, don't multiply by the intercept column, it's more
  # # efficient (when combining the kernels) to return just the scalar
  # if (identical(data_names, "(Intercept)")) {
  #   env_effect <- beta
  # } else {
  #   env_effect <- design %*% beta
  # }

  env_effect <- design %*% beta

  list(beta = beta, env_effect = env_effect)
}

# if the relationship is size-dependent, get a coefficient for that, and return
# along with the (row-vector) offsets for all size classes
get_size_effect <- function (size_dependent, prior, sizes) {

  coef <- 0

  if (size_dependent) {
    coef <- normal(prior$mean, prior$sd)
  }

  list(coef = coef, size_effect = coef * t(sizes))
}

apply_link <- function (link_function, matrix) {
  link_function(matrix)
}

# blank function for no link/transformation
none <- function (x) x

# construct kernel array for transitions with normal distributions, where each
# slice of the array represents a kernel matrix for a given site. Given a matrix
# of means at sites and sizes (n_sites x matrix_dim), a scalar standard
# deviation for the distribution, a vector of sizes to integrate over (of length
# matrix_dim), and the integration bin size, return a kernel array (n_sites x
# matrix_dim x matrix_dim) giving the transition probabilities in all cells.
distribution_kernel <- function (mean_matrix, sd, sizes, bin_width) {
  n_sites <- nrow(mean_matrix)
  matrix_dim <- ncol(mean_matrix)
  # work with the log density for computational stability
  var <- sd ^ 2
  log_constant <- -0.5 * log(2 * pi) - log(sd)
  diff <- kronecker(sizes, mean_matrix, FUN = "-")
  log_density <- log(bin_width) + log_constant - (diff ^ 2) / (2 * var)

  # normalise this in exponential-space (hope we get it the right way round!)
  log_density <- log_normalise(log_density)

  # diff was created as a matrix with dimensions (n_sites * matrix_dim) x
  # matrix_dim, so we reshape it (to n_sites x matrix_dim x matrix_dim), then
  # transpose all the matrices, to that columns sum to 1
  dim(log_density) <- c(n_sites, matrix_dim, matrix_dim)
  log_density <- aperm(log_density, c(1, 3, 2))
  exp(log_density)
}

# use greta::calculate() to sanity check the effects of different parameter
# values on kernel components
check_vals <- function (name, beta_values = NULL, coef_value = NULL) {
  beta <- betas[[name]]
  coef <- size_coefs[[name]]
  beta_vals <- rnorm(length(beta))
  coef_val <- rnorm(1)
  if (!is.null(beta_values)) beta_vals[] <- beta_values
  if (!is.null(coef_value)) coef_val[] <- coef_value
  calculate(matrices[[name]], list(beta = beta_vals, coef = coef_val))
}

# need to combine the survival components, with linear interpolation in the size
# range in which adults and seedlings overlap. do that by building a weighting
# data vector and doing a weighted sum of the two
interpolate_survival <- function (adult_survival_matrix,
                                  seedling_survival_matrix,
                                  sizes,
                                  adult_min_size = 0.182,
                                  seedling_max_size = 0.728) {

  # build weights vector
  range <- seedling_max_size - adult_min_size
  weights <- seedling_max_size / range  - (sizes / range)
  weights[sizes > seedling_max_size] <- 0
  weights[sizes <= adult_min_size] <- 1

  # apply weights
  weighted_seedling_survival_matrix <- sweep(seedling_survival_matrix, 2, weights, FUN = "*")
  weighted_adult_survival_matrix <- sweep(adult_survival_matrix, 2, 1 - weights, FUN = "*")

  # sum and return
  survival_matrix <- weighted_seedling_survival_matrix + weighted_adult_survival_matrix
  survival_matrix

}

# given an n x m x m array A, and an n x m matrix B, sweep B along the third
# dimension of A using the function FUN, to yield an array with the same
# dimensions as A . Note that this is only necessary because greta's sweep
# doesn't yet handle arrays (with more than 2 dimensions), otherwise we could
# just do: sweep(A, B, MARGIN = c(1, 3), FUN = FUN)
sweep_3 <- function(A, B, FUN = c("-", "+", "/", "*")) {
  n <- dim(A)[1]
  m <- dim(A)[2]
  Al <- aperm(A, c(2, 3, 1))
  dim(Al) <- c(m, n * m)
  Bl <- t(B)
  dim(Bl) <- n * m
  C <- sweep(Al, 2, Bl, FUN)
  dim(C) <- c(m, m, n)
  aperm(C, c(3, 1, 2))
}

# build and return a greta model for the DSDM
build_protea_model <- function (occurrence, informative_priors = FALSE) {

  # subset occurrence data to where we observe presence/absence, and subset to
  # a few observations (for now)
  occurrence %>%
    filter(!is.na(presence)) %>%
    sample_n(500) -> occ_train

  # formulae for the *environmental components* of the vital rate regressions
  formulae <- list(
    seedling_survival = ~ min_temp_july + mean_annual_precip,
    adult_survival = ~ min_temp_july + mean_annual_precip,
    growth = ~ quad(prop_fertile) + quad(winter_smd) + quad(prop_acidic) + min_temp_july + quad(summer_smd),
    flowering = ~ 1,
    seedheads = ~ min_temp_july + quad(prop_acidic),
    seeds = ~ 1,
    germination = ~ 1,
    offspring_size = ~ prop_acidic + quad(prop_fertile) + min_temp_july + quad(summer_smd)
  )

  beta_priors <- lapply(
    names(formulae),
    protea_prior,
    type = "environment",
    informative = informative_priors
  )
  names(beta_priors) <- names(formulae)

  # for each formula, get the vector of coefficients and environmental effects
  env_effects_list <- mapply(
    get_env_effect,
    formulae, beta_priors,
    MoreArgs = list(data = occ_train),
    SIMPLIFY = FALSE
  )

  betas <- lapply(env_effects_list, `[[`, "beta")
  env_effects <- lapply(env_effects_list, `[[`, "env_effect")

  # set up 100-cell kernel as in Merow et al.
  min_size <- 0
  max_size <- 5
  matrix_dim <- 5
  bounds <- min_size + (0:matrix_dim) * (max_size - min_size) / matrix_dim
  sizes <- 0.5 * (bounds[-(matrix_dim + 1)] + bounds[-1])
  bin_width <- sizes[2] - sizes[1]

  # get size effects for flowering probability, seedheads, and both survivals
  size_dependent <- list(
    seedling_survival = TRUE,
    adult_survival = TRUE,
    growth = FALSE,
    flowering = TRUE,
    seedheads = TRUE,
    seeds = FALSE,
    germination = FALSE,
    offspring_size = FALSE
  )

  size_priors <- lapply(
    names(size_dependent),
    protea_prior,
    type = "size",
    informative = informative_priors
  )

  size_effects_list <- mapply(
    get_size_effect,
    size_dependent, size_priors,
    MoreArgs = list(sizes = sizes),
    SIMPLIFY = FALSE
  )

  size_coefs <- lapply(size_effects_list, `[[`, "coef")
  size_effects <- lapply(size_effects_list, `[[`, "size_effect")

  # There's a deterministic effect of size on growth; we need to add the origin
  # size onto the (environmentally-driven) growth increment to get the
  # destination size. most efficient to do that here.
  size_effects$growth <- size_effects$growth + sizes

  # combine the size and environment effects (get the outer sum of the two) to
  # get matrices on the link scale
  link_scale_matrices <- mapply(
    kronecker,
    size_effects, env_effects,
    MoreArgs = list(FUN = "+"),
    SIMPLIFY = FALSE
  )

  # transform to parameters to the response scale to build kernel array
  link_functions <- list(
    seedling_survival = ilogit,
    adult_survival = ilogit,
    growth = none,
    flowering = ilogit,
    seedheads = exp,
    seeds = exp,
    germination = ilogit,
    offspring_size = none
  )
  # N.B. we want germination probability to be constrained 0-1, but it was not
  # estimated from a logit model so the estimates in the data are a
  # logit-transformed version of the estimates in the paper and are transformed
  # back here

  matrices <- mapply(
    apply_link,
    link_functions, link_scale_matrices,
    SIMPLIFY = FALSE
  )

  # get coefficients for variation in growth and offspring size (priors for these?)
  growth_sd <- normal(0, 10, truncation = c(0, Inf))
  offspring_size_sd <- normal(0, 10, truncation = c(0, Inf))

  # get growth and offspring size transition probabilities (integrating over
  # discretized normal distributions) to return kernel arrays
  growth_density <- distribution_kernel(
    mean_matrix = matrices$growth,
    sd = growth_sd,
    sizes = sizes,
    bin_width = bin_width
  )

  offspring_size_density <- distribution_kernel(
    mean_matrix = matrices$offspring_size,
    sd = offspring_size_sd,
    sizes = sizes,
    bin_width = bin_width
  )

  # combine adult and seedling survival rates into one
  matrices$survival <- interpolate_survival(
      adult_survival_matrix = matrices$adult_survival,
      seedling_survival_matrix = matrices$seedling_survival,
      sizes = sizes
  )

  # sweep the vector (matrix) components through the matrix (array) components
  # to build the kernel matrices (arrays)

  P_array <- sweep_3(
    growth_density,
    matrices$survival,
    FUN = "*"
  )

  # build fecundity kernel
  # environment:
  matrices$recruits <- matrices$flowering * matrices$seedheads *
    matrices$seeds * matrices$germination

  F_array <- sweep_3(
    offspring_size_density,
    matrices$recruits,
    FUN = "*"
  )

  # combine into a single kernel
  kernel_array <- P_array + F_array

  # intrinsic growth rates across locations
  stable_state <- greta.dynamics::iterate_matrix(kernel_array)
  lambda <- stable_state$lambda

  # observation model parameter
  likelihood_intercept <- normal(0, 10)

  # probability of detection under Poisson sampling

  # form the likelihood
  K <- (lambda - 1) * likelihood_intercept
  Lambda <- log1pe(K)
  p_present <- 1 - exp(-Lambda)
  distribution(occ_train$presence) <- bernoulli(p_present)

  # set up parameters to trace
  names(betas) <- paste0("beta_", names(betas))
  names(size_coefs) <- paste0("gamma_", names(size_coefs))
  params <- c(betas,
              size_coefs,
              list(growth_sd = growth_sd,
                   offspring_size_sd = offspring_size_sd,
                   likelihood_intercept = likelihood_intercept))

  # build the model
  m <- with(params, model())

  # return the model and greta arrays (used for plotting and prediction)
  list(
    model = m,
    parameters = params,
    lambda = lambda
  )

}

# do MCMC on the model
run_mcmc <- function (model_list, verbose = FALSE, ...) {

  timing <- system.time(
    draws <- mcmc(model_list$model, verbose = verbose, ...)
  )

  list(
    draws = draws,
    timing = timing
  )

}

# check
check_draws <- function (draws_list, name) {

  draws <- draws_list$draws

  pdf(file.path("figures", paste0(name, "_draws.pdf")))
  plot(draws)
  dev.off()

  list(
    r_hat = coda:::gelman.diag(draws, multivariate = FALSE),
    n_eff = coda:::effectiveSize(draws),
    summary = summary(draws)
  )

}
