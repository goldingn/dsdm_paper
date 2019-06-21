# load the parameter estimates from Merow et al. and format them as the
# parameters of normal priors
prior <- function (model = c("seedling_survival", "adult_survival", "growth",
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

  # # for intercept-only models, don't multiply by the intercept column, it's mode
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
  # diff is created as a matrix with dimensions (n_sites * matrix_dim) x matrix_dim, so we
  # reshape it (to n_sites x matrix_dim x matrix_dim)
  diff <- kronecker(sizes, mean_matrix, FUN = "-")
  dim(diff) <- c(n_sites, matrix_dim, matrix_dim)
  log_density <- log(bin_width) + log_constant - (diff ^ 2) / (2 * var)
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
sweep_3 <- function(A, B, FUN = c('-', '+', '/', '*')) {
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
build_model <- function (occurrence, informative_priors = FALSE) {

  # subset occurrence data to wheere we observe presence/absence, and subset to
  # a few observations (for now)
  occurrence %>%
    filter(!is.na(presence)) %>%
    sample_n(3) -> occ_train

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
    prior,
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
    prior,
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

  # get normalisations for the columns in the growth kernel, so we can make them
  # sum to 1 (to complete the integration)
  growth_normalisation <- apply(growth_density, c(1, 3), "sum")

  # Merow has:
  # (Ss*P)/Ps = (matrices$survival * growth_density) / growth_normalisation
  # the following (P * (Ss / Ps)) is the same but in fewer FLOPs.
  # sweep_3 is a custom function to do an array-matrix sweep on the third
  # dimension (columns of each matrix).
  P_array <- sweep_3(
    growth_density,
    matrices$survival / growth_normalisation,
    FUN = "*"
  )

  # get normalisations for the columns in the offspring size kernel, so we can
  # make them sum to 1 (to complete the integration)
  offspring_size_normalisation <- apply(offspring_size_density, c(1, 3), "sum")

  # build fecundity kernel
  # environment:
  matrices$recruits <- matrices$flowering * matrices$seedheads *
    matrices$seeds * matrices$germination

  F_array <- sweep_3(
    offspring_size_density,
    matrices$recruits / offspring_size_normalisation,
    FUN = "*"
  )

  # combine into a single kernel
  kernel_array <- P_array + F_array

  # get the intrinsic growth rate for each location
  iterations <- greta.dynamics::iterate_matrix(kernel_array)

  # form the likelihood
  log_lambda <- log(iterations$lambda)
  log_rate <- log_lambda + intercept
  p_present <- icloglog(log_rate)
  distribution(occurrence$presence) <- binomial(p_present)

  # build the model
  m <- model(intercept, beta_growth, beta_fecundity, beta_survival)

  # return the model and greta arrays (used for plotting and prediction)
  list(
    model = m,
    intercept = intercept,
    beta_growth = beta_growth,
    beta_fecundity = beta_fecundity,
    beta_survival = beta_survival
  )

}

# do MCMC on the model until converged
run_mcmc <- function (model_list) {

}

# use the model and MCMC parameter samples to visualise the fitted vital rate
# relationships (using calculate and relevant design matrices)
plot_relationships <- function(draws, model_list, occurrence) {

}

# use the model and MCMC parameter samples to make posterior predictions of
# the probability of presence in new places (using calculate and relevant design matrices)
make_predictions <- function(draws, model_list, occurrence) {

}

