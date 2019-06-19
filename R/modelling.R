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
get_env_effect <- function (formula, data) {
  design <- model.matrix(formula, data)
  beta <- normal(0, 10, dim = ncol(design))
  env_effect <- design %*% beta
  list(beta = beta, env_effect = env_effect)
}

# if the relationship is size-dependent, get a coefficient for that, and return
# along with the (row-vector) offsets for all size classes
get_size_effect <- function (size_dependent, sizes) {
  coef <- 0
  if (size_dependent) {
    coef <- normal(0, 10)
  }
  list(coef = coef, size_effect = coef * t(sizes))
}

apply_link <- function (link_function, matrix) {
  link_function(matrix)
}

# blank function for no link/transformation
none <- function (x) x

# construct and integrate over (each slice representing a kernel matrix for a
# given slice, in) a kernel array for a transition with a normal distribution.
# given a matrix of means at sites and sizes (n_sites x matrix_dim), a scalar
# standard deviation for the distribution, a vector of sizes to integrate over
# (of length matrix_dim), and the integration bin size, return a kernel array
# (n_sites x matrix_dim x matrix_dim)
integrate_kernel <- function (mean_matrix, sd, sizes, bin_width) {
  n_sites <- nrow(mean_matrix)
  matrix_dim <- ncol(mean_matrix)
  var <- sd ^ 2
  constant <- 1 / (sqrt(2 * pi * var))
  # diff is a matrix with dimensions (n_sites * matrix_dim) x matrix_dim. We'll
  # reshape it last (need to add tests for these operations as they are easy to
  # cock up!)
  diff <- kronecker(sizes, mean_matrix, FUN = "-")
  density <- bin_width * constant * exp((diff ^ 2) / (2 * var))
  dim(density) <- c(n_sites, matrix_dim, matrix_dim)
  density
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

# build and return a greta model for the DSDM
build_model <- function (occurrence) {

  # subset occurrence data to wheere we observe presence/absence, and subset to
  # a few observations (for now)
  occurrence %>%
    filter(!is.na(presence)) %>%
    sample_n(3) -> occ_train

  # formulae for the *environmental components* of the vital rate regressions
  formulae <- list(
    growth = ~ quad(prop_fertile) + quad(winter_smd) + quad(prop_acidic) + min_temp_july + quad(winter_smd),
    flowering = ~ 1,
    seedheads = ~ min_temp_july + quad(prop_acidic),
    seeds = ~ 1,
    offspring_size = ~ prop_acidic + quad(prop_fertile) + min_temp_july + quad(summer_smd),
    adult_survival = ~ min_temp_july + mean_annual_precip,
    seedling_survival = ~ min_temp_july + mean_annual_precip
  )

  # for each formula, get the vector of coefficients and environmental effects
  env_effects_list <- lapply(formulae, get_env_effect, occ_train)
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
    growth = FALSE,
    flowering = TRUE,
    seedheads = TRUE,
    seeds = FALSE,
    offspring_size = FALSE,
    adult_survival = TRUE,
    seedling_survival = TRUE
  )

  size_effects_list <- lapply(size_dependent, get_size_effect, sizes)
  size_coefs <- lapply(size_effects_list, `[[`, "coef")
  size_effects <- lapply(size_effects_list, `[[`, "size_effect")


  # There's a deterministic effect of size on growth; we need to add the origin
  # size onto the (environmentally-driven) growth increment to get the
  # destination size. most efficient to do that here.
  size_effects$growth <- size_effects$growth + sizes

  # combine the size and environment effects (get the outer sum of the two) to
  # get matrices on the link scale
  link_matrices <- mapply(kronecker,
                          size_effects, env_effects,
                          MoreArgs = list(FUN = "+"),
                          SIMPLIFY = FALSE)

  # transform to parameters to the response scale to build kernel array
  link_functions <- list(
    growth = none,
    flowering = ilogit,
    seedheads = exp,
    seeds = exp,
    offspring_size = none,
    adult_survival = ilogit,
    seedling_survival = ilogit
  )

  matrices <- mapply(apply_link,
                     link_functions, link_matrices,
                     SIMPLIFY = FALSE)

  # get coefficients for variation in growth and offspring size
  growth_sd <- normal(0, 10, truncation = c(0, Inf))
  offspring_size_sd <- normal(0, 10, truncation = c(0, Inf))

  # get growth and offspring size transition probabilities (integrating over
  # discretized normal distributions) to return kernel arrays
  growth_density <- integrate_kernel(
    mean_matrix = matrices$growth,
    sd = growth_sd,
    sizes = sizes,
    bin_width = bin_width
  )

  offspring_size_density <- integrate_kernel(
    mean_matrix = matrices$offspring_size,
    sd = offspring_size_sd,
    sizes = sizes,
    bin_width = bin_width
  )

  # need to reshape these as matrices, multiply by other kernel components,
  # sweep-divide by density sums to do integration, then reshape back again
  # phew!


  a <- matrix(1:6, 3, 2)
  b <- (1:2) / 10
  c <- kronecker(b, a, FUN = "+")
  dim(c) <- c(3, 2, 2)
  c[1, , ]



  # integrate over growth function to get growth array G



  # integrate over offspring size function to get offspring size array O

  # apply survival probabilities S to G, to get P (growth /survival subkernel)



  # build the kernel array across observation locations





  # get the intrinsic growth rate for each location
  greta.dynamics::iterate_matrix(kernel_array)

  # form the likelihood
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
# relationships
plot_relationships <- function(draws, model_list, occurrence) {

}

# use the model and MCMC parameter samples to make posterior predictions of
# the probability of presence in new places
make_predictions <- function(draws, model_list, occurrence) {

}

