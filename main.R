# this is the main script to replicate the whole analysis. It uses the drake R
# package to cache intermediate objects and skip rerunning them

# use the same RNG seed each time
set.seed(2019-06-19)

# packages ----------------------------------------------------------------

# reproducible workflow
library(future.apply)
library(drake)

# data munging and plotting
library(raster)
library(RColorBrewer)
library(dplyr)
library(maxlike)
library(pROC)
library(truncdist)

# modelling
library(greta.dynamics)

# read in data ------------------------------------------------------------

source("R/fetch_data.R")

read_all_data <- drake_plan(

  download_bbs_rasters(file_out("data/clean/bbs_covs.grd")),
  download_bbs(file_out("data/clean/bbs_occ.RDS")),
  download_maps_locations(file_out("data/clean/maps_stations.RDS")),
  download_protea(file_out("data/clean/protea_occurrence.RDS")),

)

# download BBS data chunks in parallel
old_plan <- future::plan()
future::plan(multisession)
make(read_all_data)
future::plan(old_plan)

# fit models ---------------------------------------------------------------

source("R/modelling.R")

fit_bbs_model <- drake_plan(

  # read in data
  bbs_occurrence = readRDS(file_in("data/clean/bbs_occ.RDS")),
  bbs_covs = brick(file_in("data/clean/bbs_covs.grd")),
  maps_coords = readRDS(file_in("data/clean/maps_stations.RDS")),

  # prepare the train/test split and design matrices
  bbs_formula = ~ 1 + pcMix + pcDec + pcCon + bio1 + bio6 + bio12,
  bbs_data_list = prep_data(bbs_formula, bbs_occurrence, bbs_covs, maps_coords),

  # model fitting
  bbs_model_list = build_bbs_model(bbs_data_list),
  bbs_draws_list = run_mcmc(bbs_model_list),
  bbs_draws_summary = check_draws(bbs_draws_list, "bbs"),

  # also fit the analytic version
  bbs_analytic_model_list = build_bbs_model(bbs_data_list, analytic = TRUE),
  bbs_analytic_draws_list = run_mcmc(bbs_analytic_model_list),
  bbs_analytic_draws_summary = check_draws(bbs_analytic_draws_list, "bbs"),

  # prediction
  bbs_predictions = make_bbs_predictions(bbs_model_list,
                                         bbs_draws_list,
                                         bbs_data_list),
  bbs_plots = plot_bbs_maps(bbs_predictions, bbs_data_list),
  bbs_stats = compare_bbs_predictions(bbs_predictions, bbs_data_list),

)

make(fit_bbs_model)

fit_protea_model <- drake_plan(

  protea_occurrence = readRDS(file_in("data/clean/protea_occurrence.RDS")),

  informative_model_list = build_protea_model(
    protea_occurrence,
    informative_priors = TRUE
  ),
  informative_draws_list = run_mcmc(informative_model_list, verbose = TRUE, sample = hmc(Lmin = 10, Lmax = 15)),
  informative_draws_summary = check_draws(informative_draws_list, "protea_informative"),

  naive_model_list = build_protea_model(
    protea_occurrence,
    informative_priors = FALSE
  ),
  naive_draws_list = run_mcmc(naive_model_list),
  naive_draws_summary = check_draws(naive_draws_list, "protea_naive"),

)

make(fit_protea_model)

# plot data  --------------------------------------------------------

source("R/plotting.R")

plot_protea_data <- drake_plan(
  protea_occurrence = readRDS(file_in("data/clean/protea_occurrence.RDS")),
  protea_covs_fig = plot_covs(protea_occurrence, file_out("figures/protea_covs.png")),
  protea_occ_fig = plot_occ(protea_occurrence, file_out("figures/protea_occ.png")),
  protea_estimates = readRDS(file_in("data/clean/protea_estimates.RDS")),
  protea_params_fig = plot_params(protea_estimates, file_out("figures/protea_occ.png"))
)

make(plot_protea_data)


