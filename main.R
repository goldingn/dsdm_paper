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

# modelling
library(greta.dynamics)

# read in data ------------------------------------------------------------

source("R/fetch_data.R")

read_all_data <- drake_plan(

  # bbs analysis
  bbs_covs = download_bbs_rasters(),
  tidy_covs = writeRaster(bbs_covs, file_out("data/clean/bbs_covs.grd"), overwrite = TRUE),
  bbs_occurrence = download_bbs(),
  tidy_occ = saveRDS(bbs_occurrence, file_out("data/clean/bbs_occ.RDS")),

  # protea analysis
  protea_occurrence = download_merow(),
  tidy_file = saveRDS(protea_occurrence, file_out("data/clean/protea_occurrence.RDS"))

)

# download BBS data chunks in parallel
old_plan <- future::plan()
future::plan(multisession)
make(read_all_data)
future::plan(old_plan)

# fit models ---------------------------------------------------------------

source("R/modelling.R")

fit_bbs_model <- drake_plan(
  bbs_occurrence = readRDS(file_in("data/clean/bbs_occ.RDS")),
  bbs_covs = brick(file_in("data/clean/bbs_covs.grd")),
  model_list = build_bbs_model(bbs_occurrence, bbs_covs),
  draws_list = run_mcmc(model_list),
  predictions = make_bbs_predictions(model_list, draws_list, bbs_covs),
  plots = plot_bbs_maps(predictions, model_list),
  stats = compare_bbs_predictions(predictions, model_list)
)

make(fit_bbs_model)

fit_protea_model <- drake_plan(
  protea_occurrence = readRDS(file_in("data/clean/protea_occurrence.RDS")),
  model_list = build_protea_model(protea_occurrence),
  draws_list = run_mcmc(model_list),
  relationships_fig = plot_protea_relationships(model_list, draws_list, protea_occurrence),
  predictions = make_protea_predictions(model_list, draws_list, protea_occurrence)
)

make(fit_protea_model)

# plot data  --------------------------------------------------------

source("R/plotting.R")

plot_protea_data <- drake_plan(
  protea_occurrence = readRDS(file_in("data/clean/protea_occurrence.RDS")),
  protea_covs_fig = plot_covs(protea_occurrence, file_out("figures/protea_covs.png")),
  protea_occ_fig = plot_occ(protea_occurrence, file_out("figures/protea_occ.png"))
)

make(plot_protea_data)


