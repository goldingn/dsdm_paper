# this is the main script to replicate the whole analysis. It uses the drake R
# package to cache intermediate objects and skip rerunning them

# use the same RNG seed each time
set.seed(2019-06-19)

# packages ----------------------------------------------------------------

# reproducible workflow
library(drake)

# data munging and plotting
library(raster)
library(RColorBrewer)
library(dplyr)

# read in data ------------------------------------------------------------

source("R/fetch_data.R")

read_data <- drake_plan(
  merow_file = download_merow(),
  occurrence = tidy_merow(merow_file),
  tidy_file = saveRDS(occurrence, "data/clean/occurrence.RDS")
)

make(read_data)



# plot data  --------------------------------------------------------

source("R/plotting.R")

plot_data <- drake_plan(
  occurrence = readRDS(file_in("data/clean/occurrence.RDS")),
  covs_fig = plot_covs(occurrence, file_out("figures/covs.png")),
  occ_fig = plot_occ(occurrence, file_out("figures/occ.png"))
)

make(plot_data)
