# download and unzip the RData object from Merow et al. (hosted on the Ecography data portal)
download_merow <- function () {

  # download and unzip
  url <- "http://www.ecography.org/sites/ecography.org/files/appendix/ecog-00839_appendix.zip"
  zip_file <- tempfile(tmpdir = "data/raw/")
  download.file(url, zip_file)
  unzip(zip_file, exdir = "data/raw")

  # delete the __MACOSX files and the zip file
  unlink("data/raw/__MACOSX", recursive = TRUE)
  file.remove(zip_file) -> .

  # find the path to the RData object and return that
  file_path <- list.files(
    "data/raw/ECOG-00839 Appendix/",
    pattern = ".Rdata$",
    full.names = TRUE
  )

  # read in the RData object and output a tidied up version with the objects we
  # want
  load(file_path)

  # extract the relevant parts of the cape floristic region atlas data

  # minimum July temperature - min07
  # number of summer soil moisture days - smdsum
  # number of winter soil moisture days - smdwin
  # proportion acidic soil - ph1
  # proportion high fertility soil - fert4
  # Mean fire return time - mean.tsf

  cfr.env %>%
    dplyr::select(
      lon,
      lat,
      presence = Presence_5x5,
      abundance = MaxAbun,
      abundance_class = MaxAbunClass,
      min_temp_july = min07,
      fire_return = mean.tsf,
      prop_acidic = ph1,
      prop_fertile = fert4,
      summer_smd = smdsum,
      winter_smd = smdwin,
      mean_annual_precip = map
    )

}

# get landcover and bioclim covariates for the BBS example and
# remove extranous areas
download_bbs_rasters <- function () {

  # get CarolineWren rasters from maxlike package and convert to a raster
  data(carw, envir = environment())
  rl <- lapply(carw.data$raster.data, function(x) {
    m <- matrix(x, nrow = carw.data$dim[1], ncol = carw.data$dim[2],
                byrow = TRUE)
    raster::raster(m)
  })
  rs <- raster::stack(rl)
  extent(rs) <- carw.data$ext
  names(rs) <- names(carw.data$raster.data)
  projection(rs) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  covs1 <- raster::projectRaster(rs, crs = "+init=epsg:4326")

  # bioclim layers
  covs2 <- raster::getData("worldclim", var = "bio", res = 10)
  covs2 <- raster::crop(covs2, c(-128, -66, 21, 52))

  # align them
  covs1 <- resample(covs1, covs2[[1]])
  covs2 <- mask(covs2, covs1[[1]])
  covs <- stack(covs1, covs2)

  covs
  # writeRaster(covs, file = "data/clean/bbs_covs")
}

clean_bbs <- function(file, AOU, years = 2015) {

  bbs <- read.csv(file, stringsAsFactors = FALSE)
  bbs <- bbs[bbs$Year %in% years,]

  bbs$Route_id <- paste(bbs$CountryNum,
                        bbs$StateNum,
                        bbs$Route,
                        sep = "_")

  surveys <- unique(bbs[, c("Route_id", "Year")])
  surveys$PA <- NA

  for (s in 1:nrow(surveys)) {
    idx <- bbs$Route_id == surveys$Route_id[s] &
      bbs$Year == surveys$Year[s]
    tmp <- bbs[idx,]

    if (AOU %in% tmp$AOU) {
      tmp2 <- tmp[tmp$AOU == AOU,][1,]
      count <- sum(tmp2[, paste0("Stop", 1:50)])

      if (count > 0) {
        surveys$PA[s] <- 1
      } else {
        surveys$PA[s] <- 0
      }

    } else {
      surveys$PA[s] <- 0
    }
  }

  surveys
}

# get presence-absence data for the named species from BBS (using a BBS dataset code)
download_bbs <- function (aou = 7360) {

  # download BBS data via FTP
  ftp_path <-
    "ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/50-StopData/1997ToPresent_SurveyWide/Fifty%i.zip"
  record_paths <- sprintf(ftp_path, 1:10)
  routes_path <- "ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/routes.zip"

  # combine these, download and unzip into csvs
  paths <- c(record_paths, routes_path)
  files <- file.path("data/raw", basename(paths))
  future_mapply(download.file, paths, files)
  lapply(files, unzip, exdir = "data/raw")

  # delete zip files
  lapply(files, unlink)

  # get routes in coterminous USA
  routes <- read.csv("data/raw/routes.csv")
  covs <- raster("data/clean/bbs_covs.grd")
  vals <- extract(covs[[1]], routes[, c("Longitude", "Latitude")])
  routes <- routes[!is.na(vals), ]

  # load 50-stop data for the Carolina chickadee
  files <- sprintf("data/raw/fifty%i.csv", 1:10)
  bbs_clean <- lapply(files, clean_bbs, AOU = aou)
  species <- do.call(rbind, bbs_clean)

  routes$Route_id <- paste(routes$CountryNum,
                           routes$StateNum,
                           routes$Route,
                           sep = "_")

  idx <- match(species$Route_id, routes$Route_id)
  species <- cbind(species, routes[idx, c("Longitude", "Latitude")])
  species <- species[!is.na(species$Latitude),]
  species <- species[!is.na(species$Longitude),]

  species
  # saveRDS(species, file = "data/clean/BBS_species.RDS")

}
