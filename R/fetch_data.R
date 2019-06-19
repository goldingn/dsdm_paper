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
  list.files(
    "data/raw/ECOG-00839 Appendix/",
    pattern = ".Rdata$",
    full.names = TRUE
  )

}

# read in the RData object and output a tidied up version with the objects we
# want
tidy_merow <- function (file_path) {

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
    ) -> occurrence

  occurrence

}
