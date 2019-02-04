import::from(glue, glue)
import::from(R.utils, gunzip)

years <- seq(1981, 2000)
year_urls <- glue(
  "ftp://ftp.glcf.umd.edu/glcf/GLOPEM/stow/stow/ft0001/GLCF.TSM.B7-001.00.Summed-Annual-Global/Summed_Annual/",
  "{years}_npp_latlon/{years}_npp_latlon.tif.gz"
)

my_download <- function(url, target_file, ...) {
  if (file.exists(target_file)) {
    message("File ", target_file, " already exists.")
    return(NULL)
  }
  download.file(url, target_file, ...)
}

target_dir <- "glcf_data"
dir.create(target_dir, showWarnings = FALSE)
files <- file.path(target_dir, glue("{years}_npp_latlon.tif.gz"))
.dl <- Map(my_download, year_urls, files)

all_raw_files <- list.files("glcf_data", pattern = "\\.gz$", full.names = TRUE)
.gz <- Map(gunzip, all_raw_files)
