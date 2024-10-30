## code to prepare `retention_eff.rda` dataset goes here

#reading and formating the input data
path <- dirname(rstudioapi::getSourceEditorContext()$path)
ret_eff_data <- read.table(paste0(path,"/retention_eff.csv"), header=T, sep=",", dec=".")
# ret_eff_data <- data[data$year>=start.year & data$year<=end.year,]
# ret_eff_data$date <- lubridate::mdy(ret_eff_data$date)
# ret_eff_data$doy <- lubridate::yday(ret_eff_data$date)
# ret_eff_data$year <- lubridate::year(ret_eff_data$date)
ret_eff_data$date <- as.POSIXct(ret_eff_data$date,tz="UTC")
#ret_eff_data <- ret_eff_data[ret_eff_data$year %in% 2005:2020,]
usethis::use_data(ret_eff_data, overwrite = TRUE)
