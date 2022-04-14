#' @title Generate a command to download era5 data with the era5cli tool
#'
#' @description Generates a character string you can paste into a shell
#'
#' @param vars Character vector containing names of ERA5 variables to download
#' @param startyear Numeric length 1
#' @param endyear Numeric length 1
#' @param prefix A label to begin the downloaded filenames with
#' @param target Not used
#' @param coords Vector (length 2) containing latitude and longitude
#'
#'
#' @details The location and variable name are kept in the attributes of the output.
#' @return A character vector length 1 with the command and arguments to use with the era5cli tool
#' @author Tom Shatwell
#'
#' @examples
#' \dontrun{
#' vars <- c(
#' "mean_surface_downward_long_wave_radiation_flux", # W m**-2
#' "mean_surface_downward_short_wave_radiation_flux", # W m**-2
#' "total_cloud_cover", # 0..1
#' "2m_temperature", # K
#' "surface_pressure", # Pa
#' "2m_dewpoint_temperature", # K
#' "10m_u_component_of_wind", # m s**-1
#' "10m_v_component_of_wind" # m s**-1
#' )
#'  era5cli(vars=vars,startyear=2013,endyear=2019,
#'  prefix = "era_kilpisjarvi",
#'  target = "../nc_input/era5/",
#'  coords = c(69.03, 20.8))
#' }
#'
#' @export

era5cli <- function(vars, startyear, endyear,
                    prefix, target, coords) {

  coords_lwr <- round(coords-0.006,2) # round(coords * 4)/4 - 0.125
  coords_upr <- round(coords+0.006,2)
  area <- c(coords_upr[1], coords_lwr[2], coords_lwr[1], coords_upr[2])

  cmd <- "era5cli hourly"

  args <- paste0("--variables ", paste(vars, collapse=" "),
                 " --startyear ", startyear, " --endyear ",
                 endyear,
                 " --area ", paste(area, collapse=" "),
                 " --merge --outputprefix ", prefix)
  return(paste(cmd,args))
}
