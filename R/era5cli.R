#' @title Generate a command to download era5 data with the era5cli tool
#'
#' @description Generates a character string you can paste into a shell
#'
#' @param vars Character vector containing names of ERA5 variables to download. See ERA5 longname in Details.
#' @param startyear Numeric length 1
#' @param endyear Numeric length 1
#' @param prefix A label to begin the downloaded filenames with
#' @param target Not used
#' @param coords Vector (length 2) containing latitude and longitude
#'
#'
#' @details The location and variable name are kept in the attributes of the output.
#'
#' | **ERA5 long name**                                       | **ERA5 shortname** | **Units** |
#' |----------------------------------------------------------|----------------|-------|
#' |surface_solar_radiation_downwards                         | ssrd           | J m**-2    |
#' |surface_solar_radiation_downward_clear_sky                | ssrdc          | J m**-2    |
#' |mean_surface_downward_short_wave_radiation_flux           | msdwswrf       | W m**-2    |
#' |mean_surface_downward_short_wave_radiation_flux_clear_sky | msdwswrfcs     | W m**-2    |
#' |mean_surface_downward_long_wave_radiation_flux            | msdwlwrf       | W m**-2 |
#' |surface_thermal_radiation_downwards                       | strd           | J/m^2 |
#' |2m_temperature                                            | t2m            | K     |
#' |surface_pressure                                          | ??             | Pa    |
#' |2m_dewpoint_temperature                                   | d2m            | K     |
#' |10m_u_component_of_wind                                   | u10            | m/s   |
#' |10m_v_component_of_wind                                   | v10            | m/s   |
#' |mean_total_precipitation_rate                             | mtpr             | kg m**-2 s**-1|
#' |total_precipitation                                       | ??             | m     |
#' |total_cloud_cover                                         | tcc            | 0..1     |
#' |low_cloud_cover                                           | lcc            | 0..1     |
#' |medium_cloud_cover                                        | mcc            | 0..1     |
#' |high_cloud_cover                                          | hcc            | 0..1     |
#' |cloud_base_height                                         | cbh            | m     |
#'
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
