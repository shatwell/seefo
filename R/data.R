#' Surface and bottom temperatures and ice thickness
#'
#' A dataset containing the surface and bottom temperatures as well as ice thickness from
#' historical model simulations of Lake Stechlin with the model FLake version 2.
#'
#' @format A data frame with 3653 rows and 4 variables:
#' \describe{
#'   \item{date}{date in POSIXct format}
#'   \item{Ts}{The surface temperature in degrees Celsius}
#'   \item{Tb}{The bottom temperature in degrees Celsius, depth ~23 m}
#'   \item{H_ice}{The ice thickness in meters}
#'   ...
#' }
#' @source Tom Shatwell
"Ts_Tb_ice"

#' Inflow and outflow concentration time series
#'
#' A dataset containing the inflow and outflow nutrient concentrations
#'
#' @format A data frame with 5004 rows and 7 variables:
#' \describe{
#'   \item{date}{date in POSIXct format}
#'   \item{year}{Year corresponding to date}
#'   \item{doy}{Day of year}
#'   \item{tributary}{Name of the tributary}
#'   \item{variable}{Name nutrient species}
#'   \item{value}{Measured concentration}
#'   \item{in_outlet}{Whether the tributary is an inflow or outflow}
#'   ...
#' }
#' @source Taynara Fernadez
"ret_eff_data"
