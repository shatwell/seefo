#' @title Calculate humidity from dewpoint temperature
#'
#' @description Calculate humidity from dewpoint temperature after Clausius-Claperon equations.
#'
#' @param Td Dewpoint temperature in Kelvin (numeric).
#' @param E0 Constant (numeric).
#' @param L_Rv Constant (numeric).
#' @param T0 Freezing point in Kelvin (numeric).
#'
#' @details
#' Useful for converting between units for hydrodynamic modelling.
#' To calculate saturation humidity, substitute Td = air temperature instead of dewpoint temperature.
#'
#'
#' @return The humidity in mbar.
#'
#' @references
#' Lawrence, Mark G., 2005: The relationship between relative humidity and the dewpoint temperature in moist air: A simple conversion and applications. Bull. Amer. Meteor. Soc., 86, 225-233. doi: http;//dx.doi.org/10.1175/BAMS-86-2-225)
#'
#' @author Tom Shatwell
#'
#' @examples
#' \dontrun{
#'  dewpoint_2_humidity(293)
#' }
#'
#' @export

dewpoint_2_humidity <- function(Td, E0=610.94, L_Rv = 5423, T0=273.15) {
  E0 * exp(L_Rv * ( (1/T0) - (1/Td) ) )/100 # humidity in mbar
}
