#' @title Calculates density of freshwater
#'
#' @description Calculates the density of fresh water.
#'
#' @param t Temperature in degrees Celsius (numeric)
#'
#' @details
#' Calculates the density from temperature (kg m^-3) using the formula of Millero & Poisson (1981).
#' this is the method stated in the isimip 2b protocol in July 2019
#'
#' @return numeric vector the same length as `t`
#'
#' @author Georgiy Kirillin and Tom Shatwell
#'
#' @seealso \code{\link{hmix}},  \code{\link{delta_rho}}
#'
#' @examples
#' \dontrun{
#' rho_water(0:25)
#' }
#'
#' @export

rho_water <- function(t) {
  999.842594 + (6.793952e-2 * t) - (9.095290e-3 * t^2) +
    (1.001685e-4 * t^3) - (1.120083e-6 * t^4) + (6.536336e-9 * t^5)
}

