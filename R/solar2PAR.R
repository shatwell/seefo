#' @title Convert solar radiation into PAR
#'
#' @description Converts solar (shortwave) radiation in W/m2 above the water surface into PAR (photosynthetically active radiation in umol / m2 / s) just below the water surface
#'
#'
#' @param solar Solar radiation in W m^-2 (numeric)
#' @param parfrac The fraction of solar radiation that is PAR (400 - 700 nm) (numeric)
#' @param reflect The loss due to reflection and backscatter at the surface of the water (numeric)
#' @param joule2umol Conversion factor for solar radiation to PAR with units umol per Joule (numeric)
#'
#'
#' @details This function helps to convert solar radiation usually available from weather stations into PAR.
#' Solar radiation should be in W/m2. The assumptions are that 50% of solar radiation is PAR (400 - 700 nm),
#' and 10% of that is lost at the water surface due to reflection.
#' The conversion to micromoles of photons (400-700) depends on the spectrum but the default value is typical for sunlight.
#' Units of output are umol photons / m2 / s unless you use a different conversion factor `joule2mol`.
#' The conversion is `PAR = solar * parfrac * (1 - reflect) * joule2umol`.
#'
#' @return A numeric vector the same length as `solar`.
#' @author Tom Shatwell
#'
#' @examples
#'
#' \dontrun{
#' solar2PAR(c(100, 200, 300))
#' }
#'
#' @export
#
solar2PAR <- function(solar, parfrac=0.5, reflect=0.1, joule2umol=4.56) {
  solar * parfrac * (1 - reflect) * joule2umol
}
