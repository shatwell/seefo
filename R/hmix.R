#' @title Fast calculation of mixed layer depth (upper thermocline boundary)
#'
#' @description Fast calculation of the upper thermocline boundary
#'
#' @param T Temperature data. Must be a mxn matrix with n temperature profiles in the columns
#' @param z The depths corresponding to the temperature profiles. Numeric vector of length m.
#'
#' @details
#' Calculates mixed layer depths as the depth of the maximum temperature gradient (dT/dz).
#' Designed to be as fast as possible for very large data sets, like model output.
#' It doesn't do any checks at all, and only works with regular data, so be careful.
#' Choose the right vertical resolution in the temperatures to ensure a more robust result.
#' Use h3mix or the rLakeAnalyzer functions for more robust solutions which are
#' orders of magnitude slower.
#'
#' @return A vector of mixed layer depths
#'
#' @author Georgiy Kirillin and Tom Shatwell
#'
#' @seealso \code{\link{h2mix}},  \code{\link{h3mix}},  \code{\link{delta_rho}}
#'
#' @examples
#' \dontrun{
#' z <- c(0,2,4,6,10)
#' T <-  matrix(c(15,14,10,8,7,
#' 16,13,12,11,11,
#' 22,19,18,14,13),
#' nrow=length(z))
#' hmix(T,z)
#' }
#'
#' @export

hmix <- function(T,z) z[apply(abs(diff(T)/diff(z)),2,which.max)]
