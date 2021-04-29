#' @title Fast calculation of mixed layer depth (lower thermocline boundary)
#'
#' @description Fast calculation of lower thermocline boundary
#'
#' @param T Temperature data. Must be a mxn matrix with n temperature profiles in the columns
#' @param z The depths corresponding to the temperature profiles. Numeric vector of length m.
#'
#' @details
#' Calculates mixed layer depths as the depth of the maximum curvature of the temperature profile (d2T/dz2).
#' Designed to be as fast as possible for large data sets, so doesn't do any checks at all (be careful).
#' Choose the right vertical resolution in the temperatures to ensure a more robust result.
#'
#' @return A vector of mixed layer depths
#'
#' @author Georgiy Kirillin and Tom Shatwell
#'
#' @seealso \code{\link{hmix}},  \code{\link{delta_rho}}
#'
#' @examples
#' \dontrun{
#' z <- c(0,2,4,6,10)
#' T <-  matrix(c(15,14,10,8,7,
#' 16,13,12,11,11,
#' 22,19,18,14,13),
#' nrow=length(z))
#' h2mix(T,z)
#' }
#'
#' @export

h2mix <- function(T,z) {
  d2T<-diff(T,differences = 2)
  dz2<-diff(z)^2
  d2Tdz2 <- d2T/dz2[-length(dz2)]
  return(z[apply(d2Tdz2,2,which.max)+1])
}
