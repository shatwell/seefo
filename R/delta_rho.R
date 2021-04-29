#' @title Fast calculation of mixed layer depth (depth of density threshold)
#'
#' @description Fast calculation of the depth of a set density threshold
#'
#' @param T Temperature data. Must be a mxn matrix with n temperature profiles in the columns
#' @param z The depths corresponding to the temperature profiles. Numeric vector of length m.
#' @param drho The density threshold threshold depths corresponding to the temperature profiles. Numeric vector of length m.
#'
#' @details
#' Calculates mixed layer depths as the depth of when the given density gradient `drho` is first exceeded.
#' Designed to be as fast as possible for large data sets. Does very minimal checks to remove `Inf`s and `NaN`s, so be careful.
#' Choose the right vertical resolution in the temperatures to ensure a more robust result.
#' See Wilson et al. 2020, HESS XXXX
#'
#' @return A vector of mixed layer depths
#'
#' @author Tom Shatwell
#'
#' @seealso \code{\link{h2mix}},  \code{\link{delta_rho}}
#'
#' @examples
#' \dontrun{
#' z <- c(0,2,4,6,10)
#' T <-  matrix(c(15,14,10,8,7,
#' 16,13,12,11,11,
#' 22,19,18,14,13),
#' nrow=length(z))
#' delta_rho(T,z)
#' }
#'
#' @export


delta_rho <- function(T,z,drho=0.1) {
  rho <- rho_water(T) # density
  rho_t <- rho[1,] + drho # abs density threshold
  h <- (z[-length(z)] + (-rho[-length(z),]+rho_t)  / (diff(rho)/diff(z))) * # interpolated MLDs
    (diff(rho>rho_t)==1) # select the one at the threshold
  h[!is.finite(h)]<-NA # get rid of Infs and NaNs
  h1 <- suppressWarnings(apply(h,2,max,na.rm=TRUE))
  h1[!is.finite(h1)] <- NA # filter results
  h1[h1>max(z)] <- NA
  h1[h1==0] <- NA
  return(h1)
}

