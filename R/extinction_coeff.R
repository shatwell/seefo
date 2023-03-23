#' @title Calculate the extinction coefficient from PAR
#'
#' @description Calculates the light extinction (attenuation) coefficient from measurements of PAR data according to the Lambert-Beer law
#'
#'
#' @param PAR The PAR measurements (numeric)
#' @param z Depths of PAR measurements (numeric)
#' @param plot_profile Should the data be plotted, showing the fitted extinction coefficient (logical)? Only for profiles with more than 2 measurements.
#'
#'
#' @details Extinction is best calculated from pairs of simultaneous PAR measurements at different depths with 2
#' different sensors. The function will assume this is the case if the `length(PAR)` and `length(z)` = 2.
#' In this case the function will calculate the extinction analytically from the Lambert-Beer law.
#' Sometimes only a (non-simultaneous) profile of PAR measurmentns is available, and the function will assume this is the case if `length(PAR)` and `length(z)` > 2.
#' In this case it will calculate the extinction coefficient as the slope of the regression through the logged PAR measurements.
#'
#' @return A numeric length 1.
#' @author Tom Shatwell
#'
#' @examples
#'
#' \dontrun{
#' # pair of measurements with 2 sensors
#'
#' depth <- c(0.75, 1.25)
#' light <- c(1263, 957)
#'
#' extinction_coeff(PAR=light, z=depth)
#'
#' # with a profile of measurements
#'
#' depth <- c(0.5, 1, 1.5, 2, 2.5, 3)
#' light <- c(769, 532, 370, 211, 158, 95)
#'
#' extinction_coeff(light,depth, plot_profile=TRUE)
#'
#' }
#'
#' @export
#
extinction_coeff <- function(PAR, z, plot_profile=FALSE) {
  if(length(PAR) < 2) {
    stop("Need at least 2 measurements")
  }
  if(length(z)!=length(PAR)) {
    stop("Length of z must equal the length of PAR")
  }
  if(any(z<0)) {
    warning("Negative depths were switched to positive numbers")
  }
  z <- abs(z)
  ii <- order(z)
  z <- z[ii]
  PAR <- PAR[ii]
  if(length(PAR) == 2) {
    out <- log(PAR[1]/PAR[2])/(z[2]-z[1])
  }
  if(length(PAR)>2) {
    reg <- lm(log(PAR)~z)
    out <- as.numeric(-coef(reg)[2])
    if(plot_profile) {
      plot(z~PAR, ylim=c(max(z),min(z)),log="x")
      zz <- seq(min(z),max(z), length.out=50)
      lines(exp(predict(object = reg, newdata = list(z=zz))), zz)
    }
  }
  return(out)
}
