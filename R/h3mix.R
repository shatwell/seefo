#' @title Robust calculation of mixed layer depth
#'
#' @description Robust calculation of the mixed layer depth and thermocline depth,
#' more suited for manual or irregular
#' temperature data.
#'
#' @param T Temperature profile
#' @param z The depths corresponding to the temperature profile
#' @param plot Should a plot of the profile showing the estimates be shown? (logical)
#' @param thresh The threshold for the difference between top and bottom temperature to infer stratification
#' @param therm.dz Representative thickness of the thermocline layer. 1 m works well from experience.
#' @param min.hmix I can't really remember what this is for.
#' @param min.gradient Not used
#' @param ... Arguments passed to `plot()`, when `plot=TRUE`
#'
#' @details
#' This function estimates the mixed depth and the thermocline depth.
#' It uses regressions to locate the kink in the profile at the bottom of
#' the mixed layer as the mixed layer, and the maximum gradient as the thermocline.
#' It is better suited to irregular, esp manually measured profile data.
#' It first finds the minimum curvature in the
#' profile (maximum for winter stratification) as the initial guess of the
#' border of the surface mixed
#' layer. Then it finds the thermocline as the depth of maximum T-gradient. It
#' estimates the mixed layer depth as the depth where the regression line through
#' the surface layer temperatures intersects with the regression line through the
#' thermocline temperatures. It performs some other checks, like whether
#' stratification exists surface-bottom temperature difference > `thresh`, and whether the thermocline
#' is above the surface layer, or whether the lake is completely isothermal and
#' there is no intersection (and returns `NA`). If mixed, it assumes the mixed
#' layer depth is the maximum depth. It also plots the profile if desired
#' (`plot=TRUE`), marking the surface layer and inflection point in red and the
#' thermocline in blue.
#' Surface temperature is estimated as the mean of the upper 2 m, or the top two measurements,
#' bottom temperature is estimated as the mean of the deepest 2 m layer.
#'
#'
#' @return A vector of length 2 with the mixed layer depth `hmix` and the thermocline depth `htherm`
#'
#' @author Tom Shatwell
#'
#' @seealso \code{\link{hmix}}, \code{\link{h2mix}},  \code{\link{delta_rho}}
#'
#' @examples
#' \dontrun{
#' T <- c(25.44,25.46,25.48,25.46,25.44,24.36,20.48,18.09,15.96,13.95,11.67,
#' 10.62,9.82,9.51,9.22,9.01,8.81,8.68,8.58,8.39,8.29,8.14,8.08,8.03,8.01,
#' 8.05,7.94,7.93,7.85)
#' z <- 0:28
#' h3mix(T,z,plot=TRUE)
#' }
#'
#' @export


h3mix <- function(T, z, plot=FALSE, thresh = 1,
                  therm.dz = 1, min.hmix = 1.5, min.gradient = 0.3, ...) {
  if(sum(!is.na(T))>3) {
    dTdz <- diff(T) / diff(z)
    d2Tdz2 <- diff(T,1,2) / diff(z,2)^2
    i1 <- which.min(dTdz) # index of strongest gradient (-ve for positive strat, +ve for winter strat)
    i2 <- which.min(d2Tdz2) # index of minimum curvature (bot epilimnion)
    while(z[i2]<=min.hmix) {
      d2Tdz2[i2] <- NA
      i2 <- which.min(d2Tdz2)
    }
    while(i1 <= i2) {
      dTdz[i1] <- NA
      i1 <- which.min(dTdz)
    }
    i3 <- which.max(d2Tdz2) # index of maximum curvature (top hypolimnion)
    if(min(z)<2) {
      Ts <- mean(T[z<=2], na.rm=TRUE)
    } else {
      Ts <- mean(T[1:i1], na.rm=TRUE)
      warning("Warning: Surface temp estimated as temps above greatest curvature")
    }
    Tb <- mean(T[z>=(max(z)-2)], na.rm=TRUE)
    if(Ts < 4) {
      i1 <- which.max(dTdz) # index of maximum gradient
      i2 <- which.max(d2Tdz2) # max curvature for winter stratification
      i3 <- which.min(d2Tdz2) # min curvature for winter stratification
    }
    h1 <- mean(c(z[i1], z[i1+1])) # depth of maximum gradient
    h2 <- z[i2+1] # depth of inflection

    thermo <- data.frame(z=z[z>=z[i1]-therm.dz &
                               z <= z[i1+1]+therm.dz],
                         T = T[z>=z[i1]-therm.dz &
                                 z <= z[i1+1]+therm.dz])

    regs <- lm(T[1:i2]~z[1:i2]) # regression thru surface temps
    regt <- lm(T~z, thermo) # regression thru thermocline
    dTdz.s <- coef(regs)[2] # T gradient in surface layer
    if(is.na(coef(regs)[2])) dTdz.s <- 0

    h <- (coef(regt)[1] - coef(regs)[1]) / # depth of intersection of thermocline and surface layer regression lines
      (dTdz.s - coef(regt)[2])

    if(!is.na(h) & h > z[i1+1]) h <- NA # ignore if intersection of thermocline and surface layer is below thermocline
    if(i1 < i2) h <- NA # ignore if minimum curvature is below the thermocline
    if(abs(Ts-Tb) < thresh) h <- max(z) # assume mixed to max(z) if Ts-Tb < thresh
    if(!is.na(h) & h < 0) h <- NA
    if(plot){
      plot(T, z,
           ylim=c(max(z),0), type="n", ...)
      abline(h=0:(max(z)), col="grey")
      abline(h=h)
      # abline(h=h2, lty=2)
      # if(!c(coef(regs)[2] %in% c(0,NA) | !coef(regt)[2] %in% c(0,NA))) {
      #   abline(-coef(regs)[1]/coef(regs)[2], 1/coef(regs)[2], lty=2, col="red")
      #   abline(-coef(regt)[1]/coef(regt)[2], 1/coef(regt)[2], lty=2, col="blue")
      # }
      # points(coef(regs)[2] * h + coef(regs)[1], h, pch=16, col="blue", cex=1.5)
      lines(T, z, type="o")
      lines(T[1:i2], z[1:i2], col="red", lwd=2)
      lines(T[c(i1,i1+1)], z[c(i1,i1+1)], col="blue", lwd=2)
      points(T[i2+1], z[i2+1], col="red",pch=16)
      box()
    }
  } else {
    h <- NA
    warning("Warning:need more than 3 temperature values")
  }
  out <- c(h, mean(z[c(i1,i1+1)]))
  names(out) <- c("hmix","htherm")
  return(out)
}






