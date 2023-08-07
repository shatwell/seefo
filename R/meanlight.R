#' @title Calculate mean underwater light in a water layer or column
#'
#' @description
#' Calculate mean underwater light in a water column according to the Lambert-Beer law,
#' given light at the surface, the extinction coefficient,
#' and the depth of the layer to calculate the mean for.
#'
#' @param I0 A numeric vector containing light values just beneath the water surface (ie after accounting for backscatter and reflection).
#' @param kd The light extinction coefficient in inverse units of depth (ie 1/m). Can be a single value, which is used for all `I0`, or otherwise must have the same length as `I0`.
#' @param bot The depth of the bottom of the water layer. Typically the mixed layer depth or lake depth for whole lake averages. Can be a single value, which is used for all `I0`, or otherwise must have the same length as `I0`.
#' @param areafun A function that takes elevation (+ve upwards, e.g. m a.s.l. or height above a reference datum like deepest point) as input and returns the lake area at that elevation, typically from `approxfun()`. If `NULL`, the function calculates the integral mean light without considering bathymetry. If `areafun` is given, the mean light is the volume-weighted mean in the layer.
#' @param lev The elevation of the water level (+ve upwards) relative to the same reference datum as `areafun`. Can be a single value which is used for all `I0`, or otherwise must have the same length as `I0`.
#' @param top The depth of the top of the water layer for averaging. Defaults to `0`, or the water surface. Can be a single value, which applies to all `I0`, or otherwise must have the same length as `I0`.
#' @param len The number of internal sublayers to use. Affects precision and is only used when `areafun` is given.
#'
#' @details
#' This function calculates the mean light in a water layer.
#' If you specify `areafun` and `lev`, the function numerically calculates a volume weighted mean,
#' taking the bathymetry into account, as
#' `Imean = integral(I * A * dz) / (integral(A * dz))`, where `I` is irradiance, `A` is area and `z` is depth.
#' If you don't specify `areafun`, the function returns the analytical integral of the mean light as
#' `Imean = I0 / (kd * bot) * (1 - exp(-kd * bot))`.
#' Typically you want to calculate the light in the surface mixed layer, in which case the `top` of the layer is zero and the `bot`tom of the layer is the mixed layer or thermocline depth.
#' You can also calculate mean light in a subsurface layer (e.g. DCM) by specifying `top` and `bot`.
#' The units of the mean light are the same as `I0`. The function is intended to calculate mean light for a time series of `I0`
#' values. You can supply single values for `kd`, `top`, `bot`, and `lev`, which then apply to the whole time series, or otherwise supply
#' a time series of values corresponding to `I0`. See example for how to use `areafun` and `lev` and other things.
#'
#' @return
#' A `vector` containing the vertical mean light values corresponding to `I0`.
#' Returns `NA` if part of the water column is outside the depth / elevation range defined by `areafun`.
#'
#' @author Tom Shatwell
#'
#' @seealso \code{\link{solar2PAR}}, \code{\link{extinction_coeff}}
#'
#' @examples
#' \dontrun{
#' ## example (analytical) mean without bathymetry
#'
#' # solar radiation downwelling (W/m2)
#' incomingsolar <- c(300,240,340,295,400)
#'
#' # convert to PAR just below surface (umol/m2/s)
#' PAR.surf <- seefo::solar2PAR(incomingsolar)
#'
#' # mixed layer depth (m)
#' zmix <- c(8,7,6,5.5,4)
#'
#' # light attenuation coefficient (1/m)
#' extinction <- c(0.55, 0.6, 0.48, 0.45, 0.4)
#'
#' data.frame(incomingsolar, PAR.surf, zmix, extinction)
#'
#' # calculate mean light
#' meanlight(I0 = PAR.surf, kd = extinction, bot = zmix)
#'
#' # calculate mean light using only single values for kd, bot
#' meanlight(I0 = PAR.surf, kd = 0.5, bot = 12.5)
#'
#'
#' ## example (numerical) mean with bathymetry
#' ## bathymetry
#'
#' depth <- seq(80,0,-5)
#'
#' elevation <- seq(340, 420, 5)
#'
#' area <- c(0, 32000, 78000, 132000, 186000, 313000,
#'           429000, 617000, 825000, 1034000, 1283000,
#'           1543000, 1809000, 2125000, 2453000,
#'           2915000, 3481000)
#'
#' data.frame(depth, elevation, area)
#'
#' # create function to return area
#' area_interp <- approxfun(x = elevation, y = area)
#'
#' # water levels (m) on same reference as bathymetry
#' level <- c(418, 416, 413, 409, 412)
#'
#' data.frame(incomingsolar, PAR.surf, zmix, extinction, level)
#'
#' # calculate mean light with bathymetry
#' meanlight(I0 = PAR.surf, kd = extinction, bot = zmix,
#'           areafun = area_interp, lev = level)
#'
#' # with contant water level
#' meanlight(I0 = PAR.surf, kd = extinction, bot = zmix,
#'           areafun = area_interp, lev = 420)
#'
#' # with contant everything except `I0`
#' meanlight(I0 = PAR.surf, kd = 0.5, bot = 12.5,
#'           areafun = area_interp, lev = 420)
#'
#' # mean light in a deeper layer like the thermocline from 11 to 15 m
#' meanlight(I0 = PAR.surf, kd = extinction, top = 11, bot = 15,
#'           areafun = area_interp, lev = 420)
#'
#' # works with moving thermocline
#' meanlight(I0 = PAR.surf, kd = extinction,
#'           top = c(11,11,12,14,13), bot = c(15,16,13,16,15),
#'           areafun = area_interp, lev = 420)
#' }
#'
#' @export



# hypso <- fread("raw/volume-elevation.csv")
# hypso <- fread("raw/bathy.txt")
# hypso[,elev := 338 + max(depth) - depth]
#
# Az <- approxfun(x = hypso$elev, y = hypso$area, method = "linear")


meanlight <- function(I0, kd, bot, areafun=NULL, lev=NULL, top=0, len=101) {
  if(length(kd)==1 && length(I0)>1) {
    kd <- rep(kd, length(I0))
  }
  if(length(bot)==1 && length(I0)>1) {
    bot <- rep(bot, length(I0))
  }
  if(length(top)==1 && length(I0)>1) {
    top <- rep(top, length(I0))
  }
  if(any(c(length(I0)!=length(kd),
           length(I0)!=length(bot),
           length(I0)!=length(top)))) {
    stop("kd, bot and top must be length 1 or same length as I0")
  }
  if(any(top>bot)) {
    stop("bot must be greater than top")
  }
  if(any(top!=0)) {
    I0 <- I0 * exp(-kd * top)
  }

  if(is.null(areafun)) {
    if(!is.null(lev)) {
      warning("areafun is not supplied, ignoring lev")
    }
    if(len!=101) {
      warning("ignoring len, computing analytical mean")
    }
    Im <- I0 / (kd * bot) * (1 - exp(-kd * bot))
  } else {
    if(length(lev)==1 && length(I0) > 1) {
      lev <- rep(lev, length(I0))
    }
    if(length(I0)!=length(lev)) {
      stop("lev must have length either 1 (constant water level) or have the same length as I0, kd and bot")
    }

    # zint <- sapply(bot, function(x) {seq(0,x,length.out=len+1)})
    zint <- mapply(function(x,y) {seq(x,y,length.out=len+1)}, top, bot)
    I0s <- matrix(I0, nrow=len, ncol=length(I0), byrow=TRUE)
    kds <- matrix(kd, nrow=len, ncol=length(kd), byrow=TRUE)
    levs <- matrix(lev, nrow=len, ncol=length(lev), byrow=TRUE)
    dz <- diff(zint) # zint[2,] - zint[1,]
    z <- (zint[-1,] - dz/2)
    I <- I0s * exp(-kds * z)
    A <- matrix(areafun(levs - z), nrow=len, ncol=length(lev))
    Im <- colSums(I * A * dz) / colSums(A * dz)
  }
  return(Im)
}


# meanlight(I0 = PAR.surf, kd = extinction, top = 11, bot = 15,
#           areafun = area_interp, lev = 420)
