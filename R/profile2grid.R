#' @title Interpolation of sonde profiles to a grid of fixed depths
#'
#' @description Interpolates vertical profile data (typically from a sonde or profiling system) to defined depths.
#' It uses `data.table` to more quickly and efficiently process large datasets.
#'
#' @param data A data.frame containing the profile data to interpolate, containing columns of at least datetime, depth and values.
#' @param datetime The name of the column containing the datetimes (character).
#' @param depth The name of the column containing the depths (character).
#' @param value The name of the column containing values to be interpolated (character).
#' @param outdepths The interpolation depths, which defines the desired vertical grid (numeric vector). If `NULL`, the function tries to choose suitable depths.
#' @param rule The rule defining how to interpolate outside the data range, as defined in the function \code{approxfun()} function.
#' Copied again: an integer (of length 1 or 2) describing how interpolation is to take place outside the interval [min(x), max(x)].
#' If rule is 1 then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' Use, e.g., rule = 2:1, if the left and right side extrapolation should differ.
#' @param dropNAs Should bad profiles be dropped (`TRUE`) or retained as NAs (`FALSE`) (logical).
#'
#' @details
#' The function takes as input a table in long format, with columns for a timestamp, depth, and measured values.
#' The input data is typically data of multiple profiles from a sonde with one measurement in each row.
#' Each measurement requires a timestamp (`datetime`) that uniquely associates it to a profile.
#' In principle the `datetime` column could contain anything that uniquely identifies each profile, but this is commonly a timestamp.
#' Therefore the times from raw sonde data may need to be `cut()` first. Profiles with less than 5 non-NA values are not interpolated.
#'
#' @return A `matrix` containing the interpolated profile data, nrow = the length of outdepths, and ncol=the number of profiles.
#' The attributes contain some details including the `outdepths` and unique `datetimes`.
#'
#' @author Tom Shatwell
#'
#'
#' @examples
#' \dontrun{
#'
#' depths <- c(0,2,4,6,8,10, 15, 20)
#' times <- c("2005-08-20", "2005-09-02", "2005-09-25")
#'
#' temps <- c(30.8, 30.5, 30.4, 30.3, 29.7, 28.8, 27.8, 27.3,
#'   30.5, 30.0, 29.9, 29.7, 29.1, 29.0, 27.7, 27.3,
#'   30.6, 30.6, 30.5, 30.0, 29.7, 29.2, 27.7, 27.3)
#'
#' dat <- data.frame(
#'   time = rep(times, each=length(depths)),
#'   depth=rep(depths, length(times)),
#'   temp=temps)
#'
#' out <- profile2grid(data=dat, datetime="time", depth="depth",
#'   value="temp")
#'
#' # set interpolation depths manually
#' out <- profile2grid(data=dat, datetime="time", depth="depth",
#'   value="temp", outdepths=seq(0,17,0.5))
#'
#' # works also with any identifier instead of a date
#' dat$time <- rep(c("A","B","C"), each=length(depths))
#'
#' out <- profile2grid(data=dat, datetime="time", depth="depth",
#'   value="temp")
#'
#' attributes(out)
#' }
#'
#' @export

profile2grid <- function(data, datetime, depth, value,
                         outdepths=NULL, rule=c(2,1),
                         dropNAs=FALSE) {
  # require(data.table)
  dat <- data.table::copy(data[,c(datetime,depth,value)])
  names(dat) <- c("dt", "dep", "val")
  data.table::setDT(dat, key = "dt")
  dates <- unique(dat[,dt])
  if(is.null(outdepths)) {
    nodes <- c(0.1,0.2,0.25,0.5,1,2,5,10,20,50,100)
    dz <- diff(dat[,dep])
    dz <- nodes[which.min(abs(nodes-mean(dz[dz>0], na.rm=TRUE)))]
    min_d <- min(dat[,dep], na.rm=TRUE)
    max_d <- max(dat[,dep], na.rm=TRUE)
    outdepths <- seq(max(min_d%/%dz * dz,0), max_d%/%dz * dz, dz)
  }

  # wrapper function to deal with errors
  wrapper <- function(p, fn, ...) {
    out <- NULL
    tryCatch({
      out <- fn(p, ...)
    },
    # warning=function(x) {
    #   out <<- as.character(x)
    # },
    error=function(x) {
      out <<- as.character(x)
    })
    return(out)
  }

  isValid <- function(x) {is.numeric(x)}

  func <- function(y) {
    dat[dt==y, if(sum(!is.na(val))<5) {
      stop("Not enough non-NA values to interpolate")
    } else {
      approx(x = dep, y = val, xout = ..outdepths, rule=rule)$y
    }]
  }

  out <- lapply(X=dates, FUN = wrapper, fn=func)

  retain <- which(sapply(out, isValid)) # which files were successfully read?
  drop <- which(!(1:length(dates)) %in% retain)
  out[drop] <- rep(list(rep(NA, length(outdepths))), length(drop))

  out <- matrix(unlist(out),nrow = length(outdepths))
  attributes(out)$outdepths <- outdepths
  attributes(out)$dt <- dates
  attributes(out)$value <- value
  attributes(out)$NA_profiles <- dates[drop]

  if(length(drop) > 0) { # print profiles that were not successful as warning
    warning(c("Could not interpolate ", length(drop),
              " (of ", length(dates),") profiles:\n",
              paste(dates[drop],"\t")))
    if(dropNAs) {
      out <- out[,retain]
      attirbutes(out)$dt <- dates[retain]
    }
  }

  return(out)
}
