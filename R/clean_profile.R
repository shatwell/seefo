#' @title Clean up a profile
#'
#' @description Strips unwanted values from a profile, usually read from raw data.
#'
#' @param data The profile data (`data.frame``) containing a column named either "depth" or "press"
#' and a column named "dt" containing the time stamp of the measurement
#' @param remove_lower Numeric value of the thickness of the bottom layer to truncate
#'
#' @details
#' This function removes rows in a profile that:
#' * have depth less than zero (ie in air),
#' * were taken when the sonde was pulled back up
#' * are within a certain distance of the maximum depth, given by `remove_lower`.
#' The idea of `remove_lower` is to remove the lowest values when the sonde hits the sediment.
#'
#' @return The cleaned profile data
#'
#' @author Tom Shatwell
#'
#' @seealso \code{\link{get_profile_info}}
#'
#' @examples
#' \dontrun{
#' pr_loc <- system.file("extdata", package="seefo")
#' pr_name <- file.path(pr_loc,  "profiles","20170801","YT1_20170801_1.CTM644")
#' pr1 <- read_ctm(pr_name, clean=FALSE)
#' tail(pr1)
#'
#' clean_profile(pr1)
#' tail(pr1)
#' }
#'
#' @export

clean_profile <- function(data,remove_lower=0.5) {
  depthcol <- names(data) %in% c("depth","press")
  # check increasing chronologically
  data <- data[order(data$dt),]
  # select data where depth is positive and increasing only
  # this removes measurements in air, and also when pulling the sonde up
  deepest <- rep(NA,nrow(data))
  for(i in 1:nrow(data)) {
    deepest[i] <- max(data[1:i,depthcol])
  }
  keep <- data[,depthcol] >= 0 &
   data[,depthcol] >= deepest  &
    data[,depthcol] <= (max(data[,depthcol]) - remove_lower)
  return(data[keep,])
}
