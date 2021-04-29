#' @title Calculation of leap years
#'
#' @description Determines whether a year is a leap year
#'
#' @param yr Year (numeric)
#'
#' @return Returns a vector the same length as `yr` containing `0` for non-leap years, and `1` for leap years.
#'
#' @author Tom Shatwell
#'
#' @examples
#' \dontrun{
#' leap(c(2000,2002,2004,2006,2008,2010))
#' }
#'
#' @export

leap <- function(yr) ((yr%%4)==0) - ((yr%%100)==0) + ((yr%%400)==0)
