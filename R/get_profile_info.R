#' @title Extracts information from a sonde filename
#'
#' @description Extracts information from a sonde raw data filename, like the measuring location, date, and sonde number.
#'
#' @param filename Filename to examine (character), including the file extension.
#'
#' @details
#' This function parses the filename and extracts information according to the limnophsycs naming convention.
#' It assumes that the filename has the structure `"LOC_DATE_PROFILENUMBER_OPTIONALSTUFF.SONDENUM"`.
#' It separates the name into chunks separated by `"_"`. Whole paths can be given too. In this case the path discarded from the filename.
#'
#' @return A character vector with elements containing the location, date, profile number, file extension, and possibly comments.
#' Returns `NA` if any of these are missing.
#'
#' @author Tom Shatwell
#'
#'
#' @examples
#' \dontrun{
#' get_profile_info("../mydata/lakes/YH1_20210101_1.CTM644")
#' }
#'
#' @export

# get info from profile filenames
get_profile_info <- function (filename) {
  # remove a path
  splitted    <- strsplit(x=filename, split='/')[[1]]
  filename          <- splitted [length(splitted)]
  ext <- loc <- date <- pr_no <- comment <- NA
  splitted    <- strsplit(x=filename, split='\\.')[[1]]
  l           <-length (splitted)
  if (l > 1 && sum(splitted[1:(l-1)] != ''))  {
    ext <-splitted [l]
    # the extention must be the suffix of a non-empty name
    filename <- paste(splitted[1:l-1],sep=".")
    splitted <- strsplit(x=filename, split="_")[[1]]
    l <- length(splitted)
    loc <- splitted[1]
    if(l>1) date <- splitted[2]
    if(l>2) pr_no <- splitted[3]
    if(l>3) comment <- paste(splitted[4:l],sep="_")
  }
  return(c(loc=loc,date=date,pr_no=pr_no,ext=ext,comment=comment))
}
