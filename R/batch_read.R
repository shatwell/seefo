#' @title Reads raw sonde profile data from several files and compiles the results
#'
#' @description Reads raw sonde profile data from several files and compiles them into a data.frame.
#' By default, it searches a given directory and all subdirectories for matching files for a given sonde type,
#' cleans up the data, and returns it sorted and labelled.
#'
#' @param sonde Type of sonde. One of `"ctm"`,`"bbe"`,`"ixx"`.
#' @param path The path to a directory in which to search for matching files. Used only if `fnames`== `"autosearch"`
#' @param fnames The names of the files to read. Defaults to `"autosearch"`,
#' in which case all filenames with an extension matching the sonde type will be read.
#' If `fnames` != `"autosearch"`, then `path` and `recursive` are not used.
#' @param keep names of the columns to retain in the output. Defaults to `"all"`
#' @param recursive Should subdirectories in `path` also be searched? (logical, only if `fnames`== `"autosearch"`)
#' @param clean Should the profiles be cleaned, by removing measurements in air, near the sediment, and on the upward pull?? (logical)
#' @param TZ time zone given to the output datetimes, defaults to `"UTC"` (recommended)
#'
#' @details
#' This function scans a location for raw sonde files, and matches them by their file extension (.CTM1234 for instance).
#' I then reads all the discovered raw datafiles in the given directory and returns a data.frame containing
#' the all the data sorted chronologically, also with the sampling location (which is determined from the
#' filename). If data are also located in subdirectories, set `recursive = TRUE` to search these as well.
#'
#' @return A `data.frame` containing the profile data
#'
#' @author Tom Shatwell
#'
#' @seealso \code{\link{get_profile_info}},  \code{\link{read_bbe}},  \code{\link{read_ctm}}, \code{\link{clean_profile}}
#'
#' @examples
#' \dontrun{
#' # the location of the folder containing example data
#' loc <- file.path(system.file("extdata", package="seefo"), "profiles")
#'
#' # show the contents
#' list.files(loc, recursive=TRUE)
#'
#' ctm <- batch_read(sonde="ctm", path=loc)
#' head(ctm)
#'
#' bbe <- batch_read(sonde="bbe", path=loc)
#' head(bbe)
#'
#' i72 <- batch_read(sonde="ixx", path=loc)
#' head(i72)
#'
#' rbr <- batch_read(sonde="rbr", path=paste0(loc,"/RBR"))
#' head(rbr)
#'
#' }
#'
#' @export


batch_read <- function(sonde, path, fnames="autosearch", keep="all",
                       recursive=TRUE, clean=TRUE, TZ="UTC") {
  # wrapper function to deal with errors
  wrapper <- function(p, fn, ...) {
    out <- NULL
    tryCatch({
      out <- fn(p, ...)
    }, warning=function(x) {
      out <<- as.character(x)
    }, error=function(x) {
      out <<- as.character(x)
    })
    return(out)
  }

  isValid <- function(x) {is.data.frame(x)}

  filenames <- list.files(path = path,
                          full.names = TRUE, recursive = recursive)
  if(length(filenames)==0) {
    warning("No files to read!")
  }
  if(sonde=="ctm") {
    pat <- ".CTM"
    func <- read_ctm
  } else if(sonde=="bbe") {
    pat <-  ".FP"
    func <- read_bbe
  } else if(sonde == "ixx") {
    pat <- ".i[*0-9]"
    func <- read_ixx
  } else if(sonde == "rbr") {
    pat <- ".RBR"
    func <- read_rbr
  } else {
    stop("Sonde name must be one of 'ctm' 'bbe', 'ixx', 'rbr'.")
  }

  if(fnames[1]=="autosearch") fnames <- filenames[grep(pattern = pat, x = filenames)]

  out <- lapply(X = fnames, FUN = wrapper, fn = func,
                keep = keep, clean=clean, TZ=TZ) # read all the files
  retain <- which(sapply(out, isValid)) # which files were successfully read?
  if(keep[1] == "all") FILL <- TRUE else FILL <- FALSE
  out <- data.table::rbindlist(out[retain], fill = FILL) # stack individual tables
  data.table::setkey(out,"pr_dt") # sort chronologically by profile datetime

  cat("Read", length(retain), "of", length(fnames), "files\n")
  if(length(retain) < length(fnames)) { # print file names that were not successful as warning
    warning(c("Could not read the following files:\n", paste(fnames[!fnames %in% retain],"\n")))
  }

  return(as.data.frame(out))
}

