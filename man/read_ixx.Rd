% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_ixx.R
\name{read_ixx}
\alias{read_ixx}
\title{Reads profile data from raw output file from a CTD of type ixxx}
\usage{
read_ixx(filename, keep = "all", clean = TRUE, TZ = "UTC")
}
\arguments{
\item{filename}{Name of the file to read, including the file extension}

\item{keep}{names of the columns to retain in the output. Defaults to \code{"all"}}

\item{clean}{Should the profile be cleaned? (logical)}

\item{TZ}{time zone given to the output datetimes, defaults to \code{"UTC"}(recommended)}
}
\value{
A \code{data.frame} containing the profile data
}
\description{
Reads a data from an output file generated by a ixx CTD.
}
\details{
This function first collects information from the filename using \code{get_profile_info()},
including the sampling location from the beginning of the filename,
and the sonde number from the file extension.
It begins to read data in the file 2 lines after the keyword "Datasets".
It assigns the names in \code{header} to the columns. If \code{header} is set to \code{"auto"},
it will assign hard coded names according to the sonde number, hoping that one sonde always thas the same sensors attached.
It converts the date and time columns to a POSIX datetime with the timezone \code{TZ}.
The function finally drops any columns that are not named in \code{keep}.
If \code{keep} contains names not found in the header, it will add a column of \code{NA}s for each unfound name.
}
\examples{
\dontrun{
pr_loc <- system.file("extdata", package="seefo")
pr_name <- file.path(pr_loc,  "profiles","20170801","YTN_20110321_1.i172")
pr1 <- read_ixx(pr_name)
head(pr1)
}

}
\seealso{
\code{\link{get_profile_info}},  \code{\link{read_ctm}}, \code{\link{read_bbe}}, \code{\link{batch_read}}
}
\author{
Tom Shatwell
}
