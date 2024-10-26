% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dewpoint_2_humidity.R
\name{dewpoint_2_humidity}
\alias{dewpoint_2_humidity}
\title{Calculate humidity from dewpoint temperature}
\usage{
dewpoint_2_humidity(Td, E0 = 610.94, L_Rv = 5423, T0 = 273.15)
}
\arguments{
\item{Td}{Dewpoint temperature in Kelvin (numeric).}

\item{E0}{Constant (numeric).}

\item{L_Rv}{Constant (numeric).}

\item{T0}{Freezing point in Kelvin (numeric).}
}
\value{
The humidity in mbar.
}
\description{
Calculate humidity from dewpoint temperature after Clausius-Claperon equations.
}
\details{
Useful for converting between units for hydrodynamic modelling.
To calculate saturation humidity, substitute Td = air temperature instead of dewpoint temperature.
}
\examples{
\dontrun{
 dewpoint_2_humidity(293)
}

}
\references{
Lawrence, Mark G., 2005: The relationship between relative humidity and the dewpoint temperature in moist air: A simple conversion and applications. Bull. Amer. Meteor. Soc., 86, 225-233. doi: http;//dx.doi.org/10.1175/BAMS-86-2-225)
}
\author{
Tom Shatwell
}
