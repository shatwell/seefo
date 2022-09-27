
# seefo: an R package for data analyses in lake research

<!-- badges: start -->
<!-- badges: end -->

The seefo package is a collection of useful functions and tools for common tasks in lake research in the SEEFO department at UFZ. It is a growing, collaborative project and additional contributions are welcome and encouraged. 
There are tools for calculating stratification and ice phenology, downloading and reading ERA5 reanalysis data, working with time series, like interpolating profiles to fixed depths, or interpolating over time to daily values, and calculating volume-weighted means with morphometric information. There are also functions for calculating nutrient retention in lakes or reservoirs, and for reading and working with raw sensor data for common sensors in the department.

## Installation

You can install the latest version of seefo from [gitlab](https://git.ufz.de/shatwell/seefo) with:

``` r
# install.packages("devtools")
devtools::install_gitlab(repo="shatwell/seefo", host="https://git.ufz.de/")
```
Alternatively, if you run into problems here, you can download the compressed file
`seefo_x.x.x.xxxx.tar.gz` (`x`s refer to the current version number) and install with:

``` r
install.packages("path/to/package/file/on/your/computer/seefo_v.e.r.sion.tar.gz"", type='source')
```
Don't forget to use the correct package name.

## Overview

The package contains the following functions:

`retention_eff()` calculates nutrient retention
`batch_read()` read multiple raw sonde files and compile into a data.frame.
`read_ctm()` read a ctm sonde raw datafile
`read_bbe()` read a bbe sonde raw datafile
`read_ixx()` read an ixx (e.g. i172) sonde raw datafile
`clean_profile()` cleans up a profile by removing unwanted values
`get_profile_info()` parses the raw sonde dataname to extract information (used in other functions)
`profile2grid()` interpolates multiple profiles in long format to fixed depths in wide format
`strat_phenology()`, `ice_phenology()` calculate stratification and ice phenology
`hmix()`, `h2mix()`, `h3mix()`, `delta_rho()` calculate the mixed layer depth in lakes
`era5cli()` generates a download command for using with the era5cli (command line interface)
`gather_era5_var()` reads and collates era5 reanalysis data from multiple netcdf download files for a given location
`dewpoint2humidity()` converts dewpoint temperature to humidity
`rho_water()` calculates water density
`leap()` is a helper function that determines if a year is a leap year


## Examples

### Working with raw sensor data from the limnophysics folder

The seefo package helps you read raw sensor data on the limnophysics directory. There are
currently functions to read data from our CTM, BBE and ixx (e.g. i172) sondes. Additional
functions for other sondes (e.g. EXOs or RBRs) will follow.

``` r
library(seefo)

pr_loc <- system.file("extdata", package="seefo")
pr_name <- file.path(pr_loc,  "profiles","20170801","YT1_20170801_1.CTM644")
ctm <- read_ctm(pr_name)
head(ctm)
```
You can clean up raw profiles with `clean_profile()`, which removes air values and upward haul values.
``` r
ctm_clean <- clean_profile(ctm)
plot(-press~temp,ctm, "l")
```

Probably the most useful function here is `batch_read()`, which will scan a directory and read all the 
raw datafiles contained there, including in subdirectories (if `recursive=TRUE`). It automatically cleans
the profiles, and collates the data into one data.frame.
``` r
loc <- file.path(system.file("extdata", package="seefo"), "profiles")
bbe <- batch_read(sonde="bbe", path=loc)
head(bbe)
```

From here you can use other useful functions like `profile2grid()` to convert the data from 
long format to a wide format, interpolated to given depths:
``` r
chla <- profile2grid(data=bbe, datetime="pr_dt", depth="depth", value="total_conc")

```
It returns a matrix for easy plotting, and important information like the interpolated depths
and the profile times are stored in the attributes.
This can now be easily plotted:
``` r
depths <- attributes(chla)$outdepths
times <- attributes(chla)$dt

library(plot3D)
image2D(z = t(chla), y=depths, x=times, ylim=c(max(depths),0))

```

###  Downloading and reading ERA5 reanalysis data

There are a number of functions available for downloading and reading ERA5 reanalysis data.
This is still in development.

You can also read the downloaded data easily with `gather_era5_var()`.

``` r
loc <- file.path(system.file("extdata", package="seefo"), "nc")
loc <- paste0(loc,"/") # add file separator (I should fix this)
out <- gather_era5_var(namestrings = c("era5", "cloud_base_height"),
                       varid = "cbh", lon = 5.49, lat = 51.06,
                       datapath = loc)
head(out)
```

### Calculating lake nutrient retention

The function `retention_eff()` calculates nutrient retention in lakes. 
There's more to  come, as it's under development.

### Calculating stratification information

You can use `strat_phenology()` to calculate stratification timing metrics, and `ice_phenology()` for ice metrics. 
Calculating these metrics is very tricky, and `strat_phenology()` does many checks to get it right and has been tested on over 60 lakes and works fairly reliably. You should use daily data when using it:

``` r
data(Ts_Tb_ice)
out <- strat_phenology(Ts = Ts_Tb_ice$Ts, Tb = Ts_Tb_ice$Tb, dates = Ts_Tb_ice$date)
head(out)

```
`hmix()` and `h2mix()` calculate the bottom of the epilimnion and top of the hypolimnion based on the 
density gradient. The are designed to run fast on huge modelling datasets, so they don't do any checks.
`h3mix()` is much more robust and is useful for measured and/or messy and irregular data. `delta_rho()` estimates mixed layer depth based on a density threshold and is quite robust when all else fails:

``` r
pr_loc <- system.file("extdata", package="seefo")
pr_name <- file.path(pr_loc,  "profiles","20170704","YT1_20170704_1.CTM644")
pr <- read_ctm(pr_name)
head(pr)

h3mix(T=pr$temp, z=pr$press, plot=TRUE)
````
