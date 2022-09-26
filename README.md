
# seefo

<!-- badges: start -->
<!-- badges: end -->

The seefo package is a collection of useful functions and tools for common tasks in lake research in the SEEFO department at UFZ. It is a growing, collaborative project and additional contributions are welcome and encouraged. 
These tools include: calculating ice and stratification phenology, reading raw data files from sondes, calculating water density and leap years, calculating volume-weighted mean using hyspographic information, interpolating to fixed depth grids or daily values.

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
`strat_phenology()` was the hardest function to program as ice and stratification phenology is tricky, 
but it has been tested on over 60 lakes and works reliably.:

``` r
data(Ts_Tb_ice)
out <- strat_phenology(Ts = Ts_Tb_ice$Ts, Tb = Ts_Tb_ice$Tb, dates = Ts_Tb_ice$date)
head(out)

```
`hmix()` and `h2mix()` calculate the bottom of the epilimnion and top of the hypolimnion based on the 
density gradient. The are designed to run fast on huge modelling datasets, so they don't do any checks.
`h3mix()` is much more robust and is useful for measured and/or messay and irregular data:

``` r
pr_loc <- system.file("extdata", package="seefo")
pr_name <- file.path(pr_loc,  "profiles","20170704","YT1_20170704_1.CTM644")
pr1 <- read_ctm(pr_name)
head(pr1)

h3mix(T=pr1$temp, z=pr1$press, plot=TRUE)
````
