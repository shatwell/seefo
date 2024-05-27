#' @title Import and collate ERA5 variable from different files at a location
#'
#' @description Gathers, imports and collates an ERA5 variable from different netcdf files.
#' Typically ERA5 files are downloaded containing 1 variable at several coordinates for a certain period.
#' This function will search files that match a variable name string close to a given coordinate, read the data and assemble the variable from these files in chronological order. This tidies up ERA5 data downloads eg from era5cli program.
#'
#'
#' @param namestrings list of vector containing text fragments to uniquely match the  filenames. Typically it will contain the file extension (".nc$"), variable name (e.g. "2m_temperature", and perhaps a prefix (e.g. "era5"). This is only used if files is "autosearch".
#' @param varid The name of the variable in the nc file. See ERA5 shortname in the table in Details
#' @param lon Longitude (numerical) of the location to collate data.
#' @param lat Latitude (numerical) of the location to collate data.
#' @param datapath The path to the folder in which to search for the files (character)
#' @param files The filenames to be read (character vector). The default is `"autosearch"`, in which case the files are matched with the patterns in namestrings.
#' @param recursive Should subfolders also be searched? (logical)
#' @param radius The radius (numerical vector length 2) in degrees around `[lat,lon]` that will be searched. If discovered coordinates are outside `radius[1]`, a warning is issued, if they are outside `radius[2]`, an error is issued.
#' @param pad_missings If there are fewer values than hours in the time series, should they be padded with NAs? (logical)
#'
#' @details The location and variable name are kept in the attributes of the output. Pad_missings can be useful if a variable can take on missing values, e.g. cloud base height when there is no cloud. Padding missing values can mask problems like discontinuous time series, so choose FALSE if for safety if the value cannot be NA.
#'
#' | **ERA5 long name**                                       | **ERA5 shortname** | **Units** |
#' |----------------------------------------------------------|----------------|-------|
#' |surface_solar_radiation_downwards                         | ssrd           | J m**-2    |
#' |surface_solar_radiation_downward_clear_sky                | ssrdc          | J m**-2    |
#' |mean_surface_downward_short_wave_radiation_flux           | msdwswrf       | W m**-2    |
#' |mean_surface_downward_short_wave_radiation_flux_clear_sky | msdwswrfcs     | W m**-2    |
#' |mean_surface_downward_long_wave_radiation_flux            | msdwlwrf       | W m**-2 |
#' |surface_thermal_radiation_downwards                       | strd           | J/m^2 |
#' |2m_temperature                                            | t2m            | K     |
#' |surface_pressure (surface_air_pressure)                   | sp             | Pa    |
#' |2m_dewpoint_temperature                                   | d2m            | K     |
#' |10m_u_component_of_wind                                   | u10            | m/s   |
#' |10m_v_component_of_wind                                   | v10            | m/s   |
#' |mean_total_precipitation_rate                             | mtpr           | kg m**-2 s**-1|
#' |total_precipitation                                       | tp             | m     |
#' |total_cloud_cover                                         | tcc            | 0..1     |
#' |low_cloud_cover                                           | lcc            | 0..1     |
#' |medium_cloud_cover                                        | mcc            | 0..1     |
#' |high_cloud_cover                                          | hcc            | 0..1     |
#' |cloud_base_height                                         | cbh            | m     |
#'
#' @return A data.frame containing the datetime and variable values
#' @author Tom Shatwell
#'
#' @examples
#' \dontrun{
#' ta <- gather_era5_var(namestrings = "2m_temperature", varid = "t2m",
#' datapath="../nc_input/era5_north_aral/",
#' lon=60.25, lat=46.5)
#'
#'
#' # the location of the folder containing example data
#'
#' loc <- file.path(system.file("extdata", package="seefo"), "nc")
#' loc <- paste0(loc,"/") # add file separator (I should fix this)
#'
#' # show the contents
#' list.files(loc, recursive=TRUE)
#'
#' out <- gather_era5_var(namestrings = c("era5", "cloud_base_height"),
#'                        varid = "cbh", lon = 5.49, lat = 51.06,
#'                        datapath = loc)
#' }
#'
#'
#' @export

gather_era5_var <- function(namestrings=c(".nc$"), varid, lon, lat,
                            datapath, files="autosearch", recursive=FALSE,
                            radius=c(0.5,1), pad_missings = TRUE) {
  if(files=="autosearch") {
    files <- list.files(path = datapath, recursive = recursive) # list ncdf files
    for(i in 1:length(namestrings)) {
      files <- files[grep(namestrings[i], files)]
    }
    if(length(files)==0) stop("No matching files for this variable")
  }

  if(radius[2] < radius[1]) {
    stop(paste0("Radius[2] (",radius[2],
                ") must be greater than or equal to radius[1] (",
                radius[1],")."))
  }

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

  # validity checker
  isValid <- function(x) {!is.character(unlist(x))}

  # main function to read a netcdf file close to lat/lon
  gather_one <- function(file) {
    nc <- ncdf4::nc_open(filename = paste0(datapath,file))
    nclon  <- ncdf4::ncvar_get(nc, "longitude")  # longitude [deg E]
    nclat  <- ncdf4::ncvar_get(nc, "latitude")  # latitude [deg N]
    tim  <- ncdf4::ncvar_get(nc, "time") # time in hours since tim0
    if(nc$ndims==4) {
      warning(paste0("ERA5 and ERA5T are both contained in file ",file,
                     ". Reading only ERA5 data."))
      expver <- ncdf4::ncvar_get(nc,"expver")
    }
    tim0 <- as.POSIXct(strsplit(
      ncdf4::ncatt_get( nc, varid="time")$units,
      "hours since ")[[1]][2], tz="UTC")
    dt <- tim0+tim*3600
    dist <- sqrt(outer((lon-nclon)^2, (lat-nclat)^2, FUN="+")) # distances between coords and desired point
    mindist <- min(dist)
    indlon <- which(apply(dist==mindist, 1, sum)==1)
    indlat <- which(apply(dist==mindist, 2, sum)==1)
    nearlon <- nclon[indlon]
    nearlat <- nclat[indlat]
    start <- c(indlon, indlat, 1)
    count <- c(1,1,-1)
    if(nc$ndims==4) {
      start <- c(indlon,indlat,expver[1],1)
      count <- c(1,1,1,-1)
    }
    dat <- ncdf4::ncvar_get(nc, varid, start=start, count=count)
    ncdf4::nc_close(nc)
    if(mindist > radius[2]) {
      stop(paste("Nearest location is more than", radius[2], "degrees away"))
    }
    return(data.table::data.table(dt=dt,val=dat,
                                  lon=nearlon,
                                  lat=nearlat,
                                  mindist=mindist))
  }

  out <- lapply(X = files, FUN = wrapper, fn = gather_one) # read all files
  retain <- which(sapply(out, isValid)) # which files were successfully read?
  if(length(retain)==0) stop(print(out))
  out <- data.table::rbindlist(out[retain]) # stack individual tables
  data.table::setkey(out,"dt") # sort chronologically by profile datetime
  out <- as.data.frame(out)

  cat("Read", length(retain), "of", length(files), "files:\n",
      paste(files[retain],"\n"))
  if(length(retain) < length(files)) { # print file names that were not successful as warning
    warning(c("Did not read/use the following files:\n",
              paste(files[!(1:length(files)) %in% retain],"\n")))
  }

  nearlats <- unique(out$lat) # latitudes of component data sets
  nearlons <- unique(out$lon) # longitudes of component data sets

  if(any(out$mindist > radius[1])) {
    warning(paste("At least one location is more than", radius[1], "degrees away.\n"))
  }
  if((length(nearlats)+length(nearlons)) > 2) {
    warning(paste("Coordinates from",
                  max(c(length(nearlats),length(nearlons))),
                  "different locations. Component coordinates stored in attributes.\n"))
  }

  if(nrow(out) > (difftime(out[nrow(out),"dt"],
                           out[1,"dt"], units="hour")+1)) {
    out <- stats::aggregate(val~dt, data=out, FUN=mean)
    out$dt <- as.POSIXct(out$dt,tz="UTC")
    if(nrow(out) == (difftime(out[nrow(out),"dt"],
                              out[1,"dt"], units="hour")+1)) {
      warning(paste0("More than one matching file. Duplicate values were averaged.\n",
                     paste(files, collapse=" \n")))
    }
    if(nrow(out) < (difftime(out[nrow(out),"dt"],
                             out[1,"dt"], units="hour")+1)) {
      if(pad_missings) {
        dt1 <- data.table::setDT(list(dt=seq(out[1,"dt"], out[nrow(out),"dt"], "hour"))) # hourly time sequence
        data.table::setkey(dt1,"dt")
        out <- dt1[data.table::setDT(out,key="dt")] # merge results with hour sequence, adding NAs for missings
        warning(paste0("More than one matching file. Duplicate values were averaged and missings were padded with NAs.\n",
                       paste(files, collapse=" \n")))
      } else { # end if pad_missings
        stop(paste0("More than one matching file. Any duplicate values were averaged but the time series appears not continuous.\n",
                    paste(files, collapse=" \n")))
      } # end else stop
    } # end else if nrow(out) < (difftime...)
  } # end if(nrow(out) > (difftime...))
  out <- as.data.frame(out[,c("dt","val")])
  attributes(out) <- c(attributes(out),
                       lat=lat, lon=lon,
                       varid=varid,
                       component_lats=nearlats,
                       component_lons=nearlons,
                       files=list(files[retain]))
  return(out)
}

## alternative function definition (may be less buggy than the one above)
# gather_era5_var <- function(namestrings=c(".nc$"), varid, lon, lat,
#                             datapath, files="autosearch", recursive=FALSE,
#                             radius=c(0.5,1), pad_missings = TRUE) {
#   if(files=="autosearch") {
#     files <- list.files(path = datapath, recursive = recursive) # list ncdf files
#     for(i in 1:length(namestrings)) {
#       files <- files[grep(namestrings[i], files)]
#     }
#     if(length(files)==0) stop("No matching files for this variable")
#   }
#
#   if(radius[2] < radius[1]) {
#     stop(paste0("Radius[2] (",radius[2],
#                 ") must be greater than or equal to radius[1] (",
#                 radius[1],")."))
#   }
#
#   wrapper <- function(p, fn, ...) {
#     out <- NULL
#     tryCatch({
#       out <- fn(p, ...)
#     }, warning=function(x) {
#       out <<- as.character(x)
#     }, error=function(x) {
#       out <<- as.character(x)
#     })
#     return(out)
#   }
#
#   # validity checker
#   isValid <- function(x) {!is.character(unlist(x))}
#
#   # main function to read a netcdf file close to lat/lon
#   gather_one <- function(file) {
#     nc <- ncdf4::nc_open(filename = paste0(datapath,file))
#     nclon  <- ncdf4::ncvar_get(nc, "longitude")  # longitude [deg E]
#     nclat  <- ncdf4::ncvar_get(nc, "latitude")  # latitude [deg N]
#     tim  <- ncdf4::ncvar_get(nc, "time") # time in hours since tim0
#     if(nc$ndims==4) {
#       warning(paste0("ERA5 and ERA5T are both contained in file ",file,
#                      ". Reading only ERA5 data."))
#       expver <- ncdf4::ncvar_get(nc,"expver")
#     }
#     tim0 <- as.POSIXct(strsplit(
#       ncdf4::ncatt_get( nc, varid="time")$units,
#       "hours since ")[[1]][2], tz="UTC")
#     dt <- tim0+tim*3600
#     dist <- sqrt(outer((lon-nclon)^2, (lat-nclat)^2, FUN="+")) # distances between coords and desired point
#     mindist <- min(dist)
#     indlon <- which(apply(dist==mindist, 1, sum)==1)
#     indlat <- which(apply(dist==mindist, 2, sum)==1)
#     nearlon <- nclon[indlon]
#     nearlat <- nclat[indlat]
#     start <- c(indlon, indlat, 1)
#     count <- c(1,1,-1)
#     if(nc$ndims==4) {
#       start <- c(indlon,indlat,expver[1],1)
#       count <- c(1,1,1,-1)
#     }
#     dat <- ncdf4::ncvar_get(nc, varid, start=start, count=count)
#     ncdf4::nc_close(nc)
#     if(mindist > radius[2]) {
#       stop(paste("Nearest location is more than", radius[2], "degrees away"))
#     }
#     return(data.table::data.table(dt=dt,val=dat,
#                                   lon=nearlon,
#                                   lat=nearlat,
#                                   mindist=mindist))
#   }
#
#   out <- lapply(X = files, FUN = wrapper, fn = gather_one) # read all files
#   retain <- which(sapply(out, isValid)) # which files were successfully read?
#   if(length(retain)==0) stop(print(out))
#   out <- data.table::rbindlist(out[retain]) # stack individual tables
#   data.table::setkey(out,"dt") # sort chronologically by profile datetime
#   out <- as.data.frame(out)
#
#   cat("Read", length(retain), "of", length(files), "files:\n",
#       paste(files[retain],"\n"))
#   if(length(retain) < length(files)) { # print file names that were not successful as warning
#     warning(c("Did not read/use the following files:\n",
#               paste(files[!(1:length(files)) %in% retain],"\n")))
#   }
#
#   nearlats <- unique(out$lat) # latitudes of component data sets
#   nearlons <- unique(out$lon) # longitudes of component data sets
#
#   if(any(out$mindist > radius[1])) {
#     warning(paste("At least one location is more than", radius[1], "degrees away.\n"))
#   }
#   if((length(nearlats)+length(nearlons)) > 2) {
#     warning(paste("Coordinates from",
#                   max(c(length(nearlats),length(nearlons))),
#                   "different locations. Component coordinates stored in attributes.\n"))
#   }
#
#   if(nrow(out) > (difftime(out[nrow(out),"dt"],
#                            out[1,"dt"], units="hour")+1)) {
#     out <- aggregate(val~dt, data=out, FUN=mean)
#     out$dt <- as.POSIXct(out$dt,tz="UTC")
#     if(nrow(out) == (difftime(out[nrow(out),"dt"],
#                               out[1,"dt"], units="hour")+1)) {
#       warning(paste0("More than one matching file. Duplicate values were averaged.\n",
#                      paste(files, collapse=" \n")))
#     }
#     if(nrow(out) < (difftime(out[nrow(out),"dt"],
#                              out[1,"dt"], units="hour")+1)) {
#       if(pad_missings) {
#         dt1 <- data.table::setDT(list(dt=seq(out[1,"dt"], out[nrow(out),"dt"], "hour"))) # hourly time sequence
#         data.table::setkey(dt1,"dt")
#         out <- dt1[data.table::setDT(out,key="dt")] # merge results with hour sequence, adding NAs for missings
#         warning(paste0("More than one matching file. Duplicate values were averaged and missings were padded with NAs.\n",
#                        paste(files, collapse=" \n")))
#       } else { # end if pad_missings
#         stop(paste0("More than one matching file. Any duplicate values were averaged but the time series appears not continuous.\n",
#                     paste(files, collapse=" \n")))
#       } # end else stop
#     } # end else if nrow(out) < (difftime...)
#   } # end if(nrow(out) > (difftime...))
#   out <- as.data.frame(out[,c("dt","val")])
#   attributes(out) <- c(attributes(out),
#                        lat=lat, lon=lon,
#                        varid=varid,
#                        component_lats=nearlats,
#                        component_lons=nearlons,
#                        files=list(files[retain]))
#   return(out)
# }

