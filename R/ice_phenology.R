#' @title Calculates ice phenology
#'
#' @description Calculates seasonal ice phenology including annual mean,
#' max and total duration of ice cover.
#'
#' @param H_ice A vector indicating the presence of ice cover.
#' Can be numerical or logical where `TRUE` or any number > `0` (e.g. ice thickness) indicates ice, and
#' `FALSE` or `0` indicates no ice cover.
#' @param dates The dates corresponding to the ice observations in POSIX style and same length as `H_ice`.
#' @param NH Is the lake in the northern hemisphere? (logical)
#'
#' @details
#' To be added later.
#'
#' @return A `data.frame` containing ice stratification phenology
#'
#' @author Tom Shatwell
#'
#'
#' @examples
#' \dontrun{
#' add later
#' }
#'
#' @export
# ice statistics -----------------------------------------------


# this function returns ice phenology statstics: annual mean, max and total
# this is a subset of the function above analyse_strat, but just for ice.
# NOTE: summer strat periods are allocated to the year in which the period
# starts. Winter stratification and ice periods are allocated to the year in
# which they end.
# H_ice: ice thickness time series, set to NULL if analysis not required
# dates: POSIX style date vector corresponding to rows of Ts and Tb.
# not used: dT: temperature difference between top and bottom indicating stratification.
# drho: density difference between top and bottom indicating stratificaiton [kg m^-3]
# NH: northern hemisphere? TRUE or FALSE
ice_phenology <- function(H_ice, dates, NH=TRUE) {
  the_years <- as.POSIXlt(dates)$year+1900
  yrs <- unique(the_years)
  doys <- as.POSIXlt(dates)$yday # day of the year [0..364]
  alt_doys <- doys # alternative counting from [-182 .. 182] (July 2/3) for ice in northern hemisphere or strat in southern hemisphere
  alt_doys[doys>182] <- doys[doys>182] - (365 + leap(the_years[doys>182])) # Jan 1 is day 0, correct for leap years
  alt_years <- the_years
  alt_years[alt_doys<0] <- the_years[alt_doys<0] +1 # alternative counting of years (shifted forward by half a year)

  if(NH) { # NH ice and SH stratification use alternative doy and year counts to adjust for ice and stratification events that span more than one calendar year
    ice_yrs <- alt_years
    ice_doys <- alt_doys
  } else {
    ice_yrs <- the_years
    ice_doys <- doys
  }

  ice <- H_ice > 0
  i_i_st <- diff(c(ice[1],ice))==1 # indices of ice cover onset
  i_i_en <- diff(c(ice[1],ice))==-1 #  # indices of ice cover end
  if(ice[1]) i_i_st <- c(NA, i_i_st) # if initially frozen, set first start date to NA
  if(ice[length(ice)]) i_i_en <- c(i_i_en, NA) # if frozen at end, set last thaw date to NA
  ice_st  <- dates[i_i_st] # ice start dates
  ice_en  <- dates[i_i_en] # ice end dates
  # if(sum(ice)==0)  # if there is no ice at all, set start and end to time=0

  # maximum ice thickness
  IceMax <- NULL
  for(ii in unique(the_years)) {
    Hice_maxi <- which.max(H_ice[ice_yrs == ii])
    IceMaxOut <- data.frame(year=ii,
                            HiceMax     = H_ice[ice_yrs == ii][Hice_maxi],
                            HiceMaxDay  = ice_doys[ice_yrs==ii][Hice_maxi],
                            HiceMaxDate = dates[ice_yrs == ii][Hice_maxi])
    if(sum(H_ice[ice_yrs == ii])==0) IceMaxOut[1,c("HiceMaxDay","HiceMaxDate")] <- NA
    IceMax <- rbind(IceMax, IceMaxOut)
  }

  ice_start_doys <- ice_doys[i_i_st] # day of year of start of ice cover events
  ice_end_doys <- ice_doys[i_i_en]   # day of year of end of ice cover events
  ice_event_yrs <- ice_yrs[i_i_st]   # the years assigned to each ice event

  # if there is no ice, set values to NA ...
  if(sum(ice)==0) {
    ice_start_doys <- ice_end_doys <- ice_event_yrs <- ice_st <- ice_en <- NA
  }
  ice_dur <- as.double(difftime(ice_en, ice_st, units="days")) # duration of ice periods


  # summary of ice cover events
  ice.summary <- data.frame(year = ice_event_yrs,
                            start = ice_st,
                            end = ice_en,
                            dur = ice_dur,
                            startday = ice_start_doys,
                            endday = ice_end_doys)

  ice_out <- NULL
  for(mm in unique(ice.summary$year[!is.na(ice.summary$year)])) {
    ice2 <- subset(ice.summary, year==mm)
    ice2_on <- ice2[which.max(ice2$dur),"startday"]
    ice2_off <- ice2[which.max(ice2$dur),"endday"]
    if(anyNA(ice2$dur)) ice2_on <- ice2_off <- NA

    ice_out <- rbind(ice_out,
                     data.frame(year=mm,
                                MeanIceDur=mean(ice2$dur),
                                MaxIceDur=max(ice2$dur),
                                TotIceDur=sum(ice2$dur),
                                ice_on=ice2_on,
                                ice_off=ice2_off,
                                firstfreeze=min(ice2$startday),
                                lastthaw=max(ice2$endday)))
  }

  ice_out <- ice_out[ice_out$year %in% yrs,] # trim years outside the simulation range (eg ice that forms at the end of the last year, which should be assigned to the following year outside the simulation period)

  ice_out1 <- data.frame(year=yrs, MeanIceDur=NA, MaxIceDur=NA,
                         TotIceDur=NA, IceOn=NA, IceOff=NA, FirstFreeze=NA,
                         LastThaw=NA, HiceMax=NA, HiceMaxDay=NA)
  ice_out1[match(ice_out$year, yrs),
           c("MeanIceDur","MaxIceDur","TotIceDur",
             "IceOn","IceOff","FirstFreeze","LastThaw")] <- ice_out[,-1]
  ice_out1[,c("HiceMax","HiceMaxDay")] <- IceMax[,c("HiceMax","HiceMaxDay")]

  i8 <- ice_out1$IceOff < ice_out1$IceOn
  i8[is.na(i8)]<-FALSE # this gets rid of any NAs
  if(sum(i8, na.rm=TRUE)>0) {
    ice_out1[i8, "IceOff"] <- ice_out1[i8,"IceOn"] + ice_out1[i8,"MaxIceDur"]
  }
  i9 <- ice_out1$LastThaw < ice_out1$IceOn & ice_out1$TotIceDur < 366
  i9[is.na(i9)]<-FALSE # this gets rid of any NAs
  if(sum(i9, na.rm=TRUE)>0) {
    ice_out1[i9, "LastThaw"] <- ice_out1[i9,"LastThaw"] + 365 +
      leap(ice_out1[i9,"year"])
  }

  return(ice_out1)
}
