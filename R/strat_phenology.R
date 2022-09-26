#' @title Calculates stratification phenology
#'
#' @description Calculates seasonal stratification phenology
#' including annual mean,
#' max and total duration of summer stratification and ice cover.
#'
#' @param Ts Surface temperature, as daily values (numeric vector).
#' @param Tb Bottom temperature, as daily values (numeric vector with same length as `Ts`.
#' @param H_ice Ice thickness as numeric vector of the same length as `Ts`.
#' `NULL` if ice statistics should not be calculated.
#' @param dates The dates corresponding to the temperature measurements
#' in POSIX style and same length as `Ts`.
#' @param thresh The density or temperature threshold used to define the
#' presence of stratification, ie whether the difference between Ts and Tb
#' exceeds this threshold. If `bydensity` is `TRUE`, then `thresh` must
#' be in (kg m^-3), and `Ts` and `Tb` must be in Celsius.
#' @param NH Is the lake in the northern hemisphere? (logical)
#' @param bydensity Should `thresh` be defined as a density
#' threshold (logical)?
#' If `FALSE`, a temperature threshold is assumed.
#' If `TRUE`, then `Ts` and `Tb` must be in degrees Celsius
#'
#' @details
#' These stratification phenology statistics are calculated:
#' `MaxStratDur` is the duration of the longest uninterupted stratification period,
#' `TotalStratDur` is the sum of all stratified days per year.
#' `MeanStratDur` is the average duration when more than one stratification
#' period per year exists.
#' `StratStart` and `StratEnd` are the day of the year (Jan 1 = day 0) when
#' the longest uninterupted stratification period begins and ends.
#' `StratFirst` and `StratLast` are the day of year for the earliest and latest
#' stratified days in the annual season.
#' Ice phenology is also calculated if `H_ice` is provided. The definitions of
#' the ice cover periods are analogous to those for stratification.
#'
#' The function determines that stratification exists when the difference
#' between surface (`Ts`) and bottom temperature (`Tb`) exceeds a
#' certain threshold (`thresh`). This can be defined either in terms of a
#' temperature difference or a density difference (`bydensity = TRUE`).
#' Commonly used thresholds are 1 degree Celsius or 0.1 kg/m2.
#' Density is calculated from temperature, using the formula
#' of Millero & Poisson (1981) for freshwater. Only summer (positive) stratification
#' is considered and winter (inverse) stratification is ignored because it is usually
#' transient and from experience is not easy to calculate reliably
#' (ice is perhaps a better proxy). You should choose the
#' depths of `Ts` and `Tb` to be representative for how you want to calculate
#' the stratification duration. It is important that daily values of
#' `Ts` and `Tb` are supplied, so you may need to interpolate if the
#' data resolution is lower than daily. I am planning on improving this
#' sometime in the future. Ice cover is inferred when `H_ice` is not zero,
#' so you can supply an ice thickness or simply `TRUE` or `FALSE` or `0` or `1`.
#'
#' The function defines the "summer" stratification season as 1 Jan to 31 Dec for
#' northern hemisphere lakes (`NH = TRUE`) and from 2 July until 1 July in the
#' following year for southern hemisphere lakes (`NH = FALSE`). Every stratification
#' period is assigned to the year of the season in which it began.
#' For instance, if a northern hemisphere lake stratifies in April 2015 and remains stratified
#' until January 2016, then it will be assigned to 2015, and the `StratEnd` will be the
#' day of year counting from 1 Jan 2015, so that statistics are still valid and
#' `MaxStratDur = StratEnd - StratStart`, even if stratification periods extend into
#' seasons of subsequent years. The same rules apply for ice cover, but the "winter" ice
#' season is defined as 2 July - 1 July in the northern hemisphere and 1 Jan - 31 Dec
#' in the southern hemisphere.
#'
#' @return A `data.frame` containing the stratification phenology
#'
#' @author Tom Shatwell
#'
#'
#' @examples
#' \dontrun{
#' # Load some temperature data
#'
#' data(Ts_Tb_ice)
#'
#' strat_phenology(Ts = Ts_Tb_ice$Ts, Tb = Ts_Tb_ice$Tb, dates = Ts_Tb_ice$date)
#'
#'
#' strat_phenology(Ts = Ts_Tb_ice$Ts, Tb = Ts_Tb_ice$Tb,
#'                        H_ice =  Ts_Tb_ice$H_ice, dates = Ts_Tb_ice$date)
#'
#' }
#'
#' @export

# stratification statistics -----------------------------------------------

strat_phenology <- function(Ts, Tb, H_ice=NULL, dates, thresh=1,
                            NH=TRUE, bydensity = FALSE) {
  the_years <- as.POSIXlt(dates)$year+1900
  yrs <- unique(the_years)
  doys <- as.POSIXlt(dates)$yday # day of the year [0..364]
  # alternative counting from [-182 .. 182] for ice in northern hemisphere
  # or stratification in southern hemisphere
  alt_doys <- doys
  # Jan 1 is day 0, correct for leap years
  alt_doys[doys>182] <- doys[doys>182] - (365 + leap(the_years[doys>182]))
  alt_years <- the_years
  # alternative counting of years (shifted forward by half a year)
  alt_years[alt_doys<0] <- the_years[alt_doys<0] +1

  # NH ice and SH stratification use alternative doy and year counts
  # to adjust for ice and stratification events that span more than
  # one calendar year
  if(NH) {
    ice_yrs <- alt_years
    ice_doys <- alt_doys
    strat_yrs <- the_years
    strat_doys <- doys
  } else {
    ice_yrs <- the_years
    ice_doys <- doys
    strat_yrs <- alt_years
    strat_doys <- alt_doys
  }

  if(bydensity) {
    # logical whether stratified at each time step
    s_strat <- (rho_water(t=Tb) - rho_water(t=Ts)) >= thresh & Ts > Tb
  } else {
    # logical whether stratified at each time step
    s_strat <- Ts - Tb  > thresh
  }

  # indices of stratification onset
  i_s_st <- diff(c(s_strat[1],s_strat))==1
  # indices of stratification end
  i_s_en <- diff(c(s_strat[1],s_strat))==-1
  # if stratified at beginning of simulation, make first date NA
  if(s_strat[1]) i_s_st <- c(NA, i_s_st)
  # if stratified at end of sim, set last strat date to NA
  if(s_strat[length(s_strat)]) i_s_en <- c(i_s_en, NA)
  s_start <- dates[i_s_st] # summer strat start dates
  s_end   <- dates[i_s_en] # summer strat end dates
  # if never stratifies, set to time=0
  # if(sum(s_strat)==0) s_start <- s_end <- dates[1]
  # duration of summer stratification periods
  s_dur   <- as.double(difftime(s_end, s_start, units="days"))

  a1 <- data.frame(year=strat_yrs[i_s_st],
                   start=s_start, end=s_end, dur=s_dur,
                   startday = strat_doys[i_s_st],
                   endday = strat_doys[i_s_en])
  # a1 <- subset(a1, a1$year %in% yrs)
  a1 <- a1[a1$year %in% yrs,]

  s.max <- s.mean <- s.tot <- s.on <- s.off <-
    s.first <- s.last <- yr <- NULL
  for(mm in unique(a1$year[!is.na(a1$year)])) {
    # remove NAs which are generated when the lake is
    # stratified at the satrt or end of the time series
    # a2 <- subset(a1, a1$year==mm)
    a2 <- a1[a1$year==mm,]
    ind <- which.max(a2$dur)
    # fixes issue if stratified at end of data period
    if(nrow(a2)==1) if(is.na(a2$dur)) ind <- NA
    yr <- c(yr,mm)
    s.max <- c(s.max,max(a2$dur))
    s.mean <- c(s.mean,mean(a2$dur))
    s.tot <- c(s.tot,sum(a2$dur))
    s.on <- c(s.on, as.POSIXlt(a2$start[ind])$yday)
    s.off <- c(s.off, as.POSIXlt(a2$end[ind])$yday)
    s.first <- c(s.first, min(a2$startday))
    s.last <- c(s.last, max(a2$endday))
  }

  # maximum surface temperature
  # loop thru years to find Tmax and its day of year
  TsMax <- NULL
  for(ii in unique(the_years)) {
    Ts_maxi <- which.max(Ts[strat_yrs == ii])
    TsMaxOut <- data.frame(year=ii,
                           TsMax    = Ts[strat_yrs == ii][Ts_maxi],
                           TsMaxDay = strat_doys[strat_yrs==ii][Ts_maxi],
                           TsMaxDate= dates[strat_yrs == ii][Ts_maxi]
    )

    TsMax <- rbind(TsMax, TsMaxOut)
  }

  # create empty data frame to fill with data
  # (not all years may have strat or ice)
  out <- data.frame(year=yrs, TsMax=NA, TsMaxDay=NA,
                    MaxStratDur=NA, MeanStratDur=NA, TotStratDur=NA,
                    StratStart=NA, StratEnd=NA,
                    StratFirst=NA, StratLast=NA)

  out[match(TsMax$year, yrs), c("TsMax","TsMaxDay")] <-
    TsMax[,c("TsMax","TsMaxDay")]

  out[match(yr, yrs), -1:-3] <-
    data.frame(s.max,s.mean,s.tot,s.on,s.off,s.first,s.last)


  # ice cover
  if(!is.null(H_ice)) { # only do this if ice data provided
    ice <- H_ice > 0
    # indices of ice cover onset
    i_i_st <- diff(c(ice[1],ice))==1
    # indices of ice cover end
    i_i_en <- diff(c(ice[1],ice))==-1
    # if initially frozen, set first start date to NA
    if(ice[1]) i_i_st <- c(NA, i_i_st)
    # if frozen at end, set last thaw date to NA
    if(ice[length(ice)]) i_i_en <- c(i_i_en, NA)
    # ice start dates
    ice_st  <- dates[i_i_st]
    # ice end dates
    ice_en  <- dates[i_i_en]
    # if there is no ice at all, set start and end to time=0
    # if(sum(ice)==0)

    # maximum ice thickness
    IceMax <- NULL
    for(ii in unique(the_years)) {
      Hice_maxi <- which.max(H_ice[ice_yrs == ii])
      IceMaxOut <- data.frame(year=ii,
                              HiceMax    =H_ice[ice_yrs == ii][Hice_maxi],
                              HiceMaxDay =ice_doys[ice_yrs==ii][Hice_maxi],
                              HiceMaxDate=dates[ice_yrs == ii][Hice_maxi])
      if(sum(H_ice[ice_yrs == ii])==0) {
        IceMaxOut[1,c("HiceMaxDay","HiceMaxDate")] <- NA
      }
      IceMax <- rbind(IceMax, IceMaxOut)
    }

    # day of year of start of ice cover events
    ice_start_doys <- ice_doys[i_i_st]
    # day of year of end of ice cover events
    ice_end_doys <- ice_doys[i_i_en]
    # the years assigned to each ice event
    ice_event_yrs <- ice_yrs[i_i_st]

    # if there is no ice, set values to NA ...
    if(sum(ice)==0) {
      ice_start_doys <- ice_end_doys <- ice_event_yrs <-
        ice_st <- ice_en <- NA
    }

    # duration of ice periods
    ice_dur <- as.double(difftime(ice_en, ice_st, units="days"))


    # summary of ice cover events
    ice.summary <- data.frame(year = ice_event_yrs,
                              start = ice_st,
                              end = ice_en,
                              dur = ice_dur,
                              startday = ice_start_doys,
                              endday = ice_end_doys)

    ice_out <- NULL
    for(mm in unique(ice.summary$year[!is.na(ice.summary$year)])) {
      # ice2 <- subset(ice.summary, ice.summary$year==mm)
      ice2 <- ice.summary[ice.summary$year==mm,]
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

    # trim years outside the simulation range
    # (eg ice that forms at the end of the last year,
    # which should be assigned to the following year
    # outside the simulation period)
    ice_out <- ice_out[ice_out$year %in% yrs,]

    ice_out1 <- data.frame(year=yrs,
                           MeanIceDur=NA,
                           MaxIceDur=NA,
                           TotIceDur=NA,
                           IceOn=NA,
                           IceOff=NA,
                           FirstFreeze=NA,
                           LastThaw=NA,
                           HiceMax=NA,
                           HiceMaxDay=NA)
    ice_out1[match(ice_out$year, yrs),
             c("MeanIceDur","MaxIceDur","TotIceDur",
               "IceOn","IceOff","FirstFreeze",
               "LastThaw")] <- ice_out[,-1]
    ice_out1[,c("HiceMax","HiceMaxDay")] <- IceMax[,c("HiceMax","HiceMaxDay")]

    out <- data.frame(out, ice_out1[,-1])

  }

  # adjust some exceptions where stratification or
  # ice extend longer than the cutoff period
  i6 <- out$StratEnd < out$StratStart
  i6[is.na(i6)]<-FALSE # this gets rid of any NAs
  if(sum(i6, na.rm=TRUE)>0) out[i6, "StratEnd"] <- out[i6,"StratStart"] +
    out[i6,"MaxStratDur"]
  i7 <- out$StratLast < out$StratStart & out$TotStratDur < 365
  i7[is.na(i7)]<-FALSE # this gets rid of any NAs
  if(sum(i7, na.rm=TRUE)>0) out[i7, "StratLast"] <- out[i7,"StratLast"] + 364
  i8 <- out$IceOff < out$IceOn
  i8[is.na(i8)]<-FALSE # this gets rid of any NAs
  if(sum(i8, na.rm=TRUE)>0) out[i8, "IceOff"] <- out[i8,"IceOn"] +
    out[i8,"MaxIceDur"]
  i9 <- out$LastThaw < out$IceOn & out$TotIceDur < 366
  i9[is.na(i9)]<-FALSE # this gets rid of any NAs
  if(sum(i9, na.rm=TRUE)>0) {
    out[i9, "LastThaw"] <- out[i9,"LastThaw"] + 365 +
      leap(out[i9,"year"])
  }
  return(out)
}

