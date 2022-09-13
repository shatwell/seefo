#' @title Nutrients load and retention efficiency calculator
#'
#' @description Calculates the nutrient load in tons/year using Standard Method 1 and 2 acoording to Hilde 2003 and Generalized Additive Model (GAM).
#'
#'
#' @param data Dataframe with date, variable, value, in_outlet and tributary
#' @param variable Vector with nutrient's name
#' @param value Vector with the observed concentrations in mg/l and the discharge in (m3/s)
#' @param in_outlet Vector with characters "inflow" or "outflow"
#' @param tributary Vector with characters naming the pre-dams
#' @param startyear Numeric length 1
#' @param endyear Numeric length 1
#' @param methods A vector containing one or more of the available methods ("GAM.load", "method1", "method2")
#'
#'
#' @details Needs lubridate and mgcv packages.
#' @return A dataframe with annual load and retention efficiency for each nutrient
#' @author Karsten Rinke and Taynara Fernandes
#'
#' @examples
#'
#' \dontrun{
#' methods <- c("method2", "GAM")
#' startyear=2000
#' endyear=2017
#' }
#'
#' @export
#

## --> TS these should be arguments to the function
#data provided by the user
methods <- c("GAM.load", "method1") #here the user can choose 3 out of 3 methods
start.year = 2000
end.year = 2017
## <-- TS

#define 3 functions that do the different calculations

# --> TS: this information should be in the header info above
#first the GAM. "hydrology" must be a data.frame consisting of times,year, doy, and discharge over
#the timespan calculations are done. "year", "doy","discharge", "concentration" are vectors with
#the observed concentrations
# <-- TS
load.GAM <- function(hydrology, year, doy, discharge, concentration, GOF=TRUE){
  mydata <- data.frame(year= year, doy=doy, discharge=discharge, concentration=concentration)
  myGAM <-  mgcv::gam(concentration ~ s(year)+s(doy, bs="cc")+s(discharge), data=mydata)
  dev.expl <- summary(myGAM)$dev.expl
  if(GOF){print(dev.expl)}
  predicted <- data.frame(times = hydrology$times,
                          Q     = hydrology$discharge,
                          year  = hydrology$year,
                          pred.c= predict.gam(myGAM,
                                              newdata = data.frame(year=hydrology$year,
                                                                   doy =hydrology$doy,
                                                                   discharge= hydrology$discharge)))
  if(sum(is.na(predicted$pred.c))>0){print("warning, there are NAs in GAM predictions")}

  predicted$daily.load <- predicted$pred.c*predicted$Q*3600*24 *1e-6 #in t/d
  yearly.load <- aggregate(predicted$daily.load, by=list(predicted$year),sum)
  names(yearly.load) <- c("year","load.tons")
  return(yearly.load)
}

#standard method 1
load.method1 <- function(year, discharge, concentration){

  daily.loads.measured <- discharge * concentration *3600*24*1e-6 #in t/d
  average.daily.loads  <- aggregate(daily.loads.measured, by=list(year), mean)
  names(average.daily.loads) <- c("year","mean.daily.load.tons")
  yearly.loads <- data.frame(year = average.daily.loads$year,
                             yearly.load.tons = average.daily.loads$mean.daily.load.tons*365)
  return(yearly.loads)
}

#standard method 2 (Qweighting)
load.method2 <- function(hydrology, year, discharge, concentration){
  loads.from.method1 <- load.method1(year=year, discharge=discharge, concentration=concentration)

  mean.sampled.Q <- aggregate(discharge, by=list(year), mean)
  names(mean.sampled.Q) <- c("year", "mean.sampled.Q")

  mean.overall.Q <- aggregate(hydrology$discharge, by=list(hydrology$year), mean)
  names(mean.overall.Q) <- c("year", "mean.overall.Q")

  correction.factor <- mean.overall.Q$mean.overall.Q/mean.sampled.Q$mean.sampled.Q


  yearly.from.method2 <- data.frame(year=loads.from.method1$year,
                                    yearly.load.tons=loads.from.method1$yearly.load.tons *
                                      correction.factor)
  return(yearly.from.method2)
}

# --> TS this should be done by the user before calling the function
# and should be loaded as binary into /data for the example
#reading and formating the input data
data <- read.table("data/retention_eff.csv", header=T, sep=",", dec=".")
# This should be done inside the function
data <- data[data$year>=start.year & data$year<=end.year,]
# data$date <- lubridate::mdy(data$date)
data$doy <- lubridate::yday(data$date)
data$year <- lubridate::year(data$date)
# <-- TS

# retention_eff <- function(data,method,start.year, end.year) {
  if(!method%in%c("a","b","c")) {
    stop("Method must be one of a, b, c")
  }
  my.summary.loads <- data.frame(NULL)
  my.inlets <- levels(as.factor(data$in_outlet))
  my.tributaries <- levels(as.factor(data$tributary))
  my.nutrients <- levels(as.factor(data$variable))
  my.nutrients <- my.nutrients[my.nutrients!="Q"]

  for (this.inflow in my.tributaries) {
    for(this.outlet in my.inlets){
      for (this.var in my.nutrients) {
        print(paste(this.inflow,this.var, sep=" & "))

        #select the required data to apply the load calculation methods (indexing)
        my.Q  <- data[data$tributary==this.inflow & data$variable=="Q" & data$in_outlet==this.outlet,]

        my.cq <- data[data$tributary==this.inflow & data$variable==this.var & data$in_outlet==this.outlet,]

        the.result <- data.frame(NULL)

        #we have to add the discharge Q to all measurements of the current variable in my.cq
        my.cq$Q <- my.Q$value[match(my.cq$date,my.Q$date)]
        if("GAM.load" %in% methods){
          res.GAM <- load.GAM(hydrology=data.frame(times=my.Q$date,
                                                   year=my.Q$year,
                                                   doy=my.Q$doy,
                                                   discharge=my.Q$value),
                              year=my.cq$year, doy=my.cq$doy, discharge=my.cq$Q,
                              concentration=my.cq$value, GOF=T)
          the.result <- data.frame(year = res.GAM$year, GAM.load = res.GAM$load.tons)
        }
        if("method1" %in% methods){
          res.method1 <- load.method1(year=my.cq$year, discharge=my.cq$Q, concentration=my.cq$value)
          the.result <- data.frame(year = res.method1$year, method1 = res.method1$yearly.load.tons)

        }
        if("method2" %in% methods){
          res.method2 <- load.method2(hydrology=data.frame(year=my.Q$year,
                                                           discharge=my.Q$value),
                                      year=my.cq$year, discharge=my.cq$Q, concentration=my.cq$value)
          the.result$method2 <- res.method2$yearly.load.tons[match(the.result$year, res.method2$year)]
        }

        #adding the other relevant information into the.result
        the.result$inflow   <- this.inflow
        the.result$variable <- this.var
        the.result$in_outlet <- this.outlet

        my.summary.loads <- rbind(my.summary.loads, the.result)

        my.summary.loads

      } #closing variable loop
    } #closing inflow loop
  } #closing tributary loop

  #calculating nutrient removal efficiency
  my.summary.loads$year <- as.numeric(my.summary.loads$year)
  efficiency <- data.frame(matrix(ncol = 3 + length(methods), nrow = 0))
  colnames(efficiency) <- c("year", "inflow", "variable", methods)

  for(var in levels(as.factor(my.summary.loads$variable))) {
    for(river in levels(as.factor(my.summary.loads$inflow))){
      for(yr in start.year:end.year){
        eff <- c()
        for(method in methods){
          print(paste(var,river,yr))

          relevant.data <- my.summary.loads[my.summary.loads$variable==var & my.summary.loads$inflow==river &
                                              my.summary.loads$year==yr,]

          inflow_value <- relevant.data[relevant.data$in_outlet=="inflow",][method]
          outflow_value <- relevant.data[relevant.data$in_outlet=="outflow",][method]

          eff1 <- (((inflow_value) - (outflow_value))/(inflow_value))
          eff <- append(eff, eff1)
        }
        fixos <- c(yr, river, var, eff)
        pinho <- as.data.frame(fixos)
        colnames(pinho) <- c("year", "inflow", "variable", methods)
        efficiency <- merge(efficiency, pinho, all.x = TRUE, all.y = TRUE)
        efficiency
      }
    }
  }

#   out <- data.frame(my.summary.loads,efficiency)
#   return(out)
# }
