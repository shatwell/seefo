# Save example dataset of temperature and ice (FLake model output for Lake Stechlin)

library(data.table)
Ts_Tb_ice <- read.table("C:/Users/shatwell/Documents/Modellierung/FLake/flake2_test/ascii_output/flake2_test/Stechlin/flake2TWS_EWEMBI_historical_Stechlin.rslt", skip=1, header=TRUE, dec=".")
setDT(Ts_Tb_ice)
Ts_Tb_ice <- Ts_Tb_ice[time<=365*10+1]
Ts_Tb_ice[,date:=as.POSIXct("1997-01-01", tz="UTC")+time*3600*24]
Ts_Tb_ice <- as.data.frame(Ts_Tb_ice[,.(date, Ts, Tb, H_ice)])
write.csv(Ts_Tb_ice,"data-raw/Ts_Tb_ice.csv", row.names=FALSE)
usethis::use_data(Ts_Tb_ice, overwrite = TRUE)
