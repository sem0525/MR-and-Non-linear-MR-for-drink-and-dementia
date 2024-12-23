"""
Whether one had records in cardiometabolic diseases at baseline calculated based 
on UK Biobank first occurrence or algorithmically-defined outcomes: 
myocardial infarction (Field 42000), heart failure (Field 131354), 
hypertension (Field 131286, 131288, 131290, 131292, 131294), 
type 2 diabetes mellitus (Field 130706, 130710,130712,130714)

"""
savefile = "~/project/P0_UKB_extractdata/UKB_data/UKB_RDdata/DT/"
load(paste0(savefile,"dt_00",4,"_fo.RData"))
# hear failure
dt_CVD1 <- dt_fo %>% select(
  eid,
  cvd1 =starts_with("131354"),
  cvd2 = starts_with("131286"),
  cvd3 = starts_with("131288"),
  cvd4 = starts_with("131290"),
  cvd5 = starts_with("131292"),
  cvd6 = starts_with("131294"),
  
  
)
load(paste0(savefile,"dt_00",3,"_fo.RData"))
dt_CVD2 <- dt_fo %>% select(
  eid,
  cvd7 = starts_with("130706"),
  cvd8 = starts_with("130710"),
  cvd9 = starts_with("130712"),
  cvd10 = starts_with("130714"),
  
)

dt_CVD <- merge(dt_CVD1,dt_CVD2,by="eid")


##
load(paste0(savefile,"dt",53,".RData"))
dt_CVD3 <- dt %>% select(
  eid,
  cvd11 = starts_with("42000"))
dt_CVD <- merge(dt_CVD,dt_CVD3,by="eid")

dt_CVD$min_CVD <- do.call(pmin, c(dt_CVD[,2:12], na.rm=T))
fname = "dt_CVD"
save(dt_CVD, file=paste0(savefile,fname,".RData")) #dt_CVD
