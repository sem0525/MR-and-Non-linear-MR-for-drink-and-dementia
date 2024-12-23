###############load and merge data############
#####load pakages####
library(dplyr)
library(ggplot2)
library(readxl)
library(ggsci)
library(gmodels)
library(survival)
#############-------------------load data----------------#############################################
###check the population flow chart
#load the cleaned data
####genetic data
path <- "/home/ukb/zoo_project/UKB/data_clean/"
load(paste0(path,"drink_98snp_w.RData"))  #drink_98snp_w 487412 
load(paste0(path,"drink_98snp.RData"))  #drink_98snp 487412      
load(paste0(path,"genetic_qc.RData"))  #genetic_qc 502370
load(paste0(path,"APOE.Rdata"))  #APOE 487409
####base
load(paste0(path,"base_data.Rdata"))  #base_data 502370 
load(paste0(path,"Data.RData"))  #Data 502370 #lifestyle
load(paste0(path,"sensor_date.RData"))  #sensor_date 502370 

###alcohol
load(paste0(path,"alcohol_data.RData"))  #alcohol_data 502370 
###disease
load(paste0(path,"depression.RData"))  #depression 502370 
load(paste0(path,"dementia_data.RData"))  #dementia_data 502370 
load(paste0(path,"dementia_all.RData"))  #dementia_all 502370 
load(paste0(path,"alcohol_liver_all.RData"))  #alcohol_liver_all 502370 
load(paste0(path,"FI.RData")) #FI 502370
load(paste0(path,"CVD.RData"))  #CVD 502411 


Data <- Data %>% select(-alcohlol_frequency)#502370 
data <- merge(sensor_date,base_data,by="eid") #502370 
data <- merge(data,Data,by="eid") #502370     
data <- merge(data,alcohol_data,by="eid") #502370     
data <- merge(data,depression,by="eid") #502370     
data <- merge(data,dementia_all,by="eid") #502370
data <- merge(data,alcohol_liver_all,by="eid")  #502370
data <- merge(data,dt_CVD,by="eid")  #502370

###combined fluid intelligence score
data <- merge(data,FI,by="eid",all.x=TRUE) 
table(is.na(data$FIscore)) #117230 

##genetic data
data <- merge(data,genetic_qc,by="eid")   #502370
data <- merge(data,APOE,by="eid",all.x=TRUE) #502370 #482408  ##missing 14175 
table(is.na(data$APOE))  ##missing 14175 
prop.table(table(is.na(data$APOE)))

table(is.na(data$APOE),data$ACD)  ##missing 14175 
prop.table(table(is.na(data$APOE)))
###drink_98snp
#Removing the rs1260326 and rs28929474
#Check each SNP’s information in the dbSNP (www.ncbi.nlm.hih.gov/snp)
#rs1260326 is related with fasting plasma glucose level quantitative
#rs28929474 is related with alpha-1-antitrypsin deficiency (likely-pathogenic)
#rs1229984 is related with alcohol dependence 
drink_96snp <- drink_98snp %>% select(-rs1260326,-rs28929474)
drink_96snp$drink_snp <- 1
drink_96snp_w <- drink_98snp_w %>% select(-rs1260326_w,-rs28929474_w)
####save 96 SNP data####
save_path <- "/home/zoo/project/P1_drinkMR_dementia/"
save(drink_96snp,file=paste0(save_path,"data/drink_96snp.RData"))
save(drink_96snp_w,file=paste0(save_path,"data/drink_96snp_w.RData"))


data <- merge(data,drink_96snp,by="eid",all.x=TRUE)
data <- merge(data,drink_96snp_w,by="eid",all.x=TRUE)
table(is.na(data$drink_snp))
prop.table(table(is.na(data$drink_snp)))


#####clean data
################clean data: outcome,age,edu_new,APOE4_status################
#age
data$age_group <- ifelse(data$age<=45,0,
                         ifelse(data$age<=65,1,2))
data$age_2group <- ifelse(data$age<=65,1,2)
###edu
#other: CES or equivalent, none of above  == 5
#Vocational: NVQ or HND or HNC or equivalent  == 4
#lower secondary (first stage of secondary education) == 3
#upper secondary (second/final stage of secondary educaiton) == 2
#high (college/university degree or other professional qualificaiton) == 1
data$edu_new <- data$edu
data$edu_new[data$edu==4] <- 1
data$edu_new[data$edu==3] <- 2
data$edu_new[data$edu==5] <- 3
##missing to other
data$edu_new[is.na(data$edu_new)==TRUE] <- 3 ##missing is other ##edu_new has 3771 missing
data$APOE4_status <- ifelse(data$APOE %in% c("e3e4","e4e4"),1,0)
data$APOE4 <- ifelse(data$APOE=="e4e4",2,
                     ifelse(data$APOE=="e3e4",1,0))
table(data$APOE4,data$APOE4_status )


###disease for alcoholliver
data$outcome <- data$alcohol_liver_alcoholliver
data$last_date <- data$sensor_date
data$last_date[data$outcome==1]  <- data$date_disease_first_alcoholliver[data$outcome==1]
data$last_date <- as.Date(data$last_date)
data$basedate <- as.Date(data$basedate)
data$days_alcoholliver <- data$last_date - data$basedate
data$years_alcoholliver <- round(data$days_alcoholliver/365,1)
data$years_alcoholliver <- as.numeric(data$years_alcoholliver)

###disease
data$last_date <- data$sensor_date
data$outcome <- data$ACD
data$last_date[data$outcome==1]  <- data$date_ACD_first[data$outcome==1]
data$last_date <- as.Date(data$last_date)
data$basedate <- as.Date(data$basedate)
data$days <- data$last_date - data$basedate
data$years <- round(data$days/365,1)
data$years <- as.numeric(data$years)

##CVD
table(is.na(data$min_CVD))
data$CVD <- ifelse(is.na(data$min_CVD)==F,1,0)
data$min_CVD <- as.Date(data$min_CVD)
data$CVD_before <- data$min_CVD - data$basedate
data$CVD_before <- ifelse(data$CVD_before<=0,1,0)
table(data$CVD_before)
table(is.na(data$CVD_before))
data$CVD_before[is.na(data$CVD_before)==T] <- 0
##stroke
data$Stroke <- ifelse(is.na(data$stroke)==F,1,0)
data$stroke <- as.Date(data$stroke)
data$stroke_before <- data$stroke - data$basedate
data$stroke_before <- ifelse(data$stroke_before<=0,1,0)
table(data$stroke_before)
table(is.na(data$stroke_before))
data$stroke_before[is.na(data$stroke_before)==T] <- 0

##save data
DATA_save <- data
save_path <- "/home/zoo/project/P1_drinkMR_dementia/"
save(DATA_save,file=paste0(save_path,"data/DATA_save_502370.RData"))






