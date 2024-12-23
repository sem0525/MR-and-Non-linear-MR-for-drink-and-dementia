########################Table 2: linear MR analysis ##################
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)
DATA <- data_all_analysis
DATA$IV <- DATA$wPRS_95
DATA$Exposure <- log10(DATA$perweek_alchol_drink_cal_new)
#DATA$Exposure <- DATA$perweek_alchol_drink_cal_new

##all

a = 2 #exposure: unit alcohol per week
population <- "Overall"
subgroup <- "all"
IV_name <- "weighted 95 SNP"

##orignial
##model 1
modeldata <- DATA
x="Model1"
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
re_1 <- GS_mr_one_HR(modeldata,covar_f_name,x)
r1 <- re_1 %>% select(x,HR,p_value)
##competing risk model
library(cmprsk)
x = "Competing Model"
modeldata$death <- ifelse(is.na(modeldata$death_date)==FALSE,1,0)
modeldata$outcome_new <- modeldata$outcome
modeldata$outcome_new[modeldata$death==1] <- 2
modeldata$area <- ifelse(modeldata$area=="England",1,ifelse(modeldata$area=="Scotland",2,3))
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
##descption how many death before dementia
table(modeldata$outcome,modeldata$death)
CrossTable(modeldata$sex,modeldata$death)
CrossTable(modeldata$outcome,modeldata$death)

modeldata$death_date <- as.Date( modeldata$death_date)
modeldata$death_before_dementia_years <- modeldata$death_date - modeldata$date_ACD_first
table(modeldata$death_before_dementia_years)
re_2 <- IV_HR_function_compete(modeldata,covar_f_name,x)
r2 <- re_2 %>% select(x,HR,p_value)


DATA$comorbid <- ifelse(DATA$CVD_before == 1 | DATA$stroke_before==1,1,0)
table(DATA$comorbid ,DATA$CVD_before )
table(DATA$comorbid ,DATA$stroke_before )
CrossTable(DATA$comorbid ,DATA$outcome )

modeldata <- DATA %>% filter(comorbid==1)
re_sub <- GS_mr_one_HR(modeldata,covar_f_name,x)
re_sub1 <- re_sub %>% select(x,HR,p_value)
modeldata <- DATA %>% filter(comorbid==0)
re_sub <- GS_mr_one_HR(modeldata,covar_f_name,x)
re_sub2 <- re_sub %>% select(x,HR,p_value)
re_sub1
re_sub2


###sumamry-level-analysis



