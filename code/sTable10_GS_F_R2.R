########################sTable 2: linear MR analysis R2 for GS##################
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)
DATA <- data_all_analysis

###log transfrom###
#########supplementary table R2 for PRS########
##all
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
modeldata$Exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","sex","area","array",paste0("component",1:40))
x = "Overall"
R2_overall <- R2_F_output(modeldata,x)

##Men
modeldata <- DATA %>% filter(sex==1)
modeldata$IV <- modeldata$wPRS_95
modeldata$Exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","area","array",paste0("component",1:20))
x = "Men"
R2_men <- R2_F_output(modeldata,x)

##Women
modeldata <- DATA %>% filter(sex==0)
modeldata$IV <- modeldata$wPRS_95
modeldata$Exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","area","array",paste0("component",1:20))
x = "Women"
R2_women <- R2_F_output(modeldata,x)

R2_log <- cbind(t(R2_overall),t(R2_men),t(R2_women))



#########supplementary table R2 for PRS########
##all
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
covar_f_name <- c("age","sex","area","array",paste0("component",1:40))
x = "Overall"
R2_overall <- R2_F_output(modeldata,x)

##Men
modeldata <- DATA %>% filter(sex==1)
modeldata$IV <- modeldata$wPRS_95
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
covar_f_name <- c("age","area","array",paste0("component",1:20))
x = "Men"
R2_men <- R2_F_output(modeldata,x)

##Women
modeldata <- DATA %>% filter(sex==0)
modeldata$IV <- modeldata$wPRS_95
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
covar_f_name <- c("age","area","array",paste0("component",1:20))
x = "Women"
R2_women <- R2_F_output(modeldata,x)

R2_x <- cbind(t(R2_overall),t(R2_men),t(R2_women))


All_R2 <- data.frame(rbind(R2_log,R2_x))


