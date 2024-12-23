#################sTale: PRS with assocition with X,Y and confoudners#######
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)
DATA <- data_all_analysis
######IV continous
modeldata <- DATA 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
cof <- c("Exposure","Exposure_log","townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","CVD_before","stroke_before",
         "outcome","alcohol_liver_alcoholliver","age")
Result_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  if(xvar=="age"){fun <- paste("Var ~ IV + ",paste(covar_f_name[-1],collapse=" + "))}
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$low),", ",
                     sprintf("%0.3f",res$up),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_continous)}
  Result_continous <- rbind(Result_continous,res_var)
}
Result_continous$test <- ifelse(Result_continous$`Pr(>|t|)`<(0.05/11),1,0)
Result_continous_overall <- Result_continous



##########men################
modeldata <- DATA %>% filter(sex==1)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
cof <- c("Exposure","Exposure_log","townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","CVD_before","stroke_before"
         ,"outcome","alcohol_liver_alcoholliver","age")
Result_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  if(xvar=="age"){fun <- paste("Var ~ IV + ",paste(covar_f_name[-1],collapse=" + "))}
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$low),", ",
                     sprintf("%0.3f",res$up),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_continous)}
  Result_continous <- rbind(Result_continous,res_var)
}
Result_continous$test <- ifelse(Result_continous$`Pr(>|t|)`<(0.05/11),1,0)

Result_continous_men <- Result_continous


##########women################
modeldata <- DATA %>% filter(sex==0)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
cof <- c("Exposure","Exposure_log","townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","CVD_before","stroke_before",
         "outcome","alcohol_liver_alcoholliver","age")
Result_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  if(xvar=="age"){fun <- paste("Var ~ IV + ",paste(covar_f_name[-1],collapse=" + "))}
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$low),", ",
                     sprintf("%0.3f",res$up),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_continous)}
  Result_continous <- rbind(Result_continous,res_var)
}
Result_continous$test <- ifelse(Result_continous$`Pr(>|t|)`<(0.05/11),1,0)

Result_continous_women <- Result_continous

Result_continous_blank <- Result_continous_overall[1,]
Result_continous_blank[1,] <- NA
Result_continous_blank[1,1] <- "Men"

Result_continous_all <- rbind(Result_continous_overall,Result_continous_blank,
                              Result_continous_men)
Result_continous_blank[1,1] <- "Women"

Result_continous_all <- rbind(Result_continous_all,Result_continous_blank,Result_continous_women)
