#################sTable: Q of PRS with associition exposure and outcome################
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)

DATA <- data_all_analysis

modeldata <- DATA 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
g=5
##95 SNP
modeldata$IV <- modeldata$wPRS_95
modeldata$IV_group <- cut(modeldata$IV,
                          breaks=quantile(modeldata$IV,probs= seq(0,1,0.2)))
table(modeldata$IV_group)
modeldata$group <- as.numeric(modeldata$IV_group)
g = dim(table(modeldata$IV_group))
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
cof <- c("Exposure","Exposure_log","townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","CVD_before","stroke_before",
         "outcome","alcohol_liver_alcoholliver","age")
Result_q <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV_group + ",paste(covar_f_name,collapse=" + "))
  if(xvar == "age"){
    fun <- paste("Var ~ IV_group + ",paste(covar_f_name[-1],collapse=" + "))
  }
  
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(r$coefficients[1:5,1:2])
  
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  if(xvar %in% c("Exposure","Exposure_log")){res$beta <- c("Ref",paste0(sprintf("%0.2f",res[2:5,1])," (",
                                                                        sprintf("%0.2f",res[2:5,3]),", ",
                                                                        sprintf("%0.2f",res[2:5,4]),")"        
  ))}
  if(!(xvar %in% c("Exposure","Exposure_log"))){res$beta <- c("Ref",paste0(sprintf("%0.3f",res[2:5,1])," (",
                                                                           sprintf("%0.3f",res[2:5,3]),", ",
                                                                           sprintf("%0.3f",res[2:5,4]),")"        
  ))}
  
  #p for trend
  fun <- "Var ~ group "
  for(o in covar_f_name){ fun <- paste(fun,"+",o)}
  if(xvar == "age"){
    fun <- paste("Var ~ group + ",paste(covar_f_name[-1],collapse=" + "))
  }
  
  fit1 <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit1 <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  anoval_result <- anova(fit,fit1,test="Chisq")
  p <- p_out(anoval_result$`Pr(>Chi)`[2])
  res_var <- c(xvar, res$beta,p)
  Result_q <- rbind(Result_q,res_var)
}
Result_q
names(Result_q) <- c("VAR",paste0("Q",1:g),"P for trend")

Result_q_all <- Result_q

##########men################
modeldata <- DATA %>% filter(sex==1)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
modeldata$IV_group <- cut(modeldata$IV,
                          breaks=quantile(modeldata$IV,probs= seq(0,1,0.2)))
table(modeldata$IV_group)
modeldata$group <- as.numeric(modeldata$IV_group)
g = dim(table(modeldata$IV_group))
covar_f_name <- c("age","area","array",paste0("component",1:20))
cof <- c("Exposure","Exposure_log","townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group"
         ,"CVD_before","stroke_before","outcome",
         "alcohol_liver_alcoholliver","age")
Result_q_men <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- "Var ~ IV_group "
  for(o in covar_f_name){ fun <- paste(fun,"+",o)}
  if(xvar == "age"){
    fun <- paste("Var ~ IV_group + ",paste(covar_f_name[-1],collapse=" + "))
  }
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(r$coefficients[1:5,1:2])
  
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  if(xvar %in% c("Exposure","Exposure_log")){res$beta <- c("Ref",paste0(sprintf("%0.2f",res[2:5,1])," (",
                                                                        sprintf("%0.2f",res[2:5,3]),", ",
                                                                        sprintf("%0.2f",res[2:5,4]),")"        
  ))}
  if(!(xvar %in% c("Exposure","Exposure_log"))){res$beta <- c("Ref",paste0(sprintf("%0.3f",res[2:5,1])," (",
                                                                           sprintf("%0.3f",res[2:5,3]),", ",
                                                                           sprintf("%0.3f",res[2:5,4]),")"        
  ))}
  
  #p for trend
  fun <- "Var ~ group "
  for(o in covar_f_name){ fun <- paste(fun,"+",o)}
  if(xvar == "age"){
    fun <- paste("Var ~ group + ",paste(covar_f_name[-1],collapse=" + "))
  }
  fit1 <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit1 <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  anoval_result <- anova(fit,fit1,test="Chisq")
  p <- p_out(anoval_result$`Pr(>Chi)`[2])
  res_var <- c(xvar, res$beta,p)
  Result_q_men <- rbind(Result_q_men,res_var)
}

Result_q_men
names(Result_q_men) <- c("VAR",paste0("Q",1:g),"P for trend")


##########women################
modeldata <- DATA %>% filter(sex==0)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
modeldata$IV_group <- cut(modeldata$IV,
                          breaks=quantile(modeldata$IV,probs= seq(0,1,0.2)))
table(modeldata$IV_group)
modeldata$group <- as.numeric(modeldata$IV_group)
g = dim(table(modeldata$IV_group))
covar_f_name <- c("age","area","array",paste0("component",1:20))
cof <- c("Exposure","Exposure_log","townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","APOE4","CVD_before","stroke_before",
         "outcome", "alcohol_liver_alcoholliver","age")
Result_q_women <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- "Var ~ IV_group "
  for(o in covar_f_name){ fun <- paste(fun,"+",o)}
  if(xvar == "age"){
    fun <- paste("Var ~ IV_group + ",paste(covar_f_name[-1],collapse=" + "))
  }
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(r$coefficients[1:5,1:2])
  
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  if(xvar %in% c("Exposure","Exposure_log")){res$beta <- c("Ref",paste0(sprintf("%0.2f",res[2:5,1])," (",
                                                                        sprintf("%0.2f",res[2:5,3]),", ",
                                                                        sprintf("%0.2f",res[2:5,4]),")"        
  ))}
  if(!(xvar %in% c("Exposure","Exposure_log"))){res$beta <- c("Ref",paste0(sprintf("%0.3f",res[2:5,1])," (",
                                                                           sprintf("%0.3f",res[2:5,3]),", ",
                                                                           sprintf("%0.3f",res[2:5,4]),")"        
  ))}
  
  #p for trend
  fun <- "Var ~ group "
  for(o in covar_f_name){ fun <- paste(fun,"+",o)}
  if(xvar == "age"){
    fun <- paste("Var ~ group + ",paste(covar_f_name[-1],collapse=" + "))
  }
  fit1 <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit1 <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  anoval_result <- anova(fit,fit1,test="Chisq")
  p <- p_out(anoval_result$`Pr(>Chi)`[2])
  res_var <- c(xvar, res$beta,p)
  Result_q_women <- rbind(Result_q_women,res_var)
}

Result_q_women
names(Result_q_women) <- c("VAR",paste0("Q",1:g),"P for trend")
Overall <- c("Overall",rep(NA,6))
Result_q <- rbind(Overall,Result_q)

Overall <- c("Men",rep(NA,6))
Result_q_men <- rbind(Overall,Result_q_men)

Overall <- c("Women",rep(NA,6))
Result_q_women <- rbind(Overall,Result_q_women)

Result_q_all <- rbind(Result_q,Result_q_men,Result_q_women)
