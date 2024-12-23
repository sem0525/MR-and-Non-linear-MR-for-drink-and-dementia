#################sTale 5: PRS with assocition with X,Y and confoudners#######
######IV continous

###lifelong abstainer
###############96 SNP might includeding to the PRS
## 99 SNP from GWAS lefting single 96 SNP ##
#rs2532276 is missing in the UKB
#rs1260326 is related with fasting plasma glucose 
#rs28929474 is related with α-1-antitrypsin 
##To test 96 SNP wether related with covriables##
load(paste0(project_path,"data/drink_96snp.RData")) ##drink_96snp
drink_96snp_name <- names(drink_96snp)[2:97]
cof <- c("townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","CVD_before","stroke_before",
         "outcome","alcohol_liver_alcoholliver","age")
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==0) 
## all
Result_snp <- data.frame()
for(xvar in cof){
  covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
  modeldata$Var <- modeldata[[xvar]]
  for(snpi in drink_96snp_name){
    modeldata$IV <- modeldata[[snpi]]
    if(xvar == "age") {covar_f_name <- setdiff(covar_f_name,"age")}
    fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
    fit <- glm(as.formula(fun), data=modeldata)
    if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
      fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
    }
    r <- summary(fit)
    res <- as.data.frame(t(r$coefficients[2,]))
    res$se <- res$`Std. Error`
    res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                       sprintf("%0.3f",res$se),")")
    if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
    if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
    res_var <- cbind(xvar, snpi,res)
    if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_snp)}
    Result_snp <- rbind(Result_snp,res_var)
  }
}
Result_snp_all <- Result_snp
## subgroup sex
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==0 & sex==1) 
Result_snp <- data.frame()
for(xvar in cof){
  covar_f_name <- c("age","area","array",paste0("component",1:20))
  modeldata$Var <- modeldata[[xvar]]
  for(snpi in drink_96snp_name){
    modeldata$IV <- modeldata[[snpi]]
    if(xvar == "age") {covar_f_name <- setdiff(covar_f_name,"age")}
    fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
    fit <- glm(as.formula(fun), data=modeldata)
    if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
      fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
    }
    r <- summary(fit)
    res <- as.data.frame(t(r$coefficients[2,]))
    res$se <- res$`Std. Error`
    res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                       sprintf("%0.3f",res$se),")")
    if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){res$p <- p_out(res$`Pr(>|z|)`)}
    if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver"))){res$p <- p_out(res$`Pr(>|t|)`)}
    res_var <- cbind(xvar, snpi,res)
    if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){names(res_var) <- names(Result_snp)}
    Result_snp <- rbind(Result_snp,res_var)
  }
}
Result_snp_men <- Result_snp
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==0 & sex==0) 
Result_snp <- data.frame()
for(xvar in cof){
  covar_f_name <- c("age","area","array",paste0("component",1:20))
  modeldata$Var <- modeldata[[xvar]]
  for(snpi in drink_96snp_name){
    modeldata$IV <- modeldata[[snpi]]
    if(xvar == "age") {covar_f_name <- setdiff(covar_f_name,"age")}
    fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
    fit <- glm(as.formula(fun), data=modeldata)
    if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
      fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
    }
    r <- summary(fit)
    res <- as.data.frame(t(r$coefficients[2,]))
    res$se <- res$`Std. Error`
    res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                       sprintf("%0.3f",res$se),")")
    if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){res$p <- p_out(res$`Pr(>|z|)`)}
    if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver"))){res$p <- p_out(res$`Pr(>|t|)`)}
    res_var <- cbind(xvar, snpi,res)
    if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){names(res_var) <- names(Result_snp)}
    Result_snp <- rbind(Result_snp,res_var)
  }
}
Result_snp_women <- Result_snp
### combinded the result 
Result_snp_1 <- Result_snp_all %>% filter(`Pr(>|t|)` < 0.05/(96*(6+2)))  # 2 items
Result_snp_1$group = "all"
Result_snp_2 <- Result_snp_men %>% filter(`Pr(>|t|)` < 0.05/(96*(6+2)))   # 6 items
Result_snp_2$group = "men"
Result_snp_3 <- Result_snp_women %>% filter(`Pr(>|t|)` < 0.05/(96*(6+2)))  # 0 items
Result_snp_3$group = "women"
Result_snp_combind <- rbind(Result_snp_1,Result_snp_2)
Result_snp_combind$Beta <- sprintf("%0.2f", Result_snp_combind$Estimate)
Result_snp_combind$SE <- sprintf("%0.3f", Result_snp_combind$`Std. Error`)
Result_snp_combind$t <- sprintf("%0.2f", Result_snp_combind$`t value`)
Result_snp_combind$p <- sprintf("%0.8f", Result_snp_combind$`Pr(>|t|)`)
Result_snp_combind <- Result_snp_combind %>% select(group,xvar,snpi,Beta,SE,t,p)
sTable_snp_abstainer <- Result_snp_combind 


drink_95snp_name <- drink_96snp_name[!(drink_96snp_name %in% "rs13024996")]
drink_94snp_name_men <- drink_95snp_name[!(drink_95snp_name %in% "rs561222871")]

##rs561222871 specical for men
##data$wPRS_95 <- rowSums(data[,drink_95snpw_name],na.rm=T)



##95 SNP
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==0) 
modeldata$IV <- modeldata$wPRS_95
cof <- c("townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","CVD_before","stroke_before",
         "outcome","alcohol_liver_alcoholliver","age")

Result_q_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$se <- res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$se),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_q_continous)}
  Result_q_continous <- rbind(Result_q_continous,res_var)
}
Result_q_continous_abstainer <- Result_q_continous
test_p <- 0.05/length(cof)



###subgroup into sex
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
##men
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==0,sex==1) 
modeldata$IV <- modeldata$wPRS_95
Result_q_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$se <- res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$se),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_q_continous)}
  Result_q_continous <- rbind(Result_q_continous,res_var)
}
Result_q_continous_abstainer_men <- Result_q_continous
test_p_men <- 0.05/length(cof)
##women
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==0,sex==0) 
modeldata$IV <- modeldata$wPRS_95
Result_q_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$se <- res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$se),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_q_continous)}
  Result_q_continous <- rbind(Result_q_continous,res_var)
}
Result_q_continous_abstainer_women <- Result_q_continous
test_p_women <- 0.05/length(cof)


###combind all result
table_var  <- matrix(NA,nrow=(length(cof)+2),ncol=3)
table_var[,1] <- c("Potential confounders","  Socioeconomic","  Education","  BMI",
                   "  Smoke status","  Sleep duration","  Physical activity", 
                   "Cardiometabolic disease at baseline","Stroke at baseline",
                    "Outcome", "  Dementia","  Alcohol liver disease", "  Age")
table_var[2:(length(cof)-2),2] <- Result_q_continous_abstainer$beta[1:(length(cof)-3)]
table_var[2:(length(cof)-2),3] <- Result_q_continous_abstainer$p[1:(length(cof)-3)]

table_var[length(cof):(length(cof)+2),2] <- Result_q_continous_abstainer$beta[(length(cof)-2):length(cof)]
table_var[length(cof):(length(cof)+2),3] <- Result_q_continous_abstainer$p[(length(cof)-2):length(cof)]
table_var_all <- table_var
#men
table_var  <- matrix(NA,nrow=(length(cof)+2),ncol=3)
table_var[2:(length(cof)-2),2] <- Result_q_continous_abstainer_men$beta[1:(length(cof)-3)]
table_var[2:(length(cof)-2),3] <- Result_q_continous_abstainer_men$p[1:(length(cof)-3)]
table_var[length(cof):(length(cof)+2),2] <- Result_q_continous_abstainer_men$beta[(length(cof)-2):length(cof)]
table_var[length(cof):(length(cof)+2),3] <- Result_q_continous_abstainer_men$p[(length(cof)-2):length(cof)]
table_var_men <- table_var
#women
table_var  <- matrix(NA,nrow=(length(cof)+2),ncol=3)
table_var[2:(length(cof)-2),2] <- Result_q_continous_abstainer_women$beta[1:(length(cof)-3)]
table_var[2:(length(cof)-2),3] <- Result_q_continous_abstainer_women$p[1:(length(cof)-3)]
table_var[length(cof):(length(cof)+2),2] <- Result_q_continous_abstainer_women$beta[(length(cof)-2):length(cof)]
table_var[length(cof):(length(cof)+2),3] <- Result_q_continous_abstainer_women$p[(length(cof)-2):length(cof)]
table_var_women <- table_var

sTable_PRS_abstainer <- cbind(table_var_all,table_var_men,table_var_women)


##current drinker
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==2, perweek_alchol_unit_cal_new>0) 
modeldata$exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
cof <- c("exposure","perweek_alchol_unit_cal_new","townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","CVD_before","stroke_before",
         "outcome","alcohol_liver_alcoholliver","age")

Result_q_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$se <- res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$se),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_q_continous)}
  Result_q_continous <- rbind(Result_q_continous,res_var)
}
Result_q_continous_drinker <- Result_q_continous
test_p <- 0.05/length(cof)

###subgroup into sex
covar_f_name <- c("age","area","array",paste0("component",1:20))
##men
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==2,perweek_alchol_unit_cal_new>0, sex==1) 
modeldata$exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
modeldata$IV <- modeldata$wPRS_95
Result_q_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$se <- res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$se),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_q_continous)}
  Result_q_continous <- rbind(Result_q_continous,res_var)
}
Result_q_continous_drinker_men <- Result_q_continous
test_p_men <- 0.05/length(cof)
##women
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==2,perweek_alchol_unit_cal_new>0,sex==0) 
modeldata$IV <- modeldata$wPRS_95
modeldata$exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
Result_q_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$se <- res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$se),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver","CVD_before","stroke_before")){names(res_var) <- names(Result_q_continous)}
  Result_q_continous <- rbind(Result_q_continous,res_var)
}
Result_q_continous_drinker_women <- Result_q_continous
test_p_women <- 0.05/length(cof)


###combind all result
table_var  <- matrix(NA,nrow=(length(cof)+3),ncol=3)
table_var[,1] <- c("Alcohol consumption","  Log10(unit/week)","  unit/week","Potential confounders","  Socioeconomic","  Education","  BMI",
                   "  Smoke status","  Sleep duration","  Physical activity",
                   "Cardiometabolic disease at baseline","Stroke at baseline",
                   "Outcome", "  Dementia","  Alcohol liver disease", "  Age")
table_var[2:3,2] <- Result_q_continous_drinker$beta[1:2]
table_var[2:3,3] <- Result_q_continous_drinker$p[1:2]

table_var[5:(length(cof)-1),2] <- Result_q_continous_drinker$beta[3:(length(cof)-3)]
table_var[5:(length(cof)-1),3] <- Result_q_continous_drinker$p[3:(length(cof)-3)]

table_var[(length(cof)+1):(length(cof)+3),2] <- Result_q_continous_drinker$beta[(length(cof)-2):length(cof)]
table_var[(length(cof)+1):(length(cof)+3),3] <- Result_q_continous_drinker$p[(length(cof)-2):length(cof)]
table_var_all <- table_var
#men
table_var  <- matrix(NA,nrow=(length(cof)+3),ncol=3)
table_var[2:3,2] <- Result_q_continous_drinker_men$beta[1:2]
table_var[2:3,3] <- Result_q_continous_drinker_men$p[1:2]

table_var[5:(length(cof)-1),2] <- Result_q_continous_drinker_men$beta[3:(length(cof)-3)]
table_var[5:(length(cof)-1),3] <- Result_q_continous_drinker_men$p[3:(length(cof)-3)]

table_var[(length(cof)+1):(length(cof)+3),2] <- Result_q_continous_drinker_men$beta[(length(cof)-2):length(cof)]
table_var[(length(cof)+1):(length(cof)+3),3] <- Result_q_continous_drinker_men$p[(length(cof)-2):length(cof)]
table_var_men <- table_var
#women
table_var  <- matrix(NA,nrow=(length(cof)+3),ncol=3)
table_var[2:3,2] <- Result_q_continous_drinker_women$beta[1:2]
table_var[2:3,3] <- Result_q_continous_drinker_women$p[1:2]

table_var[5:(length(cof)-1),2] <- Result_q_continous_drinker_women$beta[3:(length(cof)-3)]
table_var[5:(length(cof)-1),3] <- Result_q_continous_drinker_women$p[3:(length(cof)-3)]

table_var[(length(cof)+1):(length(cof)+3),2] <- Result_q_continous_drinker_women$beta[(length(cof)-2):length(cof)]
table_var[(length(cof)+1):(length(cof)+3),3] <- Result_q_continous_drinker_women$p[(length(cof)-2):length(cof)]
table_var_women <- table_var

sTable_PRS_drinker <- cbind(table_var_all,table_var_men,table_var_women)
