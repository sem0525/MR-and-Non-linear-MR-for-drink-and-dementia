library(survival)
library(gmodels)
######snp association test######
association_test <- function(y,x,dt){
  d <- as.data.frame(dt)
  d$var <- d[,x]
  d$conf <- d[,y]
  fun <- "conf ~ var "
  for(o in covar_f_name){ 
    fun <- paste(fun,"+",o)
  }
  fit <- summary(lm(as.formula(fun), data=d))
  var_res <- fit$coefficients[2,]
  var_Res <- c(y,x,var_res)
  var_Res
}

######snp test output save######
out_save <- function(i,keep=3,dt){
  dt[,i] <- as.numeric(dt[,i])
  dt[,i] <- sprintf(paste0("%0.",keep,"f"),dt[,i])
  dt
}

######figure function#####
hinvert_title_grob <- function(grob){
  # Swap the widths
  widths <- grob$widths
  grob$widths[1] <- widths[3]
  grob$widths[3] <- widths[1]
  grob$vp[[1]]$layout$widths[1] <- widths[3]
  grob$vp[[1]]$layout$widths[3] <- widths[1]
  # Fix the justification
  grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
  grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
  grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
  grob
}

######suplement tabe 1####
table_1_description_by_group <- function(conf,tdata,ng=1,pre=FALSE){
  tdata <- as.data.frame(tdata)
  N <- data.frame(matrix(NA,nrow=1,ncol=ng+1))
  N[1,1] <- "No of participant"
  N_conf <- list()
  for(a in 1:ng){
    N[1,(a+1)] <- paste0(dim(tdata[tdata$group==a,])[1]," (",
                         sprintf("%0.2f",dim(tdata[tdata$group==a,])[1]/dim(tdata)[1]*100), ")")
    n_conf <- data.frame()
    for(i in conf){
      tdata$var <- NA
      tdata[,dim(tdata)[2]] <- tdata[,i]
      #continous - normal
      if(i %in% c("years","age","townsend","BMI","prs","prs_w")){
        n <- matrix(NA,nrow=1,ncol=2)
        if(i=="years"){n[1,1] <- "Follow-up years, year"}
        if(i=="age"){n[1,1] <- "Age, years"}
        if(i=="townsend"){n[1,1] <- "Townsend deprivation index"}
        if(i=="BMI"){n[1,1] <- "BMI, kg/m^2"}
        if(i=="prs"){n[1,1] <- "PRS for alcohol"}
        if(i=="prs_w"){n[1,1] <- "Weight PRS for alcohol"}
        
        n[1,2] <- paste0(sprintf("%0.1f",mean(tdata$var[tdata$group==a],na.rm=TRUE)), " ± ",
                         sprintf("%0.1f",sd(tdata$var[tdata$group==a],na.rm=TRUE)))
        n_conf <- rbind(n_conf,n)
      }
      
      #continous - unnormal perweek_alchol_unit_cal_new
      if(i %in% c("edu_y","incomesocre",
                  "permonth_alchol_drink_cal","perweek_alchol_unit_cal",
                  "perweek_alchol_drink_cal_new","perweek_alchol_unit_cal_new",
                  "permonth_alchol_drink_cal","perweek_alchol_unit_cal_new",
                  "perweek_alchol_g_cal_new","perweek_alchol_g_cal"
                  ,"prs","prs_w"
      )){
        n <- matrix(NA,nrow=1,ncol=2)
        if(i=="edu_y"){n[1,1] <- "Education year, year"}
        if(i=="incomesocre"){n[1,1] <- "Income score"}
        if(i=="perweek_alchol_drink_cal"){n[1,1] <- "Alcohol consumption, drink per week"}
        if(i=="perweek_alchol_unit_cal"){n[1,1] <- "Alcohol consumption, unit per week"}
        if(i=="perweek_alchol_drink_cal_new"){n[1,1] <- "Alcohol consumption, drink per week"}
        if(i=="perweek_alchol_unit_cal_new"){n[1,1] <- "Alcohol consumption, unit per week"}
        if(i=="permonth_alchol_drink_cal"){n[1,1] <- "Alcohol consumption, drink per month"}
        if(i=="permonth_alchol_unit_cal"){n[1,1] <- "Alcohol consumption, unit per month"}
        if(i=="perweek_alchol_g_cal"){n[1,1] <- "Alcohol consumption, g per week"}
        if(i=="perweek_alchol_g_cal_new"){n[1,1] <- "Alcohol consumption, g per week"}
        if(i=="prs"){n[1,1] <- "PRS for alcohol"}
        if(i=="prs_w"){n[1,1] <- "Weight PRS for alcohol"}
        s_var <- summary(tdata$var[tdata$group==a])
        n[1,2] <- paste0(sprintf("%0.2f",s_var[3])," [",sprintf("%0.2f",s_var[2]),", ",sprintf("%0.2f",s_var[5]),"]")
        n_conf <- rbind(n_conf,n)
        
      }
      #category for multi-categroy
      if(i %in% c("age_group","age_2group","sex","area","edu","edu_new","edu_university","townsend_group_1",
                  "houseincome","smoking","smoke","BMI_group",
                  "PA_group","sleep_group",
                  "Diabetes","depression_g","APOE4_status",
                  "alcohlol_frequency","alcohol_drinker_status",
                  "family_alzheimer",
                  "Diabetes","Hypertension","depression_g","APOE4_status",
                  "alcohol_6group","alcohol_2group","ACD"
      )){
        n_dim <- dim(table(tdata$var))+1
        n <- matrix(NA,nrow=n_dim,ncol=2)
        s_table <- CrossTable(tdata$var[tdata$group==a])
        n[2:n_dim,2] <- paste0(s_table$t," (",sprintf("%0.1f",s_table$prop.row*100),")")
        
        
        if(i=="age_group"){n[,1] <- c("Age, n (%)",
                                      "  <=45 year","  (45, 65] year","  > 65 year")}
        if(i=="age_2group"){n[,1] <- c("Age, n (%)",
                                       "  <=65 year", "  > 65 year")}
        if(i=="sex"){n[,1] <- c("Gender, n (%)",
                                "  Women","  Men")}
        if(i=="area"){n[,1] <- c("Area, n (%)",
                                 "  England",
                                 "  Scotland",
                                 "  Wales")}
        if(i=="BMI_group"){n[,1] <- c("BMI , n (%)",
                                      "  <18.5  kg/m^2",
                                      "  [18.5~25.0) kg/m^2",
                                      "  [25.0~30.0) kg/m^2",
                                      "  >=30.0 kg/m^2")}
        if(i=="edu"){n[,1] <- c("Education, n (%)",
                                "  Higher",
                                "  Upper secondary",
                                "  Lower secondary",
                                "  Vocational",
                                "  Other")}
        if(i=="edu_new"){n[,1] <- c("Education, n (%)",
                                    "  Higher or vocational",
                                    "  Upper or lower secondary",
                                    "  Other")}
        if(i=="edu_university"){n[,1] <- c("College Degree, n (%)",
                                           "  Below university","  University or above")}
        if(i=="townsend_group_1"){n[,1] <- c("Socioeconomic, n (%)",
                                             "  Least deprived",
                                             "  Middle deprived",
                                             "  Most deprived)")}
        if(i=="houseincome"){n[,1] <- c("House income, n (%)",
                                        "  Less than 18,000",
                                        "  18,000 to 30,999",
                                        "  31,000 to 51,999",
                                        "  52,000 to 100,000",
                                        "  Greater than 100,000")}
        #smoking
        if(i=="smoking"){n[,1] <- c("Smoking, n (%)",
                                    "  No",
                                    "  Yes")}
        if(i=="smoke"){n[,1] <- c("Smoke status, n (%)",
                                  "  Never",
                                  "  Previous",
                                  "  Current")}
        if(i=="PA_group"){n[,1] <- c("Physical Activity, n (%)",
                                     "  Insufficient","  Sufficient","  additional")}
        if(i=="sleep_group"){n[,1] <- c("Sleep Duration, n (%)",
                                        "  <6 hour","  [6-9] hour","  >9 hour")}
        if(pre==TRUE & i=="alcohol_drinker_status"){n[,1] <- c("Aclohol drinking status, n (%)",
                                                               "  Never",
                                                               "  Previous",
                                                               "  Current")}
        if(pre==FALSE & i=="alcohol_drinker_status"){n[,1] <- c("Aclohol drinking status, n (%)",
                                                                "  Never",
                                                                "  Current")}
        if(i=="alcohlol_frequency"){n[,1] <- c("Alcohol Intake Frequency, n (%)",
                                               "  Never","  Special occasions only",
                                               "  One to three times a month","  Once or twice a week",
                                               "  Three or four times a week", "  Daily or almost daily")}
        if(i=="family_alzheimer"){n[,1] <- c("Family with Alzheimer's disease, n (%)","  Without","  With")}
        if(i=="Diabetes"){n[,1] <- c("Diabetes, n (%)","  Without","  With")}
        if(i=="Hypertension"){n[,1] <- c("Hypertension, n (%)","  Without","  With")}
        if(i=="depression_g"){n[,1] <- c("Depressive symptoms in last 2 weeks, n (%)","  Without","  With")}
        
        if(i=="APOE4_status"){n[,1] <- c("APOE e4, n (%)","  Without","  With")}
        if(i=="alcohol_2group"){n[,1] <- c("Alcohol consumption, n (%)",
                                           "  Normal (<=14 unit/week)",
                                           "  Exceed (>14 unit/week)")}
        if(i=="alcohol_6group"){n[,1] <- c("Alcohol consumption, n (%)",
                                           "  Light (<=7 unit/week)",
                                           "  Moderate ((7, 14] unit/week)",
                                           "  Heavy ((14, 21] unit/week)",
                                           "  Very heavy ((21, 28] unit/week)",
                                           "  Excessive ((28, 42] unit/week)",
                                           "  Extremely excessive (>42 unit/week)")}
        ##disease
        if(i=="ACD"){n[,1] <- c("All-cause dementia, n (%)","  Without","  With")}
        
        n_conf <- rbind(n_conf,n)
      }
      
    }
    N_conf[[a]] <- n_conf
  }
  N_conf_0 <- data.frame(N_conf[[1]])
  if(ng > 1) {for(a in 2:ng){
    subd <- N_conf[[a]]
    N_conf_0 <- cbind(N_conf_0,subd[,2])}
  }
  names(N_conf_0) <- names(N)
  N <- rbind(N,N_conf_0)
  N
}

######table 1######
table_1_incidence <- function(conf,tdata,pre=FALSE){
  tdata <- as.data.frame(tdata)
  N <- data.frame()
  N[1,1] <- "No of participant"
  N[1,2] <- dim(tdata)[1]
  N[1,3] <- paste0(dim(tdata[tdata$outcome==1,])[1]," (",
                   sprintf("%0.2f",dim(tdata[tdata$outcome==1,])[1]/dim(tdata)[1]*100), ")")
  N[1,4] <- NA
  for(i in conf){
    tdata$var <- NA
    tdata[,dim(tdata)[2]] <- tdata[,i]
    #continous - normal
    if(i %in% c("years","age","townsend","BMI","prs","prs_w")){
      n <- matrix(NA,nrow=1,ncol=4)
      if(i=="years"){n[1,1] <- "Follow-up years, year"}
      if(i=="age"){n[1,1] <- "Age, years"}
      if(i=="townsend"){n[1,1] <- "Townsend deprivation index"}
      if(i=="BMI"){n[1,1] <- "BMI, kg/m^2"}
      if(i=="prs"){n[1,1] <- "PRS for alcohol"}
      if(i=="prs_w"){n[1,1] <- "Weight PRS for alcohol"}
      
      n[1,2] <- paste0(sprintf("%0.1f",mean(tdata$var,na.rm=TRUE)), " ± ",
                       sprintf("%0.1f",sd(tdata$var,na.rm=TRUE)))
      n[1,3] <- paste0(sprintf("%0.1f",mean(tdata$var[tdata$outcome==1],na.rm=TRUE)), " ± ",
                       sprintf("%0.1f",sd(tdata$var[tdata$outcome==1],na.rm=TRUE)))
      p <- t.test(tdata$var~tdata$outcome)
      p <- as.numeric(p$p.value)
      if(p<0.001){n[1,4] <- "<0.001"}
      if(p>=0.001 & p<0.05){ n[1,4] <- sprintf("%0.3f",p)}   
      if(p>=0.05){ n[1,4] <- sprintf("%0.2f",p)}   
      
    }
    #continous - unnormal perweek_alchol_unit_cal_new
    if(i %in% c("edu_y","incomesocre",
                "permonth_alchol_drink_cal","perweek_alchol_unit_cal",
                "perweek_alchol_drink_cal_new","perweek_alchol_unit_cal_new",
                "permonth_alchol_drink_cal","perweek_alchol_unit_cal",
                "perweek_alchol_g_cal_new","perweek_alchol_g_cal"
                ,"prs","prs_w"
    )){
      n <- matrix(NA,nrow=1,ncol=4)
      if(i=="edu_y"){n[1,1] <- "Education year, year"}
      if(i=="incomesocre"){n[1,1] <- "Income score"}
      if(i=="perweek_alchol_drink_cal"){n[1,1] <- "Alcohol consumption, drink per week"}
      if(i=="perweek_alchol_unit_cal"){n[1,1] <- "Alcohol consumption, unit per week"}
      if(i=="perweek_alchol_drink_cal_new"){n[1,1] <- "Alcohol consumption, drink per week"}
      if(i=="perweek_alchol_unit_cal_new"){n[1,1] <- "Alcohol consumption, unit per week"}
      if(i=="permonth_alchol_drink_cal"){n[1,1] <- "Alcohol consumption, drink per month"}
      if(i=="permonth_alchol_unit_cal"){n[1,1] <- "Alcohol consumption, unit per month"}
      if(i=="perweek_alchol_g_cal"){n[1,1] <- "Alcohol consumption, g per week"}
      if(i=="perweek_alchol_g_cal_new"){n[1,1] <- "Alcohol consumption, g per week"}
      if(i=="prs"){n[1,1] <- "PRS for alcohol"}
      if(i=="prs_w"){n[1,1] <- "Weight PRS for alcohol"}
      s_var <- summary(tdata$var)
      n[1,2] <- paste0(sprintf("%0.2f",s_var[3])," (",sprintf("%0.2f",s_var[2]),", ",sprintf("%0.2f",s_var[5]),")")
      s_var <- summary(tdata$var[tdata$outcome==1])
      n[1,3] <- paste0(sprintf("%0.2f",s_var[3])," (",sprintf("%0.2f",s_var[2]),", ",sprintf("%0.2f",s_var[5]),")")
      
      p <- wilcox.test(tdata$var~tdata$outcome)
      p <- as.numeric(p$p.value)
      if(p<0.001){n[1,4] <- "<0.001"}
      if(p>=0.001 & p<0.05){ n[1,4] <- sprintf("%0.3f",p)}   
      if(p>=0.05){ n[1,4] <- sprintf("%0.2f",p)}   
    }
    #category for multi-categroy
    if(i %in% c("age_group","age_2group","sex","area","edu","edu_new","edu_university","townsend_group_1",
                "houseincome","smoking","smoke","BMI_group",
                "PA_group","sleep_group",
                "Diabetes","depression_g","APOE4_status",
                "alcohlol_frequency","alcohol_drinker_status",
                "family_alzheimer",
                "Diabetes","Hypertension","depression_g","APOE4_status",
                "alcohol_6group","alcohol_2group"
    )){
      
      n_dim <- dim(table(tdata$var))+1
      n_name <- names(table(tdata$var))
      
      s_table <- CrossTable(tdata$var)
      n <- matrix(NA,nrow=n_dim,ncol=5)
      n[1,1] <- i
      n[2:n_dim,1] <- n_name
      n[2:n_dim,3] <- paste0(s_table$t," (",sprintf("%0.1f",s_table$prop.row*100),")")
      
      s_table <- CrossTable(tdata$var,tdata$outcome,chisq = TRUE)
      
      n[2:n_dim,4] <- paste0(s_table$t[,2]," (",sprintf("%0.1f",s_table$prop.row[,2]*100),")")
      p <- s_table$chisq$p.value
      if(p<0.001){n[1,5] <- "<0.001"}
      if(p>=0.001 & p<0.05){ n[1,5] <- sprintf("%0.3f",p)}   
      if(p>=0.05){ n[1,5] <- sprintf("%0.2f",p)}   
      
      if(i=="age_group"){n[,2] <- c("Age, n (%)",
                                    "  <45 year","  [45, 65) year","  >= 65 year")}
      if(i=="age_2group"){n[,2] <- c("Age, n (%)",
                                     "  <=65 year", "  > 65 year")}
      if(i=="sex"){n[,2] <- c("Gender, n (%)",
                              "  Women","  Men")}
      if(i=="area"){n[,2] <- c("Area, n (%)",
                               "  England",
                               "  Scotland",
                               "  Wales")}
      if(i=="BMI_group"){n[,2] <- c("BMI , n (%)",
                                    "  <18.5  kg/m^2",
                                    "  [18.5~25.0) kg/m^2",
                                    "  [25.0~30.0) kg/m^2",
                                    "  >=30.0 kg/m^2")}
      if(i=="edu"){n[,2] <- c("Education, n (%)",
                              "  Higher",
                              "  Upper secondary",
                              "  Lower secondary",
                              "  Vocational",
                              "  Other")}
      if(i=="edu_new"){n[,2] <- c("Education, n (%)",
                                  "  Higher or vocational",
                                  "  Upper or lower secondary",
                                  "  Other")}
      if(i=="edu_university"){n[,2] <- c("College Degree, n (%)",
                                         "  Below university","  University or above")}
      if(i=="townsend_group_1"){n[,2] <- c("Socioeconomic status quintile, n (%)",
                                           "  1 (least deprived)",
                                           "  2-4",
                                           "  5 (most deprived)")}
      if(i=="houseincome"){n[,2] <- c("House income, n (%)",
                                      "  Less than 18,000",
                                      "  18,000 to 30,999",
                                      "  31,000 to 51,999",
                                      "  52,000 to 100,000",
                                      "  Greater than 100,000")}
      #smoking
      if(i=="smoking"){n[,2] <- c("Smoking, n (%)",
                                  "  No",
                                  "  Yes")}
      if(i=="smoke"){n[,2] <- c("Smoke status, n (%)",
                                "  Never",
                                "  Previous",
                                "  Current")}
      if(i=="PA_group"){n[,2] <- c("Physical Activity Level, n (%)",
                                   "  Insufficient","  Sufficient","  additional")}
      if(i=="sleep_group"){n[,2] <- c("Sleep Duration, n (%)",
                                      "  <6 hour","  6-9 hour","  >9 hour")}
      if(pre==TRUE & i=="alcohol_drinker_status"){n[,2] <- c("Aclohol drinking status, n (%)",
                                                             "  Never",
                                                             "  Previous",
                                                             "  Current")}
      if(pre==FALSE & i=="alcohol_drinker_status"){n[,2] <- c("Aclohol drinking status, n (%)",
                                                              "  Never",
                                                              "  Current")}
      if(i=="alcohlol_frequency"){n[,2] <- c("Alcohol Intake Frequency, n (%)",
                                             "  Never","  Special occasions only",
                                             "  One to three times a month","  Once or twice a week",
                                             "  Three or four times a week", "  Daily or almost daily")}
      if(i=="family_alzheimer"){n[,2] <- c("Family with Alzheimer's disease, n (%)","  Without","  With")}
      if(i=="Diabetes"){n[,2] <- c("Diabetes, n (%)","  Without","  With")}
      if(i=="Hypertension"){n[,2] <- c("Hypertension, n (%)","  Without","  With")}
      if(i=="depression_g"){n[,2] <- c("Depressive symptoms in last 2 weeks, n (%)","  Without","  With")}
      
      if(i=="APOE4_status"){n[,2] <- c("APOE e4, n (%)","  Without","  With")}
      if(i=="alcohol_2group"){n[,2] <- c("Alcohol consumption, n (%)",
                                         "  Normal (<=14 unit/week",
                                         "  Exceed (>14 unit/week)")}
      if(i=="alcohol_6group"){n[,2] <- c("Alcohol consumption, n (%)",
                                         "  Light (<=7 unit/week)",
                                         "  Moderate ((7, 14] unit/week)",
                                         "  Heavy ((14, 21] unit/week)",
                                         "  Very heavy ((21, 28] unit/week)",
                                         "  Excessive ((28, 42] unit/week)",
                                         "  Extremely excessive (>42 unit/week)")}
      
      n <- n[,c(2:5)]
    }
    
    
    N <- rbind(N,n)
  }
  N
}


######table 1######
table_1_group <- function(conf,tdata,pre=FALSE){
  tdata <- as.data.frame(tdata)
  N <- data.frame()
  N[1,1] <- "No of participant"
  N[1,2] <- dim(tdata)[1]
  N[1,3] <- paste0(dim(tdata[tdata$group==1,])[1]," (",
                   sprintf("%0.2f",dim(tdata[tdata$group==1,])[1]/dim(tdata)[1]*100), ")")
  N[1,4] <- paste0(dim(tdata[tdata$group==2,])[1]," (",
                    sprintf("%0.2f",dim(tdata[tdata$group==2,])[1]/dim(tdata)[1]*100), ")")
  N[1,5] <- NA
  for(i in conf){
    tdata$var <- NA
    tdata[,dim(tdata)[2]] <- tdata[,i]
    #continous - normal
    if(i %in% c("years","age","townsend","BMI","prs","prs_w")){
      n <- matrix(NA,nrow=1,ncol=5)
      if(i=="years"){n[1,1] <- "Follow-up years, year"}
      if(i=="age"){n[1,1] <- "Age, years"}
      if(i=="townsend"){n[1,1] <- "Townsend deprivation index"}
      if(i=="BMI"){n[1,1] <- "BMI, kg/m^2"}
      if(i=="prs"){n[1,1] <- "PRS for alcohol"}
      if(i=="prs_w"){n[1,1] <- "Weight PRS for alcohol"}
      
      n[1,2] <- paste0(sprintf("%0.1f",mean(tdata$var,na.rm=TRUE)), " ± ",
                       sprintf("%0.1f",sd(tdata$var,na.rm=TRUE)))
      n[1,3] <- paste0(sprintf("%0.1f",mean(tdata$var[tdata$group==1],na.rm=TRUE)), " ± ",
                       sprintf("%0.1f",sd(tdata$var[tdata$group==1],na.rm=TRUE)))
      n[1,4] <- paste0(sprintf("%0.1f",mean(tdata$var[tdata$group==2],na.rm=TRUE)), " ± ",
                       sprintf("%0.1f",sd(tdata$var[tdata$group==2],na.rm=TRUE)))
      
      p <- t.test(tdata$var~tdata$group)
      p <- as.numeric(p$p.value)
      if(p<0.001){n[1,5] <- "<0.001"}
      if(p>=0.001 & p<0.05){ n[1,5] <- sprintf("%0.3f",p)}   
      if(p>=0.05){ n[1,5] <- sprintf("%0.2f",p)}   
      
    }
    #continous - unnormal perweek_alchol_unit_cal_new
    if(i %in% c("edu_y","incomesocre",
                "permonth_alchol_drink_cal","perweek_alchol_unit_cal",
                "perweek_alchol_drink_cal_new","perweek_alchol_unit_cal_new",
                "permonth_alchol_drink_cal","perweek_alchol_unit_cal",
                "perweek_alchol_g_cal_new","perweek_alchol_g_cal"
                ,"prs","prs_w"
    )){
      n <- matrix(NA,nrow=1,ncol=5)
      if(i=="edu_y"){n[1,1] <- "Education year, year"}
      if(i=="incomesocre"){n[1,1] <- "Income score"}
      if(i=="perweek_alchol_drink_cal"){n[1,1] <- "Alcohol consumption, drink per week"}
      if(i=="perweek_alchol_unit_cal"){n[1,1] <- "Alcohol consumption, unit per week"}
      if(i=="perweek_alchol_drink_cal_new"){n[1,1] <- "Alcohol consumption, drink per week"}
      if(i=="perweek_alchol_unit_cal_new"){n[1,1] <- "Alcohol consumption, unit per week"}
      if(i=="permonth_alchol_drink_cal"){n[1,1] <- "Alcohol consumption, drink per month"}
      if(i=="permonth_alchol_unit_cal"){n[1,1] <- "Alcohol consumption, unit per month"}
      if(i=="perweek_alchol_g_cal"){n[1,1] <- "Alcohol consumption, g per week"}
      if(i=="perweek_alchol_g_cal_new"){n[1,1] <- "Alcohol consumption, g per week"}
      if(i=="prs"){n[1,1] <- "PRS for alcohol"}
      if(i=="prs_w"){n[1,1] <- "Weight PRS for alcohol"}
      s_var <- summary(tdata$var)
      n[1,2] <- paste0(sprintf("%0.2f",s_var[3])," (",sprintf("%0.2f",s_var[2]),", ",sprintf("%0.2f",s_var[5]),")")
      s_var <- summary(tdata$var[tdata$group==1])
      n[1,3] <- paste0(sprintf("%0.2f",s_var[3])," (",sprintf("%0.2f",s_var[2]),", ",sprintf("%0.2f",s_var[5]),")")
      s_var <- summary(tdata$var[tdata$group==2])
      n[1,4] <- paste0(sprintf("%0.2f",s_var[3])," (",sprintf("%0.2f",s_var[2]),", ",sprintf("%0.2f",s_var[5]),")")
      p <- wilcox.test(tdata$var~tdata$group)
      p <- as.numeric(p$p.value)
      if(p<0.001){n[1,5] <- "<0.001"}
      if(p>=0.001 & p<0.05){ n[1,5] <- sprintf("%0.3f",p)}   
      if(p>=0.05){ n[1,5] <- sprintf("%0.2f",p)}   
    }
    #category for multi-categroy
    if(i %in% c("age_group","age_2group","sex","area","edu","edu_new","edu_university","townsend_group_1",
                "houseincome","smoking","smoke","BMI_group",
                "PA_group","sleep_group",
                "Diabetes","depression_g","APOE4_status",
                "alcohlol_frequency","alcohol_drinker_status",
                "family_alzheimer",
                "Diabetes","Hypertension","depression_g","APOE4_status",
                "alcohol_6group","alcohol_2group","ACD"
    )){
      
      n_dim <- dim(table(tdata$var))+1
      n_name <- names(table(tdata$var))
      
      s_table <- CrossTable(tdata$var)
      n <- matrix(NA,nrow=n_dim,ncol=6)
      n[1,1] <- i
      n[2:n_dim,1] <- n_name
      n[2:n_dim,3] <- paste0(s_table$t," (",sprintf("%0.1f",s_table$prop.row*100),")")
      
      s_table <- CrossTable(tdata$var,tdata$group,chisq = TRUE)
      
      n[2:n_dim,4] <- paste0(s_table$t[,1]," (",sprintf("%0.1f",s_table$prop.col[,1]*100),")")
      n[2:n_dim,5] <- paste0(s_table$t[,2]," (",sprintf("%0.1f",s_table$prop.col[,2]*100),")")
      p <- s_table$chisq$p.value
      if(p<0.001){n[1,6] <- "<0.001"}
      if(p>=0.001 & p<0.05){ n[1,6] <- sprintf("%0.3f",p)}   
      if(p>=0.05){ n[1,6] <- sprintf("%0.2f",p)}   
      
      if(i=="age_group"){n[,2] <- c("Age, n (%)",
                                    "  <45 year","  [45, 65) year","  >= 65 year")}
      if(i=="age_2group"){n[,2] <- c("Age, n (%)",
                                     "  <=65 year", "  > 65 year")}
      if(i=="sex"){n[,2] <- c("Gender, n (%)",
                              "  Women","  Men")}
      if(i=="area"){n[,2] <- c("Area, n (%)",
                               "  England",
                               "  Scotland",
                               "  Wales")}
      if(i=="BMI_group"){n[,2] <- c("BMI , n (%)",
                                    "  <18.5  kg/m^2",
                                    "  [18.5~25.0) kg/m^2",
                                    "  [25.0~30.0) kg/m^2",
                                    "  >=30.0 kg/m^2")}
      if(i=="edu"){n[,2] <- c("Education, n (%)",
                              "  Higher",
                              "  Upper secondary",
                              "  Lower secondary",
                              "  Vocational",
                              "  Other")}
      if(i=="edu_new"){n[,2] <- c("Education, n (%)",
                                  "  Higher or vocational",
                                  "  Upper or lower secondary",
                                  "  Other")}
      if(i=="edu_university"){n[,2] <- c("College Degree, n (%)",
                                         "  Below university","  University or above")}
      if(i=="townsend_group_1"){n[,2] <- c("Socioeconomic status quintile, n (%)",
                                           "  1 (least deprived)",
                                           "  2-4",
                                           "  5 (most deprived)")}
      if(i=="houseincome"){n[,2] <- c("House income, n (%)",
                                      "  Less than 18,000",
                                      "  18,000 to 30,999",
                                      "  31,000 to 51,999",
                                      "  52,000 to 100,000",
                                      "  Greater than 100,000")}
      #smoking
      if(i=="smoking"){n[,2] <- c("Smoking, n (%)",
                                  "  No",
                                  "  Yes")}
      if(i=="smoke"){n[,2] <- c("Smoke status, n (%)",
                                "  Never",
                                "  Previous",
                                "  Current")}
      if(i=="PA_group"){n[,2] <- c("Physical Activity Level, n (%)",
                                   "  Insufficient","  Sufficient","  additional")}
      if(i=="sleep_group"){n[,2] <- c("Sleep Duration, n (%)",
                                      "  <6 hour","  6-9 hour","  >9 hour")}
      if(pre==TRUE & i=="alcohol_drinker_status"){n[,2] <- c("Aclohol drinking status, n (%)",
                                                             "  Never",
                                                             "  Previous",
                                                             "  Current")}
      if(pre==FALSE & i=="alcohol_drinker_status"){n[,2] <- c("Aclohol drinking status, n (%)",
                                                              "  Never",
                                                              "  Current")}
      if(i=="alcohlol_frequency"){n[,2] <- c("Alcohol Intake Frequency, n (%)",
                                             "  Never","  Special occasions only",
                                             "  One to three times a month","  Once or twice a week",
                                             "  Three or four times a week", "  Daily or almost daily")}
      if(i=="family_alzheimer"){n[,2] <- c("Family with Alzheimer's disease, n (%)","  Without","  With")}
      if(i=="Diabetes"){n[,2] <- c("Diabetes, n (%)","  Without","  With")}
      if(i=="Hypertension"){n[,2] <- c("Hypertension, n (%)","  Without","  With")}
      if(i=="depression_g"){n[,2] <- c("Depressive symptoms in last 2 weeks, n (%)","  Without","  With")}
      
      if(i=="APOE4_status"){n[,2] <- c("APOE e4, n (%)","  Without","  With")}
      if(i=="alcohol_2group"){n[,2] <- c("Alcohol consumption, n (%)",
                                         "  Normal (<=14 unit/week)",
                                         "  Exceed (>14 unit/week)")}
      if(i=="alcohol_6group"){n[,2] <- c("Alcohol consumption, n (%)",
                                         "  Light (<=7 unit/week)",
                                         "  Moderate ((7, 14] unit/week)",
                                         "  Heavy ((14, 21] unit/week)",
                                         "  Very heavy ((21, 28] unit/week)",
                                         "  Excessive ((28, 42] unit/week)",
                                         "  Extremely excessive (>42 unit/week)")}
      ##disease
      if(i=="ACD"){n[,2] <- c("All-cause dementia, n (%)","  Without","  With")}
      
      n <- n[,c(2:6)]
    }
    
    N <- rbind(N,n)
  }
  N
}

#####table 1: SMD####
SMD_continous_sd <- function(mean1,sd1,n1,mean2,sd2,n2,d=2){
  pooled_sd <- sqrt(((n1-1) * sd1^2 + (n2-1) * sd2^2)/(n1 + n2 - 2))
  smd <- (mean1 - mean2)/pooled_sd
  
  SE_smd <- sqrt(1/n1 + 1/n2 + smd^2/(2*(n1 +n2)))
  smd_l <- smd - 1.96 * SE_smd
  smd_u <- smd + 1.96 * SE_smd
  SMD <- sprintf("%0.*f (%0.*f, %0.*f)", d, smd, d, smd_l, d, smd_u)
  
  return(SMD)
}

SMD_continous <- function(mean1,sd1,n1,mean2,sd2,n2,d=2){
  pooled_sd <- sqrt(((n1-1) * sd1^2 + (n2-1) * sd2^2)/(n1 + n2 - 2))
  smd <- (mean1 - mean2)/pooled_sd
  abs_smd <- abs(smd)
  rounded_smd <- sprintf("%0.*f", d, abs_smd)
  
  return(rounded_smd)
}

SMD_cohens_h_sd <- function(p1,n1,p2,n2,d=2){
  h <- 2 * (asin(sqrt(p1)) - asin(sqrt(p2)))
  se_h <- sqrt(1/n1 + 1/n2)
  h_l <- h - 1.96*se_h
  h_u <- h + 1.96*se_h
  SMD <- sprintf("%0.*f (%0.*f, %0.*f)",d,h,d,h_l,d,h_u)
  return(SMD)
}

SMD_category <- function(event_count1, sample_size1, event_count2, sample_size2,d=2) {
  p1 <- event_count1 / sample_size1
  p2 <- event_count2 / sample_size2
  pooled_p <- (event_count1 + event_count2) / (sample_size1 + sample_size2)
  
  smd <- (p1 - p2) / sqrt(pooled_p * (1 - pooled_p))
  abs_smd <- abs(smd)
  rounded_smd <- sprintf("%0.*f", d, abs_smd)
  
  return(rounded_smd)
}



table_1_group_SMD <- function(conf,tdata,pre=FALSE){
  tdata <- as.data.frame(tdata)
  N <- data.frame()
  N[1,1] <- "No of participant"
  N[1,2] <- dim(tdata)[1]
  N[1,3] <- paste0(dim(tdata[tdata$group==1,])[1]," (",
                   sprintf("%0.2f",dim(tdata[tdata$group==1,])[1]/dim(tdata)[1]*100), ")")
  N[1,4] <- paste0(dim(tdata[tdata$group==2,])[1]," (",
                   sprintf("%0.2f",dim(tdata[tdata$group==2,])[1]/dim(tdata)[1]*100), ")")
  N[1,5] <- NA
  N[1,6] <- NA
  for(i in conf){
    tdata$var <- NA
    tdata[,dim(tdata)[2]] <- tdata[,i]
    #continous - normal
    if(i %in% c("years","age","townsend","BMI","prs","prs_w")){
      n <- matrix(NA,nrow=1,ncol=6)
      if(i=="years"){n[1,1] <- "Follow-up years, year"}
      if(i=="age"){n[1,1] <- "Age, years"}
      if(i=="townsend"){n[1,1] <- "Townsend deprivation index"}
      if(i=="BMI"){n[1,1] <- "BMI, kg/m^2"}
      if(i=="prs"){n[1,1] <- "PRS for alcohol"}
      if(i=="prs_w"){n[1,1] <- "Weight PRS for alcohol"}
      
      n[1,2] <- paste0(sprintf("%0.1f",mean(tdata$var,na.rm=TRUE)), " ± ",
                       sprintf("%0.1f",sd(tdata$var,na.rm=TRUE)))
      
      mean1 <- mean(tdata$var[tdata$group==1],na.rm=TRUE)
      sd1 <- sd(tdata$var[tdata$group==1],na.rm=TRUE)
      n1 <- length(tdata$var[tdata$group==1])
      mean2 <- mean(tdata$var[tdata$group==2],na.rm=TRUE)
      sd2 <- sd(tdata$var[tdata$group==2],na.rm=TRUE)
      n2 <- length(tdata$var[tdata$group==2]) 
      
      n[1,3] <- paste0(sprintf("%0.1f",mean1), " ± ",
                       sprintf("%0.1f",sd1))
      n[1,4] <- paste0(sprintf("%0.1f",mean2), " ± ",
                       sprintf("%0.1f",sd2))
      ##SMD
      SMD <- SMD_continous(mean2,sd2,n2,mean1,sd1,n1,3)
      n[1,5] <- SMD
      ##p
      p <- t.test(tdata$var~tdata$group)
      p <- as.numeric(p$p.value)
      if(p<0.001){n[1,6] <- "<0.001"}
      if(p>=0.001 & p<0.05){ n[1,6] <- sprintf("%0.3f",p)}   
      if(p>=0.05){ n[1,6] <- sprintf("%0.2f",p)}   
      
    }
    #continous - unnormal perweek_alchol_unit_cal_new
    if(i %in% c("edu_y","incomesocre",
                "permonth_alchol_drink_cal","perweek_alchol_unit_cal",
                "perweek_alchol_drink_cal_new","perweek_alchol_unit_cal_new",
                "permonth_alchol_drink_cal","perweek_alchol_unit_cal",
                "perweek_alchol_g_cal_new","perweek_alchol_g_cal"
                ,"prs","prs_w"
    )){
      n <- matrix(NA,nrow=1,ncol=6)
      if(i=="edu_y"){n[1,1] <- "Education year, year"}
      if(i=="incomesocre"){n[1,1] <- "Income score"}
      if(i=="perweek_alchol_drink_cal"){n[1,1] <- "Alcohol consumption, drink per week"}
      if(i=="perweek_alchol_unit_cal"){n[1,1] <- "Alcohol consumption, unit per week"}
      if(i=="perweek_alchol_drink_cal_new"){n[1,1] <- "Alcohol consumption, drink per week"}
      if(i=="perweek_alchol_unit_cal_new"){n[1,1] <- "Alcohol consumption, unit per week"}
      if(i=="permonth_alchol_drink_cal"){n[1,1] <- "Alcohol consumption, drink per month"}
      if(i=="permonth_alchol_unit_cal"){n[1,1] <- "Alcohol consumption, unit per month"}
      if(i=="perweek_alchol_g_cal"){n[1,1] <- "Alcohol consumption, g per week"}
      if(i=="perweek_alchol_g_cal_new"){n[1,1] <- "Alcohol consumption, g per week"}
      if(i=="prs"){n[1,1] <- "PRS for alcohol"}
      if(i=="prs_w"){n[1,1] <- "Weight PRS for alcohol"}
      s_var <- summary(tdata$var)
      n[1,2] <- paste0(sprintf("%0.2f",s_var[3])," (",sprintf("%0.2f",s_var[2]),", ",sprintf("%0.2f",s_var[5]),")")
      ##group1
      s_var <- summary(tdata$var[tdata$group==1])
      median1 <- s_var[3]
      sd1 <- s_var[5] - s_var[1]
      n1 <- length(tdata$var[tdata$group==1])
      n[1,3] <- paste0(sprintf("%0.2f",s_var[3])," (",sprintf("%0.2f",s_var[2]),", ",sprintf("%0.2f",s_var[5]),")")
      ##group2
      s_var <- summary(tdata$var[tdata$group==2])
      median2 <- s_var[3]
      sd2 <- s_var[5] - s_var[1]
      n2 <- length(tdata$var[tdata$group==1])
      n[1,4] <- paste0(sprintf("%0.2f",s_var[3])," (",sprintf("%0.2f",s_var[2]),", ",sprintf("%0.2f",s_var[5]),")")
      ##SMD
      SMD <- SMD_continous(median2,sd2,n2,median1,sd1,n1,2)
      n[1,5] <- SMD
      ##p value
      p <- wilcox.test(tdata$var~tdata$group)
      p <- as.numeric(p$p.value)
      if(is.na(p)==TRUE){n[1,6] < "NA"}
      if(is.na(p)==FALSE){
        if(p<0.001){n[1,6] <- "<0.001"}
        if(p>=0.001 & p<0.05){ n[1,6] <- sprintf("%0.3f",p)}   
        if(p>=0.05){ n[1,6] <- sprintf("%0.2f",p)} 
      }
    }
    #category for multi-categroy
    if(i %in% c("age_group","age_2group","sex","area","edu","edu_new","edu_university","townsend_group_1",
                "houseincome","smoking","smoke","BMI_group",
                "PA_group","sleep_group",
                "Diabetes","depression_g","APOE4_status",
                "alcohlol_frequency","alcohol_drinker_status",
                "family_alzheimer",
                "Diabetes","Hypertension","depression_g","APOE4_status",
                "alcohol_6group","alcohol_2group","ACD"
    )){
      
      n_dim <- dim(table(tdata$var))+1
      n_name <- names(table(tdata$var))
      
      s_table <- CrossTable(tdata$var)
      n <- matrix(NA,nrow=n_dim,ncol=7)
      n[1,1] <- i
      n[2:n_dim,1] <- n_name
      n[2:n_dim,3] <- paste0(s_table$t," (",sprintf("%0.1f",s_table$prop.row*100),")")
      
      s_table <- CrossTable(tdata$var,tdata$group,chisq = TRUE)
      
      n[2:n_dim,4] <- paste0(s_table$t[,1]," (",sprintf("%0.1f",s_table$prop.col[,1]*100),")")
      n[2:n_dim,5] <- paste0(s_table$t[,2]," (",sprintf("%0.1f",s_table$prop.col[,2]*100),")")
      #SMD
      p1 <- s_table$prop.col[,1]
      p2 <- s_table$prop.col[,2]
      SMD <- SMD_cohens_h(p2,p1,3)
      n[2:n_dim,6] <- SMD
      #3p value
      p <- s_table$chisq$p.value
      if(p<0.001){n[1,7] <- "<0.001"}
      if(p>=0.001 & p<0.05){ n[1,7] <- sprintf("%0.3f",p)}   
      if(p>=0.05){ n[1,7] <- sprintf("%0.2f",p)}   
      
      if(i=="age_group"){n[,2] <- c("Age, n (%)",
                                    "  <45 year","  [45, 65) year","  >= 65 year")}
      if(i=="age_2group"){n[,2] <- c("Age, n (%)",
                                     "  <=65 year", "  > 65 year")}
      if(i=="sex"){n[,2] <- c("Gender, n (%)",
                              "  Women","  Men")}
      if(i=="area"){n[,2] <- c("Area, n (%)",
                               "  England",
                               "  Scotland",
                               "  Wales")}
      if(i=="BMI_group"){n[,2] <- c("BMI , n (%)",
                                    "  <18.5  kg/m^2",
                                    "  [18.5~25.0) kg/m^2",
                                    "  [25.0~30.0) kg/m^2",
                                    "  >=30.0 kg/m^2")}
      if(i=="edu"){n[,2] <- c("Education, n (%)",
                              "  Higher",
                              "  Upper secondary",
                              "  Lower secondary",
                              "  Vocational",
                              "  Other")}
      if(i=="edu_new"){n[,2] <- c("Education, n (%)",
                                  "  Higher or vocational",
                                  "  Upper or lower secondary",
                                  "  Other")}
      if(i=="edu_university"){n[,2] <- c("College Degree, n (%)",
                                         "  Below university","  University or above")}
      if(i=="townsend_group_1"){n[,2] <- c("Socioeconomic status quintile, n (%)",
                                           "  1 (least deprived)",
                                           "  2-4",
                                           "  5 (most deprived)")}
      if(i=="houseincome"){n[,2] <- c("House income, n (%)",
                                      "  Less than 18,000",
                                      "  18,000 to 30,999",
                                      "  31,000 to 51,999",
                                      "  52,000 to 100,000",
                                      "  Greater than 100,000")}
      #smoking
      if(i=="smoking"){n[,2] <- c("Smoking, n (%)",
                                  "  No",
                                  "  Yes")}
      if(i=="smoke"){n[,2] <- c("Smoke status, n (%)",
                                "  Never",
                                "  Previous",
                                "  Current")}
      if(i=="PA_group"){n[,2] <- c("Physical Activity Level, n (%)",
                                   "  Insufficient","  Sufficient","  additional")}
      if(i=="sleep_group"){n[,2] <- c("Sleep Duration, n (%)",
                                      "  <6 hour","  6-9 hour","  >9 hour")}
      if(pre==TRUE & i=="alcohol_drinker_status"){n[,2] <- c("Aclohol drinking status, n (%)",
                                                             "  Never",
                                                             "  Previous",
                                                             "  Current")}
      if(pre==FALSE & i=="alcohol_drinker_status"){n[,2] <- c("Aclohol drinking status, n (%)",
                                                              "  Never",
                                                              "  Current")}
      if(i=="alcohlol_frequency"){n[,2] <- c("Alcohol Intake Frequency, n (%)",
                                             "  Never","  Special occasions only",
                                             "  One to three times a month","  Once or twice a week",
                                             "  Three or four times a week", "  Daily or almost daily")}
      if(i=="family_alzheimer"){n[,2] <- c("Family with Alzheimer's disease, n (%)","  Without","  With")}
      if(i=="Diabetes"){n[,2] <- c("Diabetes, n (%)","  Without","  With")}
      if(i=="Hypertension"){n[,2] <- c("Hypertension, n (%)","  Without","  With")}
      if(i=="depression_g"){n[,2] <- c("Depressive symptoms in last 2 weeks, n (%)","  Without","  With")}
      
      if(i=="APOE4_status"){n[,2] <- c("APOE e4, n (%)","  Without","  With")}
      if(i=="alcohol_2group"){n[,2] <- c("Alcohol consumption, n (%)",
                                         "  Normal (<=14 unit/week)",
                                         "  Exceed (>14 unit/week)")}
      if(i=="alcohol_6group"){n[,2] <- c("Alcohol consumption, n (%)",
                                         "  Light (<=7 unit/week)",
                                         "  Moderate ((7, 14] unit/week)",
                                         "  Heavy ((14, 21] unit/week)",
                                         "  Very heavy ((21, 28] unit/week)",
                                         "  Excessive ((28, 42] unit/week)",
                                         "  Extremely excessive (>42 unit/week)")}
      ##disease
      if(i=="ACD"){n[,2] <- c("All-cause dementia, n (%)","  Without","  With")}
      
      n <- n[,c(2:7)]
    }
    
    N <- rbind(N,n)
  }
  N
}


######table 2
table_2_incidence <- function(modeldata){
  inc_d <- modeldata %>% group_by(E_group) %>% 
    mutate(all=n(),
           all_t=sum(years)) %>% 
    group_by(E_group,outcome) %>%
    mutate(n1=n(),
           t1=sum(years)) %>%
    select(E_group,outcome,all,all_t,n1,t1)
  inc_d <- unique(inc_d)
  inc_d_case <- inc_d %>% filter(outcome==1)
  inc_d_case$c1 <- inc_d_case$n1/inc_d_case$t1*1000
  inc_d_case$c2 <- inc_d_case$all/inc_d_case$all_t*1000
  inc_d_case$RR <- inc_d_case$c1/inc_d_case$c2
  inc_d_case$d <- log(inc_d_case$RR)
  inc_d_case$lower <- exp(inc_d_case$d - 1.96*sqrt(1/inc_d_case$n1+1/inc_d_case$all))
  inc_d_case$upper <- exp(inc_d_case$d + 1.96*sqrt(1/inc_d_case$n1+1/inc_d_case$all))
  inc_d_case$incident <- (inc_d_case$n1*inc_d_case$t1)/(inc_d_case$all*inc_d_case$all_t)*1000
  inc_d_case$case <- paste0(inc_d_case$n1,"/",inc_d_case$all)
  inc_d_case$pyear <- paste0(sprintf("%0.2f",inc_d_case$RR),"(",
                             sprintf("%0.2f",inc_d_case$lower),", ",
                             sprintf("%0.2f",inc_d_case$upper),")")
  inc_d_case$group <- inc_d_case$E_group
  inc_d_case <- inc_d_case[,c("group","case","pyear")]
  inc_d_case <- inc_d_case[order(inc_d_case$group),]
  inc_d_case <- t(inc_d_case)
  inc_d_case
}


table_2_age_adjusted_incidence <- function(modeldata, standard_pop) {
  # 计算每个年龄组内的发病率
  incidence_data <- modeldata %>% 
    group_by(E_group) %>%
    summarize(total_cases = sum(outcome), 
              total_person_years = sum(years), 
              crude_incidence = total_cases / total_person_years * 1000)
  
  # 加入标准人口比例
  incidence_data <- merge(incidence_data, standard_pop, by = "E_group")
  
  # 计算加权发病率
  incidence_data$weighted_incidence <- incidence_data$crude_incidence * incidence_data$pop_proportion
  
  # 计算年龄调整发病率
  age_adjusted_incidence <- sum(incidence_data$weighted_incidence)
  
  return(age_adjusted_incidence)
}

########combine two plot function########
library(gtable)
library(grid)
library(jtools)
library(ggplot2)
ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  G <- g1
  G
}






######rcs figure
library(survival)
library(ggsci)
library(nnet)
library(rms)
library(gridExtra)
rcs_fig <- function(modeldata,covar_f_name,population,plotx,d=6){
  coxplot_data <- modeldata %>% select("years","outcome","Exposure",covar_f_name)
  coxplot_data <- na.omit(coxplot_data)
  coxplot_data$x <- coxplot_data$Exposure
  #Cox model
  c_d1 <- class.ind(coxplot_data$edu_new)
  c_d1 <- c_d1[,-dim(c_d1)[2]]
  
  c_d2 <- class.ind(coxplot_data$area)
  c_d2 <- c_d2[,-dim(c_d2)[2]]
  
  c_d3 <- class.ind(coxplot_data$townsend_group_1)
  c_d3 <- c_d3[,-dim(c_d3)[2]]
  
  ns <- d
  coxplot_data <- cbind(coxplot_data[,1:ns],c_d1, c_d2, c_d3) 
  names(coxplot_data)[(ns+1):(ns+2)] <- paste0("edu",1:2)
  names(coxplot_data)[(ns+5):(ns+6)] <- paste0("townsend",1:2)
  
  fun <- "Surv(years,outcome) ~ rcs(Exposure,4)"
  fun1 <- "Surv(years,outcome) ~ Exposure"
  for(o in c(4:length(names(coxplot_data)))){ #29
    fun <- paste(fun,"+",names(coxplot_data)[o])
    fun1 <- paste(fun1,"+",names(coxplot_data)[o])
  }
  dd <- datadist(coxplot_data)
  options(datadist='dd')
  m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=coxplot_data)
  m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=coxplot_data)
  dd$limits$Exposure[2] <- 1
  m0 <- update(m0)
  p_test <- anova(m0)
  p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
  p_allover <- p_test[1,3]
  
  p_allover <- p_out(p_allover)
  p_nonliner <- p_out_fig(p_nonliner)
  
  HR <- Predict(m0,Exposure, fun=exp,ref.zero=TRUE)
  #plot2
  label_text <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                           env = base::list(p_nonliner1=p_nonliner))
  g0 <- ggplot() +
    geom_line(data=HR, aes(Exposure,yhat),
              linetype="solid",size=0.5,alpha = 1,colour="skyblue4") +
    geom_ribbon(data=HR, aes(Exposure,ymin = lower, ymax = upper),
                alpha = 0.6,fill="skyblue1") +
    theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                     axis.text = element_text(color="black", size=12),
                                     axis.title=element_text(size=12),
                                     legend.position="top") +
    geom_hline(yintercept=1, linetype="dashed", color = "grey50")
  ######figure#####
  G <- g0 + labs(subtitle = label_text)
  G
}

#####output of p value####
p_out <- function(p){
  if(p>=0.01){p_value=sprintf("%0.2f",p)}
  if(p>=0.001 & p<0.01){p_value=sprintf("%0.3f",p)}
  if(p<0.001){p_value="<0.001"}
  p_value
}


p_out_fig <- function(p){
  if(p < 0.001){p <- "< 0.001"}
  if(p >= 0.001 & p < 0.01){p <- paste("=",sprintf("%0.3f",p))}
  if(p >= 0.01){p <- paste("=",sprintf("%0.2f",p))}
  p
}


#### non-linear MR function###
source(paste0(project_path,"code/SUMnlmr/haodong_functions.R"))
source(paste0(project_path,"code/SUMnlmr/mr_summarise.R"))
source(paste0(project_path,"code/SUMnlmr/piecewise_new.R"))
source(paste0(project_path,"code/SUMnlmr/nlme_summ_aes.R"))
source(paste0(project_path,"code/SUMnlmr/nlme_summ_aes_rest.R"))

library(metafor)
####non-linear MR###
sumNLMR_rank_HR <- function(DATA_m,q=50){
  mdata <- DATA_m[,c(y_time,ex_name,x_prs_name,y_name,covar_f_name)]
  #missing 972 smoke
  mdata <- na.omit(mdata)
  num <- dim(mdata)[1]
  ym <- Surv(mdata[,y_time],mdata[,y_name])
  xm <- mdata[,ex_name]
  x_prsm <- mdata[,x_prs_name]
  covar_fm <- mdata[,covar_f_name]
  summ_data <- create_nlmr_summary(y = ym,
                                   x = xm,
                                   g = x_prsm,
                                   covar = covar_fm,
                                   family = "coxph",
                                   #controlsonly=TRUE,
                                   strata_method = "ranked", 
                                   report_GR=TRUE,
                                   report_het=TRUE,
                                   extra_statistics=TRUE,
                                   q = q)
  summ_data$summary
  model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                     fig=TRUE,family="coxph",ref=min(xmean)))
  nlmr <- list(summ_data,model)
}

sumNLMR_figure_rank_HR <- function(DATA_m,puplation,plotx,ploty,mn=1,q=50){
  mdata <- DATA_m[,c(y_time,ex_name,x_prs_name,y_name,covar_f_name)]
  #missing 972 smoke
  mdata <- na.omit(mdata)
  num <- dim(mdata)[1]
  ym <- Surv(mdata[,y_time],mdata[,y_name])
  xm <- mdata[,ex_name]
  x_prsm <- mdata[,x_prs_name]
  covar_fm <- mdata[,covar_f_name]
  summ_data <- create_nlmr_summary(y = ym,
                                   x = xm,
                                   g = x_prsm,
                                   covar = covar_fm,
                                   family = "coxph",
                                   #controlsonly=TRUE,
                                   #strata_method = "ranked", 
                                   strata_method="residual",
                                   q = q)
  summ_data$summary
  model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                     fig=TRUE,family="coxph",ref=min(xmean)))
  
  
  summary(model)
  p_allover <- model$coefficients[5]
  p_quad <- model$p_tests[2]
  #p value
  p_quad0 <- p_out_fig(p_quad)
  p_allover0 <- p_out_fig(p_allover)
  
  f_data <- model$figure$data
  
  label_text <- substitute(expr=paste(puplation1," (n=",num1,"): Overall ",italic("P"),p_allover1, 
                                      ", Nonlinear ",italic("P"), p_quad1),
                           env = base::list(p_allover1=p_allover0,
                                            p_quad1=p_quad0,
                                            num1=num,
                                            puplation1=puplation))
  label_text <- substitute(expr=paste("DRnlMR test : Overall ",italic("P"),p_allover1, 
                                      ", Nonlinear ",italic("P"), p_quad1),
                           env = base::list(p_allover1=p_allover0,
                                            p_quad1=p_quad0,
                                            num1=num,
                                            puplation1=puplation,
                                            mn1=mn))
  plot_nlmr <- model$figure +
    scale_y_continuous(name=ploty,n.breaks = 10) +
    scale_x_continuous(name=plotx,expand = c(0,0), n.breaks = 10)   +
    theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                     axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
    labs(subtitle = label_text)
  plot_nlmr
}

sumNLMR_figure_res_HR <- function(DATA_m,puplation,plotx,ploty,mn=1,q=100){
  mdata <- DATA_m[,c(y_time,ex_name,x_prs_name,y_name,covar_f_name)]
  #missing 972 smoke
  mdata <- na.omit(mdata)
  num <- dim(mdata)[1]
  ym <- Surv(mdata[,y_time],mdata[,y_name])
  xm <- mdata[,ex_name]
  x_prsm <- mdata[,x_prs_name]
  covar_fm <- mdata[,covar_f_name]
  summ_data <- create_nlmr_summary(y = ym,
                                   x = xm,
                                   g = x_prsm,
                                   covar = covar_fm,
                                   family = "coxph",
                                   #controlsonly=TRUE,
                                   strata_method = "residual", 
                                   q = q)
  model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                     fig=TRUE,family="coxph",ref=min(xmean)))
  p_allover <- model$coefficients[5]
  p_quad <- model$p_tests[3]
  #p value
  p_quad0 <- p_out_fig(p_quad)
  p_allover0 <- p_out_fig(p_allover)
  
  f_data <- model$figure$data

  label_text <- substitute(expr=paste("nlMR test: Overall ",italic("P"),p_allover1, 
                                      ", Nonlinear ",italic("P"), p_quad1),
                           env = base::list(p_allover1=p_allover0,
                                            p_quad1=p_quad0,
                                            num1=num,
                                            puplation1=puplation,
                                            mn1=mn))
  plot_nlmr <- model$figure +
    scale_y_continuous(name=ploty,n.breaks = 10) +
    scale_x_continuous(name=plotx,expand = c(0,0), n.breaks = 10)   +
    theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                     axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
    labs(subtitle = label_text)
  plot_nlmr
}

sumNLMR_res_HR <- function(DATA_m,q=100){
  mdata <- DATA_m[,c(y_time,ex_name,x_prs_name,y_name,covar_f_name)]
  #missing 972 smoke
  mdata <- na.omit(mdata)
  num <- dim(mdata)[1]
  ym <- Surv(mdata[,y_time],mdata[,y_name])
  xm <- mdata[,ex_name]
  x_prsm <- mdata[,x_prs_name]
  covar_fm <- mdata[,covar_f_name]
  summ_data <- create_nlmr_summary(y = ym,
                                   x = xm,
                                   g = x_prsm,
                                   covar = covar_fm,
                                   family = "coxph",
                                   #controlsonly=TRUE,
                                   strata_method = "residual", 
                                   q = q)
  model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                     fig=TRUE,family="coxph",ref=min(xmean)))
  nlmr <- list(summ_data,model)
}

################linear MR########################
####MR HR function##
R2_F_output <- function(modeldata,x="IV"){
  ##IV-exposure association
  fun <- "Exposure ~ IV "
  for(o in covar_f_name){fun <- paste(fun,"+",o)}
  model1 <- lm(as.formula(fun), data=modeldata)
  # Extract the first stage R-squared
  model1_r2 <- summary(model1)$r.squared
  
  ##adjust R2
  fun0 <- paste("Exposure ~", paste(covar_f_name,collapse=" + "))
  model0 <- lm(as.formula(fun0), data=modeldata)
  model0_r2 <- summary(model0)$r.squared
  
  IV_R2 = (model1_r2 - model0_r2)/(1 - model0_r2)
  
  RSS1 <- deviance(model1)
  df1 <- df.residual(model1)
  
  RSS0 <- deviance(model0)
  df0 <- df.residual(model0)
  
  IV_F <- (RSS1-RSS0)/(df1-df0)/(RSS0/df0)
  
  # Calculate the number of observations (n) and the number of coefficients (k) in the first stage
  n <- nobs(model1)
  k <- length(coef(model1))
  # Calculate the first-stage F-statistic
  F_value <- (model1_r2 / (1 - model1_r2)) * ((n - k - 1) / k)
  R2_value = model1_r2
  
  ##adjusted R-squred
  adjust_r2 <- 1 - (1-model1_r2)*(n-1)/(n-k-2)
  adjust_f <- (adjust_r2 / (1 - adjust_r2)) * ((n - k - 1) / k)
  R2_value_adjsut =  sprintf("%0.2f",adjust_r2*100)
  F_value_adjsut = sprintf("%0.1f",adjust_f)
  
  R2_F <- data.frame(x,F_value_adjsut,R2_value_adjsut)
  R2_F
}


IV_HR_function <- function(modeldata,covar_f_name,x){
  ##IV-exposure association
  fun <- "Exposure ~ IV "
  for(o in covar_f_name){fun <- paste(fun,"+",o)}
  model1 <- lm(as.formula(fun), data=modeldata)
  # Extract the first stage R-squared
  model1_r2 <- summary(model1)$r.squared
  
  ##adjust R2
  fun0 <- paste("Exposure ~", paste(covar_f_name,collapse=" + "))
  model0 <- lm(as.formula(fun0), data=modeldata)
  model0_r2 <- summary(model0)$r.squared
  
  IV_R2 = (model1_r2 - model0_r2)/(1 - model0_r2)
  
  RSS1 <- deviance(model1)
  df1 <- df.residual(model1)
  
  RSS0 <- deviance(model0)
  df0 <- df.residual(model0)
  
  IV_F <- (RSS1-RSS0)/(df1-df0)/(RSS0/df0)
  
  # Calculate the number of observations (n) and the number of coefficients (k) in the first stage
  n <- nobs(model1)
  k <- length(coef(model1))
  # Calculate the first-stage F-statistic
  F_value <- (model1_r2 / (1 - model1_r2)) * ((n - k - 1) / k)
  R2_value = model1_r2
  
  
  ##adjusted R-squred
  adjust_r2 <- 1 - (1-model1_r2)*(n-1)/(n-k-2)
  adjust_f <- (adjust_r2 / (1 - adjust_r2)) * ((n - k - 1) / k)
  R2_value_adjsut = adjust_r2
  F_value_adjsut = adjust_f
  
  beta_PRS <- summary(model1)$coefficients[2,1]
  modeldata$predict <- predict(model1)
  #modeldata$predict <- modeldata$IV*beta_PRS
  ##MR - HR
  fun_out <- "Surv(years, outcome) ~ predict "
  for(o in covar_f_name){fun_out <- paste(fun_out,"+",o)}
  fit_out <- coxph(as.formula(fun_out), data=modeldata)
  mr_result <- summary(fit_out)$coefficients[1,]
  hr <- exp(mr_result[1])
  HR_low <- exp(mr_result[1]-1.96*mr_result[3])
  HR_up <- exp(mr_result[1]+1.96*mr_result[3])
  p_value <- p_out(mr_result[5])
  HR <- paste0(sprintf("%0.2f",hr), " (",sprintf("%0.2f",HR_low),"-",sprintf("%0.2f",HR_up),")")
  #MR - OR
  ##MR - HR
  fun_out <- "outcome ~ predict "
  for(o in covar_f_name){fun_out <- paste(fun_out,"+",o)}
  fit_out <- glm(as.formula(fun_out), data=modeldata,family=binomial("logit"))
  mr_result_or <- summary(fit_out)$coefficients[2,]
  or <- exp(mr_result_or[1])
  or_low <- exp(mr_result_or[1]-1.96*mr_result_or[2])
  or_up <- exp(mr_result_or[1]+1.96*mr_result_or[2])
  p_value_OR <- p_out(mr_result_or[4])
  OR <- paste0(sprintf("%0.2f",or), " (",sprintf("%0.2f",or_low),"-",sprintf("%0.2f",or_up),")")
  IV_result <- data.frame(x,IV_F,IV_R2,F_value,R2_value,F_value_adjsut,R2_value_adjsut,hr,HR_low,HR_up,HR,p_value,OR,p_value_OR)
  IV_result
}

HR_output <- function(fit_out){
  mr_result <- summary(fit_out)$coefficients[1,]
  hr <- exp(mr_result[1])
  HR_low <- exp(mr_result[1]-1.96*mr_result[3])
  HR_up <- exp(mr_result[1]+1.96*mr_result[3])
  p_value <- p_out(mr_result[5])
  HR <- paste0(sprintf("%0.2f",hr), " (",sprintf("%0.2f",HR_low),"-",sprintf("%0.2f",HR_up),")")
  HR_result <- data.frame(hr,HR_low,HR_up,HR,p_value)
  
  HR_result
}

GS_mr_one_HR <- function(DATA_m,covar_f_name,x){
  re <- IV_HR_function(DATA_m,covar_f_name,x)
  re
}

PRS_mr_one_HR <- function(DATA_m,IV_name,population="All",subgroup="None",covar_f_name,a=2){
  ##case
  rate <- sprintf("%0.1f",(dim(DATA_m[DATA_m$outcome==1,])[1]/dim(DATA_m)[1])*100)
  
  case_total <- paste0(dim(DATA_m[DATA_m$outcome==1,])[1],"/",dim(DATA_m)[1], " (",rate,")")
  if(a==1){
    #Exposure 1
    x = "Expousre 1: drink alcohol per week"
    DATA_m$Exposure <- DATA_m$perweek_alchol_drink_cal_new
    re <- IV_HR_function(DATA_m,covar_f_name,x)
  }
  if(a==2){
    #Exposure 2
    x = "Expousre 2: unit alcohol per week"
    DATA_m$Exposure <- DATA_m$perweek_alchol_unit_cal_new
    re <- IV_HR_function(DATA_m,covar_f_name,x)
  }
  if(a==3){
    #Exposure 3
    x = "Expousre 3: standar drink alcohol per week"
    DATA_m$Exposure <- DATA_m$perweek_alchol_standar_drink_cal_new
    re <- IV_HR_function(DATA_m,covar_f_name,x)
  }
  if(a==4){
    #Exposure 4
    x = "Expousre 4: log10(unit alcohol per week)"
    DATA_m$Exposure <- log10(DATA_m$perweek_alchol_unit_cal_new)
    re <- IV_HR_function(DATA_m,covar_f_name,x)
  }
  Result_PRS_one <- cbind(population,subgroup,IV_name,case_total,re)
  Result_PRS_one
}

library(cmprsk)
#####compete risk analysis
IV_HR_function_compete <- function(modeldata,covar_f_name,x){
  fun <- "Exposure ~ IV "
  for(o in covar_f_name){ 
    fun <- paste(fun,"+",o)
  }
  fit <- lm(as.formula(fun), data=modeldata)
  summary(fit)
  r <- summary(fit)
  res <- r$coefficients[2,1]
  modeldata$predict <- predict(fit)
  
  ftime <- modeldata[,"years"]
  cov <- modeldata[,c("predict",covar_f_name)]
  fstatus <- modeldata[,"outcome_new"]
  f <- crr(ftime, fstatus, cov, failcode=1,cencode=0)
  a <- summary(f)
  a$coef
  a$conf.int
  mr_result <- summary(f)$conf.int[1,]
  hr <- mr_result[1]
  HR_low <- mr_result[3]
  HR_up <- mr_result[4]
  p_value <- p_out(a$coef[1,5])
  HR <- paste0(sprintf("%0.2f",mr_result[1]), " (",sprintf("%0.2f",mr_result[3]),"-",sprintf("%0.2f",mr_result[4]),")")
  
  IV_re <- data.frame(x,hr,HR_low,HR_up,HR,p_value)
  IV_re
}

#####linear MR residual group#####
MR_analysis_residual_HR <- function(DATA_m,IV_name,covar_f_name,population,a,g=5){
  if(a==1){DATA_m$Exposure <- DATA_m$perweek_alchol_drink_cal_new}
  if(a==2){DATA_m$Exposure <- DATA_m$perweek_alchol_unit_cal_new}
  if(a==3){DATA_m$Exposure <- DATA_m$perweek_alchol_standar_drink_cal_new}
  fun <- "Exposure ~ IV "
  for(o in covar_f_name){fun <- paste(fun,"+",o)}
  fit <- lm(as.formula(fun), data=DATA_m)
  r <- summary(fit)
  res <- r$coefficients[2,1]
  DATA_m$predict <- predict(fit)
  DATA_m$res <- DATA_m$Exposure - DATA_m$predict
  if(g==5){
    DATA_m$res_group <- cut(DATA_m$res, 
                            breaks=c(min(DATA_m$res)-0.01,7,14,21,28,max(DATA_m$res)))
    subgroup_names <- c("<=7","(7,14]","(14,21]","(21,28]",">28")
  }
  if(g==2){
    DATA_m$res_group <- cut(DATA_m$res, include.lowest=TRUE,
                            breaks=c(min(DATA_m$res)-0.01,14,max(DATA_m$res)))
    subgroup_names <- c("<=14",">14")
  }
  DATA_m$res_group_g <- as.numeric(DATA_m$res_group)
  ##p trend test
  fun_out <- "Surv(years, outcome) ~ res_group_g "
  for(o in covar_f_name){ 
    fun_out <- paste(fun_out,"+",o)
  }
  fit_out <- coxph(as.formula(fun_out), data=DATA_m)
  p_trend <- summary(fit_out)$coefficients[1,5]
  p_trend <- p_out(p_trend)
  Result_mr <- data.frame()
  for(gx in 1:g){
    DATA_m_sub <- DATA_m %>% filter(DATA_m$res_group_g==gx)
    subgroup <- subgroup_names[gx]
    res_g <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
    ##compete risk event
    modeldata <- DATA_m_sub
    modeldata$death <- ifelse(is.na(modeldata$death_date)==FALSE,1,0)
    modeldata$outcome_new <- modeldata$outcome
    modeldata$outcome_new[modeldata$death==1] <- 2
    modeldata$area <- ifelse(modeldata$area=="England",1,ifelse(modeldata$area=="Scotland",2,3))
    covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
    x="weight 95 PRS"
    res_compet_w95_sub <- IV_HR_function_compete(modeldata,covar_f_name,x)
    res_g$comp_HR <- res_compet_w95_sub$HR
    res_g$comp_p <- res_compet_w95_sub$p_value
    Result_mr <- rbind(Result_mr,res_g)
  }
  
  p_trend_result <- c("p for trend",rep(NA,dim(Result_mr)[2]-2),p_trend)
  Result_mr <- rbind(Result_mr,p_trend_result)
  Result_mr
}

###subgroup#####
MR_analysis_subgroup_HR <- function(DATA_m,IV_name,covar_f_name,a){
  Result_mr <- data.frame()
  population <- "Total"
  subgroup <- "None"
  re <- PRS_mr_one_HR(DATA_m,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  ##age group
  population <- "Age"
  DATA_m_sub <- DATA_m %>% filter(age_group==0)
  subgroup <- "<45 year"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(age_group==1)
  subgroup <- "[45, 65) year"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(age_group==2)
  subgroup <- ">= 65 year"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  ##subgroup sex
  population <- "Sex"
  DATA_m_sub <- DATA_m %>% filter(sex==1)
  subgroup <- "Men"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(sex==0)
  subgroup <- "Women"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  ##Socioeconomic status group
  population <- "Socioeconomic status group"
  DATA_m_sub <- DATA_m %>% filter(townsend_group_1==1)
  subgroup <- "least deprived (Quintile 1)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(townsend_group_1==2)
  subgroup <- "middle (Quintile 2-4)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(townsend_group_1==3)
  subgroup <- "most deprived (Quintile 5)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  ##subgroup education
  population <- "Education level"
  DATA_m_sub <- DATA_m %>% filter(edu_new==1)
  subgroup <- "higher education/vocational "
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(edu_new==2)
  subgroup <- "Secondary education "
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(edu_new==3)
  subgroup <- "other"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  ##subgroup BMI
  population <- "BMI"
  DATA_m_sub <- DATA_m %>% filter(BMI_group==0)
  subgroup <- "light (<18.5 kg/m2)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(BMI_group==1)
  subgroup <- "normal (18.5-24.9 kg/m2)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(BMI_group==2)
  subgroup <- "heavy (25-29.9 kg/m2)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(BMI_group==3)
  subgroup <- "overweight (>30 kg/m2)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  ##subgroup smoke
  population <- "Smoke status"
  DATA_m_sub <- DATA_m %>% filter(smoke==0)
  subgroup <- "Non-smokers"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(smoke==1)
  subgroup <- "Ex-smokers"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(smoke==2)
  subgroup <- "Smokers"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  ##sleep group
  population <- "Sleep duration"
  DATA_m_sub <- DATA_m %>% filter(sleep_group==1)
  subgroup <- "<6 hour"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(sleep_group==2)
  subgroup <- "6-9 hour"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(sleep_group==3)
  subgroup <- "> 9 hour"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  ##subgroup physical activity group
  population <- "Physical activity group"
  DATA_m_sub <- DATA_m %>% filter(PA_group==1)
  subgroup <- "Insufficient (Tertile 1)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(PA_group==2)
  subgroup <- "sufficient (Tertile 2)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(PA_group==3)
  subgroup <- "additional (Tertile 3)"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  ##subgroup APOE
  population <- "APOE gene"
  DATA_m_sub <- DATA_m %>% filter(APOE4_status==1)
  subgroup <- "APOE e4"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  DATA_m_sub <- DATA_m %>% filter(APOE4_status==0)
  subgroup <- "Non-APOE e4"
  re <- PRS_mr_one_HR(DATA_m_sub,IV_name,population,subgroup,covar_f_name,a)
  Result_mr <- rbind(Result_mr,re)
  Result_mr
}


####2 fold summary MR###
twofold_sMR <- function(df,snp_name){
  df1 <- sample_n(df,dim(df)[1]/2)
  df2 <- df %>% filter(!(eid %in% df1$eid))
  ##with outcome
  asso_disease_snp_df1_s <- data.frame()
  for(s in snp_name){
    modeldata <- df1
    modeldata$var <- modeldata[,s]
    covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
    
    fun_out <- "Surv(years, outcome) ~ var "
    for(o in covar_f_name){ 
      fun_out <- paste(fun_out,"+",o)
    }
    fit_out <- coxph(as.formula(fun_out), data=modeldata)
    
    r <- summary(fit_out)
    res <- r$coefficients[1,]
    res0 <- c(s,res)
    asso_disease_snp_df1_s <- rbind(asso_disease_snp_df1_s,res0)
    
  }
  names(asso_disease_snp_df1_s) <- c("SNP","Beta","HR","SE","z","Pvalue")
  head(asso_disease_snp_df1_s)
  ##with exposure
  asso_exposure_snp_df2_s <- data.frame()
  for(s in snp_name){
    modeldata <- df2
    modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
    modeldata$var <- modeldata[,s]
    covar_f_name <- c("age","sex","area","smoke","array",paste0("component",1:20))
    
    fun <- "Exposure ~ var "
    for(o in covar_f_name){ fun <- paste(fun,"+",o)}
    fit <- glm(as.formula(fun), data=modeldata)
    r <- summary(fit)
    res <- r$coefficients[2,]
    res0 <- c(s,res)
    asso_exposure_snp_df2_s <- rbind(asso_exposure_snp_df2_s,res0)
  }
  names(asso_exposure_snp_df2_s) <- c("SNP","Beta","SE","z","Pvalue")
  head(asso_exposure_snp_df2_s)
  
  ##harmanazed the data
  dy <- asso_disease_snp_df1_s %>% select(SNP,Beta,SE)
  names(dy)[2:3] <- c("beta.outcome","se.outcome")
  dy$beta.outcome <- as.numeric(dy$beta.outcome)
  dy$se.outcome <- as.numeric(dy$se.outcome)
  
  dy$outcome <- "dementia"
  dy$id.outcome <- dy$outcome
  
  dx <- asso_exposure_snp_df2_s %>% select(SNP,Beta,SE) 
  names(dx)[2:3] <- c("beta.exposure","se.exposure")
  dx$beta.exposure <- as.numeric(dx$beta.exposure)
  dx$se.exposure <- as.numeric(dx$se.exposure)
  dx$exposure <- "drink_all"
  dx$id.exposure <- dx$exposure
  
  d_mr <- merge(dx,dy,by="SNP")
  d_mr$mr_keep <- TRUE
  ######two sample MR analysis --95 SNP######
  d <- d_mr
  b_out <- d$beta.outcome
  b_exp <- d$beta.exposure
  se_out <- d$se.outcome
  se_exp <- d$se.exposure
  snp <- d$SNP
  ##IVW
  res_ivw <- mr_ivw(b_exp, b_out, se_exp, se_out)
  ##MR-Egger
  res <- mr_egger_regression(b_exp, b_out, se_exp, se_out)
  res_egger_test <-  data.frame(b=res$b_i, se=res$se_i, pval=res$pval_i)
  res_egger <- data.frame(b=res$b, se=res$se, pval=res$pval, nsnp=res$nsnp, Q = res$Q, Q_df = res$Q_df, Q_pval = res$Q_pval)
  ##weighted 
  res_weimedian <- mr_weighted_median(b_exp, b_out, se_exp, se_out)
  method <- c("IVW","MR-Egger","Weighted median")
  res_mr <- rbind(res_ivw,res_egger,res_weimedian)
  res_mr <- cbind(method,res_mr)
  res_mr$hr <- exp(res_mr$b)
  res_mr$lp <- exp(res_mr$b-1.96*res_mr$se)
  res_mr$up <- exp(res_mr$b+1.96*res_mr$se)
  res_mr$HR <- paste0(sprintf("%0.2f",res_mr$hr), " (", sprintf("%0.2f",res_mr$lp),", ", sprintf("%0.2f",res_mr$up),")")
  result_list <- list(res_mr,res_egger_test)
  result_list
}


####positive control: alcohol liver disease######
nmlr_boxfigure_single <- function(data_strata){
  fdata_sum <- data_strata$figure$data
  q10_cut <- quantile(fdata_sum$x,probs=seq(0,1,0.1))
  fdata_sum$x_g <- cut(fdata_sum$x,breaks=q10_cut,include.lowest = T)
  table(fdata_sum$x_g)
  fdata_sum$hr <- exp(fdata_sum$yest)
  p_quad0 <- data_strata$p_tests[1]
  label_text <- substitute(expr=paste("Nonlinear ",italic("P"), p_quad1),
                           env = base::list(p_quad1=p_quad0))
  
  nlmr_figure <- ggplot(fdata_sum,aes(x=x_g,y=hr))+
    geom_boxplot() + 
    stat_summary(fun.y = mean, geom = "point", shape = 23, size=1) +
    
    geom_hline(yintercept =1.0,linetype='dashed',color="grey")+
    theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                     axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                     legend.position="bottom")+
    scale_y_continuous(name="pHR of Alcohol liver diseases",n.breaks = 10) +
    scale_x_discrete(name="predicted Alcohol consumption (unit/week)") +
    scale_color_lancet() +
    labs(subtitle = label_text)
  nlmr_figure 
}


nmlr_boxfigure <- function(data_strata,n,tsize=12){
  fdata_sum <- data_strata$figure$data
  q_cut <- quantile(fdata_sum$x,probs=seq(0,1,1/n))
  fdata_sum$x_g <- cut(fdata_sum$x,breaks=q_cut,include.lowest = T)
  table(fdata_sum$x_g)
  fdata_sum$x_g <- as.numeric(fdata_sum$x_g)
  for(xi in 1:n){
    fdata_sum$predx[fdata_sum$x_g==xi] <- round(mean(fdata_sum$x[fdata_sum$x_g==xi],na.rm=T),1)
  }
  fdata_sum$hr <- exp(fdata_sum$yest)
  fdata_sum$predx <- as.factor(fdata_sum$predx)
  p_quad0 <- p_out_fig(data_strata$p_tests[1])
  label_text <- substitute(expr=paste("Nonlinear ",italic("P"), p_quad1),
                           env = base::list(p_quad1=p_quad0))
  nlmr_figure <- ggplot(fdata_sum,aes(x=predx,y=hr))+
    geom_boxplot() + 
    stat_summary(fun.y = mean, geom = "point", shape = 23, size=1) +
    
    geom_hline(yintercept =1.0,linetype='dashed',color="grey")+
    theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                     axis.text = element_text(color="black", size=tsize), 
                                     axis.title=element_text(size=10),
                                     legend.position="bottom")+
    scale_y_continuous(name="pHR of Demenita",n.breaks = 10) +
    scale_x_discrete(name="predicted Alcohol consumption (unit/week)") +
    scale_color_lancet() +
    labs(subtitle = label_text)
  nlmr_figure
}


nmlr_boxfigure_summary <- function(fdata_sum,n,tsize=12){
  fdata_sum$g <- factor(seq(1,n,1))
  nlmr_figure <- ggplot()+
    geom_point(data=fdata_sum,aes(g,bx),
               size=0.5,alpha = 1,colour="skyblue4") +
    geom_errorbar(data=fdata_sum,aes(g,ymin = bx-1.96*bxse, ymax = bx+1.96*bxse),
                  width=0.2) +
    theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                     axis.text = element_text(color="black", size=tsize),
                                     axis.title=element_text(size=12),
                                     legend.position="top") +
    
    scale_y_continuous(name="Alcohol consumption (unit/week)",n.breaks = 10) +
    scale_x_discrete(name=paste(n,"Strata by non-linear MR method")) +
    scale_color_lancet()
  nlmr_figure
}



table2_MR_all_linear <- function(md,covar_f_name,population,subgroup,IV_name,a){
  modeldata <- md
  ##orignial
  modeldata$IV <- modeldata$wPRS_95
  re_1 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup,covar_f_name,a)
  r1 <- re_1 %>% select(HR,p_value)
  r1$model <- "Model 1"
  
  ##adjusted smoke
  modeldata <- md %>% filter(is.na(smoke)==FALSE)
  modeldata$IV <- modeldata$wPRS_95
  modeldata$smoke <- factor(modeldata$smoke)
  subgroup <- "Adjusted smoke"
  covar_f_name_smoke <- c(covar_f_name,"smoke")
  re_2 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup,covar_f_name_smoke,a)
  r2 <- re_2 %>% select(HR,p_value)
  r2$model <- "Model 2 adjusted smoke"
  
  ##log(expore)
  modeldata <- md 
  modeldata$IV <- modeldata$wPRS_95
  subgroup_log <- "log-exposure"
  re_3 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup_log,covar_f_name,4)
  r3 <- re_3 %>% select(HR,p_value)
  r3$model <- "Model 3 log"
  
  ##compete risk
  modeldata <- md
  modeldata$IV <- modeldata$wPRS_95
  modeldata$death <- ifelse(is.na(modeldata$death_date)==FALSE,1,0)
  modeldata$outcome_new <- modeldata$outcome
  modeldata$outcome_new[modeldata$death==1] <- 2
  modeldata$area <- ifelse(modeldata$area=="England",1,
                           ifelse(modeldata$area=="Scotland",2,3))
  x="weight 95 PRS"
  res_compet <- IV_HR_function_compete(modeldata,covar_f_name,x)
  r4 <- res_compet %>% select(HR,p_value)
  r4$model <- "Model 4 compete risk"
  
  ##residual
  modeldata <- md
  modeldata$IV <- modeldata$wPRS_95
  modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
  fun <- paste("Exposure ~ IV +",paste(covar_f_name,collapse=" + "))
  fit <- lm(as.formula(fun), data=modeldata)
  #modeldata$predict <- predict(fit)
  r <- summary(fit)
  res <- r$coefficients[2,1]
  modeldata$predict <- modeldata$IV * res
  modeldata$residual <- modeldata$Exposure -  modeldata$predict
  ##two group
  modeldata$res_2g <- ifelse(modeldata$residual<=14,1,2)
  modeldata$res_2g_group <- as.numeric(modeldata$res_2g)
  
  fun_out <- paste("Surv(years, outcome) ~ predict +",
                   paste(covar_f_name,collapse=" + "))
  d1 <- modeldata %>% filter(res_2g_group==1)
  fun_res1 <- coxph(as.formula(fun_out),data=d1)
  fun_res1_out <- HR_output(fun_res1)
  
  r5_1 <- fun_res1_out %>% select(HR,p_value)
  r5_1$model <- "Model 5 lower residual"
  
  subgroup_res <- "res<14"
  res1 <- PRS_mr_one_HR(d1,IV_name,population,subgroup_res,covar_f_name,a)
  r6_1 <- res1 %>% select(HR,p_value)
  r6_1$model <- "Model 6 lower residual"
  
  ##compete risk for residual
  modeldata <- d1
  modeldata$death <- ifelse(is.na(modeldata$death_date)==FALSE,1,0)
  modeldata$outcome_new <- modeldata$outcome
  modeldata$outcome_new[modeldata$death==1] <- 2
  modeldata$area <- ifelse(modeldata$area=="England",1,
                           ifelse(modeldata$area=="Scotland",2,3))
  x="weight 95 PRS"
  res_compet <- IV_HR_function_compete(modeldata,covar_f_name,x)
  r7_1 <- res_compet %>% select(HR,p_value)
  r7_1$model <- "Model 7 lower risidual compete risk"
  
  d2 <- modeldata %>% filter(res_2g_group==2)
  fun_res2 <- coxph(as.formula(fun_out),data=d2)
  fun_res2_out <- HR_output(fun_res2)
  r5_2 <- fun_res2_out %>% select(HR,p_value)
  r5_2$model <- "Model 5 higher residual"
  
  subgroup_res <- "res>14"
  
  res2 <- PRS_mr_one_HR(d2,IV_name,population,subgroup_res,covar_f_name,a)
  r6_2 <- res2 %>% select(HR,p_value)
  r6_2$model <- "Model 6 higher residual"
  
  ##compete risk for residual
  modeldata <- d2
  modeldata$death <- ifelse(is.na(modeldata$death_date)==FALSE,1,0)
  modeldata$outcome_new <- modeldata$outcome
  modeldata$outcome_new[modeldata$death==1] <- 2
  modeldata$area <- ifelse(modeldata$area=="England",1,
                           ifelse(modeldata$area=="Scotland",2,3))
  x="weight 95 PRS"
  res_compet <- IV_HR_function_compete(modeldata,covar_f_name,x)
  r7_2 <- res_compet %>% select(HR,p_value)
  r7_2$model <- "Model 7 higher risidual compete risk"
  
  linear_mr_all <- rbind(r1,r2,r3,r4,r5_1,r5_2,r6_1,r6_2,r7_1,r7_2)
  linear_mr_all
}

table2_MR_linear <- function(md,covar_f_name,population,subgroup,IV_name,a){
  modeldata <- md
  ##orignial
  modeldata$IV <- modeldata$wPRS_95
  re_1 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup,covar_f_name,a)
  r1 <- re_1 %>% select(HR,p_value)
  r1$model <- "Model 1"
  
  ##adjusted smoke
  modeldata <- md %>% filter(is.na(smoke)==FALSE)
  modeldata$IV <- modeldata$wPRS_95
  modeldata$smoke <- factor(modeldata$smoke)
  subgroup <- "Adjusted smoke"
  covar_f_name_smoke <- c(covar_f_name,"smoke")
  re_2 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup,covar_f_name_smoke,a)
  r2 <- re_2 %>% select(HR,p_value)
  r2$model <- "Model 2 adjusted smoke"
  
  ##log(expore)
  modeldata <- md 
  modeldata$IV <- modeldata$wPRS_95
  subgroup_log <- "log-exposure"
  re_3 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup_log,covar_f_name,4)
  r3 <- re_3 %>% select(HR,p_value)
  r3$model <- "Model 3 log"
  
  ##compete risk
  modeldata <- md
  modeldata$IV <- modeldata$wPRS_95
  modeldata$death <- ifelse(is.na(modeldata$death_date)==FALSE,1,0)
  modeldata$outcome_new <- modeldata$outcome
  modeldata$outcome_new[modeldata$death==1] <- 2
  modeldata$area <- ifelse(modeldata$area=="England",1,
                           ifelse(modeldata$area=="Scotland",2,3))
  x="weight 95 PRS"
  res_compet <- IV_HR_function_compete(modeldata,covar_f_name,x)
  r4 <- res_compet %>% select(HR,p_value)
  r4$model <- "Model 4 compete risk"
  linear_mr_all <- rbind(r1,r2,r3,r4)
  linear_mr_all
}

table2_MR_linear_HR_OR <- function(md,covar_f_name,population,subgroup,IV_name,a){
  modeldata <- md
  ##orignial
  modeldata$IV <- modeldata$wPRS_95
  re_1 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup,covar_f_name,a)
  r1 <- re_1 %>% select(HR,p_value,OR,p_value_OR)
  r1$model <- "Model 1"
  
  ##adjusted smoke
  modeldata <- md %>% filter(is.na(smoke)==FALSE)
  modeldata$IV <- modeldata$wPRS_95
  modeldata$smoke <- factor(modeldata$smoke)
  subgroup <- "Adjusted smoke"
  covar_f_name_smoke <- c(covar_f_name,"smoke")
  re_2 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup,covar_f_name_smoke,a)
  r2 <- re_2 %>% select(HR,p_value,OR,p_value_OR)
  r2$model <- "Model 2 adjusted smoke"
  
  ##log(expore)
  modeldata <- md 
  modeldata$IV <- modeldata$wPRS_95
  subgroup_log <- "log-exposure"
  re_3 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup_log,covar_f_name,4)
  r3 <- re_3 %>% select(HR,p_value,OR,p_value_OR)
  r3$model <- "Model 3 log"
  
  linear_mr_all <- rbind(r1,r2,r3)
  linear_mr_all
}

table3_res_MR_all <- function(md,covar_f_name,population,subgroup,IV_name,a){
  ##residual
  modeldata <- md
  modeldata$IV <- modeldata$wPRS_95
  fun <- paste("Exposure ~ IV +",paste(covar_f_name,collapse=" + "))
  fit <- lm(as.formula(fun), data=modeldata)
  #modeldata$predict <- predict(fit)
  r <- summary(fit)
  res <- r$coefficients[2,1]
  modeldata$predict <- modeldata$IV * res
  modeldata$residual <- modeldata$Exposure -  modeldata$predict
  ##two group
  modeldata$res_2g <- ifelse(modeldata$residual<=14,1,2)
  modeldata$res_2g_group <- as.numeric(modeldata$res_2g)
  
  fun_out <- paste("Surv(years, outcome) ~ predict +",
                   paste(covar_f_name,collapse=" + "))
  d1 <- modeldata %>% filter(res_2g_group==1)
  fun_res1 <- coxph(as.formula(fun_out),data=d1)
  fun_res1_out <- HR_output(fun_res1)
  
  r5_1 <- fun_res1_out %>% select(HR,p_value)
  r5_1$model <- "Model 5 lower residual"
  
  subgroup_res <- "res<14"
  res1 <- PRS_mr_one_HR(d1,IV_name,population,subgroup_res,covar_f_name,a)
  r6_1 <- res1 %>% select(HR,p_value)
  r6_1$model <- "Model 6 lower residual"
  
  ##compete risk for residual
  modeldata <- d1
  modeldata$death <- ifelse(is.na(modeldata$death_date)==FALSE,1,0)
  modeldata$outcome_new <- modeldata$outcome
  modeldata$outcome_new[modeldata$death==1] <- 2
  modeldata$area <- ifelse(modeldata$area=="England",1,
                           ifelse(modeldata$area=="Scotland",2,3))
  x="weight 95 PRS"
  res_compet <- IV_HR_function_compete(modeldata,covar_f_name,x)
  r7_1 <- res_compet %>% select(HR,p_value)
  r7_1$model <- "Model 7 lower risidual compete risk"
  
  d2 <- modeldata %>% filter(res_2g_group==2)
  fun_res2 <- coxph(as.formula(fun_out),data=d2)
  fun_res2_out <- HR_output(fun_res2)
  r5_2 <- fun_res2_out %>% select(HR,p_value)
  r5_2$model <- "Model 5 higher residual"
  
  subgroup_res <- "res>14"
  
  res2 <- PRS_mr_one_HR(d2,IV_name,population,subgroup_res,covar_f_name,a)
  r6_2 <- res2 %>% select(HR,p_value)
  r6_2$model <- "Model 6 higher residual"
  
  ##compete risk for residual
  modeldata <- d2
  modeldata$death <- ifelse(is.na(modeldata$death_date)==FALSE,1,0)
  modeldata$outcome_new <- modeldata$outcome
  modeldata$outcome_new[modeldata$death==1] <- 2
  modeldata$area <- ifelse(modeldata$area=="England",1,
                           ifelse(modeldata$area=="Scotland",2,3))
  x="weight 95 PRS"
  res_compet <- IV_HR_function_compete(modeldata,covar_f_name,x)
  r7_2 <- res_compet %>% select(HR,p_value)
  r7_2$model <- "Model 7 higher risidual compete risk"
  
  #combinded the result
  res_linear_mr_all <- rbind(r5_1,r5_2,r6_1,r6_2,r7_1,r7_2)
  res_linear_mr_all
}

table3_res_MR <- function(md,covar_f_name,population,subgroup,IV_name,a){
  ##residual
  modeldata <- md
  modeldata$IV <- modeldata$wPRS_95
  fun <- paste("Exposure ~ IV +",paste(covar_f_name,collapse=" + "))
  fit <- lm(as.formula(fun), data=modeldata)
  #modeldata$predict <- predict(fit)
  r <- summary(fit)
  res <- r$coefficients[2,1]
  modeldata$predict <- modeldata$IV * res
  modeldata$residual <- modeldata$Exposure -  modeldata$predict
  ##two group
  modeldata$res_2g <- ifelse(modeldata$residual<=14,1,2)
  modeldata$res_2g_group <- as.numeric(modeldata$res_2g)
  
  fun_out <- paste("Surv(years, outcome) ~ predict +",
                   paste(covar_f_name,collapse=" + "))
  d1 <- modeldata %>% filter(res_2g_group==1)
  fun_res1 <- coxph(as.formula(fun_out),data=d1)
  fun_res1_out <- HR_output(fun_res1)
  
  r5_1 <- fun_res1_out %>% select(HR,p_value)
  r5_1$model <- "Model 5 lower residual"
  
  subgroup_res <- "res<14"
  res1 <- PRS_mr_one_HR(d1,IV_name,population,subgroup_res,covar_f_name,a)
  r6_1 <- res1 %>% select(HR,p_value)
  r6_1$model <- "Model 6 lower residual"

  
  d2 <- modeldata %>% filter(res_2g_group==2)
  fun_res2 <- coxph(as.formula(fun_out),data=d2)
  fun_res2_out <- HR_output(fun_res2)
  r5_2 <- fun_res2_out %>% select(HR,p_value)
  r5_2$model <- "Model 5 higher residual"
  
  subgroup_res <- "res>14"
  
  res2 <- PRS_mr_one_HR(d2,IV_name,population,subgroup_res,covar_f_name,a)
  r6_2 <- res2 %>% select(HR,p_value)
  r6_2$model <- "Model 6 higher residual"
  
  
  #combinded the result
  res_linear_mr_all <- rbind(r5_1,r5_2,r6_1,r6_2)
  res_linear_mr_all
}

##two sample Mr figure
mr_leaveoneout_fig <- function(leaveone_mr){
    d <-leaveone_mr
    d$up <- d$b + 1.96 * d$se
    d$lo <- d$b - 1.96 * d$se
    d$tot <- 1
    d$tot[d$SNP != "All"] <- 0.01
    d$SNP <- as.character(d$SNP)
    nom <- d$SNP[d$SNP != "All"]
    nom <- nom[order(d$b)]
    d <- rbind(d, d[nrow(d),])
    d$SNP[nrow(d)-1] <- ""
    d$b[nrow(d)-1] <- NA
    d$up[nrow(d)-1] <- NA
    d$lo[nrow(d)-1] <- NA
    d$SNP <- ordered(d$SNP, levels=c("All", "", nom))
    
    f <- ggplot(d, ggplot2::aes(y=SNP, x=b)) +
      geom_vline(xintercept=0, linetype="dotted") +
      geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
      geom_point(ggplot2::aes(colour=as.factor(tot))) +
      geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
      scale_colour_manual(values=c("black", "red")) +
      scale_size_manual(values=c(0.3, 1)) +
      theme(
        legend.position="none", 
        axis.text.y=ggplot2::element_text(size=4), 
        axis.ticks.y=ggplot2::element_line(size=0),
        axis.title.x=ggplot2::element_text(size=2)) +
      labs(y="", x=paste0("MR leave-one-out sensitivity analysis for\n", "Alcohol consumption on Dementia"))+
      theme_bw()
  f
}
