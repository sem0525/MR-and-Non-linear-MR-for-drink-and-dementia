########Table 1 characteristics###########
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)
DATA <- data_all_analysis
DATA$Exposure <- DATA$perweek_alchol_unit_cal_new
##seperate the alcohol consumption 

DATA$alcohol_2group <- ifelse(DATA$perweek_alchol_unit_cal_new<=14,0,1)
table(DATA$alcohol_2group)
cut_point = 14
DATA$alcohol_2group <- cut(DATA$Exposure,
                           breaks=c(0,cut_point,max(DATA$Exposure)),
                           include.lowest = TRUE)
###funciton
table_1_group_SMD_new <- function(conf,tdata,pre=FALSE,d=2){
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
      SMD <- SMD_continous(mean2,sd2,n2,mean1,sd1,n1,d)
      n[1,5] <- SMD
      
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
      SMD <- SMD_continous(median2,sd2,n2,median1,sd1,n1,d)
      n[1,5] <- SMD
    }
    #category for multi-categroy
    if(i %in% c("age_group","age_2group","sex","area","edu","edu_new","edu_university","townsend_group_1",
                "houseincome","smoking","smoke","BMI_group",
                "PA_group","sleep_group",
                "Diabetes","depression_g","APOE4_status",
                "alcohlol_frequency","alcohol_drinker_status",
                "family_alzheimer",
                "Diabetes","Hypertension","depression_g","APOE4_status",
                "alcohol_6group","alcohol_2group","ACD","CVD_before","stroke_before"
    )){
      
      n_dim <- dim(table(tdata$var))+1
      n_name <- names(table(tdata$var))
      
      n1 <- length(tdata$var[tdata$group==1])
      n2 <- length(tdata$var[tdata$group==2])
      
      n <- matrix(NA,nrow=n_dim,ncol=6)
      n[1,1] <- i
      n[2:n_dim,1] <- n_name
      s_table <- CrossTable(tdata$var,chisq = TRUE)
      n[2:n_dim,3] <- paste0(s_table$t," (",sprintf("%0.1f",s_table$prop.row*100),")")
      s_table <- CrossTable(tdata$var,tdata$group,chisq = TRUE)
      
      t <- s_table$t
      SMD <- apply(t,1,function(x){
        SMD_category(x[1],n1,x[2],n2,d=d)
      })
      
      n[2:n_dim,4] <- paste0(s_table$t[,1]," (",sprintf("%0.1f",s_table$prop.col[,1]*100),")")
      n[2:n_dim,5] <- paste0(s_table$t[,2]," (",sprintf("%0.1f",s_table$prop.col[,2]*100),")")
      n[2:n_dim,6] <- SMD
      
      
      if(i=="age_group"){n[,2] <- c("Age, n (%)",
                                    "  ≤45 year","  (45, 65] year","  > 65 year")}
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
                                    "  ≥30.0 kg/m^2")}
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
      if(i=="townsend_group_1"){n[,2] <- c("Socioeconomic, n (%)",
                                           "  Least deprived",
                                           "  Middle deprived",
                                           "  Most deprived")}
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
      if(i=="PA_group"){n[,2] <- c("Physical Activity, n (%)",
                                   "  Insufficient","  Sufficient","  Additional")}
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
      
      if(i=="APOE4_status"){n[,2] <- c("APOE ε4, n (%)","  Without","  With")}
      if(i=="alcohol_2group"){n[,2] <- c("Alcohol consumption, n (%)",
                                         "  Safe",
                                         "  Unsafe")}
      if(i=="alcohol_6group"){n[,2] <- c("Alcohol consumption, n (%)",
                                         "  Light (≤7 unit/week)",
                                         "  Moderate ((7, 14] unit/week)",
                                         "  Heavy ((14, 21] unit/week)",
                                         "  Very heavy ((21, 28] unit/week)",
                                         "  Excessive ((28, 42] unit/week)",
                                         "  Extremely excessive (>42 unit/week)")}
      ##disease
      if(i=="ACD"){n[,2] <- c("All-cause dementia, n (%)","  Without","  With")}
      if(i=="CVD_before"){n[,2] <- c("Comorbid disease at baseline, n (%)","  NA","  Cardiometabolic disease")}
      if(i=="stroke_before"){n[,2] <- c("Comorbid disease at baseline, n (%)","  NA","  Stroke")}
      n <- n[,c(2:6)]
    }
    
    N <- rbind(N,n)
  }
  N
}


########################Table 1: Descritpitive ##################
d <- DATA

conf <-  c("years","perweek_alchol_unit_cal_new","sex",
           "age_group","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","CVD_before","stroke_before",
           "APOE4_status") 
conf <-  c("years","perweek_alchol_unit_cal_new","sex",
           "age","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","CVD_before","stroke_before",
           "APOE4_status") 


d$group <- d$sex+1
table1_result <- table_1_group_SMD_new(conf,d,FALSE)

table1_result <- table1_result[-37,]
table1_result <- table1_result[-39,]
table1_result <- table1_result[-38,]

View(table1_result)
