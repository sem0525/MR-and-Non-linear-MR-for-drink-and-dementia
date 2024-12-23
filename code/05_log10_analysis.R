##############Analysis and sub into sex group log-transform#############
###Analysis###
project_path <- "/home/zoo/project/P1_drinkMR_dementia/"
##result save date
Date <- "_drinker_20231122.csv"
###load the function
source(paste0(project_path,"code/00_functions.R"))
###load the packages
library(dplyr)
library(survival)
library(ggsci)
library(nnet)
library(rms)
library(gridExtra)
library(gmodels)
###load the data
load(paste0(project_path,"data/data_alldrink_407543.RData")) ##407983
data <- data_alldrink
######Supplement Table description the zero Drinkers######
d <- data
d$G <- d$alcohol_drinker_status
table(d$G)
d$G[d$alcohol_drinker_status==2 & 
        d$perweek_alchol_unit_cal_new==0] <- 2
d$G[d$alcohol_drinker_status==2 & 
          d$perweek_alchol_unit_cal_new >0 & 
          d$perweek_alchol_unit_cal_new <=14 ] <- 3
d$G[d$alcohol_drinker_status==2 & 
          d$perweek_alchol_unit_cal_new >14 ] <- 4
conf <-  c("years","perweek_alchol_unit_cal_new",
           "age_group","sex","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","APOE4_status","ACD"
) 
table(d$G)
##non drinker, ex-drinker
d$group <- NA
d$group[d$G==0] <- 1 
d$group[d$G==1] <- 2
table(is.na(d$group))
d1 <- d %>% filter(is.na(group)==FALSE)
re1 <- table_1_group_SMD(conf,d1,FALSE)
##dinker=0, drinkerless14
d$group <- NA
d$group[d$G==2] <- 1 
d$group[d$G==3] <- 2
table(is.na(d$group))
d1 <- d %>% filter(is.na(group)==FALSE)
re2 <- table_1_group_SMD(conf,d1,FALSE)
##dinker=0, drinker14
d$group <- NA
d$group[d$G==2] <- 1 
d$group[d$G==4] <- 2
table(is.na(d$group))
d1 <- d %>% filter(is.na(group)==FALSE)
re3 <- table_1_group_SMD(conf,d1,FALSE)
####all and current dink
d$group <- NA
d$group <- ifelse(d$alcohol_drinker_status==2,1,2)
re4 <- table_1_group_SMD(conf,d,FALSE)
###combinded
Re_sup <- data.frame(re4[,c(1,2)],re1[,c(3,4)],re4[,3], re2[,c(3,4)],re3[,4])

names(Re_sup) <- c("Characteristics","Total","Non-drinkers","Ex-drinkers","All current drinkers",
                   "0 unit/week","<=14 unit/week",">14 unit/week")

#####SupTable_all_description######
write.csv(Re_sup,paste0(project_path,"result/sTable_characterist",Date))
#################sTale: PRS with assocition with X,Y and confoudners#######
######IV continous
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==0) 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
cof <- c("townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","APOE4","outcome","alcohol_liver_alcoholliver")
Result_q_continous_abstainer <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$low),", ",
                     sprintf("%0.3f",res$up),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){names(res_var) <- names(Result_q_continous_abstainer)}
  Result_q_continous_abstainer <- rbind(Result_q_continous_abstainer,res_var)
}
Result_q_continous_abstainer
Result_q_continous_abstainer$test <- ifelse(Result_q_continous_abstainer$`Pr(>|t|)`<(0.05/9),1,0)
##men
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==0,sex==1) 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
cof <- c("townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","APOE4","outcome","alcohol_liver_alcoholliver")
Result_q_continous_abstainer_men <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$low),", ",
                     sprintf("%0.3f",res$up),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){names(res_var) <- names(Result_q_continous_abstainer_men)}
  Result_q_continous_abstainer_men <- rbind(Result_q_continous_abstainer_men,res_var)
}
Result_q_continous_abstainer_men
Result_q_continous_abstainer_men$test <- ifelse(Result_q_continous_abstainer_men$`Pr(>|t|)`<(0.05/9),1,0)
##women
modeldata <- data_alldrink %>% filter(alcohol_drinker_status==0,sex==0) 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
cof <- c("townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","APOE4","outcome","alcohol_liver_alcoholliver")
Result_q_continous_abstainer_women <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$low),", ",
                     sprintf("%0.3f",res$up),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){names(res_var) <- names(Result_q_continous_abstainer_women)}
  Result_q_continous_abstainer_women <- rbind(Result_q_continous_abstainer_women,res_var)
}
Result_q_continous_abstainer_women
Result_q_continous_abstainer_women$test <- ifelse(Result_q_continous_abstainer_women$`Pr(>|t|)`<(0.05/9),1,0)

###########################################
#######All alcohol consumption is zero
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                     data$perweek_alchol_unit_cal_new>0)
DATA <- data_all_analysis
##target exposure
DATA$Exposure <- log10(DATA$perweek_alchol_unit_cal_new)
hist(DATA$Exposure )
cut_point <- log10(14)
##seperate the alcohol consumption 
DATA$alcohol_2group <- cut(DATA$Exposure,
                           breaks=c(0,cut_point,max(DATA$Exposure)),
                           include.lowest = TRUE)
table(DATA$alcohol_2group)
##seperate the alcohol consumption 
alcohol_quintile <- quantile(DATA$Exposure,probs=seq(0,1,0.2))
DATA$alcohol_5group <- cut(DATA$Exposure,
                           breaks=alcohol_quintile,
                           include.lowest = TRUE)
table(DATA$alcohol_5group)
##q10
alcohol_quintile <- quantile(DATA$Exposure,probs=seq(0,1,0.1))
DATA$alcohol_10group <- cut(DATA$Exposure,
                            breaks=alcohol_quintile,
                            include.lowest = TRUE)
table(DATA$alcohol_10group)
##3 group
DATA$alcohol_per7_3g <- cut(DATA$Exposure,
                             breaks=c(seq(min(DATA$Exposure),max(DATA$Exposure)-0.2,cut_point),max(DATA$Exposure)),
                             include.lowest = TRUE)
table(DATA$alcohol_per7_3g)
###outcome
DATA$outcome <- DATA$ACD

#####PART One#############################################################
########################Table 1: Descritpitive ##################
d <- DATA
conf <-  c("years","perweek_alchol_unit_cal_new",
           "age_group","sex","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","APOE4_status","ACD"
) 
d$group <- as.numeric(d$alcohol_2group)
re <- table_1_group_SMD(conf,d,FALSE)
#####table1_all######
write.csv(re,paste0(project_path,"result/Table1_all",Date))
#####subgroup to sex######
conf <-  c("years","perweek_alchol_unit_cal_new",
           "age_group","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","APOE4_status","ACD") 
##men
d_men <- d %>% filter(sex==1)
re_men <- table_1_group_SMD(conf,d_men,FALSE)
##women
d_women <- d %>% filter(sex==0)
re_women <- table_1_group_SMD(conf,d_women,FALSE)
re_all <- cbind(re_men,re_women)
#####table1_sugbroup######
write.csv(re_all,paste0(project_path,"result/Table1_sub",Date))


#####PART Two#############################################################
###Figure table 2: the table of incident across alcohol ##################
modeldata <- DATA
modeldata$E_group <- modeldata$alcohol_5group
table2_all <- table_2_incidence(modeldata)
Overall <- rep(NA,5)
alcohol_quintile <- data.frame(quantile(DATA$Exposure,probs=seq(0,1,0.2)))
names(alcohol_quintile) <- "log10_E"
alcohol_quintile$exposure <- 10^alcohol_quintile$log10_E
Exposure_unit <- c(paste0("[",sprintf("%0.2f",alcohol_quintile$exposure[1]), 
                          " ,",sprintf("%0.2f",alcohol_quintile$exposure[2]),"]"
                          ),
                   paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[2]), 
                          " ,",sprintf("%0.2f",alcohol_quintile$exposure[3]),"]"
                   ),
                   paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[3]), 
                          " ,",sprintf("%0.2f",alcohol_quintile$exposure[4]),"]"
                   ),
                   paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[4]), 
                          " ,",sprintf("%0.2f",alcohol_quintile$exposure[5]),"]"
                   ),
                   paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[5]), 
                          " ,",sprintf("%0.2f",alcohol_quintile$exposure[6]),"]"
                   )
                   )

table2_all <- rbind(Overall,Exposure_unit,table2_all)
###man
modeldata <- DATA %>% filter(sex==1)
alcohol_quintile <- quantile(modeldata$Exposure,probs=seq(0,1,0.2))
modeldata$E_group <- cut(modeldata$Exposure,
                         breaks=alcohol_quintile,
                         include.lowest = TRUE)
alcohol_quintile <- data.frame(alcohol_quintile)
names(alcohol_quintile) <- "log10_E"
alcohol_quintile$exposure <- 10^alcohol_quintile$log10_E
Exposure_unit <- c(paste0("[",sprintf("%0.2f",alcohol_quintile$exposure[1]), 
                          " ,",sprintf("%0.2f",alcohol_quintile$exposure[2]),"]"
                          ),
                          paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[2]), 
                                 " ,",sprintf("%0.2f",alcohol_quintile$exposure[3]),"]"
                          ),
                          paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[3]), 
                                 " ,",sprintf("%0.2f",alcohol_quintile$exposure[4]),"]"
                          ),
                          paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[4]), 
                                 " ,",sprintf("%0.2f",alcohol_quintile$exposure[5]),"]"
                          ),
                          paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[5]), 
                                 " ,",sprintf("%0.2f",alcohol_quintile$exposure[6]),"]"
                          )
                    )

table2_men <- table_2_incidence(modeldata)
Men <- rep(NA,5)
table2_men <- rbind(Men,Exposure_unit,table2_men)

###women
modeldata <- DATA %>% filter(sex==0)
alcohol_quintile <- quantile(modeldata$Exposure,probs=seq(0,1,0.2))
modeldata$E_group <- cut(modeldata$Exposure,
                         breaks=alcohol_quintile,
                         include.lowest = TRUE)
alcohol_quintile <- data.frame(alcohol_quintile)
names(alcohol_quintile) <- "log10_E"
alcohol_quintile$exposure <- 10^alcohol_quintile$log10_E
Exposure_unit <- c(paste0("[",sprintf("%0.2f",alcohol_quintile$exposure[1]), 
                  " ,",sprintf("%0.2f",alcohol_quintile$exposure[2]),"]"
                  ),
                  paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[2]), 
                         " ,",sprintf("%0.2f",alcohol_quintile$exposure[3]),"]"
                  ),
                  paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[3]), 
                         " ,",sprintf("%0.2f",alcohol_quintile$exposure[4]),"]"
                  ),
                  paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[4]), 
                         " ,",sprintf("%0.2f",alcohol_quintile$exposure[5]),"]"
                  ),
                  paste0("(",sprintf("%0.2f",alcohol_quintile$exposure[5]), 
                         " ,",sprintf("%0.2f",alcohol_quintile$exposure[6]),"]"
                  ))

table2_women <- table_2_incidence(modeldata)
Women <- rep(NA,5)
table2_women <- rbind(Women,Exposure_unit,table2_women)
inc_d_case <- data.frame(rbind(table2_all,table2_men,table2_women))
names(inc_d_case) <- paste0("Q",1:5)
######figure table 1.2#####
write.csv(inc_d_case,paste0(project_path,"result/fig1table",Date))
########################Figure 2: HR Epidemiological Association ##################
##########overall#############
modeldata <- DATA
covar_f_name <- c("age","sex","APOE4_status","area","townsend_group_1","edu_new")
##rcs model
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
ns <- 6
coxplot_data <- cbind(coxplot_data[,1:ns],c_d1, c_d2, c_d3) 
names(coxplot_data)[(ns+1):(ns+2)] <- paste0("edu",1:2)
names(coxplot_data)[(ns+5):(ns+6)] <- paste0("townsend",1:2)
fun <- paste("Surv(years,outcome) ~ rcs(Exposure,3) + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
fun1 <- paste("Surv(years,outcome) ~ Exposure + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
dd <- datadist(coxplot_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=coxplot_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=coxplot_data)
dd$limits$Exposure[2] <- 1.084994
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
ref_line <- HR$Exposure[HR$yhat==min(HR$yhat)]
ref_exposure = 10^1.084994 #1.084994
g0 <- ggplot() +
  geom_line(data=HR, aes(Exposure,yhat),
            linetype="solid",size=0.5,alpha = 1,colour="skyblue4") +
  geom_ribbon(data=HR, aes(Exposure,ymin = lower, ymax = upper),
              alpha = 0.6,fill="skyblue1") +
  geom_vline(xintercept = ref_line,color="grey60",linetype="dashed")+
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=12),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  geom_hline(yintercept=1, linetype="dotted", color = "grey50") + 
  scale_x_continuous(name=expression("Alcohol cunsumption" ~log[10]*"(unit/week)"),
                     expand = c(0,0.03),n.breaks=10) + 
  scale_y_continuous(name="HR of Dementia in Overall",expand = c(0,0.03),
                     #limits=c(0.9,2),
                     breaks=c(1,1.1,1.2,1.3,1.4,1.5,1.6)) +
  scale_fill_jama() +
  geom_text(aes(x=-0.3, y=1.58, 
                label=paste0("Ref. ",sprintf("%0.2f",ref_exposure)," unit/week")), 
            position=position_nudge(x=0.2), hjust=0)
g0
######figure#####
G_all <- g0 + labs(subtitle = label_text)
#######rcs for overall ###########
G_all
############men##########
modeldata <- DATA %>% filter(sex==1)
covar_f_name <- c("age","APOE4_status","area","townsend_group_1","edu_new")
##rcs model
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
ns <- 5
coxplot_data <- cbind(coxplot_data[,1:ns],c_d1, c_d2, c_d3) 
names(coxplot_data)[(ns+1):(ns+2)] <- paste0("edu",1:2)
names(coxplot_data)[(ns+5):(ns+6)] <- paste0("townsend",1:2)
fun <- paste("Surv(years,outcome) ~ rcs(Exposure,3) + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
fun1 <- paste("Surv(years,outcome) ~ Exposure + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
dd <- datadist(coxplot_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=coxplot_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=coxplot_data)
dd$limits$Exposure[2] <- 1.225477
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
p_allover <- p_test[1,3]
p_allover <- p_out(p_allover)
p_nonliner <- p_out_fig(p_nonliner)
HR <- Predict(m0,Exposure, fun=exp,ref.zero=TRUE)
#plot2
label_text_men <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                         env = base::list(p_nonliner1=p_nonliner))
ref_line <- HR$Exposure[HR$yhat==min(HR$yhat)]
ref_exposure_men = 10^1.225477  #1.225477
g0_men <- ggplot() +
  geom_line(data=HR, aes(Exposure,yhat),
            linetype="solid",size=0.5,alpha = 1,colour="skyblue4") +
  geom_ribbon(data=HR, aes(Exposure,ymin = lower, ymax = upper),
              alpha = 0.6,fill="skyblue1") +
  geom_vline(xintercept = ref_line,color="grey60",linetype="dashed")+
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=12),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  geom_hline(yintercept=1, linetype="dotted", color = "grey50") + 
  scale_x_continuous(name=expression("Alcohol cunsumption" ~log[10]*"(unit/week)"),
                     expand = c(0,0.03),n.breaks=10) + 
  scale_y_continuous(name="HR of Dementia in Men",expand = c(0,0.03),
                     #limits=c(0.9,2),
                     breaks=c(1,1.2,1.4,1.6,1.8,2.0)) +
  scale_fill_jama() +
  geom_text(aes(x=-0.2, y=1.9, 
                label=paste0("Ref. ",sprintf("%0.2f",ref_exposure_men)," unit/week")), 
            position=position_nudge(x=0.2), hjust=0)
g0_men
######figure#####
G_men <- g0_men + labs(subtitle = label_text_men)
#######rcs for overall ###########
G_men
###########women##########
modeldata <- DATA %>% filter(sex==0)
covar_f_name <- c("age","APOE4_status","area","townsend_group_1","edu_new")
##rcs model
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
ns <- 5
coxplot_data <- cbind(coxplot_data[,1:ns],c_d1, c_d2, c_d3) 
names(coxplot_data)[(ns+1):(ns+2)] <- paste0("edu",1:2)
names(coxplot_data)[(ns+5):(ns+6)] <- paste0("townsend",1:2)
fun <- paste("Surv(years,outcome) ~ rcs(Exposure,3) + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
fun1 <- paste("Surv(years,outcome) ~ Exposure + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
dd <- datadist(coxplot_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=coxplot_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=coxplot_data)
dd$limits$Exposure[2] <- 0.9222726 ##need double analysis
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
p_allover <- p_test[1,3]
p_allover <- p_out(p_allover)
p_nonliner <- p_out_fig(p_nonliner)
HR <- Predict(m0,Exposure, fun=exp,ref.zero=TRUE)
#plot2
label_text_women <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                         env = base::list(p_nonliner1=p_nonliner))
ref_line <- HR$Exposure[HR$yhat==min(HR$yhat)]
ref_exposure_women = 10^0.9222726  #0.9222726
g0_women <- ggplot() +
  geom_line(data=HR, aes(Exposure,yhat),
            linetype="solid",size=0.5,alpha = 1,colour="skyblue4") +
  geom_ribbon(data=HR, aes(Exposure,ymin = lower, ymax = upper),
              alpha = 0.6,fill="skyblue1") +
  geom_vline(xintercept = ref_line,color="grey60",linetype="dashed")+
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=12),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  geom_hline(yintercept=1, linetype="dotted", color = "grey50") + 
  scale_x_continuous(name=expression("Alcohol cunsumption" ~log[10]*"(unit/week)"),
                     expand = c(0,0.03),n.breaks=10) + 
  scale_y_continuous(name="HR of Dementia in Women",expand = c(0,0.03),
                     #limits=c(0.9,2),
                     breaks=c(1,1.1,1.2,1.3,1.4,1.5,1.6)) +
  scale_fill_jama() +
  geom_text(aes(x=-0.32, y=1.65, 
                label=paste0("Ref. ",sprintf("%0.2f",ref_exposure_women)," unit/week")), 
            position=position_nudge(x=0.2), hjust=0)

######figure#####
G_women <- g0_women + labs(subtitle = label_text_women)
#######rcs for overall ###########
G_women
######figure#####
grid.arrange(G_all,G_men,G_women,ncol=3)
############
############rcs transform the x-aix############
##########overall#############
modeldata <- DATA
covar_f_name <- c("age","sex","APOE4_status","area","townsend_group_1","edu_new")
##rcs model
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
ns <- 6
coxplot_data <- cbind(coxplot_data[,1:ns],c_d1, c_d2, c_d3) 
names(coxplot_data)[(ns+1):(ns+2)] <- paste0("edu",1:2)
names(coxplot_data)[(ns+5):(ns+6)] <- paste0("townsend",1:2)
fun <- paste("Surv(years,outcome) ~ rcs(Exposure,3) + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
fun1 <- paste("Surv(years,outcome) ~ Exposure + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
dd <- datadist(coxplot_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=coxplot_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=coxplot_data)
dd$limits$Exposure[2] <- 1.084994
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
ref_line <- HR$Exposure[HR$yhat==min(HR$yhat)]
ref_exposure = 10^1.084994 #1.084994
head(HR)
HR$X <- 10^HR$Exposure
g0 <-ggplot() +
  geom_line(data=HR, aes(X,yhat),
            linetype="solid",size=0.5,alpha = 1,colour="skyblue4") +
  geom_ribbon(data=HR, aes(X,ymin = lower, ymax = upper),
              alpha = 0.6,fill="skyblue1") +
  geom_vline(xintercept = ref_exposure,color="grey60",linetype="dashed")+
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=12),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  geom_hline(yintercept=1, linetype="dotted", color = "grey50") + 
  scale_x_log10(name="Alcohol cunsumption (unit/week)",
                expand = c(0,0.03),
                breaks= round(10^seq(-0.4,2.4,0.4),1)) +
  scale_y_continuous(name="HR of Dementia in Overall",expand = c(0,0.03),
                     #limits=c(0.9,2),
                     breaks=c(1,1.1,1.2,1.3,1.4,1.5,1.6)
  ) +
  scale_fill_jama() +
  geom_text(aes(x=0.55, y=1.59, 
                label=paste0("Ref. ",sprintf("%0.1f",ref_exposure)," unit/week")), 
            position=position_nudge(x=0.2), hjust=0)
######figure#####
G_all_1 <- g0 + labs(subtitle = label_text)
#######rcs for overall ###########
G_all_1
##########men#############
modeldata <- DATA %>% filter(sex==1)
covar_f_name <- c("age","APOE4_status","area","townsend_group_1","edu_new")
##rcs model
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
ns <- 5
coxplot_data <- cbind(coxplot_data[,1:ns],c_d1, c_d2, c_d3) 
names(coxplot_data)[(ns+1):(ns+2)] <- paste0("edu",1:2)
names(coxplot_data)[(ns+5):(ns+6)] <- paste0("townsend",1:2)
fun <- paste("Surv(years,outcome) ~ rcs(Exposure,3) + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
fun1 <- paste("Surv(years,outcome) ~ Exposure + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
dd <- datadist(coxplot_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=coxplot_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=coxplot_data)
dd$limits$Exposure[2] <- 1.225477
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
p_allover <- p_test[1,3]
p_allover <- p_out(p_allover)
p_nonliner <- p_out_fig(p_nonliner)
HR <- Predict(m0,Exposure, fun=exp,ref.zero=TRUE)
#plot2
label_text_men <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                             env = base::list(p_nonliner1=p_nonliner))
ref_line <- HR$Exposure[HR$yhat==min(HR$yhat)]
ref_exposure_men = 10^1.225477  #1.225477
HR$X <- 10^HR$Exposure
g0 <-ggplot() +
  geom_line(data=HR, aes(X,yhat),
            linetype="solid",size=0.5,alpha = 1,colour="skyblue4") +
  geom_ribbon(data=HR, aes(X,ymin = lower, ymax = upper),
              alpha = 0.6,fill="skyblue1") +
  geom_vline(xintercept = ref_exposure_men,color="grey60",linetype="dashed")+
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=12),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  geom_hline(yintercept=1, linetype="dotted", color = "grey50") + 
  scale_x_log10(name="Alcohol cunsumption (unit/week)",
                expand = c(0,0.03),
                breaks= round(10^seq(-0.4,2.4,0.4),1)) + 
  scale_y_continuous(name="HR of Dementia in Men",expand = c(0,0.03),
                     #limits=c(0.9,2),
                     breaks=c(1,1.2,1.4,1.6,1.8,2.0)
  ) +
  scale_fill_jama() +
  geom_text(aes(x=0.75, y=1.92, 
                label=paste0("Ref. ",sprintf("%0.1f",ref_exposure_men)," unit/week")), 
            position=position_nudge(x=0.2), hjust=0)
######figure#####
G_men_1 <- g0 + labs(subtitle = label_text_men)
#######rcs for men ###########
G_men_1
##########women#############
modeldata <- DATA %>% filter(sex==0)
covar_f_name <- c("age","APOE4_status","area","townsend_group_1","edu_new")
##rcs model
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
ns <- 5
coxplot_data <- cbind(coxplot_data[,1:ns],c_d1, c_d2, c_d3) 
names(coxplot_data)[(ns+1):(ns+2)] <- paste0("edu",1:2)
names(coxplot_data)[(ns+5):(ns+6)] <- paste0("townsend",1:2)
fun <- paste("Surv(years,outcome) ~ rcs(Exposure,3) + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
fun1 <- paste("Surv(years,outcome) ~ Exposure + ",paste(names(coxplot_data)[4:dim(coxplot_data)[2]],collapse=" + "))
dd <- datadist(coxplot_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=coxplot_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=coxplot_data)
dd$limits$Exposure[2] <- 0.9222726 ##need double analysis
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
p_allover <- p_test[1,3]
p_allover <- p_out(p_allover)
p_nonliner <- p_out_fig(p_nonliner)
HR <- Predict(m0,Exposure, fun=exp,ref.zero=TRUE)
#plot2
label_text_women <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                               env = base::list(p_nonliner1=p_nonliner))
ref_line <- HR$Exposure[HR$yhat==min(HR$yhat)]
ref_exposure_women = 10^0.9222726  #0.9222726
HR$X <- 10^HR$Exposure
g0 <-ggplot() +
  geom_line(data=HR, aes(X,yhat),
            linetype="solid",size=0.5,alpha = 1,colour="skyblue4") +
  geom_ribbon(data=HR, aes(X,ymin = lower, ymax = upper),
              alpha = 0.6,fill="skyblue1") +
  geom_vline(xintercept = ref_exposure_women,color="grey60",linetype="dashed")+
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=12),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  geom_hline(yintercept=1, linetype="dotted", color = "grey50") + 
  scale_x_log10(name="Alcohol cunsumption (unit/week)",
                expand = c(0,0.03),
                #breaks= round(10^seq(-0.4,2.4,0.4),1)) + 
                breaks= round(10^seq(-0.4,2.0,0.4),1)) + 
  scale_y_continuous(name="HR of Dementia in Women",expand = c(0,0.03),
                     #limits=c(0.9,2),
                     breaks=c(1,1.2,1.4,1.6,1.8,2.0)
  ) +
  scale_fill_jama() +
  geom_text(aes(x=0.57, y=1.66, 
                label=paste0("Ref. ",sprintf("%0.1f",ref_exposure_women)," unit/week")), 
            position=position_nudge(x=0.2), hjust=0)
######figure#####
G_women_1 <- g0 + labs(subtitle = label_text_women)
#######rcs for women ###########
G_women_1
######figure: rcs transform the x-aix#####
grid.arrange(G_all_1,G_men_1,G_women_1,ncol=3)

#####PART Three#############################################################
######figure: PRS associated with exposure######
######oveall######
modeldata <- DATA 
modeldata$Exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
modeldata$IV <- modeldata$wPRS_95
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
modeldata$IV_group <- cut(modeldata$IV,include.lowest = T,
                          breaks=quantile(modeldata$IV,probs= seq(0,1,0.2)))
g=5
table(modeldata$IV_group)
modeldata$group <- as.numeric(modeldata$IV_group)
modeldata$Var <- modeldata$Exposure
fun <- paste("Var ~ IV_group +", paste(covar_f_name,collapse = " + "))
fit <- glm(as.formula(fun), data=modeldata)
r <- summary(fit)
res <- as.data.frame(r$coefficients[1:g,1:2])
res$low <- res$Estimate-1.96*res$`Std. Error`
res$up <- res$Estimate+1.96*res$`Std. Error`
res[1,] <- c(0,0,0,0)
res$g <- factor(paste0("Q",1:g),levels=paste0("Q",1:g))
res$group <- 1:5
res_g <- modeldata %>% group_by(group) %>% mutate(mean_PRS=mean(IV,na.rm=T)) %>% select(group,mean_PRS)
res_g <- unique(res_g)
res <- merge(res,res_g,by="group")
g1 <- ggplot(res,aes(x=mean_PRS,y=Estimate  ))+
  geom_hline(yintercept = 0,
             #linetype='dashed',
             color="grey",alpha = 0.2)+
  geom_vline(xintercept = res_g$mean_PRS,linetype='dotted',color="grey")+
  geom_point(position=position_dodge(0)) +
  geom_errorbar(aes(ymin=low,ymax=up),
                position=position_dodge(0),width=0.001)+
  #labs(subtitle = label_text) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=10),
                                   legend.position="bottom")+
  scale_y_continuous(name=expression("Coefficient of Alcohol consumption"~log[10]*"(unit/week)"),
                     #limits=c(0,6),
                     n.breaks = 5) +
  xlab("Decile of PRS") +
  scale_color_lancet() 


g2 <- ggplot(modeldata, aes(x=IV)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=10,), 
        axis.title=element_text(size=10))+
  geom_histogram(binwidth=0.02,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.5) + 
  scale_x_continuous(name="Alcohol consumption weight PRS",expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="Overall",expand = c(0,0),
                     limit=c(0,41000),
                     n.breaks=10) +
  scale_fill_jama()

plot_PRS_q5_all <- ggplot2.two_y_axis(g2,g1)
####plot_PRS_all####
grid.arrange(plot_PRS_q5_all)
#####Men#######
modeldata <- DATA %>% filter(sex==1)
modeldata$Exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
modeldata$IV <- modeldata$wPRS_95
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
modeldata$IV_group <- cut(modeldata$IV,include.lowest = T,
                          breaks=quantile(modeldata$IV,probs= seq(0,1,0.2)))
g=5
table(modeldata$IV_group)
modeldata$group <- as.numeric(modeldata$IV_group)
modeldata$Var <- modeldata$Exposure
fun <- paste("Var ~ IV_group +", paste(covar_f_name,collapse = " + "))
fit <- glm(as.formula(fun), data=modeldata)
r <- summary(fit)
res <- as.data.frame(r$coefficients[1:g,1:2])
res$low <- res$Estimate-1.96*res$`Std. Error`
res$up <- res$Estimate+1.96*res$`Std. Error`
res[1,] <- c(0,0,0,0)
res$g <- factor(paste0("Q",1:g),levels=paste0("Q",1:g))
res$group <- 1:5
res_g <- modeldata %>% group_by(group) %>% mutate(mean_PRS=mean(IV,na.rm=T)) %>% select(group,mean_PRS)
res_g <- unique(res_g)
res <- merge(res,res_g,by="group")
g1 <- ggplot(res,aes(x=mean_PRS,y=Estimate  ))+
  geom_hline(yintercept = 0,
             #linetype='dashed',
             color="grey",alpha = 0.2)+
  geom_vline(xintercept = res_g$mean_PRS,linetype='dotted',color="grey")+
  geom_point(position=position_dodge(0)) +
  geom_errorbar(aes(ymin=low,ymax=up),
                position=position_dodge(0),width=0.001)+
  #labs(subtitle = label_text) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=10),
                                   legend.position="bottom")+
  scale_y_continuous(name=expression("Coefficient of Alcohol consumption"~log[10]*"(unit/week)"),
                     #limits=c(0,6),
                     breaks=c(0.000,0.025,0.050,0.075,0.100)) +
  xlab("Decile of PRS") +
  scale_color_lancet() 


g2 <- ggplot(modeldata, aes(x=IV)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=10,), 
        axis.title=element_text(size=10))+
  geom_histogram(binwidth=0.02,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.5) + 
  scale_x_continuous(name="Alcohol consumption weight PRS",expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="Men",expand = c(0,0),
                     limit=c(0,20000),
                     n.breaks=10) +
  scale_fill_jama()

plot_PRS_q5_men <- ggplot2.two_y_axis(g2,g1)
####plot_PRS_men####
grid.arrange(plot_PRS_q5_men)
#####Women#######
modeldata <- DATA %>% filter(sex==0)
modeldata$Exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
modeldata$IV <- modeldata$wPRS_95
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
modeldata$IV_group <- cut(modeldata$IV,include.lowest = T,
                          breaks=quantile(modeldata$IV,probs= seq(0,1,0.2)))
g=5
table(modeldata$IV_group)
modeldata$group <- as.numeric(modeldata$IV_group)
modeldata$Var <- modeldata$Exposure
fun <- paste("Var ~ IV_group +", paste(covar_f_name,collapse = " + "))
fit <- glm(as.formula(fun), data=modeldata)
r <- summary(fit)
res <- as.data.frame(r$coefficients[1:g,1:2])
res$low <- res$Estimate-1.96*res$`Std. Error`
res$up <- res$Estimate+1.96*res$`Std. Error`
res[1,] <- c(0,0,0,0)
res$g <- factor(paste0("Q",1:g),levels=paste0("Q",1:g))
res$group <- 1:5
res_g <- modeldata %>% group_by(group) %>% mutate(mean_PRS=mean(IV,na.rm=T)) %>% select(group,mean_PRS)
res_g <- unique(res_g)
res <- merge(res,res_g,by="group")
g1 <- ggplot(res,aes(x=mean_PRS,y=Estimate  ))+
  geom_hline(yintercept = 0,
             #linetype='dashed',
             color="grey",alpha = 0.2)+
  geom_vline(xintercept = res_g$mean_PRS,linetype='dotted',color="grey")+
  geom_point(position=position_dodge(0)) +
  geom_errorbar(aes(ymin=low,ymax=up),
                position=position_dodge(0),width=0.001)+
  #labs(subtitle = label_text) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=10),
                                   legend.position="bottom")+
  scale_y_continuous(name=expression("Coefficient of Alcohol consumption"~log[10]*"(unit/week)"),
                     #limits=c(0,6),
                     n.breaks = 5) +
  xlab("Decile of PRS") +
  scale_color_lancet() 
g2 <- ggplot(modeldata, aes(x=IV)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=10,), 
        axis.title=element_text(size=10))+
  geom_histogram(binwidth=0.02,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.5) + 
  scale_x_continuous(name="Alcohol consumption weight PRS",expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="Women",expand = c(0,0),
                     limit=c(0,21000),
                     n.breaks=10) +
  scale_fill_jama()

plot_PRS_q5_women <- ggplot2.two_y_axis(g2,g1)
####plot_PRS_women####
grid.arrange(plot_PRS_q5_women)
###Plot_PRS##########
grid.arrange(plot_PRS_q5_all,plot_PRS_q5_men,plot_PRS_q5_women,ncol=3)

#####PART Four: IV information#############################################################
#################sTale: PRS with assocition with X,Y and confoudners#######
######IV continous
modeldata <- DATA 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
cof <- c("Exposure","Exposure_log","townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","APOE4","outcome","alcohol_liver_alcoholliver")
Result_q_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$low),", ",
                     sprintf("%0.3f",res$up),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){names(res_var) <- names(Result_q_continous)}
  Result_q_continous <- rbind(Result_q_continous,res_var)
}
Result_q_continous
Result_q_continous$test <- ifelse(Result_q_continous$`Pr(>|t|)`<(0.05/11),1,0)
##adjusted
covar_f_name_adj <- c("Exposure_log","age","sex","area","array",paste0("component",1:20))
cof <- c("townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","APOE4","outcome","alcohol_liver_alcoholliver")

Result_q_continous_adjust <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name_adj,collapse=" + "))
  
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$low),", ",
                     sprintf("%0.3f",res$up),")")
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){res$p <- p_out(res$`Pr(>|z|)`)}
  if(!(xvar %in% c("outcome","alcohol_liver_alcoholliver"))){res$p <- p_out(res$`Pr(>|t|)`)}
  res_var <- cbind(xvar, res)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){names(res_var) <- names(Result_q_continous_adjust)}
  Result_q_continous_adjust <- rbind(Result_q_continous_adjust,res_var)
}
Result_q_continous_adjust
Result_q_continous_adjust$test <- ifelse(Result_q_continous_adjust$`Pr(>|t|)`<(0.05/11),1,0)


#################sTable: Q of PRS with associition exposure and outcome################
######IV continous
modeldata <- DATA 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
g=5
##95 SNP
modeldata$IV <- modeldata$wPRS_95
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
cof <- c("Exposure","Exposure_log","townsend_group_1","edu_new","BMI_group",
         "smoke","sleep_group","PA_group","APOE4","outcome","alcohol_liver_alcoholliver")
Result_q_continous <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV + ",paste(covar_f_name,collapse=" + "))
  
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
    fit <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  r <- summary(fit)
  res <- as.data.frame(t(r$coefficients[2,]))
  res$low <- res$Estimate-1.96*res$`Std. Error`
  res$up <- res$Estimate+1.96*res$`Std. Error`
  res$beta <- paste0(sprintf("%0.3f",res$Estimate)," (",
                     sprintf("%0.3f",res$low),", ",
                     sprintf("%0.3f",res$up),")")
  res$p <- p_out(res$`Pr(>|t|)`)
 
  res_var <- cbind(xvar, res)
  Result_q_continous <- rbind(Result_q_continous,res_var)
}
Result_q_continous

###IV group######
#################sTable: Q of PRS with associition exposure and outcome################
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
         "smoke","sleep_group","PA_group","APOE4","outcome","alcohol_liver_alcoholliver")
Result_q <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- paste("Var ~ IV_group + ",paste(covar_f_name,collapse=" + "))

  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
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
  fit1 <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
    fit1 <- glm(as.formula(fun), data=modeldata,family=binomial(link="logit"))
  }
  anoval_result <- anova(fit,fit1,test="Chisq")
  p <- p_out(anoval_result$`Pr(>Chi)`[2])
  res_var <- c(xvar, res$beta,p)
  Result_q <- rbind(Result_q,res_var)
}
Result_q
names(Result_q) <- c("VAR",paste0("Q",1:g),"P for trend")

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
         "smoke","sleep_group","PA_group","APOE4","outcome","alcohol_liver_alcoholliver")
Result_q_men <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- "Var ~ IV_group "
  for(o in covar_f_name){ fun <- paste(fun,"+",o)}
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
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
  fit1 <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
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
         "smoke","sleep_group","PA_group","APOE4","outcome","alcohol_liver_alcoholliver")
Result_q_women <- data.frame()
for(xvar in cof){
  modeldata$Var <- modeldata[,xvar]
  fun <- "Var ~ IV_group "
  for(o in covar_f_name){ fun <- paste(fun,"+",o)}
  fit <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
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
  fit1 <- glm(as.formula(fun), data=modeldata)
  if(xvar %in% c("outcome","alcohol_liver_alcoholliver")){
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
######sTable_qPRS###########
write.csv(Result_q_all,paste0(project_path,"result/stable_Qprs",Date))


#####PART Five: nonlinear MR#############################################################
#########ovearall############
DATA_m <- DATA
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "Exposure"
y_name <- "outcome"
y_time <- "years"
x_name <- "Exposure"
mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 strata_method = "residual", 
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = 10)
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   powers = 1,
                                                   method = "FE",
                                                   d=1,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   #average.exposure.associations = TRUE,
                                                   ref=min(xmean),
                                                   nboot=1000,
                                                   seed=200))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
f_data <- model$figure$data
f_data$pHR <- exp(f_data$yest)
f_data$pHR_low <- exp(f_data$lci)
f_data$pHR_up <- exp(f_data$uci)
plot_nlmr <- ggplot(f_data, aes(x = x)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  
  #geom_line(aes(y = pHR_low), color = "grey90",alpha = 0.5) +
  #geom_line(aes(y = pHR_up), color = "grey90",alpha = 0.5) +
  geom_ribbon(aes(x,ymin = pHR_low, ymax = pHR_up),
              alpha = 0.6,fill="skyblue1") +
  geom_line(aes(y = pHR), color = "skyblue4") +
  theme_bw() +
  scale_y_continuous(name="pHR for Overall",n.breaks = 10) +
  scale_x_continuous(name=expression("predicted alcohol consumption"~log[10]*"(unit/week)"),
                     expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   plot.subtitle = element_text(size=8)) +
  labs(subtitle = label_text) 
#####non liear residual MR in ovearll####
plot_nlmr
###########Figure nonlinear MR - men#######
DATA_m <- DATA %>% filter(sex==1)
covar_f_name <- c("age","area","array",paste0("component",1:20))
x_prs_name <- "wPRS_95"
ex_name <- "Exposure"
y_name <- "outcome"
y_time <- "years"
x_name <- "Exposure"
mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 strata_method = "residual", 
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = 10)
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   powers = 1,
                                                   #d=1,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean),
                                                   nboot=1000,
                                                   seed=200))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
f_data <- model$figure$data
f_data$pHR <- exp(f_data$yest)
f_data$pHR_low <- exp(f_data$lci)
f_data$pHR_up <- exp(f_data$uci)
plot_nlmr_men <- ggplot(f_data, aes(x = x)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  
  #geom_line(aes(y = pHR_low), color = "grey90",alpha = 0.5) +
  #geom_line(aes(y = pHR_up), color = "grey90",alpha = 0.5) +
  geom_ribbon(aes(x,ymin = pHR_low, ymax = pHR_up),
              alpha = 0.6,fill="skyblue1") +
  geom_line(aes(y = pHR), color = "skyblue4") +
  theme_bw() +
  scale_y_continuous(name="pHR for Men",n.breaks = 10) +
  scale_x_continuous(name=expression("predicted alcohol consumption"~log[10]*"(unit/week)"),
                     expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   plot.subtitle = element_text(size=8)) +
  labs(subtitle = label_text) 
#####non liear residual MR in men####
plot_nlmr_men
###########Figure nonlinear MR - women#######
DATA_m <- DATA %>% filter(sex==0)
covar_f_name <- c("age","area","array",paste0("component",1:20))
x_prs_name <- "wPRS_95"
y_name <- "outcome"
y_time <- "years"
x_name <- "Exposure"
mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 strata_method = "residual", 
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = 10)
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   powers = 1,
                                                   method = "DL",
                                                   d=1,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   #average.exposure.associations = TRUE,
                                                   ref=min(xmean),
                                                   nboot=1000,
                                                   seed=200))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
f_data <- model$figure$data
f_data$pHR <- exp(f_data$yest)
f_data$pHR_low <- exp(f_data$lci)
f_data$pHR_up <- exp(f_data$uci)
plot_nlmr_women <- ggplot(f_data, aes(x = x)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  
  #geom_line(aes(y = pHR_low), color = "grey90",alpha = 0.5) +
  #geom_line(aes(y = pHR_up), color = "grey90",alpha = 0.5) +
  geom_ribbon(aes(x,ymin = pHR_low, ymax = pHR_up),
              alpha = 0.6,fill="skyblue1") +
  geom_line(aes(y = pHR), color = "skyblue4") +
  theme_bw() +
  scale_y_continuous(name="pHR for Women",n.breaks = 10) +
  scale_x_continuous(name=expression("predicted alcohol consumption"~log[10]*"(unit/week)"),
                     expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   plot.subtitle = element_text(size=8)) +
  labs(subtitle = label_text) 
#####non liear residual MR in women####
plot_nlmr_women
###non linear residual MR#####
grid.arrange(plot_nlmr,plot_nlmr_men,plot_nlmr_women,ncol=3)
############non linear MR transform X-aix##########
#########ovearall############
DATA_m <- DATA
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "Exposure"
y_name <- "outcome"
y_time <- "years"
x_name <- "Exposure"
mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 strata_method = "residual", 
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = 10)
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   powers = 1,
                                                   method = "FE",
                                                   d=1,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   #average.exposure.associations = TRUE,
                                                   ref=min(xmean),
                                                   nboot=1000,
                                                   seed=200))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
f_data <- model$figure$data
f_data$pHR <- exp(f_data$yest)
f_data$pHR_low <- exp(f_data$lci)
f_data$pHR_up <- exp(f_data$uci)
f_data$Xorignal <- 10^f_data$x
plot_nlmr_1 <- ggplot(f_data, aes(x = Xorignal)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  geom_ribbon(aes(Xorignal,ymin = pHR_low, ymax = pHR_up),
              alpha = 0.6,fill="skyblue1") +
  geom_line(aes(y = pHR), color = "skyblue4") +
  theme_bw() +
  scale_y_continuous(name="pHR for Overall",n.breaks = 10) +
  scale_x_continuous(name="predicted alcohol consumption (unit/week)",
                     expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   plot.subtitle = element_text(size=8)) +
  labs(subtitle = label_text) 

###########Figure nonlinear MR - men#######
DATA_m <- DATA %>% filter(sex==1)
covar_f_name <- c("age","area","array",paste0("component",1:20))
x_prs_name <- "wPRS_95"
ex_name <- "Exposure"
y_name <- "outcome"
y_time <- "years"
x_name <- "Exposure"
mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 strata_method = "residual", 
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = 10)
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   powers = 1,
                                                   #d=1,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean),
                                                   nboot=1000,
                                                   seed=200))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
f_data <- model$figure$data
f_data$pHR <- exp(f_data$yest)
f_data$pHR_low <- exp(f_data$lci)
f_data$pHR_up <- exp(f_data$uci)
f_data$Xorignal <- 10^f_data$x
plot_nlmr_men_1 <- ggplot(f_data, aes(x = Xorignal)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  geom_ribbon(aes(Xorignal,ymin = pHR_low, ymax = pHR_up),
              alpha = 0.6,fill="skyblue1") +
  geom_line(aes(y = pHR), color = "skyblue4") +
  theme_bw() +
  scale_y_continuous(name="pHR for Men",n.breaks = 10) +
  scale_x_continuous(name="predicted alcohol consumption(unit/week)",
                     expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   plot.subtitle = element_text(size=8)) +
  labs(subtitle = label_text) 
###########Figure nonlinear MR - women#######
DATA_m <- DATA %>% filter(sex==0)
covar_f_name <- c("age","area","array",paste0("component",1:20))
x_prs_name <- "wPRS_95"
y_name <- "outcome"
y_time <- "years"
x_name <- "Exposure"
mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 strata_method = "residual", 
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = 10)
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   powers = 1,
                                                   method = "DL",
                                                   d=1,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   #average.exposure.associations = TRUE,
                                                   ref=min(xmean),
                                                   nboot=1000,
                                                   seed=200))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
f_data <- model$figure$data
f_data$pHR <- exp(f_data$yest)
f_data$pHR_low <- exp(f_data$lci)
f_data$pHR_up <- exp(f_data$uci)
f_data$Xorignal <- 10^f_data$x
plot_nlmr_women_1 <- ggplot(f_data, aes(x = Xorignal)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  geom_ribbon(aes(Xorignal,ymin = pHR_low, ymax = pHR_up),
              alpha = 0.6,fill="skyblue1") +
  geom_line(aes(y = pHR), color = "skyblue4") +
  theme_bw() +
  scale_y_continuous(name="pHR for Women",n.breaks = 10) +
  scale_x_continuous(name="predicted alcohol consumption unit/week)",
                     expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   plot.subtitle = element_text(size=8)) +
  labs(subtitle = label_text) 
###non linear residual MR transform X-aix#####
grid.arrange(plot_nlmr_1,plot_nlmr_men_1,plot_nlmr_women_1,ncol=3)


#####PART Six#############################################################
########################Table 2: linear MR analysis ##################
#########supplementary table R2 for PRS########
##all
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
modeldata$Exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
x = "95w-log10E"
R2_95wprs_log10e <- R2_F_output(modeldata,x)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
x = "95w-E"
R2_95wprs_e <- R2_F_output(modeldata,x)
R2_all <- rbind(R2_95wprs_log10e,R2_95wprs_e)
R2_all$IV_R2 <- sprintf("%0.2f",R2_all$IV_R2*100)
R2_all$IV_F <- sprintf("%0.1f",R2_all$IV_F)
R2_all$R2_value_adjsut <- sprintf("%0.2f",R2_all$R2_value_adjsut*100)
R2_all$F_value_adjsut <- sprintf("%0.1f",as.numeric(R2_all$F_value_adjsut))
R2_All <- R2_all[,c(1,3,2,4)]
##Men
modeldata <- DATA %>% filter(sex==1)
modeldata$IV <- modeldata$wPRS_95
modeldata$Exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
x = "95w-log10E"
R2_95wprs_log10e_men <- R2_F_output(modeldata,x)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
covar_f_name <- c("age","area","array",paste0("component",1:20))
x = "95w-E"
R2_95wprs_e_men <- R2_F_output(modeldata,x)
R2_all <- rbind(R2_95wprs_log10e_men,R2_95wprs_e_men)
R2_all$IV_R2 <- sprintf("%0.2f",R2_all$IV_R2*100)
R2_all$IV_F <- sprintf("%0.1f",R2_all$IV_F)
R2_all$R2_value_adjsut <- sprintf("%0.2f",R2_all$R2_value_adjsut*100)
R2_all$F_value_adjsut <- sprintf("%0.1f",as.numeric(R2_all$F_value_adjsut))
R2_men <- R2_all[,c(1,3,2,4)]
##Women
modeldata <- DATA %>% filter(sex==0)
modeldata$IV <- modeldata$wPRS_95
modeldata$Exposure <- log10(modeldata$perweek_alchol_unit_cal_new)
x = "95w-log10E"
R2_95wprs_log10e_women <- R2_F_output(modeldata,x)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
covar_f_name <- c("age","area","array",paste0("component",1:20))
x = "95w-E"
R2_95wprs_e_women <- R2_F_output(modeldata,x)
R2_all <- rbind(R2_95wprs_log10e_women,R2_95wprs_e_women)
R2_all$IV_R2 <- sprintf("%0.2f",R2_all$IV_R2*100)
R2_all$IV_F <- sprintf("%0.1f",R2_all$IV_F)
R2_all$R2_value_adjsut <- sprintf("%0.2f",R2_all$R2_value_adjsut*100)
R2_all$F_value_adjsut <- sprintf("%0.1f",as.numeric(R2_all$F_value_adjsut))
R2_women <- R2_all[,c(1,3,2,4)]
R2_sTable <- data.frame(cbind(t(R2_All),t(R2_men),t(R2_women)))
names(R2_sTable) <- c("Overall","","Men","","Women","")
##########supplement table R2 for the prs#############
write.csv(R2_sTable,paste0(project_path,"result/R2_sTable",Date))
####Overall 10-fold####
library(caret)
source(paste0(project_path,"code/00_mr_function.R"))
load(paste0(project_path,"data/drink_96snp.RData")) ##drink_96snp
drink_96snp_name <- names(drink_96snp)[2:97]
drink_95snp_name <- drink_96snp_name[!(drink_96snp_name %in% "rs13024996")]
seed = 2023
md <- DATA
md$IV <- md$wPRS_95
folds <- createFolds(md$eid,k=10)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))

s =drink_96snp_name[1]
x_snp_first <- data.frame()
y_snp_first <- data.frame()

x_snp_second <- data.frame()
y_snp_second <- data.frame()
for(i in 1:10){
  set_1 <- md[folds[[i]],] 
  set_2 <-md[-folds[[i]],] 
  ##summary exposure and snp
  fun_x <- paste("Exposure ~ var +", paste(covar_f_name,collapse = " + "))
  fun_y <- paste("Surv(years, outcome) ~ var +", paste(covar_f_name,collapse = " + "))
  
  for(s in drink_96snp_name){
    modeldata <- set_1
    modeldata$var <- modeldata[,s]
    ##x
    fit_x <- lm(as.formula(fun_x),data=modeldata)
    r <- summary(fit_x)
    x_snp <- r$coefficients[2,]
    x_snp <- c(i,s,x_snp)
    names(x_snp) <- c("folds","SNP","beta","se","t","p")
    x_snp_first <- rbind(x_snp_first,x_snp)
    ##y
    fit_y <- coxph(as.formula(fun_y), data=modeldata)
    r <- summary(fit_y)
    y_snp <- r$coefficients[1,]
    y_snp <- c(i,s,y_snp)
    names(y_snp) <- c("folds","SNP","beta","HR","se","z","p")
    y_snp_second <- rbind(y_snp_second,y_snp)
  }
  for(s in drink_96snp_name){
    modeldata <- set_2
    modeldata$var <- modeldata[,s]
    ##x
    fit_x <- lm(as.formula(fun_x),data=modeldata)
    r <- summary(fit_x)
    x_snp <- r$coefficients[2,]
    x_snp <- c(i,s,x_snp)
    names(x_snp) <- c("folds","SNP","beta","se","t","p")
    x_snp_second <- rbind(x_snp_second,x_snp)
    ##y
    fit_y <- coxph(as.formula(fun_y), data=modeldata)
    r <- summary(fit_y)
    y_snp <- r$coefficients[1,]
    y_snp <- c(i,s,y_snp)
    names(y_snp) <- c("folds","SNP","beta","HR","se","z","p")
    y_snp_first <- rbind(y_snp_first,y_snp)
  }
}

names(x_snp_first) <- c("folds","SNP","xbeta","xse","t","xp")
names(x_snp_second) <- c("folds","SNP","xbeta","xse","t","xp")
names(y_snp_first) <-  c("folds","SNP","ybeta","HR","yse","z","yp")
names(y_snp_second) <-  c("folds","SNP","ybeta","HR","yse","z","yp")
########overall################
first_dataset <- merge(x_snp_first,y_snp_first,by=c("folds","SNP"))
second_dataset <- merge(x_snp_second,y_snp_second,by=c("folds","SNP"))
first_dataset_all <- first_dataset
second_dataset_all <- first_dataset
###meta-analysis##########
res_mr_sub_first <- data.frame()
for(i in 1:10){
  d <- first_dataset %>% filter(folds==i)
  b_out <-  as.numeric(d$ybeta)
  b_exp <- as.numeric(d$xbeta)
  se_out <- as.numeric(d$yse)
  se_exp <- as.numeric(d$xse)
  snp <- d$SNP
  ##IVW
  res_ivw <- mr_ivw(b_exp, b_out, se_exp, se_out)
  ##MR-Egger
  res <- mr_egger_regression(b_exp, b_out, se_exp, se_out)
  res_egger_test <-  data.frame(b=res$b_i, se=res$se_i, pval=res$pval_i)
  res_egger_test$nsnp <- NA
  res_egger_test$Q <- NA
  res_egger_test$Q_df <- NA
  res_egger_test$Q_pval <- NA
  res_egger <- data.frame(b=res$b, se=res$se, pval=res$pval, nsnp=res$nsnp, 
                          Q = res$Q, Q_df = res$Q_df, Q_pval = res$Q_pval)
  ##weighted 
  res_weimedian <- mr_weighted_median(b_exp, b_out, se_exp, se_out)
  method <- c("IVW","MR-Egger","ME-Egger_inter","Weighted median")
  res_mr <- rbind(res_ivw,res_egger,res_egger_test,res_weimedian)
  res_mr <- cbind(i,method,res_mr)
  res_mr$hr <- exp(res_mr$b)
  res_mr$lp <- exp(res_mr$b-1.96*res_mr$se)
  res_mr$up <- exp(res_mr$b+1.96*res_mr$se)
  res_mr$HR <- paste0(sprintf("%0.2f",res_mr$hr), " (", sprintf("%0.2f",res_mr$lp),", ", sprintf("%0.2f",res_mr$up),")")
  res_mr_sub_first <- rbind(res_mr_sub_first,res_mr)
}

res_mr_sub_second <- data.frame()
for(i in 1:10){
  d <- second_dataset %>% filter(folds==i)
  b_out <-  as.numeric(d$ybeta)
  b_exp <- as.numeric(d$xbeta)
  se_out <- as.numeric(d$yse)
  se_exp <- as.numeric(d$xse)
  snp <- d$SNP
  ##IVW
  res_ivw <- mr_ivw(b_exp, b_out, se_exp, se_out)
  ##MR-Egger
  res <- mr_egger_regression(b_exp, b_out, se_exp, se_out)
  res_egger_test <-  data.frame(b=res$b_i, se=res$se_i, pval=res$pval_i)
  res_egger_test$nsnp <- NA
  res_egger_test$Q <- NA
  res_egger_test$Q_df <- NA
  res_egger_test$Q_pval <- NA
  res_egger <- data.frame(b=res$b, se=res$se, pval=res$pval, nsnp=res$nsnp, 
                          Q = res$Q, Q_df = res$Q_df, Q_pval = res$Q_pval)
  ##weighted 
  res_weimedian <- mr_weighted_median(b_exp, b_out, se_exp, se_out)
  method <- c("IVW","MR-Egger","ME-Egger_inter","Weighted median")
  res_mr <- rbind(res_ivw,res_egger,res_egger_test,res_weimedian)
  res_mr <- cbind(i,method,res_mr)
  res_mr$hr <- exp(res_mr$b)
  res_mr$lp <- exp(res_mr$b-1.96*res_mr$se)
  res_mr$up <- exp(res_mr$b+1.96*res_mr$se)
  res_mr$HR <- paste0(sprintf("%0.2f",res_mr$hr), " (", sprintf("%0.2f",res_mr$lp),", ", sprintf("%0.2f",res_mr$up),")")
  res_mr_sub_second <- rbind(res_mr_sub_second,res_mr)
}
###meta-analysis
library(metafor)
method_list <- unique(res_mr_sub_first$method)
result_meta_f_1
Meta_result <- data.frame()
for(i in 1:4){
  Methods <- method_list[i]
  method_1_meta_data <- res_mr_sub_first %>% filter(method==Methods)
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="FE")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "FE"
  result_FE <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_FE)
  
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="DL")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "DL"
  result_DL <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_DL)
  
}
Meta_result_first <- Meta_result
Meta_result <- data.frame()
for(i in 1:4){
  Methods <- method_list[i]
  method_1_meta_data <- res_mr_sub_second %>% filter(method==Methods)
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="FE")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "FE"
  result_FE <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_FE)
  
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="DL")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "DL"
  result_DL <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_DL)
  
}
Meta_result_second <- Meta_result
##combinded all
res_mr_sub_all <- rbind(res_mr_sub_first,res_mr_sub_second)
Meta_result <- data.frame()
for(i in 1:4){
  Methods <- method_list[i]
  method_1_meta_data <- res_mr_sub_all %>% filter(method==Methods)
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="FE")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "FE"
  result_FE <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_FE)
  
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="DL")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "DL"
  result_DL <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_DL)
  
}
Meta_result_all <- Meta_result
###result_out_oput###
Meta_result <- Meta_result_first
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
Meta_result_first_out <- Meta_result

Meta_result <- Meta_result_second
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
Meta_result_second_out <- Meta_result


Meta_result <- Meta_result_all
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
Meta_result_all_out <- Meta_result
######overall meta-analysis for 10-fold output######
Meta_result_first_out
Meta_result_second_out
Meta_result_all_out
#####Women 10-fold###
####10-fold####
library(caret)
source(paste0(project_path,"code/00_mr_function.R"))
load(paste0(project_path,"data/drink_96snp.RData")) ##drink_96snp
drink_96snp_name <- names(drink_96snp)[2:97]
drink_95snp_name <- drink_96snp_name[!(drink_96snp_name %in% "rs13024996")]
seed = 2023
md <- DATA %>% filter(sex==0)
md$IV <- md$wPRS_95
folds <- createFolds(md$eid,k=10)
covar_f_name <- c("age","area","array",paste0("component",1:20))

x_snp_first <- data.frame()
y_snp_first <- data.frame()
x_snp_second <- data.frame()
y_snp_second <- data.frame()
for(i in 1:10){
  set_1 <- md[folds[[i]],] 
  set_2 <-md[-folds[[i]],] 
  ##summary exposure and snp
  fun_x <- paste("Exposure ~ var +", paste(covar_f_name,collapse = " + "))
  fun_y <- paste("Surv(years, outcome) ~ var +", paste(covar_f_name,collapse = " + "))
  
  for(s in drink_96snp_name){
    modeldata <- set_1
    modeldata$var <- modeldata[,s]
    ##x
    fit_x <- lm(as.formula(fun_x),data=modeldata)
    r <- summary(fit_x)
    x_snp <- r$coefficients[2,]
    x_snp <- c(i,s,x_snp)
    names(x_snp) <- c("folds","SNP","beta","se","t","p")
    x_snp_first <- rbind(x_snp_first,x_snp)
    ##y
    fit_y <- coxph(as.formula(fun_y), data=modeldata)
    r <- summary(fit_y)
    y_snp <- r$coefficients[1,]
    y_snp <- c(i,s,y_snp)
    names(y_snp) <- c("folds","SNP","beta","HR","se","z","p")
    y_snp_second <- rbind(y_snp_second,y_snp)
  }
  for(s in drink_96snp_name){
    modeldata <- set_2
    modeldata$var <- modeldata[,s]
    ##x
    fit_x <- lm(as.formula(fun_x),data=modeldata)
    r <- summary(fit_x)
    x_snp <- r$coefficients[2,]
    x_snp <- c(i,s,x_snp)
    names(x_snp) <- c("folds","SNP","beta","se","t","p")
    x_snp_second <- rbind(x_snp_second,x_snp)
    ##y
    fit_y <- coxph(as.formula(fun_y), data=modeldata)
    r <- summary(fit_y)
    y_snp <- r$coefficients[1,]
    y_snp <- c(i,s,y_snp)
    names(y_snp) <- c("folds","SNP","beta","HR","se","z","p")
    y_snp_first <- rbind(y_snp_first,y_snp)
  }
}

names(x_snp_first) <- c("folds","SNP","xbeta","xse","t","xp")
names(x_snp_second) <- c("folds","SNP","xbeta","xse","t","xp")

names(y_snp_first) <-  c("folds","SNP","ybeta","HR","yse","z","yp")
names(y_snp_second) <-  c("folds","SNP","ybeta","HR","yse","z","yp")

########Women################
first_dataset <- merge(x_snp_first,y_snp_first,by=c("folds","SNP"))
second_dataset <- merge(x_snp_second,y_snp_second,by=c("folds","SNP"))
first_dataset_women <- first_dataset
second_dataset_women <- first_dataset
###meta-analysis##########
res_mr_sub_first <- data.frame()
for(i in 1:10){
  d <- first_dataset %>% filter(folds==i)
  b_out <-  as.numeric(d$ybeta)
  b_exp <- as.numeric(d$xbeta)
  se_out <- as.numeric(d$yse)
  se_exp <- as.numeric(d$xse)
  snp <- d$SNP
  ##IVW
  res_ivw <- mr_ivw(b_exp, b_out, se_exp, se_out)
  ##MR-Egger
  res <- mr_egger_regression(b_exp, b_out, se_exp, se_out)
  res_egger_test <-  data.frame(b=res$b_i, se=res$se_i, pval=res$pval_i)
  res_egger_test$nsnp <- NA
  res_egger_test$Q <- NA
  res_egger_test$Q_df <- NA
  res_egger_test$Q_pval <- NA
  res_egger <- data.frame(b=res$b, se=res$se, pval=res$pval, nsnp=res$nsnp, 
                          Q = res$Q, Q_df = res$Q_df, Q_pval = res$Q_pval)
  ##weighted 
  res_weimedian <- mr_weighted_median(b_exp, b_out, se_exp, se_out)
  method <- c("IVW","MR-Egger","ME-Egger_inter","Weighted median")
  res_mr <- rbind(res_ivw,res_egger,res_egger_test,res_weimedian)
  res_mr <- cbind(i,method,res_mr)
  res_mr$hr <- exp(res_mr$b)
  res_mr$lp <- exp(res_mr$b-1.96*res_mr$se)
  res_mr$up <- exp(res_mr$b+1.96*res_mr$se)
  res_mr$HR <- paste0(sprintf("%0.2f",res_mr$hr), " (", sprintf("%0.2f",res_mr$lp),", ", sprintf("%0.2f",res_mr$up),")")
  res_mr_sub_first <- rbind(res_mr_sub_first,res_mr)
}

res_mr_sub_second_women <- data.frame()
for(i in 1:10){
  d <- second_dataset %>% filter(folds==i)
  b_out <-  as.numeric(d$ybeta)
  b_exp <- as.numeric(d$xbeta)
  se_out <- as.numeric(d$yse)
  se_exp <- as.numeric(d$xse)
  snp <- d$SNP
  ##IVW
  res_ivw <- mr_ivw(b_exp, b_out, se_exp, se_out)
  ##MR-Egger
  res <- mr_egger_regression(b_exp, b_out, se_exp, se_out)
  res_egger_test <-  data.frame(b=res$b_i, se=res$se_i, pval=res$pval_i)
  res_egger_test$nsnp <- NA
  res_egger_test$Q <- NA
  res_egger_test$Q_df <- NA
  res_egger_test$Q_pval <- NA
  res_egger <- data.frame(b=res$b, se=res$se, pval=res$pval, nsnp=res$nsnp, 
                          Q = res$Q, Q_df = res$Q_df, Q_pval = res$Q_pval)
  ##weighted 
  res_weimedian <- mr_weighted_median(b_exp, b_out, se_exp, se_out)
  method <- c("IVW","MR-Egger","ME-Egger_inter","Weighted median")
  res_mr <- rbind(res_ivw,res_egger,res_egger_test,res_weimedian)
  res_mr <- cbind(i,method,res_mr)
  res_mr$hr <- exp(res_mr$b)
  res_mr$lp <- exp(res_mr$b-1.96*res_mr$se)
  res_mr$up <- exp(res_mr$b+1.96*res_mr$se)
  res_mr$HR <- paste0(sprintf("%0.2f",res_mr$hr), " (", sprintf("%0.2f",res_mr$lp),", ", sprintf("%0.2f",res_mr$up),")")
  res_mr_sub_second <- rbind(res_mr_sub_second,res_mr)
}
###meta-analysis
library(metafor)
method_list <- unique(res_mr_sub_first$method)

Meta_result <- data.frame()
for(i in 1:4){
  Methods <- method_list[i]
  method_1_meta_data <- res_mr_sub_first %>% filter(method==Methods)
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="FE")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "FE"
  result_FE <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_FE)
  
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="DL")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "DL"
  result_DL <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_DL)
  
}
Meta_result_first <- Meta_result
Meta_result <- data.frame()
for(i in 1:4){
  Methods <- method_list[i]
  method_1_meta_data <- res_mr_sub_second %>% filter(method==Methods)
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="FE")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "FE"
  result_FE <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_FE)
  
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="DL")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "DL"
  result_DL <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_DL)
  
}
Meta_result_second <- Meta_result
##combinded all
res_mr_sub_all <- rbind(res_mr_sub_first,res_mr_sub_second)
Meta_result <- data.frame()
for(i in 1:4){
  Methods <- method_list[i]
  method_1_meta_data <- res_mr_sub_all %>% filter(method==Methods)
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="FE")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "FE"
  result_FE <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_FE)
  
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="DL")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "DL"
  result_DL <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_DL)
  
}
Meta_result_all <- Meta_result
###result_out_oput###
Meta_result <- Meta_result_first
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
Meta_result_first_out_women <- Meta_result

Meta_result <- Meta_result_second
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
Meta_result_second_out_women <- Meta_result


Meta_result <- Meta_result_all
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
Meta_result_all_out_women <- Meta_result
######women meta-analysis for 10-fold output######
Meta_result_first_out_women
Meta_result_second_out_women
Meta_result_all_out_women
#####Men 10-fold###
####10-fold####
library(caret)
source(paste0(project_path,"code/00_mr_function.R"))
load(paste0(project_path,"data/drink_96snp.RData")) ##drink_96snp
drink_96snp_name <- names(drink_96snp)[2:97]
drink_95snp_name <- drink_96snp_name[!(drink_96snp_name %in% "rs13024996")]
seed = 2023
md <- DATA %>% filter(sex==1)
md$IV <- md$wPRS_95
folds <- createFolds(md$eid,k=10)
covar_f_name <- c("age","area","array",paste0("component",1:20))

x_snp_first <- data.frame()
y_snp_first <- data.frame()

x_snp_second <- data.frame()
y_snp_second <- data.frame()
for(i in 1:10){
  set_1 <- md[folds[[i]],] 
  set_2 <-md[-folds[[i]],] 
  ##summary exposure and snp
  fun_x <- paste("Exposure ~ var +", paste(covar_f_name,collapse = " + "))
  fun_y <- paste("Surv(years, outcome) ~ var +", paste(covar_f_name,collapse = " + "))
  
  for(s in drink_96snp_name){
    modeldata <- set_1
    modeldata$var <- modeldata[,s]
    ##x
    fit_x <- lm(as.formula(fun_x),data=modeldata)
    r <- summary(fit_x)
    x_snp <- r$coefficients[2,]
    x_snp <- c(i,s,x_snp)
    names(x_snp) <- c("folds","SNP","beta","se","t","p")
    x_snp_first <- rbind(x_snp_first,x_snp)
    ##y
    fit_y <- coxph(as.formula(fun_y), data=modeldata)
    r <- summary(fit_y)
    y_snp <- r$coefficients[1,]
    y_snp <- c(i,s,y_snp)
    names(y_snp) <- c("folds","SNP","beta","HR","se","z","p")
    y_snp_second <- rbind(y_snp_second,y_snp)
  }
  for(s in drink_96snp_name){
    modeldata <- set_2
    modeldata$var <- modeldata[,s]
    ##x
    fit_x <- lm(as.formula(fun_x),data=modeldata)
    r <- summary(fit_x)
    x_snp <- r$coefficients[2,]
    x_snp <- c(i,s,x_snp)
    names(x_snp) <- c("folds","SNP","beta","se","t","p")
    x_snp_second <- rbind(x_snp_second,x_snp)
    ##y
    fit_y <- coxph(as.formula(fun_y), data=modeldata)
    r <- summary(fit_y)
    y_snp <- r$coefficients[1,]
    y_snp <- c(i,s,y_snp)
    names(y_snp) <- c("folds","SNP","beta","HR","se","z","p")
    y_snp_first <- rbind(y_snp_first,y_snp)
  }
}

names(x_snp_first) <- c("folds","SNP","xbeta","xse","t","xp")
names(x_snp_second) <- c("folds","SNP","xbeta","xse","t","xp")

names(y_snp_first) <-  c("folds","SNP","ybeta","HR","yse","z","yp")
names(y_snp_second) <-  c("folds","SNP","ybeta","HR","yse","z","yp")
########Men################
first_dataset <- merge(x_snp_first,y_snp_first,by=c("folds","SNP"))
second_dataset <- merge(x_snp_second,y_snp_second,by=c("folds","SNP"))
first_dataset_men <- first_dataset
second_dataset_men <- first_dataset
###meta-analysis##########
res_mr_sub_first <- data.frame()
for(i in 1:10){
  d <- first_dataset %>% filter(folds==i)
  b_out <-  as.numeric(d$ybeta)
  b_exp <- as.numeric(d$xbeta)
  se_out <- as.numeric(d$yse)
  se_exp <- as.numeric(d$xse)
  snp <- d$SNP
  ##IVW
  res_ivw <- mr_ivw(b_exp, b_out, se_exp, se_out)
  ##MR-Egger
  res <- mr_egger_regression(b_exp, b_out, se_exp, se_out)
  res_egger_test <-  data.frame(b=res$b_i, se=res$se_i, pval=res$pval_i)
  res_egger_test$nsnp <- NA
  res_egger_test$Q <- NA
  res_egger_test$Q_df <- NA
  res_egger_test$Q_pval <- NA
  res_egger <- data.frame(b=res$b, se=res$se, pval=res$pval, nsnp=res$nsnp, 
                          Q = res$Q, Q_df = res$Q_df, Q_pval = res$Q_pval)
  ##weighted 
  res_weimedian <- mr_weighted_median(b_exp, b_out, se_exp, se_out)
  method <- c("IVW","MR-Egger","ME-Egger_inter","Weighted median")
  res_mr <- rbind(res_ivw,res_egger,res_egger_test,res_weimedian)
  res_mr <- cbind(i,method,res_mr)
  res_mr$hr <- exp(res_mr$b)
  res_mr$lp <- exp(res_mr$b-1.96*res_mr$se)
  res_mr$up <- exp(res_mr$b+1.96*res_mr$se)
  res_mr$HR <- paste0(sprintf("%0.2f",res_mr$hr), " (", sprintf("%0.2f",res_mr$lp),", ", sprintf("%0.2f",res_mr$up),")")
  res_mr_sub_first <- rbind(res_mr_sub_first,res_mr)
}

res_mr_sub_second_women <- data.frame()
for(i in 1:10){
  d <- second_dataset %>% filter(folds==i)
  b_out <-  as.numeric(d$ybeta)
  b_exp <- as.numeric(d$xbeta)
  se_out <- as.numeric(d$yse)
  se_exp <- as.numeric(d$xse)
  snp <- d$SNP
  ##IVW
  res_ivw <- mr_ivw(b_exp, b_out, se_exp, se_out)
  ##MR-Egger
  res <- mr_egger_regression(b_exp, b_out, se_exp, se_out)
  res_egger_test <-  data.frame(b=res$b_i, se=res$se_i, pval=res$pval_i)
  res_egger_test$nsnp <- NA
  res_egger_test$Q <- NA
  res_egger_test$Q_df <- NA
  res_egger_test$Q_pval <- NA
  res_egger <- data.frame(b=res$b, se=res$se, pval=res$pval, nsnp=res$nsnp, 
                          Q = res$Q, Q_df = res$Q_df, Q_pval = res$Q_pval)
  ##weighted 
  res_weimedian <- mr_weighted_median(b_exp, b_out, se_exp, se_out)
  method <- c("IVW","MR-Egger","ME-Egger_inter","Weighted median")
  res_mr <- rbind(res_ivw,res_egger,res_egger_test,res_weimedian)
  res_mr <- cbind(i,method,res_mr)
  res_mr$hr <- exp(res_mr$b)
  res_mr$lp <- exp(res_mr$b-1.96*res_mr$se)
  res_mr$up <- exp(res_mr$b+1.96*res_mr$se)
  res_mr$HR <- paste0(sprintf("%0.2f",res_mr$hr), " (", sprintf("%0.2f",res_mr$lp),", ", sprintf("%0.2f",res_mr$up),")")
  res_mr_sub_second <- rbind(res_mr_sub_second,res_mr)
}
###meta-analysis
library(metafor)
method_list <- unique(res_mr_sub_first$method)

Meta_result <- data.frame()
for(i in 1:4){
  Methods <- method_list[i]
  method_1_meta_data <- res_mr_sub_first %>% filter(method==Methods)
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="FE")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "FE"
  result_FE <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_FE)
  
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="DL")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "DL"
  result_DL <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_DL)
  
}
Meta_result_first <- Meta_result
Meta_result <- data.frame()
for(i in 1:4){
  Methods <- method_list[i]
  method_1_meta_data <- res_mr_sub_second %>% filter(method==Methods)
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="FE")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "FE"
  result_FE <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_FE)
  
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="DL")
  estimate <- meta_result$b
  se <- meta_result$se
  pval <- meta_result$pval
  ci.lb <- meta_result$ci.lb
  ci.ub <- meta_result$ci.ub
  Q <- meta_result$QE
  hete_p <- meta_result$QEp
  m <- "DL"
  result_DL <- data.frame(Methods,m,estimate,se, pval, ci.lb, ci.ub, Q, hete_p)
  Meta_result <- rbind(Meta_result,result_DL)
  
}
Meta_result_second <- Meta_result
##combinded all
res_mr_sub_all <- rbind(res_mr_sub_first,res_mr_sub_second)
Meta_result <- data.frame()
for(i in 1:4){
  Methods <- method_list[i]
  method_1_meta_data <- res_mr_sub_all %>% filter(method==Methods)
  meta_result <- rma(b,sei=se,data=method_1_meta_data, method="FE")
  estimate <- meta_result$