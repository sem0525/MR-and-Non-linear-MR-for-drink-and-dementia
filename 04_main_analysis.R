##############Analysis and sub into sex group#############
###Analysis###
project_path <- "/home/zoo/project/P1_drinkMR_dementia/"
##result save date
Date <- paste0("_drinker_20231114.csv") ##313958

###load the data
#load(paste0(project_path,"data/data_all_nondrink_93585.RData"))
load(paste0(project_path,"data/data_all_analysis_313958.RData"))
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
###########################################
DATA <- data_all_analysis
##target exposure
DATA$Exposure <- DATA$perweek_alchol_unit_cal_new
##seperate the alcohol consumption 
DATA$alcohol_2group <- cut(DATA$perweek_alchol_unit_cal_new,
                           breaks=c(0,14,max(DATA$perweek_alchol_unit_cal_new)),
                           include.lowest = TRUE)
table(DATA$alcohol_2group)
##seperate the alcohol consumption 
alcohol_quintile <- quantile(DATA$perweek_alchol_unit_cal_new,probs=seq(0,1,0.2))
DATA$alcohol_5group <- cut(DATA$perweek_alchol_unit_cal_new,
                           breaks=alcohol_quintile,
                           include.lowest = TRUE)
table(DATA$alcohol_5group)
##q10
alcohol_quintile <- quantile(DATA$perweek_alchol_unit_cal_new,probs=seq(0,1,0.1))
DATA$alcohol_10group <- cut(DATA$perweek_alchol_unit_cal_new,
                            breaks=alcohol_quintile,
                            include.lowest = TRUE)
table(DATA$alcohol_10group)
##
DATA$alcohol_per7_10g <- cut(DATA$perweek_alchol_unit_cal_new,
                           breaks=c(seq(0,63,7),max(DATA$perweek_alchol_unit_cal_new)),
                           include.lowest = TRUE)
table(DATA$alcohol_per7_10g)

DATA$alcohol_per7_20g <- cut(DATA$perweek_alchol_unit_cal_new,
                             breaks=c(seq(0,139,7),max(DATA$perweek_alchol_unit_cal_new)),
                             include.lowest = TRUE)
table(DATA$alcohol_per7_20g)

#####PART One#############################################################
########################Table 1: Descritpitive ##################
d <- DATA
conf <-  c("years","perweek_alchol_unit_cal_new","alcohol_2group",
           "age_group","sex","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","APOE4_status","ACD"
) 
d$group <- as.numeric(d$alcohol_2group)
re <- table_1_group(conf,d,FALSE)
#####table1_all######
write.csv(re,paste0(project_path,"result/Table1_all",Date))
#subgroup to sex
conf <-  c("years","perweek_alchol_unit_cal_new","alcohol_2group",
           "age_group","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","APOE4_status","ACD") 
##men
d_men <- d %>% filter(sex==1)
re_men <- table_1_group(conf,d_men,FALSE)
##women
d_women <- d %>% filter(sex==0)
re_women <- table_1_group(conf,d_women,FALSE)
re_all <- cbind(re_men,re_women)
#####table1_sugbroup######
write.csv(re_all,paste0(project_path,"result/Table1_sub",Date))


#####PART Two#############################################################
###Figure table 2: the table of incident across alcohol ##################
modeldata <- DATA
modeldata$E_group <- modeldata$alcohol_10group
table2_all <- table_2_incidence(modeldata)
Overall <- rep(NA,10)
table2_all <- rbind(Overall,table2_all)
###man
modeldata <- DATA %>% filter(sex==1)
alcohol_quintile <- quantile(modeldata$perweek_alchol_unit_cal_new,probs=seq(0,1,0.1))
modeldata$E_group <- cut(modeldata$perweek_alchol_unit_cal_new,
                            breaks=alcohol_quintile,
                            include.lowest = TRUE)
table2_men <- table_2_incidence(modeldata)
Men <- rep(NA,10)
table2_men <- rbind(Men,table2_men)
###women
modeldata <- DATA %>% filter(sex==0)
alcohol_quintile <- quantile(modeldata$perweek_alchol_unit_cal_new,probs=seq(0,1,0.1))
modeldata$E_group <- cut(modeldata$perweek_alchol_unit_cal_new,
                         breaks=alcohol_quintile,
                         include.lowest = TRUE)
table2_women <- table_2_incidence(modeldata)
Women <- rep(NA,10)
table2_women <- rbind(Women,table2_women)
inc_d_case <- rbind(table2_all,table2_men,table2_women)
######figure table 1.2#####
write.csv(inc_d_case,paste0(project_path,"result/fig1table",Date))
########################Figure 2: HR Epidemiological Association ##################
modeldata <- DATA
plotx <- "Alcohol consumption (unit/week)"
population <- "Overall"
covar_f_name <- c("age","sex","APOE4_status","area","townsend_group_1","edu_new")
#,"BMI_group",
#behaviour,"smoke","sleep_group","PA_group"
rcs_fig_all <- rcs_fig(modeldata,covar_f_name,population,plotx,6)
rcs_fig_all <- rcs_fig_all  + 
  scale_x_continuous(name=plotx,expand = c(0,0),limits=c(0,250),n.breaks=10) + 
  scale_y_continuous(name=population,expand = c(0,0),limits=c(0,11),n.breaks=5) +
  scale_fill_jama()

##men
modeldata <- DATA %>% filter(sex==1)
plotx <- "Alcohol consumption (unit/week)"
population <- "Men"
covar_f_name <- c("age","APOE4_status","area","townsend_group_1","edu_new")
rcs_fig_men <- rcs_fig(modeldata,covar_f_name,population,plotx,5)
rcs_fig_men <- rcs_fig_men + 
  scale_x_continuous(name=plotx,expand = c(0,0),limits=c(0,250),n.breaks=10) + 
  scale_y_continuous(name=population,expand = c(0,0),limits=c(0,11),n.breaks=5) +
  scale_fill_jama()

##women
modeldata <- DATA %>% filter(sex==0)
plotx <- "Alcohol consumption (unit/week)"
population <- "Women"
covar_f_name <- c("age","APOE4_status","area","townsend_group_1","edu_new")
rcs_fig_women <- rcs_fig(modeldata,covar_f_name,population,plotx,5)

rcs_fig_women <- rcs_fig_women + 
  scale_x_continuous(name=plotx,expand = c(0,0),limits=c(0,150),n.breaks=10) + 
  scale_y_continuous(name=population,expand = c(0,0),limits=c(0,11),n.breaks=5) +
  scale_fill_jama()
######figure#####
grid.arrange(rcs_fig_all,rcs_fig_men,rcs_fig_women,ncol=3)


#####PART Three#############################################################
######figure: PRS associated with exposure######
modeldata <- DATA 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
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
res[1,] <- c(1,0,1,1)
res$g <- factor(paste0("Q",1:g),levels=paste0("Q",1:g))
res$group <- 1:10
res_g <- modeldata %>% group_by(group) %>% mutate(mean_PRS=mean(IV,na.rm=T)) %>% select(group,mean_PRS)
res_g <- unique(res_g)
res <- merge(res,res_g,by="group")
g1 <- ggplot(res,aes(x=mean_PRS,y=Estimate  ))+
  geom_hline(yintercept = 1.0,
             #linetype='dashed',
             color="grey",alpha = 0.2)+
  geom_vline(xintercept = res_g$mean_PRS,linetype='dotted',color="grey")+
  geom_point(position=position_dodge(0)) +
  geom_errorbar(aes(ymin=low,ymax=up),
                position=position_dodge(0),width=0.001)+
  #labs(subtitle = label_text) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   legend.position="bottom")+
  scale_y_continuous(name="Coefficient of Alcohol consumption (unit/week)",
                     limits=c(0,6),
                     n.breaks = 10) +
  xlab("Decile of PRS") +
  scale_color_lancet() 


g2 <- ggplot(modeldata, aes(x=IV)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=10,), 
        axis.title=element_text(size=10))+
  geom_histogram(binwidth=0.02,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.8) + 
  scale_x_continuous(name="Alcohol Consumption-PRS",expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="Overall",expand = c(0,0),limit=c(0,45000),n.breaks=5) +
  scale_fill_jama()

plot_PRS_q10_all <- ggplot2.two_y_axis(g2,g1)
grid.arrange(plot_PRS_q10_all)

######figure -men: PRS associated with exposure######
modeldata <- DATA %>% filter(sex==1)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
modeldata$IV <- modeldata$wPRS_95
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
modeldata$IV_group <- cut(modeldata$IV,include.lowest = T,
                          breaks=quantile(modeldata$IV,probs= seq(0,1,0.1)))
g=10
table(modeldata$group)
modeldata$group <- as.numeric(modeldata$IV_group)
modeldata$Var <- modeldata$Exposure
fun <- paste("Var ~ IV_group +", paste(covar_f_name,collapse = " + "))
fit <- glm(as.formula(fun), data=modeldata)
r <- summary(fit)
res <- as.data.frame(r$coefficients[1:g,1:2])
res$low <- res$Estimate-1.96*res$`Std. Error`
res$up <- res$Estimate+1.96*res$`Std. Error`
res[1,] <- c(1,0,1,1)
res$g <- factor(paste0("Q",1:g),levels=paste0("Q",1:g))
res$group <- 1:10
res_g <- modeldata %>% group_by(group) %>% mutate(mean_PRS=mean(IV,na.rm=T)) %>% select(group,mean_PRS)
res_g <- unique(res_g)
res <- merge(res,res_g,by="group")
g1 <- ggplot(res,aes(x=mean_PRS,y=Estimate  ))+
  geom_hline(yintercept = 1.0,
             #linetype='dashed',
             color="grey",alpha = 0.2)+
  geom_vline(xintercept = res_g$mean_PRS,linetype='dotted',color="grey")+
  geom_point(position=position_dodge(0)) +
  geom_errorbar(aes(ymin=low,ymax=up),
                position=position_dodge(0),width=0.001)+
  #labs(subtitle = label_text) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   legend.position="bottom")+
  scale_y_continuous(name="Coefficient of Alcohol consumption (unit/week)",
                     limits=c(0,8.5),
                     n.breaks = 10) +
  xlab("Decile of PRS") +
  scale_color_lancet() 


g2 <- ggplot(modeldata, aes(x=IV)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=10,), 
        axis.title=element_text(size=10))+
  geom_histogram(binwidth=0.02,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.8) + 
  scale_x_continuous(name="Alcohol Consumption-PRS",expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="Men",expand = c(0,0),limit=c(0,20000),n.breaks=) +
  scale_fill_jama()

plot_PRS_q10_men <- ggplot2.two_y_axis(g2,g1)
grid.arrange(plot_PRS_q10_men)


######figure -women: PRS associated with exposure######
modeldata <- DATA %>% filter(sex==0)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log10(modeldata$perweek_alchol_unit_cal_new)
modeldata$IV <- modeldata$wPRS_95
covar_f_name <- c("age","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
modeldata$IV_group <- cut(modeldata$IV,include.lowest = T,
                          breaks=quantile(modeldata$IV,probs= seq(0,1,0.1)))
g=10
table(modeldata$group)
modeldata$group <- as.numeric(modeldata$IV_group)
modeldata$Var <- modeldata$Exposure
fun <- paste("Var ~ IV_group +", paste(covar_f_name,collapse = " + "))
fit <- glm(as.formula(fun), data=modeldata)
r <- summary(fit)
res <- as.data.frame(r$coefficients[1:g,1:2])
res$low <- res$Estimate-1.96*res$`Std. Error`
res$up <- res$Estimate+1.96*res$`Std. Error`
res[1,] <- c(1,0,1,1)
res$g <- factor(paste0("Q",1:g),levels=paste0("Q",1:g))
res$group <- 1:10
res_g <- modeldata %>% group_by(group) %>% mutate(mean_PRS=mean(IV,na.rm=T)) %>% select(group,mean_PRS)
res_g <- unique(res_g)
res <- merge(res,res_g,by="group")
g1 <- ggplot(res,aes(x=mean_PRS,y=Estimate  ))+
  geom_hline(yintercept = 1.0,
             #linetype='dashed',
             color="grey",alpha = 0.2)+
  geom_vline(xintercept = res_g$mean_PRS,linetype='dotted',color="grey")+
  geom_point(position=position_dodge(0)) +
  geom_errorbar(aes(ymin=low,ymax=up),
                position=position_dodge(0),width=0.001)+
  #labs(subtitle = label_text) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   legend.position="bottom")+
  scale_y_continuous(name="Coefficient of Alcohol consumption (unit/week)",
                     limits=c(0,5),
                     n.breaks = 10) +
  xlab("Decile of PRS") +
  scale_color_lancet() 


g2 <- ggplot(modeldata, aes(x=IV)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=10,), 
        axis.title=element_text(size=10))+
  geom_histogram(binwidth=0.02,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.8) + 
  scale_x_continuous(name="Alcohol Consumption-PRS",expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="Women",expand = c(0,0),limit=c(0,20000),n.breaks=) +
  scale_fill_jama()

plot_PRS_q10_men <- ggplot2.two_y_axis(g2,g1)
grid.arrange(plot_PRS_q10_men)

#################sTable: Q of PRS with associition exposure and outcome################
modeldata <- DATA 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log(modeldata$perweek_alchol_unit_cal_new)

covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
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
  Result_q <- rbind(Result_q,res_var)
}

Result_q
names(Result_q) <- c("VAR",paste0("Q",1:g),"P for trend")

##########men################
modeldata <- DATA %>% filter(sex==1)
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
modeldata$Exposure_log <- log(modeldata$perweek_alchol_unit_cal_new+0.00001)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
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
modeldata$Exposure_log <- log(modeldata$perweek_alchol_unit_cal_new+0.00001)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
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

#################figure: Q of PRS with assocition exposure ################
modeldata <- DATA 
modeldata$Exposure <- modeldata$perweek_alchol_unit_cal_new
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
##95 SNP
modeldata$IV <- modeldata$wPRS_95
modeldata$IV_group <- cut(modeldata$IV,
                          breaks=quantile(modeldata$IV,probs= seq(0,1,0.1)))
modeldata$group <- as.numeric(modeldata$IV_group)
g = dim(table(modeldata$IV_group))
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
cof <- "Exposure"
modeldata$Var <- modeldata$Exposure
fun <- "Var ~ IV_group "
for(o in covar_f_name){ fun <- paste(fun,"+",o)}
fit <- glm(as.formula(fun), data=modeldata)
r <- summary(fit)
res <- as.data.frame(r$coefficients[1:g,1:2])
res$low <- res$Estimate-1.96*res$`Std. Error`
res$up <- res$Estimate+1.96*res$`Std. Error`
res[1,] <- c(1,1,1,1)
names(res) <- c("beta","std","low","up")
res$x  <- factor(c(paste0("Q",1:g)),levels=c(paste0("Q",1:g)))

modeldata$IV_group <- as.numeric(modeldata$IV_group)
fit1 <- glm(as.formula(fun), data=modeldata)

ptest <- anova(fit,fit1,test="Chisq")
p_trend <- ptest$`Pr(>Chi)`[2]
p_trend <- p_out_fig(p_trend)
###figure 1.1
label_text <- substitute(expr=paste("P for trend test: ",italic("p "),p_trend1),
                         env = base::list(p_trend1=p_trend))
fig2 <- ggplot(res,aes(x,beta))+
  geom_hline(yintercept = 1.0,linetype='dashed',color="grey")+
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(x,ymin=low,ymax=up),
                position=position_dodge(0.5),width=0.1)+
  labs(subtitle = label_text) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   legend.position="bottom")+
  scale_y_continuous(name="Coefficient of Alcohol consumption (unit/week)",n.breaks = 10) +
  scale_x_discrete(name="Decile of PRS") +
  scale_color_lancet() 
######figure2#####
fig2
##q=10 nonlinaer MR###
###########Figure nonlinear MR#######
DATA_m <- DATA
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
#DATA_m$Exposure <- log10(DATA_m$Exposure )
x_name <- "Exposure"
###q=10
mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
#ym <- mdata[,y_name]
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
q=10
####rank########
summ_data_rank <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 strata_method = "ranked", 
                                 #strata_method = "residual",
                                 strata_bound=c(0.025,0.025,0.975,0.975),
                                 #report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data_rank$strata_statistics
summ_data_rank$Heterogeneity_results

Fdata_sum <- summ_data_rank$summary
Fdata_sum$strata <- factor(1:10,levels=c(1:10))
nlmr_figure_rank_10 <- ggplot()+
  geom_point(data=Fdata_sum,aes(strata,xmean),
             size=0.5,alpha = 1,colour="skyblue4") +
  geom_errorbar(data=Fdata_sum,aes(strata,ymin = xmin, ymax = xmax),
                width=0.2) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  scale_y_continuous(name="Alcohol consumption (unit/week)",n.breaks = 10) +
  scale_x_discrete(name=paste(q,"Strata by doubly-ranked stratified method")) +
  scale_color_lancet()
########orignial data for rank######
nlmr_figure_rank_10
################predicted############
nlmr_figure_rank_10_1 <- ggplot()+
  geom_point(data=Fdata_sum,aes(strata,bx),
             size=0.5,alpha = 1,colour="skyblue4") +
  geom_errorbar(data=Fdata_sum,aes(strata,ymin = bx-1.96*bxse, ymax = bx+1.96*bxse),
                width=0.2) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  
  scale_y_continuous(name="Alcohol consumption (unit/week)",n.breaks = 10) +
  scale_x_discrete(name=paste(q,"Strata by doubly-ranked stratified method")) +
  scale_color_lancet()
########predicted data for rank######
nlmr_figure_rank_10_1
########res########
summ_data_res <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 strata_method = "residual", 
                                 strata_bound=c(0.025,0.025,0.975,0.975),
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)

Fdata_sum <- summ_data_res$summary
Fdata_sum$strata <- factor(1:10,levels=c(1:10))
nlmr_figure_res_10 <- ggplot()+
  geom_point(data=Fdata_sum,aes(strata,xmean),
             size=0.5,alpha = 1,colour="skyblue4") +
  geom_errorbar(data=Fdata_sum,aes(strata,ymin = xmin, ymax = xmax),
                width=0.2) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  scale_y_continuous(name="Alcohol consumption (unit/week)",n.breaks = 10) +
  scale_x_discrete(name=paste(q,"Strata by residual stratified method")) +
  scale_color_lancet()
########orignial data for residual######
nlmr_figure_res_10
################predicted################
nlmr_figure_res_10_1 <- ggplot()+
  geom_point(data=Fdata_sum,aes(strata,bx),
             size=0.5,alpha = 1,colour="skyblue4") +
  geom_errorbar(data=Fdata_sum,aes(strata,ymin = bx-1.96*bxse, ymax = bx+1.96*bxse),
                width=0.2) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  
  scale_y_continuous(name="Alcohol consumption (unit/week)",n.breaks = 10) +
  scale_x_discrete(name=paste(q,"Strata by residual stratified method")) +
  scale_color_lancet()
########predicted data for residual######
nlmr_figure_res_10_1

grid.arrange(nlmr_figure_rank_10,
             nlmr_figure_res_10,
             nlmr_figure_rank_10_1,
             nlmr_figure_res_10_1,nrow=2)


##########log transform #######
xm_log <- log10(mdata[,x_name])
##rank
summ_data_rank_log <- create_nlmr_summary(y = ym,
                                      x = xm_log,
                                      g = x_prsm,
                                      covar = covar_fm,
                                      family = "coxph",
                                      strata_method = "ranked", 
                                      #strata_method = "residual",
                                      strata_bound=c(0.025,0.025,0.975,0.975),
                                      #report_GR=TRUE,
                                      report_het=TRUE,
                                      extra_statistics=TRUE,
                                      q = q)
Fdata_sum <- summ_data_rank_log$summary
Fdata_sum$strata <- factor(1:10,levels=c(1:10))
nlmr_figure_rank_10log <- ggplot()+
  geom_point(data=Fdata_sum,aes(strata,xmean),
             size=0.5,alpha = 1,colour="skyblue4") +
  geom_errorbar(data=Fdata_sum,aes(strata,ymin = xmin, ymax = xmax),
                width=0.2) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  scale_y_continuous(name="Alcohol consumption (log(unit/week))",n.breaks = 10) +
  scale_x_discrete(name=paste(q,"Strata by doubly-ranked stratified method")) +
  scale_color_lancet()
########orignial data for rank######
nlmr_figure_rank_10log
############predicted############
nlmr_figure_rank_10log_1 <- ggplot()+
  geom_point(data=Fdata_sum,aes(strata,bx),
             size=0.5,alpha = 1,colour="skyblue4") +
  geom_errorbar(data=Fdata_sum,aes(strata,ymin = bx-1.96*bxse, ymax = bx+1.96*bxse),
                width=0.2) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  
  scale_y_continuous(name="Alcohol consumption (log(unit/week))",n.breaks = 10) +
  scale_x_discrete(name=paste(q,"Strata by doubly-ranked stratified method")) +
  scale_color_lancet()
########predicted data for rank######
nlmr_figure_rank_10log_1
########res########
summ_data_res_log <- create_nlmr_summary(y = ym,
                                     x = xm_log,
                                     g = x_prsm,
                                     covar = covar_fm,
                                     family = "coxph",
                                     strata_method = "residual", 
                                     strata_bound=c(0.025,0.025,0.975,0.975),
                                     report_het=TRUE,
                                     extra_statistics=TRUE,
                                     q = q)

Fdata_sum <- summ_data_res_log$summary
Fdata_sum$strata <- factor(1:10,levels=c(1:10))
nlmr_figure_res_10log <- ggplot()+
  geom_point(data=Fdata_sum,aes(strata,xmean),
             size=0.5,alpha = 1,colour="skyblue4") +
  geom_errorbar(data=Fdata_sum,aes(strata,ymin = xmin, ymax = xmax),
                width=0.2) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  scale_y_continuous(name="Alcohol consumption (log(unit/week))",n.breaks = 10) +
  scale_x_discrete(name=paste(q,"Strata by residual stratified method")) +
  scale_color_lancet()
########orignial data for residual######
nlmr_figure_res_10log
########predicted########
nlmr_figure_res_10log_1 <- ggplot()+
  geom_point(data=Fdata_sum,aes(strata,bx),
             size=0.5,alpha = 1,colour="skyblue4") +
  geom_errorbar(data=Fdata_sum,aes(strata,ymin = bx-1.96*bxse, ymax = bx+1.96*bxse),
                width=0.2) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  
  scale_y_continuous(name="Alcohol consumption (log(unit/week))",n.breaks = 10) +
  scale_x_discrete(name=paste(q,"Strata by residual stratified method")) +
  scale_color_lancet()
########predicted data for residual######
nlmr_figure_res_10log_1

grid.arrange(nlmr_figure_rank_10log,
             nlmr_figure_res_10log,
             nlmr_figure_rank_10log_1,
             nlmr_figure_res_10log_1,nrow=2)


###########Figure nonlinear MR#######
DATA_m <- DATA
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
q=100
#DATA_m$Exposure <- log10(DATA_m$Exposure )
x_name <- "Exposure"

mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
#ym <- mdata[,y_name]
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 #family="binomial",
                                 #controlsonly=TRUE,
                                 strata_method = "ranked", 
                                 #strata_method = "residual", 
                                 report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data$summary
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   #family="binomial",
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean)))
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
f_data <- model$figure$data

label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
plot_nlmr <- model$figure +
  scale_y_continuous(name="pHR for overall",n.breaks = 6) +
  scale_x_continuous(name="predicted alcohol consumption (unit/week)",expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
  labs(subtitle = label_text)
plot_nlmr
###########Figure nonlinear MR - men#######
DATA_m <- DATA %>% filter(sex==1)
covar_f_name <- c("age","area","array",paste0("component",1:20))
x_prs_name <- "wPRS_95"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
q=100
#DATA_m$Exposure <- log10(DATA_m$Exposure )
x_name <- "Exposure"
mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
#ym <- mdata[,y_name]
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 #family="binomial",
                                 #controlsonly=TRUE,
                                 strata_method = "ranked", 
                                 #strata_method = "residual", 
                                 report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data$summary
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   #family="binomial",
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean)))
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
f_data <- model$figure$data

label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
plot_nlmr_men <- model$figure +
  scale_y_continuous(name="pHR for men",n.breaks = 6) +
  scale_x_continuous(name="predicted alcohol consumption (unit/week)",expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
  labs(subtitle = label_text)
plot_nlmr_men


###########Figure nonlinear MR - women#######
DATA_m <- DATA %>% filter(sex==0)
covar_f_name <- c("age","area","array",paste0("component",1:20))
x_prs_name <- "wPRS_95"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
q=100
#DATA_m$Exposure <- log10(DATA_m$Exposure )
x_name <- "Exposure"
mdata <- DATA_m[,c(y_time,x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- Surv(mdata[,y_time],mdata[,y_name])
#ym <- mdata[,y_name]
xm <- mdata[,x_name]
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 #family="binomial",
                                 #controlsonly=TRUE,
                                 strata_method = "ranked", 
                                 #strata_method = "residual", 
                                 report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data$summary
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   #family="binomial",
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean)))
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
f_data <- model$figure$data

label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
plot_nlmr_women <- model$figure +
  scale_y_continuous(name="pHR for women",n.breaks = 6) +
  scale_x_continuous(name="predicted alcohol consumption (unit/week)",expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
  labs(subtitle = label_text)
plot_nlmr_women



###########log10: Figure nonlinear MR############
###########Figure nonlinear MR -#######
DATA_m <- DATA
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
q=10
DATA_m$Exposure <- log10(DATA_m$Exposure )
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
                                 #controlsonly=TRUE,
                                 strata_method = "ranked", 
                                 #strata_method = "residual", 
                                 #report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data$summary
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean)))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
f_data <- model$figure$data

label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
plot_nlmr_log_rank <- model$figure +
  scale_y_continuous(name="pHR for overall",n.breaks = 6) +
  scale_x_continuous(name="predicted alcohol consumption (log(unit/week))",expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
  labs(subtitle = label_text)
#######plot_nlmr_log_rank#######
plot_nlmr_log_rank
############################

summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 #controlsonly=TRUE,
                                 #strata_method = "ranked", 
                                 strata_method = "residual", 
                                 #report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data$summary
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean)))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
f_data <- model$figure$data

label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
plot_nlmr_log_res <- model$figure +
  scale_y_continuous(name="pHR for overall",n.breaks = 6) +
  scale_x_continuous(name="predicted alcohol consumption (log(unit/week))",expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
  labs(subtitle = label_text)
#######plot_nlmr_log_res#######
plot_nlmr_log_res
###########log Figure nonlinear MR -men#######
DATA_m <- DATA %>% filter(sex==1)
covar_f_name <- c("age","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
q=10
DATA_m$Exposure <- log10(DATA_m$Exposure )
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
                                 #controlsonly=TRUE,
                                 strata_method = "ranked", 
                                 #strata_method = "residual", 
                                 #report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data$summary
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean)))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
f_data <- model$figure$data

label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
plot_nlmr_log_rank_men <- model$figure +
  scale_y_continuous(name="pHR for men",n.breaks = 6) +
  scale_x_continuous(name="predicted alcohol consumption (log(unit/week))",expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
  labs(subtitle = label_text)
#######plot_nlmr_log_rank_men#######
plot_nlmr_log_rank_men
####################
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 #controlsonly=TRUE,
                                 #strata_method = "ranked", 
                                 strata_method = "residual", 
                                 #report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data$summary
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean)))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
f_data <- model$figure$data

label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
plot_nlmr_log_res_men <- model$figure +
  scale_y_continuous(name="pHR for men",n.breaks = 6) +
  scale_x_continuous(name="predicted alcohol consumption (log(unit/week))",expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
  labs(subtitle = label_text)
#######plot_nlmr_log_res_men#######
plot_nlmr_log_res_men
###########log Figure nonlinear MR -women#######
DATA_m <- DATA %>% filter(sex==0)
covar_f_name <- c("age","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
q=10
DATA_m$Exposure <- log10(DATA_m$Exposure )
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
                                 #controlsonly=TRUE,
                                 strata_method = "ranked", 
                                 #strata_method = "residual", 
                                 #report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data$summary
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean)))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
f_data <- model$figure$data

label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
plot_nlmr_log_rank_women <- model$figure +
  scale_y_continuous(name="pHR for women",n.breaks = 6) +
  scale_x_continuous(name="predicted alcohol consumption (log(unit/week))",expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
  labs(subtitle = label_text)
#######plot_nlmr_log_rank_women
plot_nlmr_log_rank_women
############################
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "coxph",
                                 #controlsonly=TRUE,
                                 #strata_method = "ranked", 
                                 strata_method = "residual", 
                                 #report_GR=TRUE,
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = q)
summ_data$summary
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   family="coxph",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean)))
summary(model)
p_allover <- model$coefficients[5]
p_quad <- model$p_tests[2]
p_Co <- model$p_tests[4]
#p value
p_quad0 <- p_out_fig(p_quad)
p_allover0 <- p_out_fig(p_allover)
p_Co <- p_out_fig(p_Co)
f_data <- model$figure$data

label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_quad1,
                                    ", Cochran Q ",italic("P"), p_Co1),
                         env = base::list(p_allover1=p_allover0,
                                          p_quad1=p_quad0,
                                          p_Co1=p_Co
                         ))
plot_nlmr_log_res_women <- model$figure +
  scale_y_continuous(name="pHR for women",n.breaks = 6) +
  scale_x_continuous(name="predicted alcohol consumption (log(unit/week))",expand = c(0,0), n.breaks = 10)   +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10)) +
  labs(subtitle = label_text)
###############plot_nlmr_log_res_women
plot_nlmr_log_res_women

###########old non-linear MR############
#################fig2 plot non-linear MR################
modeldata <- DATA 
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
##for strata for 10,100
sn <- c(10,100)
sum_DRnlmr_all <- list()
nlmr_model1_DRnlmr_all <- list()
DRnlmr_boxF_all <- list()
sum_nlmr_all <- list()
nlmr_model_nlmr_all <- list()
nlmr_boxF_all <- list()
a = 1
ts <- c(12,4)
for(n in sn){
  DRnlmr_strata <- sumNLMR_rank_HR(modeldata,n)
  sum_DRnlmr_all[[a]] <- DRnlmr_strata[[1]]
  nlmr_model1_DRnlmr_all[[a]] <- DRnlmr_strata[[2]]
  DRnlmr_boxF_all[[a]] <- nmlr_boxfigure(DRnlmr_strata[[2]],n,ts[a])
  
  nlmr_strata <-sumNLMR_res_HR(modeldata,n)
  sum_nlmr_all[[a]] <- nlmr_strata[[1]]
  nlmr_model_nlmr_all[[a]] <- nlmr_strata[[2]]
  nlmr_boxF_all[[a]] <- nmlr_boxfigure(nlmr_strata[[2]],n,ts[a])
  a = a+1
}

DRnlmr_fig_all <- list()
nlmr_fig_all <- list()
a=1
for(a in 1:2){
  DRnlmr_fig_all[[a]] <- nmlr_boxfigure_summary(sum_DRnlmr_all[[a]]$summary,sn[a],ts[a])
  nlmr_fig_all[[a]] <- nmlr_boxfigure_summary(sum_nlmr_all[[a]]$summary,sn[a],ts[a])
  a = a+ 1
}
grid.arrange(DRnlmr_fig_all[[1]],DRnlmr_boxF_all[[1]],nrow=2)
DRnlmr_fig_all[[2]] <- nmlr_boxfigure_summary(sum_DRnlmr_all[[2]]$summary,sn[2],6)
grid.arrange(DRnlmr_fig_all[[2]],DRnlmr_boxF_all[[2]],nrow=2)

#################Men fig2 plot non-linear MR################
modeldata <- DATA %>% filter(sex==1)
covar_f_name <- c("age","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
##for strata for 10,100
sn <- c(10,100)
sum_DRnlmr_men <- list()
nlmr_model1_DRnlmr_men <- list()
DRnlmr_boxF_men <- list()
sum_nlmr_men <- list()
nlmr_model_nlmr_men <- list()
nlmr_boxF_men <- list()
a = 1
ts <- c(12,4)
for(n in sn){
  DRnlmr_strata <- sumNLMR_rank_HR(modeldata,n)
  sum_DRnlmr_men[[a]] <- DRnlmr_strata[[1]]
  nlmr_model1_DRnlmr_men[[a]] <- DRnlmr_strata[[2]]
  DRnlmr_boxF_men[[a]] <- nmlr_boxfigure(DRnlmr_strata[[2]],n,ts[a])
  
  nlmr_strata <-sumNLMR_res_HR(modeldata,n)
  sum_nlmr_men[[a]] <- nlmr_strata[[1]]
  nlmr_model_nlmr_men[[a]] <- nlmr_strata[[2]]
  nlmr_boxF_men[[a]] <- nmlr_boxfigure(nlmr_strata[[2]],n,ts[a])
  a = a+1
}

DRnlmr_fig_men <- list()
nlmr_fig_men <- list()
a=1
for(a in 1:2){
  DRnlmr_fig_men[[a]] <- nmlr_boxfigure_summary(sum_DRnlmr_men[[a]]$summary,sn[a],ts[a])
  nlmr_fig_men[[a]] <- nmlr_boxfigure_summary(sum_nlmr_men[[a]]$summary,sn[a],ts[a])
  a = a+ 1
}
grid.arrange(DRnlmr_fig_men[[1]],DRnlmr_boxF_men[[1]],nrow=2)
DRnlmr_fig_men[[2]] <- nmlr_boxfigure_summary(sum_DRnlmr_men[[2]]$summary,sn[2],6)
grid.arrange(DRnlmr_fig_men[[2]],DRnlmr_boxF_men[[2]],nrow=2)

#################Women fig2 plot non-linear MR################
modeldata <- DATA %>% filter(sex==0)
covar_f_name <- c("age","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
##for strata for 10,100
sn <- c(10,100)
sum_DRnlmr_women <- list()
nlmr_model1_DRnlmr_women <- list()
DRnlmr_boxF_women <- list()
sum_nlmr_women <- list()
nlmr_model_nlmr_women <- list()
nlmr_boxF_women <- list()
a = 1
ts <- c(12,4)
for(n in sn){
  DRnlmr_strata <- sumNLMR_rank_HR(modeldata,n)
  sum_DRnlmr_women[[a]] <- DRnlmr_strata[[1]]
  nlmr_model1_DRnlmr_women[[a]] <- DRnlmr_strata[[2]]
  DRnlmr_boxF_women[[a]] <- nmlr_boxfigure(DRnlmr_strata[[2]],n,ts[a])
  
  nlmr_strata <-sumNLMR_res_HR(modeldata,n)
  sum_nlmr_women[[a]] <- nlmr_strata[[1]]
  nlmr_model_nlmr_women[[a]] <- nlmr_strata[[2]]
  nlmr_boxF_women[[a]] <- nmlr_boxfigure(nlmr_strata[[2]],n,ts[a])
  a = a+1
}

DRnlmr_fig_women <- list()
nlmr_fig_women <- list()
a=1
for(a in 1:2){
  DRnlmr_fig_women[[a]] <- nmlr_boxfigure_summary(sum_DRnlmr_women[[a]]$summary,sn[a],ts[a])
  nlmr_fig_women[[a]] <- nmlr_boxfigure_summary(sum_nlmr_women[[a]]$summary,sn[a],ts[a])
  a = a+ 1
}
grid.arrange(DRnlmr_fig_women[[1]],DRnlmr_boxF_women[[1]],nrow=2)
DRnlmr_fig_women[[2]] <- nmlr_boxfigure_summary(sum_DRnlmr_women[[2]]$summary,sn[2],6)
grid.arrange(DRnlmr_fig_women[[2]],DRnlmr_boxF_women[[2]],nrow=2)







#####PART Four#############################################################
########################Table 2: linear MR analysis ##################
##all
modeldata <- DATA
a = 2 #exposure: unit alcohol per week
population <- "Overall"
subgroup <- "all"
IV_name <- "weighted 95 SNP"
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
MR_linear_all <- table2_MR_all_linear(modeldata,covar_f_name,population,subgroup,IV_name,a)
MR_linear_all <- table2_MR_linear(modeldata,covar_f_name,population,subgroup,IV_name,a)
MR_linear_res <- table3_res_MR(modeldata,covar_f_name,population,subgroup,IV_name,a)

##R2
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
x = "95w-E"
R2_95wprs_e <- R2_F_output(modeldata,x)

modeldata$Exposure <- log(modeldata$Exposure)
x = "95w-logE"
R2_95wprs_loge <- R2_F_output(modeldata,x)

##Men
modeldata <- DATA %>% filter(sex==1)
modeldata$IV <- modeldata$wPRS_95
a = 2 #exposure: unit alcohol per week
population <- "Men"
subgroup <- "Men"
IV_name <- "weighted 95 SNP"
covar_f_name <- c("age","area","array",paste0("component",1:20))
MR_linear_men <- table2_MR_linear(modeldata,covar_f_name,population,subgroup,IV_name,a)
x = "95w-E-men"
R2_95wprs_e_men <- R2_F_output(modeldata,x)
modeldata$Exposure <- log10(modeldata$Exposure)
x = "95w-logE-men"
R2_95wprs_loge_men <- R2_F_output(modeldata,x)


##Women
ccc_women <- as.data.frame(table(DATA$Exposure[DATA$sex==0]))
modeldata <- DATA %>% filter(sex==0)
modeldata$IV <- modeldata$wPRS_95
a = 2 #exposure: unit alcohol per week
population <- "Women"
subgroup <- "Women"
IV_name <- "weighted 95 SNP"
covar_f_name <- c("age","area","array",paste0("component",1:20))
MR_linear_women <- table2_MR_linear(modeldata,covar_f_name,population,subgroup,IV_name,a)
MR_res_women <- table3_res_MR(modeldata,covar_f_name,population,subgroup,IV_name,a)

x = "95w-E-women"
R2_95wprs_e_women <- R2_F_output(modeldata,x)

modeldata$Exposure <- log10(modeldata$Exposure)
x = "95w-logE-women"
R2_95wprs_loge_women <- R2_F_output(modeldata,x)


####Sensitivity analysis: subgroup analysis########
a = 2 #exposure: unit alcohol per week
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
IV_name <- "weighted 95 SNP"
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
r_95w <- MR_analysis_subgroup_HR(modeldata,IV_name,covar_f_name,a)
####Sensitivity analysis: subgroup analysis log transform########
a = 4 #log exposure: unit alcohol per week
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
IV_name <- "weighted 95 SNP"
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
r_95w_log <- MR_analysis_subgroup_HR(modeldata,IV_name,covar_f_name,a)
r_all <- cbind(r_95w,r_95w_log)
######sTable sensitivity subgroup MR###########
write.csv(r_all,paste0(project_path,"result/stable_MR_sub",Date))
####10-fold####
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

folds[[1]]
i=1
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

###meta: x_snp_first and y_snp_firts



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
#####women 10-fold###
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
######overall meta-analysis for 10-fold output######
Meta_result_first_out_women
Meta_result_second_out_women
Meta_result_all_out_women

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
#####men 10-fold###
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
Meta_result_first_out_men <- Meta_result

Meta_result <- Meta_result_second
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
Meta_result_second_out_men <- Meta_result


Meta_result <- Meta_result_all
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
Meta_result_all_out_men <- Meta_result
######Men meta-analysis for 10-fold output######
Meta_result_first_out_men
Meta_result_second_out_men
Meta_result_all_out_men


#####PART Six#############################################################
#################################two source sample MR:HR#######
###two sample MR analysis for 95 SNP ---HR
source(paste0(project_path,"code/00_mr_function.R"))
load(paste0(project_path,"data/drink_96snp.RData")) ##drink_96snp
drink_96snp_name <- names(drink_96snp)[2:97]
drink_95snp_name <- drink_96snp_name[!(drink_96snp_name %in% "rs13024996")]
###overall
df <- DATA
###HR
asso_disease_snp <- data.frame()
for(s in drink_95snp_name){
  modeldata <- df
  modeldata$var <- modeldata[,s]
  covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
  fun_out <- paste("Surv(years, outcome) ~ var +", paste(covar_f_name,collapse=" + "))
  fit_out <- coxph(as.formula(fun_out), data=modeldata)
  r <- summary(fit_out)
  res <- r$coefficients[1,]
  res0 <- c(s,res)
  asso_disease_snp <- rbind(asso_disease_snp,res0)
}
names(asso_disease_snp) <- c("SNP","Beta","HR","SE","z","Pvalue")
asso_disease_snp_HR <- asso_disease_snp
##harmanazed the data
dy <- asso_disease_snp %>% select(SNP,Beta,SE)
names(dy)[2:3] <- c("beta.outcome","se.outcome")
dy$beta.outcome <- as.numeric(dy$beta.outcome)
dy$se.outcome <- as.numeric(dy$se.outcome)
dy$outcome <- "dementia"
dy$id.outcome <- dy$outcome
load("/home/ukb/zoo_project/UKB/data_clean/snp_info.RData")
dx_all <- snp_info %>% filter(SNP %in% drink_95snp_name)
dx_all$Beta[dx_all$EA != dx_all$refA] <- dx_all$Beta[dx_all$EA != dx_all$refA]*(-1)

dx <- dx_all %>% select(SNP,Beta,SE)
names(dx)[2:3] <- c("beta.exposure","se.exposure")
dx$exposure <- "drink_all"
dx$id.exposure <- dx$exposure
d_mr <- merge(dx,dy,by="SNP")
d_mr$mr_keep <- TRUE
##two sample MR analysis
d <- d_mr
b_out <- d$beta.outcome
b_exp <- d$beta.exposure
se_out <- d$se.outcome
se_exp <- d$se.exposure
snp <- d$SNP
mr_wald_ratio(snp, b_exp, b_out, se_exp, se_out) 

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
res_mr
res_mr_sen_HR <- res_mr %>% select(method,HR,pval)
res_mr_sen_HR$pval <- round(res_mr_sen_HR$pval,2)

res_egger_test$HR <- paste0(sprintf("%0.4f",res_egger_test$b)," (",
                            sprintf("%0.4f",res_egger_test$se),")")

res_egger_test$pval <- p_out(res_egger_test$pval)
res_egger_test$method <- "  intercept: beta (se)"
res_egger_result <- res_egger_test %>% select(method,HR,pval)

twoSMR_result <- rbind(res_mr_sen_HR[1:2,], res_egger_result, res_mr_sen_HR[3,])

###leave one out forest
leaveone_mr <- mr_leaveoneout(d)
leaveone_mr_figdata <- leaveone_mr %>% select(SNP,b,se,p)
sfig_twosample_MR1 <- mr_leaveoneout_fig(leaveone_mr_figdata)

##scater
ggplot(data=d, aes(x=beta.exposure, y=beta.outcome)) +
  geom_errorbar(aes(ymin=beta.outcome-se.outcome, 
                    ymax=beta.outcome+se.outcome),colour="grey",width=0) +
  geom_errorbarh(aes(xmin=beta.exposure-se.exposure, 
                     xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  geom_point(size=0.1) +
  geom_abline(aes(intercept=0,slope=res_mr$b[res_mr$method=="IVW"]),color="#a6cee3") +
  labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
  theme(legend.position="top", legend.direction="vertical") +
  theme_bw()
d
d1 <- d %>% filter(abs(beta.outcome)<1) 
sfig_twosample_MR2 <- ggplot(data=d1, aes(x=beta.exposure, y=beta.outcome)) +
  geom_errorbar(aes(ymin=beta.outcome-se.outcome, 
                    ymax=beta.outcome+se.outcome),colour="grey",width=0) +
  geom_errorbarh(aes(xmin=beta.exposure-se.exposure, 
                     xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  geom_point(size=0.1) +
  geom_abline(aes(intercept=0,slope=res_mr$b[res_mr$method=="IVW"]),color="#a6cee3") +
  labs(colour="MR Test", x="SNP ffect on Alcohol Consumption", y="SNP effect on Demenita") +
  theme(legend.position="top", legend.direction="vertical") +
  theme_bw()                                               

#####two sample figure####
sfig_twosample_MR1 
sfig_twosample_MR2 



