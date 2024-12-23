###Analysis###
project_path <- "/home/zoo/project/P1_drinkMR_dementia/"
##result save date
Date <- paste0("_drinker_20231031.csv") ##313958
##############Analysis and sub into sex group#############
###Analysis###
project_path <- "/home/zoo/project/P1_drinkMR_dementia/"
##result save date
Date <- "_drinker_20231114.csv"

###load the data
load(paste0(project_path,"data/data_alldrink_407543.RData")) ##407983
data <- data_alldrink
#######All alcohol consumption is zero
data_all_nondrink <- data %>% filter(!(data$alcohol_drinker_status==2 & 
                                         data$perweek_alchol_unit_cal_new>0))
table(data_all_nondrink$alcohol_drinker_status)

data_all_analysis <- data %>% filter(data$alcohol_drinker_status==2,
                                     data$perweek_alchol_unit_cal_new>0)

###load the packages
library(dplyr)
library(survival)

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

#####PART One#############################################################
########################Table 1: Descritpitive ##################
d <- DATA
conf <-  c("years","perweek_alchol_unit_cal_new","alcohol_2group",
           "age_group","sex","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","APOE4_status"
           ) 
re <- table_1_incidence(conf,d,FALSE)
#####table1_all######
write.csv(re,paste0(project_path,"result/Table1_all",Date))
###########subgrup normal alcohol group###########
DATA$group <- as.numeric(DATA$alcohol_2group)
conf <-  c("years", "perweek_alchol_unit_cal_new","age_group","sex",
           "townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","APOE4_status")
d <- DATA %>% filter(group==1)
re_1 <- table_1_incidence(conf,d,FALSE)
#exced alcohol group
d <- DATA %>% filter(group==2)
re_2 <- table_1_incidence(conf,d,FALSE)
re_all <- cbind(re_1,re_2)
#####table1_sugbroup######
write.csv(re_all,paste0(project_path,"result/Table1_sub",Date))

#####PART Two#############################################################
########################Figure 1: HR Epidemiological Association ##################
###Figure 1.1 the relationship between exposure and outcome######
modeldata <- DATA
modeldata$E_group <- modeldata$alcohol_10group
table(modeldata$E_group)
ig = dim(table(modeldata$E_group))
##model 1
covar_f_name <- c("age","sex","area","APOE4_status")
fun <- "Surv(years, outcome) ~ E_group +"
for(o in covar_f_name){fun <- paste(fun,"+",o)}
fit <- coxph(as.formula(fun), data=modeldata)
r <- summary(fit)
res_fig <- data.frame(r$coefficients[1:ig,-4])
res_fig[ig,] <- c(0,1,0,0)
names(res_fig) <- c("beta","HR","se","p")
res_fig$Model <- "Model 1"
res_fig$x  <- factor(c(paste0("Q",2:ig),"Q1"),levels=c(paste0("Q",1:ig)))
##p for trend
modeldata$E_group_g <- as.numeric(modeldata$E_group)
table(modeldata$E_group_g,modeldata$E_group)
fun <- "Surv(years, outcome) ~ E_group_g +"
for(o in covar_f_name){fun <- paste(fun,"+",o)}
fit1 <- coxph(as.formula(fun), data=modeldata)
anoval_result <- anova(fit,fit1,test="Chisq")
p_trend1 <- p_out(anoval_result$`Pr(>|Chi|)`[2])

##model 2
covar_f_name <- c("age","sex","area","APOE4_status",#"ethnicity",
                  "townsend_group_1","edu_new","BMI_group",
                  "smoke","sleep_group","PA_group")

fun <- "Surv(years, outcome) ~ E_group +"
for(o in covar_f_name){fun <- paste(fun,"+",o)}
fit <- coxph(as.formula(fun), data=modeldata)
r <- summary(fit)
res_fig1 <- data.frame(r$coefficients[1:ig,-4])
res_fig1[ig,] <- c(0,1,0,0)
names(res_fig1) <- c("beta","HR","se","p")
res_fig1$Model <- "Model 2"
res_fig1$x  <- factor(c(paste0("Q",2:ig),"Q1"),levels=c(paste0("Q",1:ig)))
##p for trend
fun <- "Surv(years, outcome) ~ E_group_g +"
for(o in covar_f_name){fun <- paste(fun,"+",o)}
fit1 <- coxph(as.formula(fun), data=modeldata)
anoval_result <- anova(fit,fit1,test="Chisq")
p_trend2 <- p_out(anoval_result$`Pr(>|Chi|)`[2])
##combinde the data
fig_data <- rbind(res_fig,res_fig1)
table(fig_data$x)
###figure 1.1
label_text <- substitute(expr=paste("P for trend test: Model 1: ",italic("p"),p_trend10,";   Model 2: ",
                                    italic("p="),p_trend20),
                         env = base::list(p_trend10=p_trend1,
                                          p_trend20=p_trend2))
fig1.1 <- ggplot(fig_data,aes(x,HR,group=Model,color=Model))+
  geom_hline(yintercept = 1.0,linetype='dashed',color="grey")+
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(x,ymin=exp(beta-1.96*se),ymax=exp(beta+1.96*se),group=Model,color=Model),
                position=position_dodge(0.5),width=0.1)+
  labs(subtitle = label_text) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   legend.position="bottom")+
  scale_y_continuous(name="HR of Dementia",n.breaks = 6, limits=c(0.7,1.3)) +
  scale_x_discrete(name="Alcohol consumption (unit/week)") +
  scale_color_lancet() 
######figure1.1#####
fig1.1
###Figure table 1.2: the table of incident across alcohol ##################
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
inc_d_case$pyear <- paste0(sprintf("%0.2f",inc_d_case$RR)," (",
                           sprintf("%0.2f",inc_d_case$lower),", ",
                           sprintf("%0.2f",inc_d_case$upper),")")
inc_d_case$group <- inc_d_case$E_group
inc_d_case <- inc_d_case[,c("group","case","pyear")]
inc_d_case <- inc_d_case[order(inc_d_case$group),]
inc_d_case <- t(inc_d_case)
######figure table 1.2#####
write.csv(inc_d_case,paste0(project_path,"result/fig1table",Date))
###eFigure: rcs plot for HR in observatioanal  ##################
modeldata <- DATA
plotx <- "Alcohol consumption (unit/week)"
covar_f_name <- c("age","sex","APOE4_status","area","smoke", 
                  "townsend_group_1","edu_new")
#,"BMI_group","smoke","sleep_group","PA_group"
coxplot_data <- modeldata %>% select("years","outcome","Exposure",covar_f_name)
coxplot_data <- na.omit(coxplot_data)
coxplot_data$x <- coxplot_data$Exposure
g <- ggplot(coxplot_data, aes(x=Exposure)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=12,), 
        axis.title=element_text(size=12))+
  geom_histogram(binwidth=10,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.4) + 
  scale_x_continuous(name=plotx,expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="Participants",expand = c(0,0),n.breaks=10) +
  scale_fill_jama()
#Cox model
c_d1 <- class.ind(coxplot_data$edu_new)
c_d1 <- c_d1[,-dim(c_d1)[2]]

c_d2 <- class.ind(coxplot_data$area)
c_d2 <- c_d2[,-dim(c_d2)[2]]

c_d3 <- class.ind(coxplot_data$townsend_group_1)
c_d3 <- c_d3[,-dim(c_d3)[2]]

c_d4 <- class.ind(coxplot_data$smoke)
c_d4 <- c_d4[,-dim(c_d4)[2]]

ns <- 6
coxplot_data <- cbind(coxplot_data[,1:ns],c_d1, c_d2, c_d3, c_d4)
names(coxplot_data)[(ns+1):(ns+2)] <- paste0("edu",1:2)
names(coxplot_data)[(ns+5):(ns+6)] <- paste0("townsend",1:2)
names(coxplot_data)[(ns+7):(ns+8)] <- paste0("smoke",1:2)

fun <- "Surv(years,outcome) ~ rcs(Exposure,3)"
fun1 <- "Surv(years,outcome) ~ Exposure"
for(o in c(4:(ns+8))){ #29
  fun <- paste(fun,"+",names(coxplot_data)[o])
  fun1 <- paste(fun1,"+",names(coxplot_data)[o])
}
dd <- datadist(coxplot_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=coxplot_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=coxplot_data)
dd$limits$Exposure[2] <- 0
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
p_allover <- p_test[1,3]

p_allover <- p_out(p_allover)
p_nonliner <- p_out(p_nonliner)

HR <- Predict(m0,Exposure, fun=exp,ref.zero=TRUE)
#plot2
label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_nonliner1),
                         env = base::list(p_allover1=p_allover,
                                          p_nonliner1=p_nonliner))
#x_text <- max(HR$Exposure)-(max(HR$Exposure)-min(HR$Exposure))*0.3
g0 <- ggplot() +
  geom_line(data=HR, aes(Exposure,yhat),
            linetype="solid",size=1,alpha = 0.3,colour="red3") +
  geom_ribbon(data=HR, aes(Exposure,ymin = lower, ymax = upper),
              alpha = 0.2,fill="red3") +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=12),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  #annotate("text", x=x_text, y=max(HR$upper), label=label_text, cex = 4) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") + 
  scale_x_continuous(name=plotx,expand = c(0,0), n.breaks = 5)  +
  scale_y_continuous(name="HR of Dementia", n.breaks = 10 )
G <- g+labs(subtitle = label_text)
plot_rcs_HR <- ggplot2.two_y_axis(G,g0)
######figure1.3#####
grid.arrange(plot_rcs_HR)

#####PART Three#############################################################
#################sTable: Q of PRS with associition exposure and outcome################
modeldata <- DATA 
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
######sTable_qPRS###########
write.csv(Result_q,paste0(project_path,"result/stable_Qprs",Date))
#################fig2 plot non-linear MR################
#load data
modeldata <- DATA 
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
plotx <- "genetic prediction of alcohol consumption (unit/week)"
y_time <- "years"
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Dementia"

qn <- c(10,100,250,500)
RD_nlMR <- list()
i=1
for(qq in qn){
  plot_x_name <- paste0(plot_x," (",qq," strata",")")
  f <-  sumNLMR_figure_rank_HR(modeldata,puplation,plot_x_name,plot_y,1,qq)
  RD_nlMR[[i]] <- f
  i = i + 1
}
##residual 
nlMR <- list()
i=1
for(qq in qn){
  plot_x_name <- paste0(plot_x," (",qq," strata",")")
  f <-  sumNLMR_figure_res_HR(modeldata,puplation,plot_x_name,plot_y,1,qq)
  nlMR[[i]] <- f
  i = i + 1
}
######figure2.1 doubleR nonMR#####
grid.arrange(RD_nlMR[[1]],RD_nlMR[[2]],RD_nlMR[[3]],RD_nlMR[[4]],nrow=2)
######figure2.2 nonMR#####
grid.arrange(nlMR[[1]],nlMR[[2]],nlMR[[3]],nlMR[[4]],nrow=2)
########box figure########
modeldata <- DATA 
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
x_prs_name <- "wPRS_95"
DRnlmr_10strata <- sumNLMR_rank_HR(modeldata,10)
sum_DRnlmr <- DRnlmr_10strata[[1]]
nlmr_model1_DRnlmr <- DRnlmr_10strata[[2]]

nlmr_10strata <-sumNLMR_res_HR(modeldata,10)
sum_nlmr <- nlmr_10strata[[1]]
nlmr_model_nlmr <- nlmr_10strata[[2]]

fdata_sum <- nlmr_model1_DRnlmr$figure$data
nmlr_boxfigure <- function(fdata_sum){
  q10_cut <- quantile(fdata_sum$x,probs=seq(0,1,0.1))
  fdata_sum$x_g <- cut(fdata_sum$x,breaks=q10_cut,include.lowest = T)
  table(fdata_sum$x_g)
  fdata_sum$hr <- exp(fdata_sum$yest)
  nlmr_figure <- ggplot(fdata_sum,aes(x=x_g,y=hr))+
    geom_boxplot() + 
    stat_summary(fun.y = mean, geom = "point", shape = 23, size=1) +
    
    geom_hline(yintercept =1.0,linetype='dashed',color="grey")+
    theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                     axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                     legend.position="bottom")+
    scale_y_continuous(name="pHR of Alcohol liver diseases",n.breaks = 10) +
    scale_x_discrete(name="predicted Alcohol consumption (unit/week)") +
    scale_color_lancet()
  nlmr_figure
}

DRnlmr_boxF <- nmlr_boxfigure(nlmr_model1_DRnlmr$figure$data)
DRnlmr_boxF
nlmr_boxF <- nmlr_boxfigure(nlmr_model_nlmr$figure$data)
nlmr_boxF
#####PART Four#############################################################
########################Table 2: MR and group residual ##################
#wprs_95snp
a = 2 #exposure: unit alcohol per week
modeldata <- DATA
population <- "Drinkers"
subgroup <- "None"
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
###weight PRS
modeldata$IV <- modeldata$wPRS_95
IV_name <- "weighted 95 SNP"
re_95w <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup,covar_f_name,a)
re_95w
###weight PRS
modeldata$IV <- modeldata$PRS_95
IV_name <- "95 SNP"
re_95 <- PRS_mr_one_HR(modeldata,IV_name,population,subgroup,covar_f_name,a)
re_95
#combinded the result
MR_linear_result<- rbind(re_95w,re_95)
######Table 2 linear MR###########
write.csv(MR_linear_result,paste0(project_path,"result/table_MR_linear_",Date))

table2 <- MR_linear_result[,c("IV_name","HR","p_value")]

##################compete risk event##############
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
modeldata$death <- ifelse(is.na(modeldata$death_date)==FALSE,1,0)
modeldata$outcome_new <- modeldata$outcome
modeldata$outcome_new[modeldata$death==1] <- 2
modeldata$area <- ifelse(modeldata$area=="England",1,
                         ifelse(modeldata$area=="Scotland",2,3))
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
x="weight 95 PRS"
res_compet_w95 <- IV_HR_function_compete(modeldata,covar_f_name,x)

##############residual in to 2 group ######################
a = 2 #exposure: unit alcohol per week
population <- "Total"
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
IV_name <- "weight 95 SNP"
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
result_residual <- MR_analysis_residual_HR(modeldata, IV_name,covar_f_name,population,a,g=2)
######Table 2 linear MR residual 2g###########
write.csv(result_residual,paste0(project_path,"result/table_MR_linear_residual2g",Date))
####Sensitivity analysis: subgroup analysis########
a = 2 #exposure: unit alcohol per week
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
IV_name <- "weighted 95 SNP"
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
r_95w <- MR_analysis_subgroup_HR(modeldata,IV_name,covar_f_name,a)
######sTable sensitivity subgroup MR###########
write.csv(r_95w,paste0(project_path,"result/stable_MR_subg",Date))
####Sensitivity analysis: subgroup analysis lgo transform########
a = 4 #log exposure: unit alcohol per week
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
IV_name <- "weighted 95 SNP"
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
r_95w_log <- MR_analysis_subgroup_HR(modeldata,IV_name,covar_f_name,a)
######sTable sensitivity subgroup MR log ###########
write.csv(r_95w_log,paste0(project_path,"result/stable_MR_subg_log",Date))

#####PART Five#############################################################
################################2-fold two sample MR with HR#######
source(paste0(project_path,"code/00_mr_function.R"))
###two sample MR analysis for 95 SNP
df <- DATA
load(paste0(project_path,"data/drink_96snp.RData")) ##drink_96snp
drink_96snp_name <- names(drink_96snp)[2:97]
drink_95snp_name <- drink_96snp_name[!(drink_96snp_name %in% "rs13024996")]

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
result_1 <- result_list
Re1_all <- data.frame()
Re2_all <- data.frame()
for(nx in 1:100){
  result_1 <- twofold_sMR(df,drink_95snp_name)
  re1 <- result_1[[1]]
  re1$n <- nx
  Re1_all <- rbind(Re1_all,re1)
  re2 <- result_1[[2]]
  re2$n <- nx
  Re2_all <- rbind(Re2_all,re2)
}

Re2_all00 <- Re2_all 
Re1_all

summary(Re1_all$hr[method=="IVW"])
summary(Re1_all$hr[method=="MR-Egger"])
summary(Re1_all$hr[method=="Weighted median"])

mean(Re1_all$hr[method=="IVW"]) - 1.96*(sd(Re1_all$hr[method=="IVW"])/100)
mean(Re1_all$hr[method=="MR-Egger"]) - 1.96*(sd(Re1_all$hr[method=="MR-Egger"])/100)
mean(Re1_all$hr[method=="Weighted median"]) - 1.96*(sd(Re1_all$hr[method=="Weighted median"])/100)


quantile(Re1_all$hr[method=="IVW"],probs=c(0.025,0.5,0.0975,1))
quantile(Re1_all$hr[method=="MR-Egger"],probs=c(0.025,0.5,0.0975,1))
quantile(Re1_all$hr[method=="Weighted median"],probs=c(0.025,0.5,0.0975,1))

####################residual for 2-fold########
#exposure: unit alcohol per week
modeldata <- DATA
modeldata$IV <- modeldata$wPRS_95
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))

fun <- "Exposure ~ IV "
for(o in covar_f_name){fun <- paste(fun,"+",o)}
fit <- lm(as.formula(fun), data=modeldata)
r <- summary(fit)
res <- r$coefficients[2,1]
modeldata$predict <- modeldata$IV * res
modeldata$res <- modeldata$Exposure - modeldata$predict
modeldata$res_group <- cut(modeldata$res, 
                        breaks=c(min(modeldata$res)-0.01,14,max(modeldata$res)))
table(modeldata$res_group)
modeldata$res_group <- as.numeric(modeldata$res_group)
###95 SNP
##res1
df <- modeldata %>% filter(res_group==1)
res1_mr_95snp <- towf_sMR(df,drink_95snp_name_wb)
##res2
df <- modeldata %>% filter(res_group==2)
res2_mr_95snp <- towf_sMR(df,drink_95snp_name_wb)


  !!!!



#####PART Six#############################################################
#################################two source sample MR:HR#######
###two sample MR analysis for 95 SNP ---HR
source(paste0(project_path,"code/00_mr_function.R"))
load(paste0(project_path,"data/drink_96snp.RData")) ##drink_96snp
drink_96snp_name <- names(drink_96snp)[2:97]
drink_95snp_name <- drink_96snp_name[!(drink_96snp_name %in% "rs13024996")]

df <- DATA
asso_disease_snp <- data.frame()
for(s in drink_95snp_name){
  modeldata <- df
  modeldata$var <- modeldata[,s]
  covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
  fun_out <- "Surv(years, outcome) ~ var "
  for(o in covar_f_name){ fun_out <- paste(fun_out,"+",o)}
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

res_mr_sen_HR <- res_mr

###leave one out forest
leaveone_mr <- mr_leaveoneout(d)
sfig_twosample_MR1 <- mr_leaveoneout_plot(leaveone_mr)

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
  labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
  theme(legend.position="top", legend.direction="vertical") +
  theme_bw()                                               

#####two sample figure####
sfig_twosample_MR1 
sfig_twosample_MR2 

#####PART Seven#############################################################
##########Positive control analysis########
###eFigure 1.1 the relationship between exposure and outcome######
modeldata <- DATA %>% filter(disease_before_alcoholliver==0)
modeldata$E_group <- modeldata$alcohol_10group
table(modeldata$E_group)
ig = dim(table(modeldata$E_group))
##alcohol livers
modeldata$years <- modeldata$years_alcoholliver
modeldata$outcome <- modeldata$alcohol_liver_alcoholliver
##model 1
covar_f_name <- c("age","sex","area","APOE4_status")
fun <- "Surv(years, outcome) ~ E_group +"
for(o in covar_f_name){fun <- paste(fun,"+",o)}
fit <- coxph(as.formula(fun), data=modeldata)
r <- summary(fit)
res_fig <- data.frame(r$coefficients[1:ig,-4])
res_fig[ig,] <- c(0,1,0,0)
names(res_fig) <- c("beta","HR","se","p")
res_fig$Model <- "Model 1"
res_fig$x  <- factor(c(paste0("Q",2:ig),"Q1"),levels=c(paste0("Q",1:ig)))
##p for trend
modeldata$E_group_g <- as.numeric(modeldata$E_group)
table(modeldata$E_group_g,modeldata$E_group)
fun <- "Surv(years, outcome) ~ E_group_g +"
for(o in covar_f_name){fun <- paste(fun,"+",o)}
fit1 <- coxph(as.formula(fun), data=modeldata)
anoval_result <- anova(fit,fit1,test="Chisq")
p_trend1 <- p_out(anoval_result$`Pr(>|Chi|)`[2])

##model 2
covar_f_name <- c("age","sex","area","APOE4_status",#"ethnicity",
                  "townsend_group_1","edu_new","BMI_group",
                  "smoke","sleep_group","PA_group")

fun <- "Surv(years, outcome) ~ E_group +"
for(o in covar_f_name){fun <- paste(fun,"+",o)}
fit <- coxph(as.formula(fun), data=modeldata)
r <- summary(fit)
res_fig1 <- data.frame(r$coefficients[1:ig,-4])
res_fig1[ig,] <- c(0,1,0,0)
names(res_fig1) <- c("beta","HR","se","p")
res_fig1$Model <- "Model 2"
res_fig1$x  <- factor(c(paste0("Q",2:ig),"Q1"),levels=c(paste0("Q",1:ig)))
##p for trend
fun <- "Surv(years, outcome) ~ E_group_g +"
for(o in covar_f_name){fun <- paste(fun,"+",o)}
fit1 <- coxph(as.formula(fun), data=modeldata)
anoval_result <- anova(fit,fit1,test="Chisq")
p_trend2 <- p_out(anoval_result$`Pr(>|Chi|)`[2])
##combinde the data
fig_data <- rbind(res_fig,res_fig1)
table(fig_data$x)
###figure 1.1
label_text <- substitute(expr=paste("P for trend test: Model 1: ",italic("p"),p_trend10,";   Model 2: ",
                                    italic("p="),p_trend20),
                         env = base::list(p_trend10=p_trend1,
                                          p_trend20=p_trend2))
fig1.1_liver <- ggplot(fig_data,aes(x,HR,group=Model,color=Model))+
  geom_hline(yintercept = 1.0,linetype='dashed',color="grey")+
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(x,ymin=exp(beta-1.96*se),ymax=exp(beta+1.96*se),group=Model,color=Model),
                position=position_dodge(0.5),width=0.1)+
  labs(subtitle = label_text) +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=10), axis.title=element_text(size=10),
                                   legend.position="bottom")+
  scale_y_continuous(name="HR of Alcohol liver diseases",n.breaks = 6, limits=c(0,60)) +
  scale_x_discrete(name="Alcohol consumption (unit/week)") +
  scale_color_lancet() 
######figure1.1#####
fig1.1_liver
###Figure table 1.2: the table of incident across alcohol ##################
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
inc_d_case$pyear <- paste0(sprintf("%0.2f",inc_d_case$RR)," (",
                           sprintf("%0.2f",inc_d_case$lower),", ",
                           sprintf("%0.2f",inc_d_case$upper),")")
inc_d_case$group <- inc_d_case$E_group
inc_d_case <- inc_d_case[,c("group","case","pyear")]
inc_d_case <- inc_d_case[order(inc_d_case$group),]
inc_d_case_liver <- t(inc_d_case)
######figure table 1.2#####
write.csv(inc_d_case_liver,paste0(project_path,"result/fig1table_liver",Date))
###eFigure 1.2: rcs plot for HR in observatioanal  ##################
plotx <- "Alcohol consumption (unit/week)"
covar_f_name <- c("age","sex","APOE4_status","area","smoke", 
                  "townsend_group_1","edu_new")
#,"BMI_group","smoke","sleep_group","PA_group"
coxplot_data <- modeldata %>% select("years","outcome","Exposure",covar_f_name)
coxplot_data <- na.omit(coxplot_data)
coxplot_data$x <- coxplot_data$Exposure
g <- ggplot(coxplot_data, aes(x=Exposure)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=12,), 
        axis.title=element_text(size=12))+
  geom_histogram(binwidth=10,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.4) + 
  scale_x_continuous(name=plotx,expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="Participants",expand = c(0,0),n.breaks=10) +
  scale_fill_jama()
#Cox model
c_d1 <- class.ind(coxplot_data$edu_new)
c_d1 <- c_d1[,-dim(c_d1)[2]]

c_d2 <- class.ind(coxplot_data$area)
c_d2 <- c_d2[,-dim(c_d2)[2]]

c_d3 <- class.ind(coxplot_data$townsend_group_1)
c_d3 <- c_d3[,-dim(c_d3)[2]]

c_d4 <- class.ind(coxplot_data$smoke)
c_d4 <- c_d4[,-dim(c_d4)[2]]

ns <- 6
coxplot_data <- cbind(coxplot_data[,1:ns],c_d1, c_d2, c_d3, c_d4)
names(coxplot_data)[(ns+1):(ns+2)] <- paste0("edu",1:2)
names(coxplot_data)[(ns+5):(ns+6)] <- paste0("townsend",1:2)
names(coxplot_data)[(ns+7):(ns+8)] <- paste0("smoke",1:2)

fun <- "Surv(years,outcome) ~ rcs(Exposure,3)"
fun1 <- "Surv(years,outcome) ~ Exposure"
for(o in c(4:(ns+8))){ #29
  fun <- paste(fun,"+",names(coxplot_data)[o])
  fun1 <- paste(fun1,"+",names(coxplot_data)[o])
}
dd <- datadist(coxplot_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=coxplot_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=coxplot_data)
dd$limits$Exposure[2] <- 0
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
p_allover <- p_test[1,3]

p_allover <- p_out(p_allover)
p_nonliner <- p_out(p_nonliner)

HR <- Predict(m0,Exposure, fun=exp,ref.zero=TRUE)
#plot2
label_text <- substitute(expr=paste("Overall ",italic("P"),p_allover1, 
                                    ", Nonlinear ",italic("P"), p_nonliner1),
                         env = base::list(p_allover1=p_allover,
                                          p_nonliner1=p_nonliner))
#x_text <- max(HR$Exposure)-(max(HR$Exposure)-min(HR$Exposure))*0.3
g0 <- ggplot() +
  geom_line(data=HR, aes(Exposure,yhat),
            linetype="solid",size=1,alpha = 0.3,colour="red3") +
  geom_ribbon(data=HR, aes(Exposure,ymin = lower, ymax = upper),
              alpha = 0.2,fill="red3") +
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA),
                                   axis.text = element_text(color="black", size=12),
                                   axis.title=element_text(size=12),
                                   legend.position="top") +
  #annotate("text", x=x_text, y=max(HR$upper), label=label_text, cex = 4) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") + 
  scale_x_continuous(name=plotx,expand = c(0,0), n.breaks = 10)  +
  scale_y_continuous(name="HR of Alcohol liver disease", n.breaks = 10 )
G <- g+labs(subtitle = label_text)
plot_rcs_HR_liver <- ggplot2.two_y_axis(G,g0)
######figure1.3#####
grid.arrange(plot_rcs_HR_liver)
#####PART Three#############################################################
#################fig2 plot non-linear MR################
#load data
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
plotx <- "genetic prediction of alcohol consumption (unit/week)"
y_time <- "years"
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "pHR of Alcohol liver diseases"

qn <- c(10,100,250,500)
RD_nlMR_liver <- list()
i=1
for(qq in qn){
  plot_x_name <- paste0(plot_x," (",qq," strata",")")
  f <-  sumNLMR_figure_rank_HR(modeldata,puplation,plot_x_name,plot_y,1,qq)
  RD_nlMR_liver[[i]] <- f
  i = i + 1
}
##residual 
nlMR_liver <- list()
i=1
for(qq in qn){
  plot_x_name <- paste0(plot_x," (",qq," strata",")")
  f <-  sumNLMR_figure_res_HR(modeldata,puplation,plot_x_name,plot_y,1,qq)
  nlMR_liver[[i]] <- f
  i = i + 1
}
######figure2.1 doubleR nonMR#####
grid.arrange(RD_nlMR_liver[[1]])
######figure2.2 nonMR#####
grid.arrange(nlMR_liver[[1]])


#####
####box figure
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
ex_name <- "perweek_alchol_unit_cal_new"
y_name <- "outcome"
y_time <- "years"
x_prs_name <- "wPRS_95"
DRnlmr_10strata_liver <- sumNLMR_rank_HR(modeldata,10)
sum_DRnlmr_liver <- DRnlmr_10strata_liver[[1]]
nlmr_model1_DRnlmr_liver <- DRnlmr_10strata_liver[[2]]

nlmr_10strata_liver <-sumNLMR_res_HR(modeldata,10)
sum_nlmr_liver <- nlmr_10strata_liver[[1]]
nlmr_model_nlmr_liver <- nlmr_10strata_liver[[2]]


DRnlmr_liver <- nmlr_boxfigure(nlmr_model1_DRnlmr_liver$figure$data)
nlmr_boxF_liver <- nmlr_boxfigure(nlmr_model_nlmr_liver$figure$data)
