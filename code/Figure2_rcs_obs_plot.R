#########Figure 2 rcs for all###
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)

#########Figure 2 rcs for all#################
DATA <- data_all_analysis
##target exposure
DATA$Exposure <- log10(DATA$perweek_alchol_unit_cal_new)
####Figure for all
modeldata <- DATA
covar_f_name <- c("area","townsend_group_1","edu_new")
                  # "BMI_group","smoke","sleep_group","PA_group","CVD_before","stroke_before")
##rcs model
coxplot_data <- modeldata %>% select("years","outcome","Exposure",
                                     "age","sex","APOE4_status",
                                     all_of(covar_f_name))
coxplot_data <- na.omit(coxplot_data)

#############figure distribution of alcohol consumption#############
f1_data <- coxplot_data
f1_data$x <- 10^coxplot_data$Exposure

#############Cox model#############
##transfored the category data to factor
ns = 6
cox_data <- coxplot_data[,1:ns]
for(a in covar_f_name){
  dummy <- class.ind(coxplot_data[[a]])
  dummy <- dummy[, -ncol(dummy)]
  colnames(dummy) <- paste0(a,1:ncol(dummy))
  cox_data <- cbind(cox_data,dummy)
}

model_var <- names(cox_data)[4:dim(cox_data)[2]]
fun <- paste("Surv(years,outcome) ~ rcs(Exposure,3) + ",paste(model_var,collapse=" + "))
fun1 <- paste("Surv(years,outcome) ~ Exposure + ",paste(model_var,collapse=" + "))
dd <- datadist(cox_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=cox_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=cox_data)
dd$limits$Exposure[2] <- 1.077285
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
p_allover <- p_test[1,3]
p_allover <- p_out(p_allover)
p_nonliner <- p_out_fig(p_nonliner)
predicted_HR <- Predict(m0, Exposure=seq(min(cox_data$Exposure), max(cox_data$Exposure), 
                                         length=dim(table(cox_data$Exposure))),fun=exp, ref.zero=TRUE)

min_HR <- which.min(predicted_HR$y) # Finds the index of the minimum HR
lowest_point <- predicted_HR$Exposure[min_HR] # Finds the corresponding Exposure value 1.077285
lowest_point
10^lowest_point

# 新加的部分：Schoenfeld残差检验
test_cox <- coxph(Surv(years,outcome) ~ rcs(Exposure,3),data = cox_data)
test_cox <- coxph(Surv(years,outcome) ~ Exposure,data = cox_data)
schoenfeld_test <- cox.zph(test_cox)
print(schoenfeld_test)
# 可视化Schoenfeld残差对时间的图
plot(schoenfeld_test)

##interaction between exposure and time
#fun_inter <- paste("Surv(years,outcome) ~ rcs(Exposure, 3) * years + ",paste(model_var,collapse=" + "))
#cph(as.formula(fun_inter), x=TRUE, y=TRUE, data=cox_data,iter.max = 100)
#
#fun_inter <- paste("Surv(years,outcome) ~ Exposure * years + ",paste(model_var,collapse=" + "))
#fun_inter <- paste("Surv(years,outcome) ~ Exposure * years + ",paste(model_var,collapse=" + "))

#cph(as.formula(fun_inter), x=TRUE, y=TRUE, data=cox_data)
#fun_inter <- paste("Surv(years,outcome) ~ rcs(Exposure, 3) + Exposure*years + ",paste(model_var,collapse=" + "))
#cph(as.formula(fun_inter), x=TRUE, y=TRUE, data=cox_data)


#############figure 2 data#############
f2_data <- predicted_HR
f2_data$x <- 10^f2_data$Exposure

ref_line <- f2_data$x[f2_data$yhat==min(f2_data$yhat)]

xbreak = c(seq(0,70,7))
xlimit = 75 
ylimit = 37000
###hist
f1 <- ggplot(f1_data, aes(x=x)) + 
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(9, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=11, color="grey50", hjust=1)
  ) + 
  geom_histogram(binwidth=2,colour="white",
                 fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.1) + 
  scale_x_continuous(name="Alcohol cunsumption (unit/week)",
                     limits=c(0,xlimit),
                     expand = c(-0.01,0.03),
                     breaks=xbreak) +
  scale_y_continuous(
    name="No of participants",
    expand = c(0,0),breaks=c(seq(5000,37000,5000))) +
  coord_cartesian(xlim=c(0, xlimit),ylim=c(0, ylimit)) 
f1

###rcs
f2_data$lower_new <- f2_data$lower
#f2_data$lower_new[f2_data$lower_new<0.95] = 0.95

f2_data$upper_new <- f2_data$upper
#f2_data$upper_new[f2_data$upper_new>1.54] = 1.54


f2 <- ggplot() +
  geom_ribbon(data=f2_data, aes(x,ymin = lower_new, ymax = upper_new),
              alpha = 0.3,fill="#6653A3") +
  geom_line(data=f2_data, aes(x,yhat),
            linetype="solid",size=0.8,alpha = 1,colour="#6653A3") +
  geom_vline(xintercept = ref_line,color="grey60",linetype="dashed")+
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(9, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=11, color="grey50", hjust=1)
  )+  
  geom_hline(yintercept=1, linetype="dashed", color = "grey50") + 
  scale_x_continuous(name="Alcohol cunsumption (unit/week)",
                     limits=c(0,xlimit),
                     expand = c(0,0),
                     breaks=xbreak)+
  scale_y_continuous(name="Hazard ratio (95% CI) for dementia",expand = c(0,0),
                     limits=c(0.95,1.54),
                     breaks=seq(1,1.5,0.1)
  ) +
  scale_fill_jama()  
f2 


########
label_text = "Overall"
label_text_nonliner <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                                  env = base::list(p_nonliner1=p_nonliner))
F2 = f2 + labs(subtitle = label_text_nonliner, title=label_text ) 
F1 = f1 + labs(subtitle = label_text_nonliner, title=label_text ) 
F2

Figur1_overall <- combind_figs(F1,F2)
grid.newpage()
grid.draw(g1_2)
grid.draw(Figur1_overall)





#########Figure 2 rcs for Men#########
DATA <- data_all_analysis %>% filter(sex==1)
##target exposure
DATA$Exposure <- log10(DATA$perweek_alchol_unit_cal_new)
####Figure for all
modeldata <- DATA
covar_f_name <- c("area","townsend_group_1","edu_new")
# "BMI_group","smoke","sleep_group","PA_group")
##rcs model
coxplot_data <- modeldata %>% select("years","outcome","Exposure",
                                     "age","APOE4_status",
                                     all_of(covar_f_name))
coxplot_data <- na.omit(coxplot_data)

#############figure distribution of alcohol consumption#############
f1_data <- coxplot_data
f1_data$x <- 10^coxplot_data$Exposure

#############Cox model#############
##transfored the category data to factor
ns = 5
cox_data <- coxplot_data[,1:ns]
for(a in covar_f_name){
  dummy <- class.ind(coxplot_data[[a]])
  dummy <- dummy[, -ncol(dummy)]
  colnames(dummy) <- paste0(a,1:ncol(dummy))
  cox_data <- cbind(cox_data,dummy)
}

model_var <- names(cox_data)[4:dim(cox_data)[2]]
fun <- paste("Surv(years,outcome) ~ rcs(Exposure,3) + ",paste(model_var,collapse=" + "))
fun1 <- paste("Surv(years,outcome) ~ Exposure + ",paste(model_var,collapse=" + "))
dd <- datadist(cox_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=cox_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=cox_data)
dd$limits$Exposure[2] <-1.220186
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
p_allover <- p_test[1,3]
p_allover <- p_out(p_allover)
p_nonliner <- p_out_fig(p_nonliner)
predicted_HR <- Predict(m0, Exposure=seq(min(cox_data$Exposure), max(cox_data$Exposure), 
                                         length=dim(table(cox_data$Exposure))),fun=exp, ref.zero=TRUE)

min_HR <- which.min(predicted_HR$y) # Finds the index of the minimum HR
lowest_point <- predicted_HR$Exposure[min_HR] # Finds the corresponding Exposure value 1.220186
lowest_point
10^lowest_point
#############figure 2 data#############
f2_data <- predicted_HR
f2_data$x <- 10^f2_data$Exposure

ref_line <- f2_data$x[f2_data$yhat==min(f2_data$yhat)]

xbreak = c(seq(0,70,7))
xlimit = 75 
ylimit = 37000

###hist
f1 <- ggplot(f1_data, aes(x=x)) + 
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(10, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=11, color="grey50", hjust=1)
  ) + 
  geom_histogram(binwidth=2,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.1) + 
  scale_x_continuous(name="Alcohol cunsumption (unit/week)",
                     limits=c(0,xlimit),
                     expand = c(-0.01,0.03),
                     breaks=xbreak) +
  scale_y_continuous(
    name="No of participants",
    expand = c(0,0),breaks=c(seq(5000,ylimit,5000))) +
  coord_cartesian(xlim=c(0, xlimit),ylim=c(0, ylimit)) 
f1

###rcs
f2_data$lower_new <- f2_data$lower
f2_data$lower_new[f2_data$lower_new<0.95] = 0.95

f2_data$upper_new <- f2_data$upper
f2_data$upper_new[f2_data$upper_new>1.54] = 1.54


f2 <- ggplot() +
  geom_ribbon(data=f2_data, aes(x,ymin = lower_new, ymax = upper_new),
              alpha = 0.3,fill="#6653A3") +
  geom_line(data=f2_data, aes(x,yhat),
            linetype="solid",size=0.8,alpha = 1,colour="#6653A3") +
  geom_vline(xintercept = ref_line,color="grey60",linetype="dashed")+
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(10, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=11, color="grey50", hjust=1)
  )+  
  geom_hline(yintercept=1, linetype="dashed", color = "grey50") + 
  scale_x_continuous(name="Alcohol cunsumption (unit/week)",
                     limits=c(0,xlimit),
                     expand = c(0,0),
                     breaks=xbreak)+
  scale_y_continuous(name="Hazard ratio (95% CI) for dementia",expand = c(0,0),
                     limits=c(0.95,1.54),
                     breaks=seq(1,1.5,0.1)
  ) +
  scale_fill_jama()  
f2 


########
label_text = "Men"
label_text_nonliner <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                                  env = base::list(p_nonliner1=p_nonliner))
F2 = f2 + labs(subtitle = label_text_nonliner, title=label_text ) 
F1 = f1 + labs(subtitle = label_text_nonliner, title=label_text ) 
F2

Figur1_men <- combind_figs(F1,F2)
grid.newpage()
grid.draw(Figur1_men)









#########Figure 2 rcs for Men#########
DATA <- data_all_analysis %>% filter(sex==0)
##target exposure
DATA$Exposure <- log10(DATA$perweek_alchol_unit_cal_new)
####Figure for all
modeldata <- DATA
covar_f_name <- c("area","townsend_group_1","edu_new")
# "BMI_group","smoke","sleep_group","PA_group")
##rcs model
coxplot_data <- modeldata %>% select("years","outcome","Exposure",
                                     "age","APOE4_status",
                                     all_of(covar_f_name))
coxplot_data <- na.omit(coxplot_data)

#############figure distribution of alcohol consumption#############
f1_data <- coxplot_data
f1_data$x <- 10^coxplot_data$Exposure

#############Cox model#############
##transfored the category data to factor
ns = 5
cox_data <- coxplot_data[,1:ns]
for(a in covar_f_name){
  dummy <- class.ind(coxplot_data[[a]])
  dummy <- dummy[, -ncol(dummy)]
  colnames(dummy) <- paste0(a,1:ncol(dummy))
  cox_data <- cbind(cox_data,dummy)
}

model_var <- names(cox_data)[4:dim(cox_data)[2]]
fun <- paste("Surv(years,outcome) ~ rcs(Exposure,3) + ",paste(model_var,collapse=" + "))
fun1 <- paste("Surv(years,outcome) ~ Exposure + ",paste(model_var,collapse=" + "))
dd <- datadist(cox_data)
options(datadist='dd')
m0 <- cph(as.formula(fun), x=TRUE, y=TRUE, data=cox_data)
m1 <- cph(as.formula(fun1), x=TRUE, y=TRUE, data=cox_data)
dd$limits$Exposure[2] <-0.9237548
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 1- pchisq( (m0$loglik[2]- m1$loglik[2]), 2)
p_allover <- p_test[1,3]
p_allover <- p_out(p_allover)
p_nonliner <- p_out_fig(p_nonliner)
predicted_HR <- Predict(m0, Exposure=seq(min(cox_data$Exposure), max(cox_data$Exposure), 
                                         length=dim(table(cox_data$Exposure))),fun=exp, ref.zero=TRUE)

min_HR <- which.min(predicted_HR$y) # Finds the index of the minimum HR
lowest_point <- predicted_HR$Exposure[min_HR] # Finds the corresponding Exposure value 0.9237548
lowest_point
10^lowest_point
#############figure 2 data#############
f2_data <- predicted_HR
f2_data$x <- 10^f2_data$Exposure

ref_line <- f2_data$x[f2_data$yhat==min(f2_data$yhat)]

xbreak = c(seq(0,70,7))
xlimit = 75 
ylimit = 37000

###hist
f1 <- ggplot(f1_data, aes(x=x)) + 
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(10, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=11, color="grey50", hjust=1)
  ) + 
  geom_histogram(binwidth=2,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.1) + 
  scale_x_continuous(name="Alcohol cunsumption (unit/week)",
                     limits=c(0,xlimit),
                     expand = c(-0.01,0.03),
                     breaks=xbreak) +
  scale_y_continuous(
    name="No of participants",
    expand = c(0,0),breaks=c(seq(5000,ylimit,5000))) +
  coord_cartesian(xlim=c(0, xlimit),ylim=c(0, ylimit)) 
f1

###rcs
f2_data$lower_new <- f2_data$lower
f2_data$lower_new[f2_data$lower_new<0.95] = 0.95

f2_data$upper_new <- f2_data$upper
f2_data$upper_new[f2_data$upper_new>1.54] = 1.54


f2 <- ggplot() +
  geom_ribbon(data=f2_data, aes(x,ymin = lower_new, ymax = upper_new),
              alpha = 0.3,fill="#6653A3") +
  geom_line(data=f2_data, aes(x,yhat),
            linetype="solid",size=0.8,alpha = 1,colour="#6653A3") +
  geom_vline(xintercept = ref_line,color="grey60",linetype="dashed")+
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(10, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=11, color="grey50", hjust=1)
  )+  
  geom_hline(yintercept=1, linetype="dashed", color = "grey50") + 
  scale_x_continuous(name="Alcohol cunsumption (unit/week)",
                     limits=c(0,xlimit),
                     expand = c(0,0),
                     breaks=xbreak)+
  scale_y_continuous(name="Hazard ratio (95% CI) for dementia",expand = c(0,0),
                     limits=c(0.95,1.54),
                     breaks=seq(1,1.5,0.1)
  ) +
  scale_fill_jama()  
f2 


########
label_text = "Women"
label_text_nonliner <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                                  env = base::list(p_nonliner1=p_nonliner))
F2 = f2 + labs(subtitle = label_text_nonliner, title=label_text ) 
F1 = f1 + labs(subtitle = label_text_nonliner, title=label_text ) 
F2

Figur1_women <- combind_figs(F1,F2)
grid.newpage()
grid.draw(Figur1_women)


