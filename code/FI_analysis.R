#### FI function
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)
DATA <- data_all_analysis
DATA$IV <- DATA$wPRS_95
DATA$Y <- DATA$FIscore
DATA$Exposure <- log10(DATA$perweek_alchol_drink_cal_new)

CrossTable(is.na(DATA$Y),DATA$outcome)
d <- DATA %>% filter(is.na(Y)==F)
CrossTable(d$sex,d$outcome)

###rcs
##target exposure
####Figure for all
modeldata <- DATA
covar_f_name <- c("area","townsend_group_1","edu_new")
# "BMI_group","smoke","sleep_group","PA_group","CVD_before","stroke_before")
##rcs model
coxplot_data <- modeldata %>% select("Y","Exposure",
                                     "age","sex","APOE4_status",
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

model_var <- names(cox_data)[3:dim(cox_data)[2]]
fun <- paste("Y ~ rcs(Exposure,3) + ",paste(model_var,collapse=" + "))
fun1 <- paste("Y ~ Exposure + ",paste(model_var,collapse=" + "))
dd <- datadist(cox_data)
options(datadist='dd')
m0 <- ols(as.formula(fun), x=TRUE, y=TRUE, data=cox_data)
m1 <- ols(as.formula(fun1), x=TRUE, y=TRUE, data=cox_data)
dd$limits$Exposure[2] <- -0.6320232
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 0.7898
p_nonliner <- p_out_fig(p_nonliner)
predicted_HR <- Predict(m0, Exposure=seq(min(cox_data$Exposure), max(cox_data$Exposure), 
                                         length=dim(table(cox_data$Exposure))),fun=exp, ref.zero=TRUE)

min_HR <- which.min(predicted_HR$y) # Finds the index of the minimum HR
lowest_point <- predicted_HR$Exposure[min_HR] # Finds the corresponding Exposure value 1.077285
lowest_point
10^lowest_point
#############figure 2 data#############
f2_data <- predicted_HR
f2_data$x <- 10^f2_data$Exposure
ref_line <- f2_data$x[f2_data$yhat==min(f2_data$yhat)]
###rcs
f2_data$lower_new <- f2_data$lower
#f2_data$lower_new[f2_data$lower_new<0.95] = 0.95
f2_data$upper_new <- f2_data$upper
#f2_data$upper_new[f2_data$upper_new>1.54] = 1.54

xlimit=70
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
  scale_y_continuous(name="Coefficient (95% CI) for Cognitive function score",expand = c(0,0),
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
F2

####### Men ######
####Figure for all
modeldata <- DATA %>% filter(sex==1)
covar_f_name <- c("area","townsend_group_1","edu_new")
# "BMI_group","smoke","sleep_group","PA_group","CVD_before","stroke_before")
##rcs model
coxplot_data <- modeldata %>% select("Y","Exposure",
                                     "age","APOE4_status",
                                     all_of(covar_f_name))
coxplot_data <- na.omit(coxplot_data)

#############figure distribution of alcohol consumption#############
f1_data <- coxplot_data
f1_data$x <- 10^coxplot_data$Exposure

#############Cox model#############
##transfored the category data to factor
ns = 4
cox_data <- coxplot_data[,1:ns]
for(a in covar_f_name){
  dummy <- class.ind(coxplot_data[[a]])
  dummy <- dummy[, -ncol(dummy)]
  colnames(dummy) <- paste0(a,1:ncol(dummy))
  cox_data <- cbind(cox_data,dummy)
}

model_var <- names(cox_data)[3:dim(cox_data)[2]]
fun <- paste("Y ~ rcs(Exposure,3) + ",paste(model_var,collapse=" + "))
fun1 <- paste("Y ~ Exposure + ",paste(model_var,collapse=" + "))
dd <- datadist(cox_data)
options(datadist='dd')
m0 <- ols(as.formula(fun), x=TRUE, y=TRUE, data=cox_data)
m1 <- ols(as.formula(fun1), x=TRUE, y=TRUE, data=cox_data)
dd$limits$Exposure[2] <- -0.6320232
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 0.5142
p_nonliner <- p_out_fig(p_nonliner)
predicted_HR <- Predict(m0, Exposure=seq(min(cox_data$Exposure), max(cox_data$Exposure), 
                                         length=dim(table(cox_data$Exposure))),fun=exp, ref.zero=TRUE)

min_HR <- which.min(predicted_HR$y) # Finds the index of the minimum HR
lowest_point <- predicted_HR$Exposure[min_HR] # Finds the corresponding Exposure value 1.077285
lowest_point
10^lowest_point
#############figure 2 data#############
f2_data <- predicted_HR
f2_data$x <- 10^f2_data$Exposure
ref_line <- f2_data$x[f2_data$yhat==min(f2_data$yhat)]
###rcs
f2_data$lower_new <- f2_data$lower
#f2_data$lower_new[f2_data$lower_new<0.95] = 0.95
f2_data$upper_new <- f2_data$upper
#f2_data$upper_new[f2_data$upper_new>1.54] = 1.54

xlimit=70
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
  scale_y_continuous(name="Coefficient (95% CI) for Cognitive function score",expand = c(0,0),
                     limits=c(0.95,1.54),
                     breaks=seq(1,1.5,0.1)
  ) +
  scale_fill_jama()  
f2 
########
label_text = "Men"
label_text_nonliner <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                                  env = base::list(p_nonliner1=p_nonliner))
F2_men = f2 + labs(subtitle = label_text_nonliner, title=label_text ) 
F2_men


####### Women ######
####Figure for all
modeldata <- DATA %>% filter(sex==0)
covar_f_name <- c("area","townsend_group_1","edu_new")
# "BMI_group","smoke","sleep_group","PA_group","CVD_before","stroke_before")
##rcs model
coxplot_data <- modeldata %>% select("Y","Exposure",
                                     "age","APOE4_status",
                                     all_of(covar_f_name))
coxplot_data <- na.omit(coxplot_data)

#############figure distribution of alcohol consumption#############
f1_data <- coxplot_data
f1_data$x <- 10^coxplot_data$Exposure

#############Cox model#############
##transfored the category data to factor
ns = 4
cox_data <- coxplot_data[,1:ns]
for(a in covar_f_name){
  dummy <- class.ind(coxplot_data[[a]])
  dummy <- dummy[, -ncol(dummy)]
  colnames(dummy) <- paste0(a,1:ncol(dummy))
  cox_data <- cbind(cox_data,dummy)
}

model_var <- names(cox_data)[3:dim(cox_data)[2]]
fun <- paste("Y ~ rcs(Exposure,3) + ",paste(model_var,collapse=" + "))
fun1 <- paste("Y ~ Exposure + ",paste(model_var,collapse=" + "))
dd <- datadist(cox_data)
options(datadist='dd')
m0 <- ols(as.formula(fun), x=TRUE, y=TRUE, data=cox_data)
m1 <- ols(as.formula(fun1), x=TRUE, y=TRUE, data=cox_data)
dd$limits$Exposure[2] <- -0.6320232
m0 <- update(m0)
p_test <- anova(m0)
p_nonliner <- 0.7160
p_nonliner <- p_out_fig(p_nonliner)
predicted_HR <- Predict(m0, Exposure=seq(min(cox_data$Exposure), max(cox_data$Exposure), 
                                         length=dim(table(cox_data$Exposure))),fun=exp, ref.zero=TRUE)

min_HR <- which.min(predicted_HR$y) # Finds the index of the minimum HR
lowest_point <- predicted_HR$Exposure[min_HR] # Finds the corresponding Exposure value 1.077285
lowest_point
10^lowest_point

#############figure 2 data#############
f2_data <- predicted_HR
f2_data$x <- 10^f2_data$Exposure
ref_line <- f2_data$x[f2_data$yhat==min(f2_data$yhat)]
###rcs
f2_data$lower_new <- f2_data$lower
f2_data$upper_new <- f2_data$upper
xlimit=70
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
  scale_y_continuous(name="Coefficient (95% CI) for Cognitive function score",expand = c(0,0),
                     limits=c(0.95,1.54),
                     breaks=seq(1,1.5,0.1)
  ) +
  scale_fill_jama()  
f2 
########
label_text = "Women"
label_text_nonliner <- substitute(expr=paste("Non-linear test: ",italic("P"), p_nonliner1),
                                  env = base::list(p_nonliner1=p_nonliner))
F2_women = f2 + labs(subtitle = label_text_nonliner, title=label_text ) 
F2_women

grid.arrange(F2,F2_men,F2_women,ncol=3)


##all
modeldata <- DATA
x="Model1"
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
###linear MR
##IV-exposure association
fun <- "Exposure ~ IV "
for(o in covar_f_name){fun <- paste(fun,"+",o)}
model1 <- lm(as.formula(fun), data=modeldata)
beta_PRS <- summary(model1)$coefficients[2,1]
modeldata$predict <- predict(model1)
fun_out <- "Y ~ predict "
for(o in covar_f_name){fun_out <- paste(fun_out,"+",o)}
fit_out <- glm(as.formula(fun_out), data=modeldata)
mr_result_all <- summary(fit_out)$coefficients[2,]

##men
modeldata <- DATA %>% filter(sex==1)
x="Model1"
covar_f_name <- c("age","area","array",paste0("component",1:20))
###linear MR
##IV-exposure association
fun <- "Exposure ~ IV "
for(o in covar_f_name){fun <- paste(fun,"+",o)}
model1 <- lm(as.formula(fun), data=modeldata)
beta_PRS <- summary(model1)$coefficients[2,1]
modeldata$predict <- predict(model1)
fun_out <- "Y ~ predict "
for(o in covar_f_name){fun_out <- paste(fun_out,"+",o)}
fit_out <- glm(as.formula(fun_out), data=modeldata)
mr_result_men <- summary(fit_out)$coefficients[2,]

##women
modeldata <- DATA %>% filter(sex==0)
x="Model1"
covar_f_name <- c("age","area","array",paste0("component",1:20))
###linear MR
##IV-exposure association
fun <- "Exposure ~ IV "
for(o in covar_f_name){fun <- paste(fun,"+",o)}
model1 <- lm(as.formula(fun), data=modeldata)
beta_PRS <- summary(model1)$coefficients[2,1]
modeldata$predict <- predict(model1)
fun_out <- "Y ~ predict "
for(o in covar_f_name){fun_out <- paste(fun_out,"+",o)}
fit_out <- glm(as.formula(fun_out), data=modeldata)
mr_result_women <- summary(fit_out)$coefficients[2,]

mr_result_FI <- rbind(mr_result_all,mr_result_men,mr_result_women)

##########nonlinear MR Figure#######
DATA_m <- modeldata
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
puplation <- "Drinkers (> 0 unit/week)"
x_prs_name <- "wPRS_95"
plot_x <- "weight 95 PRS"
plot_y <- "FI score"
ex_name <- "Exposure"
y_name <- "Y"
x_name <- "Exposure"
mdata <- DATA_m[,c(x_name,x_prs_name,y_name,covar_f_name)]
mdata <- na.omit(mdata)
num <- dim(mdata)[1]
ym <- mdata[,y_name]
xm <- mdata[,x_name] + 2
x_prsm <- mdata[,x_prs_name]
covar_fm <- mdata[,covar_f_name]
library(devtools)
install_github("jrs95/nlmr")
library(nlmr)
fp = fracpoly_mr(ym,xm,x_prsm,family="gaussian", q=5, nboot=100,fig=T,ref=min(xmean))
fp
sfp = piecewise_mr(ym,xm,x_prsm,family="gaussian", q=5, nboot=100,fig=T,ref=min(xmean))
fp
summ_data <- create_nlmr_summary(y = ym,
                                 x = xm,
                                 g = x_prsm,
                                 covar = covar_fm,
                                 family = "gaussian",
                                 strata_method = "residual", 
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = 3)
#summ_data$summary$xmean[1] <- 0.001
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   powers = 1,
                                                   method = "FE",
                                                   d=1,
                                                   family="gaussian",
                                                   ci="bootstrap_se",
                                                   average.exposure.associations = TRUE,
                                                   ref=min(xmean),
                                                   nboot=1000,
                                                   seed=200))

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
f_data$X <- 10^f_data$x
plot_nlmr <- ggplot(f_data, aes(x = X)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  geom_ribbon(aes(X,ymin = pHR_low, ymax = pHR_up),
              alpha = 0.3,fill="#6653A3") +
  geom_line(aes(y = pHR), color = "#6653A3") +
  
  scale_y_continuous(name="Hazard ratio (95% CI) of dementia",
                     expand = c(0,0),
                     limits = c(0,10,1),
                     breaks = seq(1,11,2)
  ) +
  scale_x_continuous(name=expression("Predicted alcohol consumption (unit/week)"),
                     expand = c(0,0), 
                     limits = c(1.65,35),
                     breaks = seq(5,30,5)
  )   +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(9, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=10, color="grey50", hjust=0)
  )+ 
  labs(subtitle = label_text,
       title = "Overall") 
#####non liear residual MR in ovearll####
plot_nlmr


