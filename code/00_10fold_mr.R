##########One sample summary data level analysis##########
###load the data
load(paste0(project_path,"data/drink_96snp_w.RData")) ##drink_96snp_w
load(paste0(project_path,"data/drink_96snp.RData")) ##drink_96snp
data <- data_alldrink

data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)
DATA <- data_all_analysis

##SNP information
drink_96snp_name <- names(drink_96snp)[2:97]
drink_95snp_name <- drink_96snp_name[!(drink_96snp_name %in% "rs13024996")]

###########MR analysis#######
###Overall
####10-fold####
seed = 2023
md <- DATA
md$IV <- md$wPRS_95
md$Exposure <- log10(md$perweek_alchol_unit_cal_new)
folds <- createFolds(md$eid,k=10)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))

drink_snp_1 <- data.frame()
dementia_snp_9 <- data.frame()
for(i in 1:10){
  set_drink <- md[folds[[i]],] 
  set_dementia <-md[-folds[[i]],] 
  ##summary exposure and snp
  fun_x <- paste("Exposure ~ var +", paste(covar_f_name,collapse = " + "))
  fun_y <- paste("Surv(years, outcome) ~ var +", paste(covar_f_name,collapse = " + "))
  
  for(s in drink_96snp_name){
    ##x
    modeldata <- set_drink
    modeldata$var <- modeldata[,s]
    fit_x <- lm(as.formula(fun_x),data=modeldata)
    r <- summary(fit_x)
    x_snp <- r$coefficients[2,]
    x_snp <- c(i,s,x_snp)
    names(x_snp) <- c("folds","SNP","xbeta","xse","t","xp")
    drink_snp_1 <- rbind(drink_snp_1,x_snp)
    ##y
    modeldata <- set_dementia
    modeldata$var <- modeldata[,s]
    fit_y <- coxph(as.formula(fun_y), data=modeldata)
    
    r <- summary(fit_y)
    y_snp <- r$coefficients[1,]
    y_snp <- c(i,s,y_snp)
    names(y_snp) <- c("folds","SNP","ybeta","yHR","yse","z","yp")
    dementia_snp_9 <- rbind(dementia_snp_9,y_snp)
  }
}

names(drink_snp_1) <- c("folds","SNP","beta.exposure","se.exposure","t","xp")
names(dementia_snp_9) <-  c("folds","SNP","beta.outcome","HR","se.outcome","z","yp")

drink_snp_1$id.exposure <- "alcohol_log10"
dementia_snp_9$id.outcome <- "dementia"
########overall################
first_dataset <- merge(drink_snp_1,dementia_snp_9,by=c("folds","SNP"))
first_dataset_all <- first_dataset
save(first_dataset_all,file=paste0(project_path,"data/Beta_dataset_all_seed2023.RData"))
###meta-analysis for MR##########
res_mr_sub_first <- data.frame()
for(i in 1:10){
  d <- first_dataset %>% filter(folds==i)
  b_out <-  as.numeric(d$beta.outcome)
  b_exp <- as.numeric(d$beta.exposure)
  se_out <- as.numeric(d$se.outcome)
  se_exp <- as.numeric(d$se.exposure)
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
method_list <- c("IVW","MR-Egger","ME-Egger_inter","Weighted median")
##MR
res_mr_sub_all <- res_mr_sub_first
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
Meta_result <- Meta_result_all
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
##scatter plot
library(ggplot2)
names(first_dataset)
f_dt <- first_dataset
str(f_dt)
for(i in c(3:6,8:12)){
  f_dt[[i]] <- as.numeric(f_dt[[i]])
}
##save data
scate_data_all <- f_dt
Meta_result_all <- Meta_result


############Men#########
####10-fold####
seed = 2023
md <- DATA  %>% filter(sex==1)
md$IV <- md$wPRS_95
md$Exposure <- log10(md$perweek_alchol_unit_cal_new)
folds <- createFolds(md$eid,k=10)
covar_f_name <- c("age","area","array",paste0("component",1:20))

drink_snp_1 <- data.frame()
dementia_snp_9 <- data.frame()
for(i in 1:10){
  set_drink <- md[folds[[i]],] 
  set_dementia <-md[-folds[[i]],] 
  ##summary exposure and snp
  fun_x <- paste("Exposure ~ var +", paste(covar_f_name,collapse = " + "))
  fun_y <- paste("Surv(years, outcome) ~ var +", paste(covar_f_name,collapse = " + "))
  
  for(s in drink_96snp_name){
    ##x
    modeldata <- set_drink
    modeldata$var <- modeldata[,s]
    fit_x <- lm(as.formula(fun_x),data=modeldata)
    r <- summary(fit_x)
    x_snp <- r$coefficients[2,]
    x_snp <- c(i,s,x_snp)
    names(x_snp) <- c("folds","SNP","xbeta","xse","t","xp")
    drink_snp_1 <- rbind(drink_snp_1,x_snp)
    ##y
    modeldata <- set_dementia
    modeldata$var <- modeldata[,s]
    fit_y <- coxph(as.formula(fun_y), data=modeldata)
    
    r <- summary(fit_y)
    y_snp <- r$coefficients[1,]
    y_snp <- c(i,s,y_snp)
    names(y_snp) <- c("folds","SNP","ybeta","yHR","yse","z","yp")
    dementia_snp_9 <- rbind(dementia_snp_9,y_snp)
  }
}

names(drink_snp_1) <- c("folds","SNP","beta.exposure","se.exposure","t","xp")
names(dementia_snp_9) <-  c("folds","SNP","beta.outcome","HR","se.outcome","z","yp")

drink_snp_1$id.exposure <- "alcohol_log10"
dementia_snp_9$id.outcome <- "dementia"
########MEn################
first_dataset <- merge(drink_snp_1,dementia_snp_9,by=c("folds","SNP"))
first_dataset_men <- first_dataset
save(first_dataset_men,file=paste0(project_path,"data/Beta_dataset_men_seed2023.RData"))
###meta-analysis for MR##########
res_mr_sub_first <- data.frame()
for(i in 1:10){
  d <- first_dataset %>% filter(folds==i)
  b_out <-  as.numeric(d$beta.outcome)
  b_exp <- as.numeric(d$beta.exposure)
  se_out <- as.numeric(d$se.outcome)
  se_exp <- as.numeric(d$se.exposure)
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
method_list <- c("IVW","MR-Egger","ME-Egger_inter","Weighted median")
##MR
res_mr_sub_all <- res_mr_sub_first
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
Meta_result_men <- Meta_result
###result_out_oput###
Meta_result <- Meta_result_men
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
##scatter plot
library(ggplot2)
names(first_dataset)
f_dt <- first_dataset
str(f_dt)
for(i in c(3:6,8:12)){
  f_dt[[i]] <- as.numeric(f_dt[[i]])
}
##save data
scate_data_men <- f_dt
Meta_result_men <- Meta_result


############Women#########
####10-fold####
seed = 2023
md <- DATA  %>% filter(sex==0)
md$IV <- md$wPRS_95
md$Exposure <- log10(md$perweek_alchol_unit_cal_new)
folds <- createFolds(md$eid,k=10)
covar_f_name <- c("age","area","array",paste0("component",1:20))

drink_snp_1 <- data.frame()
dementia_snp_9 <- data.frame()
for(i in 1:10){
  set_drink <- md[folds[[i]],] 
  set_dementia <-md[-folds[[i]],] 
  ##summary exposure and snp
  fun_x <- paste("Exposure ~ var +", paste(covar_f_name,collapse = " + "))
  fun_y <- paste("Surv(years, outcome) ~ var +", paste(covar_f_name,collapse = " + "))
  
  for(s in drink_96snp_name){
    ##x
    modeldata <- set_drink
    modeldata$var <- modeldata[,s]
    fit_x <- lm(as.formula(fun_x),data=modeldata)
    r <- summary(fit_x)
    x_snp <- r$coefficients[2,]
    x_snp <- c(i,s,x_snp)
    names(x_snp) <- c("folds","SNP","xbeta","xse","t","xp")
    drink_snp_1 <- rbind(drink_snp_1,x_snp)
    ##y
    modeldata <- set_dementia
    modeldata$var <- modeldata[,s]
    fit_y <- coxph(as.formula(fun_y), data=modeldata)
    
    r <- summary(fit_y)
    y_snp <- r$coefficients[1,]
    y_snp <- c(i,s,y_snp)
    names(y_snp) <- c("folds","SNP","ybeta","yHR","yse","z","yp")
    dementia_snp_9 <- rbind(dementia_snp_9,y_snp)
  }
}

names(drink_snp_1) <- c("folds","SNP","beta.exposure","se.exposure","t","xp")
names(dementia_snp_9) <-  c("folds","SNP","beta.outcome","HR","se.outcome","z","yp")

drink_snp_1$id.exposure <- "alcohol_log10"
dementia_snp_9$id.outcome <- "dementia"
########MEn################
first_dataset <- merge(drink_snp_1,dementia_snp_9,by=c("folds","SNP"))
first_dataset_women <- first_dataset
save(first_dataset_women,file=paste0(project_path,"data/Beta_dataset_women_seed2023.RData"))
###meta-analysis for MR##########
res_mr_sub_first <- data.frame()
for(i in 1:10){
  d <- first_dataset %>% filter(folds==i)
  b_out <-  as.numeric(d$beta.outcome)
  b_exp <- as.numeric(d$beta.exposure)
  se_out <- as.numeric(d$se.outcome)
  se_exp <- as.numeric(d$se.exposure)
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
method_list <- c("IVW","MR-Egger","ME-Egger_inter","Weighted median")
##MR
res_mr_sub_all <- res_mr_sub_first
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
Meta_result_women <- Meta_result
###result_out_oput###
Meta_result <- Meta_result_women
Meta_result$hr <- exp(Meta_result$estimate)
Meta_result$hr_l <- exp(Meta_result$ci.lb)
Meta_result$hr_u <- exp(Meta_result$ci.ub)
Meta_result$HR <- paste0(sprintf("%0.2f",Meta_result$hr)," (",
                         sprintf("%0.2f",Meta_result$hr_l),",",
                         sprintf("%0.2f",Meta_result$hr_u),")")
Meta_result$p <- unlist(lapply(Meta_result$pval,p_out))
##scatter plot
library(ggplot2)
names(first_dataset)
f_dt <- first_dataset
str(f_dt)
for(i in c(3:6,8:12)){
  f_dt[[i]] <- as.numeric(f_dt[[i]])
}
##save data
scate_data_women <- f_dt
Meta_result_women <- Meta_result



'
f_dt0 <- f_dt %>% filter(folds==1)
mrres <- Meta_result %>% filter(m=="FE")
mrres <- mrres[-3,]
mrres$b <- c(0,-0.001820217,0)

ggplot(data=f_dt, aes(x=beta.exposure, y=beta.outcome)) +
  geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
  geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  geom_point(aes(x=beta.exposure, y=beta.outcome,group=folds,colour=folds)) +
  scale_y_continuous(limits = c(-1,1)) +
  scale_x_continuous(limits = c(-0.1,0.2)) +
  geom_abline(aes(intercept=0,slope=0.650107923),colour="#a6cee3") + 

  geom_abline(data=mrres, aes(intercept=b, slope=estimate, colour=Methods), show.legend=TRUE) +
  scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", 
                               "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) #+
  labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
  theme(legend.position="top", legend.direction="vertical") +
  guides(colour=guide_legend(ncol=2)

f_dt$folds <- as.factor(f_dt$folds )
ggplot(data=f_dt, aes(x=beta.exposure, y=beta.outcome)) +
  geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
  geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  geom_point(aes(x=beta.exposure, y=beta.outcome, group=folds, colour=folds)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-0.1, 0.2)) +
  #geom_abline(aes(intercept=0, slope=0.650107923), colour="#a6cee3") +
  geom_abline(data=mrres, aes(intercept=b, slope=estimate, colour=Methods), show.legend=TRUE) +
  scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", 
                               "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#add8e6")) +
  labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
  theme(legend.position="top", legend.direction="vertical") +
  guides(colour=guide_legend(ncol=2))

'
