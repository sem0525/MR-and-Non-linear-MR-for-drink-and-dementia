##########nonlinear MR Figure#######
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)
DATA <- data_all_analysis
DATA$Exposure <- log10(DATA$perweek_alchol_unit_cal_new)


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
                                 strata_method = "ranked", 
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
f_data$X <- 10^f_data$x
plot_nlmr <- ggplot(f_data, aes(x = X)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  geom_ribbon(aes(X,ymin = pHR_low, ymax = pHR_up),
              alpha = 0.3,fill="#6653A3") +
  geom_line(aes(y = pHR), color = "#6653A3") +
  
  scale_y_continuous(name="Hazard ratio (95% CI)",
                     expand = c(0,0),
                     limits = c(0,10,1),
                     breaks = seq(1,40,2)
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
                                 strata_method = "ranked", 
                                 report_het=TRUE,
                                 extra_statistics=TRUE,
                                 q = 10)
model <- with(summ_data$summary, frac_poly_summ_mr(by, bx, byse, bxse, xmean, 
                                                   fig=TRUE,
                                                   powers = 1,
                                                   #d=1,
                                                   family="ranked",
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
f_data$X <- 10^f_data$x
plot_nlmr_men <- ggplot(f_data, aes(x = X)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  geom_ribbon(aes(X,ymin = pHR_low, ymax = pHR_up),
              alpha = 0.3,fill="#6653A3") +
  geom_line(aes(y = pHR), color = "#6653A3") +
  scale_y_continuous(name="Hazard ratio (95% CI)",
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
       title = "Men")
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
                                 strata_method = "ranked", 
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
f_data$X <- 10^f_data$x
f_data$pHR_up_new <- f_data$pHR_up
f_data$pHR_up_new[f_data$pHR_up_new>10] <- 10

plot_nlmr_women <- ggplot(f_data, aes(x = X)) +
  geom_hline(aes(yintercept = 1),linetype="dotted", color = "grey60") +
  
  geom_ribbon(aes(X,ymin = pHR_low, ymax = pHR_up_new),
              alpha = 0.3,fill="#6653A3") +
  geom_line(aes(y = pHR), color = "#6653A3") +
  scale_y_continuous(name="Hazard ratio (95% CI)",
                     expand = c(0,0),
                     limits = c(0,10,1),
                     breaks = seq(1,11,2)
  ) +
  scale_x_continuous(name=expression("Predicted alcohol consumption (unit/week)"),
                     expand = c(0,0), 
                     breaks = seq(5,50,5)
  )   +
  
  
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(9, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=10, color="grey50", hjust=0)
  )+ 
  labs(subtitle = label_text,
       title = "Women")
#####non liear residual MR in women####
plot_nlmr_women


plot_nlmr_ranked_log <- plot_nlmr
plot_nlmr_ranked_log_men <- plot_nlmr_men
plot_nlmr_ranked_log_women <- plot_nlmr_women

