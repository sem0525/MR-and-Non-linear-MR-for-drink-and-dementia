#########supplment figure 1 Alcohol-GS with alcohol consumption
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)
DATA <- data_all_analysis
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
  xlab("Alcohol-GS") +
  scale_color_lancet() +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(9, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=10, color="grey50", hjust=0)
  )+
  scale_fill_jama()


g2 <- ggplot(modeldata, aes(x=IV)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=10,), 
        axis.title=element_text(size=10))+
  geom_histogram(binwidth=0.02,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.2) + 
  scale_x_continuous(name="Alcohol consumption weight PRS",expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="No of participants",
                     expand = c(0,0),
                     limit=c(0,41000),
                     breaks = c(2000,seq(10000,50000,10000))
                     ) +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(9, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=10, color="grey50", hjust=0)
  )+
  scale_fill_jama()

########
label_text = "Overall"

g2 = g2 + labs(title=label_text ) 
g1 = g1 + labs(title=label_text ) 

plot_PRS_q5_all <- ggplot2.two_y_axis(g1,g2)
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
  xlab("Alcohol-GS") +
  scale_color_lancet()  +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(9, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=10, color="grey50", hjust=0)
  )+
  scale_fill_jama()



g2 <- ggplot(modeldata, aes(x=IV)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=10,), 
        axis.title=element_text(size=10))+
  geom_histogram(binwidth=0.02,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.2) + 
  scale_x_continuous(name="Alcohol consumption weight PRS",expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="No of participants",expand = c(0,0),
                     limit=c(0,41000),
                     breaks = c(2000,seq(10000,50000,10000))) +
  scale_fill_jama()

########
label_text = "Men"

g2 = g2 + labs(title=label_text ) 
g1 = g1 + labs(title=label_text ) 

plot_PRS_q5_men <- ggplot2.two_y_axis(g1,g2)
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
                     breaks=c(0.000,0.025,0.050,0.075,0.100)) +
  xlab("Alcohol-GS") +
  scale_color_lancet() +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text = element_text(color=RColorBrewer::brewer.pal(9, "Greys")[7], size=11), 
        axis.title = element_text(size=12),
        plot.title= element_text(face="bold", size=12.5, color="grey20", hjust=0.5),
        plot.subtitle=element_text(face="bold", size=10, color="grey50", hjust=0)
  )+
  scale_fill_jama()
g2 <- ggplot(modeldata, aes(x=IV)) + ylab("Count") +
  theme_classic() %+replace%
  theme(panel.background = element_rect(fill = NA), 
        axis.text= element_text(color="black", size=10,), 
        axis.title=element_text(size=10))+
  geom_histogram(binwidth=0.02,colour="white",
                 #fill=RColorBrewer::brewer.pal(4, "Greys")[4],
                 alpha = 0.2) + 
  scale_x_continuous(name="Alcohol consumption weight PRS",expand = c(0,0),n.breaks=10) + 
  scale_y_continuous(name="No of participants",expand = c(0,0),
                     limit=c(0,41000),
                     breaks = c(2000,seq(10000,50000,10000))) +
  scale_fill_jama()
########
label_text = "Women"

g2 = g2 + labs(title=label_text ) 
g1 = g1 + labs(title=label_text ) 
plot_PRS_q5_women <- ggplot2.two_y_axis(g1,g2)
####plot_PRS_women####
grid.arrange(plot_PRS_q5_women)
###Plot_PRS##########
grid.arrange(plot_PRS_q5_all,plot_PRS_q5_men,plot_PRS_q5_women,ncol=3)

