###Figure table 2: the table of incident across alcohol ##################
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
alcohol_q10 <- quantile(DATA$Exposure,probs=seq(0,1,0.1))
DATA$alcohol_10group <- cut(DATA$Exposure,
                            breaks=alcohol_quintile,
                            include.lowest = TRUE)
table(DATA$alcohol_10group)

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
