######Supplement Table description the zero Drinkers######
d <- data
conf <-  c("years","perweek_alchol_unit_cal_new","sex",
           "age_group","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","CVD_before","stroke_before",
           "APOE4_status") 
conf <-  c("years","perweek_alchol_unit_cal_new","sex",
           "age","townsend_group_1","edu_new","BMI_group",
           "smoke","sleep_group","PA_group","CVD_before","stroke_before",
           "APOE4_status") 
d$G <- d$alcohol_drinker_status
table(d$G)
d$G[d$alcohol_drinker_status==2 & 
      d$perweek_alchol_unit_cal_new==0] <- 2
d$G[d$alcohol_drinker_status==2 & 
      d$perweek_alchol_unit_cal_new >0 & 
      d$perweek_alchol_unit_cal_new <=14 ] <- 3
d$G[d$alcohol_drinker_status==2 & 
      d$perweek_alchol_unit_cal_new >14 ] <- 4
table(d$G)
##Total
d$group <- 1
a1 <- table_1_description_by_group(conf,d,ng=dim(table(d$group)))

##non drinker, ex-drinker, current drinker
d$group <- d$G + 1
d$group[d$G>2] <- 3
table(d$group)
a2 <- table_1_description_by_group(conf,d,ng=dim(table(d$group)))

##current drinker: 0, <=14, >14 unite/week
d1 <- d %>% filter(G>1)
table(d1$G)
d1$group <- d1$G - 1
a3 <- table_1_description_by_group(conf,d1,ng=dim(table(d1$group)))

###combinded
Re_sup <- data.frame(a1[,c(1,2)],a2[,-1],a3[,-1])
names(Re_sup) <- c("Characteristics","Total","Non-drinkers","Ex-drinkers","All current drinkers",
                   "0 unit/week","<=14 unit/week",">14 unit/week")
Re_sup <- Re_sup[-38,]
