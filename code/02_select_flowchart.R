#####select tha analysis population
####load the combinded data
project_path <- "/home/zoo/project/P1_drinkMR_dementia/"
load(paste0(project_path,"data/DATA_save_502370.RData")) ##DATA_save
#####load pakages####
library(dplyr)
library(ggplot2)
library(readxl)
library(ggsci)
library(gmodels)
library(survival)
library(gmodels)

###select the data accroding to different criteria
data <- DATA_save
################Step 1: choose the data had the genetic ################
#Due to the study aim is to expolore the genetic causal relationship
table(data$drink_snp) #487163
table(is.na(data$drink_snp)) #15207
prop.table(table(is.na(data$drink_snp))) #0.03027052
table(is.na(data$drink_snp),data$ACD) #426
CrossTable(is.na(data$drink_snp),data$ACD,chisq = T,digits=4)
####
data <- data %>% filter(drink_snp==1)

###########Step 2: genetic data quality control###########
#QC1,22019: Sex chromosome aneuploidy marker.This indicates samples which were identified as 
#putatively carrying sex chromosome configurations that are not either XX or XY.
table(is.na(data$qc1)) #651 
prop.table(table(is.na(data$qc1))) #0.001336308   
#QC2,22027: Outliers for heterozygosity or missing rate
#Indicates samples identified as outliers in heterozygosity and missing rates, w
#hich implies that the genotypes for these samples are of poor quality.
table(is.na(data$qc2)) #0
prop.table(table(is.na(data$qc2))) 
#Inconsistent reported and genetic sex 
table(data$sex == data$genectic_sex) #367
prop.table(table(data$sex == data$genectic_sex)) #0.0007533413 
###
data_qc_del <- data %>% filter(is.na(qc1)==FALSE | is.na(qc2)==FALSE| 
                                 sex!=genectic_sex) #missing=651
table(data_qc_del$ACD)
prop.table(table(data_qc_del$ACD))
dim(data_qc_del)
####
data <- data %>% filter(is.na(qc1)==TRUE,is.na(qc2)==TRUE,
                        sex==genectic_sex) #missing=651
table(data$ACD)
prop.table(table(data$ACD))
##genetic_ethnic,22006: Genetic ethnic grouping
CrossTable(data$ethnic_new,is.na(data$genetic_ethnic))
table(is.na(data$genetic_ethnic)) #73890 
table(is.na(data$ethnic_new)) #393597
table(is.na(data$ethnicity)) #393597
table(data$ethnicity)
data$ethnicity[is.na(data$ethnicity)==TRUE] <- 1
table(data$ethnicity,is.na(data$genetic_ethnic))
table(data$ethnicity[is.na(data$genetic_ethnic)==TRUE],data$ACD[is.na(data$genetic_ethnic)==TRUE])
prop.table(table(data$ethnicity[is.na(data$genetic_ethnic)==TRUE & data$ethnicity==1],data$ACD[is.na(data$genetic_ethnic)==TRUE & data$ethnicity==1])) #0.158
prop.table(table(data$ethnicity[is.na(data$genetic_ethnic)==TRUE & data$ethnicity==0],data$ACD[is.na(data$genetic_ethnic)==TRUE & data$ethnicity==0])) #0.158
prop.table(table(data$ACD[is.na(data$genetic_ethnic)==TRUE])) #0.158
table(data$ACD[is.na(data$genetic_ethnic)==TRUE])
prop.table(table(data$ACD[is.na(data$genetic_ethnic)==FALSE])) #0.158
table(data$ACD[is.na(data$genetic_ethnic)==FALSE])
#Indicates samples who self-identified as 'White British' according to (Field 21000) and have very similar genetic 
#ancestry based on a principal components analysis of the genotypes.
table(is.na(data$genetic_ethnic)) #78,278
table(data$genetic_ethnic)
####
data <- data %>% filter(is.na(data$genetic_ethnic)==FALSE)
###########Step 4: baseline with ACD ###########
table(data$ACD_before)  #215
prop.table(table(data$ACD_before)) #0.0004454823  
table(data$ACD)
prop.table(table(data$ACD_before[data$ACD==1])) #0.02338149  
data <- data %>% filter(ACD_before==0)  #482408     
table(data$ACD)
prop.table(table(data$ACD))
###########Step 5: clean alcohol status missing ###########
table(is.na(data$alcohol_drinker_status)) #1654
prop.table(table(is.na(data$alcohol_drinker_status))) #0.0039242
table(is.na(data$alcohol_drinker_status),data$ACD) #1654
prop.table(table(is.na(data$alcohol_drinker_status),data$ACD))
####
data <- data %>% filter(is.na(alcohol_drinker_status)==FALSE) #407726
dim(data)  #407543    
table(data$ACD)
prop.table(table(data$ACD))
data_alldrink <- data


############lifelong abstainers: test the pleitropy effect of the SNP########
load(paste0(project_path,"data/drink_96snp_w.RData")) ##drink_96snp_w
load(paste0(project_path,"data/drink_96snp.RData")) ##drink_96snp
source(paste0(project_path,"code/00_functions.R"))
table(data$alcohol_drinker_status)
table(data$alcohol_drinker_status,data$ACD)
dt <- data %>% filter(alcohol_drinker_status==0)
covar_f_name <- c("age","sex","area","array",paste0("component",1:20))
#96 SNP: drink_96snp
drink_96snp_name <- names(drink_96snp)[2:97]
conf <- c("BMI_group","PA_group","sleep_group",
          "smoke","edu_new","townsend_group_1","edu","APOE4")
Test_var <- data.frame()
for(f in conf){
  for(a in drink_96snp_name){
    test_var <- association_test(f,a,dt)
    Test_var <- rbind(Test_var,test_var)
  }
}
names(Test_var) <- c("Confounders","SNP","beta","se","t","p")
Test_var$p <- as.numeric(Test_var$p)

b_p <- 0.05/(96*7)
Test_var_sig <- Test_var 
Test_var_sig$include <- ifelse(Test_var_sig$p < b_p, 0, 1)
table(Test_var_sig$include )
Test_var_sig <- out_save(3,2,Test_var_sig) #beta
Test_var_sig <- out_save(4,3,Test_var_sig) #se
Test_var_sig <- out_save(5,2,Test_var_sig) #t
Test_var_sig <- out_save(6,6,Test_var_sig) #p
#rs13024996 related to education level need to delete
#rs1229984 related to smoke keep the rs1229984
############
write.csv(Test_var_sig,file=paste0(project_path,"result/sT_test_var_sig.csv"))

############calculate the prs and wprs############
data <- data_alldrink
data <- as.data.frame(data)
drink_96snp_name <- names(drink_96snp)[2:97]
drink_96snpw_name <- names(drink_96snp_w)[2:97]
##96 PRS
data$PRS_96 <- rowSums(data[,drink_96snp_name],na.rm=T)
data$wPRS_96 <- rowSums(data[,drink_96snpw_name],na.rm=T)
plot(data$PRS_96,data$wPRS_96)
##95 SNP without rs13024996
drink_95snp_name <- drink_96snp_name[!(drink_96snp_name %in% "rs13024996")]
length(drink_95snp_name)
drink_95snpw_name <- drink_96snpw_name[!(drink_96snpw_name %in% "rs13024996_w")]
length(drink_95snpw_name)
data$PRS_95 <- rowSums(data[,drink_95snp_name],na.rm=T)
data$wPRS_95 <- rowSums(data[,drink_95snpw_name],na.rm=T)
##rs561222871 specical for men
drink_94snp_name_men <- drink_95snp_name[!(drink_95snp_name %in% "rs561222871")]
drink_94snpw_name_men <- drink_95snpw_name[!(drink_95snpw_name %in% "rs561222871_w")]

data$wPRS_94_men <- rowSums(data[,drink_94snp_name_men],na.rm=T)
data$PRS_94_men <- rowSums(data[,drink_94snpw_name_men],na.rm=T)


##describe the wPRS
summary(data$wPRS_95[data$alcohol_drinker_status==0])
summary(data$wPRS_95[data$alcohol_drinker_status==1])
summary(data$wPRS_95[data$alcohol_drinker_status==2])
summary(data$wPRS_95[data$alcohol_drinker_status==2 & data$perweek_alchol_unit_cal_new==0])

#######save data with all drink #######
data_alldrink <- data
save(data_alldrink,file=paste0(project_path,"data/data_alldrink_407543.RData"))



