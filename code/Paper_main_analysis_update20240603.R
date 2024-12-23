##########One sample summary data level analysis##########
###Setting enviroment###
project_path <- "/home/zoo/project/P1_drinkMR_dementia/"
today <- Sys.Date()
Date <- paste0("_drinker_",today) ##313958

### Load the data
load(paste0(project_path,"data/data_alldrink_407543.RData")) ##407983

### Load the packages
library(dplyr)
library(survival)
library(ggsci)
library(nnet)
library(rms)
library(gridExtra)
library(gmodels)
library(caret)
library(dplyr)
source(paste0(project_path,"code/00_functions.R"))
source(paste0(project_path,"code/00_mr_function.R"))


data <- data_alldrink
data_all_analysis <- data_alldrink %>% filter(data$alcohol_drinker_status==2,
                                              data$perweek_alchol_unit_cal_new>0)
## Check for the comments
DATA <- data_all_analysis
DATA$Exposure <- DATA$perweek_alchol_unit_cal_new
table(is.na(DATA$lossfollow_date))
prop.table(table(is.na(DATA$lossfollow_date)))
DATA$alcohol_21 <- ifelse(DATA$perweek_alchol_unit_cal_new<21,0,1)
CrossTable(DATA$alcohol_21,DATA$sex)

######Supplement Table 1 description all participant with alcohol and genes######
source(paste0(project_path,"code/sTable1_description_by_group.R"))
write.csv(Re_sup,paste0(project_path,"result/sTable1_characterist",Date,".csv"))

######Table 1 description current participant by sex######
source(paste0(project_path,"code/Table1_description_by_group.R"))
write.csv(table1_result,paste0(project_path,"result/Table1_characterist",Date,".csv"))

######Supplement Table 4&5&6 continous Alcohol-GS with varibales in abstainer and drinekrs#####
source(paste0(project_path,"code/sTable4&5&6_GS_test_abstainer_drinker.R"))
write.csv(sTable_snp_abstainer,paste0(project_path,"result/sTable4_snp_abstainer",Date,".csv"))
write.csv(sTable_PRS_abstainer,paste0(project_path,"result/sTable5_PRS_abstainer",Date,".csv"))
write.csv(sTable_PRS_drinker,paste0(project_path,"result/sTable6_PRS_drinker",Date,".csv"))

######Supplement Table 7-9 qunitle Alcohol-GS with varibales in drinekrs#####
######sTable_qPRS###########
source(paste0(project_path,"code/sTable_Q5_GS_variables_drinkers.R"))
write.csv(Result_q_all,paste0(project_path,"result/stable_Qprs",Date))

#######sTale: PRS with assocition with X,Y and confoudners#######
source(paste0(project_path,"code/sTable_GS_variables.R"))
write.csv(Result_continous_all,paste0(project_path,"result/stable_GS_variable_drinkers",Date))


#######COX model: non-linear#############
######Figure 1 cox regression of aclohol consumption with dementia #####
#source(paste0(project_path,"code/Figure2_incidence_table.R"))
#write.csv(sTable_PRS_drinker,paste0(project_path,"result/sTable6_PRS_drinker",Date,".csv"))
##figure
source(paste0(project_path,"code/functions/combind_figs.R"))
source(paste0(project_path,"code/Figure2_rcs_obs_plot.R"))
grid.arrange(Figur1_overall,Figur1_men,Figur1_women,ncol=3)
#subgroup
##safe drinkers
#source(paste0(project_path,"code/Figure2_rcs_obs_plot_safe.R"))
#grid.arrange(Figur1_all_safe,Figur1_men_safe,Figur1_women_safe,ncol=3)
##unsafe drinkers
#source(paste0(project_path,"code/Figure2_rcs_obs_plot_unsafe.R"))
#grid.arrange(Figur1_all_unsafe,Figur1_men_unsafe,Figur1_women_unsafe,ncol=3)
###combind plot
#grid.arrange(Figur1_all_safe,Figur1_men_safe,Figur1_women_safe,Figur1_all_unsafe,Figur1_men_unsafe,Figur1_women_unsafe,ncol=3)


######sFigure: non-linear MR############
source(paste0(project_path,"code/sFigure1_GS.R"))
grid.arrange(plot_PRS_q5_all,plot_PRS_q5_men,plot_PRS_q5_women,ncol=3)
grid.arrange(plot_PRS_q5_all,plot_PRS_q5_men,plot_PRS_q5_women,ncol=1)


#######Residual non-linear MR#############
###Figure2 non linear residual MR#####
source(paste0(project_path,"code/Figure3_nonlinearMR.R"))
grid.arrange(plot_nlmr,plot_nlmr_men,plot_nlmr_women,ncol=3)
###sFigure non linear residual MR#####
source(paste0(project_path,"code/sFigure2_nlMR_res.R"))
grid.arrange(plot_nlmr,plot_nlmr_men,plot_nlmr_women,ncol=3)
###sFigure non linear residual MR with positive control#####
source(paste0(project_path,"code/sFigure2_nlMR_res_Pos_control.R"))
grid.arrange(plot_nlmr_res_log_pos,plot_nlmr_res_log_pos_men,
             plot_nlmr_res_log_pos_women,ncol=3)
###sFigure non linear residual MR with negative control#####
source(paste0(project_path,"code/sFigure2_nlMR_res_Neg_control.R"))
grid.arrange(plot_nlmr_res_log_neg,plot_nlmr_res_log_neg_men,
             plot_nlmr_res_log_neg_women,ncol=3)

#######Rank non-linear MR#############
###sFigure non linear ranked MR#####
source(paste0(project_path,"code/sFigure_nonlinearMR_ranked.R"))
grid.arrange(plot_nlmr_ranked_log,plot_nlmr_ranked_log_men,plot_nlmr_ranked_log_women,ncol=3)
###sFigure non linear ranked MR with positive control#####
source(paste0(project_path,"code/sFigure2_nlMR_rank_Pos_control.R"))
grid.arrange(plot_nlmr_ranked_log_pos,plot_nlmr_ranked_log_pos_men,
             plot_nlmr_ranked_log_pos_women,ncol=3)
###sFigure non linear ranked MR with negative control#####
source(paste0(project_path,"code/sFigure2_nlMR_rank_Neg_control.R"))
grid.arrange(plot_nlmr_ranked_log_neg,plot_nlmr_ranked_log_neg_men,
             plot_nlmr_ranked_log_neg_women,ncol=3)


##############Linear MR#####################
##R2 in method to instruduce
source(paste0(project_path,"code/sTable10_GS_F_R2.R"))
write.csv(All_R2,paste0(project_path,"result/stable_All_R2",Date))
##


