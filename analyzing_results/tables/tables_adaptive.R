parent_dir <- '/Users/susanany/phd_research/sequential_voi/code/'
source_name <- paste(parent_dir, 'simulation_case/src/functions.R', sep = '')
source(source_name)
today <- format(Sys.time(), '%d%m%y')
load(paste(parent_dir, 'data/saturation_data/cases.RData', sep = ''))

indexes_for_voi <- new_partition_indexes(labels = labels_array, seed = 1234, size_of_partition = 300)
indexes_for_survey <- new_complement_indexes_balanced(labels = labels_array, indexes_for_voi, seed = 1234, size_of_partition = 100)

#load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_nov/list_adaptive_all_061122.RData")
load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_dec/list_mc_281122154747_all.RData")
df <- summarize_results(list_global, indexes_for_survey, labels_array, surveys = 7)$df

#%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Table 6: Summary of stop #
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

df_stop_all <- df %>% filter(binary_stop == 1)
df_stop_all %>% count(survey_1)
df_stop_all %>% count(survey_2)
df_stop_all %>% count(survey_3)
df_stop_all %>% count(survey_4)
df_stop_all %>% count(survey_5)
df_stop_all %>% count(survey_6)
df_stop_all %>% count(vector_stop_time)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Table 7: Summary of continue #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

df_continue_all <- df %>% filter(binary_stop == 0)
df_continue_all %>% count(survey_1)
df_continue_all %>% count(survey_2)
df_continue_all %>% count(survey_3)
df_continue_all %>% count(survey_4)
df_continue_all %>% count(survey_5)
df_continue_all %>% count(survey_6)
df_continue_all %>% count(survey_7)
# Surveys summary for the realizations where we continue injection. At time $t=5$ we acquire higher quality seismic for all the 56 realizations and decide continue injection in all of them. More surveys are done at later time steps and injection continues until the end of the project.

#%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Table 8: False negatives #
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

df_fn_all <- df_continue_all %>% filter(df_continue_all$scenario != df_continue_all$binary_stop)

#%%%%%%%%%%%%%%%%%%%%%%%%%#
# Table 9: False positive #
#%%%%%%%%%%%%%%%%%%%%%%%%%#

df_fp_all <- df_stop_all %>% filter(df_stop_all$scenario != df_stop_all$binary_stop)

file_out <- paste('/Users/susanany/phd_research/sequential_voi/code/analyzing_results/tables_for_paper/tables_adaptive_', today, '.RData', sep = '')
save.image(file = file_out)
