parent_dir <- '/Users/susanany/phd_research/sequential_voi/code/'
source_name <- paste(parent_dir, 'simulation_case/src/functions.R', sep = '')
source(source_name)
today <- format(Sys.time(), '%d%m%y')
load(paste(parent_dir, 'data/saturation_data/cases.RData', sep = ''))

indexes_for_voi <- new_partition_indexes(labels = labels_array, seed = 1234, size_of_partition = 300)
indexes_for_survey <- new_complement_indexes_balanced(labels = labels_array, indexes_for_voi, seed = 1234, size_of_partition = 100)

# one type of data
# "0: false_positive" "1: true_negative"  "2: true_positive"  "3: false_negative"
cm_tab <- matrix(data = NA, nrow = 5, ncol = 5)
#%%%%%%%%%%%%#
# lower_csem #
#%%%%%%%%%%%%#

load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_dec/list_mc_011222103223_lower_csem.RData")
#load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_nov/list_adaptive_lower_csem_061122.RData")
res <- summarize_results(list_global, indexes_for_survey, labels_array, surveys = 7)
cm <- res$summary
acc <- colSums(cm %>% filter(cat %in% c(1, 2)))[3]
cm
acc
cm_tab[1, ] <- c(0, cm$n, acc)

#%%%%%%%%%%%%%#
# higher_csem #
#%%%%%%%%%%%%%#

load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_dec/list_mc_011222104929_higher_csem.RData")
#load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_nov/list_adaptive_higher_csem_061122.RData")
res <- summarize_results(list_global, indexes_for_survey, labels_array, surveys = 7)
cm <- res$summary
acc <- colSums(cm %>% filter(cat %in% c(1, 2)))[3]
cm
acc
cm_tab[2, ] <- c(0, cm$n, acc)

#%%%%%%%%%%%%%%%#
# lower_seismic #
#%%%%%%%%%%%%%%%#

load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_dec/list_mc_011222112951_lower_seismic.RData")
#load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_nov/list_adaptive_lower_seismic_061122.RData")
res <- summarize_results(list_global, indexes_for_survey, labels_array, surveys = 7)
cm <- res$summary
acc <- colSums(cm %>% filter(cat %in% c(1, 2)))[3]
cm
acc
cm_tab[3, ] <- c(cm$n, acc)

#%%%%%%%%%%%%%%%#
# higher_seismic #
#%%%%%%%%%%%%%%%#

load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_dec/list_mc_011222115110_higher_seismic.RData")
#load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_nov/list_adaptive_higher_seismic_061122.RData")
res <- summarize_results(list_global, indexes_for_survey, labels_array, surveys = 7)
cm <- res$summary
acc <- colSums(cm %>% filter(cat %in% c(1, 2)))[3]
cm
acc
cm_tab[4, ] <- c(cm$n, acc)

#%%%%%#
# all #
#%%%%%#

load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_dec/list_mc_281122154747_all.RData")
#load("/Users/susanany/phd_research/sequential_voi/code/analyzing_results/results_nov/list_adaptive_all_061122.RData")
res <- summarize_results(list_global, indexes_for_survey, labels_array, surveys = 7)
cm <- res$summary
acc <- colSums(cm %>% filter(cat %in% c(1, 2)))[3]
cm
acc
cm_tab[5, ] <- c(cm$n, acc)

#%%%%%%#
# save #
#%%%%%%#

cm_tab
colnames(cm_tab) <- c("0: false_positive", "1: true_negative", "2: true_positive",  "3: false_negative", "accuracy")
rownames(cm_tab) <- c('lower_csem', 'higher_csem', 'lower_seismic', 'higher_seismic', 'all')
file_out <- paste('/Users/susanany/phd_research/sequential_voi/code/analyzing_results/tables_for_paper/table_accuracy_', today, '.RData', sep = '')
save(cm_tab, file = file_out)

