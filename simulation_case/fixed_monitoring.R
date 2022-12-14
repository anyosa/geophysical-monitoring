parent_dir <- '/.../code/'

source_name <- paste(parent_dir, 'simulation_case/src/functions.R', sep = '')
source(source_name)
load(paste(parent_dir, 'data/cases.RData', sep = ''))

indexes_for_voi <- new_partition_indexes(labels = labels_array, seed = 1234, size_of_partition = 300)
indexes_for_survey <- new_complement_indexes_balanced(labels = labels_array, indexes_for_voi, seed = 1234, size_of_partition = 100)

path_out <- paste(parent_dir, 'simulation_case/output/', sep = '')

today <- format(Sys.time(), '%d%m%y')

parent_directory <- parent_dir

###########

list_fixed_plans <- list()
list_fixed_plans[[1]] <- list('survey_vector' = c('lower_seismic', 'lower_csem'),
                              'time_vector' = c(5, 14))
list_fixed_plans[[2]] <- list('survey_vector' = c('higher_seismic', 'higher_csem'),
                              'time_vector' = c(5, 14))
list_fixed_plans[[3]] <- list('survey_vector' = c('lower_seismic', 'lower_csem', 'lower_seismic', 'lower_csem'),
                              'time_vector' = c(5, 8, 11, 14))
list_fixed_plans[[4]] <- list('survey_vector' = c('higher_seismic', 'higher_csem', 'higher_seismic', 'higher_csem'),
                              'time_vector' = c(5, 8, 11, 14))

times_for_fixed_plans <- c()
for(plan in 1:length(list_fixed_plans)){
survey_vector <- list_fixed_plans[[plan]]$survey_vector
time_vector <- list_fixed_plans[[plan]]$time_vector
file_out <- paste(path_out, 'prints_fixed_plan', plan, '_', today, '.txt', sep = '')
list_global <- list()
start <- proc.time()
for(k in 1:100){
  print(k)
  list_global[[k]] <- fixed_monitoring_plan(indexes_for_voi, index_survey = indexes_for_survey[k], survey_vector, time_vector, labels_array, parent_directory)
}
end <- proc.time()
print(end-start)
times_for_fixed_plans[plan] <- (end-start)[3]
file_out <- paste(path_out, 'list_fixed_',  plan, '_', today, '.RData', sep = '')
save(list_global, file = file_out)
}