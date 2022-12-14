parent_directory <- '/.../repository/'

source_name <- paste(parent_directory, 'simulation_case/src/functions.R', sep = '')
source(source_name)
today <- format(Sys.time(), '%d%m%y%H%M%S')
load(paste(parent_directory, 'data/cases.RData', sep = ''))

indexes_for_voi <- new_partition_indexes(labels = labels_array, seed = 1234, size_of_partition = 300)
indexes_for_survey <- new_complement_indexes_balanced(labels = labels_array, indexes_for_voi, seed = 1234, size_of_partition = 100)

path_out <- paste(parent_directory, 'simulation_case/output/', sep = '')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

data_vector <- c('lower_csem', 'higher_csem', 'lower_seismic', 'higher_seismic')
print(data_vector)

list_prices <- list()
list_prices[[1]] <- data.frame(lower_csem = 0.1)
list_prices[[2]] <- data.frame(higher_csem = 0.15)
list_prices[[3]] <- data.frame(lower_seismic = 0.6)
list_prices[[4]] <- data.frame(higher_seismic = 0.65)
list_prices[[5]] <- data.frame(lower_csem = 0.1, higher_csem = 0.15, lower_seismic = 0.6, higher_seismic = 0.65)

df_price <- data.frame(lower_csem = 0.1, higher_csem = 0.15, lower_seismic = 0.6, higher_seismic = 0.65)

file_out <-  paste(path_out, 'prints_local_', today, '.txt', sep = '')
list_global <- list()

start <- proc.time()
capture.output(
  for (k in 1:10){
    print(k)
    list_global[[k]] <- update_from_survey_multiple_data_types(indexes_for_voi, index_survey = indexes_for_survey[k], data_vector, df_price, labels_array, parent_directory, number_of_cores = 6, starting_time = 5)
  }, file = file_out)
end <- proc.time()
print(end-start)

file_out <- paste(path_out, 'list_local_', today, '.RData', sep = '')
save(list_global, file = file_out)