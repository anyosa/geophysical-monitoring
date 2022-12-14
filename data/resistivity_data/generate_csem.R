# CSEM resistivity data by Archies law
# From Jos paper: Value of information of seismic amplitude and CSEM resistivity

#library(R.matlab)
today <- format(Sys.Date(), "%d%m%y")
dir.create(file.path('/Users/susanany/phd_research/sequential_voi/code/data/resistivity_data/', today))

load("/Users/susanany/phd_research/sequential_voi/code/data/resistivity_data/porosityData.RData")
load("/Users/susanany/phd_research/sequential_voi/code/data/resistivity_data/cases.RData")

#generate_resistivity_data <- function(saturation_values, porosity_values){
#    number_of_realizations <- nrow(saturation_array)
#    number_of_cells <- ncol(saturation_array)
#  
#    gamma_parameter <- 0.78
#    alpha_parameter <- 2.0  # before 25.02.22: -1.31
#    beta_parameter <- 2.0 # before 25.02.22: -0.14
#    resistivity_array <- matrix(data = NA, nrow = number_of_realizations, ncol = number_of_cells)
#
#    for (i in 1:number_of_realizations) {
#        saturation <- saturation_array[i, ]
#        saturation <- ifelse(saturation < 0.05, 0.05, saturation)
#        saturation <- ifelse(saturation > 0.95, 0.95, saturation)
#        porosity <- porosity_array[i, ]
#        resistivity <- gamma_parameter * saturation**alpha_parameter * porosity**beta_parameter
#        resistivity_array[i, ] <- resistivity
#    }
#    return(resistivity_array)
#}

generate_resistivity_data <- function(saturation_values, porosity_values){
  
  number_of_realizations <- nrow(saturation_array)
  number_of_cells <- ncol(saturation_array)
  a <- 1
  R_brine <- 2/11
  m <- 2.5 # cem.exp 
  n <- 2 # sat.exp
  
  resistivity_array <- matrix(data = NA, nrow = number_of_realizations, ncol = number_of_cells)
  
  for (i in 1:number_of_realizations) {
    saturation <- saturation_array[i, ]
    saturation <- ifelse(saturation < 0.05, 0.05, saturation)
    saturation <- ifelse(saturation > 0.95, 0.95, saturation)
    porosity <- porosity_array[i, ]
    log_resistivity <- log(a*R_brine/((porosity**m)*((1-saturation)**n)))
    resistivity_array[i, ] <- log_resistivity
  }
  return(resistivity_array)
}



# modus 1: only 3 cells are aggregated. modus 2: 9 cells are aggregated.
reduce_data_resolution <- function(resistivity_array, neighborhood_size, n_easting, n_northing, sail_line, modus){
    number_of_realizations <- nrow(resistivity_array)
    number_of_cells <- ncol(resistivity_array)
    
    number_of_sites <- as.integer(n_northing/neighborhood_size)
    #line <- as.integer(n_easting/2)
    reduced_resolution_array <- matrix(data = NA, nrow = number_of_realizations, ncol = number_of_sites)
    
    if(modus == '1'){
      cells_to_aggregate <- neighborhood_size*number_of_sites
      for (i in 1:number_of_realizations) {
        values_in_line <- t(matrix(resistivity_array[i, ], n_easting, n_northing))[, sail_line]
        aggregated_values <- rowMeans(matrix(values_in_line[1:cells_to_aggregate], ncol = neighborhood_size, byrow = TRUE))
        reduced_resolution_array[i, ] <- aggregated_values
      }
    }
    if(modus == '2'){
      cells_to_aggregate <- neighborhood_size*number_of_sites
      from <- sail_line - 1
      to <- sail_line + 1
      for (i in 1:number_of_realizations) {
        matrix_values_in_line <- t(matrix(resistivity_array[i, ], n_easting, n_northing))[, from:to]
        values_in_line <- rowMeans(matrix_values_in_line)
        aggregated_values <- rowMeans(matrix(values_in_line[1:cells_to_aggregate], ncol = neighborhood_size, byrow = TRUE))
        reduced_resolution_array[i, ] <- aggregated_values
      }
    }
    if(modus == '3'){
      cells_to_aggregate <- neighborhood_size*number_of_sites
      from <- sail_line - 2
      to <- sail_line + 2
      for (i in 1:number_of_realizations) {
        matrix_values_in_line <- t(matrix(resistivity_array[i, ], n_easting, n_northing))[, from:to]
        values_in_line <- rowMeans(matrix_values_in_line)
        aggregated_values <- rowMeans(matrix(values_in_line[1:cells_to_aggregate], ncol = neighborhood_size, byrow = TRUE))
        reduced_resolution_array[i, ] <- aggregated_values
      }
    }
    if(modus == '4'){
      cells_to_aggregate <- neighborhood_size*number_of_sites
      from <- sail_line - 2
      to <- sail_line + 2
      for (i in 1:number_of_realizations) {
        matrix_values_in_line <- t(matrix(resistivity_array[i, ], n_easting, n_northing))[, from:to]
        values_in_line <- rowMeans(matrix_values_in_line)
        aggregated_values <- rowMeans(matrix(values_in_line[1:cells_to_aggregate], ncol = neighborhood_size, byrow = TRUE))
        reduced_resolution_array[i, ] <- aggregated_values
      }
    }
  
  return(reduced_resolution_array)
}

# 30.06.22 changed neighborhood_size from 3 to 5
# 14.11.22 increased from-to to cover 5 horizontal cells
#for(k in 1:2){
#    modus <- k
modus <- 2
for (t in c(2,5,8,11,14,17,20,23)) {
    path_in <- '/Users/susanany/phd_research/sequential_voi/code/data/saturation_data/'
    file_in <- paste(path_in, 'Sat_time', t, '_mod.RData', sep = '')
    load(file_in)# array 1000x324
    resistivity_array <- generate_resistivity_data(saturation_values = saturation_array, porosity_values = porosity_array)
    reduced_resolution_array <- reduce_data_resolution(resistivity_array, neighborhood_size = 2, 18, 18, 15, modus)
    path_out <- paste('/Users/susanany/phd_research/sequential_voi/code/data/resistivity_data/', today, '/', sep = '')
    file_out <- paste(path_out, 'resistivity', 9, '_time', t, '.RData', sep = '')
    save(reduced_resolution_array, file = file_out)
}
#}


#%%%%%%%%%%%%%#
#### Cases ####
#%%%%%%%%%%%%%#

path_in <- '/Users/susanany/phd_research/co2_dynamic/code/preprocessing_data/resistivity_data/source_data/'
file_in <- paste(path_in, 'cases.mat', sep = '')
labels_array <- readMat(file_in)$cases
path_out <- '/Users/susanany/phd_research/co2_dynamic/code/preprocessing_data/resistivity_data/output_data/'
file_out <- paste(path_out, 'cases.RData', sep = '')
save(labels_array, file = file_out)

#%%%%%%%%%%%%%%%%#
#### Porosity ####
#%%%%%%%%%%%%%%%%#

path_in <- '/Users/susanany/phd_research/co2_dynamic/code/preprocessing_data/resistivity_data/source_data/'
file_in <- paste(path_in, 'porosityData.mat', sep = '')
porosity_array <- readMat(file_in)[[1]]
path_out <- '/Users/susanany/phd_research/co2_dynamic/code/preprocessing_data/resistivity_data/output_data/'
file_out <- paste(path_out, 'porosityData.RData', sep = '')
save(porosity_array, file = file_out)
