library(dplyr)
library(parallel)

select_locations_by_variance <- function(dataset, locations_to_keep){
  number_of_features <- ncol(dataset) - 2
  dataset <- dataset[ , 1:number_of_features]
  vars <- as.numeric(unlist(lapply(dataset, var)))
  locs <- seq(1, length(vars))
  df_locs <- data.frame(locs = locs,
                        vars = vars)
  df_sorted <- df_locs[order(df_locs$vars, decreasing = TRUE), ]
  indexes <- (df_locs[order(df_locs$vars, decreasing = TRUE), ]$locs)[1:locations_to_keep]
  return(list(indexes = indexes, df_sorted = df_sorted))
}

get_alphas <- function(w1, df_vs){
  w0 <- 1 - w1
  df_vs$w_x<- ifelse(df_vs$scenario == 0, w0, w1)
  df_vs$alphas <- df_vs$v_value*df_vs$w_x 
  return(df_vs)
}

ESS <- function(w1, df_vs){
  w0 <- 1 - w1
  df_vs$w_x<- ifelse(df_vs$scenario == 0, w0, w1)
  df_vs$alphas <- df_vs$v_value*df_vs$w_x 
  ESS <- 1/(sum(df_vs$alphas**2)) 
  return(ESS)
}

fixed_monitoring_plan <- function(indexes_for_voi, index_survey, survey_vector, time_vector, labels_array, parent_directory){
  years <-time_vector
  time_steps <- length(years)
  stop_time <- 0
  list_weights <- list()
  list_v_values <- list()
  survey_data_list <- list()
  for (step in 1:time_steps){
    particles_set <- list()
    print(paste('length list particles=', length(particles_set)))
    survey_set <- list()
    # start block
    t = years[step]
    print(t)
    data_s <- survey_vector[step]
    dataset <- load_realizations(data_s, t, labels_array, parent_directory)
    particles_set[[data_s]] <- new_select_realizations(dataset, indexes_for_voi) 
    survey_set[[data_s]] <- new_select_realizations(dataset, index_survey)
    if(step==1){
      w1 <- 0.25
      df_vs <- new_set_initial_weights(particles_set[[data_s]])
    }
    print(paste('Current w1 =',round(w1,3)))
    print(paste('We gather data:', data_s))
    survey_data_list[[step]] <- data_s
    print(paste('Index of data_y is index =', index_survey))
    sampled_data_y = survey_set[[data_s]] %>% filter(index==index_survey)
    number_of_features <- ncol(sampled_data_y) - 2
    updated_values <- assimilate_data(particles_set[[data_s]], sampled_data_y, w1, df_vs, inverse_of_matrix_R(data_type = data_s, number_of_features = number_of_features))
    w1 <- updated_values$w1
    df_vs <- updated_values$posterior_vs
    chosen_alternative = make_decision(t, w1)
    print(paste('Decision is at time t =', t, 'is',chosen_alternative))
    print(paste('Updated value of w1 =',round(w1,3)))
    list_weights[[step]] <- w1
    list_v_values[[step]] <- df_vs
    if(chosen_alternative == 'stop'){
      stop_time <- t
      print(paste('We choose to stop at time step t =',t))
      break
    }
  }
  #
  return(list(w_values = list_weights, v_values = list_v_values, stop_time = stop_time, surveys = survey_data_list))
}

get_expval <- function(vector_p_x, vector_index, t){
  expval_a_0 <- value_function(x = 0, a = 0, t)*(1 - vector_p_x) + value_function(x = 1, a = 0, t)*vector_p_x
  expval_a_1 <- value_function(x = 0, a = 1, t)*(1 - vector_p_x) + value_function(x = 1, a = 1, t)*vector_p_x
  df_ <- data.frame(p_x = vector_p_x, index = vector_index, ev_a_0 = expval_a_0, ev_a_1 = expval_a_1)
  return(df_)
}

calculate_voi <- function(df_expval, w1, df_vs){
  w0 <- 1-w1
  df_vs$w_x <- ifelse(df_vs$scenario == 0, w0, w1)
  df_merged <- merge(df_vs, df_expval, by = c('index')) %>% arrange(scenario, index)
  df_merged$alphas <- df_merged$v_value*df_merged$w_x
  prior_value <- pmax(weighted.mean(df_merged$ev_a_0, df_merged$alphas), weighted.mean(df_merged$ev_a_1, df_merged$alphas))
  posterior_value <- weighted.mean(pmax(df_merged$ev_a_0, df_merged$ev_a_1), df_merged$alphas)
  voi <- posterior_value - prior_value
  return(voi)
}

inverse_of_matrix_R <- function(data_type, number_of_features){
  if(data_type == 'lower_seismic'){ # R0
    c <- 0.3 # variance_factor
    N <- number_of_features
    inv_R <- solve(c*(0.06**2)*diag(N))
  }
  if(data_type == 'higher_seismic'){ 
    c <- 0.3 # variance_factor
    N <- number_of_features/2
    R0_part <- (0.06**2)*diag(N)
    G_part <- (0.17**2)*diag(N)
    R0_and_G_part <- (-0.7*0.06*0.17)*diag(N)
    joined_matrix <- rbind(cbind(R0_part, R0_and_G_part), cbind(R0_and_G_part, G_part))
    inv_R <- solve(c*joined_matrix)  
  }
  if(data_type == 'lower_csem'){ 
    inv_R <- solve((0.35**2)*diag(number_of_features)) 
  }
  if(data_type == 'higher_csem'){ 
    inv_R <- solve((0.35**2)*diag(number_of_features)) 
  }
  return (inv_R)
}


load_realizations <- function(type_of_data, t, labels_array, parent_directory){
  labels_array <- labels_array
  if(type_of_data == 'lower_csem'){
    path_in <- paste(parent_directory, '/data/resistivity_data/141122/', sep = '')
    file_in <- paste(path_in, 'resistivity2_time', t, '.RData', sep = '')
    load(file_in)
    dataset <- data.frame(csem = reduced_resolution_array, scenario = labels_array, index = seq(1, length(labels_array))) 
  }
  if(type_of_data == 'higher_csem'){
    path_in <- paste(parent_directory, '/data/resistivity_data/141122/', sep = '')
    file_in <- paste(path_in, 'resistivity9_time', t, '.RData', sep = '')
    load(file_in)
    dataset <- data.frame(csem = reduced_resolution_array, scenario = labels_array, index = seq(1, length(labels_array))) 
  }
  if(type_of_data == 'lower_seismic'){
    path_in <- paste(parent_directory, '/data/seismic_data/281122/', sep = '') 
    file_in <- paste(path_in, 'R0_time', t, '.RData', sep = '')
    load(file_in) 
    dataset <- data.frame(seismic = seismic_array, scenario = labels_array, index = seq(1, length(labels_array))) 
  }
  if(type_of_data == 'higher_seismic'){
    path_in <- paste(parent_directory, '/data/seismic_data/281122/', sep = '')
    file_in <- paste(path_in, 'R0_G_time', t, '.RData', sep = '')
    load(file_in)
    dataset <- data.frame(seismic = seismic_array, scenario = labels_array, index = seq(1, length(labels_array))) 
  }
  return(dataset)
}

value_function <- function(x,a,t){
    if(x == 0 && a == 0){v=-25+0.6*t}
    if(x == 1 && a == 0){v=-27-0.6*t}
    if(x == 0 && a == 1){v=-10}
    if(x == 1 && a == 1){v=-42}
    return(v)
}

value_of_information <- function(df_probabilities, t, w1, df_vs){
  
  vector_p_x <- df_probabilities$vector_p_x
  vector_index <- df_probabilities$vector_index
  
  expval_a_0 <- value_function(x = 0, a = 0, t)*(1 - vector_p_x) + value_function(x = 1, a = 0, t)*vector_p_x
  expval_a_1 <- value_function(x = 0, a = 1, t)*(1 - vector_p_x) + value_function(x = 1, a = 1, t)*vector_p_x
  
  w0 <- 1-w1
  df_vs$w_x <- ifelse(df_vs$scenario == 0, w0, w1)
  df_ev <- data.frame(expval_a_0, expval_a_1, vector_p_x, vector_index)
  colnames(df_ev) <- c('ev0', 'ev1', 'p_x', 'index')
  df_merged <- merge(df_vs, df_ev, by = c('index'))
  df_merged$alphas <- df_merged$v_value*df_merged$w_x
  prior_value <- pmax(weighted.mean(df_merged$ev0, df_merged$alphas), weighted.mean(df_merged$ev1, df_merged$alphas))
  posterior_value <- weighted.mean(pmax(df_merged$ev0, df_merged$ev1), df_merged$alphas)
  voi <- posterior_value - prior_value
  return(voi)
}


get_voi <- function(particles, w1, df_vs, inv_R, t, number_of_cores){
  get_results <- function(j){
    p <- nrow(inv_R) 
    probability_x <- w1 
    data_y_j = particles[j, 1:p] 
    df_without_j = particles[-j, ] 
    df_phi = new_gaussian_measurement_model(data_y_j, inv_R, df_without_j) #returns n-1 weights
    df_phi_with_weights = new_reweigh_particle_weights(df_phi, df_vs, j, particles)
    df_u_values = new_compute_unnormalized_joint_probability(probability_x, df_phi_with_weights)
    probability_x = new_classification_probabilities_given_y(df_u_values)$w_1
    particle_index <- particles$index[j]
    return(list(p_x = probability_x, index = particle_index))
  }
  n <- nrow(particles)
  j <- seq(1, n)
  numCores <- number_of_cores
  results <- mclapply(j, get_results, mc.cores = numCores)
  vector_p_x <- sapply(j, function(j) results[[j]]$p_x)
  vector_index <- sapply(j, function(j) results[[j]]$index)
  df_probabilities <- data.frame(vector_p_x, vector_index)
  #vector_p_x <- sapply(j, function(j) results[[j]]$p_x)
  #vector_index <- sapply(j, function(j) results[[j]]$index)
  #df_ev <- get_expval(vector_p_x, vector_index, t)
  #voi <- calculate_voi(df_ev, w1, df_vs)
  voi <- value_of_information(df_probabilities, t, w1, df_vs)
  return(voi)
}

partition_indexes <- function(labels, seed, size_of_partition){
  df <- data.frame(scenario = labels_array,
                   index = seq(1, nrow(labels_array)))
  df_0 <- df %>% filter(scenario == 0)
  df_1 <- df %>% filter(scenario == 1) 
  size_of_sample <- as.integer(size_of_partition/2) #150
  set.seed(seed)
  sample_ids_0 <- sample(df_0$index, size_of_sample, replace = FALSE)
  set.seed(seed)
  sample_ids_1 <- sample(df_1$index, size_of_sample, replace = FALSE)
  return(list(index_0 = sort(sample_ids_0), index_1 = sort(sample_ids_1)))
}

split_realisations <- function(data_frame, seed, size_of_partition){
  df_scenario_0 <- data_frame %>% filter(scenario == 0)
  df_scenario_1 <- data_frame %>% filter(scenario == 1) 
  size_of_sample <- size_of_partition/2 #150
  set.seed(seed)
  sample_ids_0 <- sample(df_scenario_0$index, size_of_sample, replace = FALSE)
  #sample_scenario_0 <- sample.int(n = nrow(df_scenario_0), size = size_of_sample, replace = FALSE)
  set.seed(seed)
  sample_ids_1 <- sample(df_scenario_1$index, size_of_sample, replace = FALSE)
  #sample_scenario_1 <- sample.int(n = nrow(df_scenario_1), size = size_of_sample, replace = FALSE)
  part_scenario_0 <- df_scenario_0 %>% filter(index %in% c(sample_ids_0))
  part_scenario_1 <- df_scenario_1 %>% filter(index %in% c(sample_ids_1))
  df = data.frame(rbind(part_scenario_0, part_scenario_1), row.names = NULL)
  return(df)
}

split_realisations_2 <- function(data_frame, seed, size_of_partition){
  df_scenario_0 <- data_frame %>% filter(scenario == 0)
  df_scenario_1 <- data_frame %>% filter(scenario == 1) 
  size_of_sample <- size_of_partition/2 #150
  set.seed(seed)
  sample_ids_0 <- sample(df_scenario_0$index, size_of_sample, replace = FALSE)
  #sample_scenario_0 <- sample.int(n = nrow(df_scenario_0), size = size_of_sample, replace = FALSE)
  set.seed(seed)
  sample_ids_1 <- sample(df_scenario_1$index, size_of_sample, replace = FALSE)
  #sample_scenario_1 <- sample.int(n = nrow(df_scenario_1), size = size_of_sample, replace = FALSE)
  part_scenario_0 <- df_scenario_0 %>% filter(index %in% c(sample_ids_0))
  part_scenario_1 <- df_scenario_1 %>% filter(index %in% c(sample_ids_1))
  df = data.frame(rbind(part_scenario_0, part_scenario_1), row.names = NULL)
  return(df)
}


select_realizations <- function(df, list_of_indexes){
  sample_ids_0 <- list_of_indexes$index_0
  sample_ids_1 <- list_of_indexes$index_1
  part_scenario_0 <- df %>% filter(index %in% c(sample_ids_0))
  part_scenario_1 <- df %>% filter(index %in% c(sample_ids_1))
  df = data.frame(rbind(part_scenario_0, part_scenario_1), row.names = NULL)
  return(df)
}

complement_indexes <- function(labels, list){
  index_voi <- sort(c(unlist(list$index_0), unlist(list$index_1)))
  df <- data.frame(scenario = labels_array,
                   index = seq(1, nrow(labels_array)))
  df_voi <- df %>% filter(index %in% c(index_voi))
  df_survey <- anti_join(df, df_voi, by = 'index')
  index_survey <- df_survey$index
  return(index_survey)
}

new_gaussian_measurement_model <- function(data_y_j, inv_R, df_without_j){
  p <- length(data_y_j) # number of features
  data_y_j <- data_y_j[, 1:p]
  particles_to_evaluate <- df_without_j[, 1:p]
  n <- nrow(particles_to_evaluate)
  list_phis <- list()
  for(i in 1:n){
    diff <- as.numeric(data_y_j) - as.numeric(particles_to_evaluate[i, ])
    diff <- matrix(diff, nrow = 1)
    phi <- exp(-0.5* diff%*%inv_R%*%t(diff))
    list_phis <- append(list_phis, phi)
  }
  
  df_phis <- data.frame(phis = unlist(list_phis),
                        scenario = df_without_j$scenario,
                        index = df_without_j$index)
  # weights used to be inside df_phis now they are added in the function add_weights
  return (df_phis)
}

new_compute_unnormalized_joint_probability <- function(probability_x, df_phis_w){ # need weights of particles here
  df_scenario_0 <- df_phis_w %>% filter(scenario == 0)
  df_scenario_1 <- df_phis_w %>% filter(scenario == 1) 
  
  u_values_x0 <- (1 - probability_x) * df_scenario_0$v_values * df_scenario_0$phis
  u_values_x1 <- probability_x * df_scenario_1$v_values * df_scenario_1$phis
  
  
  df_u_x0 <- data.frame(u_values = u_values_x0, 
                        scenario = df_scenario_0$scenario,
                        index = df_scenario_0$index, 
                        v_values = df_scenario_0$v_values)
  df_u_x1 <- data.frame(u_values = u_values_x1, 
                        scenario = df_scenario_1$scenario,
                        index = df_scenario_1$index,
                        v_values = df_scenario_1$v_values)
  df_u <- rbind(df_u_x0, df_u_x1)
  df_u <- df_u %>% arrange(scenario, index)
  return(df_u)
}

new_compute_posterior_weights <- function(df_u){
  
  df_scenario_0 <- df_u %>% filter(scenario == 0)
  df_scenario_1 <- df_u %>% filter(scenario == 1) 
  u_x0 <- df_scenario_0$u_values
  u_x1 <- df_scenario_1$u_values
  v_values_x0 = u_x0/sum(u_x0)
  v_values_x1 = u_x1/sum(u_x1)
  df_v_x0 <- data.frame(index = df_scenario_0$index, scenario = df_scenario_0$scenario, v_value = v_values_x0)
  df_v_x1 <- data.frame(index = df_scenario_1$index, scenario = df_scenario_1$scenario, v_value = v_values_x1)
  df_v <- rbind(df_v_x0, df_v_x1)
  df_v <- df_v %>% arrange(scenario, index)
  return (df_v)
}

new_classification_probabilities_given_y <- function(df_u){
  df_scenario_0 <- df_u %>% filter(scenario == 0)
  df_scenario_1 <- df_u %>% filter(scenario == 1) 
  u_x0 <- df_scenario_0$u_values
  u_x1 <- df_scenario_1$u_values
  total <- sum(u_x0) + sum(u_x1)
  w_0 <- sum(u_x0)/total
  w_1 <- sum(u_x1)/total
  return (list('w_0' = w_0, 'w_1' = w_1))
}

find_site <- function(location_in_grid, n_easting, n_northing){
  number_of_sites <- n_easting*n_northing
  coords <- expand.grid(1:n_easting, n_northing:1)
  df_temp <- data.frame(sites = seq(1, number_of_sites),
                        easting = coords$Var1,
                        northing = coords$Var2)
  index_site <- filter(df_temp, easting == location_in_grid[1] & northing == location_in_grid[2])$sites
  return(index_site)
}


set_initial_weights <- function(indexes_for_voi){
  B_0 <- length(indexes_for_voi$index_0)
  B_1 <- length(indexes_for_voi$index_1)
  weights_0 <- rep(1/B_0, B_0)
  weights_1 <- rep(1/B_1, B_1)
  df_vs <- data.frame(v_value = c(weights_0, weights_1),
                      index_particle = c(indexes_for_voi$index_0, indexes_for_voi$index_1),
                      scenario = c(rep(0, B_0), rep(1, B_1)))
  return(df_vs)
}


new_set_initial_weights <- function(particles_set){
  indexes_for_voi <- particles_set$index
  scenario_for_voi <- particles_set$scenario
  
  B_0 <- (particles_set %>% filter(scenario ==0) %>% count(scenario))$n
  B_1 <- (particles_set %>% filter(scenario ==1) %>% count(scenario))$n
  weights_0 <- 1/B_0
  weights_1 <- 1/B_1
  df_vs <- data.frame(index = indexes_for_voi,
                      scenario = scenario_for_voi)
  df_vs$v_value <- ifelse(df_vs$scenario == 0, weights_0, weights_1)
  return(df_vs)
}


make_decision = function(t, w1){
  w0 = 1-w1
  Exp_v_stop=value_function(x=0,a=0,t)*(w0)+value_function(x=1,a=0,t)*w1
  Exp_v_cont=value_function(x=0,a=1,t)*(w0)+value_function(x=1,a=1,t)*w1
  res = ifelse(Exp_v_stop>Exp_v_cont,'stop','continue')
  return(res)
}


assimilate_data <- function(particles, sampled_data_y, w1, df_vs, inv_R){
    print(paste('To assimilate data we used w1 =', round(w1,3)))
    p <- ncol(particles) - 2
    probability_x <- w1 
    data_y = sampled_data_y[ , 1:p] 
    df_phi = new_gaussian_measurement_model(data_y, inv_R, particles) 
    df_phi_with_weights = add_particle_weights(df_phi, df_vs)
    df_u_values = new_compute_unnormalized_joint_probability(probability_x, df_phi_with_weights)
    df_v_values = new_compute_posterior_weights(df_u_values) 
    probability_x = new_classification_probabilities_given_y(df_u_values)$w_1
    return(list(w1 = probability_x, posterior_vs = df_v_values, u_values = df_u_values))
}

compute_data_alternative <- function(results_t, df_price){
  voi <- data.frame(results_t)
  price <- df_price
  df <- data.frame(t(rbind(voi, price)))
  colnames(df) <- c('voi', 'price')
  df$diff <- df$voi - df$price
  max_diff <- max(df$diff)
  print(df)
  if(max_diff <= 0){
    chosen_data <- 'none'
    difference_value <- 0
    voi_value <- NA
  } else{
    chosen_data <- rownames(df %>% filter(diff == max_diff))
    difference_value <- max_diff
    voi_value <- (df %>% filter(diff == max_diff))$voi}
  return(list(chosen_data = chosen_data, difference_value = difference_value, voi_value =  voi_value))
} 


update_from_survey_multiple_data_types <- function(indexes_for_voi, index_survey, data_vector, df_price, labels_array, parent_directory, number_of_cores, starting_time){
  start <- proc.time()
  years <- seq(from = starting_time, to = 23, by = 3)
  print(paste('Max. number of surveys:', length(years)))
  time_steps <- length(years)
  stop_time <- 0
  list_weights <- list()
  list_v_values <- list()
  list_ESS <- list()
  chosen_data_list <- list()
  D <- length(data_vector)
  data_alternative_list <- list()
  difference_value_list <- list()
  voi_value_list <- list()
  all_voi <- list() 
  for (step in 1:time_steps){
    particles_set <- list()
    survey_set <- list()
    results_t <- list()
    # start block
    t = years[step]
    for(d in 1:D){
      data_d <- data_vector[d]
      dataset <- load_realizations(data_d, t, labels_array, parent_directory)
      number_of_features <- ncol(dataset) - 2
      particles_set[[data_d]] <- new_select_realizations(dataset, indexes_for_voi) 
      survey_set[[data_d]] <- new_select_realizations(dataset, index_survey)
      if(step == 1){ # t==2
        w1 <- 0.25
        df_vs <- new_set_initial_weights(particles_set[[data_d]])
      }
      results_t[[data_d]] <-  get_voi(particles_set[[data_d]], w1, df_vs, inverse_of_matrix_R(data_type = data_d, number_of_features = number_of_features), t, number_of_cores)
    }
    all_voi[[step]] <- unlist(results_t) 
    print(paste('Current w1 =',round(w1,3)))
    data_alternative_list<- compute_data_alternative(results_t, df_price) # return 5 options: 4 data_types and 'none', a df is printed
    chosen_data_list[[step]] <- data_alternative_list$chosen_data
    difference_value_list[[step]] <- data_alternative_list$difference_value
    voi_value_list[[step]] <- data_alternative_list$voi_value
    if(difference_value_list[[step]]>0){
      data_s <- chosen_data_list[[step]]
      print(paste('We gather data:', data_s))
      print(paste('Index of data_y is index =', index_survey))
      sampled_data_y <- survey_set[[data_s]] %>% filter(index==index_survey)
      #print(sampled_data_y)
      #print(class(sampled_data_y))
      number_of_features <- ncol(sampled_data_y) - 2
      updated_values <- assimilate_data(particles_set[[data_s]], sampled_data_y, w1, df_vs, inverse_of_matrix_R(data_type = data_s, number_of_features = number_of_features))
      w1 <- updated_values$w1
      df_vs <- updated_values$posterior_vs
      chosen_alternative <- make_decision(t, w1)
      print(paste('Decision is at time t =', t, 'is',chosen_alternative))
      print(paste('Updated value of w1 =',round(w1,3)))
      list_weights[[step]] <- w1
      list_v_values[[step]] <- df_vs
      list_ESS[[step]] <- ESS(w1, df_vs)
      if(chosen_alternative == 'stop'){
        stop_time <- t
        print(paste('We choose to stop at time step t =',t))
        break
      }
    } else{print(paste('No updates at time t =',t))
      list_weights[[step]] <- w1
      list_v_values[[step]] <- df_vs
      list_ESS[[step]] <- ESS(w1, df_vs)
    }
    #
  }
  end <- proc.time()
  total_time <- (end-start)[3]
  return(list(w_values = list_weights, v_values = list_v_values, stop_time = stop_time, voi = voi_value_list, surveys = chosen_data_list, all_voi = all_voi, elapsed_time = total_time, ESS = list_ESS))
}


summarize_results <- function(list_global, indexes_for_survey, labels_array, surveys){
  n <- nrow(labels_array)
  df_indexes <- data.frame(scenario = labels_array, index = seq(1, n)) %>% filter(index %in% indexes_for_survey)
  m <- length(list_global)
  j <- 1:m
  vector_stop_time <- sapply(j, function(j) list_global[[j]]$stop_time)
  w_values <- sapply(j, function(j) unlist(list_global[[j]]$w_values))
  steps <- surveys
  matrix_w <- matrix(data = NA, nrow = m, ncol = steps) #10x8
  colnames(matrix_w) <- paste(rep('w', steps), 1:steps, sep = '_')
  surveys <- sapply(j, function(j) unlist(list_global[[j]]$surveys))
  matrix_surveys <- matrix(data = NA, nrow = m, ncol = steps) #10x8
  colnames(matrix_surveys) <- paste(rep('survey', steps), 1:steps, sep = '_')
  for(k in 1:m){
    matrix_w[k, 1:length(w_values[[k]])] = w_values[[k]]
    matrix_surveys[k, 1:length(surveys[[k]])] = surveys[[k]]
    }
  df <- cbind(df_indexes, round(matrix_w, 5), matrix_surveys, vector_stop_time)
  df$binary_stop = ifelse(df$vector_stop_time>0, 1, 0)
  df$cat <- 2*df$scenario-df$binary_stop + 1
  results <- df %>% count(cat)
  meanings = c('0: false_positive', '1: true_negative', '2: true_positive', '3: false_negative')
  results$prop <- round(results$n/m, 4)
  return(list(meanings = meanings, summary = results, df = df))
}

complement_indexes_balanced <- function(labels, list, seed, size_of_partition){
  index_voi <- sort(c(unlist(list$index_0), unlist(list$index_1)))
  df <- data.frame(scenario = labels_array,
                   index = seq(1, nrow(labels_array)))
  df_voi <- df %>% filter(index %in% c(index_voi))
  df_survey <- anti_join(df, df_voi, by = 'index')
  df_survey_0 <- df_survey %>% filter(scenario == 0)
  df_survey_1 <- df_survey %>% filter(scenario == 1)
  size_of_sample <- as.integer(size_of_partition/2) #50
  set.seed(seed)
  sample_ids_0 <- sample(df_survey_0$index, size_of_sample, replace = FALSE)
  set.seed(seed)
  sample_ids_1 <- sample(df_survey_1$index, size_of_sample, replace = FALSE)
  index_survey <- c(sort(sample_ids_0), sort(sample_ids_1))
  return(index_survey)
}


new_reweigh_particle_weights <- function(df_phi, df_vs, j, particles){
  df_merged <- merge(df_phi, df_vs, by=c('scenario', 'index')) # size 299
  index_data_y_j <- particles[j, ]$index 
  df <- df_merged %>% filter(!index %in% c(index_data_y_j)) # %>% arrange(scenario, index)
  df_0 <- df %>% filter(scenario == 0)
  df_1 <- df %>% filter(scenario == 1)
  sum_0 <- sum(df_0$v_value)
  sum_1 <- sum(df_1$v_value)
  total <- nrow(df_merged)
  new_v_values <- c()
  for(i in 1:total){
    new_v_values[i] <- ifelse(df$scenario[i] == 0, df$v_value[i]/sum_0, df$v_value[i]/sum_1)
  }
  df_ <- data.frame(phis = df$phis,
                    scenario = df$scenario,
                    index = df$index,
                    v_values = new_v_values)
  df_ <- df_ %>% arrange(scenario, index)
  return(df_)
}


add_particle_weights <- function(df_phi, df_vs){
  df_merged <- merge(df_phi, df_vs, by=c('scenario', 'index'))
  df_ <- data.frame(phis = df_merged$phis,
                    scenario = df_merged$scenario,
                    index = df_merged$index,
                    v_values = df_merged$v_value)
  df_ <- df_ %>% arrange(scenario, index)
  return(df_)
}

new_partition_indexes <- function(labels, seed, size_of_partition){
  df <- data.frame(scenario = labels_array,
                   index = seq(1, nrow(labels_array)))
  df_0 <- df %>% filter(scenario == 0)
  df_1 <- df %>% filter(scenario == 1) 
  size_of_sample <- as.integer(size_of_partition/2) #150
  set.seed(seed)
  sample_ids_0 <- sample(df_0$index, size_of_sample, replace = FALSE)
  set.seed(seed)
  sample_ids_1 <- sample(df_1$index, size_of_sample, replace = FALSE)
  vec <- sort(c(sample_ids_0, sample_ids_1))
  return(vec)
}

new_complement_indexes_balanced <- function(labels, vector_of_indexes, seed, size_of_partition){
  index_voi <- vector_of_indexes
  df <- data.frame(scenario = labels_array,
                   index = seq(1, nrow(labels_array)))
  df_voi <- df %>% filter(index %in% c(index_voi))
  df_survey <- anti_join(df, df_voi, by = 'index')
  df_survey_0 <- df_survey %>% filter(scenario == 0)
  df_survey_1 <- df_survey %>% filter(scenario == 1)
  size_of_sample <- as.integer(size_of_partition/2) 
  set.seed(seed)
  sample_ids_0 <- sample(df_survey_0$index, size_of_sample, replace = FALSE)
  set.seed(seed)
  sample_ids_1 <- sample(df_survey_1$index, size_of_sample, replace = FALSE)
  index_survey <- sort(c(sample_ids_0, sample_ids_1))
  return(index_survey)
}

new_select_realizations <- function(df, vector_of_indexes){
  df_ <- df %>% filter(index %in% vector_of_indexes) %>% arrange(scenario, index)
  return(df_)
}


get_labels_for_indexes <- function(labels_array, vector_of_indexes){
  total <- length(vector_of_indexes)
  df <- data.frame(scenario = labels_array, index = seq(1, nrow(labels_array))) %>% filter(index %in% vector_of_indexes)
  if(sum(df$index == vector_of_indexes) == total){res <- df$scenario}
  return(res)
}