# Smoking Cessation Markov model
# Howard Thom 14-June-2019

# Load necessary libraries
# If not installed use the following line first

#' Markov smoking model
#'
#' @param duration Numeric, the number of model cycles. Defaults to 10. This is a 
#' combination of time and cycle length (so 5 years with a cycle length of 6 months would give
#' n.cycles).
#' @param no_samples Numeric, the number of samples to take. Defaults to 10000. 
#' @param discount Numeric, defaults to 1.035. Discounting to apply to costs and qalys
#' @return Output
#' @export
#' @examples 
#' 
markov_seabbs <- function(no_samples = 10000, duration = 100, discount = 1.035) {

  future::plan(future::multiprocess())
  # Transitions -------------------------------------------------------------
  # 1. Specify transition matrices for each intervention
  # Baseline - Soc
  soc_transition <- function() {
    tmp <- rbind(VGAM::rdiric(1, c(88,12)),
                 VGAM::rdiric(1, c(8,92)))
    
    colnames(tmp) <- c("Smoking", "Not smoking")
    rownames(tmp) <- c("Smoking", "Not smoking")
    
    return(tmp)
  }
  
  # Intervention - Soc with website
  # Depends on Soc
  soc_with_website_transition <- function(baseline = NULL) {
    baseline[1, ] <- VGAM::rdiric(1,c(85,15))
    
    return(baseline)
  }
  
  
  ## Test
  #soc_trans_sample <- soc_transition()
 # soc_trans_sample
  
  #soc_with_website_trans_sample <- soc_with_website_transition(soc_trans_sample)
 # soc_with_website_trans_sample
  
  #Set up transition list
  transitions_list <- list(soc_transition, 
                           soc_with_website_transition)
  
  names(transitions_list) <- c("SoC", "Soc with Website")
  
  # Qualies -----------------------------------------------------------------
  # 2. Specify qaly costs per intervention (random sampling)
  
  
  qalys <- function() {
    qaly <- function(samples = 1) {
      smoking <- rnorm(1, mean = 0.95,sd = 0.01) / 2
      not_smoking <- 1 / 2
      
      out <- c(smoking, not_smoking)
      names(out) <- c("Smoking", "Not smoking")
      
      return(out)
      
    }
    
    soc <- qaly()
    soc_with_website <- soc
    
    out <- list(soc, soc_with_website)
    names(out) <- list("SoC", "Soc with Website")
    
    return(out)
  }
  
 # qalys()
  
  
  # Costs -------------------------------------------------------------------
  # 3. Specify costs per intervention (random sampling)
  
  intervention_costs <- function(samples = 1) {
    soc <- 0
    soc_with_website <- 50
    
    out <- c(soc, soc_with_website)
    names(out) <- c("SoC", "Soc with Website")
    
    return(out)
  }
  
 # intervention_costs()
  
  state_costs <- function(samples = 1) {
    state_cost <- function(samples = 1) {
      smoking <- 0
      not_smoking <- 0
      
      out <- c(smoking, not_smoking)
      names(out) <- c("Smoking", "Not smoking")
      
      return(out)
      
    }
    
    soc <- state_cost()
    soc_with_website <- soc
    
    out <- list(soc, soc_with_website)
    names(out) <- list("SoC", "Soc with Website")
    return(out)
  }
  
 # state_costs()
  
  
  
  # Cohort ------------------------------------------------------------------
  #4. Define cohort
  
  cohorts <- function() {
    
    cohort <- function() {
      smoking <- 1
      not_smoking <- 0
      
      out <- matrix(c(smoking, not_smoking), ncol = 2)
      colnames(out) <- c("Smoking", "Not smoking")
      
      return(out)
    }
    
    soc <- cohort()
    soc_with_website <- soc
    
    out <- list(soc, soc_with_website)
    names(out) <- list("SoC", "Soc with Website")
    
    return(out)
  }
  
  #cohorts()
  
  # Set up single sample -----------------------------------------------------
  #5. Set up single sample
  
  ## Set up sampling function
  single_sample <- function(transitions = NULL, state_costs = NULL, 
                            intervention_costs = NULL, cohorts = NULL, 
                            qalys = NULL) {
    
    #sample baseline transition matrix
    baseline <- transitions[[1]]()
    
    #sample all interventions depending on baseline
    interventions <- purrr::map(2:length(transitions), ~ transitions[[.]](baseline))
    
    #update transitions as a single sample
    transitions[[1]] <- baseline
    transitions[-1] <- interventions
    
    sample <- tibble::tibble(intervention = names(transitions), 
                             transition = transitions, 
                             state_cost = state_costs(), 
                             intervention_cost = intervention_costs(),
                             cohort = cohorts(),
                             qalys = qalys()
    )
    
  }
  
  # List of all transitions (named), all subsequent transitions from the first must take the first transition as
  # an inp
  
  #test_sample <- single_sample(transitions = transitions_list,
 #                              state_costs = state_costs,
 #                              intervention_costs = intervention_costs,
 #                              cohorts = cohorts, qalys = qalys)
  
  
  
  # Run the model ----------------------------------------------------------
  
  run_markov <- function(transition = NULL, cohort = NULL, state_cost = NULL, 
                         intervention_cost = NULL, qalys = NULL, duration = NULL,
                         discount = NULL) {
    
    ## Preallocate
    sim <- matrix(NA, nrow = duration, ncol = nrow(transition))
    colnames(sim) <- colnames(transition)
    
    ## Assign initial pop
    sim[1, ] <- cohort
    
    ##Loop over the rest of the model 
    for (i in 2:duration) {
      sim[i, ] <- sim[i - 1, ] %*% transition
    }
    
    ## Discounting
    discounting <-  (1 / 1.035)^(0:(duration - 1))
    
    ##Total costs per cycle
    total_costs_cycle <- (sim %*% state_cost) * discounting
    
    ##Total QALYs per cycle
    discounted_qalys <- (sim %*% qalys) * discounting
    total_qalys <- sum(discounted_qalys)
    
    ## Overall costs
    total_costs <- sum(total_costs_cycle) + intervention_cost
    
    out <- tibble::tibble(total_costs = total_costs, total_qalys = total_qalys)
    
    return(out)
  }

  
  ## Test model run
#  test_sim <- run_markov(transition = test_sample$transition[[1]],
 #                        cohort = test_sample$cohort[[1]],
 #                        state_cost = test_sample$state_cost[[1]], 
 #                        intervention_cost = test_sample$intervention_cost[[1]], 
 #                        qalys = test_sample$qalys[[1]], 
 #                        duration = 10,
 #                        discount = 1.035)
  
  
  
  
  # Analyse samples --------------------------------------------------------
  
  analyse_model <- function(results) {
    
    ## Work out incremental costs
    
    ## Work out mean costs
    mean_costs <- results %>% 
      dplyr::group_by(intervention) %>% 
      dplyr::summarise(
        mean_costs = mean(total_costs),
        mean_qalys = mean(total_qalys)
      ) %>% 
      dplyr::ungroup()
  }
  
  
  # Generate samples --------------------------------------------------------
  
  samples <- furrr::future_map_dfr(1:no_samples, ~ single_sample(transitions = transitions_list,
                                                          state_costs = state_costs,
                                                          intervention_costs = intervention_costs,
                                                          cohorts = cohorts, qalys = qalys), .id = "sample")
  
  results <- furrr::future_map_dfr(1:nrow(samples), 
                            ~ run_markov(transition = samples$transition[[.]],
                                         cohort = samples$cohort[[.]],
                                         state_cost = samples$state_cost[[.]], 
                                         intervention_cost = samples$intervention_cost[[.]], 
                                         qalys = samples$qalys[[.]], 
                                         duration = duration,
                                         discount = discount))
  
  combined <- dplyr::bind_cols(samples, results)
  
  # Analyse model -----------------------------------------------------------
  
  sum <- analyse_model(combined)
  
}
