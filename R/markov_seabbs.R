library(tibble)


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
sco_with_website_transition <- function(baseline) {
  baseline[1, ] <- VGAM::rdiric(n.samples,c(85,15))
  
  return(baseline)
}


## Test
soc_trans_sample <- soc_transition()
soc_trans_sample

soc_with_website_trans_sample <- sco_with_website_transition(soc_trans_sample)
soc_with_website_trans_sample


# Qualies -----------------------------------------------------------------
# 2. Specify qualy costs per intervention (random sampling)


qualy <- function(samples = 1) {
  smoking <- rnorm(1, mean = 0.95,sd = 0.01) / 2
  not_smoking <- 1 / 2
  
  out <- c(smoking, not_smoking)
  names(out) <- c("Smoking", "Not smoking")
  
  return(out)
  
}

qualy()


# Costs -------------------------------------------------------------------
# 3. Specify costs per intervention (random sampling)

costs <- function(samples = 1) {
  soc <- 0
  soc_with_website <- 50
  
  out <- c(soc, soc_with_website)
  names(out) <- c("SoC", "Soc with Website")
  
  return(out)
}

costs()



# Cohort ------------------------------------------------------------------
#4. Define cohort

cohort <- function() {
  smoking <- 1
  not_smoking <- 0
  
  out <- c(smoking, not_smoking)
  names(out) <- c("Smoking", "Not smoking")
  
  return(out)
}

cohort()


# Set up and sampling -----------------------------------------------------
#5. Set up sampling

# Smoking Cessation Markov model
# Howard Thom 14-June-2019

# Load necessary libraries
# If not installed use the following line first


#' Markov smoking model
#'
#' @param n.cycles Numeric, the number of model cycles. Defaults to 10. This is a 
#' combination of time and cycle length (so 5 years with a cycle length of 6 months would give
#' n.cycles).
#' @param n.samples Numeric, the number of samples to take. Defaults to 10000. 
#' @return Output
#' @export
#'
markov_seabbs <- function(n.cycles = 100, n.samples = 10000) {
  set.seed(14143)
  

  #############################################################################
  ## Input parameters #########################################################
  #############################################################################
  
  # The transition matrix is a 2x2 matrix
  # Rows sum to 1
  # Top left entry is transition probability from smoking to smoking
  # Top right is transition probability from smoking to not smoking
  # Bottom left is transition probability from not smoking to smoking
  # Bottom right is transition probability from not smoking to not smoking
  
  # There is one transition matrix for each reatment option and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  transition.matrices<-array(dim=c(n.treatments,n.samples,n.states,n.states),
                             dimnames=list(treatment.names,NULL,state.names,state.names))
  
  # First the transition matrix for Standard of Care with website
  # Transitions from smoking 
  transition.matrices["SoC with website",,"Smoking",]<-VGAM::rdiric(n.samples,c(85,15))
  # Transitions from not smoking
  transition.matrices["SoC with website",,"Not smoking",]<-VGAM::rdiric(n.samples,c(8,92))
  
  # Second the transition matrix for Standard of Care
  # Transitions from smoking 
  transition.matrices["SoC",,"Smoking",]<-VGAM::rdiric(n.samples,c(88,12))
  # Transitions from not smoking
  # These should be the same as the transition probabilities from not smoking for SoC with website
  # as the website has no impact on probability of relapse
  transition.matrices["SoC",,"Not smoking",]<-transition.matrices["SoC with website",,"Not smoking",]
  
  # Now define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state.qalys<-array(dim=c(n.samples, n.states),dimnames=list(NULL,state.names))
  
  # QALY associated with 1-year in the smoking state is Normal(mean=0.95, SD=0.01)
  # Divide by 2 as cycle length is 6 months
  state.qalys[,"Smoking"]<-rnorm(n.samples,mean=0.95,sd=0.01)/2
  
  # QALY associated with 1-year in the not smoking state is 1 (no uncertainty)
  # So all PSA samples have the same value
  # Again divide by 2 as cycle length is 6 months
  state.qalys[,"Not smoking"]<-1/2
  
  # And finally define the state costs
  # These are all zero as the only cost is a one-off subscription fee of ?50
  # to the smoking cessation website
  state.costs<-array(0,dim=c(n.samples, n.states),dimnames=list(NULL,state.names))
  
  # Define the treatment costs
  # One for each PSA sample and each treatment
  # Treatment costs are actually fixed but this allows flexibility if we
  # want to include uncertainty/randomness in the cost
  treatment.costs<-array(dim=c(n.treatments,n.samples),dimnames=list(treatment.names,NULL))
  
  # Cost of the smoking cessation website is a one-off subscription fee of ?50
  treatment.costs["SoC with website",]<-50
  # Zero cost for standard of care
  treatment.costs["SoC",]<-0
  
  #############################################################################
  ## Simulation ###############################################################
  #############################################################################
  
  # Build an array to store the cohort vector at each cycle
  # Each cohort vector has 2 (=n.states) elements: probability of being in smoking state,
  # and probability of being in the not smoking state
  # There is one cohort vector for each treatment, for each PSA sample, for each cycle.
  cohort.vectors<-array(dim=c(n.treatments,n.samples,n.cycles,n.states),
                        dimnames=list(treatment.names,NULL,NULL,state.names))
  
  # Assume that everyone starts in the smoking state no matter the treatment
  cohort.vectors[,,1,"Smoking"]<-1
  cohort.vectors[,,1,"Not smoking"]<-0
  
  # Build an array to store the costs and QALYs accrued per cycle
  # One for each treatment, for each PSA sample, for each cycle
  # These will be filled in below in the main model code
  # Then discounted and summed to contribute to total costs and total QALYs
  cycle.costs<-array(dim=c(n.treatments,n.samples,n.cycles),
                     dimnames=list(treatment.names,NULL,NULL))
  cycle.qalys<-array(dim=c(n.treatments,n.samples,n.cycles),
                     dimnames=list(treatment.names,NULL,NULL))
  
  # Build arrays to store the total costs and total QALYs
  # There is one for each treatment and each PSA sample
  # These are filled in below using cycle.costs, 
  # treatment.costs, and cycle.qalys
  total.costs<-array(dim=c(n.treatments,n.samples),
                     dimnames=list(treatment.names,NULL))
  total.qalys<-array(dim=c(n.treatments,n.samples),
                     dimnames=list(treatment.names,NULL))
  
  
  # The remainder of the cohort.vectors will be filled in by Markov updating below
  
  # Main model code
  # Loop over the treatment options
  for(i.treatment in 1:n.treatments)
  {
    # Loop over the PSA samples
    for(i.sample in 1:n.samples)
    {
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n.cycles
      for(i.cycle in 2:n.cycles)
      {
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i.e. pi_j = pi_(j-1)*P
        cohort.vectors[i.treatment,i.sample,i.cycle,]<-
          cohort.vectors[i.treatment,i.sample,i.cycle-1,]%*%
          transition.matrices[i.treatment,i.sample,,]
      }
      
      # Now use the cohort vectors to calculate the 
      # total costs for each cycle
      cycle.costs[i.treatment,i.sample,]<-
        cohort.vectors[i.treatment,i.sample,,]%*%state.costs[i.sample,]
      # And total QALYs for each cycle
      cycle.qalys[i.treatment,i.sample,]<-
        cohort.vectors[i.treatment,i.sample,,]%*%state.qalys[i.sample,]
      
      # Combine the cycle.costs and treatment.costs to get total costs
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total.costs[i.treatment,i.sample]<-treatment.costs[i.treatment,i.sample]+
        cycle.costs[i.treatment,i.sample,]%*%
        (1/1.035)^rep(c(0:(n.cycles/2-1)),each=2)
      
      # Combine the cycle.qalys to get total qalys
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total.qalys[i.treatment,i.sample]<-cycle.qalys[i.treatment,i.sample,]%*%
        (1/1.035)^rep(c(0:(n.cycles/2-1)),each=2)
    }
  }
  
  #############################################################################
  ## Analysis of results ######################################################
  #############################################################################
  output <- list()
  # Average costs
  # These are ?50 on the website and 0 on standard of care as there are no
  # costs other than the website subscription cost
  output$average.costs<-rowMeans(total.costs)
  # Average effects (in QALY units)
  # These are slightly higher on the website as higher probability of 
  # quitting smoking
  output$average.effects<-rowMeans(total.qalys)
  
  # Incremental costs and effects relative to standard of care
  # No uncertainty in the costs as the website cost is fixed at ?50
  output$incremental.costs<-total.costs["SoC with website",]-total.costs["SoC",]
  # In some samples the website leads to higher QALYs but in others it is negative
  # There is uncertainty as to whether the website is an improvement over SoC
  output$incremental.effects<-total.qalys["SoC with website",]-total.qalys["SoC",]
  
  # The ICER comparing Standard of care with website to standard of care
  # This is much lower than the ?20,000 willingness-to-pay threshold indicating
  # good value for money
  output$ICER<-mean(output$incremental.costs)/mean(output$incremental.effects)
  
  # Incremental net benefit at the ?20,000 willingness-to-pay
  # Sometimes positive (website more cost-effective) and sometimes negative (SoC more cost-effective)
  # Need to look at averages and consider probabilities of cost-effectiveness
  output$incremental.net.benefit<-20000*output$incremental.effects-output$incremental.costs
  
  # Average incremental net benefit
  # This is positive indicating cost-effectiveness at the ?20,000 threshold
  output$average.inb<-mean(output$incremental.net.benefit)
  
  # Probability cost-effective
  # This is the proportion of samples for which the incremental net benefit is positive
  # It is clost to 72%, representing good degree of certainty
  # in recommendation to adopt the smoking cessation website
  output$probability.cost.effective<-sum(output$incremental.net.benefit>0)/n.samples
  
  # Now use the BCEA package to analyse the results...
  output
}
