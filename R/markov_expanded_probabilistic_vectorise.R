# Smoking Cessation Markov model
# Howard Thom 14-June-2019

# Edited 7-November-2019 to use 10 health states instead of 2.

# Load necessary libraries
# If not installed use the following line first
# install.packages("VGAM")
 
#' Reduced dimensions in markov smoking probabilistic model
#'
#' @return Output
#' @export
markov_expanded_vectorise <- function() {
  set.seed(14143)
  
  # Define the number and names of treatments
  # These are Standard of Care with website
  # and Standard of Care without website
  n.treatments<-2
  treatment.names<-c("SoC with website","SoC")
  
  # Define the number and names of states of the model
  # Any number is allowed and states are named "State 1", "State 2", etc.
  n.states<-10
  state.names<-paste("State",c(1:n.states))
  
  # Define the number of cycles
  # This is 10 as the time horizon is 5 years and cycle length is 6 months
  # The code will work for any even n.cycles (need to change the discounting code if
  # an odd number of cycles is desired)
  
  n.cycles<-100
  
  # Define simulation parameters
  # This is the number of PSA samples to use
  n.samples<-25000
  
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
  
  # Now define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state.qalys<-array(dim=c(n.samples, n.states),dimnames=list(NULL,state.names))
  
  # And finally define the state costs
  # There is one for each PSA sample and each state
  # Store in an NA array and then fill in below
  state.costs<-array(dim=c(n.samples, n.states),dimnames=list(NULL,state.names))
  
  
  # Define transition matrices, state utilities and costs 
  for(i.state in 1:n.states)
  {
    # Use dirichlet distributions to define transition matrices
    # This is for illustration only; in practice this would be estimated on some data
    transition.matrices["SoC",,paste("State",c(i.state)),]<-VGAM::rdiric(n.samples,sample(c(50:100),n.states))
    transition.matrices["SoC with website",,paste("State",i.state),]<-VGAM::rdiric(n.samples,sample(c(50:100),n.states))
    
    # State utilities
    # Anything between 0 and 1
    # Divide by 2 as cycle length is 6 months
    state.qalys[,paste("State",i.state)]<-runif(n.samples,min=0,max=1)/2
    
    # State costs
    # Assumed normal with sd small enough to avoid negative values
    state.costs[,paste("State",i.state)]<-rnorm(n.samples,mean=100,sd=10)
    
  }
  
  
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
  # Each cohort vector has n.states elements: probability of being in each state,
  # There is one cohort vector for each treatment, for each PSA sample, for each cycle.
  cohort.vectors<-array(0,dim=c(n.treatments,n.samples,n.cycles,n.states),
                        dimnames=list(treatment.names,NULL,NULL,state.names))
  
  # Assume that everyone starts in first state
  cohort.vectors[,,1,"State 1"]<-1
  #cohort.vectors[,,1,paste("State",c(2:10))]<-0
  
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
  
  #i.treatment <- 1
  #i.sample <- 1
  #i.cycle <- 2
  
  disc_vec <- (1/1.035)^rep(c(0:(n.cycles/2-1)),each=2)
  
  # The remainder of the cohort.vectors will be filled in by Markov updating below
  
  
  
  # Main model code
  # Loop over the treatment options
  
  for(i.treatment in 1:n.treatments)
  {
    # Pre-extract treatment matrix
    transition.matrices_tr <- transition.matrices[i.treatment,,,]
    
    
     # transition.matrices_tr_sample <- transition.matrices_tr[i.sample,,]
      
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n.cycles
      for(i.cycle in 2:n.cycles)
      {
        # Avoiding matrix multiplication as it cannot be vectorised.
        # Instead looping over states and vectorising over samples
        for(i.state in 1:n.states)
        {
            # Occupancy probability for each state
            # Was initialised to 0
            cohort.vectors[i.treatment, ,i.cycle,i.state]<- 
              rowSums(cohort.vectors[i.treatment, ,i.cycle-1, ] * transition.matrices_tr[,,i.state])
        }
    }
    for(i.sample in 1:n.samples)
    {
      cohort.vectors_tr_sample <- cohort.vectors[i.treatment,i.sample,,]
      
      # Now use the cohort vectors to calculate the 
      # total costs for each cycle
      cycle.costs[i.treatment,i.sample,]<-
        cohort.vectors_tr_sample%*%state.costs[i.sample,]
      # And total QALYs for each cycle
      cycle.qalys[i.treatment,i.sample,]<-
        cohort.vectors_tr_sample%*%state.qalys[i.sample,]
      
      # Combine the cycle.costs and treatment.costs to get total costs
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total.costs[i.treatment,i.sample]<-treatment.costs[i.treatment,i.sample] +
        cycle.costs[i.treatment,i.sample,]%*%
        disc_vec
      
      # Combine the cycle.qalys to get total qalys
      # Apply the discount factor 
      # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total.qalys[i.treatment,i.sample]<-cycle.qalys[i.treatment,i.sample,]%*%
        disc_vec
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

system.time(output<-markov_expanded_vectorise())
