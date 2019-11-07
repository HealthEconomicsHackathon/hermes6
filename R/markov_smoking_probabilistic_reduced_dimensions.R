# Smoking Cessation Markov model
# Howard Thom 14-June-2019

# Load necessary libraries
# If not installed use the following line first
# install.packages("VGAM")

#' Reduced dimensions in markov smoking probabilistic model
#'
#' @return Output
#' @export
markov_reduced_dimensions <- function() {
  set.seed(14143)
  
  # Define the number and names of treatments
  # These are Standard of Care with website
  # and Standard of Care without website
  n.treatments<-2
  treatment.names<-c("SoC with website","SoC")
  
  # Define the number and names of states of the model
  # This is two and they are "Smoking" and "Not smoking"
  n.states<-2
  state.names<-c("Smoking","Not smoking")
  
  # Define the number of cycles
  # This is 10 as the time horizon is 5 years and cycle length is 6 months
  # The code will work for any even n.cycles (need to change the discounting code if
  # an odd number of cycles is desired)
  
  n.cycles<-100
  
  # Define simulation parameters
  # This is the number of PSA samples to use
  n.samples<-10000
  
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
  state.qalys[,"Smoking"]<-stats::rnorm(n.samples,mean=0.95,sd=0.01)/2
  
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
  
  #i.treatment <- 1
  #i.sample <- 1
  #i.cycle <- 2
  
  disc_vec <- (1/1.035)^rep(c(0:(n.cycles/2-1)),each=2)
  
  # The remainder of the cohort.vectors will be filled in by Markov updating below
  
  # Main model code
  # Loop over the treatment options
  
  for(i.treatment in 1:n.treatments)
  {
    transition.matrices_tr <- transition.matrices[i.treatment,,,]
    
    # Loop over the PSA samples
    for(i.sample in 1:n.samples)
    {
      
      transition.matrices_tr_sample <- transition.matrices_tr[i.sample,,]
      
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n.cycles
      for(i.cycle in 2:n.cycles)
      {
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i.e. pi_j = pi_(j-1)*P
        cohort.vectors[i.treatment, i.sample,i.cycle,]<-
          cohort.vectors[i.treatment, i.sample,i.cycle-1,] %*%
          transition.matrices_tr_sample
      }
      
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
markov_reduced_dimensions()
# ICER 746.4855 (benchmark for future changes)

# Test Time for the basic function
system.time(markov_reduced_dimensions())
#  
# user  system elapsed 
# 7.16    0.03    7.41

# $average.costs
# SoC with website              SoC 
# 50                0 

# $average.effects
# SoC with website              SoC 
# 23.78637         23.71939 



###########################################################
# try try
###########################################################

# Function to run the analysis based on treatment arms (scenarios)
markov_reduced_dimensions2 <- function(scenario = "UoC_Web") {
  
# Define cycles and samples (simulation times)
n.cycles<-100
n.samples<-10000

# function of define probabilites that depend on treatment arms
treatment_arm <- function(treatment = "UoC") {
  # Probabilites of remaining in certain states
  if(treatment == "UoC") {
    p.remain_smoke <- VGAM::rdiric(n.samples,c(88,12))[, 1]
    cost.oneoff <-rep(0, times=n.samples)
    }
  if(treatment == "UoC_Web") {
    p.remain_smoke <- VGAM::rdiric(n.samples,c(85,15))[, 1]
    cost.oneoff <-rep(50, times=n.samples)
    }
  temp_list <- list("p.remain_smoke" = p.remain_smoke, "cost.oneoff" = cost.oneoff)
  return(temp_list) 
}

######################################################
# Define all variables (objects) that have 
# uncertanties based on treatment arms
#######################################################
#####!!!!!!!!!!!!!!!!!!!
p.remain_smoke <- treatment_arm (treatment = scenario) [["p.remain_smoke"]]
cost.oneoff <- treatment_arm (treatment = scenario) [["cost.oneoff"]]

p.remain_nonsmoke <-VGAM::rdiric(n.samples,c(8,92))[, 1] 


# QALY associated with 1-year in the smoking state is Normal(mean=0.95, SD=0.01)
# Divide by 2 as cycle length is 6 months
qaly.smoke <-rnorm(n.samples,mean=0.95,sd=0.01)/2
qaly.notsmoke <- rep(1/2, times=n.samples)

# State Cost
cost.smoke <-rep(0, times=n.samples)
cost.notsmoke <- rep(0, times=n.samples)

#Discount
disc_vec <- (1/1.035)^rep(c(0:(n.cycles/2-1)),each=2)


# Create matrices for probability of staying in each states cross each cycle (colomns) and simulations (rows)
# Number of matrices created = number of states
m.notsmoke <- m.smoke <- matrix(0, nrow =  n.samples, ncol = n.cycles + 1, # starts from 0 and run n.cycle, so needs n.cycle + 1 cycles
                                dimnames = list(1:n.samples, 0:n.cycles))  #initialize state matrices
# at cycle 0, all people start from smoking in all treatment arms, none from not-smoking (define column 0 for base population)
# matrices for smoking state
m.smoke [, 1] <- 1       
# matrices for smoking state  
m.notsmoke [, 1] <- 0


# Run cycles for all simulations at once 
for (cycle in 1:n.cycles) {
    m.smoke[, cycle + 1] <-  m.smoke[, cycle] * p.remain_smoke +
                             m.notsmoke [,cycle]* (1-p.remain_nonsmoke)  # calculate the prop of smoke at time t + 1
    m.notsmoke[,cycle + 1] <-  m.notsmoke[, cycle] * p.remain_nonsmoke +
                               m.smoke[, cycle]* (1-p.remain_smoke)
}  

# Distcount and Sum up

TC <- TE <- vector("numeric", length = n.samples) 
TE <- (qaly.smoke * m.smoke[,-1] + qaly.notsmoke* m.notsmoke[,-(n.cycles)]) %*% disc_vec   # total discounted cost for all simulations
TC <- (cost.smoke * m.smoke[,-1] + cost.notsmoke* m.notsmoke[,-(n.cycles)]) %*% disc_vec + cost.oneoff    # total discounted QALYs for all simulations


# Mean of QALY and Costs and output
results <- list("TE" = mean(TE), "TC" = mean(TC))
return(results) 
}

# return list of results
results_UoC <- markov_reduced_dimensions2 (scenario = "UoC")
results_UoC_Web <- markov_reduced_dimensions2 (scenario = "UoC_Web")

ICER <- (results_UoC_Web [["TC"]] - results_UoC [["TC"]]) / (results_UoC_Web [["TE"]] - results_UoC [["TE"]])


results_UoC [["TE"]]

system.time(markov_reduced_dimensions2()) 

#    user  system elapsed 
#    0.1     0.0     0.1



