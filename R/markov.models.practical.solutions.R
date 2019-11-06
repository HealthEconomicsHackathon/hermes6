# Solution to Markov models practical
# Howard Thom 8-July-2019

###############################################################################################################
###############################################################################################################
###############################################################################################################

# Exercise 1.  Understanding the Markov update loop
# Run the file "markov.smoking.probabilistic.X.R" (version number is X, use the latest available)
# Can do this using the source() command
source("markov.smoking.probabilistic.2.R")

# a)
# The dimensions of the array 
# These are n.treatments (number of treatments are 2), n.samples (the number of PSA samples), 
# n.states (the number of states), and n.states (the number of states)
dim(transition.matrices)
# Look at 1st PSA sample of the transition matrix for SoC with website
# First row is probability of moving from smoking to smoking (i.e. staying in smoking state)
# and to not smoking. Second row is probability of moving from not smoking to
# smoking and to not smoking (i.e. staying in the not smoking state)
transition.matrices["SoC with website",1,,]

# Similarly the transition probability matrix for SoC
transition.matrices["SoC",1,,]

# b) 
# The first PSA sample of transition probabilities from Smoking on SoC with website are given by
transition.matrices["SoC with website",1,"Smoking",]
# The average transition probabilities from Smoking state over all PSA samples for SoC with website are given by
colMeans(transition.matrices["SoC with website",,"Smoking",])
# The average transition probabilities from Smoking state over all PSA samples for SoC are given by
# Higher probability of remaining in the smoking state if on SoC than on SoC with website
colMeans(transition.matrices["SoC",,"Smoking",])
# The transition probabilities from not smoking are the same for SoC with website and SoC alone
colMeans(transition.matrices["SoC with website",,"Not smoking",])
colMeans(transition.matrices["SoC",,"Not smoking",])

# c) Investigating the Markov update step
# Proportion on SoC with wesbite in first PSA sample starting in smoking and not smoking states is given by
cohort.vectors["SoC with website",1,1,]
# The first PSA transition matrix for SoC with website is given by
transition.matrices["SoC with website",1,,]
# We can therefore perform Matrix multiplication to see the proportions in these states on the second cycle
cohort.vectors["SoC with website",1,1,]%*%transition.matrices["SoC with website",1,,]
# Confirm that this matches the second cycle for SoC with website for the first PSA sample
# It does (ot at least should!)
cohort.vectors["SoC with website",1,2,]

###############################################################################################################
###############################################################################################################
###############################################################################################################

# Exercise 2. Using BCEA to analyse the results
# The total costs and QALYs are n.treatments x n.samples but BCEA expects n.samples x n.treatments
# Will therefore need to use the t() transpose function to convert to format expected by BCEA.
dim(total.costs)
dim(total.qalys)

# a)
# Build a BCEA object using the transposed QALYs and costs.
smoking.bcea<-bcea(e=t(total.qalys),c=t(total.costs),ref=1, interventions=treatment.names)

# b)
# The EIB is the expected incremental net benefit of SoC with wesbite compared to SoC. As it is positive
# we favour SoC with website at the supplied willingness-to-pay threshold of £20,000.
# CEAC is the probability that SoC with website has a higher net beneft, which is >70% so good certainty SoC
# with website is the best.
# And the ICER is about £3500 so also confident SoC with website is good value for money at £20,000 willingness-to-pay
summary(smoking.bcea,wtp=20000)

# c)
# The CE-plane looks very extreme with only very small differences in effectiveness and a constant difference in costs
ceplane.plot(smoking.bcea, wtp=20000)
# Is this correct? Yes, as the costs of the website are fixed at £50
rowMeans(total.costs)
# And there is only a small difference in effects
rowMeans(total.qalys)

# d)
# The CEAC plot is the probability that SoC with website (the reference) has the greatest net benefit
# over a range of willingness-to-pay thresholds.
# We see that this converges to around 75% as the willingness-to-pay increases
ceac.plot(smoking.bcea)

###############################################################################################################
###############################################################################################################
###############################################################################################################

# Exercise 3. Using the info.rank() function from BCEA
# This generates the proportion of total EVPI (value of information, or the monetary
# value of decision uncertainty)
# Construct the list with a matrix of input parameters and a set of names of the parameters
input.parameters<-list()
# This consist of 6 transition probabilies and the state qalys for smoking
# State QALYs for not smoking and the costs for smoking and not smoking are fixed so not included in the sensitivity analysis.
input.parameters$parameters<-c("Probs from smoking to smoking SoC with website","Probs from smoking to not smoking SoC with website",
                               "Probs from smoking to smoking SoC","Probs from smoking to not smoking SoC",
                               "Probs from not smoking to smoking","Probs from not smoking to not smoking","QALYs smoking")
# Build a matrix to store each of the uncertain input parameters
# Again 6 transition probabilities and the state qalys for smoking
input.parameters$mat<-matrix(NA, nrow=n.samples, ncol=7)
# Name the columns after the parameters
colnames(input.parameters$mat)<-input.parameters$parameters
# Fill in the matrix with the parameter samples used in the Markov model
# Transition parameters from Not smoking are the same for SoC as for SoC with website so only need it once.
input.parameters$mat[,"Probs from smoking to smoking SoC with website"]<-transition.matrices["SoC with website",,"Smoking","Smoking"]
input.parameters$mat[,"Probs from smoking to not smoking SoC with website"]<-transition.matrices["SoC with website",,"Smoking","Not smoking"]
input.parameters$mat[,"Probs from smoking to smoking SoC"]<-transition.matrices["SoC",,"Smoking","Smoking"]
input.parameters$mat[,"Probs from smoking to not smoking SoC"]<-transition.matrices["SoC",,"Smoking","Not smoking"]
input.parameters$mat[,"Probs from not smoking to smoking"]<-transition.matrices["SoC with website",,"Not smoking","Smoking"]
input.parameters$mat[,"Probs from not smoking to not smoking"]<-transition.matrices["SoC with website",,"Not smoking","Not smoking"]
# Only the state QALYs for smoking are uncertain
input.parameters$mat[,"QALYs smoking"]<-state.qalys[,"Smoking"]

# Finally, use the info.rank() function to generate the Bayesian decision analysis
# equivalent of deterministic sensitivity analyses
# We see that the transition probabilities that vary by intervention are most influential.
# The QALYs associated with smoking have moderate influence but much less than the transition probabilities.
# The transition probabilities that are the same for both interventions (i.e. those from the not smoking state)
# have limited impact on the decision
info.rank(input.parameters$parameters,input.parameters$mat,smoking.bcea)


