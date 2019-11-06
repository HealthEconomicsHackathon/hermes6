source("markov_smoking_probabilistic_original.R")
source("markov_smoking_probabilistic_reduced_dimensions.R")

system.time({markov()})
system.time({markov_reduced_dimensions()})
