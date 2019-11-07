library(VGAM)
library(tidyverse)
library(parallel)
library(foreach)

files.sources <- file.path(getwd(), "R", list.files(file.path(getwd(), "R")))
sapply(files.sources, source)

n.states <- 10
n.cycles <- 100
n.samples <- 10000
n.cores <- parallel::detectCores() - 2

s1 <- system.time(markov(n.cycles = n.cycles, n.samples = n.samples)) 
s2 <- system.time(markov_reduced_dimensions(n.cycles = n.cycles, n.samples = n.samples))

## Vectorising
s3 <- system.time(markov_vectorisetx(n.cycles = n.cycles, n.samples = n.samples))
s5 <- system.time(markov_vectorisetxsmp(n.cycles = n.cycles, n.samples = n.samples))

s4 <- system.time(markov_expanded_vectorisesmp(n.states = n.states, n.cycles = n.cycles, n.samples = n.samples))

## Parallelising 
s6 <- system.time(markov_expanded_parallisesmp_foreach(n.states = n.states, n.cycles = n.cycles, n.samples = n.samples, n.cores = 2))
s7 <- system.time(markov_expanded_parallisesmp_furrr(n.states = n.states, n.cycles = n.cycles, n.samples = n.samples))

