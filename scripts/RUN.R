devtools::install_github('nathanvan/parallelsugar')


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

stime_baseline <- system.time(markov(n.cycles = n.cycles, n.samples = n.samples)) 
stime_reduceddim <- system.time(markov_reduced_dimensions(n.cycles = n.cycles, n.samples = n.samples))

## Vectorising
stime_vector_tx <- system.time(markov_vectorisetx(n.cycles = n.cycles, n.samples = n.samples))
stime_vector_txsmp <- system.time(markov_vectorisetxsmp(n.cycles = n.cycles, n.samples = n.samples))

stime_vector_smp <- system.time(markov_expanded_vectorisesmp(n.states = n.states, n.cycles = n.cycles, n.samples = n.samples))

## Parallelising 
stime_parallel_foreach <- system.time(markov_expanded_parallisesmp_foreach(n.states = n.states, n.cycles = n.cycles, n.samples = n.samples, n.cores = 2))
stime_parallel_furrr <- system.time(markov_expanded_parallisesmp_furrr(n.states = n.states, n.cycles = n.cycles, n.samples = n.samples))
stime_parallel_mclapply <- system.time(markov_expanded_parallisesmp_mclapply(n.states = n.states, n.cycles = n.cycles, n.samples = n.samples))
