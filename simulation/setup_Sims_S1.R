library("MASS")
library("PROreg")
library("splines")
library("statmod")
library("rstan")

phi <- 0.05
alpha <- 0

source("simulation/runSimsS1_lp.R")
saveRDS(out, file="simulation/results/S1/S1_p005a0.rds")

phi <- 0.05
alpha <- 2

source("simulation/runSimsS1_lp.R")
saveRDS(out, file="simulation/results/S1/S1_p005a2.rds")


phi <- 0.05
alpha <- 4

source("simulation/runSimsS1_lp.R")
saveRDS(out, file="simulation/results/S1/S1_p005a4.rds")


phi <- 0.5
alpha <- 0

source("simulation/runSimsS1_lp.R")
saveRDS(out, file="simulation/results/S1/S1_p05a0.rds")

phi <- 0.5
alpha <- 2

source("simulation/runSimsS1_lp.R")
saveRDS(out, file="simulation/results/S1/S1_p05a2.rds")

phi <- 0.5
alpha <- 4

source("simulation/runSimsS1_lp.R")
saveRDS(out, file="simulation/results/S1/S1_p05a4.rds")

phi <- 1
alpha <- 0

source("simulation/runSimsS1_lp.R")
saveRDS(out, file="simulation/results/S1/S1_p1a0.rds")

phi <- 1
alpha <- 2

source("simulation/runSimsS1_lp.R")
saveRDS(out, file="simulation/results/S1/S1_p1a2.rds")

phi <- 1
alpha <- 4

source("simulation/runSimsS1_lp.R")
saveRDS(out, file="simulation/results/S1/S1_p1a4.rds")