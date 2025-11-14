rm(list=ls())

library("MASS")
library("PROreg")
library("splines")
library("statmod")
library("rstan")

phi <- 0.05
alpha <- 0

source("simulation/runSimsS2.R")
saveRDS(out, file="simulation/results/S2/S2_p005a_0.rds")

phi <- 0.05
alpha <- -1.5

source("simulation/runSimsS2.R")
saveRDS(out, file="simulation/results/S2/S2_p005a_15.rds")


phi <- 0.05
alpha <- -3

source("simulation/runSimsS2.R")
saveRDS(out, file="simulation/results/S2/S2_p005a_3.rds")


phi <- 0.5
alpha <- 0

source("simulation/runSimsS2.R")
saveRDS(out, file="simulation/results/S2/S2_p05a_0.rds")

phi <- 0.5
alpha <- -1.5

source("simulation/runSimsS2.R")
saveRDS(out, file="simulation/results/S2/S2_p05a_15.rds")

phi <- 0.5
alpha <- -3

source("simulation/runSimsS2.R")
saveRDS(out, file="simulation/results/S2/S2_p05a_3.rds")

phi <- 1
alpha <- 0

source("simulation/runSimsS2.R")
saveRDS(out, file="simulation/results/S2/S2_p1a_0.rds")

phi <- 1
alpha <- -1.5

source("simulation/runSimsS2.R")
saveRDS(out, file="simulation/results/S2/S2_p1a_15.rds")

phi <- 1
alpha <- -3

source("simulation/runSimsS2.R")
saveRDS(out, file="simulation/results/S2/S2_p1a_3.rds")