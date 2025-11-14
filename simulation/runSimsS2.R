set.seed(8)

source("simulation/JMbb_datasets.R") # simulation of the data sets function
source("simulation/Simulation_function.R")  # here the simulation functions 

#####################################################
# SET THE VALUES OF THE PARAMETERS THAT DONT CHANGE #
#####################################################

n <- 500 # number of subjects
K <- 4 # number of planned repeated measurements per subject, per outcome
t.max <- 6.5 # maximum follow-up time
cens <- 0.1 # the random censoring

# Fixed in scenario 2

betas <- c("Intercept" = 2.17 , "time" = -0.07) 
ntrial <- 8
D <- matrix(c(2.45,0,0,0.07),ncol=2)

nu <- 1.75 
gamma <- c(-1.6) 

#############
# SIMULATE  #
#############

nSims <- 100
out <- simulation.function(nSims = nSims, n = n, K = K, t.max = t.max, betas = betas,
                           phi = phi, ntrial = ntrial, alpha = alpha, nu = nu,
                           gamma = gamma, D = D, cens = cens) 

