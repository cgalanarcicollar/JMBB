library("MASS")
library("statmod")
library("rstan")
library("adaptMCMC")
library("emdbook")
library("pracma")

source("simulation/JMbb_datasets.R")
source("dyn_pred/dynpred_funcL.R")


set.seed(902)

n <- 500 # number of subjects
t.max <- 6.5 # maximum follow-up time
nu <- 1.2 # shape for the Weibull baseline hazard
gamma <- c(-3.5) 

cens <- 0.1 # the random censoring

betas <- c("Intercept" = -0.19 , "time" = 0.03) 
ntrial <- 24
D <- matrix(c(1.6,0,0,0.05),ncol=2)

phi <- 0.05
alpha <- 2


# Simulate jm with PRO longitudinal response 

jm_data <- data_simulation(n = n, K = 4, t.max = t.max, betas = betas, 
                           phi = phi, ntrial = ntrial, alpha = alpha,
                           nu = nu, gamma = gamma, D = D, cens = cens)

jm_data$longitudinal <- cbind(jm_data$longitudinal, m = ntrial)

long_data <- jm_data$longitudinal
surv_data <- jm_data$survival

# Perform dynamic predictions : 

# First divide training and test samples 

test_id <- sample(1:n, 2) #two test id 

test <- long_data[(long_data$id %in% test_id),]
test_1 <- surv_data[(surv_data$id %in% test_id),]

train <- long_data[!(long_data$id %in% test_id),]
train$idx<-match(train$id, unique(train$id))
train_1 <-  surv_data[!(surv_data$id %in% test_id),]
train_1$idx<-match(train_1$id, unique(train_1$id))


# Fit the model to train data set 

# GAUS LEGENDRE QUADRATRUE  
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes   # nodes
wk <- glq$weights # weights
k <- length(xk)   # K-points

dataset_stan <- list(N= nrow(train) ,
                     n=nrow(train_1), 
                     nTrial = rep(ntrial, nrow(train)),
                     m = ntrial,
                     y=train$y,
                     times=train$time, 
                     ID=as.numeric(train$idx),
                     Time=train_1$Time,
                     status=train_1$event,
                     L=train_1$L,
                     K= k, xk= xk, wk= wk)


JMBB_train  <- stan(file = "simulation/JML_lp.stan",
                  data = dataset_stan,pars=c("betas","Alpha","Sigma","phi","gamma","nu"),
                  init = 0, 
                  chains = 3,
                  iter = 2000,
                  cores = getOption("mc.cores",3))

saveRDS(JMBB_train, file="dyn_pred/JMBB_train.rds")

# The id data to which we aim to perform the dynamic prediction

Long <- test[which(test$id==test_id[2]),]
surv <- test_1[which(test_1$id==test_id[2]),]

data_id <- list(Long=Long,surv=surv)

update_JM <- list()

for(i in 1:nrow(data_id$Long)){
  dd <- data_id
  dd$Long <- dd$Long[1:i,]
  update_JM[[i]] <- update_bi(data=dd, fit=JMBB_train, k=15, iter=1000)
}

saveRDS(update_JM, file="dyn_pred/update_JM.rds")


# Update longitudinal
n.meas <- nrow(data_id$Long)

y_pred <- list()
y_mean <- y_upper <- y_lower <- matrix(NA,n.meas,n.meas)

y_pred_gen <- list()
y_mean_gen <- y_upper_gen <- y_lower_gen <- matrix(NA,n.meas,n.meas)


for (tt in 1:n.meas) {
  beta0 <-  update_JM[[tt]]$betas.1
  beta1 <-  update_JM[[tt]]$betas.2
  bi0 <- update_JM[[tt]]$bi.1
  bi1 <-  update_JM[[tt]]$bi.2
  y_pred[[tt]] <- matrix(NA, ncol = length(1:tt), nrow=1000)
  y_pred_gen[[tt]] <- matrix(NA, ncol = length(1:tt), nrow=1000)
  
  
  for (j in 1:tt) {
    lp <- (beta0+bi0)+(beta1+bi1)*data_id$Long$time[j]
    p <- exp(lp)/(1+exp(lp))
    y_pred[[tt]][,j] <- data$Long$m[1]*p
    
    lp_gen <- (beta0)+(beta1)*data_id$Long$time[j]
    p_gen <- exp(lp_gen)/(1+exp(lp_gen))
    y_pred_gen[[tt]][,j] <- data_id$Long$m[1]*p_gen
  }
  
  y_mean[tt,1:tt] <- apply(y_pred[[tt]],2,mean)
  y_upper[tt,1:tt] <- apply(y_pred[[tt]],2,quantile,probs=0.975)
  y_lower[tt,1:tt] <- apply(y_pred[[tt]],2,quantile,probs=0.025)
  
  y_mean_gen[tt,1:tt] <- apply(y_pred_gen[[tt]],2,mean)
  y_upper_gen[tt,1:tt] <- apply(y_pred_gen[[tt]],2,quantile,probs=0.975)
  y_lower_gen[tt,1:tt] <- apply(y_pred_gen[[tt]],2,quantile,probs=0.025)
}



# Update survival
m <- data_id$Long$m[1]
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes   # nodes
wk <- glq$weights # weights
K <- length(xk)   # K-points


Surv_cond_mean <- Surv_cond_upper <- Surv_cond_lower <- matrix(NA, ncol = n.meas, nrow = 10)

for (tt in 1:n.meas) {
  beta0 <-  update_JM[[tt]]$betas.1
  beta1 <-  update_JM[[tt]]$betas.2
  bi0 <- update_JM[[tt]]$bi.1
  bi1 <-  update_JM[[tt]]$bi.2
  Alpha <- update_JM[[tt]]$Alpha
  nu <- update_JM[[tt]]$nu
  gamma <- update_JM[[tt]]$gamma
  
  t.grid <- seq(data_id$Long$time[tt],data_id$Long$time[tt]+3,len = 10)
  Surv <- matrix(NA,ncol = 1000, nrow = length(t.grid))
  for (i in 1:length(t.grid)) {
    
    
    cumHazK <- matrix(NA, ncol = 1000, nrow = K)
    for(k in 1:K){
      cumHazK[k,] =  nu * (t.grid[i]/2*(xk[k]+1))^(nu-1) *
        exp( gamma + Alpha  * exp((beta0 + bi0) + (beta1 + bi1)*t.grid[i]/2*(xk[k]+1))/(1+exp((beta0 + bi0) + (beta1 + bi1)*t.grid[i]/2*(xk[k]+1))))
    }
    
    Surv[i,] = exp(-t.grid[i]/ 2 * dot(wk, cumHazK))
    
  }
  Surv_cond <-t( t(Surv)/Surv[1,])
  
  Surv_cond_mean[,tt] <- apply(Surv_cond,1,mean)
  Surv_cond_upper[,tt] <- apply(Surv_cond,1,quantile, probs = 0.975)
  Surv_cond_lower[,tt] <- apply(Surv_cond,1,quantile, probs = 0.025)
}



#################################
#   EXAMPLE PLOTS OF FIGURE 1   #
#################################

#1st meas 
plot.new()
par(mar = c(5, 4, 4, 4) + 0.1)


plot.window(xlim=c(0,data_id$Long$time[n.meas]+3), ylim=c(0,ntrial+1))
points(data_id$Long$time[1:1],data_id$Long$y[1:1],pch=8)
axis(1)
axis(2, col.axis="black")
box()

plot.window(xlim=c(0,data_id$Long$time[n.meas]+3), ylim=c(0,1))
lines(seq(data_id$Long$time[1],data_id$Long$time[1]+3,len = 10), Surv_cond_mean[,2] )
lines(seq(data_id$Long$time[1],data_id$Long$time[1]+3,len = 10), Surv_cond_upper[,2], col="red", lty = 2)
lines(seq(data_id$Long$time[1],data_id$Long$time[1]+3,len = 10), Surv_cond_lower[,2], col="red", lty = 2)
abline(v=data_id$Long$time[1],lty = 2)
axis(4)

title("Patient X")
mtext("Survival probability", side = 4, las=3, line=2, col="black")
mtext("PRO", side = 2, las=3, line=2.5, col="black")
mtext("years", side = 1, line=2, col="black")



#2nd meas
plot.new()
par(mar = c(5, 4, 4, 4) + 0.1)


plot.window(xlim=c(0,data_id$Long$time[n.meas]+3), ylim=c(0,ntrial+1))
points(data_id$Long$time[1:2],data_id$Long$y[1:2],pch=8)
points(data_id$Long$time[1:2],y_mean[2,1:2],type = "l",col="black")
lines(data_id$Long$time[1:2],y_mean_gen[2,1:2],lty = 2,col="black")
axis(1)
axis(2, col.axis="black")
box()

plot.window(xlim=c(0,data_id$Long$time[n.meas]+3), ylim=c(0,1))
lines(seq(data_id$Long$time[2],data_id$Long$time[2]+3,len = 10), Surv_cond_mean[,2] )
lines(seq(data_id$Long$time[2],data_id$Long$time[2]+3,len = 10), Surv_cond_upper[,2], col="red", lty = 2)
lines(seq(data_id$Long$time[2],data_id$Long$time[2]+3,len = 10), Surv_cond_lower[,2], col="red", lty = 2)
abline(v=data_id$Long$time[2],lty = 2)
axis(4)

title("Patient X")
mtext("Survival probability", side = 4, las=3, line=2, col="black")
mtext("PRO", side = 2, las=3, line=2.5, col="black")
mtext("years", side = 1, line=2, col="black")


# 3rd meas

plot.new()
par(mar = c(5, 4, 4, 4) + 0.1)

plot.window(xlim=c(0,data_id$Long$time[n.meas]+3), ylim=c(0,ntrial+1))
points(data_id$Long$time[1:3],data_id$Long$y[1:3],pch=8)
points(data_id$Long$time[1:3],y_mean[3,1:3],type = "l",col="black")
lines(data_id$Long$time[1:3],y_mean_gen[3,1:3],lty = 2,col="black")
axis(1)
axis(2, col.axis="black")
box()

plot.window(xlim=c(0,data_id$Long$time[n.meas]+3), ylim=c(0,1))
lines(seq(data_id$Long$time[3],data_id$Long$time[3]+3,len = 10), Surv_cond_mean[,3] )
lines(seq(data_id$Long$time[3],data_id$Long$time[3]+3,len = 10), Surv_cond_upper[,3], col="red", lty = 2)
lines(seq(data_id$Long$time[3],data_id$Long$time[3]+3,len = 10), Surv_cond_lower[,3], col="red", lty = 2)
abline(v=data_id$Long$time[3],lty = 2)
axis(4)

title("Patient X")
mtext("Survival probability", side = 4, las=3, line=2, col="black")
mtext("PRO", side = 2, las=3, line=2.5, col="black")
mtext("years", side = 1, line=2, col="black")
