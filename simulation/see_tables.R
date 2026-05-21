library(xtable)


# --- CP function
proportion <- function(b, se, true, level = 0.95, df = Inf) {
  qtile <- level + (1 - level)/2
  lower <- b - qt(qtile, df = df) * se
  upper <- b + qt(qtile, df = df) * se
  cp <- mean(true >= lower & true <= upper)
  return(cp)
}

# --- Main processing function (for Alpha and Beta_1) ---

process_file_alpha <- function(file_path) {
  out <- readRDS(file_path)
  alpha <- out$parameters[[1]] # Real association parameter
  
  # ---- JMBB Alpha ----
  Alpha_JM <- out$result_JM$surv$Alpha
  est_JM <- mean(Alpha_JM[, 1], na.rm = TRUE)
  esd_JM <- sqrt(var(Alpha_JM[, 1], na.rm = TRUE))
  asd_JM <- mean(Alpha_JM[, 3], na.rm = TRUE)
  bias_JM <- if (alpha == 0) (est_JM - alpha) else (est_JM - alpha) / alpha
  CP_JM <- mean(Alpha_JM[, 4] < alpha & alpha < Alpha_JM[, 8], na.rm = TRUE) * 100
  
  # ---- TSBB Alpha ----
  Alpha_TS <- out$result_TS$survival$Alpha
  est_TS <- mean(Alpha_TS[, 1], na.rm = TRUE)
  esd_TS <- sqrt(var(Alpha_TS[, 1], na.rm = TRUE))
  asd_TS <- mean(Alpha_TS[, 2], na.rm = TRUE)
  bias_TS <- if (alpha == 0) (est_TS - alpha) else (est_TS - alpha) / alpha
  CP_TS <- proportion(Alpha_TS[, 1], Alpha_TS[, 2], alpha) * 100
  
  # ---- Combine para Alpha ----
  res <- data.frame(
    File = basename(file_path),
    Parameter = "Alpha",
    Model = c("JM", "TS"),
    Bias = c(bias_JM, bias_TS),
    ESD = c(esd_JM, esd_TS),
    ASD = c(asd_JM, asd_TS),
    CP = c(CP_JM, CP_TS)
  )
  
  return(res)
}

process_file_beta1 <- function(file_path) {
  out <- readRDS(file_path)
  beta1 <- out$parameters[[3]][2]  # Real longitudinal slope parameter 
  
  # ---- JMBB Beta1 ----
  Beta_JM <- out$result_JM$long$beta1
  est_JM <- mean(Beta_JM[, 1], na.rm = TRUE)
  esd_JM <- sqrt(var(Beta_JM[, 1], na.rm = TRUE))
  asd_JM <- mean(Beta_JM[, 3], na.rm = TRUE)
  bias_JM <- (est_JM - beta1) / beta1
  CP_JM <- mean(Beta_JM[, 4] < beta1 & beta1 < Beta_JM[, 8], na.rm = TRUE) * 100
  
  # ---- TSBB Beta1 ----
  Beta_TS <- out$result_TS$longitudinal$beta1
  est_TS <- mean(Beta_TS[, 1], na.rm = TRUE)
  esd_TS <- sqrt(var(Beta_TS[, 1], na.rm = TRUE))
  asd_TS <- mean(Beta_TS[, 2], na.rm = TRUE)
  bias_TS <- (est_TS - beta1) / beta1
  CP_TS <- proportion(Beta_TS[, 1], Beta_TS[, 2], beta1) * 100
  
  # ---- Combine Beta1 ----
  res <- data.frame(
    File = basename(file_path),
    Parameter = "Beta1",
    Model = c("JM", "TS"),
    Bias = c(bias_JM, bias_TS),
    ESD = c(esd_JM, esd_TS),
    ASD = c(asd_JM, asd_TS),
    CP = c(CP_JM, CP_TS)
  )
  
  return(res)
}


#################
#      S1       #
#################

S1_path <- "simulation/results/S1/"

S1_rds_files <- list.files(S1_path, pattern = "\\.rds$", full.names = TRUE)

S1_results_alpha <- do.call(rbind, lapply(S1_rds_files, process_file_alpha))
S1_results_beta1 <- do.call(rbind, lapply(S1_rds_files, process_file_beta1))

S1_results <- rbind(S1_results_alpha, S1_results_beta1)



#################
#      S2       #
#################

S2_path <- "simulation/results/S2/"

S2_rds_files <- list.files(S2_path, pattern = "\\.rds$", full.names = TRUE)

S2_results_alpha <- do.call(rbind, lapply(S2_rds_files, process_file_alpha))
S2_results_beta1 <- do.call(rbind, lapply(S2_rds_files, process_file_beta1))

S2_results <- rbind(S2_results_alpha, S2_results_beta1)

# ======== #
#  TABLES  #
# ======== #


# Table 4 --- Beta1 Results (S1) 
S1_beta1_table <- xtable(S1_results[S1_results$Parameter == "Beta1",], 
                         digits = c(0, 0, 0, 0, 2, 3, 3, 1),  
                         caption = "Summary for Beta1 parameter (S1)",
                         label = "tab:summary_S1_Beta1")
print(S1_beta1_table, include.rownames = FALSE)


# Table 5 ---  Beta1 Results (S2)
S2_beta1_table <- xtable(S2_results[S2_results$Parameter == "Beta1",], 
                         digits = c(0, 0, 0, 0, 2, 3, 3, 1),  
                         caption = "Summary for Beta1 parameter (S2)",
                         label = "tab:summary_S2_Beta1")
print(S2_beta1_table, include.rownames = FALSE)


# Table 6 ---  Alpha Results (S1)
S1_alpha_table <- xtable(S1_results[S1_results$Parameter == "Alpha",], 
                         digits = c(0, 0, 0, 0, 2, 3, 3, 1),  
                         caption = "Summary for Alpha parameter (S1)",
                         label = "tab:summary_S1_alpha")

print(S1_alpha_table, include.rownames = FALSE)

# Table 7 ---  Alpha Results (S2)
S2_alpha_table <- xtable(S2_results[S2_results$Parameter == "Alpha",], 
                         digits = c(0, 0, 0, 0, 2, 3, 3, 1),  
                         caption = "Summary for Alpha parameter (S2)",
                         label = "tab:summary_S2_alpha")

print(S2_alpha_table, include.rownames = FALSE)


