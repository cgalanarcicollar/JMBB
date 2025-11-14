library(xtable)

#################
#      S1       #
#################

S1_path <- "simulation/results/S1/"

S1_rds_files <- list.files(S1_path, pattern = "\\.rds$", full.names = TRUE)

# --- CP function
proportion <- function(b, se, true, level = 0.95, df = Inf) {
  qtile <- level + (1 - level)/2
  lower <- b - qt(qtile, df = df) * se
  upper <- b + qt(qtile, df = df) * se
  cp <- mean(true >= lower & true <= upper)
  return(cp)
}

# --- Main processing function ---
process_file <- function(file_path) {
  out <- readRDS(file_path)
  alpha <- out$parameters[[1]]
  
  # ---- JMBB  ----
  Alpha_JM <- out$result_JM$surv$Alpha
  est_JM <- mean(Alpha_JM[, 1], na.rm = TRUE)
  esd_JM <- sqrt(var(Alpha_JM[, 1], na.rm = TRUE))
  asd_JM <- mean(Alpha_JM[, 3], na.rm = TRUE)
  bias_JM <- if (alpha == 0) (est_JM - alpha) else (est_JM - alpha) / alpha
  CP_JM <- mean(Alpha_JM[, 4] < alpha & alpha < Alpha_JM[, 8], na.rm = TRUE) * 100
  
  # ---- TSBB ----
  Alpha_TS <- out$result_TS$survival$Alpha
  est_TS <- mean(Alpha_TS[, 1], na.rm = TRUE)
  esd_TS <- sqrt(var(Alpha_TS[, 1], na.rm = TRUE))
  asd_TS <- mean(Alpha_TS[, 2], na.rm = TRUE)
  bias_TS <- if (alpha == 0) (est_TS - alpha) else (est_TS - alpha) / alpha
  CP_TS <- proportion(Alpha_TS[, 1], Alpha_TS[, 2], alpha) * 100
  
  # ---- Combine  ----
  res <- data.frame(
    File = basename(file_path),
    Model = c("JM", "TS"),
    Bias = c(bias_JM, bias_TS),
    ESD = c(esd_JM, esd_TS),
    ASD = c(asd_JM, asd_TS),
    CP = c(CP_JM, CP_TS)
  )
  
  return(res)
}


S1_results <- do.call(rbind, lapply(S1_rds_files, process_file))

# ======== #
#  TABLE 4 #
# ======== #

S1_latex_table <- xtable(S1_results, 
                      digits = c(0, 0, 0, 3, 3, 3, 1),
                      caption = "Summary of Bias, ESD, ASD, and CP for all simulations (S1)",
                      label = "tab:summary_S1")


S1_latex_table

#################
#      S2       #
#################

S2_path <- "simulation/results/S2/"

S2_rds_files <- list.files(S2_path, pattern = "\\.rds$", full.names = TRUE)

S2_results <- do.call(rbind, lapply(S2_rds_files, process_file))

# ======== #
#  TABLE 5 #
# ======== #

S2_latex_table <- xtable(S2_results, 
                      digits = c(0, 0, 0, 3, 3, 3, 1),
                      caption = "Summary of Bias, ESD, ASD, and CP for all simulations (S2)",
                      label = "tab:summary_S2")

S2_latex_table
