library(rstan)
library(dplyr)
library(stringr)
library(xtable)


files <- list.files("case_study/results", pattern = "\\.RData$", full.names = TRUE)


extract_model_summary <- function(file_path) {
  load(file_path)  
  
  obj <- rstan::summary(long_rstan_LT)
  sm <- obj$summary
  
  # Longitudinal parameters
  beta1   <- sm[2, c(1, 4, 8)]
  Phi     <- sm[4, c(1, 4, 8)]
  sigmau  <- sm[5, c(1, 4, 8)]
  sigmav  <- sm[8, c(1, 4, 8)]
  
  # Survival parameter
  Alpha   <- sm[3, c(1, 4, 8)]
  
  # Extract subscale name from filename
  shortname <- str_replace(basename(file_path), ".RData", "")
  shortname <- str_replace(shortname, "_.*$", "")  
  
  # Combine results
  data.frame(
    Subscale = shortname,
    beta1_mean = beta1[1], beta1_LCL = beta1[2], beta1_UCL = beta1[3],
    Phi_mean = Phi[1], Phi_LCL = Phi[2], Phi_UCL = Phi[3],
    sigma_b0_mean = sigmau[1], sigma_b0_LCL = sigmau[2], sigma_b0_UCL = sigmau[3],
    sigma_b1_mean = sigmav[1], sigma_b1_LCL = sigmav[2], sigma_b1_UCL = sigmav[3],
    alpha_mean = Alpha[1], alpha_LCL = Alpha[2], alpha_UCL = Alpha[3]
  )
}



results_df <- do.call(rbind, lapply(files, extract_model_summary))

fmt <- function(m, l, u) sprintf("%.2f (%.2f, %.2f)", m, l, u)


tab_order <- c("PF", "RP", "BP", "GH", "VT", "SF", "RE", "MH", "SYMP", "IMP", "ACT")

# Reorder the data frame
results_df <- results_df %>%
  mutate(Subscale = factor(Subscale, levels = tab_order)) %>%
  arrange(Subscale) %>%
  mutate(Subscale = as.character(Subscale))

rownames(results_df) <- NULL

# ======== #
#  TABLE 1 #
# ======== #

latex_tab <- xtable(results_df)

