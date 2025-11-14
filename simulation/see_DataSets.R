library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(xtable)



# =================================================================
# Function to get baseline measurement for first simulated data set
# =================================================================

get_baseline_data <- function(file_path) {
  result_list <- readRDS(file_path)
  data_long <- result_list$Data[[1]]$longitudinal
  
  if ("meas" %in% names(data_long)) {
    baseline_data <- data_long[data_long$meas == "t0", ]
  } else {
    baseline_data <- data_long[!duplicated(data_long$id), ]
  }
  baseline_data
}

phi_vals <- c("005", "05", "1")



# S1
baseline_S1 <- lapply(phi_vals, function(phi) {
  get_baseline_data(paste0("simulation/results/S1/S1_p", phi, "a0.rds"))
})
names(baseline_S1) <- phi_vals

max_y_S1 <- max(sapply(baseline_S1, function(df) {
  ggplot_build(
    ggplot(df, aes(x = y)) + geom_histogram(binwidth = 1)
  )$data[[1]]$count
}))

plots_S1 <- lapply(phi_vals, function(phi) {
  ggplot(baseline_S1[[phi]], aes(x = y)) +
    geom_histogram(binwidth = 1, fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
    scale_x_continuous(breaks = seq(0, 24, by = 1)) +
    ylim(0, max_y_S1) +
    theme_bw() +
    labs(
      x = "Score", y = "frequency") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
})

# S2
baseline_S2 <- lapply(phi_vals, function(phi) {
  get_baseline_data(paste0("simulation/results/S2/S2_p", phi, "a_0.rds"))
})
names(baseline_S2) <- phi_vals

max_y_S2 <- max(sapply(baseline_S2, function(df) {
  ggplot_build(
    ggplot(df, aes(x = y)) + geom_histogram(binwidth = 1)
  )$data[[1]]$count
}))

plots_S2 <- lapply(phi_vals, function(phi) {
  ggplot(baseline_S2[[phi]], aes(x = y)) +
    geom_histogram(binwidth = 1, fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
    scale_x_continuous(breaks = seq(0, 8, by = 1)) +
    ylim(0, max_y_S2) +
    theme_bw() +
    labs(
      x = "Score", y = "frequency") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
})

# ======== #
# FIGURE 2 #
# ======== #

plots_S1[[1]] + plots_S1[[2]] + plots_S1[[3]]


# ======== #
# FIGURE 3 #
# ======== #

plots_S2[[1]] + plots_S2[[2]] + plots_S2[[3]]


# =========================================================
# Function to compute frequency / measurement statistics
# =========================================================
result_freq <- function(file_path) {
  result_list <- readRDS(file_path)
  
  n_iter <- length(result_list$Data)
  events <- numeric(n_iter)
  record_meas <- numeric(n_iter)
  freq_summary <- matrix(0, ncol = 4, nrow = n_iter)
  
  for (i in seq_len(n_iter)) {
    data_long <- result_list$Data[[i]]$longitudinal
    data_surv <- result_list$Data[[i]]$survival
    
    freq_i <- tabulate(data_long$id)
    tab_counts <- table(factor(freq_i, levels = 1:4))
    
    freq_summary[i, ] <- as.numeric(tab_counts)
    record_meas[i] <- nrow(data_long)
    events[i] <- sum(data_surv$event)
  }
  
  meas_table <- data.frame(
    n_meas = 1:4,
    freq = colMeans(freq_summary)
  )
  
  bar_plot <- ggplot(meas_table, aes(x = n_meas, y = freq)) +
    geom_bar(stat = 'identity', fill = "#3399FF", color = "#e9ecef", alpha = 0.9) +
    scale_x_continuous(breaks = 1:4) +
    ylim(0, 420) +
    theme_bw() +
    labs(x = "measurement", y = "frequency") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
  
  list(
    file = basename(file_path),
    plot = bar_plot,
    mean_events = mean(events),
    mean_record_meas = mean(record_meas),
    freq1 = mean(freq_summary[, 1])
  )
}

# ======================================== #
# Define all combinations for subscenarios #
# ======================================== #


# Scenario 1: alpha = 0, 2, 4  -> files like "S1_p005a0.rds"
alpha_vals_S1 <- c("0", "2", "4")
files_S1 <- expand.grid(phi = phi_vals, alpha = alpha_vals_S1) %>%
  mutate(path = paste0("simulation/results/S1/S1_p", phi, "a", alpha, ".rds"))

# Scenario 2: alpha = 0, 1.5, 3  -> files like "S2_p005a_0.rds"
alpha_vals_S2 <- c("0", "15", "3")
files_S2 <- expand.grid(phi = phi_vals, alpha = alpha_vals_S2) %>%
  mutate(path = paste0("simulation/results/S2/S2_p", phi, "a_", alpha, ".rds"))

# Combine 
all_files <- bind_rows(
  files_S1 %>% mutate(scenario = "S1"),
  files_S2 %>% mutate(scenario = "S2")
)

# load results

all_results <- lapply(all_files$path, result_freq)
names(all_results) <- paste0(all_files$scenario, "_p", all_files$phi, "_a", all_files$alpha)


# Events and measurements summary

summary_table <- data.frame(
  Scenario = rep(c("Scenario 1", "Scenario 2"), each = 3),
  Alpha = rep(c("weak", "moderate", "strong"), times = 2),
  Events = NA,
  Meas = NA,
  Freq1 = NA
)

# S1
for (i in seq_along(alpha_vals_S1)) {

  subset_names <- paste0("S1_p", phi_vals, "_a", alpha_vals_S1[i])
  subset_results <- all_results[subset_names]
  
  # compute mean across phi
  mean_events <- mean(sapply(subset_results, `[[`, "mean_events"))
  mean_record <- mean(sapply(subset_results, `[[`, "mean_record_meas"))
  mean_freq1  <- mean(sapply(subset_results, `[[`, "freq1"))
  
  summary_table[i, c("Events", "Meas", "Freq1")] <- c(mean_events, mean_record, mean_freq1)
}

# S2 
for (i in seq_along(alpha_vals_S2)) {
  subset_names <- paste0("S2_p", phi_vals, "_a", alpha_vals_S2[i])
  subset_results <- all_results[subset_names]
  
  # compute mean across phi
  mean_events <- mean(sapply(subset_results, `[[`, "mean_events"))
  mean_record <- mean(sapply(subset_results, `[[`, "mean_record_meas"))
  mean_freq1  <- mean(sapply(subset_results, `[[`, "freq1"))
  
  summary_table[i + 3, c("Events", "Meas", "Freq1")] <- c(mean_events, mean_record, mean_freq1)
}

# ======== #
#  TABLE 3 #
# ======== #

print(
  xtable(
    summary_table,
    include.rownames = FALSE,
    comment = FALSE,
    sanitize.text.function = identity
  ),
  include.rownames = FALSE
)



# Barplots for same phi and variations in association parameter

# ======== #
# FIGURE 4 #
# ======== #

final_plot_S1 <- all_results$S1_p005_a0$plot +
  all_results$S1_p005_a2$plot +
  all_results$S1_p005_a4$plot +
  plot_layout(guides = "collect")

final_plot_S1

# ======== #
# FIGURE 5 #
# ======== #

final_plot_S2 <- all_results$S2_p005_a0$plot +
  all_results$S2_p005_a15$plot +
  all_results$S2_p005_a3$plot +
  plot_layout(guides = "collect")

final_plot_S2
