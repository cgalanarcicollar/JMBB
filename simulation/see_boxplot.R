library(patchwork)
library(dplyr)
library(ggplot2)

make_file_grid <- function(base_path, prefix, phi_suffixes, phi_values,
                           alpha_suffixes, alpha_values) {
  
  # Ensure lengths match
  stopifnot(length(phi_suffixes) == length(phi_values))
  stopifnot(length(alpha_suffixes) == length(alpha_values))
  
  # expand grid of suffixes
  grid <- expand.grid(phi_suf = phi_suffixes, alpha_suf = alpha_suffixes, stringsAsFactors = FALSE)
  
  # map suffix -> numeric
  phi_map <- setNames(phi_values, phi_suffixes)
  alpha_map <- setNames(alpha_values, alpha_suffixes)
  
  # build file paths and numeric values
  grid$file <- file.path(base_path, sprintf("%s_p%s%s.rds", prefix, grid$phi_suf, grid$alpha_suf))
  grid$phi_numeric <- unname(phi_map[grid$phi_suf])
  grid$alpha_numeric <- unname(alpha_map[grid$alpha_suf])
  
  # return grid with file, phi_numeric, alpha_numeric
  grid
}

# --- load results for alpha ---
make_plot_data_alpha <- function(grid_df, jm_name = "JMBB", ts_name = "TSBB") {
  rows <- list()
  for (i in seq_len(nrow(grid_df))) {
    f <- grid_df$file[i]
    if (!file.exists(f)) stop("File not found: ", f)
    
    out <- readRDS(f)
    
    # numeric values from different scenarios
    phi <- grid_df$phi_numeric[i]
    alpha <- grid_df$alpha_numeric[i]
    
    rows[[length(rows) + 1]] <- cbind(out$result_JM$surv$Alpha[,1], alpha, phi, jm_name)
    rows[[length(rows) + 1]] <- cbind(out$result_TS$surv$Alpha[,1], alpha, phi, ts_name)
  }
  
  df <- as.data.frame(do.call(rbind, rows))
  colnames(df) <- c("est", "alpha", "phi", "model")
  
  df <- df %>%
    mutate(
      est   = as.numeric(est),
      alpha = as.numeric(alpha),
      phi   = as.factor(phi),           
      model = factor(model, levels = c(jm_name, ts_name))
    )
  
  df$bias <- df$est - df$alpha
  df$bias[df$alpha != 0] <- df$bias[df$alpha != 0] / df$alpha[df$alpha != 0]
  df
}

# --- load results for beta1 ---
make_plot_data_beta1 <- function(grid_df, jm_name = "JMBB", ts_name = "TSBB") {
  rows <- list()
  for (i in seq_len(nrow(grid_df))) {
    f <- grid_df$file[i]
    if (!file.exists(f)) stop("File not found: ", f)
    
    out <- readRDS(f)
    
    # numeric values from different scenarios
    phi <- grid_df$phi_numeric[i]
    alpha <- grid_df$alpha_numeric[i]
    
    # Extract true beta1 value from out$parameters[[3]][2]
    true_beta1 <- out$parameters[[3]][2]
    
    # Extract estimated beta1 from the correct locations
    jm_beta1 <- out$result_JM$long$beta1[,1]
    ts_beta1 <- out$result_TS$longitudinal$beta1[,1]
    
    # Calculate relative bias for each estimate
    jm_bias <- (jm_beta1 - true_beta1) / true_beta1
    ts_bias <- (ts_beta1 - true_beta1) / true_beta1
    
    rows[[length(rows) + 1]] <- cbind(jm_bias, alpha, phi, jm_name)
    rows[[length(rows) + 1]] <- cbind(ts_bias, alpha, phi, ts_name)
  }
  
  df <- as.data.frame(do.call(rbind, rows))
  colnames(df) <- c("bias", "alpha", "phi", "model")
  
  df <- df %>%
    mutate(
      bias  = as.numeric(bias),
      alpha = as.numeric(alpha),
      phi   = as.factor(phi),           
      model = factor(model, levels = c(jm_name, ts_name))
    )
  
  df
}

# --- box plot for alpha ---
plot_bias_alpha <- function(df, alpha_levels, alpha_labels, show_y_label = TRUE) {
  p <- df %>%
    mutate(
      Alpha = factor(alpha, levels = alpha_levels, labels = alpha_labels)
    ) %>%
    ggplot(aes(x = phi, y = bias, fill = model)) +
    geom_boxplot() +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    scale_fill_discrete(
      breaks = c("JMBB", "TSBB"),
      type = c(JMBB = "grey74", TSBB = "white")
    ) +
    facet_grid(cols = vars(Alpha), labeller = label_parsed) +
    xlab(expression(Phi)) +
    geom_hline(yintercept = 0, col = "red")
  
  if (show_y_label) {
    p <- p + ylab(expression("Relative bias of " * alpha))
  } else {
    p <- p + ylab(NULL)
  }
  
  return(p)
}

# --- box plot for beta1 ---
plot_bias_beta1 <- function(df, alpha_levels, alpha_labels, show_y_label = TRUE) {
  p <- df %>%
    mutate(
      Alpha = factor(alpha, levels = alpha_levels, labels = alpha_labels)
    ) %>%
    ggplot(aes(x = phi, y = bias, fill = model)) +
    geom_boxplot() +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    scale_fill_discrete(
      breaks = c("JMBB", "TSBB"),
      type = c(JMBB = "grey74", TSBB = "white")
    ) +
    facet_grid(cols = vars(Alpha), labeller = label_parsed) +
    xlab(expression(Phi)) +
    geom_hline(yintercept = 0, col = "red")
  
  if (show_y_label) {
    p <- p + ylab(expression("Relative bias of " * beta[1]))
  } else {
    p <- p + ylab(NULL)
  }
  
  return(p)
}

add_title <- function(plot, title) {
  plot + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = 14))
}

########################
#         S1           #
########################

phi_suffixes_s1 <- c("005", "05", "1")          
phi_values_s1   <- c(0.05, 0.5, 1)              

alpha_suffixes_s1 <- c("a0", "a2", "a4")
alpha_values_s1   <- c(0, 2, 4)

s1_grid <- make_file_grid(
  base_path = "simulation/results/S1",
  prefix = "S1",
  phi_suffixes = phi_suffixes_s1,
  phi_values   = phi_values_s1,
  alpha_suffixes = alpha_suffixes_s1,
  alpha_values   = alpha_values_s1
)

s1_alpha_data <- make_plot_data_alpha(s1_grid)
s1_beta1_data <- make_plot_data_beta1(s1_grid)

q1_alpha <- plot_bias_alpha(
  s1_alpha_data,
  alpha_levels = c("0", "2", "4"),
  alpha_labels = c(
    expression(paste(alpha, " = 0")),
    expression(paste(alpha, " = 2")),
    expression(paste(alpha, " = 4"))
  ),
  show_y_label = TRUE
)

q1_beta1 <- plot_bias_beta1(
  s1_beta1_data,
  alpha_levels = c("0", "2", "4"),
  alpha_labels = c(
    expression(paste(alpha, " = 0")),
    expression(paste(alpha, " = 2")),
    expression(paste(alpha, " = 4"))
  ),
  show_y_label = TRUE  # Cambiado a TRUE para que muestre el label
)

########################
#         S2           #
########################

phi_suffixes_s2 <- c("005", "05", "1")       
phi_values_s2   <- c(0.05, 0.5, 1)

alpha_suffixes_s2 <- c("a_0", "a_15", "a_3")
alpha_values_s2   <- c(0, -1.5, -3)

s2_grid <- make_file_grid(
  base_path = "simulation/results/S2",
  prefix = "S2",
  phi_suffixes = phi_suffixes_s2,
  phi_values   = phi_values_s2,
  alpha_suffixes = alpha_suffixes_s2,
  alpha_values   = alpha_values_s2
)

s2_alpha_data <- make_plot_data_alpha(s2_grid)
s2_beta1_data <- make_plot_data_beta1(s2_grid)

q2_alpha <- plot_bias_alpha(
  s2_alpha_data,
  alpha_levels = c("0", "-1.5", "-3"),
  alpha_labels = c(
    expression(paste(alpha, " = 0")),
    expression(paste(alpha, " = -1.5")),
    expression(paste(alpha, " = -3"))
  ),
  show_y_label = FALSE
)

q2_beta1 <- plot_bias_beta1(
  s2_beta1_data,
  alpha_levels = c("0", "-1.5", "-3"),
  alpha_labels = c(
    expression(paste(alpha, " = 0")),
    expression(paste(alpha, " = -1.5")),
    expression(paste(alpha, " = -3"))
  ),
  show_y_label = FALSE
)

# Add titles to each plot
q1_alpha <- add_title(q1_alpha, "Scenario 1")
q2_alpha <- add_title(q2_alpha, "Scenario 2")
q1_beta1 <- add_title(q1_beta1, "Scenario 1")
q2_beta1 <- add_title(q2_beta1, "Scenario 2")

# ============== #
# FIGURE 6 and 7 #
# ============== #

# Arrange alpha plots
final_boxplot_alpha <- (q1_alpha | q2_alpha) + plot_layout(guides = "collect")
final_boxplot_alpha

# Arrange beta1 plots
final_boxplot_beta1 <- (q1_beta1 | q2_beta1) + plot_layout(guides = "collect")
final_boxplot_beta1
