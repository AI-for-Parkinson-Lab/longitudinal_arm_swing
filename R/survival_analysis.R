# vignette("timedep")
# install.packages("ggplot2")

library("survival")
library(jsonlite)
library(dplyr)
library(broom)
library(tidyr)
library(readr)

# Set measure to "median" or "95p", and med_state to "med" or "no_med"
measure <- "95p"
med_state <- "no_med"
var_colname <- ifelse(measure == "95p", "X95p", measure)

# Load the data
path_to_df <- "C:\\Users\\erik_\\Documents\\PhD\\data\\longitudinal_arm_swing\\prepared_data\\measures\\filtered_gait\\very_long_gait_segments"
filename_df <- sprintf("\\coxph_ctv_%s_%s.csv", measure, med_state)
full_path <- file.path(path_to_df, filename_df)
df = read.csv(full_path)

# + %s_delta_pct_2_week_lag"
cox_base_formula <- as.formula(sprintf("Surv(start, stop, event) ~ %s", var_colname))
cox_spline_formula <- as.formula(sprintf("Surv(start, stop, event) ~ pspline(%s) ", var_colname))
  
# Fit a basic time-varying CoxPH model
cox_base = coxph(cox_base_formula, data = df)

# Check linearity assumption
cox_spline <- coxph(
  cox_spline_formula,
  data = df
)

mart_res <- residuals(cox_base, type = "martingale")

plot(df[[var_colname]], mart_res,
     xlab = paste(tools::toTitleCase(measure), "range of motion"),
     ylab = "Martingale Residuals",
     main = "")

lines(lowess(df[[var_colname]], mart_res), col = "red")

# Based on the martingale residuals, we don't need to add spline complexity to the model for both median and 95p.

# Check proportional hazards assumption
schoenfeld_res <- cox.zph(cox_base)

label <- paste(tools::toTitleCase(gsub("_", " ", measure)), "range of motion")
plot(schoenfeld_res, var = 1, ylab = label, ylim = c(-0.25, 0.45))

# Based on the schoenfeld residuals, we should add time-interactions effects
# Add time-interaction effects: 95p (abs: log-linear, delta:  none), median (abs: none, delta: None)
tmax_global <- max(df$stop)

# Scaling
# scaled_values <- scale(df[[var_colname]])
# df[[var_colname]] <- as.numeric(scaled_values)

# + %s_delta_pct_2_week_lag
formula_time <- if (measure == '95p') {
  as.formula(sprintf(
    "Surv(start, stop, event) ~ %s + tt(%s)",
    var_colname, var_colname
  ))
} else {
  as.formula(sprintf(
    "Surv(start, stop, event) ~ %s + tt(%s)",
    var_colname, var_colname
  ))
}
cox_time <- coxph(
  formula = formula_time,
  data = df,
  tt = function(x, t, ...) {
    if (measure == "95p") {
      return(x * log(t + 1))
    } else {
      return(x * log(t + 1))
    }
  }
)

# Examine performance
summary(cox_time)
AIC(cox_time)

# Predict survival probabilities
# We cannot use the usual survfit() model because it does not work with time-varying coefficients, so
# we should calculate hazards manualy.

# Baseline cumulative hazard
base_haz <- basehaz(cox_base, centered = FALSE)

# Linear predictor function with rescaling
linear_predictor <- function(row, t, measure, tmax) {

  x <- row[[var_colname]]
  delta <- row[[paste0(var_colname, "_delta_pct_2_week_lag")]]
  
  main_coef <- coef(cox_time)[[var_colname]]
  # delta_coef <- coef(cox_time)[[paste0(var_colname, "_delta_pct_2_week_lag")]]
  tt_coef <- coef(cox_time)[[sprintf("tt(%s)", var_colname)]]
  
  time_effect <- if (measure == "95p") {
    x * log(t + 1)
  } else {
    x * log(t + 1)
  }
  
  # + delta_coef
  eta <- main_coef * x + tt_coef * time_effect
  return(eta)
}

# Step 3: Compute survival probabilities for one subject
compute_surv_prob <- function(row, base_haz_df, measure, tmax, times) {
  surv_probs <- numeric(length(times))
  
  for (i in seq_along(times)) {
    t <- times[i]
    idx <- max(which(base_haz_df$time <= t))
    if (length(idx) == 0 || is.infinite(idx)) {
      H0_t <- 0  # no baseline hazard yet
    } else {
      H0_t <- base_haz_df$hazard[idx]
    }
    eta_t <- linear_predictor(row, t, measure, tmax)
    surv_probs[i] <- exp(-H0_t * exp(eta_t))
  }
  
  return(surv_probs)
}

# Define the fixed weeks you want columns for (every 2 weeks)
fixed_weeks <- seq(2, max(df$stop), by = 2)
results <- list()

for (subject_id in unique(df$id)) {
  subject_data <- df[df$id == subject_id, ]
  if (nrow(subject_data) == 0) next
  
  latest_row <- subject_data[which.max(subject_data$stop), ]
  tmax_subj <- max(subject_data$stop) 
  
  weeks_for_subj <- fixed_weeks[fixed_weeks <= tmax_subj]
  
  surv_probs <- compute_surv_prob(latest_row, base_haz, measure, tmax_subj, weeks_for_subj)
  
  # Pad NA for weeks > subject's last time
  surv_vec <- rep(NA_real_, length(fixed_weeks))
  surv_vec[fixed_weeks <= tmax_subj] <- surv_probs
  
  results[[as.character(subject_id)]] <- surv_vec
}

results_df <- do.call(rbind, results) %>% as.data.frame()
colnames(results_df) <- paste0("week", fixed_weeks)
results_df <- cbind(id = names(results), results_df)

write.csv(results_df, file = sprintf("%s/cox_survival_predictions_tt_%s.csv", path_to_df, measure), row.names = FALSE)

# Repeat, but for the method excluding time-varying coefficients

# Predict survival per subject using their full longitudinal data
ids <- unique(df$id)

output_list <- lapply(ids, function(subj_id) {
  subject_df <- df %>% filter(id == subj_id)
  
  if (nrow(subject_df) < 1) return(NULL)  # skip empty subjects
  latest_row <- subject_df[which.max(subject_df$stop), , drop = FALSE]
  
  surv_fit <- survfit(cox_base, newdata = latest_row)
  
  # Interpolate survival at fixed time points
  surv_interp <- approx(surv_fit$time, surv_fit$surv, xout = fixed_weeks, rule = 2)$y
  
  data.frame(id = subj_id, setNames(as.list(surv_interp), paste0("week", fixed_weeks)))
})

output_df <- bind_rows(output_list)

output_file <- file.path(path_to_df, sprintf("cox_survival_predictions_base_%s.csv", measure))
write_csv(output_df, output_file)


