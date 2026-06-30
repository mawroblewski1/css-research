# Statistical Analysis of Reproducibility Predictors
# Data: reprstitcherOutputFeb8.csv (output of reprstitcher.py)
#
# To run the full script from the RStudio console:
#   source("reprstitcherAnalysis_clean.R", echo = TRUE)
#
# Required packages (install once from the R console if needed):
#   install.packages("ggplot2")
#   install.packages("rgl")

library(ggplot2)  # 2D plotting
library(rgl)      # 3D plotting (used in exploratory section at bottom)

# =============================================================================
# SECTION 1: DATA LOADING AND FILTERING
# =============================================================================

# Load the stitched output (predictor metrics + ground truth reproducibility)
all_data <- read.csv("dataFiles/reprstitcherOutputFeb8.csv")
cat("Total rows (all build pairs):", nrow(all_data), "\n")

# Keep only the one sampled build pair per repository
# (reprpredictor.py selects one commit pair per repo; others are retained in
# the CSV for reference but excluded from analysis to avoid pseudoreplication)
all_data <- subset(all_data, sampled == "True")
cat("Rows after keeping one sample per repo:", nrow(all_data), "\n")

# Keep only Maven projects
# (availability >= 99 indicates an unsupported build system)
all_data <- subset(all_data, availability < 99)
cat("Rows after keeping Maven projects only:", nrow(all_data), "\n")

# Keep only projects where the repo was successfully downloaded and parsed
# (availability >= 2 indicates a download or parsing error)
maven_data <- subset(all_data, availability < 2)
cat("Rows after keeping successfully parsed projects:", nrow(maven_data), "\n")

# Subset: projects where the parent POM file was found (availability == 0)
# and further split by whether the parent POM had any listed dependencies
parent_pom_data       <- subset(maven_data, availability < 1)
parent_pom_with_deps  <- subset(maven_data, deps_parent_total > 0)
parent_pom_zero_deps  <- subset(maven_data, deps_parent_total == 0)

# Note: all projects in maven_data had at least one dependency project-wide
# (no projects with deps_proj_total == 0 were found in this dataset)

# =============================================================================
# SECTION 2: DESCRIPTIVE STATISTICS
# =============================================================================

cat("\n--- Reproducibility by dependency presence ---\n")

# Overall mean reproducibility (ground truth)
cat("Overall mean reproducibility:", mean(maven_data$gt_repr), "\n")

# Projects with zero vs. nonzero project-wide dependencies
# Note: no projects had zero project-wide dependencies in this dataset
cat("Mean reproducibility, zero project-wide deps:",
    mean(maven_data$gt_repr[maven_data$deps_proj_total == 0]), "\n")
cat("Mean reproducibility, nonzero project-wide deps:",
    mean(maven_data$gt_repr[maven_data$deps_proj_total > 0]), "\n")

# Projects with zero vs. nonzero parent POM dependencies
# (distinct from whether a parent POM file was found at all)
cat("Mean reproducibility, zero parent POM deps:",
    mean(maven_data$gt_repr[maven_data$deps_parent_total == 0]), "\n")
cat("Mean reproducibility, nonzero parent POM deps:",
    mean(maven_data$gt_repr[maven_data$deps_parent_total > 0]), "\n")

cat("\n--- Reproducibility by project-wide versioning rate ---\n")

# Compute the project-wide dependency versioning rate:
# proportion of all dependencies (across all POM files in the project)
# that have an explicit version attribute
versioning_rate <- maven_data$proj_vers / maven_data$deps_proj_total

# Reproducibility at versioning rate extremes and in between
cat("Mean reproducibility, no versioning (rate = 0):",
    mean(maven_data$gt_repr[maven_data$proj_vers == 0]), "\n")
cat("Mean reproducibility, full versioning (rate = 1):",
    mean(maven_data$gt_repr[maven_data$proj_vers == maven_data$deps_proj_total]), "\n")
cat("Mean reproducibility, partial versioning (0 < rate < 1):",
    mean(maven_data$gt_repr[
      (maven_data$proj_vers != maven_data$deps_proj_total) &
      (maven_data$proj_vers != 0)
    ]), "\n")
cat("Mean reproducibility, low versioning (rate < 25%):",
    mean(maven_data$gt_repr[maven_data$proj_vers < 0.25 * maven_data$deps_proj_total]), "\n")

# How many projects have a very low versioning rate (< 2%)?
cat("Number of projects with versioning rate < 2%:",
    sum(versioning_rate < 0.02), "\n")

cat("\n--- Parent POM subset reproducibility ---\n")
cat("Mean reproducibility, parent POM with deps:", mean(parent_pom_with_deps$gt_repr), "\n")
cat("Mean reproducibility, parent POM zero deps:", mean(parent_pom_zero_deps$gt_repr), "\n")

# =============================================================================
# SECTION 3: MODEL FITTING — CANDIDATE PREDICTORS
# =============================================================================
# For each candidate predictor:
#   1. Apply a transformation where appropriate (log, sqrt, or ratio)
#   2. Fit a quasibinomial GLM with polynomial basis expansion
#      - Quasibinomial chosen over binomial to account for overdispersion
#        in the outcome (gt_repr is a ratio of matches to attempts,
#        not a strict 0/1 outcome)
#      - Polynomial degree chosen by visual inspection of the fitted curve
#   3. Plot the fitted curve against the raw data
#   4. Print coefficient table with p-values

# Helper function: fit, plot, and summarize a quasibinomial polynomial model
fit_and_plot <- function(x_vals, y_vals, degree, x_label, title_label) {
  df <- data.frame(x = x_vals, y = y_vals)
  model <- glm(y ~ poly(x, degree), data = df, family = quasibinomial)
  df$predicted <- predict(model, type = "response")

  p <- ggplot(df, aes(x, y)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = predicted), color = "blue", size = 1) +
    labs(title = title_label,
         x = x_label,
         y = "Predicted probability of reproducibility") +
    theme_minimal()
  print(p)

  cat("\n--- Model summary:", title_label, "---\n")
  print(summary(model))

  return(model)
}

# Companion helper: same as fit_and_plot, but excludes observations at the
# extremes (x == 0 or x == 1) before fitting. Intended for rate/proportion
# variables (e.g. versioning rate), where many projects pile up at the
# extremes (0% or 100% versioned) and may be disproportionately influencing
# the fitted curve. Labeled clearly in the plot title and a console message.
fit_and_plot_trimmed <- function(x_vals, y_vals, degree, x_label, title_label) {
  n_total <- length(x_vals)
  keep <- (x_vals != 0) & (x_vals != 1)
  n_excluded <- n_total - sum(keep)

  cat("\n[Trimmed variant of:", title_label, "]\n")
  cat("Excluding", n_excluded, "of", n_total,
      "observations at the extremes (x = 0 or x = 1)\n")

  trimmed_title <- paste0(title_label, " (extremes excluded, n=", sum(keep), ")")

  fit_and_plot(
    x_vals  = x_vals[keep],
    y_vals  = y_vals[keep],
    degree  = degree,
    x_label = x_label,
    title_label = trimmed_title
  )
}

# --- 3.1 Project-wide dependency versioning rate ---
# Predictor: proportion of all dependencies across all POM files
# in the project that have an explicit version attribute
# Transformation: none (already a proportion in [0, 1])
# Polynomial degree: 4 (quartic), captures the hump-shaped pattern

model_versioning_rate <- fit_and_plot(
  x_vals  = maven_data$proj_vers / maven_data$deps_proj_total,
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Project-wide dependency versioning rate",
  title_label = "Quartic quasibinomial GLM: Project-wide versioning rate"
)

# Companion plot: same model, excluding projects at the extremes (0% or
# 100% versioned). Many projects cluster at these extremes, and this
# trimmed view checks whether the hump-shaped relationship still holds
# among projects with genuinely intermediate versioning practices.
model_versioning_rate_trimmed <- fit_and_plot_trimmed(
  x_vals  = maven_data$proj_vers / maven_data$deps_proj_total,
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Project-wide dependency versioning rate",
  title_label = "Quartic quasibinomial GLM: Project-wide versioning rate"
)

# Alternative fit with degree 6 — for comparison
# (Higher degree captures more curvature but risks overfitting at n=50)
model_versioning_rate_deg6 <- fit_and_plot(
  x_vals  = maven_data$proj_vers / maven_data$deps_proj_total,
  y_vals  = maven_data$gt_repr,
  degree  = 6,
  x_label = "Project-wide dependency versioning rate",
  title_label = "Degree-6 quasibinomial GLM: Project-wide versioning rate (alternative)"
)

# Companion plot: degree-6 model, extremes excluded
model_versioning_rate_deg6_trimmed <- fit_and_plot_trimmed(
  x_vals  = maven_data$proj_vers / maven_data$deps_proj_total,
  y_vals  = maven_data$gt_repr,
  degree  = 6,
  x_label = "Project-wide dependency versioning rate",
  title_label = "Degree-6 quasibinomial GLM: Project-wide versioning rate (alternative)"
)

# --- 3.2 Log of raw count of versioned dependencies (project-wide) ---
# Predictor: total number of dependencies with explicit versions across all POM files
# Transformation: log (to reduce right skew)
# Polynomial degree: 4

model_log_versioned_count <- fit_and_plot(
  x_vals  = log(maven_data$proj_vers),
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Log of versioned dependency count (project-wide)",
  title_label = "Quartic quasibinomial GLM: Log versioned dependency count"
)

# --- 3.3 Commit age ---
# Predictor: time elapsed since the failing commit was created
# Note: 'age' in the data is stored as a negative number of days
# (negative of days elapsed), so we negate it for interpretability.
# The dataset is restricted to a ~40-day window from mining.
# Transformation: none (linear scale, small range)
# Polynomial degree: 2 (quadratic)
# NOTE: In the original script this model was accidentally fit on the
# log versioned count dataframe instead of the age dataframe. Fixed here.

model_age <- fit_and_plot(
  x_vals  = -maven_data$age,  # convert to positive days elapsed
  y_vals  = maven_data$gt_repr,
  degree  = 2,
  x_label = "Days since failing commit was created",
  title_label = "Quadratic quasibinomial GLM: Commit age"
)

# --- 3.4 Log of number of POM files in the project ---
# Predictor: total number of Maven POM files found across the repository
# (a proxy for project complexity / multi-module structure)
# Transformation: log
# Polynomial degree: 2
# Note: standard binomial used here (not quasibinomial) — revisit if
# overdispersion is apparent from the summary

model_log_pom_count <- fit_and_plot(
  x_vals  = log(maven_data$num_poms),
  y_vals  = maven_data$gt_repr,
  degree  = 2,
  x_label = "Log of number of POM files in repository",
  title_label = "Quadratic binomial GLM: Log POM file count"
)
# Override family to quasibinomial for consistency with other models if needed:
# model_log_pom_count <- glm(gt_repr ~ poly(log(num_poms), 2),
#                            data = maven_data, family = quasibinomial)

# --- 3.5 Parent POM dependency versioning rate ---
# Predictor: proportion of dependencies in the parent POM file specifically
# that have an explicit version attribute
# Subset: only projects where the parent POM had at least one dependency
# Transformation: none (already a proportion in [0, 1])
# Polynomial degree: 4

model_parent_versioning_rate <- fit_and_plot(
  x_vals  = parent_pom_with_deps$parent_vers / parent_pom_with_deps$deps_parent_total,
  y_vals  = parent_pom_with_deps$gt_repr,
  degree  = 4,
  x_label = "Parent POM dependency versioning rate",
  title_label = "Quartic quasibinomial GLM: Parent POM versioning rate"
)

# Companion plot: same model, excluding projects at the extremes (0% or
# 100% of parent POM dependencies versioned)
model_parent_versioning_rate_trimmed <- fit_and_plot_trimmed(
  x_vals  = parent_pom_with_deps$parent_vers / parent_pom_with_deps$deps_parent_total,
  y_vals  = parent_pom_with_deps$gt_repr,
  degree  = 4,
  x_label = "Parent POM dependency versioning rate",
  title_label = "Quartic quasibinomial GLM: Parent POM versioning rate"
)

# --- 3.6 Square root of versioned dependency count in parent POM ---
# Predictor: raw count of explicitly versioned dependencies in the parent POM
# Subset: same as 3.5 (parent POM with at least one dependency)
# Transformation: square root (to reduce right skew)
# Polynomial degree: 4

model_parent_versioned_count <- fit_and_plot(
  x_vals  = sqrt(parent_pom_with_deps$parent_vers),
  y_vals  = parent_pom_with_deps$gt_repr,
  degree  = 4,
  x_label = "Square root of versioned dependency count (parent POM)",
  title_label = "Quartic quasibinomial GLM: Parent POM versioned count (sqrt)"
)

# =============================================================================
# SECTION 4: EXPLORATORY 3D PLOTS
# =============================================================================
# These were used during exploratory analysis to look for multivariate patterns.
# They are not intended as thesis figures.

# Project-wide versioning rate vs. log total dependency count vs. reproducibility
plot3d(
  maven_data$proj_vers / maven_data$deps_proj_total,
  log(maven_data$deps_proj_total),
  maven_data$gt_repr,
  xlab = "Versioning rate", ylab = "Log total deps", zlab = "Reproducibility",
  col = "blue", size = 5
)

# Same plot, color-coded by whether parent POM was found
# (blue = parent POM found, red = parent POM missing)
plot3d(
  maven_data$proj_vers / maven_data$deps_proj_total,
  log(maven_data$deps_proj_total),
  maven_data$gt_repr,
  xlab = "Versioning rate", ylab = "Log total deps", zlab = "Reproducibility",
  col = ifelse(maven_data$availability == 1, "red", "blue"), size = 5
)

# Log versioned count vs. log total dependency count vs. reproducibility
plot3d(
  log(maven_data$proj_vers),
  log(maven_data$deps_proj_total),
  maven_data$gt_repr,
  xlab = "Log versioned deps", ylab = "Log total deps", zlab = "Reproducibility",
  col = "blue", size = 5
)

# =============================================================================
# SECTION 5: POSTER/PRESENTATION VERSIONS OF KEY PLOTS
# =============================================================================
# Large font sizes for use in UHP/URC poster presentation.
# These use the same model as Section 3.1 (degree-4 versioning rate).

poster_theme <- theme_minimal() +
  theme(
    axis.title.x = element_text(size = 32, margin = margin(t = 30)),
    axis.title.y = element_text(size = 32, margin = margin(r = 30)),
    axis.text.x  = element_text(size = 28),
    axis.text.y  = element_text(size = 28)
  )

# Versioning rate — poster version
versioning_rate_df <- data.frame(
  x         = maven_data$proj_vers / maven_data$deps_proj_total,
  y         = maven_data$gt_repr,
  predicted = predict(model_versioning_rate, type = "response")
)

ggplot(versioning_rate_df, aes(x, y)) +
  geom_point(alpha = 0.5, size = 6) +
  geom_line(aes(y = predicted), color = "blue", size = 2) +
  labs(x = "% of dependencies with versions",
       y = "% of reproduction attempts matching") +
  poster_theme

# Commit age — poster version (scatter only, no fitted curve)
age_df <- data.frame(
  x = -maven_data$age,  # positive days elapsed
  y = maven_data$gt_repr
)

ggplot(age_df, aes(x, y)) +
  geom_point(alpha = 0.5, size = 6) +
  labs(x = "Time since commit was created (days)",
       y = "% of reproduction attempts matching") +
  poster_theme
