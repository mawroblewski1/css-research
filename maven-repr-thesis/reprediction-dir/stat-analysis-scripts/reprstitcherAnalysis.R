# Statistical Analysis of Reproducibility Predictors
# Data: reprstitcherOutputFeb8.csv (output of reprstitcher.py)
#
# To run the full script from the RStudio console:
#   source("reprstitcherAnalysis_clean.R", echo = TRUE)
#
# Required packages (install once from the R console if needed):
#   install.packages("ggplot2")
#   install.packages("rgl")

# =============================================================================
# ANALYSIS OPTIONS — edit these before sourcing
# =============================================================================

# Set to TRUE to overlay a loess smoother (red curve) on every 2D plot,
# alongside the quasibinomial GLM fit (blue curve). Useful for diagnosing
# whether the GLM curve is being driven by influential points at the edges
# of the data. Set to FALSE for the standard GLM-only plots.
show_loess <- FALSE

# =============================================================================

library(ggplot2)  # 2D plotting
library(rgl)      # 3D plotting (used in exploratory section at bottom)

# Ensure rgl renders as widgets in the RStudio Viewer pane
options(rgl.printRglwidget = TRUE)

# =============================================================================
# SECTION 1: DATA LOADING AND FILTERING
# =============================================================================

# Load the stitched output (predictor metrics + ground truth reproducibility)
all_data <- read.csv("reprstitcherOutputFeb8.csv")
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

# Helper function: fit, plot, and summarize a quasibinomial polynomial model.
# If show_loess is TRUE (set at the top of the script), a loess smoother
# is overlaid in red for comparison with the GLM fit in blue.
fit_and_plot <- function(x_vals, y_vals, degree, x_label, title_label) {
  df <- data.frame(x = x_vals, y = y_vals)
  model <- glm(y ~ poly(x, degree), data = df, family = quasibinomial)
  df$predicted <- predict(model, type = "response")

  p <- ggplot(df, aes(x, y)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = predicted), color = "blue", linewidth = 1) +
    labs(title = title_label,
         x = x_label,
         y = "Predicted probability of reproducibility") +
    theme_minimal()

  if (show_loess) {
    p <- p + geom_smooth(method = "loess", color = "red", se = FALSE,
                         linewidth = 1, linetype = "dashed")
  }

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

# Companion plot: excludes the single leftmost point (smallest log versioned
# count). With n=50, one extreme point can disproportionately pull a quartic
# fit upward or downward near that edge. This checks whether the fitted
# curve's shape near the left edge depends on that one observation.
log_versioned_x <- log(maven_data$proj_vers)
leftmost_x <- min(log_versioned_x)
cat("\n[Trimmed variant: Log versioned dependency count]\n")
cat("Excluding leftmost point(s) at x =", leftmost_x,
    "(raw versioned count =", maven_data$proj_vers[which.min(log_versioned_x)], ")\n")
cat("Corresponding y (reproducibility) value(s):",
    maven_data$gt_repr[log_versioned_x == leftmost_x], "\n")

keep_log_versioned <- log_versioned_x != leftmost_x
model_log_versioned_count_trimmed <- fit_and_plot(
  x_vals  = log_versioned_x[keep_log_versioned],
  y_vals  = maven_data$gt_repr[keep_log_versioned],
  degree  = 4,
  x_label = "Log of versioned dependency count (project-wide)",
  title_label = paste0("Quartic quasibinomial GLM: Log versioned dependency count ",
                        "(leftmost point excluded, n=", sum(keep_log_versioned), ")")
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

# --- 3.7 Log of total project-wide dependency count ---
# Predictor: total number of dependencies declared across all POM files,
# including duplicates across modules, regardless of whether they are versioned.
# This tests whether the sheer scale of a project's dependency footprint
# predicts reproducibility, independent of versioning practice.
# Transformation: log (to reduce right skew)
# Polynomial degree: 4

model_log_total_deps <- fit_and_plot(
  x_vals  = log(maven_data$deps_proj_total),
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Log of total dependency count (project-wide, including duplicates)",
  title_label = "Quartic quasibinomial GLM: Log total dependency count"
)

# Companion plot: exclude the leftmost point (smallest total dependency count),
# as a single extreme point may disproportionately influence the quartic fit.
log_total_x <- log(maven_data$deps_proj_total)
leftmost_total_x <- min(log_total_x)
cat("\n[Trimmed variant: Log total dependency count]\n")
cat("Excluding leftmost point(s) at x =", leftmost_total_x,
    "(raw total count =", maven_data$deps_proj_total[which.min(log_total_x)], ")\n")
cat("Corresponding y (reproducibility) value(s):",
    maven_data$gt_repr[log_total_x == leftmost_total_x], "\n")

keep_log_total <- log_total_x != leftmost_total_x
model_log_total_deps_trimmed <- fit_and_plot(
  x_vals      = log_total_x[keep_log_total],
  y_vals      = maven_data$gt_repr[keep_log_total],
  degree      = 4,
  x_label     = "Log of total dependency count (project-wide, including duplicates)",
  title_label = paste0("Quartic quasibinomial GLM: Log total dependency count ",
                       "(leftmost point excluded, n=", sum(keep_log_total), ")")
)

# --- 3.8 Log of total file count in the repository ---
# Predictor: total number of files in the repository, a general proxy for
# project size independent of Maven structure.
# Transformation: log (to reduce right skew)
# Polynomial degree: 4

model_log_num_files <- fit_and_plot(
  x_vals  = log(maven_data$num_files),
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Log of total file count in repository",
  title_label = "Quartic quasibinomial GLM: Log total file count"
)

# --- 3.9 Log of unique named dependency count (project-wide) ---
# Predictor: number of distinct dependencies by name (groupId/artifactId),
# deduplicating across all POM files in the project. Contrasts with
# deps_proj_total (which counts duplicates) to isolate the effect of
# the number of distinct external dependencies rather than total declarations.
# Transformation: log (to reduce right skew)
# Polynomial degree: 4

model_log_proj_unique <- fit_and_plot(
  x_vals  = log(maven_data$proj_unique),
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Log of unique named dependency count (project-wide)",
  title_label = "Quartic quasibinomial GLM: Log unique dependency count"
)

# --- 3.10 Dependency redundancy ratio (proj_unique / deps_proj_total) ---
# Predictor: ratio of unique named dependencies to total dependency declarations.
# A ratio near 1 means few duplicates across modules (each dependency declared
# roughly once); a ratio near 0 means heavy duplication across submodules.
# Transformation: none (already a ratio in [0, 1])
# Polynomial degree: 4
# Trimmed companions: exclude extremes at 0 and 1

model_redundancy_ratio <- fit_and_plot(
  x_vals  = maven_data$proj_unique / maven_data$deps_proj_total,
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Dependency redundancy ratio (unique / total declarations)",
  title_label = "Quartic quasibinomial GLM: Dependency redundancy ratio"
)

model_redundancy_ratio_trimmed <- fit_and_plot_trimmed(
  x_vals  = maven_data$proj_unique / maven_data$deps_proj_total,
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Dependency redundancy ratio (unique / total declarations)",
  title_label = "Quartic quasibinomial GLM: Dependency redundancy ratio"
)

# --- 3.11 Versioning rate among unique dependencies (proj_vers / proj_unique) ---
# Predictor: proportion of uniquely named dependencies that are explicitly
# versioned somewhere in the project. Strips out the effect of multi-module
# duplication present in PRVRATE (proj_vers / deps_proj_total).
# Transformation: none (ratio in [0, 1])
# Polynomial degree: 4
# Trimmed companions: exclude extremes at 0 and 1

model_unique_versioning_rate <- fit_and_plot(
  x_vals  = maven_data$proj_vers / maven_data$proj_unique,
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Versioning rate among unique dependencies (versioned / unique)",
  title_label = "Quartic quasibinomial GLM: Unique dependency versioning rate"
)

model_unique_versioning_rate_trimmed <- fit_and_plot_trimmed(
  x_vals  = maven_data$proj_vers / maven_data$proj_unique,
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Versioning rate among unique dependencies (versioned / unique)",
  title_label = "Quartic quasibinomial GLM: Unique dependency versioning rate"
)

# --- 3.12 POM density (num_poms / num_files) ---
# Predictor: proportion of all files in the repository that are POM files.
# A proxy for how Maven-centric the project structure is — a high ratio
# suggests a heavily multi-module Maven project.
# Transformation: none (ratio, but unlikely to hit exactly 0 or 1)
# Polynomial degree: 4

model_pom_density <- fit_and_plot(
  x_vals  = maven_data$num_poms / maven_data$num_files,
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "POM file density (POM files / total files)",
  title_label = "Quartic quasibinomial GLM: POM file density"
)

# --- 3.13 Average dependencies per POM file (deps_proj_total / num_poms) ---
# Predictor: mean number of dependency declarations per POM file across the
# project. Captures dependency density at the file level rather than total
# scale — distinguishes a project with many dependencies spread across many
# modules from one with many dependencies concentrated in few files.
# Transformation: log (to reduce right skew)
# Polynomial degree: 4

model_log_deps_per_pom <- fit_and_plot(
  x_vals  = log(maven_data$deps_proj_total / maven_data$num_poms),
  y_vals  = maven_data$gt_repr,
  degree  = 4,
  x_label = "Log of average dependency count per POM file",
  title_label = "Quartic quasibinomial GLM: Log average deps per POM"
)

# =============================================================================
# SECTION 3 SUMMARY: ALL MODEL SUMMARIES AND MINIMUM P-VALUES
# =============================================================================
# Prints every model summary in one block, and extracts the minimum p-value
# among non-intercept coefficients for each model.
# This section is designed to be robust to source() truncation — all key
# numbers appear together at the end of Section 3 regardless of earlier output.

# Helper: extract minimum non-intercept p-value from a glm summary
min_pval <- function(model) {
  coef_table <- summary(model)$coefficients
  # Exclude intercept (row 1), extract p-values (column 4)
  pvals <- coef_table[-1, 4]
  return(min(pvals))
}

# List of all 2D models with descriptive labels
model_list <- list(
  list(label = "3.1  PRVRATE, degree-4 (full data)",
       model = model_versioning_rate),
  list(label = "3.1T PRVRATE, degree-4 (extremes excluded)",
       model = model_versioning_rate_trimmed),
  list(label = "3.1b PRVRATE, degree-6 (full data)",
       model = model_versioning_rate_deg6),
  list(label = "3.1bT PRVRATE, degree-6 (extremes excluded)",
       model = model_versioning_rate_deg6_trimmed),
  list(label = "3.2  Log versioned dependency count (full data)",
       model = model_log_versioned_count),
  list(label = "3.2T Log versioned dependency count (leftmost excluded)",
       model = model_log_versioned_count_trimmed),
  list(label = "3.3  Commit age",
       model = model_age),
  list(label = "3.4  Log POM file count",
       model = model_log_pom_count),
  list(label = "3.5  Parent POM versioning rate (full data)",
       model = model_parent_versioning_rate),
  list(label = "3.5T Parent POM versioning rate (extremes excluded)",
       model = model_parent_versioning_rate_trimmed),
  list(label = "3.6  Parent POM versioned count (sqrt)",
       model = model_parent_versioned_count),
  list(label = "3.7  Log total dependency count (full data)",
       model = model_log_total_deps),
  list(label = "3.7T Log total dependency count (leftmost excluded)",
       model = model_log_total_deps_trimmed),
  list(label = "3.8  Log total file count",
       model = model_log_num_files),
  list(label = "3.9  Log unique dependency count",
       model = model_log_proj_unique),
  list(label = "3.10  Dependency redundancy ratio (full data)",
       model = model_redundancy_ratio),
  list(label = "3.10T Dependency redundancy ratio (extremes excluded)",
       model = model_redundancy_ratio_trimmed),
  list(label = "3.11  Unique dependency versioning rate (full data)",
       model = model_unique_versioning_rate),
  list(label = "3.11T Unique dependency versioning rate (extremes excluded)",
       model = model_unique_versioning_rate_trimmed),
  list(label = "3.12  POM file density",
       model = model_pom_density),
  list(label = "3.13  Log average deps per POM file",
       model = model_log_deps_per_pom)
)

cat("\n")
cat("=============================================================================\n")
cat("SECTION 3 SUMMARY: MINIMUM NON-INTERCEPT P-VALUE PER MODEL\n")
cat("=============================================================================\n")
cat(sprintf("%-55s %s\n", "Model", "Min p-value (non-intercept)"))
cat(strrep("-", 75), "\n")
for (m in model_list) {
  pval <- tryCatch(min_pval(m$model), error = function(e) NA)
  cat(sprintf("%-55s %.4f\n", m$label, pval))
}
cat(strrep("=", 75), "\n")

cat("\n")
cat("=============================================================================\n")
cat("SECTION 3 SUMMARY: FULL MODEL SUMMARIES\n")
cat("=============================================================================\n")
for (m in model_list) {
  cat("\n---", m$label, "---\n")
  print(summary(m$model))
}

# =============================================================================
# SECTION 4: EXPLORATORY 3D PLOTS
# =============================================================================
# The heuristic for these plots:
#   X axis — structural vulnerability: how exposed is the project to other
#             programmers' changes? (versioning rate or dependency count)
#   Y axis — time: the medium through which collective effects accumulate.
#             Other programmers update libraries over time; longer elapsed
#             time means more opportunity for the ecosystem to shift.
#   Z axis — reproducibility: the observable outcome.
#
# If the collective effect is real, projects with higher structural
# vulnerability should show steeper reproducibility decline with age
# than low-vulnerability projects. A flat surface along the age axis
# is consistent with the 1D age finding (no significant effect in 40 days).
#
# Note: rgl renders interactive 3D windows. Use rgl::snapshot3d() to
# capture static images for the thesis from a good viewing angle.
#
# The fitted surfaces use a 2D quasibinomial GLM with degree-2 polynomial
# terms and their interaction. With n=50 this is statistically fragile;
# treat the surfaces as illustrative rather than confirmatory.

# Helper: fit a 3D quasibinomial GLM (model fitting only, no plotting)
# This is called first so model objects always exist even if rgl plotting fails.
fit_3d <- function(x_vals, y_vals, z_vals, title_label) {
  cat("\n--- 3D model fit:", title_label, "---\n")
  df3d <- data.frame(x = x_vals, y = y_vals, z = z_vals)
  model_3d <- glm(z ~ poly(x, 2) * poly(y, 2), data = df3d, family = quasibinomial)
  print(summary(model_3d))
  return(model_3d)
}

# Helper: render a 3D scatter + fitted surface for an already-fitted model.
# Call open3d() before calling this function.
plot_3d <- function(model_3d, x_vals, y_vals, z_vals,
                    xlab, ylab, zlab, title_label) {
  tryCatch({
    plot3d(x_vals, y_vals, z_vals,
           xlab = xlab, ylab = ylab, zlab = zlab,
           col = "black", size = 5, alpha = 0.7)
    title3d(main = title_label, cex = 0.9)

    x_grid <- seq(min(x_vals), max(x_vals), length.out = 40)
    y_grid <- seq(min(y_vals), max(y_vals), length.out = 40)
    grid_df <- expand.grid(x = x_grid, y = y_grid)
    grid_df$z_pred <- predict(model_3d, newdata = grid_df, type = "response")
    z_matrix <- matrix(grid_df$z_pred, nrow = 40, ncol = 40)
    surface3d(x_grid, y_grid, z_matrix, col = "blue", alpha = 0.4, back = "lines")
    rglwidget()
  }, error = function(e) {
    cat("Note: rgl plotting failed for", title_label, "—", conditionMessage(e), "\n")
  })
}

# --- 4.1 Existing plots (from original exploratory analysis) ---

# Project-wide versioning rate vs. log total dependency count vs. reproducibility
open3d()
plot3d(
  maven_data$proj_vers / maven_data$deps_proj_total,
  log(maven_data$deps_proj_total),
  maven_data$gt_repr,
  xlab = "Versioning rate", ylab = "Log total deps", zlab = "Reproducibility",
  col = "blue", size = 5
)
rglwidget()

# Same plot, color-coded by whether parent POM was found
# (blue = parent POM found, red = parent POM missing)
open3d()
plot3d(
  maven_data$proj_vers / maven_data$deps_proj_total,
  log(maven_data$deps_proj_total),
  maven_data$gt_repr,
  xlab = "Versioning rate", ylab = "Log total deps", zlab = "Reproducibility",
  col = ifelse(maven_data$availability == 1, "red", "blue"), size = 5
)
rglwidget()

# Log versioned count vs. log total dependency count vs. reproducibility
open3d()
plot3d(
  log(maven_data$proj_vers),
  log(maven_data$deps_proj_total),
  maven_data$gt_repr,
  xlab = "Log versioned deps", ylab = "Log total deps", zlab = "Reproducibility",
  col = "blue", size = 5
)
rglwidget()

# --- 4.2 New plots: structural vulnerability × time × reproducibility ---
# Models are fitted first (fit_3d), then plotted (plot_3d).
# This ensures model objects exist for the summary block even if rgl fails.

# Project-wide versioning rate × commit age × reproducibility
model_3d_prvrate_age <- fit_3d(
  x_vals = maven_data$proj_vers / maven_data$deps_proj_total,
  y_vals = -maven_data$age,
  z_vals = maven_data$gt_repr,
  title_label = "Versioning rate x Commit age x Reproducibility"
)
open3d()
plot_3d(model_3d_prvrate_age,
  x_vals = maven_data$proj_vers / maven_data$deps_proj_total,
  y_vals = -maven_data$age,
  z_vals = maven_data$gt_repr,
  xlab = "Project-wide versioning rate",
  ylab = "Days since failing commit",
  zlab = "Reproducibility",
  title_label = "Versioning rate x Commit age x Reproducibility"
)

# Log total dependency count × commit age × reproducibility
model_3d_totaldeps_age <- fit_3d(
  x_vals = log(maven_data$deps_proj_total),
  y_vals = -maven_data$age,
  z_vals = maven_data$gt_repr,
  title_label = "Log total dependency count x Commit age x Reproducibility"
)
open3d()
plot_3d(model_3d_totaldeps_age,
  x_vals = log(maven_data$deps_proj_total),
  y_vals = -maven_data$age,
  z_vals = maven_data$gt_repr,
  xlab = "Log total dependency count (project-wide)",
  ylab = "Days since failing commit",
  zlab = "Reproducibility",
  title_label = "Log total dependency count x Commit age x Reproducibility"
)

# Log versioned dependency count × commit age × reproducibility
model_3d_versioned_age <- fit_3d(
  x_vals = log(maven_data$proj_vers),
  y_vals = -maven_data$age,
  z_vals = maven_data$gt_repr,
  title_label = "Log versioned dependency count x Commit age x Reproducibility"
)
open3d()
plot_3d(model_3d_versioned_age,
  x_vals = log(maven_data$proj_vers),
  y_vals = -maven_data$age,
  z_vals = maven_data$gt_repr,
  xlab = "Log versioned dependency count (project-wide)",
  ylab = "Days since failing commit",
  zlab = "Reproducibility",
  title_label = "Log versioned dependency count x Commit age x Reproducibility"
)

# Log total dependency count × versioning rate × reproducibility
model_3d_totaldeps_prvrate <- fit_3d(
  x_vals = log(maven_data$deps_proj_total),
  y_vals = maven_data$proj_vers / maven_data$deps_proj_total,
  z_vals = maven_data$gt_repr,
  title_label = "Log total dependency count x Versioning rate x Reproducibility"
)
open3d()
plot_3d(model_3d_totaldeps_prvrate,
  x_vals = log(maven_data$deps_proj_total),
  y_vals = maven_data$proj_vers / maven_data$deps_proj_total,
  z_vals = maven_data$gt_repr,
  xlab = "Log total dependency count (project-wide)",
  ylab = "Project-wide versioning rate",
  zlab = "Reproducibility",
  title_label = "Log total dependency count x Versioning rate x Reproducibility"
)

# Number of POM files × versioning rate × reproducibility
model_3d_numpoms_prvrate <- fit_3d(
  x_vals = maven_data$num_poms,
  y_vals = maven_data$proj_vers / maven_data$deps_proj_total,
  z_vals = maven_data$gt_repr,
  title_label = "Number of POM files x Versioning rate x Reproducibility"
)
open3d()
plot_3d(model_3d_numpoms_prvrate,
  x_vals = maven_data$num_poms,
  y_vals = maven_data$proj_vers / maven_data$deps_proj_total,
  z_vals = maven_data$gt_repr,
  xlab = "Number of POM files",
  ylab = "Project-wide versioning rate",
  zlab = "Reproducibility",
  title_label = "Number of POM files x Versioning rate x Reproducibility"
)

# =============================================================================
# SECTION 4 SUMMARY: 3D MODEL SUMMARIES AND MINIMUM P-VALUES
# =============================================================================
# Placed here so all 3D model objects are guaranteed to exist.

model_list_3d <- list(
  list(label = "4.2a PRVRATE x Commit age",
       model = model_3d_prvrate_age),
  list(label = "4.2b Log total dependency count x Commit age",
       model = model_3d_totaldeps_age),
  list(label = "4.2c Log versioned dependency count x Commit age",
       model = model_3d_versioned_age),
  list(label = "4.2d Log total dependency count x PRVRATE",
       model = model_3d_totaldeps_prvrate),
  list(label = "4.2e Number of POM files x PRVRATE",
       model = model_3d_numpoms_prvrate)
)

cat("\n")
cat("=============================================================================\n")
cat("SECTION 4 SUMMARY: 3D MODEL MINIMUM P-VALUES\n")
cat("=============================================================================\n")
cat(sprintf("%-55s %s\n", "Model", "Min p-value (non-intercept)"))
cat(strrep("-", 75), "\n")
for (m in model_list_3d) {
  pval <- tryCatch(min_pval(m$model), error = function(e) NA)
  cat(sprintf("%-55s %.4f\n", m$label, pval))
}
cat(strrep("=", 75), "\n")

cat("\n")
cat("=============================================================================\n")
cat("SECTION 4 SUMMARY: FULL 3D MODEL SUMMARIES\n")
cat("=============================================================================\n")
for (m in model_list_3d) {
  cat("\n---", m$label, "---\n")
  print(summary(m$model))
}

# =============================================================================
# SECTION 5: POSTER/PRESENTATION VERSIONS OF KEY PLOTS
# =============================================================================
# The original poster plots (OPP 1 and OPP 2) are preserved untouched at the
# end of this section (5.4 and 5.5). New variants with standard analysis
# styling and toggleable loess are provided first (5.1, 5.2, 5.3).

poster_theme <- theme_minimal() +
  theme(
    axis.title.x = element_text(size = 32, margin = margin(t = 30)),
    axis.title.y = element_text(size = 32, margin = margin(r = 30)),
    axis.text.x  = element_text(size = 28),
    axis.text.y  = element_text(size = 28)
  )

# Precompute data frames used across multiple plots below
versioning_rate_df <- data.frame(
  x         = maven_data$proj_vers / maven_data$deps_proj_total,
  y         = maven_data$gt_repr,
  predicted = predict(model_versioning_rate, type = "response")
)

prvrate_all  <- maven_data$proj_vers / maven_data$deps_proj_total
prvrate_keep <- (prvrate_all != 0) & (prvrate_all != 1)
versioning_rate_trimmed_df <- data.frame(
  x         = prvrate_all[prvrate_keep],
  y         = maven_data$gt_repr[prvrate_keep],
  predicted = predict(model_versioning_rate_trimmed, type = "response")
)

age_df <- data.frame(
  x = -maven_data$age,
  y = maven_data$gt_repr
)

# --- 5.1 OPP 1 variant: PRVRATE, full data, standard styling, toggleable loess ---
p_51 <- ggplot(versioning_rate_df, aes(x, y)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = predicted), color = "blue", linewidth = 1) +
  labs(title = "Quartic quasibinomial GLM: Project-wide versioning rate",
       x = "Project-wide dependency versioning rate",
       y = "Predicted probability of reproducibility") +
  theme_minimal()
if (show_loess) {
  p_51 <- p_51 + geom_smooth(method = "loess", color = "red", se = FALSE,
                              linewidth = 1, linetype = "dashed")
}
print(p_51)

# --- 5.2 OPP 1 variant: PRVRATE, extremes excluded, standard styling, toggleable loess ---
p_52 <- ggplot(versioning_rate_trimmed_df, aes(x, y)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = predicted), color = "blue", linewidth = 1) +
  labs(title = "Quartic quasibinomial GLM: Project-wide versioning rate (extremes excluded)",
       x = "Project-wide dependency versioning rate",
       y = "Predicted probability of reproducibility") +
  theme_minimal()
if (show_loess) {
  p_52 <- p_52 + geom_smooth(method = "loess", color = "red", se = FALSE,
                              linewidth = 1, linetype = "dashed")
}
print(p_52)

# --- 5.3 OPP 2 variant: commit age, standard styling, GLM curve, toggleable loess ---
age_df$predicted <- predict(model_age, type = "response")

p_53 <- ggplot(age_df, aes(x, y)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = predicted), color = "blue", linewidth = 1) +
  labs(title = "Commit age vs. reproducibility",
       x = "Days since failing commit was created",
       y = "Predicted probability of reproducibility") +
  theme_minimal()
if (show_loess) {
  p_53 <- p_53 + geom_smooth(method = "loess", color = "red", se = FALSE,
                              linewidth = 1, linetype = "dashed")
}
print(p_53)

# --- 5.4 OPP 1 original — untouched ---
ggplot(versioning_rate_df, aes(x, y)) +
  geom_point(alpha = 0.5, size = 6) +
  geom_line(aes(y = predicted), color = "blue", linewidth = 2) +
  labs(x = "% of dependencies with versions",
       y = "% of reproduction attempts matching") +
  poster_theme

# --- 5.5 OPP 2 original — untouched ---
ggplot(age_df, aes(x, y)) +
  geom_point(alpha = 0.5, size = 6) +
  labs(x = "Time since commit was created (days)",
       y = "% of reproduction attempts matching") +
  poster_theme

