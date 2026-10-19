################################################################################
################################################################################

# Fixel-wise Bayes Factor Analysis

# Author: Jackson M. Lee
# Last updated: 19.10.2026

# Description:
# This script computes a Bayes factor (BF) per fixel from a t-statistic map.
# The output is a log10(BF) NIfTI map that can be visualised in MRview.
# The script can be used with any fixel metric (e.g., FDC, FC, or logFC).

# Before running:
#   1. Convert the t-statistic map from .mif to .nii/.nii.gz using mrconvert.
#   2. Set the input file path and output directory below.
#   3. Enter the sample sizes used in the original statistical analysis.
#   4. Set the contrast name and fixel metric for file naming.
#   5. Check the direction of the alternative hypothesis against the group ordering and sign convention of the original t-statistic map.
#   6. Ensure the input image elements correspond to the intended fixels.

# Important:
#   - The sample sizes must match the original analysis.
#   - c(0, Inf) specifies a positive-direction alternative; c(-Inf, 0) specifies a negative-direction alternative.
#   - NULL specifies a two-sided alternative.
#   - Positive log10(BF) favours H1; negative log10(BF) favours H0.
#   - The direction of H1 must be verified from the original t-statistics.
#   - The default prior is the medium Cauchy prior (r = 0.707).
#   - Sensitivity analyses compare the results using narrower and wider priors.
#   - Evidence categories follow the thresholds described by Lee & Wagenmakers (2018); cite the relevant source in any publication.
#   - This script assumes the t-statistic map has the correct dimensions and element ordering for the intended fixel-wise analysis.

# Outputs:
#   - A log10(BF) NIfTI map.
#   - A table counting fixels in each evidence category.
#   - Sensitivity analyses assessing the robustness of the results.

################################################################################
################################################################################

################################################################################
# Step 0: Set analysis parameters and load libraries
################################################################################

library(BayesFactor)
library(RNifti)

metric <- "FDC"             # e.g., FDC, FC, log-FC
contrast <- "EP_vs_HC"      # e.g., EP_vs_HC

input_file <- "path/to/t_statistic_map.nii.gz"
output_dir <- "path/to/output"

# Sample sizes from the original statistical analysis (chang these)
n1 <- 52
n2 <- 99

# Direction of the alternative hypothesis:
alternative <- c(0, Inf)

# Prior settings
primary_rscale <- "medium"  # Default Cauchy scale: approximately 0.707
sensitivity_rscales <- c(0.5, 1.0)

################################################################################
# Step 1: Load and prepare the t-statistic map
################################################################################

if (!file.exists(input_file)) {
  stop("Input file not found. Please check input_file.")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

tnii <- readNifti(input_file)
tvals <- as.numeric(tnii)

################################################################################
# Step 2: Compute Bayes factor per fixel
################################################################################

# BayesFactor returns log(BF) when extractBF(..., logbf = TRUE) is used.
# Dividing by log(10) converts natural logarithms to log10(BF).

calculate_log10bf <- function(tvals, rscale) {
  
  sapply(tvals, function(tval) {
    
    if (!is.finite(tval)) return(NA_real_)
    
    bf <- if (is.null(alternative)) {
      ttest.tstat(
        t = tval, n1 = n1, n2 = n2, rscale = rscale
      )
    } else {
      ttest.tstat(
        t = tval, n1 = n1, n2 = n2,
        nullInterval = alternative, rscale = rscale
      )
    }
    
    as.numeric(extractBF(bf, logbf = TRUE)$bf[1]) / log(10)
    
  })
}

# Primary analysis
log10bf <- calculate_log10bf(tvals, primary_rscale)

summary(log10bf)

################################################################################
# Step 3: Categorise evidence strength (Lee & Wagenmakers, 2018)
################################################################################

# Evidence thresholds expressed as log10(BF).
# These correspond approximately to BF values of 1, 3, 10, 30, and 100.

breaks_log10 <- c(
  -Inf, -2, -1.4771213, -1, -0.4771213, 0, 0.4771213, 1, 1.4771213, 2, Inf
)

labels_bf <- c(
  "Extreme H0", "Very strong H0", "Strong H0", "Moderate H0",
  "Anecdotal H0", "Anecdotal H1", "Moderate H1", "Strong H1",
  "Very strong H1", "Extreme H1"
)

categories <- cut(
  log10bf,
  breaks = breaks_log10,
  labels = labels_bf,
  include.lowest = TRUE
)

# Count fixels in each evidence category.
# Percentages are calculated over valid Bayes factors only.

fixel_counts <- as.data.frame(table(categories))
fixel_counts$Percentage <- round(
  100 * fixel_counts$Freq / sum(fixel_counts$Freq), 2
)

fixel_counts

################################################################################
# Step 4: Export log10(BF) map as NIfTI
################################################################################

# Check that the number of Bayes factors matches the input image.
if (length(log10bf) != length(tnii)) {
  stop("Mismatch between t-statistics and Bayes factors.")
}

# Reconstruct the array using the original image dimensions and reference.
log10bf_array <- array(log10bf, dim = dim(tnii))
log10bf_nii <- asNifti(log10bf_array, reference = tnii)

output_file <- file.path(
  output_dir,
  paste0("log10BF_", contrast, "_", metric, ".nii.gz")
)

writeNifti(log10bf_nii, file = output_file, version = 2)

################################################################################
# Step 5: Prior sensitivity analysis
################################################################################

# Recalculate log10(BF) using narrower and wider prior scales.
# All other settings remain unchanged.

log10bf_narrow <- calculate_log10bf(tvals, sensitivity_rscales[1])
log10bf_wide   <- calculate_log10bf(tvals, sensitivity_rscales[2])

################################################################################
# 5.1 Directional consistency
################################################################################

# Percentage of fixels for which the sign of log10(BF) is unchanged.
# This measures evidence direction (H0 vs H1), not the direction of the underlying group difference.

direction_narrow <- sign(log10bf_narrow) == sign(log10bf)
direction_wide   <- sign(log10bf_wide) == sign(log10bf)

directional_consistency <- data.frame(
  Comparison = c("Primary vs narrow", "Primary vs wide"),
  Percent_consistent = round(c(
    mean(direction_narrow, na.rm = TRUE) * 100,
    mean(direction_wide, na.rm = TRUE) * 100
  ), 2)
)

directional_consistency

################################################################################
# 5.2 Categorical stability
################################################################################

# Percentage of fixels remaining in the same evidence category.

categories_narrow <- cut(
  log10bf_narrow, breaks = breaks_log10,
  labels = labels_bf, include.lowest = TRUE
)

categories_wide <- cut(
  log10bf_wide, breaks = breaks_log10,
  labels = labels_bf, include.lowest = TRUE
)

categorical_stability <- data.frame(
  Comparison = c("Primary vs narrow", "Primary vs wide"),
  Percent_same_category = round(c(
    mean(categories_narrow == categories, na.rm = TRUE) * 100,
    mean(categories_wide == categories, na.rm = TRUE) * 100
  ), 2)
)

categorical_stability

################################################################################
# 5.3 Direction and magnitude of category transitions
################################################################################

# Category ranks:
# 1 = Extreme H0 ... 10 = Extreme H1

# Positive transitions indicate movement towards H1; negative transitions indicate movement towards H0.

rank_primary <- as.numeric(categories)
rank_narrow  <- as.numeric(categories_narrow)
rank_wide    <- as.numeric(categories_wide)

transition_narrow <- rank_narrow - rank_primary
transition_wide   <- rank_wide - rank_primary

transition_summary <- data.frame(
  Comparison = c("Narrow vs primary", "Wide vs primary"),
  Mean_transition = round(c(
    mean(transition_narrow, na.rm = TRUE),
    mean(transition_wide, na.rm = TRUE)
  ), 3),
  Median_transition = round(c(
    median(transition_narrow, na.rm = TRUE),
    median(transition_wide, na.rm = TRUE)
  ), 3),
  Mean_absolute_transition = round(c(
    mean(abs(transition_narrow), na.rm = TRUE),
    mean(abs(transition_wide), na.rm = TRUE)
  ), 3)
)

transition_summary

################################################################################
################################################################################
