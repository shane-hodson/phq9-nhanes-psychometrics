# ==============================================================================
# Stage 3: Dimensionality analysis
# Project: PHQ-9 psychometric evaluation using NHANES 2017–March 2020
#
# Current implementation scope:
#   - reconstruct the complete-PHQ-9 analytic sample with SEQN retained
#   - create the locked 50:50 development and validation split
#   - validate both split samples
#   - estimate and validate the development-sample polychoric matrix
#   - run the prespecified single-core ordinal parallel analysis
#   - reconstruct and validate the 95th-percentile factor recommendation
#   - save the recoverable split, polychoric and parallel-analysis objects
#
# Parallel analysis supports three factors.
# EFA and CFA are not yet run in this script.
# ==============================================================================

# 1. Install any missing packages ----------------------------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "here",
  "tibble",
  "psych"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}


# 2. Load packages --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(here)
library(tibble)
library(psych)

# 3. Define locked Stage 3 constants -------------------------------------------

phq9_items <- c(
  "DPQ010",
  "DPQ020",
  "DPQ030",
  "DPQ040",
  "DPQ050",
  "DPQ060",
  "DPQ070",
  "DPQ080",
  "DPQ090"
)

split_seed <- 20260723L
expected_complete_n <- 8276L
expected_development_n <- 4138L
expected_validation_n <- 4138L

# 4. Define input and output paths ---------------------------------------------

stage2a_data_path <- here(
  "data",
  "processed",
  "phq9_stage2a_data.rds"
)

stage3_split_path <- here(
  "data",
  "processed",
  "phq9_stage3_split.rds"
)

split_summary_path <- here(
  "tables",
  "phq9_stage3_split_summary.csv"
)

if (!file.exists(stage2a_data_path)) {
  stop(
    "The Stage 2A processed data file was not found at: ",
    stage2a_data_path
  )
}

# 5. Load and validate the source object ---------------------------------------

phq9_stage2a_data <- readRDS(stage2a_data_path)

required_variables <- c(
  "SEQN",
  phq9_items
)

missing_variables <- setdiff(
  required_variables,
  names(phq9_stage2a_data)
)

if (length(missing_variables) > 0) {
  stop(
    "The following required variables are missing: ",
    paste(missing_variables, collapse = ", ")
  )
}

if (anyNA(phq9_stage2a_data$SEQN)) {
  stop("SEQN contains missing values in the Stage 2A source object.")
}

if (anyDuplicated(phq9_stage2a_data$SEQN) > 0L) {
  stop("SEQN contains duplicate values in the Stage 2A source object.")
}

# 6. Reconstruct the complete-PHQ-9 Stage 3 sample -----------------------------

phq9_stage3_data <- phq9_stage2a_data |>
  filter(
    if_all(
      all_of(phq9_items),
      ~ !is.na(.x) & .x %in% 0:3
    )
  ) |>
  select(
    SEQN,
    all_of(phq9_items)
  ) |>
  arrange(SEQN)

if (nrow(phq9_stage3_data) != expected_complete_n) {
  stop(
    "Unexpected complete-PHQ-9 sample size. Expected ",
    expected_complete_n,
    " but found ",
    nrow(phq9_stage3_data),
    "."
  )
}

if (anyNA(phq9_stage3_data$SEQN)) {
  stop("SEQN contains missing values in the complete-PHQ-9 sample.")
}

if (anyDuplicated(phq9_stage3_data$SEQN) > 0L) {
  stop("SEQN contains duplicate values in the complete-PHQ-9 sample.")
}

complete_item_matrix <- as.matrix(
  phq9_stage3_data[phq9_items]
)

if (anyNA(complete_item_matrix)) {
  stop("The complete-PHQ-9 item matrix contains missing values.")
}

if (min(complete_item_matrix) < 0 || max(complete_item_matrix) > 3) {
  stop("The complete-PHQ-9 item matrix contains responses outside 0 to 3.")
}

if (any(complete_item_matrix != floor(complete_item_matrix))) {
  stop("The complete-PHQ-9 item matrix contains non-integer responses.")
}

# 7. Create the locked development and validation split ------------------------

set.seed(split_seed)

development_rows <- sample.int(
  n = nrow(phq9_stage3_data),
  size = expected_development_n,
  replace = FALSE
)

phq9_stage3_split <- phq9_stage3_data
phq9_stage3_split$stage3_sample <- "validation"
phq9_stage3_split$stage3_sample[development_rows] <- "development"

phq9_stage3_split <- phq9_stage3_split |>
  relocate(stage3_sample, .after = SEQN)

actual_split_counts <- table(phq9_stage3_split$stage3_sample)

if (as.integer(actual_split_counts["development"]) != expected_development_n) {
  stop("Unexpected development-sample size.")
}

if (as.integer(actual_split_counts["validation"]) != expected_validation_n) {
  stop("Unexpected validation-sample size.")
}

# 8. Define split-sample validation summary ------------------------------------

validate_split_sample <- function(data, sample_name) {
  sample_data <- data |>
    filter(stage3_sample == sample_name)

  item_matrix <- as.matrix(sample_data[phq9_items])
  total_scores <- rowSums(item_matrix)

  tibble(
    stage3_sample = sample_name,
    participant_n = nrow(sample_data),
    unique_seqn_n = n_distinct(sample_data$SEQN),
    duplicated_seqn_n = sum(duplicated(sample_data$SEQN)),
    missing_seqn_n = sum(is.na(sample_data$SEQN)),
    missing_item_value_n = sum(is.na(item_matrix)),
    minimum_item_response = min(item_matrix),
    maximum_item_response = max(item_matrix),

    noninteger_item_value_n = sum(item_matrix != floor(item_matrix)),
    phq9_total_mean = mean(total_scores),
    phq9_total_sd = sd(total_scores),
    phq9_total_median = median(total_scores),
    phq9_total_minimum = min(total_scores),
    phq9_total_maximum = max(total_scores)
  )
}

split_summary <- bind_rows(
  validate_split_sample(
    phq9_stage3_split,
    "development"
  ),
  validate_split_sample(
    phq9_stage3_split,
    "validation"
  )
)

if (any(split_summary$participant_n != 4138L)) {
  stop("At least one split sample does not contain 4,138 participants.")
}

if (any(split_summary$participant_n != split_summary$unique_seqn_n)) {
  stop("At least one split sample contains non-unique SEQN values.")
}

if (any(split_summary$duplicated_seqn_n != 0L)) {
  stop("At least one split sample contains duplicated SEQN values.")
}

if (any(split_summary$missing_seqn_n != 0L)) {
  stop("At least one split sample contains missing SEQN values.")
}

if (any(split_summary$missing_item_value_n != 0L)) {
  stop("At least one split sample contains missing PHQ-9 item values.")
}

if (
  any(split_summary$minimum_item_response != 0) ||
    any(split_summary$maximum_item_response != 3)
) {
  stop("At least one split sample has an unexpected item-response range.")
}

if (any(split_summary$noninteger_item_value_n != 0L)) {
  stop("At least one split sample contains non-integer item responses.")
}

# 9. Check item-category availability and floor effects ------------------------

item_category_frequencies <- phq9_stage3_split |>
  pivot_longer(
    cols = all_of(phq9_items),
    names_to = "variable",
    values_to = "response"
  ) |>
  count(
    stage3_sample,
    variable,
    response,
    name = "frequency"
  ) |>
  complete(
    stage3_sample,
    variable,
    response = 0:3,
    fill = list(frequency = 0L)
  ) |>
  arrange(
    stage3_sample,
    match(variable, phq9_items),
    response
  )

absent_categories <- item_category_frequencies |>
  filter(frequency == 0L)

if (nrow(absent_categories) > 0L) {
  print(absent_categories)

  stop(
    "At least one PHQ-9 response category is absent from a split sample."
  )
}

item_floor_summary <- item_category_frequencies |>
  group_by(
    stage3_sample,
    variable
  ) |>
  summarise(
    participant_n = sum(frequency),
    floor_frequency = frequency[response == 0],
    floor_percentage = 100 * floor_frequency / participant_n,
    .groups = "drop"
  ) |>
  arrange(
    stage3_sample,
    match(variable, phq9_items)
  )

# 10. Prepare development-sample item data -------------------------------------

development_items <- phq9_stage3_split |>
  filter(stage3_sample == "development") |>
  select(all_of(phq9_items)) |>
  as.data.frame()

# 11. Estimate the development polychoric correlation matrix -------------------

polychoric_warnings <- character()

development_poly <- withCallingHandlers(
  psych::polychoric(
    development_items,
    correct = 0.5,
    smooth = FALSE,
    global = TRUE,
    progress = FALSE
  ),
  warning = function(warning_condition) {
    polychoric_warnings <<- c(
      polychoric_warnings,
      conditionMessage(warning_condition)
    )

    invokeRestart("muffleWarning")
  }
)

development_poly_matrix <- development_poly$rho
development_poly_thresholds <- development_poly$tau

if (length(polychoric_warnings) > 0L) {
  stop(
    "Development polychoric estimation produced warning(s): ",
    paste(unique(polychoric_warnings), collapse = " | ")
  )
}

if (!identical(dim(development_poly_matrix), c(9L, 9L))) {
  stop("The development polychoric matrix is not 9 by 9.")
}

if (
  !identical(rownames(development_poly_matrix), phq9_items) ||
  !identical(colnames(development_poly_matrix), phq9_items)
) {
  stop("The development polychoric matrix has an unexpected item order.")
}

if (any(!is.finite(development_poly_matrix))) {
  stop("The development polychoric matrix contains non-finite estimates.")
}

if (
  !isTRUE(
    all.equal(
      development_poly_matrix,
      t(development_poly_matrix),
      tolerance = 1e-12
    )
  )
) {
  stop("The development polychoric matrix is not symmetric.")
}

if (any(abs(diag(development_poly_matrix) - 1) > 1e-12)) {
  stop("The development polychoric matrix does not have a unit diagonal.")
}

if (
  min(development_poly_matrix) < -1 ||
  max(development_poly_matrix) > 1
) {
  stop("The development polychoric matrix contains values outside -1 to 1.")
}

if (!identical(dim(development_poly_thresholds), c(9L, 3L))) {
  stop("The development threshold matrix is not 9 by 3.")
}

if (any(!is.finite(development_poly_thresholds))) {
  stop("The development threshold matrix contains non-finite estimates.")
}

if (!identical(rownames(development_poly_thresholds), phq9_items)) {
  stop("The development threshold matrix has an unexpected item order.")
}

if (!identical(colnames(development_poly_thresholds), c("1", "2", "3"))) {
  stop("The development threshold matrix has unexpected column names.")
}

development_thresholds_increasing <- apply(
  development_poly_thresholds,
  1,
  function(item_thresholds) {
    all(diff(item_thresholds) > 0)
  }
)

if (!all(development_thresholds_increasing)) {
  stop(
    "Thresholds are not strictly increasing for: ",
    paste(
      phq9_items[!development_thresholds_increasing],
      collapse = ", "
    )
  )
}

development_poly_eigenvalues <- eigen(
  development_poly_matrix,
  symmetric = TRUE,
  only.values = TRUE
)$values

development_poly_minimum_eigenvalue <- min(
  development_poly_eigenvalues
)

if (development_poly_minimum_eigenvalue <= 0) {
  stop(
    "The development polychoric matrix is not positive definite. ",
    "Minimum eigenvalue: ",
    development_poly_minimum_eigenvalue
  )
}

development_poly_off_diagonal <- development_poly_matrix[
  lower.tri(development_poly_matrix)
]

development_poly_diagnostics <- tibble(
  sample_n = nrow(development_items),
  item_n = ncol(development_items),
  threshold_n_per_item = ncol(development_poly_thresholds),
  warning_n = length(polychoric_warnings),
  minimum_correlation = min(development_poly_off_diagonal),
  maximum_correlation = max(development_poly_off_diagonal),
  minimum_eigenvalue = development_poly_minimum_eigenvalue,
  positive_definite = development_poly_minimum_eigenvalue > 0,
  smoothing_requested = FALSE,
  continuity_correction = 0.5,
  global_thresholds = TRUE
)

# 12. Run the prespecified ordinal parallel analysis ---------------------------

parallel_analysis_seed <- 20260724L
parallel_analysis_iterations <- 1000L
parallel_analysis_quantile <- 0.95

parallel_analysis_warnings <- character()

parallel_results <- local({
  previous_mc_cores <- getOption("mc.cores")

  on.exit(
    options(mc.cores = previous_mc_cores),
    add = TRUE
  )

  options(mc.cores = 1L)
  set.seed(parallel_analysis_seed)

  withCallingHandlers(
    psych::fa.parallel(
      x = development_items,
      fm = "minres",
      fa = "fa",
      nfactors = 1,
      n.iter = parallel_analysis_iterations,
      SMC = FALSE,
      sim = FALSE,
      quant = parallel_analysis_quantile,
      cor = "poly",
      correct = 0.5,
      plot = FALSE
    ),
    warning = function(warning_condition) {
      parallel_analysis_warnings <<- c(
        parallel_analysis_warnings,
        conditionMessage(warning_condition)
      )

      invokeRestart("muffleWarning")
    }
  )
})

# 13. Validate and reconstruct the parallel-analysis recommendation ------------

parallel_factor_columns <- paste0(
  "F",
  seq_along(phq9_items)
)

if (length(parallel_analysis_warnings) > 0L) {
  stop(
    "Parallel analysis produced warning(s): ",
    paste(
      unique(parallel_analysis_warnings),
      collapse = " | "
    )
  )
}

if (!is.matrix(parallel_results$values)) {
  stop("The parallel-analysis simulation results are not a matrix.")
}

if (nrow(parallel_results$values) != parallel_analysis_iterations) {
  stop(
    "The parallel-analysis simulation matrix has an unexpected ",
    "number of iterations."
  )
}

if (
  !all(
    parallel_factor_columns %in%
    colnames(parallel_results$values)
  )
) {
  stop(
    "The parallel-analysis simulation matrix does not contain ",
    "the expected F1 to F9 columns."
  )
}

parallel_observed_eigenvalues <- as.numeric(
  parallel_results$fa.values
)

if (
  length(parallel_observed_eigenvalues) !=
  length(phq9_items)
) {
  stop(
    "The parallel analysis returned an unexpected number of ",
    "observed factor eigenvalues."
  )
}

if (any(!is.finite(parallel_observed_eigenvalues))) {
  stop(
    "The observed parallel-analysis factor eigenvalues contain ",
    "non-finite values."
  )
}

parallel_simulated_factor_eigenvalues <- parallel_results$values[
  ,
  parallel_factor_columns,
  drop = FALSE
]

if (
  any(
    !is.finite(
      parallel_simulated_factor_eigenvalues
    )
  )
) {
  stop(
    "The simulated parallel-analysis factor eigenvalues contain ",
    "non-finite values."
  )
}

parallel_95th_thresholds <- apply(
  parallel_simulated_factor_eigenvalues,
  2,
  stats::quantile,
  probs = parallel_analysis_quantile,
  na.rm = TRUE,
  names = FALSE
)

if (
  length(parallel_95th_thresholds) !=
  length(phq9_items)
) {
  stop(
    "The reconstructed parallel-analysis thresholds have an ",
    "unexpected length."
  )
}

if (any(!is.finite(parallel_95th_thresholds))) {
  stop(
    "The reconstructed parallel-analysis thresholds contain ",
    "non-finite values."
  )
}

parallel_exceeds_threshold <- (
  parallel_observed_eigenvalues >
    parallel_95th_thresholds
)

parallel_custom_nfact <- if (
  all(parallel_exceeds_threshold)
) {
  length(parallel_exceeds_threshold)
} else {
  which(!parallel_exceeds_threshold)[1] - 1L
}

parallel_psych_nfact <- as.integer(
  parallel_results$nfact
)

if (
  length(parallel_psych_nfact) != 1L ||
  is.na(parallel_psych_nfact)
) {
  stop(
    "psych::fa.parallel() did not return one valid factor ",
    "recommendation."
  )
}

if (
  !identical(
    as.integer(parallel_custom_nfact),
    parallel_psych_nfact
  )
) {
  stop(
    "The reconstructed leading-consecutive recommendation does ",
    "not agree with psych::fa.parallel(). Custom recommendation: ",
    parallel_custom_nfact,
    "; psych recommendation: ",
    parallel_psych_nfact,
    "."
  )
}

parallel_comparison <- tibble(
  factor_number = seq_along(phq9_items),
  observed_eigenvalue = parallel_observed_eigenvalues,
  threshold_95th = parallel_95th_thresholds,
  exceeds_threshold = parallel_exceeds_threshold
)

parallel_analysis_diagnostics <- tibble(
  sample_n = nrow(development_items),
  item_n = ncol(development_items),
  seed = parallel_analysis_seed,
  iterations = parallel_analysis_iterations,
  quantile = parallel_analysis_quantile,
  warning_n = length(parallel_analysis_warnings),
  mc_cores = 1L,
  custom_recommended_factors = parallel_custom_nfact,
  psych_recommended_factors = parallel_psych_nfact,
  recommendations_agree = identical(
    as.integer(parallel_custom_nfact),
    parallel_psych_nfact
  )
)

# 14. Define the validated Stage 3 split object --------------------------------

stage3_split_object <- list(
  metadata = list(
    source_file = "data/processed/phq9_stage2a_data.rds",
    split_seed = split_seed,
    complete_sample_n = expected_complete_n,
    development_sample_n = expected_development_n,
    validation_sample_n = expected_validation_n,
    phq9_items = phq9_items,
    r_version = R.version.string,
    package_versions = c(
      dplyr = as.character(packageVersion("dplyr")),
      tidyr = as.character(packageVersion("tidyr")),
      readr = as.character(packageVersion("readr")),
      here = as.character(packageVersion("here")),
      tibble = as.character(packageVersion("tibble")),
      psych = as.character(packageVersion("psych"))
    )
  ),
  data = phq9_stage3_split,
  split_summary = split_summary,
  item_category_frequencies = item_category_frequencies,
  item_floor_summary = item_floor_summary,
  development_polychoric = list(
    correlation_matrix = development_poly_matrix,
    thresholds = development_poly_thresholds,
    eigenvalues = development_poly_eigenvalues,
    diagnostics = development_poly_diagnostics,
    warnings = unique(polychoric_warnings),
    estimation_settings = list(
      correct = 0.5,
      smooth = FALSE,
      global = TRUE,
      progress = FALSE
    )
  ),
  development_parallel_analysis = list(
    observed_factor_eigenvalues = parallel_observed_eigenvalues,
    simulated_factor_eigenvalues = (
      parallel_simulated_factor_eigenvalues
    ),
    threshold_95th = parallel_95th_thresholds,
    comparison = parallel_comparison,
    recommended_factors = parallel_custom_nfact,
    psych_recommended_factors = parallel_psych_nfact,
    diagnostics = parallel_analysis_diagnostics,
    warnings = unique(parallel_analysis_warnings),
    psych_call = parallel_results$Call,
    estimation_settings = list(
      seed = parallel_analysis_seed,
      iterations = parallel_analysis_iterations,
      fm = "minres",
      fa = "fa",
      nfactors = 1L,
      SMC = FALSE,
      sim = FALSE,
      quantile = parallel_analysis_quantile,
      correlation = "poly",
      continuity_correction = 0.5,
      plot = FALSE,
      mc_cores = 1L
    )
  )
)

# 15. Calculate PHQ-9 total-score distribution --------------------------------

phq9_total_distribution <- phq9_stage3_split |>
  mutate(phq9_total = rowSums(pick(all_of(phq9_items)))) |>
  count(stage3_sample, phq9_total, name = "frequency") |>
  complete(
    stage3_sample,
    phq9_total = 0:27,
    fill = list(frequency = 0L)
  ) |>
  group_by(stage3_sample) |>
  mutate(percentage = 100 * frequency / sum(frequency)) |>
  ungroup() |>
  arrange(stage3_sample, phq9_total)

stage3_split_object$phq9_total_distribution <- phq9_total_distribution

# 16. Save validated Stage 3 outputs -------------------------------------------

saveRDS(
  stage3_split_object,
  stage3_split_path
)

write_csv(
  split_summary,
  split_summary_path
)
