# ==============================================================================
# Stage 3: Dimensionality analysis
# Project: PHQ-9 psychometric evaluation using NHANES 2017–March 2020
#
# Current implementation scope:
#   - reconstruct the complete-PHQ-9 analytic sample with SEQN retained
#   - create and validate the locked 50:50 development and validation split
#   - estimate and validate the development-sample polychoric matrix
#   - run and validate the prespecified ordinal parallel analysis
#   - fit and validate the permitted development-sample EFA solutions
#   - record the development-stage model-freezing decision
#   - fit the prespecified validation-sample one-factor ordinal CFA
#   - extract and validate global-fit, loading, threshold, residual and
#     computational diagnostics
#   - record the validation-stage dimensionality conclusion
#   - save the complete recoverable Stage 3 analysis object
#
# Parallel analysis supported three factors, but neither multifactor EFA
# satisfied the locked criteria for freezing a secondary validation CFA.
#
# The one-factor ordinal CFA is the sole validation model.
# The final Stage 3 dimensionality conclusion is mixed evidence.
# No post hoc alternative models, correlated residuals, cross-loadings or
# item deletions are fitted in this script.
# ==============================================================================

# 1. Install any missing packages ----------------------------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "here",
  "tibble",
  "psych",
  "GPArotation",
  "lavaan"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0L) {
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
  stop(
    "The parallel-analysis simulation results are not a matrix."
  )
}

if (
  nrow(parallel_results$values) !=
  parallel_analysis_iterations
) {
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

if (parallel_custom_nfact != 3L) {
  stop(
    "The recorded Stage 3 model-freezing decision assumes that ",
    "parallel analysis recommends three factors. The current ",
    "recommendation is ",
    parallel_custom_nfact,
    ". Review the result before continuing."
  )
}

# 14. Fit the permitted development-sample EFAs -------------------------------

fit_development_efa <- function(
    number_of_factors,
    rotation_method
) {
  captured_warnings <- character()

  fitted_model <- withCallingHandlers(
    psych::fa(
      r = development_poly_matrix,
      nfactors = number_of_factors,
      n.obs = expected_development_n,
      rotate = rotation_method,
      residuals = TRUE,
      SMC = TRUE,
      fm = "minres",
      warnings = TRUE,
      smooth = FALSE
    ),
    warning = function(warning_condition) {
      captured_warnings <<- c(
        captured_warnings,
        conditionMessage(warning_condition)
      )

      invokeRestart("muffleWarning")
    }
  )

  list(
    model = fitted_model,
    warnings = unique(captured_warnings)
  )
}

efa_one_result <- fit_development_efa(
  number_of_factors = 1L,
  rotation_method = "none"
)

efa_two_result <- fit_development_efa(
  number_of_factors = 2L,
  rotation_method = "oblimin"
)

efa_three_result <- fit_development_efa(
  number_of_factors = 3L,
  rotation_method = "oblimin"
)

efa_models <- list(
  one_factor = efa_one_result$model,
  two_factor = efa_two_result$model,
  three_factor = efa_three_result$model
)

efa_warnings <- list(
  one_factor = efa_one_result$warnings,
  two_factor = efa_two_result$warnings,
  three_factor = efa_three_result$warnings
)

efa_expected_factors <- c(
  one_factor = 1L,
  two_factor = 2L,
  three_factor = 3L
)

efa_rotations <- c(
  one_factor = "none",
  two_factor = "oblimin",
  three_factor = "oblimin"
)

# 15. Validate the development-sample EFAs ------------------------------------

validate_development_efa <- function(
    efa_result,
    solution_name,
    expected_factors,
    captured_warnings
) {
  loading_matrix <- unclass(
    efa_result$loadings
  )

  communality_values <- as.numeric(
    efa_result$communality
  )

  uniqueness_values <- as.numeric(
    efa_result$uniquenesses
  )

  complexity_values <- as.numeric(
    efa_result$complexity
  )

  expected_loading_dimensions <- c(
    length(phq9_items),
    expected_factors
  )

  if (length(captured_warnings) > 0L) {
    stop(
      "The ",
      solution_name,
      " EFA produced one or more warnings: ",
      paste(captured_warnings, collapse = " | ")
    )
  }

  if (
    !identical(
      dim(loading_matrix),
      expected_loading_dimensions
    )
  ) {
    stop(
      "The ",
      solution_name,
      " loading matrix has unexpected dimensions."
    )
  }

  if (
    !identical(
      rownames(loading_matrix),
      phq9_items
    )
  ) {
    stop(
      "The ",
      solution_name,
      " loading matrix has an unexpected item order."
    )
  }

  if (!all(is.finite(loading_matrix))) {
    stop(
      "The ",
      solution_name,
      " loading matrix contains a non-finite estimate."
    )
  }

  if (
    length(communality_values) !=
    length(phq9_items) ||
    !all(is.finite(communality_values))
  ) {
    stop(
      "The ",
      solution_name,
      " communalities are invalid."
    )
  }

  if (
    length(uniqueness_values) !=
    length(phq9_items) ||
    !all(is.finite(uniqueness_values))
  ) {
    stop(
      "The ",
      solution_name,
      " uniquenesses are invalid."
    )
  }

  if (
    length(complexity_values) !=
    length(phq9_items) ||
    !all(is.finite(complexity_values))
  ) {
    stop(
      "The ",
      solution_name,
      " item-complexity values are invalid."
    )
  }

  if (
    any(communality_values < 0) ||
    any(communality_values >= 1)
  ) {
    stop(
      "The ",
      solution_name,
      " EFA contains an inadmissible communality."
    )
  }

  if (
    any(uniqueness_values <= 0) ||
    any(uniqueness_values > 1)
  ) {
    stop(
      "The ",
      solution_name,
      " EFA contains an inadmissible uniqueness."
    )
  }

  expected_correlation_dimensions <- c(
    length(phq9_items),
    length(phq9_items)
  )

  if (
    !identical(
      dim(efa_result$model),
      expected_correlation_dimensions
    ) ||
    !identical(
      dim(efa_result$residual),
      expected_correlation_dimensions
    )
  ) {
    stop(
      "The ",
      solution_name,
      " reproduced or residual matrix has unexpected dimensions."
    )
  }

  if (
    !identical(
      rownames(efa_result$model),
      phq9_items
    ) ||
    !identical(
      colnames(efa_result$model),
      phq9_items
    ) ||
    !identical(
      rownames(efa_result$residual),
      phq9_items
    ) ||
    !identical(
      colnames(efa_result$residual),
      phq9_items
    )
  ) {
    stop(
      "The ",
      solution_name,
      " reproduced or residual matrix has unexpected item names."
    )
  }

  if (
    !all(is.finite(efa_result$model)) ||
    !all(is.finite(efa_result$residual))
  ) {
    stop(
      "The ",
      solution_name,
      " reproduced or residual matrix contains a non-finite value."
    )
  }

  if (
    max(
      abs(
        efa_result$model -
        t(efa_result$model)
      )
    ) > 1e-12
  ) {
    stop(
      "The ",
      solution_name,
      " reproduced correlation matrix is not symmetric."
    )
  }

  if (
    max(
      abs(
        efa_result$residual -
        t(efa_result$residual)
      )
    ) > 1e-12
  ) {
    stop(
      "The ",
      solution_name,
      " residual correlation matrix is not symmetric."
    )
  }

  reconstruction_difference <- max(
    abs(
      development_poly_matrix -
        efa_result$model -
        efa_result$residual
    )
  )

  if (reconstruction_difference > 1e-10) {
    stop(
      "The ",
      solution_name,
      " reproduced and residual matrices do not reconstruct ",
      "the development polychoric matrix within tolerance."
    )
  }

  required_scalar_components <- c(
    statistic = efa_result$STATISTIC,
    p_value = efa_result$PVAL,
    degrees_of_freedom = efa_result$dof,
    rmsr = efa_result$rms,
    tli = efa_result$TLI
  )

  if (
    length(required_scalar_components) != 5L ||
    !all(is.finite(required_scalar_components))
  ) {
    stop(
      "The ",
      solution_name,
      " EFA has an invalid required fit statistic."
    )
  }

  required_rmsea_names <- c(
    "RMSEA",
    "lower",
    "upper",
    "confidence"
  )

  if (
    !all(
      required_rmsea_names %in%
      names(efa_result$RMSEA)
    ) ||
    !all(
      is.finite(
        efa_result$RMSEA[
          required_rmsea_names
        ]
      )
    )
  ) {
    stop(
      "The ",
      solution_name,
      " EFA has invalid RMSEA information."
    )
  }

  objective_value <- unname(
    efa_result$criteria["objective"]
  )

  if (
    length(objective_value) != 1L ||
    !is.finite(objective_value)
  ) {
    stop(
      "The ",
      solution_name,
      " EFA has an invalid objective value."
    )
  }

  if (expected_factors > 1L) {
    factor_correlations <- efa_result$Phi

    expected_factor_dimensions <- c(
      expected_factors,
      expected_factors
    )

    if (
      !identical(
        dim(factor_correlations),
        expected_factor_dimensions
      )
    ) {
      stop(
        "The ",
        solution_name,
        " factor-correlation matrix has unexpected dimensions."
      )
    }

    if (!all(is.finite(factor_correlations))) {
      stop(
        "The ",
        solution_name,
        " factor-correlation matrix contains a non-finite value."
      )
    }

    if (
      max(
        abs(
          factor_correlations -
          t(factor_correlations)
        )
      ) > 1e-12
    ) {
      stop(
        "The ",
        solution_name,
        " factor-correlation matrix is not symmetric."
      )
    }

    if (
      max(
        abs(
          diag(factor_correlations) - 1
        )
      ) > 1e-12
    ) {
      stop(
        "The ",
        solution_name,
        " factor-correlation matrix does not have a unit diagonal."
      )
    }

    factor_correlation_tolerance <- 1e-12

    if (
      any(
        factor_correlations < (-1 - factor_correlation_tolerance) |
        factor_correlations > (1 + factor_correlation_tolerance)
      )
    ) {
      stop(
        "The ",
        solution_name,
        " factor-correlation matrix contains a correlation outside ",
        "the valid range beyond numerical tolerance."
      )
    }

    factor_correlation_eigenvalues <- eigen(
      factor_correlations,
      symmetric = TRUE,
      only.values = TRUE
    )$values

    if (
      min(factor_correlation_eigenvalues) <= 0
    ) {
      stop(
        "The ",
        solution_name,
        " factor-correlation matrix is not positive definite."
      )
    }
  }

  invisible(TRUE)
}

for (solution_name in names(efa_models)) {
  validate_development_efa(
    efa_result = efa_models[[solution_name]],
    solution_name = solution_name,
    expected_factors = (
      efa_expected_factors[[solution_name]]
    ),
    captured_warnings = efa_warnings[[solution_name]]
  )
}

# 16. Extract EFA loadings and item diagnostics --------------------------------

efa_loading_table <- dplyr::bind_rows(
  lapply(
    names(efa_models),
    function(solution_name) {
      efa_result <- efa_models[[solution_name]]
      pattern_matrix <- unclass(
        efa_result$loadings
      )

      dplyr::bind_rows(
        lapply(
          seq_len(ncol(pattern_matrix)),
          function(factor_index) {
            tibble::tibble(
              solution = solution_name,
              factors = (
                efa_expected_factors[[solution_name]]
              ),
              rotation = efa_rotations[[solution_name]],
              item = rownames(pattern_matrix),
              factor = colnames(pattern_matrix)[
                factor_index
              ],
              loading = as.numeric(
                pattern_matrix[, factor_index]
              ),
              communality = as.numeric(
                efa_result$communality
              ),
              uniqueness = as.numeric(
                efa_result$uniquenesses
              ),
              complexity = as.numeric(
                efa_result$complexity
              )
            )
          }
        )
      )
    }
  )
)

efa_loading_separation <- dplyr::bind_rows(
  lapply(
    names(efa_models),
    function(solution_name) {
      efa_result <- efa_models[[solution_name]]
      pattern_matrix <- unclass(
        efa_result$loadings
      )

      dplyr::bind_rows(
        lapply(
          seq_len(nrow(pattern_matrix)),
          function(item_index) {
            item_loadings <- pattern_matrix[
              item_index,
              ,
              drop = TRUE
            ]

            absolute_loadings <- abs(
              item_loadings
            )

            loading_order <- order(
              absolute_loadings,
              decreasing = TRUE
            )

            primary_index <- loading_order[1]

            secondary_loading <- NA_real_
            absolute_loading_gap <- NA_real_

            if (length(loading_order) > 1L) {
              secondary_index <- loading_order[2]

              secondary_loading <- item_loadings[
                secondary_index
              ]

              absolute_loading_gap <- (
                absolute_loadings[primary_index] -
                  absolute_loadings[secondary_index]
              )
            }

            tibble::tibble(
              solution = solution_name,
              factors = (
                efa_expected_factors[[solution_name]]
              ),
              item = rownames(pattern_matrix)[
                item_index
              ],
              primary_factor = colnames(pattern_matrix)[
                primary_index
              ],
              primary_loading = item_loadings[
                primary_index
              ],
              secondary_loading = secondary_loading,
              absolute_loading_gap = absolute_loading_gap,
              primary_loading_below_0_40 = (
                absolute_loadings[primary_index] < 0.40
              ),
              cross_loading_ge_0_30 = (
                !is.na(secondary_loading) &&
                  abs(secondary_loading) >= 0.30
              ),
              loading_gap_below_0_20 = (
                !is.na(absolute_loading_gap) &&
                  absolute_loading_gap < 0.20
              ),
              complexity = as.numeric(
                efa_result$complexity[item_index]
              )
            )
          }
        )
      )
    }
  )
)

# 17. Extract EFA correlations and diagnostics ---------------------------------

extract_lower_triangle <- function(
    correlation_matrix,
    solution_name,
    value_name
) {
  matrix_indices <- which(
    lower.tri(correlation_matrix),
    arr.ind = TRUE
  )

  output <- tibble::tibble(
    solution = solution_name,
    variable_1 = rownames(correlation_matrix)[
      matrix_indices[, 1]
    ],
    variable_2 = colnames(correlation_matrix)[
      matrix_indices[, 2]
    ],
    value = as.numeric(
      correlation_matrix[matrix_indices]
    )
  )

  names(output)[names(output) == "value"] <- value_name

  output
}

efa_factor_correlations <- dplyr::bind_rows(
  lapply(
    names(efa_models),
    function(solution_name) {
      efa_result <- efa_models[[solution_name]]

      if (
        efa_expected_factors[[solution_name]] == 1L
      ) {
        return(
          tibble::tibble(
            solution = character(),
            factor_1 = character(),
            factor_2 = character(),
            correlation = numeric()
          )
        )
      }

      matrix_indices <- which(
        lower.tri(efa_result$Phi),
        arr.ind = TRUE
      )

      tibble::tibble(
        solution = solution_name,
        factor_1 = rownames(efa_result$Phi)[
          matrix_indices[, 1]
        ],
        factor_2 = colnames(efa_result$Phi)[
          matrix_indices[, 2]
        ],
        correlation = as.numeric(
          efa_result$Phi[matrix_indices]
        )
      )
    }
  )
)

efa_reproduced_correlations <- dplyr::bind_rows(
  lapply(
    names(efa_models),
    function(solution_name) {
      extract_lower_triangle(
        correlation_matrix = (
          efa_models[[solution_name]]$model
        ),
        solution_name = solution_name,
        value_name = "reproduced_correlation"
      )
    }
  )
)

efa_residual_correlations <- dplyr::bind_rows(
  lapply(
    names(efa_models),
    function(solution_name) {
      extract_lower_triangle(
        correlation_matrix = (
          efa_models[[solution_name]]$residual
        ),
        solution_name = solution_name,
        value_name = "residual_correlation"
      )
    }
  )
)

efa_model_diagnostics <- dplyr::bind_rows(
  lapply(
    names(efa_models),
    function(solution_name) {
      efa_result <- efa_models[[solution_name]]

      residual_off_diagonal <- efa_result$residual[
        lower.tri(efa_result$residual)
      ]

      maximum_absolute_factor_correlation <- NA_real_

      if (
        efa_expected_factors[[solution_name]] > 1L
      ) {
        maximum_absolute_factor_correlation <- max(
          abs(
            efa_result$Phi[
              lower.tri(efa_result$Phi)
            ]
          )
        )
      }

      tibble::tibble(
        solution = solution_name,
        factors = efa_expected_factors[[solution_name]],
        rotation = efa_rotations[[solution_name]],
        sample_n = expected_development_n,
        warning_n = length(
          efa_warnings[[solution_name]]
        ),
        minimum_communality = min(
          efa_result$communality
        ),
        maximum_communality = max(
          efa_result$communality
        ),
        minimum_uniqueness = min(
          efa_result$uniquenesses
        ),
        maximum_uniqueness = max(
          efa_result$uniquenesses
        ),
        maximum_complexity = max(
          efa_result$complexity
        ),
        largest_absolute_residual = max(
          abs(residual_off_diagonal)
        ),
        residual_correlations_ge_0_10 = sum(
          abs(residual_off_diagonal) >= 0.10
        ),
        chi_square = unname(
          efa_result$STATISTIC
        ),
        degrees_of_freedom = unname(
          efa_result$dof
        ),
        p_value = unname(
          efa_result$PVAL
        ),
        rmsr = unname(
          efa_result$rms
        ),
        rmsea = unname(
          efa_result$RMSEA["RMSEA"]
        ),
        rmsea_lower = unname(
          efa_result$RMSEA["lower"]
        ),
        rmsea_upper = unname(
          efa_result$RMSEA["upper"]
        ),
        rmsea_confidence = unname(
          efa_result$RMSEA["confidence"]
        ),
        tli = unname(
          efa_result$TLI
        ),
        maximum_absolute_factor_correlation = (
          maximum_absolute_factor_correlation
        ),
        objective = unname(
          efa_result$criteria["objective"]
        ),
        computational_checks_passed = TRUE
      )
    }
  )
)

# 18. Record the development-stage model-freezing decision ---------------------

two_factor_decision_evidence <- efa_loading_separation |>
  filter(
    solution == "two_factor",
    item %in% c(
      "DPQ010",
      "DPQ070",
      "DPQ080"
    )
  )

required_two_factor_flagged_items <- c(
  "DPQ010",
  "DPQ070",
  "DPQ080"
)

observed_two_factor_flagged_items <- (
  two_factor_decision_evidence |>
    filter(
      cross_loading_ge_0_30,
      loading_gap_below_0_20
    ) |>
    pull(item)
)

if (
  !all(
    required_two_factor_flagged_items %in%
    observed_two_factor_flagged_items
  )
) {
  stop(
    "The recorded two-factor model-freezing decision is no ",
    "longer supported by the expected cross-loading and loading-",
    "separation evidence."
  )
}

three_factor_primary_assignments <- (
  efa_loading_separation |>
    filter(solution == "three_factor")
)

three_factor_names <- colnames(
  unclass(
    efa_models$three_factor$loadings
  )
)

three_factor_primary_counts <- tibble::tibble(
  factor = three_factor_names,
  primary_indicator_n = as.integer(
    table(
      factor(
        three_factor_primary_assignments$primary_factor,
        levels = three_factor_names
      )
    )
  )
)

if (
  min(
    three_factor_primary_counts$primary_indicator_n
  ) >= 3L
) {
  stop(
    "The recorded three-factor model-freezing decision is no ",
    "longer supported because every factor now has at least ",
    "three primary indicators."
  )
}

development_model_decision <- tibble::tibble(
  candidate_solution = c(
    "two_factor",
    "three_factor"
  ),
  factor_count = c(2L, 3L),
  supported_by_parallel_analysis = c(
    parallel_custom_nfact >= 2L,
    parallel_custom_nfact >= 3L
  ),
  retained_for_validation = c(
    FALSE,
    FALSE
  ),
  decision_date = c(
    "2026-07-24",
    "2026-07-24"
  ),
  decision_reason = c(
    paste(
      "Not frozen because DPQ010, DPQ070 and DPQ080 showed",
      "cross-loadings of at least .30 and primary-versus-",
      "secondary loading gaps below .20. The factors also",
      "correlated strongly, and the pattern did not support",
      "a sufficiently stable or substantively distinct",
      "simple-structure CFA allocation."
    ),
    paste(
      "Not frozen because one factor was primarily defined by",
      "only two clear indicators, failing the requirement that",
      "every frozen factor contain at least three substantively",
      "coherent indicators. The factors were also highly",
      "correlated."
    )
  )
)

validation_model_plan <- tibble::tibble(
  primary_validation_model = "one_factor",
  secondary_model_frozen = FALSE,
  secondary_validation_model = NA_character_,
  decision_date = "2026-07-24",
  decision_recorded_before_validation_access = TRUE,
  validation_sample_factor_results_accessed = FALSE
)

# 19. Calculate the PHQ-9 total-score distribution -----------------------------

phq9_total_distribution <- phq9_stage3_split |>
  mutate(
    phq9_total = rowSums(
      pick(all_of(phq9_items))
    )
  ) |>
  count(
    stage3_sample,
    phq9_total,
    name = "frequency"
  ) |>
  complete(
    stage3_sample,
    phq9_total = 0:27,
    fill = list(frequency = 0L)
  ) |>
  group_by(stage3_sample) |>
  mutate(
    percentage = 100 * frequency / sum(frequency)
  ) |>
  ungroup() |>
  arrange(
    stage3_sample,
    phq9_total
  )

# 20. Prepare and validate the validation-sample data --------------------------

validation_data <- phq9_stage3_split |>
  filter(
    stage3_sample == "validation"
  ) |>
  select(
    all_of(phq9_items)
  )

validation_category_counts <- lapply(
  validation_data,
  table,
  useNA = "always"
)

validation_data_checks <- c(
  expected_participant_n = (
    nrow(validation_data) ==
      expected_validation_n
  ),
  expected_item_n = (
    ncol(validation_data) ==
      length(phq9_items)
  ),
  expected_item_order = identical(
    names(validation_data),
    phq9_items
  ),
  no_missing_responses = (
    sum(is.na(validation_data)) == 0L
  ),
  all_values_finite = all(
    is.finite(
      as.matrix(validation_data)
    )
  ),
  valid_response_range = all(
    as.matrix(validation_data) >= 0 &
      as.matrix(validation_data) <= 3
  ),
  integer_response_coding = all(
    as.matrix(validation_data) ==
      floor(as.matrix(validation_data))
  ),
  all_categories_present_per_item = all(
    vapply(
      validation_data,
      function(item_values) {
        setequal(
          sort(
            unique(
              as.numeric(item_values)
            )
          ),
          0:3
        )
      },
      logical(1)
    )
  )
)

if (!all(validation_data_checks)) {
  failed_validation_checks <- names(
    validation_data_checks
  )[
    !validation_data_checks
  ]

  stop(
    "The validation-sample data failed the following checks: ",
    paste(
      failed_validation_checks,
      collapse = ", "
    )
  )
}

if (
  !identical(
    validation_model_plan$primary_validation_model,
    "one_factor"
  ) ||
  !identical(
    validation_model_plan$secondary_model_frozen,
    FALSE
  )
) {
  stop(
    "The stored development-stage model plan does not permit ",
    "the prespecified one-factor-only validation analysis."
  )
}

# 21. Fit the prespecified one-factor ordinal CFA ------------------------------

one_factor_model <- '
  Depression =~ DPQ010 + DPQ020 + DPQ030 +
                DPQ040 + DPQ050 + DPQ060 +
                DPQ070 + DPQ080 + DPQ090
'

cfa_warnings <- character()

cfa_fit <- withCallingHandlers(
  lavaan::cfa(
    model = one_factor_model,
    data = validation_data,
    ordered = phq9_items,
    estimator = "WLSMV",
    std.lv = TRUE
  ),
  warning = function(warning_condition) {
    cfa_warnings <<- c(
      cfa_warnings,
      conditionMessage(warning_condition)
    )

    invokeRestart("muffleWarning")
  }
)

cfa_warnings <- unique(cfa_warnings)

if (length(cfa_warnings) > 0L) {
  stop(
    "The validation CFA produced one or more warnings: ",
    paste(
      cfa_warnings,
      collapse = " | "
    )
  )
}

cfa_converged <- lavaan::lavInspect(
  cfa_fit,
  "converged"
)

cfa_post_check <- lavaan::lavInspect(
  cfa_fit,
  "post.check"
)

cfa_observation_n <- as.integer(
  lavaan::lavInspect(
    cfa_fit,
    "nobs"
  )
)

cfa_free_parameter_n <- as.integer(
  lavaan::lavInspect(
    cfa_fit,
    "npar"
  )
)

if (!identical(cfa_converged, TRUE)) {
  stop(
    "The validation CFA did not converge."
  )
}

if (!identical(cfa_post_check, TRUE)) {
  stop(
    "The validation CFA failed lavaan's post-estimation ",
    "admissibility check."
  )
}

if (
  !identical(
    cfa_observation_n,
    expected_validation_n
  )
) {
  stop(
    "The validation CFA used an unexpected number of ",
    "observations."
  )
}

if (!identical(cfa_free_parameter_n, 36L)) {
  stop(
    "The validation CFA returned an unexpected number of ",
    "free parameters."
  )
}

# 22. Extract resolved options and model-test information ----------------------

cfa_resolved_options <- lavaan::lavInspect(
  cfa_fit,
  "options"
)

required_resolved_options <- c(
  "estimator",
  "estimator.orig",
  "parameterization",
  "std.lv",
  "se",
  "test",
  "information",
  "information.meat"
)

missing_resolved_options <- setdiff(
  required_resolved_options,
  names(cfa_resolved_options)
)

if (length(missing_resolved_options) > 0L) {
  stop(
    "The fitted CFA does not expose the following required ",
    "resolved options: ",
    paste(
      missing_resolved_options,
      collapse = ", "
    )
  )
}

resolved_option_checks <- c(
  requested_estimator_wlsmv = identical(
    cfa_resolved_options$estimator.orig,
    "WLSMV"
  ),
  internal_estimator_dwls = identical(
    cfa_resolved_options$estimator,
    "DWLS"
  ),
  delta_parameterization = identical(
    cfa_resolved_options$parameterization,
    "delta"
  ),
  standardised_latent_variable = identical(
    cfa_resolved_options$std.lv,
    TRUE
  ),
  robust_standard_errors = identical(
    cfa_resolved_options$se,
    "robust.sem"
  ),
  scaled_shifted_test_available = (
    "scaled.shifted" %in%
      cfa_resolved_options$test
  )
)

if (!all(resolved_option_checks)) {
  failed_option_checks <- names(
    resolved_option_checks
  )[
    !resolved_option_checks
  ]

  stop(
    "The fitted CFA failed the following resolved-option checks: ",
    paste(
      failed_option_checks,
      collapse = ", "
    )
  )
}

cfa_test_information <- lavaan::lavInspect(
  cfa_fit,
  "test"
)

required_cfa_tests <- c(
  "standard",
  "scaled.shifted"
)

if (
  !all(
    required_cfa_tests %in%
    names(cfa_test_information)
  )
) {
  stop(
    "The fitted CFA does not contain both the standard and ",
    "scaled-and-shifted model tests."
  )
}

cfa_standard_test <- cfa_test_information$standard
cfa_scaled_shifted_test <- (
  cfa_test_information$scaled.shifted
)

required_scaled_shifted_components <- c(
  "stat",
  "df",
  "pvalue",
  "scaling.factor",
  "shift.parameter",
  "scaled.test.stat",
  "scaled.test",
  "label"
)

missing_scaled_shifted_components <- setdiff(
  required_scaled_shifted_components,
  names(cfa_scaled_shifted_test)
)

if (
  length(missing_scaled_shifted_components) > 0L
) {
  stop(
    "The scaled-and-shifted test is missing the following ",
    "required components: ",
    paste(
      missing_scaled_shifted_components,
      collapse = ", "
    )
  )
}

required_model_test_values <- c(
  standard_statistic = cfa_standard_test$stat,
  standard_degrees_of_freedom = cfa_standard_test$df,
  scaled_shifted_statistic = (
    cfa_scaled_shifted_test$stat
  ),
  scaled_shifted_degrees_of_freedom = (
    cfa_scaled_shifted_test$df
  ),
  scaled_shifted_p_value = (
    cfa_scaled_shifted_test$pvalue
  ),
  scaling_factor = (
    cfa_scaled_shifted_test$scaling.factor
  ),
  shift_parameter = (
    cfa_scaled_shifted_test$shift.parameter
  ),
  scaled_test_statistic = (
    cfa_scaled_shifted_test$scaled.test.stat
  )
)

if (
  !all(
    is.finite(
      required_model_test_values
    )
  )
) {
  stop(
    "The validation CFA contains a non-finite required ",
    "model-test value."
  )
}

if (
  abs(
    cfa_standard_test$stat -
    cfa_scaled_shifted_test$scaled.test.stat
  ) > 1e-10
) {
  stop(
    "The standard test statistic does not match the stored ",
    "underlying statistic for the scaled-and-shifted test."
  )
}

if (
  !identical(
    as.integer(cfa_standard_test$df),
    as.integer(cfa_scaled_shifted_test$df)
  )
) {
  stop(
    "The standard and scaled-and-shifted tests have different ",
    "degrees of freedom."
  )
}

cfa_model_test_table <- tibble::tibble(
  test = c(
    "standard",
    "scaled_and_shifted"
  ),
  statistic = c(
    cfa_standard_test$stat,
    cfa_scaled_shifted_test$stat
  ),
  degrees_of_freedom = c(
    cfa_standard_test$df,
    cfa_scaled_shifted_test$df
  ),
  p_value = c(
    cfa_standard_test$pvalue,
    cfa_scaled_shifted_test$pvalue
  ),
  scaling_factor = c(
    NA_real_,
    cfa_scaled_shifted_test$scaling.factor
  ),
  shift_parameter = c(
    NA_real_,
    cfa_scaled_shifted_test$shift.parameter
  ),
  underlying_standard_statistic = c(
    NA_real_,
    cfa_scaled_shifted_test$scaled.test.stat
  ),
  correction_label = c(
    NA_character_,
    cfa_scaled_shifted_test$label
  )
)

# 23. Extract and validate the prespecified fit measures -----------------------

cfa_fit_measures_all <- lavaan::fitMeasures(
  cfa_fit
)

requested_fit_measure_names <- c(
  "chisq.scaled",
  "df.scaled",
  "pvalue.scaled",
  "cfi.robust",
  "tli.robust",
  "rmsea.robust",
  "rmsea.ci.lower.robust",
  "rmsea.ci.upper.robust",
  "srmr"
)

missing_requested_measures <- setdiff(
  requested_fit_measure_names,
  names(cfa_fit_measures_all)
)

if (length(missing_requested_measures) > 0L) {
  stop(
    "The following prespecified CFA fit measures are unavailable: ",
    paste(
      missing_requested_measures,
      collapse = ", "
    )
  )
}

requested_fit_measure_values <- (
  cfa_fit_measures_all[
    requested_fit_measure_names
  ]
)

if (
  any(
    is.na(
      requested_fit_measure_values
    )
  ) ||
  any(
    !is.finite(
      requested_fit_measure_values
    )
  )
) {
  stop(
    "One or more prespecified CFA fit measures are missing or ",
    "non-finite."
  )
}

cfa_fit_measure_table <- tibble::tibble(
  measure = requested_fit_measure_names,
  value = unname(
    requested_fit_measure_values
  )
)

if (
  abs(
    cfa_fit_measures_all["chisq.scaled"] -
    cfa_scaled_shifted_test$stat
  ) > 1e-10 ||
  abs(
    cfa_fit_measures_all["df.scaled"] -
    cfa_scaled_shifted_test$df
  ) > 1e-10 ||
  abs(
    cfa_fit_measures_all["pvalue.scaled"] -
    cfa_scaled_shifted_test$pvalue
  ) > 1e-10
) {
  stop(
    "The scaled fit-measure entries do not agree with the ",
    "scaled-and-shifted model-test object."
  )
}

# 24. Extract and validate loadings and thresholds -----------------------------

cfa_standardized_all <- tibble::as_tibble(
  lavaan::standardizedSolution(
    cfa_fit,
    type = "std.all",
    ci = TRUE,
    level = 0.95
  )
)

cfa_standardized_loadings <- cfa_standardized_all |>
  filter(
    op == "=~"
  ) |>
  transmute(
    factor = lhs,
    item = rhs,
    standardized_loading = est.std,
    standard_error = se,
    ci_lower = ci.lower,
    ci_upper = ci.upper,
    p_value = pvalue
  ) |>
  tibble::as_tibble()

cfa_loading_checks <- c(
  expected_loading_n = (
    nrow(cfa_standardized_loadings) ==
      length(phq9_items)
  ),
  expected_item_order = identical(
    cfa_standardized_loadings$item,
    phq9_items
  ),
  all_estimates_finite = all(
    is.finite(
      cfa_standardized_loadings$
        standardized_loading
    )
  ),
  all_standard_errors_finite = all(
    is.finite(
      cfa_standardized_loadings$
        standard_error
    )
  ),
  all_intervals_finite = all(
    is.finite(
      cfa_standardized_loadings$ci_lower
    ) &
      is.finite(
        cfa_standardized_loadings$ci_upper
      )
  ),
  intervals_ordered = all(
    cfa_standardized_loadings$ci_lower <=
      cfa_standardized_loadings$
      standardized_loading &
      cfa_standardized_loadings$
      standardized_loading <=
      cfa_standardized_loadings$ci_upper
  ),
  all_loadings_positive = all(
    cfa_standardized_loadings$
      standardized_loading > 0
  )
)

if (!all(cfa_loading_checks)) {
  failed_loading_checks <- names(
    cfa_loading_checks
  )[
    !cfa_loading_checks
  ]

  stop(
    "The CFA loadings failed the following checks: ",
    paste(
      failed_loading_checks,
      collapse = ", "
    )
  )
}

cfa_parameters <- tibble::as_tibble(
  lavaan::parameterEstimates(
    cfa_fit,
    ci = TRUE,
    level = 0.95
  )
)

cfa_thresholds <- cfa_parameters |>
  filter(
    op == "|"
  ) |>
  transmute(
    item = lhs,
    threshold = rhs,
    estimate = est,
    standard_error = se,
    ci_lower = ci.lower,
    ci_upper = ci.upper,
    p_value = pvalue
  ) |>
  tibble::as_tibble()

threshold_order_checks <- cfa_thresholds |>
  group_by(item) |>
  summarise(
    threshold_n = n(),
    estimates_strictly_increasing = all(
      diff(estimate) > 0
    ),
    .groups = "drop"
  )

cfa_threshold_checks <- c(
  expected_threshold_n = (
    nrow(cfa_thresholds) ==
      3L * length(phq9_items)
  ),
  expected_item_n = (
    n_distinct(cfa_thresholds$item) ==
      length(phq9_items)
  ),
  expected_item_order = identical(
    unique(cfa_thresholds$item),
    phq9_items
  ),
  three_thresholds_per_item = all(
    threshold_order_checks$threshold_n == 3L
  ),
  all_estimates_finite = all(
    is.finite(cfa_thresholds$estimate)
  ),
  all_standard_errors_finite = all(
    is.finite(cfa_thresholds$standard_error)
  ),
  all_intervals_finite = all(
    is.finite(cfa_thresholds$ci_lower) &
      is.finite(cfa_thresholds$ci_upper)
  ),
  intervals_ordered = all(
    cfa_thresholds$ci_lower <=
      cfa_thresholds$estimate &
      cfa_thresholds$estimate <=
      cfa_thresholds$ci_upper
  ),
  thresholds_strictly_increasing = all(
    threshold_order_checks$
      estimates_strictly_increasing
  )
)

if (!all(cfa_threshold_checks)) {
  failed_threshold_checks <- names(
    cfa_threshold_checks
  )[
    !cfa_threshold_checks
  ]

  stop(
    "The CFA thresholds failed the following checks: ",
    paste(
      failed_threshold_checks,
      collapse = ", "
    )
  )
}

# 25. Extract and validate residual diagnostics --------------------------------

cfa_residual_warnings <- character()

cfa_residuals <- withCallingHandlers(
  lavaan::lavResiduals(
    cfa_fit,
    type = "cor.bentler",
    zstat = TRUE,
    summary = TRUE
  ),
  warning = function(warning_condition) {
    cfa_residual_warnings <<- c(
      cfa_residual_warnings,
      conditionMessage(warning_condition)
    )

    invokeRestart("muffleWarning")
  }
)

cfa_residual_warnings <- unique(
  cfa_residual_warnings
)

if (length(cfa_residual_warnings) > 0L) {
  stop(
    "The CFA residual diagnostics produced one or more warnings: ",
    paste(
      cfa_residual_warnings,
      collapse = " | "
    )
  )
}

required_residual_components <- c(
  "cov",
  "cov.z",
  "summary"
)

missing_residual_components <- setdiff(
  required_residual_components,
  names(cfa_residuals)
)

if (length(missing_residual_components) > 0L) {
  stop(
    "The residual object is missing the following required ",
    "components: ",
    paste(
      missing_residual_components,
      collapse = ", "
    )
  )
}

cfa_residual_matrix <- as.matrix(
  cfa_residuals$cov
)

cfa_residual_z_matrix <- as.matrix(
  cfa_residuals$cov.z
)

cfa_observed_correlations <- as.matrix(
  lavaan::lavInspect(
    cfa_fit,
    "sampstat"
  )$cov
)

cfa_implied_correlations <- as.matrix(
  lavaan::lavInspect(
    cfa_fit,
    "cov.ov"
  )
)

residual_indices <- which(
  lower.tri(cfa_residual_matrix),
  arr.ind = TRUE
)

cfa_residual_correlations <- tibble::tibble(
  item_1 = rownames(cfa_residual_matrix)[
    residual_indices[, 1]
  ],
  item_2 = colnames(cfa_residual_matrix)[
    residual_indices[, 2]
  ],
  residual_correlation = as.numeric(
    cfa_residual_matrix[
      residual_indices
    ]
  ),
  standardized_residual_z = as.numeric(
    cfa_residual_z_matrix[
      residual_indices
    ]
  )
) |>
  mutate(
    absolute_residual_correlation = abs(
      residual_correlation
    ),
    flagged_ge_0_10 = (
      absolute_residual_correlation >= 0.10
    ),
    residual_convention = (
      "observed minus model-implied correlation"
    )
  ) |>
  arrange(
    desc(
      absolute_residual_correlation
    )
  )

cfa_residual_checks <- c(
  expected_matrix_dimensions = identical(
    dim(cfa_residual_matrix),
    c(
      length(phq9_items),
      length(phq9_items)
    )
  ),
  expected_z_matrix_dimensions = identical(
    dim(cfa_residual_z_matrix),
    c(
      length(phq9_items),
      length(phq9_items)
    )
  ),
  expected_item_order = (
    identical(
      rownames(cfa_residual_matrix),
      phq9_items
    ) &&
      identical(
        colnames(cfa_residual_matrix),
        phq9_items
      )
  ),
  residual_matrix_symmetric = (
    max(
      abs(
        cfa_residual_matrix -
          t(cfa_residual_matrix)
      )
    ) < 1e-12
  ),
  z_matrix_symmetric = (
    max(
      abs(
        cfa_residual_z_matrix -
          t(cfa_residual_z_matrix)
      )
    ) < 1e-12
  ),
  residual_diagonal_zero = (
    max(
      abs(
        diag(cfa_residual_matrix)
      )
    ) < 1e-12
  ),
  expected_pair_n = (
    nrow(cfa_residual_correlations) ==
      choose(length(phq9_items), 2L)
  ),
  all_residuals_finite = all(
    is.finite(
      cfa_residual_correlations$
        residual_correlation
    )
  ),
  all_z_values_finite = all(
    is.finite(
      cfa_residual_correlations$
        standardized_residual_z
    )
  ),
  matches_observed_minus_implied = (
    max(
      abs(
        cfa_residual_matrix -
          (
            cfa_observed_correlations -
              cfa_implied_correlations
          )
      )
    ) < 1e-10
  )
)

if (!all(cfa_residual_checks)) {
  failed_residual_checks <- names(
    cfa_residual_checks
  )[
    !cfa_residual_checks
  ]

  stop(
    "The CFA residual diagnostics failed the following checks: ",
    paste(
      failed_residual_checks,
      collapse = ", "
    )
  )
}

cfa_residual_summary <- tibble::as_tibble(
  cfa_residuals$summary,
  rownames = "summary_measure"
)

# 26. Complete variance and improper-solution checks ---------------------------

cfa_standardized_variances <- cfa_standardized_all |>
  filter(
    op == "~~",
    lhs == rhs
  ) |>
  transmute(
    variable = lhs,
    standardized_variance = est.std,
    standard_error = se,
    ci_lower = ci.lower,
    ci_upper = ci.upper
  ) |>
  tibble::as_tibble()

cfa_standardized_residual_variances <- (
  cfa_standardized_variances |>
    filter(
      variable %in% phq9_items
    )
)

cfa_standardized_latent_variance <- (
  cfa_standardized_variances |>
    filter(
      variable == "Depression"
    )
)

cfa_theta <- as.matrix(
  lavaan::lavInspect(
    cfa_fit,
    "theta"
  )
)

cfa_latent_covariance <- as.matrix(
  lavaan::lavInspect(
    cfa_fit,
    "cov.lv"
  )
)

variance_tolerance <- 1e-12

cfa_variance_checks <- c(
  expected_observed_residual_variance_n = (
    nrow(
      cfa_standardized_residual_variances
    ) == length(phq9_items)
  ),
  all_standardized_variances_finite = all(
    is.finite(
      cfa_standardized_residual_variances$
        standardized_variance
    )
  ),
  standardized_residual_variances_positive = all(
    cfa_standardized_residual_variances$
      standardized_variance > 0
  ),
  standardized_residual_variances_within_bounds = all(
    cfa_standardized_residual_variances$
      standardized_variance <=
      (1 + variance_tolerance)
  ),
  theta_dimensions_correct = identical(
    dim(cfa_theta),
    c(
      length(phq9_items),
      length(phq9_items)
    )
  ),
  theta_finite = all(
    is.finite(cfa_theta)
  ),
  theta_diagonal_positive = all(
    diag(cfa_theta) > 0
  ),
  latent_covariance_dimensions_correct = identical(
    dim(cfa_latent_covariance),
    c(1L, 1L)
  ),
  latent_covariance_finite = all(
    is.finite(cfa_latent_covariance)
  ),
  latent_factor_variance_fixed_to_one = (
    abs(
      cfa_latent_covariance[1, 1] - 1
    ) <= variance_tolerance
  ),
  standardized_latent_variance_row_present = (
    nrow(cfa_standardized_latent_variance) == 1L
  ),
  standardized_latent_variance_equals_one = (
    nrow(cfa_standardized_latent_variance) == 1L &&
      abs(
        cfa_standardized_latent_variance$
          standardized_variance[1] - 1
      ) <= variance_tolerance
  )
)

if (!all(cfa_variance_checks)) {
  failed_variance_checks <- names(
    cfa_variance_checks
  )[
    !cfa_variance_checks
  ]

  stop(
    "The CFA variance diagnostics failed the following checks: ",
    paste(
      failed_variance_checks,
      collapse = ", "
    )
  )
}

cfa_residual_variance_summary <- tibble::tibble(
  residual_variance_n = nrow(
    cfa_standardized_residual_variances
  ),
  minimum_standardized_residual_variance = min(
    cfa_standardized_residual_variances$
      standardized_variance
  ),
  maximum_standardized_residual_variance = max(
    cfa_standardized_residual_variances$
      standardized_variance
  ),
  all_finite = all(
    is.finite(
      cfa_standardized_residual_variances$
        standardized_variance
    )
  ),
  all_positive = all(
    cfa_standardized_residual_variances$
      standardized_variance > 0
  ),
  all_within_admissible_bounds = all(
    cfa_standardized_residual_variances$
      standardized_variance <=
      (1 + variance_tolerance)
  ),
  heywood_type_result_present = FALSE
)

cfa_computational_diagnostics <- tibble::tibble(
  model = "one_factor",
  sample = "validation",
  sample_n = cfa_observation_n,
  free_parameter_n = cfa_free_parameter_n,
  requested_estimator = (
    cfa_resolved_options$estimator.orig
  ),
  internal_estimator = (
    cfa_resolved_options$estimator
  ),
  parameterization = (
    cfa_resolved_options$parameterization
  ),
  std_lv = cfa_resolved_options$std.lv,
  standard_error_method = (
    cfa_resolved_options$se
  ),
  test_methods = paste(
    cfa_resolved_options$test,
    collapse = "; "
  ),
  information = paste(
    cfa_resolved_options$information,
    collapse = "; "
  ),
  information_meat = paste(
    cfa_resolved_options$information.meat,
    collapse = "; "
  ),
  fitting_warning_n = length(cfa_warnings),
  residual_warning_n = length(
    cfa_residual_warnings
  ),
  converged = cfa_converged,
  post_estimation_check_passed = (
    cfa_post_check
  ),
  loadings_checked = all(
    cfa_loading_checks
  ),
  thresholds_checked = all(
    cfa_threshold_checks
  ),
  residuals_checked = all(
    cfa_residual_checks
  ),
  variances_checked = all(
    cfa_variance_checks
  ),
  improper_solution_present = FALSE
)

# 27. Record the validation-stage dimensionality conclusion --------------------

validation_dimensionality_conclusion <- tibble::tibble(
  conclusion_category = "mixed_evidence",
  conclusion_label = "Mixed evidence",
  decision_date = "2026-07-26",
  primary_validation_model = "one_factor",
  secondary_cfa_fitted = FALSE,
  minimum_standardized_loading = min(
    cfa_standardized_loadings$
      standardized_loading
  ),
  maximum_standardized_loading = max(
    cfa_standardized_loadings$
      standardized_loading
  ),
  robust_cfi = unname(
    cfa_fit_measures_all["cfi.robust"]
  ),
  robust_tli = unname(
    cfa_fit_measures_all["tli.robust"]
  ),
  robust_rmsea = unname(
    cfa_fit_measures_all["rmsea.robust"]
  ),
  robust_rmsea_lower = unname(
    cfa_fit_measures_all[
      "rmsea.ci.lower.robust"
    ]
  ),
  robust_rmsea_upper = unname(
    cfa_fit_measures_all[
      "rmsea.ci.upper.robust"
    ]
  ),
  srmr = unname(
    cfa_fit_measures_all["srmr"]
  ),
  largest_absolute_residual_correlation = max(
    cfa_residual_correlations$
      absolute_residual_correlation
  ),
  residual_correlations_ge_0_10 = sum(
    cfa_residual_correlations$
      flagged_ge_0_10
  ),
  interpretation = paste(
    "A substantial general PHQ-9 factor was supported by",
    "consistently strong item loadings. However, the",
    "prespecified one-factor ordinal CFA showed notable global",
    "misfit, indicating that a single factor did not fully",
    "reproduce the relationships among the nine items. The",
    "evidence was therefore classified as mixed rather than as",
    "proof of strict unidimensionality. Development-sample",
    "multifactor solutions did not provide a stable,",
    "substantively coherent simple structure eligible for",
    "confirmatory validation."
  )
)

validation_analysis_boundaries <- tibble::tibble(
  secondary_cfa_fitted = FALSE,
  post_hoc_correlated_residuals_added = FALSE,
  cfa_cross_loadings_added = FALSE,
  items_deleted = FALSE,
  bifactor_model_fitted = FALSE,
  alternative_model_frozen_retrospectively = FALSE
)

# 28. Define the complete validated Stage 3 object -----------------------------

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
      dplyr = as.character(
        packageVersion("dplyr")
      ),
      tidyr = as.character(
        packageVersion("tidyr")
      ),
      readr = as.character(
        packageVersion("readr")
      ),
      here = as.character(
        packageVersion("here")
      ),
      tibble = as.character(
        packageVersion("tibble")
      ),
      psych = as.character(
        packageVersion("psych")
      ),
      GPArotation = as.character(
        packageVersion("GPArotation")
      ),
      lavaan = as.character(
        packageVersion("lavaan")
      )
    )
  ),
  data = phq9_stage3_split,
  split_summary = split_summary,
  item_category_frequencies = item_category_frequencies,
  item_floor_summary = item_floor_summary,
  phq9_total_distribution = phq9_total_distribution,
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
    observed_factor_eigenvalues = (
      parallel_observed_eigenvalues
    ),
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
  ),
  development_efa = list(
    models = efa_models,
    warnings = efa_warnings,
    loadings = efa_loading_table,
    loading_separation = efa_loading_separation,
    factor_correlations = efa_factor_correlations,
    reproduced_correlations = (
      efa_reproduced_correlations
    ),
    residual_correlations = (
      efa_residual_correlations
    ),
    diagnostics = efa_model_diagnostics,
    three_factor_primary_counts = (
      three_factor_primary_counts
    ),
    two_factor_decision_evidence = (
      two_factor_decision_evidence
    ),
    estimation_settings = list(
      correlation_matrix = (
        "explicit development polychoric matrix"
      ),
      sample_n = expected_development_n,
      fm = "minres",
      SMC = TRUE,
      smooth = FALSE,
      residuals = TRUE,
      one_factor_rotation = "none",
      multifactor_rotation = "oblimin"
    )
  ),
  development_model_freezing_decision = list(
    candidate_assessment = development_model_decision,
    validation_model_plan = validation_model_plan,
    secondary_model_frozen = FALSE,
    validation_model = "one_factor_only",
    decision_date = "2026-07-24",
    recorded_before_validation_access = TRUE
  ),
  validation_cfa = list(
    validation_data_checks = validation_data_checks,
    category_counts = validation_category_counts,
    model_syntax = one_factor_model,
    fitted_model = cfa_fit,
    fitting_warnings = cfa_warnings,
    resolved_options = cfa_resolved_options,
    resolved_option_checks = resolved_option_checks,
    complete_test_information = (
      cfa_test_information
    ),
    model_test_table = cfa_model_test_table,
    complete_fit_measures = (
      cfa_fit_measures_all
    ),
    requested_fit_measures = (
      cfa_fit_measure_table
    ),
    complete_standardized_solution = (
      cfa_standardized_all
    ),
    standardized_loadings = (
      cfa_standardized_loadings
    ),
    loading_checks = cfa_loading_checks,
    thresholds = cfa_thresholds,
    threshold_order_checks = (
      threshold_order_checks
    ),
    threshold_checks = cfa_threshold_checks,
    complete_residual_object = cfa_residuals,
    residual_summary = cfa_residual_summary,
    observed_correlations = (
      cfa_observed_correlations
    ),
    implied_correlations = (
      cfa_implied_correlations
    ),
    residual_correlations = (
      cfa_residual_correlations
    ),
    residual_checks = cfa_residual_checks,
    residual_convention = (
      "observed minus model-implied correlation"
    ),
    standardized_residual_variances = (
      cfa_standardized_residual_variances
    ),
    residual_variance_summary = (
      cfa_residual_variance_summary
    ),
    theta = cfa_theta,
    latent_covariance = cfa_latent_covariance,
    latent_factor_identification = (
      "Factor variance fixed to 1 under std.lv = TRUE"
    ),
    variance_checks = cfa_variance_checks,
    computational_diagnostics = (
      cfa_computational_diagnostics
    ),
    conclusion = (
      validation_dimensionality_conclusion
    ),
    analysis_boundaries = (
      validation_analysis_boundaries
    ),
    validation_access = list(
      validation_results_accessed = TRUE,
      access_date = "2026-07-26",
      primary_model_fitted = TRUE,
      secondary_model_fitted = FALSE
    )
  )
)

# 29. Save the complete validated Stage 3 outputs ------------------------------

saveRDS(
  stage3_split_object,
  stage3_split_path
)

write_csv(
  split_summary,
  split_summary_path
)
