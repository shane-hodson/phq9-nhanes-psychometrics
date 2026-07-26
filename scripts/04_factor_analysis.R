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
#   - fit and validate the permitted development-sample EFA solutions
#   - extract EFA loading, factor-correlation, residual and fit diagnostics
#   - record the development-stage model-freezing decision
#   - save the recoverable split, polychoric, parallel-analysis and EFA objects
#
# Parallel analysis supports three factors.
# Neither multifactor EFA satisfied the locked criteria for freezing a
# secondary validation CFA.
# No validation-sample CFA is run in the current implementation.
# ==============================================================================

# 1. Install any missing packages ----------------------------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "here",
  "tibble",
  "psych",
  "GPArotation"
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

# 20. Define the validated Stage 3 object --------------------------------------

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
  )
)

# 21. Save the validated Stage 3 outputs ---------------------------------------

saveRDS(
  stage3_split_object,
  stage3_split_path
)

write_csv(
  split_summary,
  split_summary_path
)
