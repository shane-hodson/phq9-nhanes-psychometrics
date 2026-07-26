# ==============================================================================
# Stage 3: Dimensionality-analysis tables and figures
# Project: PHQ-9 psychometric evaluation using NHANES 2017–March 2020
#
# This script:
#   - reads the validated Stage 3 processed-data object
#   - exports the prespecified parallel-analysis, EFA, model-decision and CFA
#     tables
#   - creates the prespecified parallel-analysis, EFA-loading and CFA-loading
#     figures
#
# It does not refit the development EFA or validation CFA models.
# ==============================================================================

# 1. Install any missing packages ----------------------------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "here",
  "tibble"
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
library(ggplot2)
library(here)
library(tibble)

# 3. Define input and output paths ---------------------------------------------

stage3_object_path <- here(
  "data",
  "processed",
  "phq9_stage3_split.rds"
)

parallel_analysis_table_path <- here(
  "tables",
  "phq9_parallel_analysis_results.csv"
)

efa_loadings_table_path <- here(
  "tables",
  "phq9_efa_loadings.csv"
)

efa_factor_correlations_table_path <- here(
  "tables",
  "phq9_efa_factor_correlations.csv"
)

efa_diagnostics_table_path <- here(
  "tables",
  "phq9_efa_diagnostics.csv"
)

stage3_model_decision_table_path <- here(
  "tables",
  "phq9_stage3_model_decision.csv"
)

cfa_fit_measures_table_path <- here(
  "tables",
  "phq9_cfa_fit_measures.csv"
)

cfa_standardized_loadings_table_path <- here(
  "tables",
  "phq9_cfa_standardized_loadings.csv"
)

cfa_thresholds_table_path <- here(
  "tables",
  "phq9_cfa_thresholds.csv"
)

cfa_residual_correlations_table_path <- here(
  "tables",
  "phq9_cfa_residual_correlations.csv"
)

cfa_computational_diagnostics_table_path <- here(
  "tables",
  "phq9_cfa_computational_diagnostics.csv"
)

parallel_analysis_figure_path <- here(
  "figures",
  "phq9_parallel_analysis.png"
)

efa_loadings_figure_path <- here(
  "figures",
  "phq9_efa_loadings.png"
)

cfa_standardized_loadings_figure_path <- here(
  "figures",
  "phq9_cfa_standardized_loadings.png"
)

if (!file.exists(stage3_object_path)) {
  stop(
    "The validated Stage 3 object was not found at: ",
    stage3_object_path
  )
}

dir.create(
  here("tables"),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  here("figures"),
  showWarnings = FALSE,
  recursive = TRUE
)

# 4. Load and validate the Stage 3 object --------------------------------------

stage3_output_object <- readRDS(
  stage3_object_path
)

required_top_level_components <- c(
  "metadata",
  "development_parallel_analysis",
  "development_efa",
  "development_model_freezing_decision",
  "validation_cfa"
)

missing_top_level_components <- setdiff(
  required_top_level_components,
  names(stage3_output_object)
)

if (length(missing_top_level_components) > 0L) {
  stop(
    "The Stage 3 object is missing the following required components: ",
    paste(
      missing_top_level_components,
      collapse = ", "
    )
  )
}

parallel_comparison <- (
  stage3_output_object$
    development_parallel_analysis$
    comparison
)

efa_loadings <- (
  stage3_output_object$
    development_efa$
    loadings
)

efa_factor_correlations <- (
  stage3_output_object$
    development_efa$
    factor_correlations
)

efa_diagnostics <- (
  stage3_output_object$
    development_efa$
    diagnostics
)

development_model_decision <- (
  stage3_output_object$
    development_model_freezing_decision$
    candidate_assessment
)

validation_model_plan <- (
  stage3_output_object$
    development_model_freezing_decision$
    validation_model_plan
)

cfa_fit_measures <- (
  stage3_output_object$
    validation_cfa$
    requested_fit_measures
)

cfa_model_tests <- (
  stage3_output_object$
    validation_cfa$
    model_test_table
)

cfa_standardized_loadings <- (
  stage3_output_object$
    validation_cfa$
    standardized_loadings
)

cfa_thresholds <- (
  stage3_output_object$
    validation_cfa$
    thresholds
)

cfa_residual_correlations <- (
  stage3_output_object$
    validation_cfa$
    residual_correlations
)

cfa_residual_variance_summary <- (
  stage3_output_object$
    validation_cfa$
    residual_variance_summary
)

cfa_computational_diagnostics <- (
  stage3_output_object$
    validation_cfa$
    computational_diagnostics
)

dimensionality_conclusion <- (
  stage3_output_object$
    validation_cfa$
    conclusion
)

validation_analysis_boundaries <- (
  stage3_output_object$
    validation_cfa$
    analysis_boundaries
)

stage3_source_checks <- c(
  parallel_comparison_dimensions = identical(
    dim(parallel_comparison),
    c(9L, 4L)
  ),
  parallel_factor_numbers_correct = identical(
    parallel_comparison$factor_number,
    1:9
  ),
  parallel_recommendation_three = identical(
    sum(
      parallel_comparison$exceeds_threshold
    ),
    3L
  ),
  efa_loading_dimensions = identical(
    dim(efa_loadings),
    c(54L, 9L)
  ),
  efa_solutions_complete = setequal(
    unique(efa_loadings$solution),
    c(
      "one_factor",
      "two_factor",
      "three_factor"
    )
  ),
  efa_factor_correlation_dimensions = identical(
    dim(efa_factor_correlations),
    c(4L, 4L)
  ),
  efa_diagnostic_dimensions = identical(
    dim(efa_diagnostics),
    c(3L, 24L)
  ),
  all_efa_checks_passed = all(
    efa_diagnostics$
      computational_checks_passed
  ),
  two_development_candidates = (
    nrow(development_model_decision) == 2L
  ),
  no_candidate_retained = !any(
    development_model_decision$
      retained_for_validation
  ),
  one_factor_primary_model = identical(
    validation_model_plan$
      primary_validation_model,
    "one_factor"
  ),
  no_secondary_model_frozen = identical(
    validation_model_plan$
      secondary_model_frozen,
    FALSE
  ),
  nine_cfa_fit_measures = (
    nrow(cfa_fit_measures) == 9L
  ),
  two_cfa_model_tests = (
    nrow(cfa_model_tests) == 2L
  ),
  scaled_shifted_test_present = (
    sum(
      cfa_model_tests$test ==
        "scaled_and_shifted"
    ) == 1L
  ),
  nine_cfa_loadings = (
    nrow(cfa_standardized_loadings) == 9L
  ),
  twenty_seven_thresholds = (
    nrow(cfa_thresholds) == 27L
  ),
  thirty_six_residual_pairs = (
    nrow(cfa_residual_correlations) == 36L
  ),
  residuals_sorted_descending = all(
    diff(
      cfa_residual_correlations$
        absolute_residual_correlation
    ) <= 0
  ),
  one_residual_pair_flagged = (
    sum(
      cfa_residual_correlations$
        flagged_ge_0_10
    ) == 1L
  ),
  residual_variance_summary_present = (
    nrow(cfa_residual_variance_summary) == 1L
  ),
  cfa_diagnostics_present = (
    nrow(cfa_computational_diagnostics) == 1L
  ),
  mixed_evidence_recorded = identical(
    dimensionality_conclusion$
      conclusion_category,
    "mixed_evidence"
  ),
  no_secondary_cfa_fitted = identical(
    dimensionality_conclusion$
      secondary_cfa_fitted,
    FALSE
  ),
  all_analysis_boundaries_respected = !any(
    unlist(
      validation_analysis_boundaries,
      use.names = FALSE
    )
  ),
  lavaan_version_available = (
    "lavaan" %in%
      names(
        stage3_output_object$
          metadata$
          package_versions
      )
  )
)

if (!all(stage3_source_checks)) {
  failed_source_checks <- names(
    stage3_source_checks
  )[
    !stage3_source_checks
  ]

  stop(
    "The Stage 3 export sources failed the following checks: ",
    paste(
      failed_source_checks,
      collapse = ", "
    )
  )
}

# 5. Prepare and validate the public Stage 3 tables ----------------------------

parallel_analysis_export <- parallel_comparison |>
  transmute(
    sample = "development",
    sample_n = as.integer(
      stage3_output_object$
        metadata$
        development_sample_n
    ),
    factor_number = factor_number,
    observed_factor_eigenvalue = unname(
      observed_eigenvalue
    ),
    simulated_95th_percentile = unname(
      threshold_95th
    ),
    observed_exceeds_95th_percentile = as.logical(
      exceeds_threshold
    ),
    recommended_factor_count = as.integer(
      stage3_output_object$
        development_parallel_analysis$
        recommended_factors
    ),
    simulation_iterations = as.integer(
      stage3_output_object$
        development_parallel_analysis$
        estimation_settings$
        iterations
    ),
    recommendation_rule = paste(
      "Leading consecutive observed factor eigenvalues",
      "exceeding the simulated 95th-percentile values"
    )
  )

efa_loadings_export <- efa_loadings |>
  group_by(
    solution,
    item
  ) |>
  mutate(
    absolute_loading = abs(loading),
    loading_rank_within_item = min_rank(
      desc(absolute_loading)
    ),
    primary_loading = (
      absolute_loading ==
        max(absolute_loading)
    ),
    loading_ge_0_30 = (
      absolute_loading >= 0.30
    ),
    substantial_secondary_loading = (
      !primary_loading &
        absolute_loading >= 0.30
    )
  ) |>
  ungroup() |>
  arrange(
    factors,
    item,
    factor
  ) |>
  select(
    solution,
    factors,
    rotation,
    item,
    factor,
    loading,
    absolute_loading,
    loading_rank_within_item,
    primary_loading,
    loading_ge_0_30,
    substantial_secondary_loading,
    communality,
    uniqueness,
    complexity
  )

efa_factor_correlations_export <- (
  efa_factor_correlations |>
    mutate(
      absolute_correlation = abs(
        correlation
      )
    ) |>
    arrange(
      solution,
      factor_1,
      factor_2
    )
)

efa_diagnostics_export <- efa_diagnostics |>
  arrange(factors)

development_decision_export <- (
  development_model_decision |>
    transmute(
      decision_stage = (
        "development_model_freezing"
      ),
      model = candidate_solution,
      factor_count = factor_count,
      is_primary_validation_model = FALSE,
      supported_by_parallel_analysis = as.logical(
        supported_by_parallel_analysis
      ),
      retained_for_validation = (
        retained_for_validation
      ),
      secondary_cfa_fitted = NA,
      conclusion_category = NA_character_,
      conclusion_label = NA_character_,
      decision_date = decision_date,
      decision_reason = decision_reason
    )
)

validation_plan_export <- validation_model_plan |>
  transmute(
    decision_stage = "validation_model_plan",
    model = primary_validation_model,
    factor_count = 1L,
    is_primary_validation_model = TRUE,
    supported_by_parallel_analysis = NA,
    retained_for_validation = TRUE,
    secondary_cfa_fitted = FALSE,
    conclusion_category = NA_character_,
    conclusion_label = NA_character_,
    decision_date = decision_date,
    decision_reason = paste(
      "Prespecified compulsory one-factor validation model.",
      "No secondary multifactor model was frozen."
    )
  )

validation_conclusion_export <- (
  dimensionality_conclusion |>
    transmute(
      decision_stage = (
        "validation_conclusion"
      ),
      model = primary_validation_model,
      factor_count = 1L,
      is_primary_validation_model = TRUE,
      supported_by_parallel_analysis = NA,
      retained_for_validation = TRUE,
      secondary_cfa_fitted = (
        secondary_cfa_fitted
      ),
      conclusion_category = (
        conclusion_category
      ),
      conclusion_label = conclusion_label,
      decision_date = decision_date,
      decision_reason = interpretation
    )
)

stage3_model_decision_export <- bind_rows(
  development_decision_export,
  validation_plan_export,
  validation_conclusion_export
)

cfa_fit_measures_export <- cfa_fit_measures |>
  mutate(
    measure_label = case_when(
      measure == "chisq.scaled" ~
        "Scaled-and-shifted chi-square",
      measure == "df.scaled" ~
        "Scaled-and-shifted degrees of freedom",
      measure == "pvalue.scaled" ~
        "Scaled-and-shifted p-value",
      measure == "cfi.robust" ~
        "Robust CFI",
      measure == "tli.robust" ~
        "Robust TLI",
      measure == "rmsea.robust" ~
        "Robust RMSEA",
      measure == "rmsea.ci.lower.robust" ~
        "Robust RMSEA 90% CI lower bound",
      measure == "rmsea.ci.upper.robust" ~
        "Robust RMSEA 90% CI upper bound",
      measure == "srmr" ~
        "SRMR",
      TRUE ~ measure
    ),
    measure_family = case_when(
      measure %in% c(
        "chisq.scaled",
        "df.scaled",
        "pvalue.scaled"
      ) ~ "model_test",
      measure %in% c(
        "cfi.robust",
        "tli.robust"
      ) ~ "incremental_fit",
      measure %in% c(
        "rmsea.robust",
        "rmsea.ci.lower.robust",
        "rmsea.ci.upper.robust"
      ) ~ "error_of_approximation",
      measure == "srmr" ~
        "standardized_residual",
      TRUE ~ "other"
    ),
    confidence_level = if_else(
      measure %in% c(
        "rmsea.ci.lower.robust",
        "rmsea.ci.upper.robust"
      ),
      0.90,
      NA_real_
    ),
    display_order = match(
      measure,
      c(
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
    )
  ) |>
  arrange(display_order) |>
  select(
    display_order,
    measure,
    measure_label,
    measure_family,
    value,
    confidence_level
  )

cfa_standardized_loadings_export <- (
  cfa_standardized_loadings |>
    mutate(
      confidence_interval_type = (
        "95% Wald-type"
      ),
      confidence_level = 0.95
    ) |>
    arrange(
      match(
        item,
        stage3_output_object$
          metadata$
          phq9_items
      )
    )
)

cfa_thresholds_export <- cfa_thresholds |>
  mutate(
    confidence_interval_type = (
      "95% Wald-type"
    ),
    confidence_level = 0.95
  ) |>
  arrange(
    match(
      item,
      stage3_output_object$
        metadata$
        phq9_items
    ),
    threshold
  )

cfa_residual_correlations_export <- (
  cfa_residual_correlations |>
    mutate(
      pair_rank = row_number()
    ) |>
    select(
      pair_rank,
      everything()
    )
)

standard_model_test <- cfa_model_tests |>
  filter(
    test == "standard"
  )

scaled_shifted_model_test <- cfa_model_tests |>
  filter(
    test == "scaled_and_shifted"
  )

if (
  nrow(standard_model_test) != 1L ||
  nrow(scaled_shifted_model_test) != 1L
) {
  stop(
    "The stored CFA model-test table does not contain exactly ",
    "one standard and one scaled-and-shifted test."
  )
}

cfa_computational_diagnostics_export <- (
  cfa_computational_diagnostics |>
    bind_cols(
      tibble(
        standard_test_statistic = (
          standard_model_test$statistic
        ),
        standard_test_degrees_of_freedom = (
          standard_model_test$
            degrees_of_freedom
        ),
        standard_test_p_value = (
          standard_model_test$p_value
        ),
        scaled_shifted_test_statistic = (
          scaled_shifted_model_test$
            statistic
        ),
        scaled_shifted_degrees_of_freedom = (
          scaled_shifted_model_test$
            degrees_of_freedom
        ),
        scaled_shifted_p_value = (
          scaled_shifted_model_test$
            p_value
        ),
        scaling_factor = (
          scaled_shifted_model_test$
            scaling_factor
        ),
        shift_parameter = (
          scaled_shifted_model_test$
            shift_parameter
        ),
        underlying_standard_statistic = (
          scaled_shifted_model_test$
            underlying_standard_statistic
        ),
        correction_label = (
          scaled_shifted_model_test$
            correction_label
        ),
        residual_convention = (
          stage3_output_object$
            validation_cfa$
            residual_convention
        ),
        standardized_residual_variance_n = (
          cfa_residual_variance_summary$
            residual_variance_n
        ),
        minimum_standardized_residual_variance = (
          cfa_residual_variance_summary$
            minimum_standardized_residual_variance
        ),
        maximum_standardized_residual_variance = (
          cfa_residual_variance_summary$
            maximum_standardized_residual_variance
        ),
        residual_variances_all_finite = (
          cfa_residual_variance_summary$
            all_finite
        ),
        residual_variances_all_positive = (
          cfa_residual_variance_summary$
            all_positive
        ),
        residual_variances_within_bounds = (
          cfa_residual_variance_summary$
            all_within_admissible_bounds
        ),
        heywood_type_result_present = (
          cfa_residual_variance_summary$
            heywood_type_result_present
        ),
        latent_factor_identification = (
          stage3_output_object$
            validation_cfa$
            latent_factor_identification
        ),
        conclusion_category = (
          dimensionality_conclusion$
            conclusion_category
        ),
        conclusion_label = (
          dimensionality_conclusion$
            conclusion_label
        ),
        lavaan_version = unname(
          stage3_output_object$
            metadata$
            package_versions[
              "lavaan"
            ]
        )
      )
    ) |>
    bind_cols(
      validation_analysis_boundaries
    )
)

stage3_export_table_checks <- c(
  parallel_table_has_nine_rows = (
    nrow(parallel_analysis_export) == 9L
  ),
  parallel_table_recommends_three = (
    unique(
      parallel_analysis_export$
        recommended_factor_count
    ) == 3L
  ),
  efa_loading_table_has_fifty_four_rows = (
    nrow(efa_loadings_export) == 54L
  ),
  one_primary_loading_per_item_solution = all(
    efa_loadings_export |>
      group_by(
        solution,
        item
      ) |>
      summarise(
        primary_loading_n = sum(
          primary_loading
        ),
        .groups = "drop"
      ) |>
      pull(primary_loading_n) == 1L
  ),
  efa_factor_correlation_table_has_four_rows = (
    nrow(
      efa_factor_correlations_export
    ) == 4L
  ),
  efa_diagnostic_table_has_three_rows = (
    nrow(efa_diagnostics_export) == 3L
  ),
  model_decision_table_has_four_rows = (
    nrow(stage3_model_decision_export) == 4L
  ),
  no_multifactor_candidate_retained = !any(
    stage3_model_decision_export$
      retained_for_validation[
        stage3_model_decision_export$
          decision_stage ==
          "development_model_freezing"
      ]
  ),
  mixed_evidence_conclusion_present = (
    sum(
      stage3_model_decision_export$
        conclusion_category ==
        "mixed_evidence",
      na.rm = TRUE
    ) == 1L
  ),
  cfa_fit_measure_table_has_nine_rows = (
    nrow(cfa_fit_measures_export) == 9L
  ),
  all_cfa_fit_values_finite = all(
    is.finite(
      cfa_fit_measures_export$value
    )
  ),
  cfa_loading_table_has_nine_rows = (
    nrow(
      cfa_standardized_loadings_export
    ) == 9L
  ),
  cfa_threshold_table_has_twenty_seven_rows = (
    nrow(cfa_thresholds_export) == 27L
  ),
  cfa_residual_table_has_thirty_six_rows = (
    nrow(
      cfa_residual_correlations_export
    ) == 36L
  ),
  cfa_residual_table_sorted = all(
    diff(
      cfa_residual_correlations_export$
        absolute_residual_correlation
    ) <= 0
  ),
  cfa_residual_table_has_one_flag = (
    sum(
      cfa_residual_correlations_export$
        flagged_ge_0_10
    ) == 1L
  ),
  computational_diagnostic_table_has_one_row = (
    nrow(
      cfa_computational_diagnostics_export
    ) == 1L
  ),
  scaled_shifted_correction_complete = all(
    is.finite(
      c(
        cfa_computational_diagnostics_export$
          scaled_shifted_test_statistic,
        cfa_computational_diagnostics_export$
          scaling_factor,
        cfa_computational_diagnostics_export$
          shift_parameter
      )
    )
  ),
  residual_variance_bounds_recorded = (
    cfa_computational_diagnostics_export$
      minimum_standardized_residual_variance >
      0 &&
      cfa_computational_diagnostics_export$
      maximum_standardized_residual_variance <=
      1
  ),
  no_heywood_type_result = !(
    cfa_computational_diagnostics_export$
      heywood_type_result_present
  ),
  all_post_hoc_boundaries_false = !any(
    unlist(
      validation_analysis_boundaries,
      use.names = FALSE
    )
  )
)

if (!all(stage3_export_table_checks)) {
  failed_export_table_checks <- names(
    stage3_export_table_checks
  )[
    !stage3_export_table_checks
  ]

  stop(
    "The prepared Stage 3 tables failed the following checks: ",
    paste(
      failed_export_table_checks,
      collapse = ", "
    )
  )
}

# 6. Prepare and validate the Stage 3 figures ----------------------------------

phq9_item_order <- (
  stage3_output_object$
    metadata$
    phq9_items
)

# Parallel-analysis figure data

parallel_long_data <- (
  parallel_analysis_export |>
    select(
      factor_number,
      observed_factor_eigenvalue,
      simulated_95th_percentile
    ) |>
    pivot_longer(
      cols = c(
        observed_factor_eigenvalue,
        simulated_95th_percentile
      ),
      names_to = "series",
      values_to = "factor_eigenvalue"
    ) |>
    mutate(
      series = recode(
        series,
        observed_factor_eigenvalue = (
          "Observed factor eigenvalue"
        ),
        simulated_95th_percentile = (
          "Simulated 95th percentile"
        )
      )
    )
)

parallel_figure_data <- bind_rows(
  parallel_long_data |>
    mutate(
      display_panel = "All factors"
    ),
  parallel_long_data |>
    filter(
      factor_number >= 2L
    ) |>
    mutate(
      display_panel = (
        "Factors 2–9, expanded scale"
      )
    )
) |>
  mutate(
    display_panel = factor(
      display_panel,
      levels = c(
        "All factors",
        "Factors 2–9, expanded scale"
      )
    ),
    series = factor(
      series,
      levels = c(
        "Observed factor eigenvalue",
        "Simulated 95th percentile"
      )
    )
  )

parallel_analysis_figure <- ggplot(
  parallel_figure_data,
  aes(
    x = factor_number,
    y = factor_eigenvalue,
    group = series,
    linetype = series,
    shape = series
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.3
  ) +
  facet_wrap(
    vars(display_panel),
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = 1:9
  ) +
  labs(
    title = "Ordinal parallel analysis of the PHQ-9",
    subtitle = paste(
      "Development sample, n =",
      stage3_output_object$
        metadata$
        development_sample_n
    ),
    x = "Factor number",
    y = "Factor eigenvalue",
    linetype = NULL,
    shape = NULL,
    caption = paste(
      "Three factors were retained because the first three consecutive",
      "observed factor eigenvalues exceeded the corresponding simulated",
      "95th-percentile values. The lower panel enlarges Factors 2–9.",
      sep = "\n"
    )
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust = 0,
      margin = margin(
        t = 8
      )
    ),
    plot.margin = margin(
      t = 10,
      r = 12,
      b = 12,
      l = 12
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

# EFA-loading figure data

efa_panel_levels <- c(
  "One-factor solution: MR1",
  "Two-factor solution: MR1",
  "Two-factor solution: MR2",
  "Three-factor solution: MR1",
  "Three-factor solution: MR2",
  "Three-factor solution: MR3"
)

efa_figure_data <- efa_loadings_export |>
  mutate(
    item = factor(
      item,
      levels = rev(phq9_item_order)
    ),
    solution_label = case_when(
      solution == "one_factor" ~
        "One-factor solution",
      solution == "two_factor" ~
        "Two-factor solution",
      solution == "three_factor" ~
        "Three-factor solution",
      TRUE ~ solution
    ),
    panel_label = paste(
      solution_label,
      factor,
      sep = ": "
    ),
    panel_label = factor(
      panel_label,
      levels = efa_panel_levels
    )
  )

efa_loadings_figure <- ggplot(
  efa_figure_data,
  aes(
    x = item,
    y = loading
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = c(
      -0.30,
      0.30
    ),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_col(
    width = 0.72
  ) +
  coord_flip() +
  facet_wrap(
    vars(panel_label),
    ncol = 3
  ) +
  scale_y_continuous(
    limits = c(
      -1,
      1
    ),
    breaks = seq(
      -1,
      1,
      by = 0.5
    )
  ) +
  labs(
    title = "Development-sample exploratory factor loadings",
    subtitle = paste(
      "MINRES extraction using the audited polychoric matrix;",
      "multifactor solutions use oblimin rotation"
    ),
    x = NULL,
    y = "Pattern loading",
    caption = paste(
      "Dashed reference lines mark absolute loadings of .30.",
      "The multifactor patterns were exploratory and were not",
      "carried forward as validation models."
    )
  ) +
  theme_minimal(
    base_size = 10
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    strip.text = element_text(
      face = "bold",
      size = 9
    )
  )

# CFA-loading figure data

cfa_loading_figure_data <- (
  cfa_standardized_loadings_export |>
    mutate(
      item = factor(
        item,
        levels = rev(phq9_item_order)
      )
    )
)

cfa_standardized_loadings_figure <- ggplot(
  cfa_loading_figure_data,
  aes(
    x = standardized_loading,
    y = item
  )
) +
  geom_errorbar(
    aes(
      xmin = ci_lower,
      xmax = ci_upper
    ),
    orientation = "y",
    width = 0.15
  ) +
  geom_point(
    size = 2.7
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.1
    )
  ) +
  labs(
    title = "One-factor ordinal CFA loadings",
    subtitle = paste(
      "Validation sample, n =",
      stage3_output_object$
        metadata$
        validation_sample_n
    ),
    x = "Fully standardised loading",
    y = NULL,
    caption = paste(
      "Error bars show 95% Wald-type confidence intervals.",
      "Strong loadings support a substantial general factor,",
      "but do not establish strict unidimensionality.",
      sep = "\n"
    )
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.caption = element_text(
      hjust = 0
    )
  )
stage3_figure_checks <- c(
  parallel_figure_data_has_expected_rows = (
    nrow(parallel_figure_data) == 34L
  ),
  parallel_figure_has_two_panels = (
    n_distinct(
      parallel_figure_data$
        display_panel
    ) == 2L
  ),
  efa_figure_has_fifty_four_rows = (
    nrow(efa_figure_data) == 54L
  ),
  efa_figure_has_six_panels = (
    n_distinct(
      efa_figure_data$
        panel_label
    ) == 6L
  ),
  cfa_figure_has_nine_items = (
    nrow(cfa_loading_figure_data) == 9L
  ),
  cfa_intervals_are_ordered = all(
    cfa_loading_figure_data$ci_lower <=
      cfa_loading_figure_data$
      standardized_loading &
      cfa_loading_figure_data$
      standardized_loading <=
      cfa_loading_figure_data$ci_upper
  ),
  parallel_figure_is_ggplot = inherits(
    parallel_analysis_figure,
    "ggplot"
  ),
  efa_figure_is_ggplot = inherits(
    efa_loadings_figure,
    "ggplot"
  ),
  cfa_figure_is_ggplot = inherits(
    cfa_standardized_loadings_figure,
    "ggplot"
  )
)

if (!all(stage3_figure_checks)) {
  failed_figure_checks <- names(
    stage3_figure_checks
  )[
    !stage3_figure_checks
  ]

  stop(
    "The prepared Stage 3 figures failed the following checks: ",
    paste(
      failed_figure_checks,
      collapse = ", "
    )
  )
}

# 7. Export and validate the Stage 3 outputs -----------------------------------

write_csv(
  parallel_analysis_export,
  parallel_analysis_table_path
)

write_csv(
  efa_loadings_export,
  efa_loadings_table_path
)

write_csv(
  efa_factor_correlations_export,
  efa_factor_correlations_table_path
)

write_csv(
  efa_diagnostics_export,
  efa_diagnostics_table_path
)

write_csv(
  stage3_model_decision_export,
  stage3_model_decision_table_path
)

write_csv(
  cfa_fit_measures_export,
  cfa_fit_measures_table_path
)

write_csv(
  cfa_standardized_loadings_export,
  cfa_standardized_loadings_table_path
)

write_csv(
  cfa_thresholds_export,
  cfa_thresholds_table_path
)

write_csv(
  cfa_residual_correlations_export,
  cfa_residual_correlations_table_path
)

write_csv(
  cfa_computational_diagnostics_export,
  cfa_computational_diagnostics_table_path
)

ggsave(
  filename = parallel_analysis_figure_path,
  plot = parallel_analysis_figure,
  width = 8.5,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = efa_loadings_figure_path,
  plot = efa_loadings_figure,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = cfa_standardized_loadings_figure_path,
  plot = cfa_standardized_loadings_figure,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

stage3_table_paths <- c(
  parallel_analysis_table_path,
  efa_loadings_table_path,
  efa_factor_correlations_table_path,
  efa_diagnostics_table_path,
  stage3_model_decision_table_path,
  cfa_fit_measures_table_path,
  cfa_standardized_loadings_table_path,
  cfa_thresholds_table_path,
  cfa_residual_correlations_table_path,
  cfa_computational_diagnostics_table_path
)

stage3_figure_paths <- c(
  parallel_analysis_figure_path,
  efa_loadings_figure_path,
  cfa_standardized_loadings_figure_path
)

stage3_output_paths <- c(
  stage3_table_paths,
  stage3_figure_paths
)

stage3_output_file_information <- file.info(
  stage3_output_paths
)

stage3_written_output_checks <- c(
  ten_tables_exist = (
    sum(
      file.exists(stage3_table_paths)
    ) == 10L
  ),
  three_figures_exist = (
    sum(
      file.exists(stage3_figure_paths)
    ) == 3L
  ),
  all_output_files_exist = all(
    file.exists(stage3_output_paths)
  ),
  all_output_files_nonempty = all(
    is.finite(
      stage3_output_file_information$size
    ) &
      stage3_output_file_information$size > 0
  ),
  all_table_extensions_csv = all(
    tools::file_ext(
      stage3_table_paths
    ) == "csv"
  ),
  all_figure_extensions_png = all(
    tools::file_ext(
      stage3_figure_paths
    ) == "png"
  )
)

if (!all(stage3_written_output_checks)) {
  failed_written_output_checks <- names(
    stage3_written_output_checks
  )[
    !stage3_written_output_checks
  ]

  stop(
    "The Stage 3 output files failed the following checks: ",
    paste(
      failed_written_output_checks,
      collapse = ", "
    )
  )
}

cat(
  "\nStage 3 tables exported successfully:\n",
  paste0(
    "- ",
    stage3_table_paths,
    collapse = "\n"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Stage 3 figures exported successfully:\n",
  paste0(
    "- ",
    stage3_figure_paths,
    collapse = "\n"
  ),
  "\n",
  sep = ""
)
