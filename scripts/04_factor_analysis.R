# ==============================================================================
# Stage 3: Dimensionality analysis
# Project: PHQ-9 psychometric evaluation using NHANES 2017–March 2020
#
# Current implementation scope:
#   - reconstruct the complete-PHQ-9 analytic sample with SEQN retained
#   - create the locked 50:50 development and validation split
#   - validate both split samples
#   - save the recoverable split object
#
# No polychoric analysis, parallel analysis, EFA or CFA is run in this block.
# ==============================================================================

# 1. Install any missing packages ----------------------------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "here",
  "tibble"
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

# 10. Define the validated Stage 3 split object --------------------------------

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
      tibble = as.character(packageVersion("tibble"))
    )
  ),

  data = phq9_stage3_split,
  split_summary = split_summary,
  item_category_frequencies = item_category_frequencies,
  item_floor_summary = item_floor_summary
)

# 10a. Calculate PHQ-9 total-score distribution -------------------------------

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

# 11. Save validated Stage 3 split outputs -------------------------------------

saveRDS(
  stage3_split_object,
  stage3_split_path
)

write_csv(
  split_summary,
  split_summary_path
)
