# Analysis Plan

**Repository:** `phq9-nhanes-psychometrics`
**Analysis plan version:** 1.7
**Last updated:** 26 July 2026
**Current stage:** The pre-specified Stage 3 development and validation analyses are complete. No secondary multifactor CFA was frozen or fitted. The validation-sample one-factor ordinal CFA produced mixed evidence: all nine items loaded strongly on a substantial general factor, but global fit indicated notable departure from a complete one-factor representation. Stage 3 output export, report integration and final repository review remain outstanding.

## 1. Project title

**Psychometric Evaluation of the PHQ-9 in U.S. Adults Using NHANES 2017–March 2020**

**Subtitle:** *Reliability, Dimensionality, Functional Difficulty, and Clinical Interpretation*

## 2. Dataset and cycle

The project uses the publicly available **National Health and Nutrition Examination Survey 2017–March 2020 pre-pandemic release**.

The primary files are:

- `P_DPQ.xpt`: Mental Health—Depression Screener
- `P_DEMO.xpt`: demographic and survey-design variables

Additional questionnaire files may be imported only after a separate external-variable audit.

The completed 2017–2018 cycle and the partial 2019–March 2020 cycle are treated as one combined pre-pandemic release. The partial 2019–March 2020 data will not be analysed independently.

## 3. Current analytic sample definitions

Stage 1 established the following sample characteristics:

- `P_DPQ` records: 8,965
- `P_DEMO` records: 15,560
- adults aged 18 years or older in `P_DEMO`: 9,693
- `P_DPQ` records successfully linked to `P_DEMO`: 8,965
- duplicate `SEQN` values: 0 in both files
- unmatched `P_DPQ` records: 0

### Complete-PHQ-9 sample

The primary complete-PHQ-9 sample is defined as adults who:

- have valid responses from 0 to 3 on all nine PHQ-9 items;
- are linked to `P_DEMO` using `SEQN`;
- have valid survey-design variables available for weighted analyses.

This sample contains:

- **8,276 adults with complete PHQ-9 data**

This sample is used for:

- weighted PHQ-9 total-score descriptives;
- weighted symptom-severity band estimates;
- reliability and item analysis;
- polychoric correlations;
- dimensionality analysis.

### Symptom-positive complete-PHQ-9 sample

Among participants with complete PHQ-9 data:

- **5,522 endorsed at least one PHQ-9 symptom**
- **2,754 scored zero on all nine PHQ-9 items**

### Primary functional-difficulty sample

The primary functional-difficulty sample is defined as adults who:

- have complete PHQ-9 data;
- endorsed at least one PHQ-9 symptom;
- have a valid `DPQ100` response.

This sample contains:

- **5,517 adults with complete PHQ-9 data, at least one endorsed symptom and valid `DPQ100`**

Five symptom-positive participants with complete PHQ-9 data did not have a valid `DPQ100` response.

Analysis-specific sample sizes will be reported for:

1. item-level descriptive analyses;
2. complete-PHQ-9 reliability and dimensionality analyses;
3. functional-difficulty analyses;
4. any optional external-variable analyses.

All major exclusions and analysis-specific sample sizes will be reported.

## 4. PHQ-9 variables and items

The nine PHQ-9 symptom items are:

| Variable | Symptom content |
|---|---|
| `DPQ010` | Little interest or pleasure |
| `DPQ020` | Depressed mood or hopelessness |
| `DPQ030` | Sleep disturbance |
| `DPQ040` | Tiredness or low energy |
| `DPQ050` | Appetite disturbance |
| `DPQ060` | Negative self-evaluation |
| `DPQ070` | Concentration difficulties |
| `DPQ080` | Psychomotor slowing or agitation |
| `DPQ090` | Thoughts of death or self-harm |

Responses are scored:

- 0: not at all
- 1: several days
- 2: more than half the days
- 3: nearly every day

A complete PHQ-9 total score ranges from 0 to 27.

`DPQ100` measures reported difficulty with work, home responsibilities and relationships. It is not included in the PHQ-9 total score and is treated as the primary external functional variable.

## 5. DPQ100 routing decision

Stage 1 found that:

- 2,754 participants had valid zero responses on all nine PHQ-9 items;
- all 2,754 had `DPQ100` recorded as system missing;
- 5,522 participants had complete PHQ-9 data and endorsed at least one symptom;
- 5,517 of these participants had a valid `DPQ100` response;
- five symptom-positive participants did not have a valid `DPQ100` response.

The observed pattern is consistent with structural questionnaire routing among participants reporting no PHQ-9 symptoms.

Structurally skipped `DPQ100` values will not be recoded as “not difficult at all.”

The primary functional-difficulty analysis will use participants with:

- complete PHQ-9 data;
- at least one endorsed PHQ-9 symptom;
- valid `DPQ100`.

For extended descriptive tables and figures, structurally skipped responses are shown as a separate category:

> No PHQ-9 symptoms / DPQ100 structurally skipped

The five symptom-positive participants without a valid `DPQ100` response are also shown separately where complete sample accounting is required.

Neither structural skips nor missing responses are treated as observed `DPQ100` response categories.

## 6. External variables

### Required for version 1

- `DPQ100`: functional difficulty

### Secondary candidates

- `HUQ010`: self-rated general health
- `HUQ090`: contact with a mental-health professional during the previous year
- `SLD012`: usual weekday or workday sleep duration
- `SLQ050`: reported trouble sleeping to a health professional
- `SLQ120`: frequency of excessive daytime sleepiness

Only `DPQ100` is required for version 1.

A secondary variable will be included only when:

- its coding and questionnaire routing are sufficiently clear;
- the available sample remains adequate;
- it addresses a defined question about PHQ-9 symptom scores or reported functioning;
- conceptual overlap with PHQ-9 item content is acknowledged.

Sleep variables will be treated as exploratory because sleep disturbance and fatigue are already represented within the PHQ-9.

## 7. Research questions

1. What are the weighted response distributions, endorsement rates and missingness patterns of the nine PHQ-9 items?
2. What is the weighted distribution of PHQ-9 total scores and conventional symptom-severity bands?
3. How consistently do the nine PHQ-9 items measure depressive symptom severity?
4. Are the relationships among the nine ordered PHQ-9 items reasonably compatible with a single general depressive-symptom factor?
5. How are PHQ-9 scores associated with reported functional difficulty?
6. What conclusions about symptom severity and functioning are supported, and which diagnostic or treatment-related conclusions are not supported?

The Stage 3 dimensionality conclusion will concern compatibility with a one-factor representation. It will not be interpreted as proof of strict unidimensionality, construct validity, diagnostic accuracy or treatment need.

## 8. Exclusion and missing-data rules

- Restrict the main sample to participants aged 18 years or older.
- Link files using `SEQN`.
- Check for duplicate and unmatched identifiers before analysis.
- Treat refusal, “don’t know” and system-missing responses as missing values.
- Never treat non-substantive response codes such as 7 or 9 as symptom scores.
- Use all available valid responses for item-specific descriptive analyses.
- Calculate the primary PHQ-9 total only when all nine items have valid scores from 0 to 3.
- Do not prorate incomplete PHQ-9 scores.
- Do not impute missing PHQ-9 items in version 1.
- Use complete PHQ-9 cases for reliability, polychoric-correlation and dimensionality analyses.
- Do not recode structurally skipped `DPQ100` values as zero or as “not difficult at all.”
- Use participants with complete PHQ-9 data, at least one endorsed symptom and valid `DPQ100` for the primary functional-difficulty analysis.
- Use analysis-specific complete cases for secondary external variables.
- Report the number and percentage excluded under each major rule.
- Compare the retained complete-PHQ-9 sample with excluded participants on basic demographic variables where feasible.

## 9. Completed Stage 1 analyses

Stage 1 review and validation are complete.

Completed work includes:

- programmatic download of `P_DPQ.xpt` and `P_DEMO.xpt`;
- source-file inspection;
- `SEQN` linkage verification;
- PHQ-9 variable dictionary;
- response-code audit;
- missingness audit;
- `DPQ100` routing audit;
- sample-flow output;
- survey-design variable audit;
- local recoding of non-substantive questionnaire codes as missing.

Stage 1 established that:

- `P_DPQ` contains 8,965 unique respondents;
- all `P_DPQ` respondents link successfully to `P_DEMO`;
- 8,276 respondents have valid responses to all nine PHQ-9 items;
- 5,522 complete-PHQ-9 respondents endorsed at least one symptom;
- 5,517 respondents have complete PHQ-9 data, at least one endorsed symptom and valid `DPQ100`;
- 2,754 symptom-free respondents have structurally skipped `DPQ100`;
- non-substantive questionnaire responses must be treated as missing before analysis;
- all required survey-design variables are available.

## 10. Completed Stage 2A analyses

Stage 2A weighted descriptive analysis and validation are complete.

### 10.1 Survey design

The NHANES complex survey-design object was created using:

- examination weight: `WTMECPRP`
- masked variance stratum: `SDMVSTRA`
- masked variance primary sampling unit: `SDMVPSU`
- PSU identifiers nested within strata

The design contains 49 nested PSUs across 24 strata.

Analysis-specific samples were defined as survey domains after construction of the full survey-design object.

### 10.2 Completed descriptive outputs

Stage 2A produced:

- weighted PHQ-9 item-response proportions;
- weighted PHQ-9 total-score summaries and distribution;
- weighted conventional symptom-severity band estimates;
- weighted `DPQ100` descriptives among routed respondents;
- an extended functional-difficulty table including structural skips and missing responses as separate categories;
- exported descriptive tables;
- five exported descriptive figures;
- an updated Quarto report.

Stage 2A retained structurally skipped `DPQ100` responses as a separate descriptive category rather than treating them as observed functional-difficulty responses.

### 10.3 Selected Stage 2A results

In the complete-PHQ-9 sample of 8,276 adults:

- weighted mean PHQ-9 total: 3.18
- 95% confidence interval for the weighted mean: 3.06–3.31
- weighted standard deviation: 4.09
- weighted minimal band: 74.92%
- weighted mild band: 16.54%
- weighted moderate band: 5.95%
- weighted moderately severe band: 1.81%
- weighted severe band: 0.78%

Among the 5,517 routed respondents with valid `DPQ100`:

- not difficult at all: 74.64%
- somewhat difficult: 20.79%
- very difficult: 3.40%
- extremely difficult: 1.17%

These results are descriptive. Stage 2A did not estimate reliability, test dimensionality, model the relationship between PHQ-9 scores and functional difficulty, or draw diagnostic conclusions.

## 11. Completed Stage 2B: Reliability and item analysis

Stage 2B used the complete-PHQ-9 sample of 8,276 participants.

The completed analysis included:

- item floor effects and response-category sparsity;
- corrected item-total correlations with Fisher’s z confidence intervals;
- estimation and validation of the polychoric correlation matrix;
- ordinal coefficient alpha;
- McDonald’s omega total from a one-factor congeneric reliability model;
- reliability estimates after removing each item individually;
- validation of the omega uniqueness estimates;
- 2,000 participant-level nonparametric bootstrap resamples;
- percentile bootstrap confidence intervals and computational diagnostics;
- exported reliability and item-analysis tables;
- exported corrected item-total, polychoric-correlation and reliability figures;
- an updated Quarto report section.

Ordinal coefficient alpha was 0.923, with a 95% percentile bootstrap confidence interval of 0.918–0.927.

McDonald’s omega total was 0.923, with a 95% percentile bootstrap confidence interval of 0.919–0.927.

All nine omega uniqueness estimates were finite, non-negative and no greater than 1, with values ranging from 0.236 to 0.525.

Removing any individual item lowered both ordinal alpha and omega total.

All 2,000 bootstrap resamples completed successfully. No resampled polychoric matrix required smoothing, and no bootstrap replicate produced a computational warning.

The estimation method used for each reliability coefficient was documented in the analysis script and Quarto report.

Item-removal statistics were treated as diagnostic information and were not used as an automatic basis for recommending deletion of PHQ-9 items.

Reliability coefficients were reported as evidence of internal consistency, not as evidence that the PHQ-9 is valid or unidimensional.

Stage 2B did not include:

- exploratory or confirmatory factor analysis;
- measurement invariance;
- functional-difficulty regression modelling;
- diagnostic interpretation;
- treatment-need claims;
- substantive clinical interpretation.

## 12. Locked Stage 3 dimensionality-analysis specification

### 12.1 Primary dimensionality question

Stage 3 will evaluate:

> Are the relationships among the nine ordered PHQ-9 items reasonably compatible with a single general depressive-symptom factor?

The conclusion will concern compatibility with a one-factor representation.

A favourable result will not be described as proof that:

- the PHQ-9 is strictly unidimensional;
- the PHQ-9 has established construct validity in this sample;
- participants meet diagnostic criteria for a depressive disorder;
- any participant requires treatment;
- a particular clinical decision is warranted.

### 12.2 Stage 3 analysis sample

Stage 3 will use the existing complete-PHQ-9 analytic sample of 8,276 participants.

Eligibility requires:

- valid integer responses from 0 to 3 on all nine PHQ-9 items;
- successful linkage through `SEQN`;
- no missing PHQ-9 item response.

`SEQN` will be retained for reproducible sample assignment but excluded from all factor models.

The Stage 3 factor analyses will be unweighted. Their parameters will be described as psychometric estimates for the complete-PHQ-9 analytic sample rather than nationally representative latent-variable estimates.

NHANES weights, strata and primary sampling units will not enter the version 1 factor models.

### 12.3 Fixed development and validation split

The complete-PHQ-9 sample will be divided into a fixed 50:50 development and validation split using:

```r
set.seed(20260723)
```

The intended sample sizes are:

- development sample: `n = 4,138`
- validation sample: `n = 4,138`

The split will be a simple random participant-level allocation.

The assignment will be saved with `SEQN` in a local processed-data object so that the exact split remains recoverable.

Before any factor model is fitted, the following will be checked separately in the development and validation samples:

- participant count;
- duplicate `SEQN` values;
- missingness;
- valid item-response range;
- integer response coding;
- presence of all four response categories for every PHQ-9 item;
- item-category frequencies;
- item floor effects;
- PHQ-9 total-score summaries and distribution.

The split will not be rerandomised after its characteristics are inspected.

An absent response category, unexpected missingness, incorrect sample size or other serious data problem will stop execution and require a documented analysis-plan amendment rather than use of another seed.

Validation-sample factor results will not be inspected during the exploratory development stage.

### 12.4 Development-sample polychoric matrix

The development-sample polychoric correlation matrix will be estimated explicitly using:

```r
development_poly <- psych::polychoric(
  development_items,
  correct = 0.5,
  smooth = FALSE,
  global = TRUE,
  progress = FALSE
)
```

The following properties will be checked:

- finite estimates;
- expected dimensions;
- expected item order;
- symmetry;
- unit diagonal;
- correlations within the valid range;
- positive definiteness;
- minimum eigenvalue;
- whether any warning was produced.

Automatic matrix smoothing will not be used for the primary observed development matrix.

A non-positive-definite matrix, invalid estimate or unexplained warning will stop Stage 3 execution for methodological review.

The explicitly audited development matrix will be used for the exploratory factor analyses.

### 12.5 Ordinal parallel analysis

Ordinal factor parallel analysis will be conducted in the development sample using a separate fixed seed:

```r
set.seed(20260724)
```

The analysis will be run in a controlled single-core scope:

```r
parallel_results <- local({
  old_options <- options(mc.cores = 1L)
  on.exit(options(old_options), add = TRUE)

  set.seed(20260724)

  psych::fa.parallel(
    x = development_items,
    fm = "minres",
    fa = "fa",
    nfactors = 1,
    n.iter = 1000,
    SMC = FALSE,
    sim = FALSE,
    quant = 0.95,
    cor = "poly",
    correct = 0.5,
    plot = FALSE
  )
})
```

The locked settings are:

- `fm = "minres"`
- `fa = "fa"`
- `nfactors = 1`
- `n.iter = 1000`
- `SMC = FALSE`
- `sim = FALSE`
- `quant = 0.95`
- `cor = "poly"`
- `correct = 0.5`
- `plot = FALSE`
- `mc.cores = 1L`
- fixed seed `20260724`

With `sim = FALSE`, the null reference distributions will be generated by independently resampling each PHQ-9 item column with replacement. This preserves each item’s marginal category distribution while disrupting inter-item associations.

The primary retention comparison will use factor eigenvalues rather than principal-component eigenvalues.

Each observed factor eigenvalue will be compared with the corresponding 95th percentile of the resampled factor-eigenvalue distribution.

The complete returned object will be saved locally as an excluded processed `.rds` file.

The R version and installed `psych` package version will be recorded.

Because `fa.parallel()` does not expose every matrix-control argument used in the explicit `polychoric()` call, the report will not claim that every internally generated null matrix was estimated with `smooth = FALSE`. The observed development-sample polychoric matrix itself will remain explicitly unsmoothed and separately audited.

### 12.6 Parallel-analysis extraction and retention rule

The percentile thresholds will be reconstructed from the replicate-level factor columns in `parallel_results$values`.

The mean resampled eigenvalues stored in `parallel_results$fa.simr` will not be used as substitutes for the requested 95th-percentile thresholds.

The expected factor columns will be defined as:

```r
factor_columns <- paste0("F", seq_along(phq9_items))
```

The returned object will be checked before extraction:

```r
stopifnot(
  length(parallel_results$fa.values) == length(phq9_items),
  nrow(parallel_results$values) == 1000L,
  all(factor_columns %in% colnames(parallel_results$values))
)
```

The replicate-level factor eigenvalues and their 95th percentiles will be extracted as follows:

```r
resampled_factor_eigenvalues <- parallel_results$values[
  ,
  factor_columns,
  drop = FALSE
]

resampled_q95 <- apply(
  resampled_factor_eigenvalues,
  2,
  stats::quantile,
  probs = 0.95,
  na.rm = TRUE
)

parallel_comparison <- tibble::tibble(
  factor_number = seq_along(phq9_items),
  observed_eigenvalue = parallel_results$fa.values,
  resampled_q95 = unname(resampled_q95),
  retained = observed_eigenvalue > resampled_q95
)
```

The recorded custom recommendation will be the number of leading consecutive factors for which:

```text
observed factor eigenvalue > corresponding resampled 95th-percentile eigenvalue
```

Both sides of every comparison will be exported.

The custom leading-consecutive recommendation will be checked against `parallel_results$nfact`.

A disagreement between the custom result and `parallel_results$nfact` will stop execution rather than be silently resolved.

A non-consecutive crossing pattern, unstable result or recommendation above three factors will be treated as an interpretive ambiguity. It will not justify unrestricted expansion of the model set.

Parallel analysis will be treated as one source of evidence rather than an automatic factor-retention rule because large samples can support minor factors that have limited substantive importance.

### 12.7 Limited exploratory factor analysis

Development-sample exploratory factor analysis will include:

1. a one-factor MINRES solution;
2. a two-factor MINRES solution with oblimin rotation;
3. a three-factor MINRES solution with oblimin rotation only when parallel analysis supports at least three factors.

No four-factor or larger exploratory solution will be fitted in version 1.

The exploratory analyses will use:

- the explicitly estimated development-sample polychoric matrix;
- `fm = "minres"`;
- no rotation for the one-factor solution;
- oblimin rotation for multifactor solutions;
- `n.obs = 4138`;
- the prespecified PHQ-9 item order.

The following will be reported for each fitted solution where applicable:

- pattern-matrix loadings;
- factor correlations;
- communalities;
- uniquenesses;
- item complexity;
- reproduced correlations;
- residual correlations;
- root mean square residual information;
- variance-accounted-for information;
- convergence and warning information.

The following will be treated as descriptive diagnostic flags rather than automatic deletion or retention rules:

- primary loading below `.40`;
- cross-loading of `.30` or greater;
- difference below `.20` between the two largest loadings;
- unusually large residual correlations;
- a factor with fewer than three substantively coherent indicators;
- a factor defined primarily by one unusual item;
- a factor apparently driven by sparse response categories;
- a factor apparently driven by wording or content overlap without broader theoretical coherence.

No PHQ-9 item will be deleted during Stage 3.

Exploratory results will not be described as confirmatory evidence.

### 12.8 Development-stage model-freezing decision

The one-factor CFA is compulsory and remains the primary validation model regardless of the exploratory results.

At most one alternative multifactor CFA may be frozen before validation results are inspected.

The permitted alternative is:

- one correlated two-factor CFA; or
- one correlated three-factor CFA.

Both alternatives will not be carried forward as competing post hoc models.

A multifactor candidate may be frozen only when all of the following conditions are met:

- the factor count is supported by the development evidence;
- the exploratory loading pattern is stable and substantively interpretable;
- every factor has at least three coherent indicators;
- the structure is psychologically meaningful;
- the structure is not primarily driven by one unusual item;
- the structure is not primarily driven by sparse categories;
- the structure is not merely a wording-overlap artefact;
- the proposed factors are allowed to correlate;
- no CFA cross-loadings are introduced;
- no correlated residuals are introduced;
- the exact item-to-factor allocation is written into the analysis record before validation results are accessed.

The model-freezing record will state:

- whether an alternative model was retained;
- the exploratory evidence supporting the decision;
- the exact item allocation;
- the theoretical interpretation of each factor;
- the reasons the candidate met or failed the freezing criteria;
- the date on which the decision was recorded.

When no multifactor solution satisfies these conditions, only the prespecified one-factor CFA will be fitted in the validation sample.

#### Recorded development-stage decision — 24 July 2026

Ordinal parallel analysis recommended three factors in the development sample. The custom leading-consecutive calculation agreed with `psych::fa.parallel()`. However, the second and third retained factor eigenvalues were small in absolute magnitude, at `.273` and `.111`, respectively, and the exploratory multifactor solutions did not satisfy the prespecified criteria for confirmatory freezing.

The three-factor solution was not eligible for freezing because one factor was primarily defined by only two clear indicators, `DPQ070` and `DPQ080`. This failed the requirement that every frozen factor contain at least three substantively coherent indicators. The three-factor correlations were also high, ranging from `.693` to `.807`, indicating limited separation between the proposed dimensions.

The two-factor solution improved the residual and approximate-fit statistics relative to the one-factor EFA, but its loading pattern was not sufficiently stable or distinct for a simple-structure CFA. `DPQ010`, `DPQ070` and `DPQ080` each had cross-loadings above `.30`, with primary-versus-secondary absolute loading gaps of only `.075`, `.099` and `.115`, respectively. The two factors correlated at `.780`. Freezing this model would therefore require several substantial exploratory cross-loadings to be fixed to zero despite weak empirical separation and insufficiently distinct substantive interpretation.

The exploratory evidence is interpreted as showing a strong general PHQ-9 factor accompanied by reproducible residual multidimensionality that does not organise into a sufficiently coherent simple structure for separate confirmatory validation.

No secondary multifactor CFA was frozen.

The validation stage will contain only the prespecified one-factor ordinal CFA. This decision was recorded before the validation-sample factor structure or CFA results were accessed. Validation results will not be used to retrospectively freeze another model, add correlated residuals or cross-loadings, delete items, or introduce bifactor or other post hoc models.

No bifactor model will be fitted in version 1.

### 12.9 Primary validation-sample CFA

The prespecified primary validation model is:

```r
one_factor_model <- '
  Depression =~ DPQ010 + DPQ020 + DPQ030 +
                DPQ040 + DPQ050 + DPQ060 +
                DPQ070 + DPQ080 + DPQ090
'
```

The model will be estimated using:

```r
cfa_fit <- lavaan::cfa(
  model = one_factor_model,
  data = validation_data,
  ordered = phq9_items,
  estimator = "WLSMV",
  std.lv = TRUE
)
```

The nine PHQ-9 indicators will be declared ordered.

WLSMV is selected for the four-category ordinal indicators. Parameter estimation uses diagonally weighted least squares, with robust standard errors and a mean- and variance-adjusted test statistic.

The default categorical-data parameterisation used by the installed `lavaan` version will be documented rather than changed without a methodological reason.

The R version and installed `lavaan` version will be recorded.

### 12.10 Role of any frozen alternative CFA

No alternative CFA was frozen in the development stage. The one-factor CFA is therefore the sole validation model. The remainder of this subsection documents the prespecified role that would have applied had an eligible alternative been frozen.

The one-factor CFA is the primary validation model.

Any frozen multifactor CFA will be treated as a secondary explanatory model intended to clarify coherent departures from the one-factor representation.

The secondary model will be evaluated using the same categories of evidence as the primary model:

- global-fit indices;
- fully standardised loadings;
- 95% Wald-type confidence intervals for standardised loadings;
- factor correlations;
- threshold diagnostics;
- residual-correlation diagnostics;
- convergence and computational checks;
- substantive interpretability.

A multifactor model will not replace the primary model solely because it produces numerically better fit.

Interpretation of the alternative model will depend on whether any improvement is:

- substantial rather than trivial;
- consistent across global and local evidence;
- supported by a coherent loading structure;
- supported by interpretable factor correlations;
- consistent with the frozen development-sample rationale.

The one-factor compatibility conclusion will be based primarily on the one-factor model itself.

The alternative model will be used to explain departures from the primary model rather than retroactively redefine the primary research question.

No naïve subtraction of WLSMV chi-square statistics will be used as the primary model-selection rule.

No formal chi-square difference test is required for version 1. Any later use of a formal robust difference-testing procedure would require a documented methodological amendment and an implementation appropriate to the fitted estimator.

### 12.11 CFA fit-measure extraction

Fit measures will be extracted using:

```r
cfa_fit_measures <- lavaan::fitMeasures(cfa_fit)
```

The complete named `fitMeasures()` vector will be retained locally for auditability.

The exported fit table will report, where available:

- mean- and variance-adjusted or scaled model chi-square statistic;
- model degrees of freedom;
- model-test p-value;
- robust CFI;
- robust TLI;
- robust RMSEA;
- robust RMSEA 90% confidence interval;
- SRMR.

Any prespecified robust measure that is unavailable or returned as `NA` will be reported and investigated rather than silently replaced with a different statistic.

Global fit will be interpreted as a combined pattern.

No single fit index or conventional cutoff will independently determine whether the model passes or fails.

### 12.12 Fully standardised loadings

The fully standardised solution will be extracted using:

```r
cfa_standardized <- lavaan::standardizedSolution(
  cfa_fit,
  type = "std.all",
  ci = TRUE,
  level = 0.95
)
```

Rows with `op == "=~"` will be used to report:

- fully standardised factor loadings;
- standard errors;
- 95% Wald-type confidence intervals.

These intervals will not be described as bootstrap confidence intervals.

Loading magnitude, direction, uncertainty and consistency across items will be considered together.

No item will be removed solely because of a comparatively low loading.

### 12.13 Threshold extraction

Parameter estimates will be extracted using:

```r
cfa_parameters <- lavaan::parameterEstimates(
  cfa_fit,
  ci = TRUE,
  level = 0.95
)

cfa_thresholds <- cfa_parameters |>
  dplyr::filter(op == "|")
```

Thresholds will be checked within each item for:

- finite estimates;
- expected number of thresholds;
- monotonically increasing order;
- unexplained extreme values;
- missing or duplicated estimates.

An invalid threshold pattern or unexplained estimation problem will stop interpretation pending review.

### 12.14 Residual and local-fit diagnostics

Residual diagnostics will be extracted using:

```r
cfa_residuals <- lavaan::lavResiduals(
  cfa_fit,
  type = "cor.bentler",
  zstat = TRUE,
  summary = TRUE
)
```

The analysis will distinguish between:

- residual-correlation values;
- z-standardised residual statistics;
- summary measures such as SRMR.

The absolute `.10` flag applies to residual-correlation values, not to z-standardised residual statistics.

An absolute residual correlation of `.10` or greater will be treated as:

- a descriptive flag;
- evidence of a potentially notable local area of misfit;
- a reason to inspect the substantive content of the item pair;
- not a model pass/fail threshold;
- not an automatic reason to add a correlated residual;
- not an automatic reason to delete an item.

The largest residual-correlation pairs will be reported.

The analysis will examine whether flagged residuals:

- are isolated;
- cluster around a particular item;
- form a coherent content-based pattern;
- correspond to the frozen exploratory interpretation;
- suggest a meaningful limitation of the one-factor representation.

Modification indices will not be used as an automatic model-building procedure.

No post hoc correlated residuals will be added in version 1.

### 12.15 Computational and improper-solution checks

For every validation model, the following will be documented:

- convergence status;
- estimator;
- parameterisation;
- number of observations;
- number of free parameters;
- model-test information;
- warnings;
- successful post-estimation checks;
- finite parameter estimates;
- threshold ordering;
- standardised loading direction and magnitude;
- residual variances;
- latent-factor variances;
- latent-factor correlations where applicable;
- presence or absence of improper estimates.

Potential improper solutions include:

- non-finite estimates;
- negative residual variances;
- standardised estimates outside an interpretable range;
- factor correlations at or above an absolute value of 1;
- non-positive-definite latent covariance matrices;
- unexplained computational warnings.

An improper or non-converged model will not be substantively interpreted until the problem has been reviewed.

### 12.16 Dimensionality interpretation categories

The final Stage 3 conclusion will use one of three formulations.

#### Compatible with a one-factor representation

This conclusion may be used when:

- global and local evidence is broadly supportive;
- standardised loadings are substantively meaningful;
- there is no serious coherent pattern of multidimensionality;
- residual misfit is limited or substantively minor;
- any alternative model does not reveal a clearly superior and coherent structure that materially changes interpretation.

#### Mixed evidence

This conclusion may be used when:

- a strong general factor is evident;
- most loadings are meaningful;
- one or more global-fit indices indicate notable departure;
- local residual patterns indicate some coherent misfit;
- exploratory or secondary-model evidence suggests limited multidimensional structure;
- the evidence does not justify a simple claim of either adequate unidimensionality or clear multidimensionality.

#### Not adequately represented by one factor

This conclusion may be used when:

- several forms of evidence converge on substantial misfit;
- local residuals form a coherent pattern;
- exploratory evidence supports an interpretable multifactor structure;
- a frozen alternative model provides substantially and consistently better explanation;
- the departure is not merely statistical sensitivity caused by the large sample.

Even a favourable result will not be described as proving strict unidimensionality.

### 12.17 Stage 3 interpretation boundaries

Stage 3 will distinguish:

- internal consistency from dimensionality;
- dimensionality from validity;
- a strong general factor from strict unidimensionality;
- sample-based psychometric estimates from population-weighted estimates;
- symptom measurement from diagnosis;
- statistical model fit from clinical utility;
- exploratory findings from validation evidence.

Stage 3 will not make claims about:

- diagnostic accuracy;
- sensitivity or specificity;
- treatment need;
- suicide-risk assessment;
- causal mechanisms;
- longitudinal change;
- clinical decision-making.

### 12.18 Planned Stage 3 outputs

Local processed objects will include:

```text
data/processed/phq9_stage3_split.rds
data/processed/phq9_parallel_analysis_object.rds
```

These processed objects will remain excluded from the public repository unless a later data-publication decision is made.

Planned public tables include:

```text
tables/phq9_stage3_split_summary.csv
tables/phq9_parallel_analysis_results.csv
tables/phq9_efa_loadings.csv
tables/phq9_efa_factor_correlations.csv
tables/phq9_efa_diagnostics.csv
tables/phq9_stage3_model_decision.csv
tables/phq9_cfa_fit_measures.csv
tables/phq9_cfa_standardized_loadings.csv
tables/phq9_cfa_thresholds.csv
tables/phq9_cfa_residual_correlations.csv
tables/phq9_cfa_computational_diagnostics.csv
```

Planned figures include:

```text
figures/phq9_parallel_analysis.png
figures/phq9_efa_loadings.png
figures/phq9_cfa_standardized_loadings.png
```

Only necessary, interpretable outputs will be included in the Quarto report.

The session-information record will be updated after Stage 3 execution.

### 12.19 Stage 3 execution order

Stage 3 will proceed in the following order:

1. commit analysis plan version 1.5;
2. implement the fixed development and validation split;
3. validate both split samples;
4. estimate and validate the development polychoric matrix;
5. run the locked ordinal parallel analysis;
6. extract and validate the 95th-percentile comparison;
7. fit the permitted development-sample EFA solutions;
8. record the alternative-model freezing decision;
9. access and fit the primary validation-sample one-factor CFA;
10. fit one frozen secondary CFA only when permitted by the development-stage decision;
11. extract global-fit, loading, threshold, residual and computational diagnostics;
12. export tables and figures;
13. update the Quarto report;
14. complete Stage 3 methodological, academic-voice and repository review.

At step 8, no alternative multifactor model satisfied the freezing criteria. Step 10 is therefore not applicable, and the validation stage will contain only the prespecified one-factor ordinal CFA.

No parallel analysis, EFA, CFA or model comparison will be run before analysis plan version 1.5 is committed.

## 13. Later planned analyses

### 13.1 Functional-difficulty association

The relationship between PHQ-9 symptom severity and `DPQ100` will be examined after Stage 3 using:

- weighted PHQ-9 estimates across functional-difficulty categories;
- visualisation of the score–difficulty relationship;
- survey-weighted ordinal logistic regression with `DPQ100` as the ordered outcome and PHQ-9 total score as the main predictor;
- assessment of the proportional-odds assumption;
- effect estimates with confidence intervals.

If the proportional-odds assumption is not supported, a less restrictive survey-weighted model will be selected and documented before inferential results are interpreted.

Results will be interpreted as evidence about the relationship between symptom scores and reported functional impairment, rather than evidence of diagnostic accuracy.

### 13.2 Limitations and interpretation

The discussion will address:

- the cross-sectional design;
- self-reported symptoms;
- missingness and complete-case selection;
- the absence of an independent diagnostic interview;
- the use of a population rather than clinical sample;
- limits to generalisation outside the United States;
- the distinction between symptom severity and psychiatric diagnosis;
- limitations of conventional severity thresholds;
- interpretation of PHQ-9 item 9;
- the distinction between statistical associations and clinical decisions.

## 14. Weighted analyses and psychometric modelling

### Weighted analyses

Population-description and functional-association analyses use the NHANES complex survey design:

- weight: `WTMECPRP`
- strata: `SDMVSTRA`
- cluster: `SDMVPSU`

These analyses are implemented using the `survey` and/or `srvyr` packages.

Weighted estimates describe the U.S. civilian non-institutionalised adult population represented by the combined 2017–March 2020 pre-pandemic release.

### Psychometric analyses

Reliability, polychoric correlations and factor analyses use the complete-PHQ-9 analytic sample.

Unless survey weights are incorporated into a specific psychometric model, resulting parameters will be described as estimates for the analytic sample rather than nationally representative latent-variable parameters.

The report will distinguish:

- survey-weighted population descriptions;
- sample-based psychometric estimates.

Survey-weighted latent-variable modelling is outside the version 1 scope.

## 15. Clinical interpretation

PHQ-9 scores will be interpreted as **self-reported indicators of recent depressive symptom severity**. They will not be presented as psychiatric diagnoses.

The dataset does not establish whether a participant met diagnostic criteria for a depressive disorder, required treatment or would benefit from a particular intervention.

`DPQ090`, corresponding to PHQ-9 item 9, concerns thoughts of death or self-harm. Endorsement is clinically important, but item 9 is not a standalone suicide-risk assessment. It does not establish intent, planning, imminence or overall suicide risk.

## 16. Version 1 completion criteria

Version 1 will contain:

1. a documented data and codebook audit;
2. a transparent sample-flow table;
3. weighted item-response descriptives;
4. weighted PHQ-9 total-score and severity-band descriptives;
5. item floor-effect and category-sparsity assessment;
6. corrected item-total correlations;
7. ordinal coefficient alpha;
8. McDonald’s omega total;
9. a polychoric correlation matrix;
10. a development and validation dimensionality analysis centred on the prespecified one-factor ordinal CFA;
11. an analysis of PHQ-9 scores in relation to functional difficulty;
12. a discussion of interpretation and limitations.

## 17. Not included in version 1

Version 1 will not include:

- diagnostic-accuracy analysis;
- sensitivity or specificity estimates;
- claims about whether participants need treatment;
- claims that PHQ-9 scores establish a depressive disorder;
- measurement-invariance testing;
- survey-weighted latent-variable modelling;
- bifactor modelling;
- post hoc item deletion;
- post hoc correlated residuals;
- unrestricted exploratory model searching;
- separate analysis of the partial 2019–March 2020 cycle.

These exclusions reflect the information available in the dataset and the defined scope of the project.

## 18. Optional extensions

Potential extensions after version 1 is complete include:

- PHQ-8 sensitivity analysis excluding item 9;
- subgroup descriptives;
- measurement invariance across sex or broad age groups;
- secondary associations with general health;
- secondary associations with mental-health professional contact;
- exploratory sleep-related analyses;
- survey-weighted reliability estimates;
- survey-weighted latent-variable modelling;
- bifactor analysis only if later theoretically and methodologically justified;
- replication using another NHANES cycle.

Measurement invariance and bifactor modelling are extensions and are not required for version 1.

## 19. Repository structure

```text
phq9-nhanes-psychometrics/
├── README.md
├── ANALYSIS_PLAN.md
├── phq9-nhanes-psychometrics.Rproj
├── data/
│   ├── raw/
│   └── processed/
├── scripts/
│   ├── 00_download_data.R
│   ├── 01_clean_data.R
│   ├── 02_descriptives.R
│   ├── 03_reliability.R
│   ├── 04_factor_analysis.R
│   └── 05_figures_tables.R
├── report/
│   ├── report.qmd
│   └── report.html
├── figures/
├── tables/
├── references.bib
├── session-info.txt
├── .gitignore
└── LICENSE
```

## 20. Execution checklist

### Stage 1: data audit

- [x] Create the GitHub repository and local RStudio project.
- [x] Add the documented directory structure.
- [x] Save the analysis plan.
- [x] Create the initial `README.md`.
- [x] Create the six numbered R scripts.
- [x] Create and render the initial `report/report.qmd`.
- [x] Download `P_DPQ.xpt` and `P_DEMO.xpt` programmatically.
- [x] Import and inspect both files.
- [x] Verify file linkage using `SEQN`.
- [x] Create a PHQ-9 variable dictionary.
- [x] Audit all item response codes.
- [x] Audit missingness and the `DPQ100` response pattern.
- [x] Construct and document the eligible analytic samples.
- [x] Audit the NHANES survey-design variables.
- [x] Produce the initial sample-flow and missingness outputs.
- [x] Complete Stage 1 review and validation.

### Stage 2A: weighted descriptives

- [x] Create the NHANES survey-design object.
- [x] Calculate weighted PHQ-9 item-response proportions.
- [x] Describe the weighted PHQ-9 total-score distribution.
- [x] Calculate the weighted symptom-severity band distribution.
- [x] Describe `DPQ100` among routed respondents.
- [x] Create the extended functional-difficulty category including structural skips.
- [x] Export descriptive tables and figures.
- [x] Update the Quarto report.
- [x] Complete the academic-voice and formatting check.
- [x] Complete Stage 2A review and validation.

### Stage 2B: reliability and item analysis

- [x] Assess item floor effects and response-category sparsity.
- [x] Calculate corrected item-total correlations.
- [x] Estimate the polychoric correlation matrix.
- [x] Estimate ordinal coefficient alpha.
- [x] Estimate McDonald’s omega.
- [x] Examine reliability estimates after removing individual items where methodologically informative.
- [x] Calculate uncertainty intervals where feasible.
- [x] Export Stage 2B tables and any necessary figures.
- [x] Update the Quarto report.
- [x] Complete the academic-voice, methodological-language and formatting check.
- [x] Complete Stage 2B review and validation.

### Stage 3: dimensionality analysis

- [x] Lock the Stage 3 dimensionality-analysis specification before analysis.
- [x] Create the fixed-seed development and validation split.
- [x] Validate sample size, response coding and category availability in both split samples.
- [x] Estimate and validate the unsmoothed development-sample polychoric matrix.
- [x] Conduct the controlled single-core ordinal parallel analysis.
- [x] Reconstruct the resampled 95th-percentile factor-eigenvalue thresholds.
- [x] Check the custom leading-consecutive recommendation against `parallel_results$nfact`.
- [x] Conduct the permitted development-sample exploratory factor analyses.
- [x] Record the alternative-model freezing decision before accessing validation results; no eligible alternative CFA was frozen.
- [x] Fit the prespecified one-factor ordinal CFA in the validation sample.
- [x] Determine that no secondary multifactor CFA will be fitted because no candidate met the freezing criteria.
- [x] Extract and validate the prespecified robust global-fit measures.
- [x] Report fully standardised factor loadings with 95% Wald-type confidence intervals.
- [x] Extract and validate all item thresholds.
- [x] Document the estimator, parameterisation, model-test correction, convergence and improper-solution diagnostics.
- [x] Examine residual correlations using the descriptive absolute `.10` flag.
- [x] Confirm that all standardised residual variances are finite, positive and within admissible bounds.
- [x] Record the validation-stage dimensionality conclusion as mixed evidence.
- [x] Avoid post hoc correlated residuals, cross-loadings, item deletion, bifactor modelling and retrospective alternative-model selection.
- [ ] Export Stage 3 tables and figures.
- [ ] Update the session-information record.
- [ ] Update the Quarto report.
- [ ] Complete the academic-voice, methodological-language and formatting check.
- [ ] Complete Stage 3 review and validation.

Stage 2A and Stage 2B review and validation are complete.

The Stage 3 dimensionality-analysis specification was locked before analysis on 23 July 2026. The fixed split, split validation, development polychoric matrix, ordinal parallel analysis, permitted EFAs, model-freezing decision and prespecified validation-sample one-factor ordinal CFA are complete. No secondary multifactor CFA was frozen or fitted.

The validation CFA supported a substantial general PHQ-9 factor, with fully standardised loadings from `.707` to `.873`. However, global fit showed notable departure from the one-factor representation: robust CFI was `.924`, robust TLI was `.898`, robust RMSEA was `.127`, with a 90% confidence interval from `.114` to `.140`, and SRMR was `.048`. One residual correlation met the prespecified absolute `.10` descriptive flag. The Stage 3 dimensionality conclusion was therefore recorded as **mixed evidence**, rather than as establishing strict unidimensionality.

Stage 3 table and figure export, session-information updating, Quarto integration, academic-voice review and final repository validation remain outstanding.