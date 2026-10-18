library(DatabaseConnector)
library(CohortMethod)
library(FeatureExtraction)
library(Eunomia)
library(dplyr)

# 1. Connect to SQLite database
db_path <- file.path("data", "eunomia_cdm.sqlite")
connectionDetails <- createConnectionDetails(
  dbms = "sqlite",
  server = db_path
)

# 2. Exclude treatment & comparator drugs from covariate matrix
covariateSettings <- createDefaultCovariateSettings(
  excludedCovariateConceptIds = c(1118084, 1124300),
  addDescendantsToExclude = TRUE
)

cmDataArgs <- createGetDbCohortMethodDataArgs(
  washoutPeriod = 180,
  firstExposureOnly = TRUE,
  removeDuplicateSubjects = "remove all",
  covariateSettings = covariateSettings
)

# 3. Extract baseline data from OMOP CDM
cohortMethodData <- getDbCohortMethodData(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = "main",
  targetId = 1,          # Celecoxib
  comparatorId = 2,      # Diclofenac
  outcomeIds = 3,        # GI Bleed
  exposureDatabaseSchema = "main",
  exposureTable = "cohort",
  outcomeDatabaseSchema = "main",
  outcomeTable = "cohort",
  getDbCohortMethodDataArgs = cmDataArgs
)

# 4. Define Study Population (On-treatment: Day 1 to Day 30 post-exposure)
studyPopArgs <- createCreateStudyPopulationArgs(
  removeSubjectsWithPriorOutcome = TRUE,
  minDaysAtRisk = 1,
  riskWindowStart = 1,
  startAnchor = "cohort start",
  riskWindowEnd = 30,
  endAnchor = "cohort end"
)

studyPop <- createStudyPopulation(
  cohortMethodData = cohortMethodData,
  outcomeId = 3,
  createStudyPopulationArgs = studyPopArgs
)

# 5. Fit Propensity Score Model (L1 Regularized Logistic Regression)
psArgs <- createCreatePsArgs(
  prior = createPrior("laplace", exclude = c(0), useCrossValidation = TRUE)
)

ps <- createPs(
  cohortMethodData = cohortMethodData,
  population = studyPop,
  createPsArgs = psArgs
)

# Save Propensity Score Distribution Plot
plotPs(ps, fileName = "output/ps_distribution.png")

# 6. Perform 1:1 Matching on Propensity Score
matchArgs <- createMatchOnPsArgs(
  caliper = 0.2,
  caliperScale = "standardized",
  maxRatio = 1
)

matchedPop <- matchOnPs(
  population = ps,
  matchOnPsArgs = matchArgs
)

# 7. Compute Covariate Balance (Standardized Mean Difference)
balance <- computeCovariateBalance(
  population = matchedPop,
  cohortMethodData = cohortMethodData
)

# Save the Covariate Balance Scatter Plot
plotCovariateBalanceScatterPlot(balance, fileName = "output/covariate_balance.png")

# 8. Fit Stratified Cox Proportional Hazards Model
outcomeModelArgs <- createFitOutcomeModelArgs(
  modelType = "cox",
  stratified = TRUE
)

outcomeModel <- fitOutcomeModel(
  population = matchedPop,
  cohortMethodData = cohortMethodData,
  fitOutcomeModelArgs = outcomeModelArgs
)

# Print Treatment Effect (Hazard Ratio)
print(summary(outcomeModel))

# Plot Kaplan-Meier Survival Curves
plotKaplanMeier(
  population = matchedPop,
  includeZero = FALSE,
  fileName = "output/kaplan_meier.png"
)

library(ggplot2)

# ----------------------------------------------------------------------
# 9. Vertical Bar Chart: Patient Counts (Pre-matching vs Post-matching)
# ----------------------------------------------------------------------
# Extract sample sizes from the study population and matched population
pre_target  <- sum(studyPop$treatment == 1)
pre_comp    <- sum(studyPop$treatment == 0)
post_target <- sum(matchedPop$treatment == 1)
post_comp   <- sum(matchedPop$treatment == 0)

count_data <- data.frame(
  Cohort = factor(
    rep(c("Celecoxib (Target)", "Diclofenac (Comparator)"), 2),
    levels = c("Celecoxib (Target)", "Diclofenac (Comparator)")
  ),
  Stage = factor(
    c("Before Matching", "Before Matching", "After Matching (1:1)", "After Matching (1:1)"),
    levels = c("Before Matching", "After Matching (1:1)")
  ),
  Patients = c(pre_target, pre_comp, post_target, post_comp)
)

bar_plot <- ggplot(count_data, aes(x = Cohort, y = Patients, fill = Cohort)) +
  geom_bar(stat = "identity", width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(Patients)), vjust = -0.4, size = 4.2, fontface = "bold") +
  facet_wrap(~ Stage) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Celecoxib (Target)" = "#E15759", "Diclofenac (Comparator)" = "#4E79A7")) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Patient Population Before and After 1:1 PS Matching",
    subtitle = "New-user design with 180-day washout period",
    x = NULL,
    y = "Number of Patients"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank()
  )

# Save to output folder
ggsave("output/patient_counts_barchart.png", plot = bar_plot, width = 7.5, height = 5, dpi = 300)


# ----------------------------------------------------------------------
# 10. Bonus OHDSI Diagram: Attrition Flowchart
# ----------------------------------------------------------------------
# Generates the standard publication attrition box diagram
drawAttritionDiagram(
  studyPop,
  targetLabel = "Celecoxib (Target)",
  comparatorLabel = "Diclofenac (Comparator)",
  fileName = "output/attrition_diagram.png"
)