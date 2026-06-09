source("scripts/config.R")
source("scripts/preprocessing.R")
source("scripts/pmdd_diagnosis.R")
source("scripts/visualization.R")
source("scripts/menstrualcycleR_functions.R")
source("scripts/pmdd_scoring.R")

# ============================================================================
# STEP 2: RUN PMDD ALGORITHM ON SCALED DATA
# ============================================================================

pmdd_results <- data_combined_filtered %>%
  group_by(participant_id) %>%
  summarise(
    pmdd_result = list(check_pmdd(pick(everything()))),
    n_days_analyzed = n(),
    .groups = "drop"
  ) %>%
  mutate(
    algorithmic_pmdd     = map_int(pmdd_result, ~ .x$pmdd_diagnosis %||% NA_integer_),
    valid_cycles         = map_int(pmdd_result, ~ .x$valid_cycles),
    pmdd_positive_cycles = map_int(pmdd_result, ~ .x$pmdd_positive_cycles)
  ) %>%
  select(-pmdd_result) %>%
  left_join(
    profile_combined %>%
      select(participant_id, age, PMDD, MDD, GAD, PTSD, BD, ADHD, BPD, med, diag),
    by = "participant_id"
  ) %>%
  rename(clinical_pmdd = PMDD)

results_one_cycle  <- pmdd_results %>% filter(valid_cycles == 1)
results_two_cycles <- pmdd_results %>% filter(valid_cycles >= 2)

ids_one_cycle <- results_one_cycle %>% pull(participant_id)
ids_two_cycle <- results_two_cycles %>% pull(participant_id)

data_one_cycle <- data_combined_filtered %>% 
  filter(participant_id %in% ids_one_cycle) %>% 
  mutate(id = participant_id) %>% 
  mutate(menses = if_else(bleeding %in% c(1,2), 1, 0)) %>% 
  mutate(ovtoday = 0) %>% 
  select(-depression, -mania, -comments, -hours)
data_two_cycle <- data_combined_filtered %>% 
  filter(participant_id %in% ids_two_cycle) %>% 
  mutate(id = participant_id) %>% 
  mutate(menses = if_else(bleeding %in% c(1,2), 1, 0)) %>% 
  mutate(ovtoday = 0) %>% 
  select(-depression, -mania, -comments, -hours)

# ============================================================================
# PACTS and GAMM
# ============================================================================
# 0. Setup: Apply PACTS scaling to get your final dataset
data_two_pacts <- pacts_scaling(data = data_two_cycle, id = id, date = date, menses = menses, ovtoday = ovtoday) %>% 
  mutate(bleeding = replace_na(bleeding, 0))

# 1. Run the batch function 
batch_results <- fit_and_predict_all_symptoms(data_two_pacts, symptoms_list)

# 2. Extract the pieces you need from the results
my_fitted_models    <- batch_results$models
gamm_predictions_df <- batch_results$predictions_df

# 3. Calculate your PMDD scores!
final_pmdd_scores <- compute_all_scores(data_two_pacts, gamm_predictions_df)

test <- results_two_cycles %>% 
  left_join(final_pmdd_scores, by = "participant_id")
# ============================================================================
# OTHER
# ============================================================================
plot_score_distributions(test)

test <- data_combined_filtered %>% 
  filter(participant_id == "4gHpF0976Jhpc8gfycO0HqwJSom1")

create_pmdd_calendar(test, "4gHpF0976Jhpc8gfycO0HqwJSom1")

test <- data_two_cycle %>% 
  group_by(participant_id) %>% 
  summarise(data = n()) %>%    # Use n() to count rows
  arrange(desc(data))          # Tell desc() which column to sort

# GAMM MODEL first fits model around baseline popn INCLUDING NON PMDD PEOPLE. Then it tries to skew the results per participant.