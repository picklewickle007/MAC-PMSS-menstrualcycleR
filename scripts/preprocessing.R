source("scripts/config.R")

# ============================================================================
# PROFILE PROCESSING
# ============================================================================

profile_android <- profile_android_ori %>% mutate(source = "android")
profile_ios <- profile_ios_ori %>% mutate(source = "ios")
profile_combined_ori <- full_join(profile_android, profile_ios)

combine_unique <- function(x) {
  all_lists <- map(x, ~{
    parsed <- str_replace_all(.x, "\\[|\\]|'", "") %>%
      str_split(",") %>%
      unlist() %>%
      str_trim()
    parsed[parsed != ""]
  })
  
  all_items <- unlist(all_lists)
  unique_items <- unique(all_items[all_items != "" & !is.na(all_items)])
  
  if (length(unique_items) == 0) {
    return("['']")
  }
  paste0("['", paste(unique_items, collapse = "', '"), "']")
}

profile_combined <- profile_combined_ori %>% 
  select(-query_done, -date) %>% 
  group_by(participant_id) %>% 
  summarise(
    age = last(age),
    medication = combine_unique(medication),
    diagnosis = combine_unique(diagnosis),
    source = last(source),
    .groups = 'drop'
  ) %>%  
  # Diagnostic flags
  mutate(
    MDD  = ifelse(str_detect(diagnosis, regex("depression", ignore_case = TRUE)), 1, 0),
    PMDD = ifelse(str_detect(diagnosis, regex("dysphoric", ignore_case = TRUE)), 1, 0),
    GAD  = ifelse(str_detect(diagnosis, regex("anxiety", ignore_case = TRUE)), 1, 0),
    PTSD = ifelse(str_detect(diagnosis, regex("PTSD", ignore_case = TRUE)), 1, 0),
    BD   = ifelse(str_detect(diagnosis, regex("bipolar", ignore_case = TRUE)), 1, 0),
    ADHD = ifelse(str_detect(diagnosis, regex("ADHD", ignore_case = TRUE)), 1, 0),
    BPD  = ifelse(str_detect(diagnosis, regex("borderline", ignore_case = TRUE)), 1, 0)
  ) %>%
  # Medication & diagnosis summary flags
  mutate(
    med = ifelse(
      medication == "['']" |
        medication == "['', '']" |
        str_detect(medication, regex("none", ignore_case = TRUE)),
      0, 1
    ),
    diag = ifelse(
      diagnosis == "['']" |
        diagnosis == "['No previous diagnosis']",
      0, 1
    )
  )

# ============================================================================
# DATA PROCESSING
# ============================================================================

data_android <- data_android_ori %>% mutate(source = "android")
data_ios <- data_ios_ori %>% mutate(source = "ios")

data_combined_ori <- full_join(data_android, data_ios)

data_combined <- data_combined_ori %>%
  mutate(across(all_of(symptoms_list), ~ . + 1)) %>% # add 1 to all symptom values
  mutate(total = rowSums(across(all_of(symptoms_list)))) %>%
  mutate(comments = replace_na(comments, "")) %>% 
  filter(if_all(all_of(symptoms_list), ~ . >= 1)) %>% # filter any symptom severities less than 1
  mutate(date = as.Date(date))

data_combined_filtered <- data_combined %>% #keep streaks of contin day greater than minimum amount 
  arrange(participant_id, date) %>% 
  group_by(participant_id) %>% 
  mutate(
    day_diff = as.numeric(date - lag(date)),
    new_block = is.na(day_diff) | day_diff > 1,
    block_id = cumsum(new_block)
  ) %>% 
  group_by(participant_id, block_id) %>% 
  mutate(block_length = n()) %>% 
  ungroup() %>% 
  filter(block_length >= MIN_CONTINUOUS_DAYS) %>% 
  select(-day_diff, -new_block, -block_id, -block_length) 
