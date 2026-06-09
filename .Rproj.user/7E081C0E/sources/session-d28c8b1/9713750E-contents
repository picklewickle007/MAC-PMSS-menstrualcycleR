source("scripts/config.R")

# ============================================================================
# PMDD ALGORITHM
# ============================================================================

test <- data_combined_filtered %>% 
  filter(participant_id == "Cp5m74wnhtT3UTYEp3wxzKOCwCY2")

find_bleeding_start_days <- function(df) {
  
  min_gap <- 10
  period_starts <- c()
  current_sequence <- c()
  current_sequence_index <- NA
  last_period_start <- -Inf
  in_sequence <- FALSE
  

  for(i in 1:nrow(df)) {
    bleeding_val <- df$bleeding[i]
    
    if(!in_sequence && (bleeding_val == 1 | bleeding_val == 2)) {
      in_sequence <- TRUE
      current_sequence_index <- i
      current_sequence <- c(bleeding_val)
    } else if (in_sequence && (bleeding_val == 1 | bleeding_val == 2)) {
      current_sequence <- c(current_sequence, bleeding_val)
    } else if (in_sequence) {
  
      if(length(current_sequence) >= 2 && any(current_sequence == 2) && (current_sequence_index - last_period_start >= min_gap)) {
        period_starts <- c(period_starts, current_sequence_index)
        last_period_start <- current_sequence_index
      }
      in_sequence <- FALSE
      current_sequence <- c()
    }
  }
  
  if(in_sequence && length(current_sequence) >= 2 && any(current_sequence == 2) && (current_sequence_index - last_period_start >= min_gap)) {
    period_starts <- c(period_starts, current_sequence_index)
  }
  
  return(period_starts)
}

diagnosis_PMDD <- function(df_before, df_after) {
  # Dynamically pick 'max' or 'min' based on your configuration
  f_after <- match.fun(FOLLICULAR_METHOD)
  
  # Calculate metrics for all symptoms at once
  symptom_stats <- map_dfr(symptoms_list, ~ {
    b_vals <- df_before[[.x]]
    a_vals <- df_after[[.x]]
    
    list(
      max_before = max(b_vals, na.rm = TRUE),
      severity_after = f_after(a_vals, na.rm = TRUE),
      is_key = .x %in% c("lability", "anger", "depressed", "anxiety")
    )
  })
  
  # Apply PMDD criteria using vector math
  criteria <- symptom_stats %>%
    filter(max_before >= 4) %>%
    mutate(
      # The core PMDD rule: "lack of improvement"
      # If extreme (6), it must drop below 3; otherwise, below 4.
      no_improvement = if_else(max_before == 6, 
                               severity_after >= 3, 
                               severity_after >= 4)
    )
  
  # Return 1 if criteria met, else 0
  as.integer(
    nrow(criteria) >= 5 &&               # At least 5 severe symptoms present
      any(criteria$is_key) &&            # At least one key symptom
      sum(!criteria$no_improvement) >= 5 # CHANGED: At least 5 symptoms DID improve (dropped to baseline)
  )
}

check_pmdd <- function(df) {
  starts <- find_bleeding_start_days(df)
  
  # Generate all valid cycle segments at once
  cycles <- map(starts, ~ {
    if (.x - 7 < 1 || .x + 8 > nrow(df)) return(NULL)
    list(before = df[(.x-7):(.x-1), ], after = df[(.x+2):(.x+8), ])
  }) %>% compact()
  
  # Evaluate cycles
  results <- map_lgl(cycles, ~ diagnosis_PMDD(.x$before, .x$after))
  
  list(
    pmdd_diagnosis       = as.integer(any(results)),
    valid_cycles         = length(results),
    pmdd_positive_cycles = sum(results)
  )
}

