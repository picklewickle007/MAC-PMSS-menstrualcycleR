source("scripts/config.R")

# ============================================================================
# CUSTOM FUNCTIONS
# ============================================================================
# Calculates Area Under the Curve using the trapezoidal rule (replaces MESS::auc)
calculate_auc <- function(x, y) {
  # Remove NAs to prevent calculation errors
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]
  
  # Ensure x is in sequential order
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  
  # If we don't have at least 2 points, we can't calculate an area
  if (length(x) < 2) return(NA_real_)
  
  # Trapezoidal calculation
  dx <- diff(x)
  mean_y <- (y[-1] + y[-length(y)]) / 2
  
  return(sum(dx * mean_y))
}

# ============================================================================
# PMDD INDEX SCORING
# ============================================================================
# All three scores use the same windows as the PMDD algorithm:
#   - Luteal window: 7 days before period start (idx-7 to idx-1)
#   - Follicular window: days +2 to +8 after period start (idx+2 to idx+8)
#
# Score 1: Mean Luteal Severity
#   - Average symptom severity across the luteal window, averaged across cycles
#
# Score 2: Cyclical Severity (Mean Luteal-Follicular Difference)
#   - Average drop in symptom severity from luteal to follicular, across cycles
#
# Score 3: DSM-Weighted Composite
#   - Weighted average of luteal-follicular drop: 60% core, 40% secondary
#   - Only computed for cycles where at least 1 core symptom peaked >= 4
#
# Score 4:
# ============================================================================

CORE_SYMPTOMS      <- c("lability", "anger", "depressed", "anxiety")
SECONDARY_SYMPTOMS <- setdiff(symptoms_list, CORE_SYMPTOMS)

extract_windows <- function(df, idx, need_follicular = TRUE) {
  
  if(idx - 7 < 1) return(NULL)
  luteal_vals <- df[(idx - 7):(idx - 1), symptoms_list, drop = FALSE]
  
  if(!need_follicular) {
    return(list(luteal_vals = luteal_vals))
  }
  
  if(idx +8 > nrow(df)) return(NULL)
  follicular_vals <- df[(idx+2):(idx+8), symptoms_list, drop = FALSE]
  
  return(list(
    luteal_vals = luteal_vals,
    follicular_vals = follicular_vals
  ))
                        
}

# ============================================================================
# Score 1: Mean Luteal Severity
# ============================================================================

compute_score1 <- function(df) {
  
  bleeding_starts <- find_bleeding_start_days(df)
  
  if (length(bleeding_starts) == 0) {
    return(list(score1 = NA_real_, n_cycles_scored = 0L, cycle_scores = tibble()))
  }
  
  # Map directly to a combined Data Frame, automatically dropping NULLs
  cycle_scores_df <- map_dfr(bleeding_starts, ~ {
    windows <- extract_windows(df, .x, need_follicular = FALSE)
    if (is.null(windows)) return(NULL)
    
    L_per_symptom <- colMeans(windows$luteal_vals, na.rm = FALSE)
    cycle_mean_L  <- mean(L_per_symptom)
    
    # Safely construct the row
    bind_cols(
      tibble(cycle_start_idx = .x, cycle_mean_L = cycle_mean_L),
      as_tibble_row(L_per_symptom, .name_repair = ~ paste0("L_", .x))
    )
  })
  
  # Catch the case where windows were extracted but all returned NULL due to NAs
  if (nrow(cycle_scores_df) == 0) {
    return(list(score1 = NA_real_, n_cycles_scored = 0L, cycle_scores = tibble()))
  }
  
  score1 <- mean(cycle_scores_df$cycle_mean_L, na.rm = TRUE)
  
  return(list(
    score1          = round(score1, 3),
    n_cycles_scored = nrow(cycle_scores_df),
    cycle_scores    = cycle_scores_df
  ))
}

# ============================================================================
# Score 2: Cyclical Severity
# ============================================================================

compute_score2 <- function(df) {
  bleeding_starts <- find_bleeding_start_days(df)
  
  if (length(bleeding_starts) == 0) {
    return(list(score2 = NA_real_, cycle_deltas = tibble()))
  }
  
  cycle_deltas_df <- map_dfr(bleeding_starts, ~ {
    windows <- extract_windows(df, .x, need_follicular = TRUE)
    if (is.null(windows)) return(NULL)
    
    L_per_symptom <- colMeans(windows$luteal_vals,     na.rm = FALSE)
    F_per_symptom <- colMeans(windows$follicular_vals, na.rm = FALSE)
    
    delta_per_symptom <- L_per_symptom - F_per_symptom
    cycle_avg_delta   <- mean(delta_per_symptom)
    
    bind_cols(
      tibble(cycle_start_idx = .x, cycle_avg_delta = cycle_avg_delta),
      as_tibble_row(delta_per_symptom, .name_repair = ~ paste0("delta_", .x))
    )
  })
  
  if (nrow(cycle_deltas_df) == 0) {
    return(list(score2 = NA_real_, cycle_deltas = tibble()))
  }
  
  score2 <- mean(cycle_deltas_df$cycle_avg_delta, na.rm = TRUE)
  
  return(list(
    score2       = round(score2, 3),
    cycle_deltas = cycle_deltas_df
  ))
}

# ============================================================================
# Score 3: DSM-Weighted Severity
# ============================================================================

compute_score3 <- function(df) {
  bleeding_starts <- find_bleeding_start_days(df)
  
  if (length(bleeding_starts) == 0) {
    return(list(score3 = NA_real_, n_cycles_threshold_met = 0L, cycle_composites = tibble()))
  }
  
  cycle_composites_df <- map_dfr(bleeding_starts, ~ {
    windows <- extract_windows(df, .x, need_follicular = TRUE)
    if (is.null(windows)) return(NULL)
    
    # Safe subsetting with drop = FALSE to prevent crashes if 1 core symptom
    core_luteal_vals <- windows$luteal_vals[, CORE_SYMPTOMS, drop = FALSE]
    core_max_luteal  <- apply(core_luteal_vals, 2, max, na.rm = TRUE)
    
    # Skip cycle if no core symptom peaked >= 4
    if (!any(core_max_luteal >= 4)) return(NULL)
    
    L_per_symptom <- colMeans(windows$luteal_vals,     na.rm = FALSE)
    F_per_symptom <- colMeans(windows$follicular_vals, na.rm = FALSE)
    
    delta_per_symptom   <- L_per_symptom - F_per_symptom
    avg_core_delta      <- mean(delta_per_symptom[CORE_SYMPTOMS])
    avg_secondary_delta <- mean(delta_per_symptom[SECONDARY_SYMPTOMS])
    
    cycle_composite <- 0.6 * avg_core_delta + 0.4 * avg_secondary_delta
    
    tibble(
      cycle_start_idx     = .x,
      avg_core_delta      = avg_core_delta,
      avg_secondary_delta = avg_secondary_delta,
      cycle_composite     = cycle_composite
    )
  })
  
  if (nrow(cycle_composites_df) == 0) {
    return(list(score3 = NA_real_, n_cycles_threshold_met = 0L, cycle_composites = tibble()))
  }
  
  score3 <- mean(cycle_composites_df$cycle_composite, na.rm = TRUE)
  
  return(list(
    score3                 = round(score3, 3),
    n_cycles_threshold_met = nrow(cycle_composites_df),
    cycle_composites       = cycle_composites_df
  ))
}

# ============================================================================
# Score 4: Composite Perimenstrual Worsening Score (GAMM)
# ============================================================================
# Focuses on worsening around menses
# Daily composite score = mean of all MAC-PMSS symptoms
# Score = Mean severity during perimenstrual window − Mean severity during baseline

compute_score4 <- function(gamm_df) {
  if (is.null(gamm_df) || nrow(gamm_df) == 0) return(list(score4 = NA_real_))
  
  # Peri window: Days -7 to -1 AND Days 4 to 10
  peri_window <- gamm_df %>% filter((cycle_day >= -7 & cycle_day <= -1) | (cycle_day >= 4 & cycle_day <= 10))
  
  # Baseline window: The rest of the cycle (everything NOT in the peri window)
  baseline_window <- gamm_df %>% filter(!((cycle_day >= -7 & cycle_day <= -1) | (cycle_day >= 4 & cycle_day <= 10)))
  
  if (nrow(peri_window) == 0 | nrow(baseline_window) == 0) return(list(score4 = NA_real_))
  
  # Calculate daily composites (mean of all symptoms for that specific day)
  peri_composites <- rowMeans(peri_window[, symptoms_list, drop = FALSE], na.rm = TRUE)
  baseline_composites <- rowMeans(baseline_window[, symptoms_list, drop = FALSE], na.rm = TRUE)
  
  score4 <- mean(peri_composites, na.rm = TRUE) - mean(baseline_composites, na.rm = TRUE)
  
  return(list(score4 = round(score4, 3)))
}

# ============================================================================
# Score 5: AUC-based Cyclic Worsening Score (GAMM)
# ============================================================================
# Measures excess burden. 
# Positive score = symptoms accumulate more strongly during perimenstrual window compared to the rest of the cycle.

compute_score5 <- function(gamm_df) {
  if (is.null(gamm_df) || nrow(gamm_df) == 0) return(list(score5 = NA_real_))
  
  # Peri window: Days -7 to -1 AND Days 4 to 10
  peri_window <- gamm_df %>% filter((cycle_day >= -7 & cycle_day <= -1) | (cycle_day >= 4 & cycle_day <= 10))
  
  # Baseline window: The rest of the cycle (everything NOT in the peri window)
  baseline_window <- gamm_df %>% filter(!((cycle_day >= -7 & cycle_day <= -1) | (cycle_day >= 4 & cycle_day <= 10)))
  
  if (nrow(peri_window) == 0 | nrow(baseline_window) == 0) return(list(score5 = NA_real_))
  
  # Calculate AUC for each symptom
  auc_totals <- sapply(symptoms_list, function(symp) {
    baseline <- mean(baseline_window[[symp]], na.rm = TRUE)
    excess   <- peri_window[[symp]] - baseline
    
    # MESS::auc calculates the mathematical area under the curve
    calculate_auc(peri_window$cycle_day, excess)
  })
  
  score5 <- sum(auc_totals, na.rm = TRUE)
  
  return(list(score5 = round(score5, 3)))
}

# ============================================================================
# Score 6: DSM-based Weighted AUC-based Cyclic Worsening Score (GAMM)
# ============================================================================
# Weighted global score. Core symptoms = 60% total weight, Secondary = 40%

compute_score6 <- function(gamm_df) {
  if (is.null(gamm_df) || nrow(gamm_df) == 0) return(list(score6 = NA_real_))
  
  # Peri window: Days -7 to -1 AND Days 4 to 10
  peri_window <- gamm_df %>% filter((cycle_day >= -7 & cycle_day <= -1) | (cycle_day >= 4 & cycle_day <= 10))
  
  # Baseline window: The rest of the cycle (everything NOT in the peri window)
  baseline_window <- gamm_df %>% filter(!((cycle_day >= -7 & cycle_day <= -1) | (cycle_day >= 4 & cycle_day <= 10)))
  
  if (nrow(peri_window) == 0 | nrow(baseline_window) == 0) return(list(score6 = NA_real_))
  
  # Calculate base AUCs
  auc_vals <- sapply(symptoms_list, function(symp) {
    baseline <- mean(baseline_window[[symp]], na.rm = TRUE)
    excess   <- peri_window[[symp]] - baseline
    calculate_auc(peri_window$cycle_day, excess)
  })
  
  # Dynamically calculate weights based on symptom lists to ensure they always equal 1.0
  # e.g., 0.60 / 4 core symptoms = 0.15 each
  weight_core <- 0.60 / length(CORE_SYMPTOMS)
  weight_sec  <- 0.40 / length(SECONDARY_SYMPTOMS)
  
  weighted_auc <- sum(auc_vals[CORE_SYMPTOMS] * weight_core, na.rm = TRUE) + 
    sum(auc_vals[SECONDARY_SYMPTOMS] * weight_sec, na.rm = TRUE)
  
  return(list(score6 = round(weighted_auc, 3)))
}

# ============================================================================
# Apply all scores across all participants
# ============================================================================

# Note: gamm_predictions_df must contain columns: participant_id, cycle_day, 
# and all symptom columns calculated from predict(gamm_model$gam, type="link")

compute_all_scores <- function(data_filtered, gamm_predictions_df = NULL) {
  data_filtered %>%
    group_by(participant_id) %>%
    group_modify(~ {
      
      # 1. Compute raw data scores (Scores 1-3)
      s1 <- compute_score1(.x)
      s2 <- compute_score2(.x)
      s3 <- compute_score3(.x)
      
      # 2. Extract this specific participant's GAMM predictions (if they exist)
      current_id <- .y$participant_id
      
      if (!is.null(gamm_predictions_df)) {
        participant_gamm <- gamm_predictions_df %>% filter(participant_id == current_id)
      } else {
        participant_gamm <- NULL
      }
      
      # 3. Compute GAMM data scores (Scores 4-6)
      s4 <- compute_score4(participant_gamm)
      s5 <- compute_score5(participant_gamm)
      s6 <- compute_score6(participant_gamm)
      
      # 4. Bind together
      tibble(
        score1 = s1$score1,
        score2 = s2$score2,
        score3 = s3$score3,
        score4 = s4$score4,
        score5 = s5$score5,
        score6 = s6$score6
      )
    }) %>%
    ungroup()
}
