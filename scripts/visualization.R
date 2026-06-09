source("scripts/config.R")

# ============================================================================
# SYMPTOM TRACKING CALENDAR (heatmap per participant)
# ============================================================================

create_pmdd_calendar <- function(participant_data, participant_id) {
  
  severity_levels <- c("1", "2", "3", "4", "5", "6",
                       "Bleeding 0", "Bleeding 1", "Bleeding 2")
  
  cal_data <- participant_data %>%
    mutate(
      date        = as.Date(date),
      month       = floor_date(date, "month"),
      day         = day(date),
      month_label = format(month, "%B %Y")
    ) %>%
    arrange(month, date) %>%
    mutate(
      month_label = factor(month_label,
                           levels = unique(month_label[order(month)]))
    )
  
  complete_cal_data <- cal_data %>%
    group_by(month, month_label) %>%
    summarise(days_in_this_month = lubridate::days_in_month(first(month)),
              .groups = 'drop') %>%
    rowwise() %>%
    mutate(day = list(1:days_in_this_month)) %>%
    unnest(day) %>%
    select(-days_in_this_month) %>%
    left_join(cal_data, by = c("month", "month_label", "day"))
  
  for (symp in symptoms_list) {
    if (!symp %in% names(complete_cal_data)) complete_cal_data[[symp]] <- NA_real_
  }
  
  symptom_long <- complete_cal_data %>%
    select(day, month, month_label, all_of(symptoms_list)) %>%
    pivot_longer(cols = all_of(symptoms_list),
                 names_to = "symptom", values_to = "severity") %>%
    mutate(
      symptom      = factor(symptom, levels = rev(symptoms_list)),
      severity_cat = case_when(
        is.na(severity) ~ NA_character_,
        severity == 1   ~ "1",
        severity == 2   ~ "2",
        severity == 3   ~ "3",
        severity == 4   ~ "4",
        severity == 5   ~ "5",
        severity == 6   ~ "6",
        TRUE            ~ NA_character_
      )
    )
  
  bleeding_long <- complete_cal_data %>%
    select(day, month, month_label, bleeding) %>%
    mutate(
      symptom      = "Bleeding",
      severity     = bleeding,
      severity_cat = case_when(
        is.na(bleeding) ~ NA_character_,
        bleeding == 2   ~ "Bleeding 2",
        bleeding == 1   ~ "Bleeding 1",
        bleeding == 0   ~ "Bleeding 0",
        TRUE            ~ NA_character_
      )
    ) %>%
    select(day, month_label, symptom, severity, severity_cat)
  
  all_data <- bind_rows(symptom_long, bleeding_long) %>%
    mutate(
      symptom      = factor(symptom, levels = c(rev(symptoms_list), "Bleeding")),
      severity_cat = factor(severity_cat, levels = severity_levels)
    )
  
  color_palette <- c(
    "1"          = "#5a7d5a",
    "2"          = "#c4d68c",
    "3"          = "#f4d679",
    "4"          = "#e8a66a",
    "5"          = "#d9534f",
    "6"          = "#6C3BAA",
    "Bleeding 2" = "#666666",
    "Bleeding 1" = "#cccccc",
    "Bleeding 0" = "#f5f5f5"
  )
  
  p <- ggplot(all_data, aes(x = day, y = symptom)) +
    geom_tile(aes(fill = severity_cat), color = "white", linewidth = 0.5) +
    scale_fill_manual(values   = color_palette,
                      breaks   = severity_levels,
                      na.value = "grey90") +
    scale_x_continuous(breaks = 1:31, expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    facet_wrap(~month_label, ncol = 1) +
    labs(
      title = paste("Symptom Tracking Calendar - Participant", participant_id),
      x     = "Day of Month",
      y     = NULL,
      fill  = "Severity"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      panel.grid      = element_blank(),
      strip.text      = element_text(face = "bold", size = 12, hjust = 0),
      axis.text.y     = element_text(size = 10),
      axis.text.x     = element_text(size = 8),
      legend.position = "right",
      plot.title      = element_text(face = "bold", size = 14),
      panel.spacing   = unit(1.5, "lines")
    )
  
  return(p)
}

# ============================================================================
# SCORE DISTRIBUTION HISTOGRAMS
# ============================================================================

plot_score_distributions <- function(results_df) {
  
  # Added colors and labels for Scores 4, 5, and 6
  score_colours <- c(
    "Score 1\n(Mean Luteal)"       = "#5a7d5a",
    "Score 2\n(Cyclical Severity)" = "#e8a66a",
    "Score 3\n(DSM-Weighted)"      = "#6C3BAA",
    "Score 4\n(GAMM Composite)"    = "#2C7FB8",
    "Score 5\n(GAMM AUC)"          = "#C51B8A",
    "Score 6\n(GAMM DSM-AUC)"      = "#1D91C0"
  )
  
  scores_long <- results_df %>%
    # Added score4, score5, score6 to select
    select(participant_id, algorithmic_pmdd, clinical_pmdd,
           score1, score2, score3, score4, score5, score6) %>%
    pivot_longer(
      # Added score4, score5, score6 to cols
      cols      = c(score1, score2, score3, score4, score5, score6),
      names_to  = "score_name",
      values_to = "score_value"
    ) %>%
    filter(!is.na(score_value)) %>%
    mutate(
      score_label = case_when(
        score_name == "score1" ~ "Score 1\n(Mean Luteal)",
        score_name == "score2" ~ "Score 2\n(Cyclical Severity)",
        score_name == "score3" ~ "Score 3\n(DSM-Weighted)",
        score_name == "score4" ~ "Score 4\n(GAMM Composite)",
        score_name == "score5" ~ "Score 5\n(GAMM AUC)",
        score_name == "score6" ~ "Score 6\n(GAMM DSM-AUC)"
      ),
      score_label = factor(score_label, levels = names(score_colours))
    )
  
  score_stats <- scores_long %>%
    group_by(score_label) %>%
    summarise(
      mean_val = mean(score_value, na.rm = TRUE),
      sd_val   = sd(score_value,   na.rm = TRUE),
      n        = n(),
      .groups  = "drop"
    )
  
  cat("\n=== PMDD SCORE STATISTICS ===\n")
  print(score_stats %>% mutate(mean_val = round(mean_val, 3), sd_val = round(sd_val, 3)))
  cat("=============================\n\n")
  
  # Plot 1: overall distribution
  p1 <- ggplot(scores_long, aes(x = score_value, fill = score_label)) +
    geom_histogram(bins = 30, colour = "white", linewidth = 0.3) +
    geom_vline(
      data = score_stats,
      aes(xintercept = mean_val),
      linetype = "dashed", colour = "grey30", linewidth = 0.7
    ) +
    geom_text(
      data = score_stats,
      aes(x = mean_val, y = Inf,
          label = paste0("M = ", round(mean_val, 2), "\nSD = ", round(sd_val, 2))),
      vjust = 1.5, hjust = -0.1, size = 3.5,
      fontface = "bold", colour = "grey30"
    ) +
    scale_fill_manual(values = score_colours) +
    facet_wrap(~score_label, scales = "free", ncol = 3) +
    labs(
      title    = "Distribution of PMDD Index Scores",
      subtitle = "Dashed line = mean. M = mean, SD = standard deviation.",
      x        = "Score Value",
      y        = "Number of Participants"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      legend.position  = "none",
      strip.text       = element_text(face = "bold", size = 11),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank()
    )
  
  # Plot 2: by algorithmic PMDD
  algo_stats <- scores_long %>%
    filter(!is.na(algorithmic_pmdd)) %>%
    mutate(algorithmic_pmdd = factor(algorithmic_pmdd,
                                     levels = c(0, 1),
                                     labels = c("No PMDD", "PMDD"))) %>%
    group_by(score_label, algorithmic_pmdd) %>%
    summarise(
      mean_val = mean(score_value, na.rm = TRUE),
      sd_val   = sd(score_value,   na.rm = TRUE),
      .groups  = "drop"
    )
  
  # Wilcoxon rank-sum p-values: PMDD vs No PMDD for each score (algorithmic)
  algo_pvals <- scores_long %>%
    filter(!is.na(algorithmic_pmdd)) %>%
    group_by(score_label) %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(score_value ~ algorithmic_pmdd)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = case_when(
        is.na(p_value)  ~ "p = NA",
        p_value < 0.001 ~ "p < 0.001",
        TRUE            ~ paste0("p = ", round(p_value, 3))
      ),
      stars = case_when(
        is.na(p_value) ~ "ns",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  p2 <- scores_long %>%
    filter(!is.na(algorithmic_pmdd)) %>%
    mutate(algorithmic_pmdd = factor(algorithmic_pmdd,
                                     levels = c(0, 1),
                                     labels = c("No PMDD", "PMDD"))) %>%
    ggplot(aes(x = score_value, fill = algorithmic_pmdd)) +
    geom_histogram(bins = 30, colour = "white", linewidth = 0.3,
                   position = "identity", alpha = 0.7) +
    geom_vline(
      data = algo_stats,
      aes(xintercept = mean_val, colour = algorithmic_pmdd),
      linetype = "dashed", linewidth = 0.7
    ) +
    scale_fill_manual(values = c("No PMDD" = "#a8c8a8", "PMDD" = "#d9534f")) +
    scale_colour_manual(values = c("No PMDD" = "#5a7d5a", "PMDD" = "#8b0000"),
                        guide = "none") +
    facet_wrap(~score_label, scales = "free", ncol = 3) +
    labs(
      title    = "Score Distributions by Algorithmic PMDD Diagnosis",
      subtitle = paste0(
        "Wilcoxon rank-sum test (PMDD vs No PMDD) —\n",
        paste(paste0(algo_pvals$score_label, ": ", algo_pvals$p_label), collapse = "  |  ")
      ),
      x        = "Score Value",
      y        = "Number of Participants",
      fill     = NULL
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      legend.position  = "bottom",
      strip.text       = element_text(face = "bold", size = 11),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank()
    ) +
    geom_text(
      data        = algo_pvals,
      aes(x = Inf, y = Inf, label = stars),
      hjust       = 1.3,
      vjust       = 1.5,
      size        = 7,
      fontface    = "bold",
      colour      = "black",
      inherit.aes = FALSE
    )
  
  cat("\n=== ALGORITHMIC PMDD GROUP STATISTICS ===\n")
  print(algo_stats %>% mutate(mean_val = round(mean_val, 3), sd_val = round(sd_val, 3)))
  cat("\n--- Wilcoxon Rank-Sum Test (PMDD vs No PMDD) ---\n")
  print(algo_pvals %>% select(score_label, p_value, p_label))
  cat("=========================================\n\n")
  
  # Plot 3: by clinical PMDD
  clin_stats <- scores_long %>%
    filter(!is.na(clinical_pmdd)) %>%
    mutate(clinical_pmdd = factor(clinical_pmdd,
                                  levels = c(0, 1),
                                  labels = c("No PMDD", "PMDD"))) %>%
    group_by(score_label, clinical_pmdd) %>%
    summarise(
      mean_val = mean(score_value, na.rm = TRUE),
      sd_val   = sd(score_value,   na.rm = TRUE),
      .groups  = "drop"
    )
  
  # Wilcoxon rank-sum p-values: PMDD vs No PMDD for each score (clinical)
  clin_pvals <- scores_long %>%
    filter(!is.na(clinical_pmdd)) %>%
    group_by(score_label) %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(score_value ~ clinical_pmdd)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = case_when(
        is.na(p_value)  ~ "p = NA",
        p_value < 0.001 ~ "p < 0.001",
        TRUE            ~ paste0("p = ", round(p_value, 3))
      ),
      stars = case_when(
        is.na(p_value)  ~ "ns",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  p3 <- scores_long %>%
    filter(!is.na(clinical_pmdd)) %>%
    mutate(clinical_pmdd = factor(clinical_pmdd,
                                  levels = c(0, 1),
                                  labels = c("No PMDD", "PMDD"))) %>%
    ggplot(aes(x = score_value, fill = clinical_pmdd)) +
    geom_histogram(bins = 30, colour = "white", linewidth = 0.3,
                   position = "identity", alpha = 0.7) +
    geom_vline(
      data = clin_stats,
      aes(xintercept = mean_val, colour = clinical_pmdd),
      linetype = "dashed", linewidth = 0.7
    ) +
    scale_fill_manual(values = c("No PMDD" = "#a8c8a8", "PMDD" = "#d9534f")) +
    scale_colour_manual(values = c("No PMDD" = "#5a7d5a", "PMDD" = "#8b0000"),
                        guide = "none") +
    facet_wrap(~score_label, scales = "free", ncol = 3) +
    labs(
      title    = "Score Distributions by Clinical PMDD Diagnosis",
      subtitle = paste0(
        "Wilcoxon rank-sum test (PMDD vs No PMDD) —\n",
        paste(paste0(clin_pvals$score_label, ": ", clin_pvals$p_label), collapse = "  |  ")
      ),
      x        = "Score Value",
      y        = "Number of Participants",
      fill     = NULL
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      legend.position  = "bottom",
      strip.text       = element_text(face = "bold", size = 11),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank()
    ) +
    geom_text(
      data        = clin_pvals,
      aes(x = Inf, y = Inf, label = stars),
      hjust       = 1.3,
      vjust       = 1.5,
      size        = 7,
      fontface    = "bold",
      colour      = "black",
      inherit.aes = FALSE
    )
  
  cat("\n=== CLINICAL PMDD GROUP STATISTICS ===\n")
  print(clin_stats %>% mutate(mean_val = round(mean_val, 3), sd_val = round(sd_val, 3)))
  cat("\n--- Wilcoxon Rank-Sum Test (PMDD vs No PMDD) ---\n")
  print(clin_pvals %>% select(score_label, p_value, p_label))
  cat("======================================\n\n")
  
  print(p1)
  print(p2)
  print(p3)
}

# ============================================================================
# WINDOW SYMPTOM SEVERITY PLOT (Mean ± SD)
# ============================================================================

plot_window_symptoms <- function(data_combined_filtered) {
  
  luteal_rows     <- list()
  follicular_rows <- list()
  
  for (pid in unique(data_combined_filtered$participant_id)) {
    df <- data_combined_filtered %>% filter(participant_id == pid)
    bleeding_starts <- find_bleeding_start_days(df)
    
    for (idx in bleeding_starts) {
      if (idx - 7 < 1)        next
      if (idx + 8 > nrow(df)) next
      
      df_luteal <- df[(idx - 7):(idx - 1), symptoms_list]
      df_folli  <- df[(idx + 2):(idx + 8), symptoms_list]
      
      if (any(is.na(df_luteal))) next
      if (any(is.na(df_folli)))  next
      
      luteal_rows[[length(luteal_rows) + 1]] <- df_luteal %>%
        summarise(across(everything(), mean)) %>%
        mutate(participant_id = pid)
      
      follicular_rows[[length(follicular_rows) + 1]] <- df_folli %>%
        summarise(across(everything(), mean)) %>%
        mutate(participant_id = pid)
    }
  }
  
  if (length(luteal_rows) == 0 || length(follicular_rows) == 0) {
    message("No valid cycles found for window plot")
    return(invisible(NULL))
  }
  
  summarise_window <- function(rows) {
    bind_rows(rows) %>%
      select(all_of(symptoms_list)) %>%
      pivot_longer(everything(), names_to = "symptom", values_to = "value") %>%
      group_by(symptom) %>%
      summarise(
        mean_val = mean(value, na.rm = TRUE),
        sd_val   = sd(value,   na.rm = TRUE),
        .groups  = "drop"
      ) %>%
      mutate(symptom = factor(symptom, levels = rev(symptoms_list)))
  }
  
  luteal_summary     <- summarise_window(luteal_rows)
  follicular_summary <- summarise_window(follicular_rows)
  
  make_plot <- function(summary_df, title) {
    ggplot(summary_df, aes(x = mean_val, y = symptom)) +
      # 1. Use the actual SD, no pmax manipulation
      geom_errorbarh(
        aes(xmin = mean_val - sd_val,
            xmax = mean_val + sd_val),
        height = 0.3, colour = "#a8c8f0", linewidth = 1
      ) +
      geom_point(colour = "#1a3a6b", size = 3) +
      # 2. Move text to the end of the error bar so it doesn't overlap the lines
      geom_text(
        aes(x = mean_val + sd_val, label = round(mean_val, 2)),
        hjust = -0.3, vjust = -0.5, size = 3, colour = "#1a3a6b"
      ) +
      geom_vline(xintercept = 3, colour = "red", linewidth = 0.8) +
      # 3. Set the breaks, but handle the limits safely in coord_cartesian
      scale_x_continuous(breaks = 1:6) +
      coord_cartesian(xlim = c(1, 6)) + 
      labs(title = title, x = "Mean ± SD", y = "Symptoms") +
      theme_minimal(base_family = "Helvetica") +
      theme(
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 13),
        axis.text        = element_text(size = 10)
      )
  }
  
  p_luteal     <- make_plot(luteal_summary,     "Luteal Window — Symptom Severity (Mean ± SD)")
  p_follicular <- make_plot(follicular_summary, "Follicular Window — Symptom Severity (Mean ± SD)")
  
  print(p_luteal)
  print(p_follicular)
  
  invisible(list(luteal = p_luteal, follicular = p_follicular))
}

# ============================================================================
# RESPONSE DISTRIBUTION PLOT
# ============================================================================

plot_response_distribution <- function(data, participant_ids = NULL,
                                       window = c("luteal", "follicular"),
                                       title = NULL) {
  window <- match.arg(window)
  
  window_rows <- list()
  
  df_all <- if (!is.null(participant_ids)) {
    data %>% filter(participant_id %in% participant_ids)
  } else {
    data
  }
  
  for (pid in unique(df_all$participant_id)) {
    df <- df_all %>% filter(participant_id == pid)
    bleeding_starts <- find_bleeding_start_days(df)
    
    for (idx in bleeding_starts) {
      if (idx - 7 < 1)        next
      if (idx + 8 > nrow(df)) next
      
      df_window <- if (window == "luteal") {
        df[(idx - 7):(idx - 1), symptoms_list]
      } else {
        df[(idx + 2):(idx + 8), symptoms_list]
      }
      
      if (any(is.na(df_window))) next
      
      window_rows[[length(window_rows) + 1]] <- df_window
    }
  }
  
  if (length(window_rows) == 0) {
    message("No valid cycles found")
    return(invisible(NULL))
  }
  
  plot_data <- bind_rows(window_rows) %>%
    pivot_longer(everything(), names_to = "symptom", values_to = "response") %>%
    filter(!is.na(response)) %>%
    mutate(
      response = factor(round(response), levels = 1:6),
      symptom  = factor(symptom, levels = rev(symptoms_list))
    ) %>%
    group_by(symptom, response) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(symptom) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup()
  
  window_label <- if (window == "luteal") "Luteal" else "Follicular"
  plot_title   <- if (!is.null(title)) title else paste(window_label, "Window — Response Distribution")
  
  ggplot(plot_data, aes(x = percentage, y = symptom, fill = response)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.7) +
    geom_text(
      data     = plot_data,
      aes(label = ifelse(percentage >= 8, paste0(round(percentage, 1), "%"), "")),
      position = position_stack(vjust = 0.5, reverse = TRUE),
      size     = 5.5,
      color    = "black"
    ) +
    geom_vline(xintercept = 50, colour = "red", linetype = "dashed", linewidth = 0.7) +
    scale_fill_manual(
      values = c(
        "1" = "#deeef7",
        "2" = "#b8d9ee",
        "3" = "#7ab8d9",
        "4" = "#3d8fbf",
        "5" = "#1a5f8a",
        "6" = "#0a2f4a"
      ),
      name = "Response (1–6)"
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = c("0", "25", "50", "75", "100")
    ) +
    labs(title = plot_title, x = "Percentage (%)", y = "Symptom") +
    theme_minimal(base_family = "Helvetica") +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title         = element_text(face = "bold", size = 13),
      axis.text          = element_text(size = 11),
      legend.position    = "right"
    )
}

# ============================================================================
# GENERATE ALL FOUR RESPONSE DISTRIBUTION PLOTS
# ============================================================================

plot_all_response_distributions <- function(data_combined_filtered, two_cycle_ids) {
  p1 <- plot_response_distribution(
    data_combined_filtered,
    participant_ids = NULL,
    window          = "luteal",
    title           = "All Participants — Luteal Window Response Distribution"
  )
  
  p2 <- plot_response_distribution(
    data_combined_filtered,
    participant_ids = NULL,
    window          = "follicular",
    title           = "All Participants — Follicular Window Response Distribution"
  )
  
  p3 <- plot_response_distribution(
    data_combined_filtered,
    participant_ids = two_cycle_ids,
    window          = "luteal",
    title           = "2 Evaluable Cycles — Luteal Window Response Distribution"
  )
  
  p4 <- plot_response_distribution(
    data_combined_filtered,
    participant_ids = two_cycle_ids,
    window          = "follicular",
    title           = "2 Evaluable Cycles — Follicular Window Response Distribution"
  )
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
  invisible(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

# ============================================================================
# RESPONSE FREQUENCY TABLE
# ============================================================================

make_response_table <- function(data, participant_ids = NULL, window = c("luteal", "follicular")) {
  window <- match.arg(window)
  
  window_rows <- list()
  
  df_all <- if (!is.null(participant_ids)) {
    data %>% filter(participant_id %in% participant_ids)
  } else {
    data
  }
  
  for (pid in unique(df_all$participant_id)) {
    df <- df_all %>% filter(participant_id == pid)
    bleeding_starts <- find_bleeding_start_days(df)
    
    for (idx in bleeding_starts) {
      if (idx - 7 < 1)        next
      if (idx + 8 > nrow(df)) next
      
      df_window <- if (window == "luteal") {
        df[(idx - 7):(idx - 1), symptoms_list]
      } else {
        df[(idx + 2):(idx + 8), symptoms_list]
      }
      
      if (any(is.na(df_window))) next
      
      window_rows[[length(window_rows) + 1]] <- df_window
    }
  }
  
  if (length(window_rows) == 0) {
    message("No valid cycles found")
    return(invisible(NULL))
  }
  
  bind_rows(window_rows) %>%
    pivot_longer(everything(), names_to = "symptom", values_to = "response") %>%
    filter(!is.na(response)) %>%
    mutate(
      response = factor(round(response), levels = 1:6),
      symptom  = factor(symptom, levels = symptoms_list)
    ) %>%
    group_by(symptom, response) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(symptom) %>%
    mutate(percentage = round(n / sum(n) * 100, 1)) %>%
    ungroup() %>%
    mutate(cell = paste0(percentage, "% (", n, ")")) %>%
    select(symptom, response, cell) %>%
    pivot_wider(names_from = response, values_from = cell)
}

print_all_response_tables <- function(data_combined_filtered, two_cycle_ids) {
  cat("\n=== ALL PARTICIPANTS — LUTEAL WINDOW ===\n")
  make_response_table(data_combined_filtered, participant_ids = NULL, window = "luteal") %>% print()
  
  cat("\n=== ALL PARTICIPANTS — FOLLICULAR WINDOW ===\n")
  make_response_table(data_combined_filtered, participant_ids = NULL, window = "follicular") %>% print()
  
  cat("\n=== 2 EVALUABLE CYCLES — LUTEAL WINDOW ===\n")
  make_response_table(data_combined_filtered, participant_ids = two_cycle_ids, window = "luteal") %>% print()
  
  cat("\n=== 2 EVALUABLE CYCLES — FOLLICULAR WINDOW ===\n")
  make_response_table(data_combined_filtered, participant_ids = two_cycle_ids, window = "follicular") %>% print()
}

# ============================================================================
# YES/NO CLINICAL SIGNIFICANCE PLOT
# ============================================================================

plot_yes_no_distribution <- function(data, participant_ids = NULL, title = NULL) {
  
  df_all <- if (!is.null(participant_ids)) {
    data %>% filter(participant_id %in% participant_ids)
  } else {
    data
  }
  
  plot_data <- df_all %>%
    select(all_of(symptoms_list)) %>%
    pivot_longer(everything(), names_to = "symptom", values_to = "response") %>%
    filter(!is.na(response)) %>%
    mutate(
      response = ifelse(response >= 4, "yes", "no"),
      response = factor(response, levels = c("no", "yes")),
      symptom  = factor(symptom, levels = rev(symptoms_list))
    ) %>%
    group_by(symptom, response) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(symptom) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup()
  
  plot_title <- if (!is.null(title)) title else "% Yes/No per symptom"
  
  ggplot(plot_data, aes(x = percentage, y = symptom, fill = response)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.7) +
    geom_text(
      aes(label = ifelse(percentage >= 8, paste0(round(percentage, 1), "%"), "")),
      position = position_stack(vjust = 0.5, reverse = TRUE),
      size     = 5.5,
      color    = "black"
    ) +
    scale_fill_manual(
      values = c("yes" = "#4472C4", "no" = "#D9D9D9"),
      name   = "Response"
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = c("0", "25", "50", "75", "100")
    ) +
    labs(title = plot_title, x = "Percentage (%)", y = "Symptom") +
    theme_minimal(base_family = "Helvetica") +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title         = element_text(face = "bold", size = 13),
      axis.text          = element_text(size = 10),
      legend.position    = "right"
    )
}

plot_all_yes_no <- function(data_combined_filtered, two_cycle_ids, one_cycle_ids) {
  p1 <- plot_yes_no_distribution(
    data_combined_filtered,
    participant_ids = NULL,
    title           = "% Yes/No per symptom (all participants)"
  )
  
  p2 <- plot_yes_no_distribution(
    data_combined_filtered,
    participant_ids = two_cycle_ids,
    title           = "% Yes/No per symptom (2 evaluable cycles)"
  )
  
  p3 <- plot_yes_no_distribution(
    data_combined_filtered,
    participant_ids = one_cycle_ids,
    title           = "% Yes/No per symptom (1 evaluable cycle)"
  )
  
  print(p1)
  print(p2)
  print(p3)
  
  invisible(list(p1 = p1, p2 = p2, p3 = p3))
}

source("scripts/config.R")

# ============================================================================
# ROC CURVE ANALYSIS
# Compares Scores 1-6 against:
#   (A) clinical_pmdd  — self-reported diagnosis
#   (B) algorithmic_pmdd — check_pmdd() algorithm
# Requires the pROC package
# ============================================================================

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

# ============================================================================
# HELPER: compute ROC and extract AUC + CI for one score vs one outcome
# ============================================================================

compute_roc <- function(data, score_col, outcome_col) {
  df <- data %>%
    select(score = all_of(score_col), outcome = all_of(outcome_col)) %>%
    filter(!is.na(score) & !is.na(outcome))
  
  if (nrow(df) == 0 || length(unique(df$outcome)) < 2) return(NULL)
  
  roc_obj <- roc(df$outcome, df$score,
                 levels    = c(0, 1),
                 direction = "<",   # higher score = more likely PMDD
                 quiet     = TRUE)
  
  ci_obj  <- ci.auc(roc_obj, conf.level = 0.95)
  
  list(
    roc     = roc_obj,
    auc     = as.numeric(auc(roc_obj)),
    ci_low  = ci_obj[1],
    ci_high = ci_obj[3],
    n       = nrow(df),
    score   = score_col,
    outcome = outcome_col
  )
}

# ============================================================================
# HELPER: plot ROC curves on one panel dynamically
# ============================================================================

plot_roc_panel <- function(roc_list, outcome_label) {
  
  # Colour palette expanded for 6 scores
  cols <- c(
    score1 = "#5a7d5a", 
    score2 = "#e8a66a", 
    score3 = "#6C3BAA",
    score4 = "#1F78B4", # Blue
    score5 = "#E31A1C", # Red
    score6 = "#B15928"  # Copper
  )
  
  # Function to build tidy dataframe for a single ROC result
  roc_to_df <- function(roc_result, score_name) {
    # Generate a clean label like "Score 1 (AUC = 0.850)"
    clean_name <- gsub("score", "Score ", score_name)
    label <- sprintf("%s (AUC = %.3f)", clean_name, roc_result$auc)
    
    tibble(
      fpr          = 1 - roc_result$roc$specificities,
      tpr          = roc_result$roc$sensitivities,
      score_id     = score_name,
      legend_label = label
    )
  }
  
  # Bind all available ROC objects into one dataframe
  roc_df <- map_dfr(names(roc_list), ~ roc_to_df(roc_list[[.x]], .x))
  
  # Extract the dynamic color mapping for the legend
  mapping_df <- unique(roc_df[, c("score_id", "legend_label")])
  color_mapping <- setNames(cols[mapping_df$score_id], mapping_df$legend_label)
  
  ggplot(roc_df, aes(x = fpr, y = tpr, colour = legend_label)) +
    geom_line(linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey60", linewidth = 0.7) +
    scale_colour_manual(values = color_mapping) +
    scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(
      title    = outcome_label,
      x        = "False Positive Rate (1 - Specificity)",
      y        = "True Positive Rate (Sensitivity)",
      colour   = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )
}

# ============================================================================
# RUN ROC ANALYSIS
# ============================================================================

run_roc_analysis <- function(results_df, dataset_label = "") {
  
  # Updated to include all 6 scores
  scores   <- paste0("score", 1:6) 
  outcomes <- c("clinical_pmdd", "algorithmic_pmdd")
  
  # Print AUC summary table
  cat("\n=== ROC AUC SUMMARY ===\n")
  
  for (outcome in outcomes) {
    
    cat("\nOutcome:", outcome, "\n")
    cat(sprintf("  %-10s  AUC    95%% CI            N\n", "Score"))
    cat("  ", strrep("-", 45), "\n", sep = "")
    
    roc_list <- list()
    
    for (score in scores) {
      r <- compute_roc(results_df, score, outcome)
      if (is.null(r)) {
        cat(sprintf("  %-10s  insufficient data\n", score))
        next
      }
      roc_list[[score]] <- r
      cat(sprintf("  %-10s  %.3f  [%.3f – %.3f]  n=%d\n",
                  score, r$auc, r$ci_low, r$ci_high, r$n))
    }
    
    # Plot as long as we have at least one valid ROC curve to show
    if (length(roc_list) > 0) {
      outcome_label <- ifelse(outcome == "clinical_pmdd",
                              "Clinical PMDD (self-reported)",
                              "Algorithmic PMDD (check_pmdd)")
      plot_title <- paste(outcome_label, "-", dataset_label)
      
      # Pass the entire list at once
      p <- plot_roc_panel(roc_list, plot_title)
      print(p)
    }
  }
}

