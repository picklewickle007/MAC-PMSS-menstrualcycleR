# ==========================================
# FUNCTION 1: Fit the GAMM Model
# ==========================================
fit_cyclic_gamm <- function(data, target_symptom) {
  
  # 1. Prepare and transform the data
  data$symptom_log <- log(data[[target_symptom]])
  data$id <- as.factor(data$id)
  
  # 2. Filter for complete cases
  selected_vars <- c("cyclic_time_impute", "symptom_log", "id")
  datSX <- data[complete.cases(data[selected_vars]), ]
  
  # 3. Fit the model
  gamm_model <- gamm(
    symptom_log ~ 
      s(cyclic_time_impute, bs = "cc") + 
      s(id, bs = 're') + 
      s(cyclic_time_impute, id, bs=c("re", "cc")), 
    knots = list(cyclic_time_impute = c(-1, 1)), 
    data = datSX,
    method = 'REML'
  )
  
  # Return both the model and the cleaned data (needed for plotting)
  return(list(
    model = gamm_model, 
    clean_data = datSX, 
    symptom_name = target_symptom
  ))
}

# ==========================================
# FUNCTION 2: Plot the GAMM Model
# ==========================================
plot_cyclic_gamm <- function(gamm_output) {
  
  # Extract components from the fitted output
  gamm_model <- gamm_output$model
  datSX <- gamm_output$clean_data
  symptom_name <- gamm_output$symptom_name
  
  # 1. Create plot data using the first valid ID as a placeholder
  plotdat <- expand.grid(
    cyclic_time_impute = seq(-1, 1, by = 0.05),
    id = datSX$id[1] 
  )
  
  # 2. Generate Predictions
  preds <- predict(
    gamm_model$gam, 
    newdata = plotdat, 
    type = "link", 
    se.fit = TRUE,
    exclude = c("s(id)", "s(cyclic_time_impute,id)") 
  )
  
  # 3. Reverse the log(x+1) transformation
  plotdat$estimate  <- exp(preds$fit)
  plotdat$conf.low  <- exp(preds$fit - (1.96 * preds$se.fit)) 
  plotdat$conf.high <- exp(preds$fit + (1.96 * preds$se.fit)) 
  
  # 4. Build the plot (Changed color = "white" to NA to avoid vertical lines)
  gamplot <- ggplot(plotdat, aes(x = cyclic_time_impute, y = estimate)) +
    geom_rect(xmin = -0, xmax = 0.08, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.2, color = NA) +
    geom_rect(xmin = 0.92, xmax = 1, ymin = -Inf, ymax = Inf, fill = "grey87", alpha = 0.2, color = NA) +
    geom_rect(xmin = -1, xmax = -0.92, ymin = -Inf, ymax = Inf, fill = "grey87", alpha = 0.2, color = NA) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "lightgrey", alpha = 0.3) +
    geom_line(linewidth = 1, color = "black", show.legend = FALSE) +
    scale_x_continuous(
      limits = c(-1, 1), 
      breaks = seq(-1, 1, by = 0.50),
      labels = c("Ovulation", "50%L", "Menses Onset", "50%F", "Ovulation")
    ) + 
    # Dynamically capitalize the symptom name for the y-axis
    labs(x = "Menstrual Cycle Phase", y = tools::toTitleCase(symptom_name)) +
    theme_minimal()
  
  return(gamplot)
}

# ==========================================
# FUNCTION 3: Batch Fit Models & Generate Predictions
# ==========================================
fit_and_predict_all_symptoms <- function(data, symptoms_list) {
  
  # 1. Setup the master predictions grid for all participants
  # We extract unique IDs directly from your raw data
  unique_ids <- unique(as.character(data$id)) 
  
  predictions_df <- expand.grid(
    cycle_day = seq(-14, 13, by = 1), 
    participant_id = unique_ids
  )
  
  # Standardize the days to match the GAMM's -1 to 1 scale
  # and duplicate participant_id to 'id' so the GAMM recognizes it
  predictions_df$cyclic_time_impute <- predictions_df$cycle_day / 14
  predictions_df$id <- as.factor(predictions_df$participant_id)
  
  # 2. Setup an empty list to hold the model outputs
  fitted_models_list <- list()
  
  # 3. Loop through each symptom
  for (symp in symptoms_list) {
    message(paste("Fitting GAMM for:", symp, "..."))
    
    # Run your existing Function 1
    gamm_output <- fit_cyclic_gamm(data, target_symptom = symp)
    
    # Save the output to the list (so you can still plot them later!)
    fitted_models_list[[symp]] <- gamm_output
    
    # Generate individualized predictions for this specific symptom
    preds <- predict(
      gamm_output$model$gam, 
      newdata = predictions_df, 
      type = "link", 
      se.fit = FALSE
    )
    
    # Back-transform from log(x+1) and attach it as a new column
    predictions_df[[symp]] <- exp(preds) 
  }
  
  # 4. Clean up the temporary 'id' column
  predictions_df$id <- NULL
  
  # Return BOTH the list of models and the master prediction dataframe
  return(list(
    models = fitted_models_list,
    predictions_df = predictions_df
  ))
}