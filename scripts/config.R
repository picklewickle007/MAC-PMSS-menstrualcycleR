# ============================================================================
# GITHUB AND PACKAGE LOADING
# ============================================================================
if (exists("config_loaded", envir = .GlobalEnv)) {
  message("Config already loaded")
} else {
  config_loaded <- TRUE
  
  if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
  }
  pacman::p_load(tidyverse, pROC) # install packages 
  
  if(!requireNamespace("menstrualcycleR", quietly = TRUE)) {
    install.packages("menstrualcycleR")
  }
  library(menstrualcycleR)
  
  if(!requireNamespace("mgcv", quietly = TRUE)) {
    install.packages("mgcv")
  }
  library(mgcv)
  
  if(!requireNamespace("marginaleffects", quietly = TRUE)) {
    install.packages("marginaleffects")
  }
  library(marginaleffects)

  
  GITHUB_TOKEN <- readLines("~/Documents/github_token.txt")[1]
  GITHUB_USER  <- "bengid07"
  REPO_NAME    <- "MACPMSS_complete_dataset"
  local_path   <- "MACPMSS_complete_dataset"
  
  if (!dir.exists(local_path)) {
    repo_url <- paste0("https://", GITHUB_TOKEN, "@github.com/", GITHUB_USER, "/", REPO_NAME, ".git")
    system(paste("git clone", repo_url, local_path))
  }
  
  system(paste("git -C", local_path, "pull"))
  
  data_android_ori    <- read_csv(file.path(local_path, "firebase_data_android-monitorization.csv"), show_col_types = FALSE)
  data_ios_ori        <- read_csv(file.path(local_path, "firebase_data_ios-monitorization.csv"), show_col_types = FALSE)
  profile_android_ori <- read_csv(file.path(local_path, "firebase_data_android-profile.csv"), show_col_types = FALSE)
  profile_ios_ori     <- read_csv(file.path(local_path, "firebase_data_ios-profile.csv"), show_col_types = FALSE)
}

# ============================================================================
# CONFIGURATION PARAMETERS
# ============================================================================
symptoms_list <- c("lability", "anger", "depressed", "anxiety", "activity",
                   "concentration", "lethargic", "appetite", "sleep",
                   "overwhelm", "physical")

MIN_CONTINUOUS_DAYS <- 16
FOLLICULAR_METHOD   <- "max"

