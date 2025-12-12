####---
###     MR Analysis Based on GWAS-derived Exposure Factors (Loop Version) ----------------------------------------
####---

### Outcome processed
rm(list = ls())
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(ggplot2)
library(openxlsx)

# Set fixed working directory path
target_dir <- "E:/20250303 TB NHANES/20251209 IJS revise/TSMR"
setwd(target_dir)

target_outcome <- c("TB") #### Replace as needed

# Automatically get ZIP file and process path
file_list <- list.files(path = "data/outcome", pattern = "\\.gz$", full.names = TRUE)
zip_file <- file_list[2]
stopifnot("No gz file found" = !is.na(zip_file))

library(tools)
zip_name <- file_path_sans_ext(zip_file)
file_name <- basename(zip_name)
print(file_name)

# Create new folder
dir.create(file_name, showWarnings = FALSE)
setwd(file_name)

# Read outcome data
outcome_dat <- readRDS("outcome_dat.RDS")

# Exposure factors folder
target_folder <- "E:/20250303 TB NHANES/20251209 IJS revise/TSMR/data/exposure"
setwd(target_folder)

# Get all RDS files
rds_file_list <- list.files(path = target_folder, 
                            pattern = "\\.RDS$", 
                            full.names = TRUE, 
                            recursive = TRUE)

print(paste("Found", length(rds_file_list), "RDS files"))

# Set loop range
start_index <- 1
end_index <- 8

# Ensure indices are within valid range
start_index <- max(1, min(start_index, length(rds_file_list)))
end_index <- min(length(rds_file_list), max(end_index, start_index))

# Create log directory
log_dir <- file.path(target_dir, file_name, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, append = TRUE, split = TRUE)

cat("=== MR Analysis Started ===\n")
cat("Start time:", format(Sys.time()), "\n")
cat("Analysis range: Exposure factor", start_index, "to", end_index, "\n")

# Initialize data frames to store Egger and PRESSO results across all exposures
egger_intercept_df <- data.frame()
presso_global_df <- data.frame()
presso_corrected_df <- data.frame()

# Smart SNP filtering function
smart_snp_filtering <- function(dat_clean, target_min = 10, target_max = 100) {
  # Define p-value threshold sequence (only stricter to stricter, less is less)
  pval_thresholds <- c(5e-8,5e-9, 5e-10,  # Stricter
                       5e-11,5e-12, 5e-13,  
                       5e-14,5e-15, 5e-16,  
                       5e-17,5e-18, 5e-19,  
                       5e-20,5e-21, 5e-22,  
                       5e-23,5e-24, 5e-25,  
                       5e-26,5e-27, 5e-28,  
                       5e-29,5e-30, 5e-31,  
                       5e-32,5e-33, 5e-34)   
  best_filtered <- NULL
  best_threshold <- 5e-8
  best_count <- 0
  
  # First try default threshold
  dat_filtered <- dat_clean %>% filter(pval.exposure < 5e-8)
  current_count <- nrow(dat_filtered)
  
  cat("Initial threshold 5e-8 yields SNP count:", current_count, "\n")
  
  # If initial count is within target range, return directly
  if (current_count >= target_min & current_count <= target_max) {
    return(list(data = dat_filtered, threshold = 5e-8, note = "Using default threshold"))
  }
  
  # If count is too high, try stricter thresholds
  if (current_count > target_max) {
    cat("SNP count too high, trying stricter thresholds...\n")
    for (pval_thresh in pval_thresholds[pval_thresholds < 5e-8]) {
      temp_filtered <- dat_clean %>% filter(pval.exposure < pval_thresh)
      temp_count <- nrow(temp_filtered)
      cat("Threshold", format(pval_thresh, scientific = TRUE), ":", temp_count, "SNPs\n")
      
      if (temp_count >= target_min & temp_count <= target_max) {
        best_filtered <- temp_filtered
        best_threshold <- pval_thresh
        best_count <- temp_count
        break
      } else if (temp_count <= target_max & temp_count > best_count) {
        best_filtered <- temp_filtered
        best_threshold <- pval_thresh
        best_count <- temp_count
      }
    }
  }
  
  # If count is too low, try more lenient thresholds (only when MR-PRESSO is needed)
  if (current_count < target_min) {
    cat("SNP count insufficient, trying more lenient thresholds...\n")
    for (pval_thresh in pval_thresholds[pval_thresholds > 5e-8]) {
      temp_filtered <- dat_clean %>% filter(pval.exposure < pval_thresh)
      temp_count <- nrow(temp_filtered)
      cat("Threshold", format(pval_thresh, scientific = TRUE), ":", temp_count, "SNPs\n")
      
      if (temp_count >= target_min & temp_count <= target_max) {
        best_filtered <- temp_filtered
        best_threshold <- pval_thresh
        best_count <- temp_count
        break
      } else if (temp_count >= target_min & (best_count == 0 | temp_count < best_count)) {
        best_filtered <- temp_filtered
        best_threshold <- pval_thresh
        best_count <- temp_count
      }
      
      # Stop if upper limit reached
      if (temp_count > target_max * 2) break
    }
  }
  
  # If no ideal range found, use the closest value to target
  if (is.null(best_filtered)) {
    if (current_count > 0) {
      best_filtered <- dat_filtered
      best_threshold <- 5e-8
      best_count <- current_count
      cat("Using default threshold, SNP count:", best_count, "\n")
    } else {
      # If even lenient thresholds don't yield enough SNPs, use the most lenient
      most_lenient <- max(pval_thresholds)
      best_filtered <- dat_clean %>% filter(pval.exposure < most_lenient)
      best_threshold <- most_lenient
      best_count <- nrow(best_filtered)
      cat("Using most lenient threshold", format(most_lenient, scientific = TRUE), 
          ", SNP count:", best_count, "\n")
    }
  }
  
  note <- paste("Optimized threshold:", format(best_threshold, scientific = TRUE),
                "SNP count:", best_count)
  
  return(list(data = best_filtered, threshold = best_threshold, note = note))
}

# Main loop
for(i in start_index:end_index) {
  targetrds <- rds_file_list[i]
  exposure_name <- tools::file_path_sans_ext(basename(targetrds))
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Processing exposure factor", i, ":", exposure_name, "\n")
  cat("Start time:", format(Sys.time()), "\n")
  
  result <- tryCatch({
    # Read exposure data
    exposure_dat <- readRDS(targetrds)
    cat("Successfully read exposure data, containing", nrow(exposure_dat), "SNPs\n")
    
    # Check common SNPs
    common_snps <- intersect(exposure_dat$SNP, outcome_dat$SNP)
    cat("Common SNP count:", length(common_snps), "\n")
    
    if(length(common_snps) < 3) {
      warning("Insufficient common SNPs, skipping this exposure factor")
      next
    }
    
    # Data harmonization
    harmonised_dat <- harmonise_data(
      exposure_dat = exposure_dat, 
      outcome_dat = outcome_dat, 
      action = 2
    )
    cat("After harmonization, retained", nrow(harmonised_dat), "SNPs\n")
    
    # Data cleaning
    dat_clean <- na.omit(harmonised_dat[,-c(18:20)])
    dat_clean <- subset(dat_clean, se.exposure > 1e-10 & se.outcome > 1e-10)
    
    # Smart SNP filtering
    cat("Performing smart SNP filtering...\n")
    filtering_result <- smart_snp_filtering(dat_clean, target_min = 10, target_max = 100)
    dat_filtered <- filtering_result$data
    final_threshold <- filtering_result$threshold
    cat(filtering_result$note, "\n")
    
    num_snps <- nrow(dat_filtered)
    if(num_snps < 3) {
      stop("Insufficient SNPs for analysis")
    }
    
    # Create result folder
    result_folder <- file.path(target_dir, file_name, basename(dirname(targetrds)))
    dir.create(result_folder, recursive = TRUE, showWarnings = FALSE)
    
    # Skip MR-Egger and MR-PRESSO if SNPs < 3
    if (num_snps < 3) {
      cat("Skipping", exposure_name, ": SNPs < 3, cannot perform MR-Egger and MR-PRESSO\n")
      next
    }
    
    ### 1. MR-Egger intercept test (with error handling)
    pleiotropy_test <- tryCatch({
      mr_pleiotropy_test(dat_filtered)
    }, error = function(e) {
      warning(paste("MR-Egger failed for", exposure_name, ":", e$message))
      return(list(egger_intercept = NA, se = NA, pval = NA))
    })
    
    # Store Egger results
    egger_intercept_df <- rbind(egger_intercept_df, 
                                data.frame(
                                  id.exposure = exposure_name,
                                  egger_intercept = pleiotropy_test$egger_intercept,
                                  se = pleiotropy_test$se,
                                  pval = pleiotropy_test$pval,
                                  stringsAsFactors = FALSE
                                ))
    
    ### 2. MR-PRESSO analysis (only run when SNPs ≥ 10)
    presso_result <- NULL
    if (num_snps >= 10) {
      cat("Running MR-PRESSO analysis...\n")
      presso_result <- tryCatch({
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          data = dat_filtered,
          NbDistribution = 1000,
          SignifThreshold = 0.05
        )
      }, error = function(e) {
        warning(paste("MR-PRESSO failed for", exposure_name, ":", e$message))
        return(NULL)
      })
      
      # Store PRESSO results
      if (!is.null(presso_result)) {
        # Global test results
        global_test <- presso_result$`MR-PRESSO results`$`Global Test`
        presso_global_df <- rbind(presso_global_df, 
                                  data.frame(
                                    id.exposure = exposure_name,
                                    RSSobs = global_test$RSSobs,
                                    Pvalue = global_test$Pvalue,
                                    stringsAsFactors = FALSE
                                  ))
        
        # Corrected MR results
        corrected_res <- presso_result$`Main MR results`[2, ]
        presso_corrected_df <- rbind(presso_corrected_df, 
                                     data.frame(
                                       id.exposure = exposure_name,
                                       Causal_Estimate = ifelse(is.null(corrected_res$`Causal Estimate`), 
                                                                NA, corrected_res$`Causal Estimate`),
                                       Sd = ifelse(is.null(corrected_res$Sd), NA, corrected_res$Sd),
                                       T_stat = ifelse(is.null(corrected_res$`T-stat`), NA, corrected_res$`T-stat`),
                                       P_value = ifelse(is.null(corrected_res$`P-value`), NA, corrected_res$`P-value`),
                                       stringsAsFactors = FALSE
                                     ))
      } else {
        # Record PRESSO failure
        presso_global_df <- rbind(presso_global_df, 
                                  data.frame(
                                    id.exposure = exposure_name,
                                    RSSobs = NA,
                                    Pvalue = NA,
                                    stringsAsFactors = FALSE
                                  ))
        presso_corrected_df <- rbind(presso_corrected_df, 
                                     data.frame(
                                       id.exposure = exposure_name,
                                       Causal_Estimate = NA,
                                       Sd = NA,
                                       T_stat = NA,
                                       P_value = NA,
                                       stringsAsFactors = FALSE
                                     ))
      }
    } else {
      cat("Skipping MR-PRESSO: SNPs < 10\n")
      # Record NA for PRESSO when not run
      presso_global_df <- rbind(presso_global_df, 
                                data.frame(
                                  id.exposure = exposure_name,
                                  RSSobs = NA,
                                  Pvalue = NA,
                                  stringsAsFactors = FALSE
                                ))
      presso_corrected_df <- rbind(presso_corrected_df, 
                                   data.frame(
                                     id.exposure = exposure_name,
                                     Causal_Estimate = NA,
                                     Sd = NA,
                                     T_stat = NA,
                                     P_value = NA,
                                     stringsAsFactors = FALSE
                                   ))
    }
    
    # MR analysis
    cat("Running MR analysis...\n")
    mr_result <- mr(dat_filtered)
    mr_result <- generate_odds_ratios(mr_result)
    
    # Other tests
    heterogeneity <- mr_heterogeneity(dat_filtered)
    pleiotropy <- mr_pleiotropy_test(dat_filtered)
    singlesnp <- mr_singlesnp(dat_filtered)
    leaveoneout <- mr_leaveoneout(dat_filtered)
    
    # Save filter information
    filter_info <- data.frame(
      exposure = exposure_name,
      initial_snps = nrow(harmonised_dat),
      cleaned_snps = nrow(dat_clean),
      final_snps = nrow(dat_filtered),
      pval_threshold = final_threshold,
      mr_presso_performed = !is.null(presso_result),
      timestamp = Sys.time()
    )
    
    # Save RDS files
    saveRDS(harmonised_dat, file.path(result_folder, "harmonised_data.RDS"))
    saveRDS(dat_filtered, file.path(result_folder, "filtered_data.RDS"))
    saveRDS(mr_result, file.path(result_folder, "mr_results.RDS"))
    saveRDS(filter_info, file.path(result_folder, "filter_info.RDS"))
    
    if(!is.null(presso_result)) {
      saveRDS(presso_result, file.path(result_folder, "mr_presso_results.RDS"))
    }
    
    # Create three separate Excel files as requested
    
    ## File 1: Main MR Results
    main_results_list <- list(
      "Main Results" = mr_result
    )
    write.xlsx(main_results_list, file.path(result_folder, "MR_Main_Results.xlsx"))
    
    ## File 2: Sensitivity Tests
    sensitivity_list <- list(
      "Heterogeneity" = heterogeneity,
      "Pleiotropy" = pleiotropy,
      "LeaveOneOut" = leaveoneout,
      "SingleSNP" = singlesnp
    )
    write.xlsx(sensitivity_list, file.path(result_folder, "MR_Sensitivity_Tests.xlsx"))
    
    ## File 3: Egger and PRESSO Results
    # Prepare Egger results for this exposure
    egger_current <- egger_intercept_df[egger_intercept_df$id.exposure == exposure_name, ]
    # Prepare PRESSO results for this exposure
    presso_global_current <- presso_global_df[presso_global_df$id.exposure == exposure_name, ]
    presso_corrected_current <- presso_corrected_df[presso_corrected_df$id.exposure == exposure_name, ]
    
    egger_presso_list <- list(
      "Egger_Intercept" = egger_current,
      "PRESSO_Global" = presso_global_current,
      "PRESSO_Corrected" = presso_corrected_current
    )
    write.xlsx(egger_presso_list, file.path(result_folder, "MR_Egger_PRESSO.xlsx"))
    
    cat("Analysis for exposure factor", i, "completed\n")
    cat("Results saved to:", result_folder, "\n")
    
    list(status = "success", exposure = exposure_name, snp_count = nrow(dat_filtered),
         threshold = final_threshold, presso_performed = !is.null(presso_result))
    
  }, error = function(e) {
    error_msg <- paste("Error:", e$message)
    cat(error_msg, "\n")
    cat("Exposure factor", i, "processing failed:", exposure_name, "-", error_msg, "\n")
    list(status = "error", exposure = exposure_name, error = e$message)
  })
  
  cat("End time:", format(Sys.time()), "\n")
  cat("Processing status:", ifelse(result$status == "success", "Success", "Failed"), "\n")
  
  if(result$status == "success") {
    cat("SNPs used:", result$snp_count, "\n")
    cat("Final p-value threshold:", format(result$threshold, scientific = TRUE), "\n")
    cat("MR-PRESSO performed:", ifelse(result$presso_performed, "Yes", "No"), "\n")
  } else {
    cat("Error message:", result$error, "\n")
  }
  
  cat(rep("=", 60), "\n\n", sep = "")
}

# Save combined Egger and PRESSO results across all exposures
if (nrow(egger_intercept_df) > 0) {
  combined_egger_presso <- list(
    "Egger_Intercept" = egger_intercept_df,
    "PRESSO_Global" = presso_global_df,
    "PRESSO_Corrected" = presso_corrected_df
  )
  write.xlsx(combined_egger_presso, file.path(target_dir, file_name, "Combined_Egger_PRESSO_Results.xlsx"))
}

# Close log file
cat("\n=== MR Analysis Completed ===\n")
cat("End time:", format(Sys.time()), "\n")
sink()

# Restore console output
sink()

cat("Analysis completed!\n")
cat("Log file saved to:", log_file, "\n")
cat("Total exposure factors processed:", (end_index - start_index + 1), "\n")