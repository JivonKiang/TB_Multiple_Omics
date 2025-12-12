####---
###     基于GWAS获得的暴露因素的分析（循环版本）----------------------------------------
####---

### outcome 处理好了
rm(list = ls())
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(ggplot2)
library(openxlsx)

# 设置固定工作目录路径
target_dir <- "E:/20250303 TB NHANES/20251209 IJS revise/TSMR"
setwd(target_dir)

target_outcome <- c("TB") ####根据需求替换

# 自动获取ZIP文件并处理路径
file_list <- list.files(path = "data/outcome", pattern = "\\.gz$", full.names = TRUE)
zip_file <- file_list[1]
stopifnot("未找到gz文件" = !is.na(zip_file))

library(tools)
zip_name <- file_path_sans_ext(zip_file)
file_name <- basename(zip_name)
print(file_name)

# 创建新文件夹
dir.create(file_name, showWarnings = FALSE)
setwd(file_name)

# 读取结局数据
outcome_dat <- readRDS("outcome_dat.RDS")

# 暴露因素文件夹
target_folder <- "E:/20250303 TB NHANES/20251209 IJS revise/TSMR/data/exposure"
setwd(target_folder)

# 获取所有RDS文件
rds_file_list <- list.files(path = target_folder, 
                            pattern = "\\.RDS$", 
                            full.names = TRUE, 
                            recursive = TRUE)

print(paste("找到", length(rds_file_list), "个RDS文件"))

# 设置循环范围
start_index <- 7
end_index <- 7

# 确保索引在有效范围内
start_index <- max(1, min(start_index, length(rds_file_list)))
end_index <- min(length(rds_file_list), max(end_index, start_index))

# 创建日志文件
log_dir <- file.path(target_dir, file_name, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, append = TRUE, split = TRUE)

cat("=== MR分析开始 ===\n")
cat("开始时间:", format(Sys.time()), "\n")
cat("分析范围: 第", start_index, "到第", end_index, "个暴露因素\n")

# 智能SNP筛选函数
smart_snp_filtering <- function(dat_clean, target_min = 10, target_max = 100) {
  # 定义p值阈值序列（从严格到宽松）
  pval_thresholds <- c(5e-8, 1e-8, 5e-9, 1e-9, 5e-10,  # 更严格
                       1e-7, 5e-7, 1e-6, 5e-6, 1e-5)   # 更宽松
  
  best_filtered <- NULL
  best_threshold <- 5e-8
  best_count <- 0
  
  # 首先尝试默认阈值
  dat_filtered <- dat_clean %>% filter(pval.exposure < 5e-8)
  current_count <- nrow(dat_filtered)
  
  cat("初始阈值5e-8得到的SNP数量:", current_count, "\n")
  
  # 如果初始数量在目标范围内，直接返回
  if (current_count >= target_min & current_count <= target_max) {
    return(list(data = dat_filtered, threshold = 5e-8, note = "使用默认阈值"))
  }
  
  # 如果数量过多，尝试更严格的阈值
  if (current_count > target_max) {
    cat("SNP数量过多，尝试更严格的阈值...\n")
    for (pval_thresh in pval_thresholds[pval_thresholds < 5e-8]) {
      temp_filtered <- dat_clean %>% filter(pval.exposure < pval_thresh)
      temp_count <- nrow(temp_filtered)
      cat("阈值", format(pval_thresh, scientific = TRUE), ": ", temp_count, "个SNP\n")
      
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
  
  # 如果数量过少，尝试更宽松的阈值（仅在MR-PRESSO需要时）
  if (current_count < target_min) {
    cat("SNP数量不足，尝试更宽松的阈值...\n")
    for (pval_thresh in pval_thresholds[pval_thresholds > 5e-8]) {
      temp_filtered <- dat_clean %>% filter(pval.exposure < pval_thresh)
      temp_count <- nrow(temp_filtered)
      cat("阈值", format(pval_thresh, scientific = TRUE), ": ", temp_count, "个SNP\n")
      
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
      
      # 如果已经达到上限，停止尝试
      if (temp_count > target_max * 2) break
    }
  }
  
  # 如果没有找到理想范围，使用最接近目标的值
  if (is.null(best_filtered)) {
    if (current_count > 0) {
      best_filtered <- dat_filtered
      best_threshold <- 5e-8
      best_count <- current_count
      cat("使用默认阈值，SNP数量:", best_count, "\n")
    } else {
      # 如果连宽松阈值都没有足够SNP，使用最宽松的阈值
      most_lenient <- max(pval_thresholds)
      best_filtered <- dat_clean %>% filter(pval.exposure < most_lenient)
      best_threshold <- most_lenient
      best_count <- nrow(best_filtered)
      cat("使用最宽松阈值", format(most_lenient, scientific = TRUE), 
          "，SNP数量:", best_count, "\n")
    }
  }
  
  note <- paste("优化后阈值:", format(best_threshold, scientific = TRUE),
                "SNP数量:", best_count)
  
  return(list(data = best_filtered, threshold = best_threshold, note = note))
}

# 优化的MR-PRESSO函数
run_optimized_presso <- function(dat_filtered, max_snps = 100) {
  snp_count <- nrow(dat_filtered)
  
  # 如果SNP数量不足，直接跳过
  if (snp_count < 10) {
    cat("SNP数量不足10个，跳过MR-PRESSO分析\n")
    return(NULL)
  }
  
  # 如果SNP数量过多，进行随机抽样
  if (snp_count > max_snps) {
    cat("SNP数量", snp_count, "大于", max_snps, "个，进行随机抽样...\n")
    set.seed(123) # 保证可重复性
    sampled_indices <- sample(1:snp_count, max_snps)
    dat_sampled <- dat_filtered[sampled_indices, ]
    cat("抽样后SNP数量:", nrow(dat_sampled), "\n")
  } else {
    dat_sampled <- dat_filtered
  }
  
  # 动态设置重抽样次数
  nb_dist <- ifelse(nrow(dat_sampled) > 50, 1000, 500)
  cat("设置NbDistribution为:", nb_dist, "\n")
  
  # 运行MR-PRESSO
  presso <- tryCatch({
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure", 
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      data = dat_sampled,
      SignifThreshold = 0.05,
      NbDistribution = nb_dist,
      seed = 12345
    )
  }, error = function(e) {
    cat("MR-PRESSO分析出错：", e$message, "\n")
    return(NULL)
  })
  
  return(presso)
}

# 主循环
for(i in start_index:end_index) {
  targetrds <- rds_file_list[i]
  exposure_name <- tools::file_path_sans_ext(basename(targetrds))
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("正在处理第", i, "个暴露因素:", exposure_name, "\n")
  cat("开始时间:", format(Sys.time()), "\n")
  
  result <- tryCatch({
    # 读取暴露数据
    exposure_dat <- readRDS(targetrds)
    cat("成功读取暴露数据，包含", nrow(exposure_dat), "个SNP\n")
    
    # 检查共有SNP
    common_snps <- intersect(exposure_dat$SNP, outcome_dat$SNP)
    cat("共有SNP数量:", length(common_snps), "\n")
    
    if(length(common_snps) < 3) {
      warning("共有SNP数量不足，跳过该暴露因素")
      next
    }
    
    # 数据协调
    harmonised_dat <- harmonise_data(
      exposure_dat = exposure_dat, 
      outcome_dat = outcome_dat, 
      action = 2
    )
    cat("协调后保留", nrow(harmonised_dat), "个SNP\n")
    
    # 数据清理
    dat_clean <- na.omit(harmonised_dat[,-c(18:20)])
    dat_clean <- subset(dat_clean, se.exposure > 1e-10 & se.outcome > 1e-10)
    
    # 智能SNP筛选
    cat("执行智能SNP筛选...\n")
    filtering_result <- smart_snp_filtering(dat_clean, target_min = 10, target_max = 100)
    dat_filtered <- filtering_result$data
    final_threshold <- filtering_result$threshold
    cat(filtering_result$note, "\n")
    
    if(nrow(dat_filtered) < 3) {
      stop("SNP数量不足，无法进行分析")
    }
    
    # 创建结果文件夹
    result_folder <- file.path(target_dir, file_name, basename(dirname(targetrds)))
    dir.create(result_folder, recursive = TRUE, showWarnings = FALSE)
    
    # 优化的MR-PRESSO分析
    cat("运行优化的MR-PRESSO分析...\n")
    presso_result <- run_optimized_presso(dat_filtered, max_snps = 100)
    
    if(!is.null(presso_result)) {
      cat("MR-PRESSO全局检验p值:", ifelse(!is.null(presso_result$`Global Test`$Pvalue), 
                                     presso_result$`Global Test`$Pvalue, "NA"), "\n")
    }
    
    # MR分析
    cat("运行MR分析...\n")
    mr_result <- mr(dat_filtered)
    mr_result <- generate_odds_ratios(mr_result)
    
    # 其他检验
    heterogeneity <- mr_heterogeneity(dat_filtered)
    pleiotropy <- mr_pleiotropy_test(dat_filtered)
    singlesnp <- mr_singlesnp(dat_filtered)
    leaveoneout <- mr_leaveoneout(dat_filtered)
    
    # 保存筛选信息
    filter_info <- data.frame(
      exposure = exposure_name,
      initial_snps = nrow(harmonised_dat),
      cleaned_snps = nrow(dat_clean),
      final_snps = nrow(dat_filtered),
      pval_threshold = final_threshold,
      mr_presso_performed = !is.null(presso_result),
      timestamp = Sys.time()
    )
    
    # 保存结果
    saveRDS(harmonised_dat, file.path(result_folder, "harmonised_data.RDS"))
    saveRDS(dat_filtered, file.path(result_folder, "filtered_data.RDS"))
    saveRDS(mr_result, file.path(result_folder, "mr_results.RDS"))
    saveRDS(filter_info, file.path(result_folder, "filter_info.RDS"))
    
    if(!is.null(presso_result)) {
      saveRDS(presso_result, file.path(result_folder, "mr_presso_results.RDS"))
    }
    
    # 创建详细的Excel报告
    wb <- createWorkbook()
    
    # 添加筛选信息
    addWorksheet(wb, "Filter Info")
    writeData(wb, "Filter Info", filter_info)
    
    addWorksheet(wb, "Main Results")
    writeData(wb, "Main Results", mr_result)
    
    addWorksheet(wb, "Single SNP")
    writeData(wb, "Single SNP", singlesnp)
    
    addWorksheet(wb, "Leave-One-Out")
    writeData(wb, "Leave-One-Out", leaveoneout)
    
    addWorksheet(wb, "Heterogeneity")
    writeData(wb, "Heterogeneity", heterogeneity)
    
    addWorksheet(wb, "Pleiotropy")
    writeData(wb, "Pleiotropy", pleiotropy)
    
    if(!is.null(presso_result)) {
      addWorksheet(wb, "MR-PRESSO")
      # 处理MR-PRESSO结果格式
      presso_summary <- data.frame(
        Test = c("Global Test", "Outlier Test"),
        Pvalue = c(ifelse(!is.null(presso_result$`Global Test`$Pvalue), 
                          presso_result$`Global Test`$Pvalue, NA),
                   ifelse(!is.null(presso_result$`Outlier Test`$Pvalue), 
                          presso_result$`Outlier Test`$Pvalue, NA))
      )
      writeData(wb, "MR-PRESSO", presso_summary)
    }
    
    saveWorkbook(wb, file.path(result_folder, "MR_Results.xlsx"), overwrite = TRUE)
    
    cat("第", i, "个暴露因素分析完成\n")
    cat("结果保存至:", result_folder, "\n")
    
    list(status = "success", exposure = exposure_name, snp_count = nrow(dat_filtered),
         threshold = final_threshold, presso_performed = !is.null(presso_result))
    
  }, error = function(e) {
    error_msg <- paste("错误:", e$message)
    cat(error_msg, "\n")
    cat("第", i, "个暴露因素处理失败:", exposure_name, "-", error_msg, "\n")
    list(status = "error", exposure = exposure_name, error = e$message)
  })
  
  cat("结束时间:", format(Sys.time()), "\n")
  cat("处理状态:", ifelse(result$status == "success", "成功", "失败"), "\n")
  
  if(result$status == "success") {
    cat("使用SNP数量:", result$snp_count, "\n")
    cat("最终p值阈值:", format(result$threshold, scientific = TRUE), "\n")
    cat("MR-PRESSO执行:", ifelse(result$presso_performed, "是", "否"), "\n")
  } else {
    cat("错误信息:", result$error, "\n")
  }
  
  cat(rep("=", 60), "\n\n", sep = "")
}

# 关闭日志文件
cat("\n=== MR分析结束 ===\n")
cat("结束时间:", format(Sys.time()), "\n")
sink()

# 恢复控制台输出
sink()

cat("分析完成！\n")
cat("日志文件保存至:", log_file, "\n")
cat("总共处理了", (end_index - start_index + 1), "个暴露因素\n")