# 合并重复的库加载（避免重复加载）
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

##TwoSampleMR包
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR",force = TRUE)

##MRInstruments包
#remotes::install_github("MRCIEU/MRInstruments")

#library(devtools)
#install_github("qingyuanzhao/mr.raps")
####---
###    准备暴露因素数据集----------------------------------------
####---
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(ieugwasr)
library(MRInstruments)
library(mr.raps)
library(MendelianRandomization)
library(vroom)
rm(list = ls())

## Fenn 通过表格下载数据
## GWAS 通过https://gwas.mrcieu.ac.uk/ 查询，然后下载vcf格式文件

target_outcome <- c("TB") ####根据需求替换

# 设置固定工作目录路径（根据你的需求修改此处）
target_dir <- "E:/20250303 TB NHANES/20251209 IJS revise/TSMR"  # 注意使用正斜杠或双反斜杠
setwd(target_dir)

# 自动获取ZIP文件并处理路径
# 然后查找文件
file_list <- list.files(path = "data/outcome", pattern = "\\.gz$", full.names = TRUE)

zip_file <- file_list[1]          # 取当前目录第一个ZIP文件

stopifnot("未找到gz文件" = !is.na(zip_file))          # 确保文件存在

library(tools)
zip_name <- file_path_sans_ext(zip_file)               # 自动提取文件夹名称

print(zip_name)

file_name <- basename(zip_name)
print(file_name)

# 1. 创建新文件夹
dir.create(file_name)                                   # 新建空文件夹

# 结局因素的处理
# 安装必要包
#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

# 读取 VCF 文件
library(VariantAnnotation)
outcome_dat <- readVcf(zip_file, genome = "hg19")

# 读取 tsv.gz 文件
library(readr)
#outcome_dat <- read_tsv(zip_file)

# 2. 设置新工作目录
setwd(file_name)                                        # 切换到解压文件夹
message("当前工作目录已设置为：", getwd())

#head(outcome_dat,10)

library(dplyr)

library(gwasglue)
outcome_raw <- gwasglue::gwasvcf_to_TwoSampleMR(outcome_dat, type = "outcome")  # 自动转换为TwoSampleMR格式

# 重命名核心列
outcome_dat <- outcome_raw

outcome_dat$id.outcome <- zip_name

saveRDS(outcome_dat,"outcome_dat.RDS")

###    proteomic_qtls的分析开始了----------------------------------------
####---

### outcome 处理好了
rm(list = ls())
outcome_dat <- readRDS("outcome_dat.RDS")
#data("drug_interactions")
data("proteomic_qtls")

# 所有分析
exposure_dat <- proteomic_qtls

print(exposure_dat)

colnames(exposure_dat)

# 3. 检查每个基因-组织的SNP数量
snp_counts <- exposure_dat %>%
  group_by(analyte, SNP
  ) %>%
  summarise(N_SNPs = n())

# 在调用format_data前生成唯一ID
exposure_dat$unique_id <- paste(exposure_dat$analyte#,
                                #exposure_dat$gene,
                                #exposure_dat$location, 
                                #exposure_dat$annotation, 
                                #sep = "_"
                                )

exposure_dat <- format_data(
  exposure_dat,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  id_col = "unique_id",  # ID列 设置为 phenotype
  eaf_col = "eaf",     # 假设原始数据包含效应等位基因频率
  gene_col = "gene"#,
  #phenotype_col = "phenotype",
  # samplesize_col = "n"  # 样本量列（如有）
)

head(exposure_dat)

head(exposure_dat[, c("SNP", "id.exposure", "beta.exposure")])

saveRDS(exposure_dat,"exposure_dat protein.RDS")

### harmonise ---
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
saveRDS(dat,"dat protein.RDS")

library(MRPRESSO)
# 检查有效SNP数量
cat("Available SNPs:", nrow(dat), "\n")
if (nrow(dat) < 10) {
  warning("SNP数量不足，MR-PRESSO无法运行！建议放宽筛选条件。")
} else {
  # 运行MR-PRESSO的代码
  # MR-PRESSO全局检验（需安装MRPRESSO包）
  presso <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                                BetaExposure = "beta.exposure",
                                SdOutcome = "se.outcome",
                                SdExposure = "se.exposure",
                                data = dat,
                                SignifThreshold = 0.05,  # 默认0.05，放宽至0.1
                                NbDistribution = 1000)# 默认1000
  
  cat("MR-PRESSO global test p-value:", presso$`Global Test`$Pvalue, "\n")
}

MR_PRESSO <- as.data.frame(presso)

### MR 分析--
# MR 分析
result <- mr(dat#, method_list = mr_method_list[-16,]$obj
)#去掉mr raps

# 计算OR值
result <- generate_odds_ratios(result)

head(result)

# 异质性检验
heterogeneity <- mr_heterogeneity(dat)
#若Q检验P < 0.05：存在显著异质性，提示潜在水平多效性。

# 水平多效性检验
pleiotropy <- mr_pleiotropy_test(dat)

# 散点图
#p1 <- mr_scatter_plot(result, dat)

# 森林图
singlesnp <- mr_singlesnp(dat)
#p2 <- mr_forest_plot(result_single)

# 留一图
leaveoneout <- mr_leaveoneout(dat)
#p3 <- mr_leaveoneout_plot(result_loo)

# 漏斗图
#result_single <- mr_singlesnp(dat)
#p4 <- mr_funnel_plot(result_single)

# 安装并加载所需包（如果尚未安装）
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)

# 创建 Excel 工作簿对象
wb <- createWorkbook()

# 添加第一个工作表：主要结果（result）
addWorksheet(wb, sheetName = "Main Results")
writeData(wb, sheet = 1, x = result, startCol = 1, startRow = 1)

# 添加第二个工作表：单SNP分析（result_single）
addWorksheet(wb, sheetName = "Single SNP")
writeData(wb, sheet = 2, x = singlesnp, startCol = 1, startRow = 1)

# 添加第三个工作表：留一法分析（result_loo）
addWorksheet(wb, sheetName = "Leave-One-Out")
writeData(wb, sheet = 3, x = leaveoneout, startCol = 1, startRow = 1)

# 添加第4个工作表：MRPRESSO
addWorksheet(wb, sheetName = "MRPRESSO")
writeData(wb, sheet = 4, x = MR_PRESSO, startCol = 1, startRow = 1)

# 添加第5个工作表：heterogeneity
addWorksheet(wb, sheetName = "heterogeneity")
writeData(wb, sheet = 5, x = heterogeneity, startCol = 1, startRow = 1)

# 添加第6个工作表：pleiotropy
addWorksheet(wb, sheetName = "pleiotropy")
writeData(wb, sheet = 6, x = pleiotropy, startCol = 1, startRow = 1)

# 保存Excel文件到当前工作目录
saveWorkbook(wb, 
             file = "MR_protein_Results.xlsx",
             overwrite = TRUE)

# 确认输出完成
message("Results saved to: ", file.path(getwd(), "MR_protein_Results.xlsx"))