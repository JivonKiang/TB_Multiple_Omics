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
#library(ieugwasr)
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
file_list <- list.files(path = "data/exposure", pattern = "\\.gz$", full.names = TRUE)

zip_file <- file_list[8]          # 取当前目录第一个ZIP文件

stopifnot("未找到gz文件" = !is.na(zip_file))          # 确保文件存在

library(tools)
zip_name <- file_path_sans_ext(zip_file)               # 自动提取文件夹名称

print(zip_name)

file_name <- c("data/exposure")
print(file_name)

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

# 1. 创建新文件夹
dir.create(zip_name)                                   # 新建空文件夹

# 2. 设置新工作目录
setwd(zip_name)                                        # 切换到解压文件夹
message("当前工作目录已设置为：", getwd())

#head(outcome_dat,10)

library(dplyr)

library(gwasglue)
outcome_raw <- gwasglue::gwasvcf_to_TwoSampleMR(outcome_dat, type = "exposure")  # 自动转换为TwoSampleMR格式

exposure_dat <- outcome_raw

colnames(exposure_dat)

# 在调用format_data前生成唯一ID
exposure_dat <- format_data(
  exposure_dat,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  id_col = "id.exposure",  # ID列 设置为 phenotype
  phenotype_col = "exposure",
  eaf_col = "eaf.exposure"#,     # 假设原始数据包含效应等位基因频率
  #gene_col = "gene"#,
  
  # samplesize_col = "n"  # 样本量列（如有）
)

# 重命名核心列
exposure_dat$id.exposure <- basename(zip_name)

saveRDS(exposure_dat,"exposure_dat.RDS")

