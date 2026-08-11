# 加载必要的包
library(tidyverse)
library(broom)
setwd("C:/Users/HP/Desktop/locus/")
abundance_data<-read.table("1.txt",header = T,row.names = 1)
phenotype_data<-read.table("2.txt",header=T,row.names = 1)


# 数量性状模型分析函数（单个表型，所有菌属）
quantitative_analysis_single <- function(abundance, phenotype, min_positive_samples = 5) {
  pheno_name <- colnames(phenotype)
  pheno_values <- phenotype[, 1]
  cat(sprintf("\n=== 分析表型: %s ===\n", pheno_name))
    results <- data.frame(
    taxon = character(),
    n_positive = numeric(),
    p_value = numeric(),
    effect_size = numeric(),
    r_squared = numeric(),
    significant = logical(),
    stringsAsFactors = FALSE
  )
  
  # 对每个菌属进行分析
  for(taxon in colnames(abundance)) {
    taxon_abundance <- abundance[, taxon]
    positive_samples <- which(taxon_abundance > 0)
    n_positive <- length(positive_samples)
    
    # 只在有足够阳性样本时进行分析
    if(n_positive >= min_positive_samples) {
      # 对数转换
      log_abundance <- log(taxon_abundance[positive_samples] + 1e-10)
      pheno_subset <- pheno_values[positive_samples]
      
      # 线性回归
      model <- lm(pheno_subset ~ log_abundance)
      model_summary <- summary(model)
      
      # 提取结果
      coef_df <- tidy(model)
      p_val <- coef_df$p.value[2]
      effect <- coef_df$estimate[2]
      r2 <- model_summary$r.squared
      
      results <- rbind(results, data.frame(
        taxon = taxon,
        n_positive = n_positive,
        p_value = p_val,
        effect_size = effect,
        r_squared = r2,
        significant = p_val < 0.05
      ))
    } else {
      # 记录样本不足的菌属
      results <- rbind(results, data.frame(
        taxon = taxon,
        n_positive = n_positive,
        p_value = NA,
        effect_size = NA,
        r_squared = NA,
        significant = FALSE
      ))
    }
  }
  
  # 按p值排序
  results <- results[order(results$p_value, na.last = TRUE), ]
  return(results)
}

# 运行分析
cat("\n=== 开始数量性状模型分析 ===\n")
analysis_results <- quantitative_analysis_single(abundance_data, phenotype_data)

# 直接输出所有结果
cat("\n=== 所有菌属回归分析结果 ===\n")
print(analysis_results)
write.table(analysis_results,file = "Rongda-SKIN.FAT-weight.txt")
