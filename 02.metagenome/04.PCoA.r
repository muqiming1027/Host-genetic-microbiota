#!/usr/bin/env Rscript

# 加载 vegan 包（如果未安装，请先执行 install.packages("vegan")）
library(vegan)

# 读取数据：第一列为行名（Contigs），其余为样品列
data <- read.table("result_reordered.txt", header = TRUE, row.names = 1, 
                   check.names = FALSE, stringsAsFactors = FALSE)

# 转置：使行为样品、列为基因（contigs）
data_t <- t(data)

# 计算 Bray-Curtis 距离矩阵
dist_bc <- vegdist(data_t, method = "bray")

# 执行 PCoA（cmdscale），返回所有可能的坐标轴
pcoa <- cmdscale(dist_bc, eig = TRUE, k = nrow(data_t) - 1)

# 提取样品坐标（行名即为样品名）
coords <- pcoa$points

# 保存到 PCoA.txt（制表符分隔，保留行名）
write.table(coords, file = "PCoA.txt", sep = "\t", quote = FALSE, col.names = NA)

# 可选：打印解释方差比例到屏幕
eig <- pcoa$eig
var_explained <- eig / sum(eig)
cat("Variance explained by each axis:\n")
print(var_explained)
