library(clusterProfiler)
library(stringr)
library(openxlsx)
library(org.Hs.eg.db)
##01.准备文件
remove(list = ls()) #清除 Global Environment
getwd()  #查看当前工作路径
setwd("C:/Users/HP/Desktop/locus")  #设置需要的工作路径
list.files()  #查看当前工作目录下的文件
data = read.xlsx("DEG_edgeR_symbol.xlsx",sheet= "DEG_edgeR_symbol",sep=',') #导入数据
head(data)
##02.筛选差异基因
#数据处理-差异基因筛选#
vector = abs(data$log2FoldChange) > 1 & data$FDR < 0.01 & data$gene_name !="" ##abs绝对值;通常logFC> 1和PValue< 0.05条件进行筛选；data$gene_name != ""表示gene_name不为空白
#data$gene_name<-str_to_title(data$gene_name)#用stringr将基因名称的第一个字母大写（小鼠首字母为大写）
data_sgni= data[vector,]#筛选差异基因
head(data_sgni)
##03.差异基因ID转换
#已有OrgDb的常见物种富集分析#
#BiocManager::install("org.Hs.eg.db") 
library(org.Hs.eg.db)
#基因ID转换#
keytypes(org.Hs.eg.db) #查看所有可转化类型
entrezid_all = mapIds(x = org.Hs.eg.db,  #id转换的比对基因组的物种，这处示例为人
                      keys = data_sgni$gene_name, #将输入的gene_name列进行数据转换
                      keytype = "SYMBOL", #输入数据的类型
                      column = "ENTREZID") #输出数据的类型
entrezid_all = na.omit(entrezid_all)  #na省略entrezid_all中不是一一对应的数据情况
entrezid_all = data.frame(entrezid_all) #设置数据格式；数据格式还可选择vector，matrix，array，list，data.frame
head(entrezid_all)
##04.KEGG和GO富集分析
#GO富集分析#
go_enrich = enrichGO(gene = entrezid_all[,1], #表示前景基因，即待富集的基因列表;[,1]表示对entrezid_all数据集的第1列进行处理
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", #输入数据的类型
                     ont = "ALL", #可以输入CC\MF\BP\ALL
                     pAdjustMethod = "fdr", # 指定多重假设检验矫正的方法，选项包含 "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                     pvalueCutoff = 0.05,  #指定 p 值阈值（指定 1 以输出全部），默认为0.05
                     qvalueCutoff = 0.05, #指定 q 值阈值（指定 1 以输出全部），默认0.2
                     readable = T) # 是否将gene ID转换到 gene symbol
go_enrich  = data.frame(go_enrich) #将GO_enrich导成数据框格式
write.csv(go_enrich,'C:/Users/HP/Desktop/locus/go_enrich.csv') #数据导出#
#KEGG富集分析#
KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], #即待富集的基因列表
                         keyType = "kegg",
                         organism= "human",  #hsa，可根据你自己要研究的物种更改，可在https://www.kegg.jp/brite/br08611中寻找
                         pAdjustMethod = "fdr", # 指定多重假设检验矫正的方法，选项包含 "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                         pvalueCutoff = 1,  #指定 p 值阈值（指定 1 以输出全部），默认为0.05
                         qvalueCutoff =1 #指定 q 值阈值（指定 1 以输出全部），默认0.2
                         )  
KEGG_enrich  = data.frame(KEGG_enrich)
write.csv(KEGG_enrich,'C:/Users/HP/Desktop/locus/KEGG_enrich.csv') #数据导出
##作图
#install.packages("ggplot2")
pdf(file="myplot.pdf",bg="white",width=10,height=12)
library(ggplot2)
go_enrich = read.xlsx("GO.xlsx",sheet= "Sheet1",sep=',') 
#go_enrich$term <- paste(go_enrich$ID, go_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
#go_enrich$term <- factor(go_enrich$term, levels = go_enrich$term,ordered = T)
go_enrich$term <- factor(go_enrich$Description, levels = go_enrich$Description,ordered = T)
#纵向柱状图#
ggplot(go_enrich,
       aes(x=Description,y=Count, fill=ONTOLOGY)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Count") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme_bw()+
  theme(axis.text=element_text(size=15,face = "bold"),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))
dev.off()
#气泡图#
pdf(file="myplot.pdf",bg="white",width=12,height=12)
library(ggplot2)
KEGG_enrich = read.xlsx("KEGG.xlsx",sheet= "Sheet1",sep=',') 
ggplot(KEGG_enrich,
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  #facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Ratio",y="KEGG term",title="KEGG Enrichment")+
  theme_bw()
dev.off()
