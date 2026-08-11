##When you caculate Fst, you just need one step:
vcftools --vcf /home/qmmou/merge.vcf/gc.snp.flt.recode.vcf --weir-fst-pop pop1.list --weir-fst-pop pop2.list --fst-window-size 10000 --fst-window-step 5000 --out pop1_VS_pop2
##or you can caculate Fst for each SNP with one step:
vcftools --vcf /home/qmmou/merge.vcf/gc.snp.flt.recode.vcf --weir-fst-pop pop1.list --weir-fst-pop pop2.list --out pop1_VS_pop2.snp
##Then you can use script drawing.
perl /home/qmmou/code/fst.result.sub.piture.pl chr1.fst
##R
R
mydata=read.table("R1.R7.Fst.txt.picture", header=T)
colnames(mydata)=c("CHR","CHR","BP","P") 
mydata=mydata[,c(1,3,4)]
#png(file="myplot.png",width=15,height=5)
#source("/home/zkzhou/demo_data/Fst/change.Ylabname.FST.source_v2.r")
source("/home/qmmou/code/change.Ylabname.FST.source_v3.r")
GAPIT.Manhattan(mydata,name.of.trait=paste("R1_vs_R7",sep = ""))
#dev.off()
q("no")
