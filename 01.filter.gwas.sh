vcftools --gzvcf merge.raw.vcf.gz \
--recode-INFO-all --remove-indels --recode --keep sample.txt \
--maf 0.05 --max-maf 0.99 --max-missing 0.9 \
--min-meanDP 3 --max-meanDP 30 \
--out gc.snp.flt --min-alleles 2 --max-alleles 2
##GATK的VariantsToTable将VCF文件中每个变体的指定字段提取到以制表符分隔的表中
/home/qmmou/anaconda3/envs/java-v8/bin/java -Xmx100g -jar /home/qmmou/software/GenomeAnalysisTK-3.5.jar \
                 -R /home/dpliu/index/PBH.v1.5.fasta \
                 -T VariantsToTable \
                 -V gc.snp.flt.recode.vcf \
                 -F CHROM -F POS -F REF -F ALT -F QUAL -F AC -GF GT \
                 -o gc.snp.results.table
##生成tfam文件（300：有多少个体写多少）（perl语言见附件）
perl /home/dpliu/meat_pop/code/VCFtable2plink.tped.tfam.pl gc.snp.results.table 789 gc.snp.tped gc.snp.tfam
##文件格式转化，转化成二进制
/home/hfliu/mambaforge-pypy3/envs/eqtl/bin/plink --tfile  gc.snp --recode12 --make-bed --out  gc300.snp --allow-extra-chr --chr-set 31
##对SNP进行质控
/home/hfliu/mambaforge-pypy3/envs/eqtl/bin/plink --bfile  gc300.snp --output-missing-genotype 0 --recode 12 transpose --out  gc300.snp --allow-extra-chr --chr-set 31 --geno 0.1 --maf 0.05
##构建亲缘关系矩阵
/home/qmmou/software/emmax-kin -v  -h -s -d 10 gc300.snp
##PCA分析
vcftools --vcf gc.snp.flt.recode.vcf --out sample.snp.plink --plink
/home/hfliu/mambaforge-pypy3/envs/eqtl/bin/plink --allow-extra-chr --chr-set 31 --file sample.snp.plink --out sample.plink.pca --pca
