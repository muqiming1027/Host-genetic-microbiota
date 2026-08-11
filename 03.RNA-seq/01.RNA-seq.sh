#进行比对
file=$(sed -n "1,5p" sample.txt)
for d in $file
do
hisat2 -p 64 -x /home/qmmou/index/PBH.v1.5 \
-1 $d'_R1_clean.fastq.gz' -2 $d'_R2_clean.fastq.gz' \
-S $d'.sam'
#排序sam文件
samtools sort -@ 4 -o $d'_out.bam' $d'.sam'
#将sam文件转换成bam文件，定量分析
samtools view -h $d'_out.bam' | htseq-count  -r pos -t gene -i gene \
- /home/dpliu/index/ref_IASCAAS_PekingDuck_PBH1.5_top_level12.gff3 \
--stranded=no >$d'_out.readscount.txt'
rm -r $d'.sam'
done
