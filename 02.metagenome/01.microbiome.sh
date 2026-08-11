#!/bin/sh
file=$(sed -n "1,789p" sample.txt)
for i in $file
do
/home/qmmou/anaconda3/envs/kaiju/bin/fastp -i /home/qmmou/metagenome/raw.data/${i}_1.fastq.gz -o ${i}_1.cutadapt.fastq.gz \
-I /home/qmmou/metagenome/raw.data/${i}_2.fastq.gz -O ${i}_2.cutadapt.fastq.gz \
-u 20 \
--length_required 30 \
--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
#rm -r /home/qmmou/metagenome2/fastq.gz/$d'_1.fastq.gz'
#rm -r /home/qmmou/metagenome2/fastq.gz/$d'_2.fastq.gz'
printf "Step 1: Finish removing adapters  at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> ${i}.calling.log
 
#######################step2.clean-host###################
# 运行比对程序
bwa mem -t 16 /home/dpliu/index/PBH.v1.5.fasta ./${i}_1.cutadapt.fastq.gz ./${i}_2.cutadapt.fastq.gz > ${i}_mapped_and_unmapped.sam
printf "Step 2: Finish mapping to PBH.v1.5  at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $i'.calling.log'
#将sam文件转化为bam文件
samtools view -@ 10 -bS ${i}_mapped_and_unmapped.sam > ${i}_mapped_and_unmapped.bam
printf "Step 3: Finish transforming the SAM to BAM  at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $i'.calling.log'
#去除比对到宿主基因组两条链的基因组
/usr/bin/samtools view -b -@ 10 -f 12 -F 256 ${i}_mapped_and_unmapped.bam > ${i}_bothEndsmapped.bam
printf "Step 4: Finish removing the mapped reads  at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $i'.calling.log'
#排序bam文件
/usr/bin/samtools sort -n ${i}_bothEndsmapped.bam -o ${i}_bothEndsmapped_sorted.bam
printf "Step 5: Finish sorting the BAM  at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $i'.calling.log'
#将bam文件转化成fastq文件
/home/qmmou/anaconda3/envs/kaiju/bin/bedtools bamtofastq -i ${i}_bothEndsmapped_sorted.bam -fq ${i}_r1.fastq -fq2 ${i}_r2.fastq
#压缩fastq文件
gzip ${i}_r1.fastq
gzip ${i}_r2.fastq
#删除中间文件
rm -r ${i}_mapped_and_unmapped.sam
rm -r ${i}_mapped_and_unmapped.bam
rm -r ${i}_bothEndsmapped.bam
rm -r ${i}_bothEndsmapped_sorted.bam
rm -r ${i}_1.cutadapt.fastq.gz ${i}_2.cutadapt.fastq.gz
printf "Step 2: Finish cleaning the host  at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $i'.calling.log'
done

#######################step3.kraken2###################
file=$(sed -n "1,789p" sample.txt)
for i in $file
do
/home/qmmou/anaconda3/envs/kaiju/bin/kraken2 --db /home/qmmou/index/kraken2/standard \
--gzip-compressed \
./${i}_r1.fastq.gz ./${i}_r2.fastq.gz \
--report ${i}_kraken2-metaphlan.mpa.report \
--use-mpa-style \
--output ${i}_kraken2-metaphlan.mpa.txt \
--threads 60
rm -r ${i}_kraken2-metaphlan.mpa.txt
done
