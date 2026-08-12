file=$(sed -n "1,789p" sample.txt)
for d in $file
do
#BWA align reads1
bwa aln -t 64 /home/dpliu/index/PBH.v1.5.fasta \
               /home/dataset/DNA-seq/Z2_Z5_LC_MAL/$d'_clean_1.fq.gz' >$d'_1.sai'
#BWA align reads2
bwa aln -t 64 /home/dpliu/index/PBH.v1.5.fasta \
               /home/dataset/DNA-seq/Z2_Z5_LC_MAL/$d'_clean_2.fq.gz' >$d'_2.sai'
#BWA generate sam file
bwa sampe /home/dpliu/index/PBH.v1.5.fasta \
                $d'_1.sai' \
                $d'_2.sai' \
               /home/dataset/DNA-seq/Z2_Z5_LC_MAL/$d'_clean_1.fq.gz' \
               /home/dataset/DNA-seq/Z2_Z5_LC_MAL/$d'_clean_2.fq.gz' \
		>$d'.sam'
#Convert sam to bam
samtools view -bS $d'.sam' > $d'.bam'
#Bam file sort
samtools sort -o $d'.sorted.bam' -T $d'.sorasted.tmp' -@ 8 -O bam $d'.bam'
#Basic statistic for bam file
samtools flagstat $d'.sorted.bam' >$d'.flagstat.txt'
# REMOVE process files
rm $d'_1.sai'
rm $d'_2.sai'
rm $d'.sam'
rm $d'.bam'
# Step 1: Remove unmapped and multihit reads
samtools view -h -F 4 -b $d'.sorted.bam' >$d'.mapped.sorted.bam'
samtools view -bq 1 $d'.mapped.sorted.bam' > $d'.sorted.uniqe.bam'
samtools index $d'.sorted.uniqe.bam'
printf "Step 1: remove unmapped and multihit reads finished at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $d'.calling.log'
# Step 2: Add bam reads group
/home/hfliu/mambaforge-pypy3/envs/java8/bin/java -jar -Xmx100g /home/dpliu/software/AddOrReplaceReadGroups.jar \
                I=$d'.sorted.uniqe.bam' \
                O=$d'.sorted.uniqe.rg.bam' \
                LB=$d \
		PL=illumina \
                PU=IAS \
                SM=$d \
                VALIDATION_STRINGENCY=SILENT
samtools index $d'.sorted.uniqe.rg.bam'
printf "Step 2: add bam reads group finished at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $d'.calling.log'
# Step 3: Mark Duplicates
/home/hfliu/mambaforge-pypy3/envs/java8/bin/java -jar -Xmx4g /home/dpliu/software/MarkDuplicates.jar \
                I=$d'.sorted.uniqe.rg.bam' \
                O=$d'.sorted.uniqe.rg.dedup.bam' \
		M=$d'.dedup.metrics' \
                VALIDATION_STRINGENCY=SILENT
samtools index $d'.sorted.uniqe.rg.dedup.bam'
printf "Step 3: Mark Duplicates finished at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $d'.calling.log'
# Step 4: Create INDEL position(realigne target)
/home/hfliu/mambaforge-pypy3/envs/java8/bin/java -Xmx100g -jar /home/dpliu/software/GenomeAnalysisTK-3.5.jar \
     -R /home/dpliu/index/PBH.v1.5.fasta \
     -T RealignerTargetCreator \
     -o $d'.realn.intervals' \
     -I $d'.sorted.uniqe.rg.dedup.bam'
printf "Step 4: Create INDEL position(realigne target) finished at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $d'.calling.log'
## Step 5: INDEL realigner
/home/hfliu/mambaforge-pypy3/envs/java8/bin/java -Xmx100g -jar /home/dpliu/software/GenomeAnalysisTK-3.5.jar \
     -R /home/dpliu/index/PBH.v1.5.fasta \
     -T IndelRealigner \
     -targetIntervals $d'.realn.intervals' \
     -o $d'.sorted.uniqe.rg.dedup.realn.bam' \
     -I $d'.sorted.uniqe.rg.dedup.bam'
printf "Step 5: INDEL realigner finished at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $d'.calling.log'
## Step 6: calling
/home/hfliu/mambaforge-pypy3/envs/java8/bin/java -Xmx200g -jar /home/dpliu/software/GenomeAnalysisTK-3.5.jar \
     -R /home/dpliu/index/PBH.v1.5.fasta \
     -T HaplotypeCaller -nct 32 \
     -I $d'.sorted.uniqe.rg.dedup.realn.bam' \
     -o $d'.gvcf' \
     --genotyping_mode DISCOVERY \
     -stand_emit_conf 30 \
     -stand_call_conf 30 \
     -ERC GVCF \
     -variant_index_type LINEAR \
     -variant_index_parameter 128000
## Step 7: file compression
bgzip $d'.gvcf'
tabix -p vcf $d'.gvcf.gz'
rm $d'.dedup.metrics'
rm $d'.sorted.bam'
rm $d'.mapped.sorted.bam'
rm $d'.realn.intervals'
rm $d'.sorted.uniqe.bam'
rm $d'.sorted.uniqe.bam.bai'
rm $d'.sorted.uniqe.rg.bam'
rm $d'.sorted.uniqe.rg.bam.bai'
rm $d'.sorted.uniqe.rg.dedup.bam'
rm $d'.sorted.uniqe.rg.dedup.bam.bai'
done
