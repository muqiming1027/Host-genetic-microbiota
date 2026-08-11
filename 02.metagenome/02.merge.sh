file=$(cat sample.txt)
for i in $file
do
awk '{print $1"#  "$2}' ${i}_kraken2-metaphlan.mpa.report > ${i}.report && mv ${i}.report ${i}_kraken2.report
done
printf "Finsh add_#"


file=$(cat sample.txt)
for i in $file
do
#awk '{print $1"#  "$2}' ${i}_kraken2.report > ${i}.report && mv ${i}.report ${i}_kraken2.report
awk '{print $1}' ${i}_kraken2.report >> microbes1.txt
done
sort microbes1.txt|uniq > microbes.txt
rm microbes1.txt 
printf "Finish grep microbes"



touch merge.txt
file=$(cat microbes.txt)
samfile=$(cat sample.txt)
for i in $file
do
#awk '$1=="kkkkkkkkkkkkkkkkkkkkkkkkkkk"' *_kraken2.report >> 1.txt || echo '-' >> 1.txt
#cat sed2.sh|sed "s/kkkkkkkkkkkkkkkkkkkkkkkkkkk/${i}/g" > sed3.sh
#sh sed3.sh
for d in $samfile
do
grep -w ${i} ${d}_kraken2.report >> 1.txt || echo '-' >> 1.txt
#rm sed3.sh
done
awk '{print $NF}' 1.txt > 2.txt
paste merge.txt 2.txt > merge1.txt && mv merge1.txt merge.txt
rm 1.txt 2.txt
done
printf "<<<< Finish!!! >>>>"
 
