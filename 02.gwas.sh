 #!/bin/bash
for (( d=1; d<=5312; d++ ))
do
m=$(sed -n "${d}p" ./pheno.list.txt)
n=$[$d+2]
cut -f1-2,$n ./data-gwas.txt > $m'.txt'
/home/qmmou/software/emmax -v -d 10 -t ./gc300.snp  -p $m'.txt' -k ./gc300.snp.hIBS.kinf -c ./sex.txt -o 'gwas.'$m
cut -d" " -f1-4 ./gc300.snp.tped > gc300.snp.headcol
paste gc300.snp.headcol 'gwas.'$m'.ps' >'gwas.'$m'.result.txt'
awk '{$3=null;$5=null;$6=null;print}' 'gwas.'$m'.result.txt' >'gwas.'$m'.assoc.txt'
sed -i '1i\CHR SNPID BP Pi P' 'gwas.'$m'.assoc.txt'
perl /home/qmmou/code/emmaxfilter.pl 'gwas.'$m'.assoc.txt' 'gwas.'$m'.assoc.plot.txt'
#awk '{print $1,$3,$5}' 'gwas.'$m'.assoc.txt' > 'gwas.'$m'.assoc.plot.txt'
#rm 'gwas.'$m'.result.txt'
#rm 'gwas.'$m'.ps'
#rm 'gwas.'$m'.assoc.txt'
echo "source(\"/home/qmmou/code/MANHATTAN_QQ_png.r\")
mydata=read.table(\"gwas."$m".assoc.plot.txt\", header=TRUE)
GAPIT.Manhattan(mydata,name.of.trait=\""$m"\")
mydata=mydata[order(as.numeric(mydata[,3]),decreasing = FALSE),]
myoutdata=mydata[1:2001,]
myoutdata[,3] <-  -log10(myoutdata[,3])
write.table(myoutdata,file=\"gwas."$m".assoc.t2000.txt\", col.names=T, row.names=F, quote=F,sep=\" \")" >plot."$m".r
perl -e 'system ("R <plot.'$m'.r --vanilla")' &
done
##GWAS分析结束，获得manhattan图和qq图
