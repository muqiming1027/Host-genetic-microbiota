##### my locuszoom for GWAS causitive gene plot ######

##### Write by zhkzhou@126.com 2014-07-16
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("image.plot")
library(spectralGP)
#install.packages("spectralGP")
`GAPIT.Pruning` <-
function(values,DPP=100){
if(length(values)<=DPP)return(c(1:length(values)))
values=sqrt(values)  #This shift the weight a little bit to the low building.
theMin=min(values)
theMax=max(values)
range=theMax-theMin
interval=range/DPP
ladder=round(values/interval)
ladder2=c(ladder[-1],0)
keep=ladder-ladder2
index=which(keep>0)
return(index)
}#end of GAPIT.Pruning

`GAPIT.Manhattan` <-
function(GI.MP = NULL, name.of.trait = "Trait",plot.type2 = "Chromosomewise",
DPP=10000,cutOff=0.01,band=5,seqQTN=NULL){
if(is.null(GI.MP)) return
borrowSlot=5     ######################## column number +1 ;                         
GI.MP[,borrowSlot]=0 #Inicial as 0
if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1

#Eeep QTN with NA p values (set it to 1)
index=which(GI.MP[,borrowSlot]==1  & is.na(GI.MP[,3]))
GI.MP[index,3]=1
GI.MP=matrix(as.numeric(as.matrix(GI.MP) ) ,nrow(GI.MP),ncol(GI.MP))

#Remove all SNPs that do not have a choromosome, bp position and p value(NA)
GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
GI.MP <- GI.MP[!is.na(GI.MP[,3]),]
GI.MP <- GI.MP[!is.na(GI.MP[,4]),]
#Remove all SNPs that have P values between 0 and 1 (not na etc)
GI.MP <- GI.MP[GI.MP[,3]>0,]
GI.MP <- GI.MP[GI.MP[,3]<=1,]
GI.MP <- GI.MP[GI.MP[,1]!=0,]

numMarker=nrow(GI.MP)
bonferroniCutOff=-log10(cutOff/numMarker)

#Replace P the -log10 of the P-values
GI.MP[,3] <-  -log10(GI.MP[,3])

y.lim <- as.integer(ceiling(max(GI.MP[,3])))
y.lim = y.lim+2
print("The max -logP vlaue is")
print(y.lim)

chm.to.analyze <- unique(GI.MP[,1])
chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
numCHR= length(chm.to.analyze)

#Chromosomewise plot
if(plot.type2 == "Chromosomewise")
{
  pdf(paste(name.of.trait,".Chromosomewise.pdf" ,sep = ""), width = 10)
  par(mar = c(1,5,17,3), lab = c(8,5,7))
  for(i in 1:numCHR)
  {

    subset=GI.MP[GI.MP[,1]==chm.to.analyze[i],]
  	y.lim <- as.integer(ceiling(max(subset[,3])))  #set upper for each chr
  	if(length(subset)>4){
      x <- as.numeric(subset[,2])/10^(6)
      y <- as.numeric(subset[,3])
      z <- as.numeric(subset[,4])
    }else{
      x <- as.numeric(subset[2])/10^(6)
      y <- as.numeric(subset[3])
      z <- as.numeric(subset[4])      
    }

    order=order(y,decreasing = TRUE)
    y=y[order]
    x=x[order]
    z=z[order]

    index=GAPIT.Pruning(y,DPP)
   	x=x[index]
  	y=y[index]
  	z=z[index]
  	#print(z)
    y.lim = y.lim+1
    par(bty="l", lwd=1.5)
    par(mar = c(15,6,5,2),mgp=c(3.2, 1.2, 0))  ## mgp: adjust the distnace of label and axis  
    plot(y~x,type="p", cex.lab=1.5, cex.axis=1.5,
    ylim=c(0,y.lim), xlim = c(min(x), max(x)),    
    #col="black", lwd=0.5,bg= ifelse(x >= 40, "red", ifelse(x <40 & x >= 30, "blue", ifelse((x >= 20 & x < 30), "chocolate", ifelse((x >= 10 & x < 20), "darkgreen", "darkorchid4")))),
    col="black", lwd=0.5, bg= ifelse(z >= 0.8, "red", ifelse(z <0.8 & z >= 0.6, "blue", ifelse((z >= 0.4 & z < 0.6), "orange", ifelse((z >= 0.2 & z < 0.4), "darkgreen", "darkorchid4")))),
    pch= ifelse(z >= 1, 23, 21),
    lwd.ticks=1.5,
    xlab =paste("Chromosome ",chm.to.analyze[i]," (Mb)",sep="") ,
    ylab=expression(-log[10]~italic(P)),
    xaxs="i", yaxs="i",  cex=1.5 ,
    main =paste("GWAS on ",name.of.trait,sep=""),cex.main=2.5, font.main = 1)
    #abline(v=37.375,col="dimgray",lty=2, lwd=1.5)
    #abline(v=38.66,col="dimgray",lty=2, lwd=1.5)
    legend("topright", c(c(expression(italic(R)^2))), cex=1.0,bty="n")
      par( mar=c(27,4,3.5,1),lwd=0.5,cex=0.7 )
      #legend("top", c(c(expression(italic(R)^2))), cex=1.0,bty="n")
      image_plot( zlim=c(0.0,1,0.2), lwd=0.5 ,col=c("darkorchid4","darkgreen","orange","blue","red") , legend.only=TRUE, axis.args = list(cex.axis = 0.5), axes=FALSE )
  }

  #legend(0,0, c(c(expression(italic(R)^2))), cex=1.0,bty="n")
  print("Manhattan-Plot.Chromosomewise finished!")
  dev.off()
}
}
