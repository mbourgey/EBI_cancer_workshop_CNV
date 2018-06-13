## set environement
cd $HOME/ebicancerworkshop2018/CNV/SNParray

R

library(ASCAT)

ascat.bc = ascat.loadData("tumor2.LRR.tsv","tumor2.BAF.tsv","normal2.LRR.tsv","normal2.BAF.tsv")

ascat.plotRawData(ascat.bc)

ascat.seg = ascat.aspcf(ascat.bc)

ascat.plotSegmentedData(ascat.seg)

ascat.output = ascat.runAscat(ascat.seg)

params.estimate=data.frame(Sample=names(ascat.output$aberrantcellfraction),Aberrant_cell_fraction=round(ascat.output$aberrantcellfraction,2),Ploidy=round(ascat.output$ploidy,2))
write.table(params.estimate,"sample.Param_estimate.tsv",sep="\t",quote=F,col.names=T,row.names=F)

CNA=rep(".",dim(ascat.output$segments)[1])
CNA[rowSums(ascat.output$segments[,5:6]) > round(ascat.output$ploidy)]="DUP"
CNA[rowSums(ascat.output$segments[,5:6]) < round(ascat.output$ploidy)]="DEL"
output.table=data.frame(ascat.output$segments,CNA=CNA)
write.table(output.table[output.table$CNA != ".",],"sample_CNVcalls.tsv",quote=F,sep="\t",col.names=T,row.names=F)
q(save="yes")
