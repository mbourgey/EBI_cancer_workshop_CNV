
## set environement

#launch docker
docker run --privileged -v /tmp:/tmp --network host \
  -it -w $PWD -v $HOME:$HOME --user $UID:$GROUPS \
  -v /etc/group:/etc/group  -v /etc/passwd:/etc/passwd \
  -e DISPLAY=$DISPLAY -v /etc/fonts/:/etc/fonts/  c3genomics/genpipes:0.8


cd $HOME/ebicancerworkshop2020/CNV/SNParray


module load mugqic/R_Bioconductor/3.6.0_3.9




R


library(ASCAT)


ascat.bc = ascat.loadData(
  "C0053/tumor/tumor2.LRR.tsv", 
  "C0053/tumor/tumor2.BAF.tsv", 
  "C0053/normal/normal2.LRR.tsv",
  "C0053/normal/normal2.BAF.tsv")


ascat.plotRawData(ascat.bc)


ascat.seg = ascat.aspcf(ascat.bc)


ascat.plotSegmentedData(ascat.seg)


ascat.output = ascat.runAscat(ascat.seg)


params.estimate=data.frame(
   Sample=names(ascat.output$aberrantcellfraction),
   Aberrant_cell_fraction=round(ascat.output$aberrantcellfraction,2),
   Ploidy=round(ascat.output$ploidy,2)
)

write.table(
   params.estimate,
   "sample.Param_estimate.tsv",
   sep="\t",
   quote=F,
   col.names=T,
   row.names=F
)


CNA=rep(".",dim(ascat.output$segments)[1])
CNA[rowSums(ascat.output$segments[,5:6]) > round(ascat.output$ploidy)]="DUP"
CNA[rowSums(ascat.output$segments[,5:6]) < round(ascat.output$ploidy)]="DEL"
output.table=data.frame(ascat.output$segments,CNA=CNA)

write.table(output.table[output.table$CNA != ".",],"sample_CNVcalls.tsv",quote=F,sep="\t",col.names=T,row.names=F)
q(save="yes")

exit

