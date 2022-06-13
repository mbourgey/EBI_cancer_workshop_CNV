## set environement

#launch docker
docker run --privileged -v /tmp:/tmp --network host -it \
    -w $PWD -v $HOME:$HOME -v /etc/fonts/:/etc/fonts/ \
    -v $HOME/cvmfs_caches/:/cvmfs-cache/ c3genomics/genpipes:v2.1.0

module purge 

export REF=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/
export COURSE=/home/training/ebicancerworkshop2022

mkdir -p $COURSE/CNV/NGS

cd $COURSE/CNV/NGS

 
module load mugqic/samtools/1.4.1 mugqic/R_Bioconductor/3.6.0_3.9 mugqic/python/2.7.14

  
## sequenza preprocessing step 1 - bam 2 seqz format
###Already preprocessed
mkdir -p sequenza

sequenza-utils bam2seqz \
-n C0053/normal/normal_chr2_60Mb.bam \
-t C0053/tumor/tumor_chr2_60Mb.bam \
--fasta ${REF}/genome/Homo_sapiens.GRCh37.fa  \
-gc ${REF}/annotations/Homo_sapiens.GRCh37.gc50Base.txt.gz \
-q 20 \
-N 20 \
-C 2:106000000-166000000 | gzip > \
sequenza/C0053.seqz.gz 

 
## sequenza preprocessing step 2 - seqz binning 500bp
sequenza-utils seqz_binning \
 -w 500 \
 -s sequenza/C0053.seqz.gz | gzip > \
 sequenza/C0053.seqz.bin500.gz

 
#zless -S sequenza/C0053.seqz.gz

 
R

 
library("sequenza")

 
data.file = "sequenza/C0053.seqz.bin500.gz"
seqzdata = sequenza.extract(data.file)

 
CP.example = sequenza.fit(seqzdata)

 
sequenza.results(sequenza.extract = seqzdata, 
 cp.table = CP.example, 
 sample.id = "C0053", 
 out.dir="sequenza/results")

 
q("yes")

 
exit
