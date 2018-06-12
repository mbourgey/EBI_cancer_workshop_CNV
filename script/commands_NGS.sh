## set environement
export SEQUENZA_UTILS=/home/training/R/x86_64-pc-linux-gnu-library/3.4/sequenza/exec/sequenza-utils.py
export REF=/home/training/ebicancerworkshop2018/reference

cd $HOME/ebicancerworkshop2018/CNV/NGS

## sequenza preprocessing step 1 - bam 2 seqz format
mkdir -p sequenza

${SEQUENZA_UTILS} bam2seqz \
 -n normal/normal.sorted.dup.recal.bam \
 -t tumor/tumor.sorted.dup.recal.bam \
 --fasta ${REF}/Homo_sapiens.GRCh37.fa  \
 -gc ${REF}/Homo_sapiens.GRCh37.gc50Base.txt.gz \
 -q 20 \
 -N 20 \
 -C 2:106000000-166000000 | gzip > \
 sequenza/C0053.seqz.gz 


## sequenza preprocessing step 2 - seqz binning 500bp
${SEQUENZA_UTILS} seqz-binning \
 -w 500 \
 -s sequenza/C0053.seqz.gz | gzip > \
 sequenza/C0053.seqz.bin500.gz

zless -S sequenza/C0053.seqz.gz

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

evince sequenza/results/C0053_model_fit.pdf

evince sequenza/results/C0053_genome_view.pdf

igv &

