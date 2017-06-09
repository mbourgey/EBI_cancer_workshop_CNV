## set environement
export SEQUENZA_UTILS=/usr/local/lib/R/site-library/sequenza/exec/sequenza-utils.py
export REF=/home/training/ebicancerworkshop2017/reference

cd $HOME/ebicancerworkshop2017/CNV/NGS

## sequenza preprocessing step 1 - bam 2 seqz format
mkdir -p sequenza

${SEQUENZA_UTILS} bam2seqz \
 -n C0053/normal/normal_chr2_60Mb.bam \
 -t C0053/tumor/tumor_chr2_60Mb.bam \
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

