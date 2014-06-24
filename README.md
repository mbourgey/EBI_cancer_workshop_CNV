CNVworkshopEBI
==============

hands-on for NGS/SNParray CNV call trainning

# Introduction
This workshop will show you how to proceed when to detect Copy Number Variants (CNV) using either Next Generation Sequencing (NGS) or SNParray data 

The initial structure of your folders should look like this:
```
/lb/project/mugqic/hgen698_X
|-- reference/               # genome and annotation files
|-- SNParray/                # BAF and logR files per sample
|-- alignment/               # bam files generated previously
    `-- SAMPLE               # One sample directory     

```

# Data
in this hands-on two types of data are available: SNParray and NGS data

**NGS data**

approximately 1M of Illumina paired-end 100 reads covering the region _chr19 : 50500375 - 52501256_ from a whole genome experiment for germline and matched tumor samples.

For each sample:

	1 alignment file (.bam generated in the previous hands-on) on the human reference séquence

**SNParray data**

Illumina data for more than 650000 SNPs located genome wide for germline and matched tumor samples

For each sample:

	1 file containing the LogRatio (LRR) mesurement of each probe: the signal intensities compared to a collection of reference hybridizations

	1 file containing the B Allele Frequency (BAF) mesurement of each probe: the proportion of the total allele signal (A + B) explained by a single allele (A) 



**What are the advantages and limitations of using each type of technology ?**
[solution](solutions/1.dataDiff.md)

# NGS data analysis

We will start from two alignment files of the tumor and germline samples form the same individual and we will apply a manual read depth procedure to estimate the tumoral copy number.

![read Depth](images/readDeapth.jpg "readDepth")

**What are the main steps to proceed this analysis ?**
[solution](solutions/2.NgsAnalysisSummary.md)


## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI and Louis Letourneau from MUGQIC, who I would like to thank and acknowledge.