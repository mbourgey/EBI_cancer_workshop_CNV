CNVworkshopEBI
==============

hands-on for NGS/SNParray CNV call trainning

# Introduction
This workshop will show you how to proceed when to detect Copy Number Variants (CNV) using either Next Generation Sequencing (NGS) or SNParray data 

The initial structure of your folders should look like this:
```
UKworkshop/
|-- src/                     # scripts
|-- NGS_CNV/                 # Analysis folder for NGS data
|-- SNP_CNV/                 # Analysis folder for SNParray data
|-- SNParray/                # BAF and logR files per sample
|-- alignment/               # bam files generated previously
    `-- normal               # Normal sample directory
    `-- tumor                # Tumor sample directory

```

# Data
in this hands-on two types of data are available: SNParray and NGS data

**NGS data**

approximately 1M of Illumina paired-end 100 reads covering the region _chr19 : 50500375 - 52501256_ from a whole genome experiment for germline and matched tumor samples.

For each sample:

	1 alignment file (.bam generated in the previous hands-on) on the human reference sequence

**SNParray data**

Illumina data for more than 650000 SNPs located genome wide for germline and matched tumor samples

For each sample:

	1 file containing the LogRatio (LRR) mesurement of each probe: the signal intensities compared to a collection of reference hybridizations
	1 file containing the B Allele Frequency (BAF) mesurement of each probe: the proportion of the total allele signal (A + B) explained by a single allele (A) 



**What are the advantages and limitations of using each type of technology ?**
[solution](solutions/1.dataDiff.md)

# NGS data analysis

We will start from two alignment files of the tumor and germline samples form the same individual and we will apply a manual read depth procedure to estimate the tumoral copy number.

![read Depth](./images/readDepth.jpg "readDepth")

**What are the steps to proceed this analysis ?**
[solution](solutions/2.NgsAnalysisSummary.md)


For your information, here is a non-exhaustive list of available softwares for calling CNV using whole genome NGS Data:


| Tool | URL | Language | Input | Comments |
|:--:|:-:|:------:|:---:|:------:|
| SegSeq | http://www.broad.mit.edu/cancer/pub/solexa_copy_numbers/ | Matlab | Aligned read positions | Detecting CNV breakpoints using massively parallel sequence data |
| CNV-seq | http://tiger.dbs.nus.edu.sg/cnv-seq | Perl, R | Aligned read positions | Identifying CNVs using the difference of observed copy number ratios |
| RDXplorer | http://rdxplorer.sourceforge.net  | Python, Shell | BAM | Detecting CNVs through event-wise testing algorithm on normalized read depth of coverage |
| BIC-seq | http://compbio.med.harvard.edu/Supplements/PNAS11.html | Perl, R | BAM | Using the Bayesian information criterion to detect CNVs based on uniquely mapped reads |
| CNAsega | http://www.compbio.group.cam.ac.uk/software/cnaseg | R | BAM | Using flowcell-to-flowcell variability in cancer and control samples to reduce false positives |
| cn.MOPS | http://www.bioinf.jku.at/software/cnmops/ | R | BAM/read count matrices | Modelling of read depths across samples at each genomic position using mixture Poisson model |
| JointSLMb | http://nar.oxfordjournals.org/content/suppl/2011/02/16/gkr068.DC1/JointSLM_R_Package.zip | R | SAM/BAM | Population-based approach to detect common CNVs using read depth data |
| ReadDepth | http://code.google.com/p/readdepth | R | BED files | Using breakpoints to increase the resolution of CNV detection from low-coverage reads |
| rSW-seqa | http://compbio.med.harvard.edu/Supplements/BMCBioinfo10-2.html | C | Aligned read positions | Identifying CNVs by comparing matched tumor and control sample |
| CNVnator | http://sv.gersteinlab.org | C++ | BAM | Using mean-shift approach and performing multiple-bandwidth partitioning and GC correction |
| CNVnorma | http://www.precancer.leeds.ac.uk/cnanorm | R | Aligned read positions | Identifying contamination level with normal cells |
| CMDS | https://dsgweb.wustl.edu/qunyuan/software/cmds | C, R | Aligned read positions | Discovering CNVs from multiple samples |
| mrCaNaVar | http://mrcanavar.sourceforge.net | C | SAM | A tool to detect large segmental duplications and insertions |
| cnvHMM | http://genome.wustl.edu/software/cnvhmm | C | Consensus sequence from SAMtools | Using HMM to detect CNV |
| _**PopSV**_ | NA | R | BAM | Use a population of controls sample to detect individual CNV |
| _**DNACRD**_ | NA | R | bin read count | Correct for GC and mappability and use a sample specific outlier detection approach | 


# SNParray data analysis

SNParray analysis are very similar to NGS data analysis while incorporating the additional information bring by the SNP: the BAF. In this analysis we will start from one LRR signal file and one BAF signal file for each of the germline and matched tumor samples from an individul.

Many software are avaiable for doing CNV call from SNParray. Here is a non-exhaustive list of  softaware that could be used:

1. Proprietary softwares
  1. [GenomeStudio/CNVpartition](http://support.illumina.com/array/array_software/genomestudio/downloads.ilmn) - Illumina
  2. [Genotyping Console/Birdsuite](http://www.broadinstitute.org/science/programs/medical-and-population-genetics/birdsuite/birdsuite-0) - Affymetrix
2. Affymetrix oriented softwares
  1. [Genome Alteration Detection Algorithm (GADA)](https://r-forge.r-project.org/R/?group_id=915)
  2. [Cokgen](http://mendel.gene.cwru.edu/laframboiselab/software.php)
3. Commercial softwares
  1. [Partek Genomics Suite](http://www.partek.com/)
  2. [Golden Helix SNP](http://www.goldenhelix.com/SNP_Variation/CNV_Analysis_Package/index.html)
4. Freely available general software
  1. [PennCNV](http://www.openbioinformatics.org/penncnv/)
  2. [QuantiSNP](https://sites.google.com/site/quantisnp/)
5. Freely available cancer oriented software
  1. [**Allele-Specific Copy number Analysis of Tumors (ASCAT)**](http://heim.ifi.uio.no/bioinf/Projects/ASCAT/)
  2. [OncoSNP](https://sites.google.com/site/oncosnp/)


**What are the major cancer factors that could biais a CNV analysis  ?**

[solution](solutions/5.cancerChallenge.md)


Here is a video showing how the ASCAT software works:

<a href="http://www.youtube.com/watch?feature=player_embedded&v=V6UdEEeVBz8
" target="_blank"><img src="http://img.youtube.com/vi/V6UdEEeVBz8/0.jpg" 
alt="ASCAT CNV" width="240" height="180" border="10" /></a>


[local Shortened version] (videos/ASCAT_Algorithm_BioDiscovery.com.mp4)


**What are the steps to proceed this analysis ?**

[solution](solutions/6.SnpAnalysisSummary.md)



## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI and Louis Letourneau from MUGQIC. This tutorial use materails form the [BioDicovery.com website](http://www.biodiscovery.com/video-library/), the [Alkan, Coe & Eichler's review article](http://www.nature.com/nrg/journal/v12/n5/full/nrg2958.html) and the [Zhao _et al._ artcile](http://www.biomedcentral.com/1471-2105/14/S11/S1). I would like to thank and acknowledge all of them.


