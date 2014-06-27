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
[solution](#technology-differences)

# NGS data analysis

We will start from two alignment files of the tumor and germline samples form the same individual and we will apply a manual read depth procedure to estimate the tumoral copy number.

![read Depth](./images/readDepth.jpg "readDepth")

**What are the steps to proceed this analysis ?**
[solution](#ngs-data-analysis-for-cnv-detection)



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

# Analysis

## Bin count

We use the software [BVAtools](https://bitbucket.org/mugqic/bvatools/downloads) to run this sep of the analysis.

### Environment setup

Open a terminal and set-up the environment in order to be able running BVAtools:

```
export BVATOOLS_JAR=/home/training/Applications/bvatools-1.3/bvatools-1.3-full.jar 
```

### BVAtools overview

Take a look at the BVAtools usage:

```
java7 -jar $BVATOOLS_JAR

java7 -jar $BVATOOLS_JAR bincounter
```

Which parameters should we used ? 
[solution](#bincounter-parameters)


<div class="pagebreak"></div>

# Solutions

## Technology differences

### SNP array
#### Pros
SNP array data offers many advantages in terms of cost, throughput,  turn around time and processing. But mainly the major advantage of this type of approach is the key advantage of SNP arrays is the use of SNP BAF in addtion to LRR which really increase the CNV detection sensitivity and a better identification of the event type.


#### Cons
SNP array are limited to detecting copy-number differences of sequences present in the reference assembly used to design the probes, provide no information on the location of duplicated copies. A major limitation of SNP array is  the size range of alterations detected. In order to provide good CNV call a minimum nbumber of probes is required which allows only the detection of large event and disable the resolution of breakpoints at the single-base-pair level.


### NGS
#### Pros
The most important benefits of NGS technologies are a genome-wide analysis without a priori information, the specificity and linear dynamic range response of NGS data offer many advantages for estimation of copy number. Additionally NGS allows fine detection of small event and several methods of analysis coming from the CGH arrays domain are available.

#### Cons
NGS data are limited  by the cost, the turn around time and the processing complexity. Moreover the major limitations of NGS approach for CNV are the use of uniform reads distribution assumption (False) and the inherente noise in the data introduced by the quality of the reference sequence during the alignment step.

### NGS (CGH) vs. SNP-array data representation
![NGS (CGH) vs. SNP-array data representation](images/CGH_SNParray.jpg "CGH_SNParray")


### Conclusion
When studying CNV, using a combined approach of SNParray and NGS is actually the best alternative.


<div class="pagebreak"></div>


## NGS data analysis for CNV detection 

  Starting from a fastq files (1-9) or from BAM files (4-9):

#### 1. Sequence trimming

  Removing low quality bases 

#### 2. Genome alignment

  Locate individual reads to the reference genome

#### 3. Alignment refinement

  Realign around INDELs, remove duplicates, etc...

#### 4. Bin creation and count 


  The choice of the bin size has a major importance (moving or overlaping windows)

  [step 4 - analysis](#bin-count)

#### 5. Bin count correction (optional)
  1. GC content
  2. Mappability

**It won't be done today: it needs a whole genome data set**

#### 6. Computing LRR

  Normalize count between tumor and germaline then estimate the logRatio in each bin

#### 7. LRR signal smoothing (optional)

  Reduce single point outliers impact before analysis

#### 8. LRR signal segmentation

  Copy number aberrations (CNA) occur in contiguous regions of the chromosome that often cover multiple bins up to whole chromosome arms or chromosomes. The segmentation split the chromosomes into regions of equal copy number that accounts for the noise in the data

#### 9. CNV calling from segments

   Use any method that can differenciate the copy number of each segments 

   Non-exhaustive available methods: 
  * Thresholds
  * Sd deviation
  * Poisson distribution
  * Z-score
  * Event-Wise Testing
  * Etc...

[step 6 to 9  - analysis](../analysis/DNAcopy.md)



<div class="pagebreak"></div>


## Bincounter parameters
| Parameter | explanation | Sugested value  |
| ------------- |:-------------:| -----:|
| --minMapQ | Filter out reads a mapping quality lower or equal at  | 35 |
| --gc      | Compute GC% per window |  _not neccessary_ |
| --ref | Path to the indexed Reference Genome |  _not neccessary_ |
| --refbam | Path to the germline sample bam | ../alignment/normal/normal.sorted.dup.bam |
| --bam | Path to the tumor sample bam | ../alignment/tumor/tumor.sorted.dup.bam |
| --norm | unit for count normalization | chr |
| --windows | bin size | 1000 | 


### bincoutner analysis
Feel free to test different sets of parameters.

Sugested command line

```
cd /home/training/EBI_CNV_workshop2014/NGS_CNV
java7 -jar $BVATOOLS_JAR bincounter --minMapQ 35 --refbam ../alignment/normal/normal.sorted.dup.bam --bam ../alignment/tumor/tumor.sorted.dup.bam --norm chr --windows 1000 > sample_binCount_1kb.tsv
```

The output file should looks like:
```
$ head sample_binCount_1kb.tsv 
chr	start	end	sample_raw	ref_raw	sample_normalized	ref_normalized	ln(sample/ref)
chr1	0	199	0	0	0.0	0.0	NaN
chr1	200	399	0	0	0.0	0.0	NaN
chr1	400	599	0	0	0.0	0.0	NaN
chr1	600	799	0	0	0.0	0.0	NaN
chr1	800	999	0	0	0.0	0.0	NaN
chr1	1000	1199	0	0	0.0	0.0	NaN
chr1	1200	1399	0	0	0.0	0.0	NaN
chr1	1400	1599	0	0	0.0	0.0	NaN
chr1	1600	1799	0	0	0.0	0.0	NaN
```




## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI and Louis Letourneau from MUGQIC. This tutorial use materails form the [BioDicovery.com website](http://www.biodiscovery.com/video-library/), the [Alkan, Coe & Eichler's review article](http://www.nature.com/nrg/journal/v12/n5/full/nrg2958.html) and the [Zhao _et al._ artcile](http://www.biomedcentral.com/1471-2105/14/S11/S1). I would like to thank and acknowledge all of them.


