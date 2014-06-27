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
 
 
-------------------------------- 
 
 

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
[solution](Technology differences)

 
 
 
----------------------------------
 
 
# NGS data analysis

We will start from two alignment files of the tumor and germline samples form the same individual and we will apply a manual read depth procedure to estimate the tumoral copy number.

![read Depth](./images/readDepth.jpg "readDepth")

**What are the steps to proceed this analysis ?**
[solution](NGS data analysis for CNV detection)



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


 
 
 
--------------------------------------------------------
 

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

[solution](Cancer Challenges)


Here is a video showing how the ASCAT software works:

<a href="http://www.youtube.com/watch?feature=player_embedded&v=V6UdEEeVBz8
" target="_blank"><img src="http://img.youtube.com/vi/V6UdEEeVBz8/0.jpg" 
alt="ASCAT CNV" width="240" height="180" border="10" /></a>


[local Shortened version] (videos/ASCAT_Algorithm_BioDiscovery.com.mp4)


**What are the steps to proceed this analysis ?**

[solution](SNP data analysis for CNV detection)




-------------------------------------------------------------



# Analysis



-----------------------------------------------------------------

## 1. Bin count

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
[solution](Bincounter parameters)


 

--------------------------------------------------------------



## 2. DNAcopy

In order to give you a good understanding of the process of CNV calling from NGS data, we will manually do the CNV calls in R. To do this we will apply the 
[DNAcopy](http://www.ncbi.nlm.nih.gov/pubmed/15475419) 
approach orginally developped for CGH array to the NGS data.

### Environment setup

DNAcopy is an 
[R](http://cran.r-project.org/) 
package and it can be downloaded directly from the 
[Bioconductor packe page](http://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html).

### DNAcopy overview

Take a look at the DNAcopy available function:

[DNAcopy manual](http://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf)


**Which functions should we used ?**  

[solution](DNACopy explanation)


### Analysis

If necessary take a look at these basic R commands:

[Surviving in R](http://personality-project.org/r/r.commands.html)


To call CNV, we will apply a smoothing step, followed by a circulary binary segmentation and retain the segement that shows a sufficient ratio deviation (sugested: _[0.5;2]_ less than a 1/2 deviation in both direction) for a sufficient number of consecutive bins (sugested: _10 bins_ )

Try to do the analysis on your own:

1. Open R 
2. Load the DNAcopy package in R 
3. Load the bin count data in R 
4. Clean count data (remove uncovered region) 
5. Normalize count 
6. Compute LogRatio 
7. Create a copy number alteration object from the data 
8. Smooth the LRR signal 
9. Generate segments 
10. Call CNA segments 
11. Save the CNA 

[solution](DNAcopy analysis code)

** What can we see from the result file ?**

[solution](DNAcopy results examination)

To confirm these result, it is recommended to proceed to addtional validation steps:

1. Visual inspection of calls
2. Lab validation

 
------------------------------------------------------------- 
 

## 3. ASCAT

Now that we have some knowledge in CNV calling, we will run a software that could automate the CNV calling for us. To do this we will use the 
[ASCAT](http://www.pnas.org/content/107/39/16910.abstract) 
approach.

### Environment setup

ASCAT is an 
[R](http://cran.r-project.org/) 
set of script that can be downloaded directly from the 
[ASCAT web site](http://heim.ifi.uio.no/bioinf/Projects/ASCAT/).

### ASCAT overview

Unfortunately there is no real documentation for this software.

Here are the available function in the main scripts:
[Analysis](ASCAT functions)


### Analysis

If necessary take a look at these basic R commands:

[Surviving in R](http://personality-project.org/r/r.commands.html)


To call CNV, we will apply the segmentation and the allele specific copy number estimation step and then ectract the results and plot the data

Try to do the analysis on your own:

1. Go to the working directory  and Open R 
2. Load the ASCAT in R from the folder ../src 
3. Load the data into an ASCAT object 
4. Plot the raw data 
5. Apply the segmentation 
6. Plot the segments 
7. Compute allele specific copy number 
8. Save aberant cell fraction and ploidy in a file 
9. Call CNA segments 
10. Save the CNA 

[solution](ASCAT analysis code)



**What can we see from result files ?**

[solution](ASCAT results examination)




----------------------------------



## 4. ASCAT functions


**ASCAT version 2.2, 11/12/2012**

### ascat.loadData

#### Usage

ascat.loadData = function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, Germline_BAF_file = NULL, chrs = c(1:22,"X","Y"), gender = NULL, sexchromosomes = c("X","Y"))

#### exlpanation

Function to read in SNP array data

Input:
  - tumor LogR file
  - tumor BAF file
  - germline LogR file
  - germline BAF file

Output:
  - ASCAT object containing:
    1. Tumor_LogR data matrix
    2. Tumor_BAF data matrix
    3. Tumor_LogR_segmented: placeholder, NULL
    4. Tumor_BAF_segmented: placeholder, NULL
    5. Germline_LogR data matrix
    6. Germline_BAF data matrix
    7. SNPpos: position of all SNPs
    8. ch: a list containing vectors with the indices for each chromosome (e.g. Tumor_LogR[ch[[13]],] will output the Tumor_LogR data of chromosome 13
    9. chr: a list containing vectors with the indices for each distinct part that can be segmented separately (e.g. chromosome arm, stretch of DNA between gaps in the array design)
    10. gender: a vector of gender for each cases ("XX" or "XY"). Default = NULL: all female ("XX")

### ascat.plotRawData

#### Usage

ascat.plotRawData = function(ASCATobj, tumorfiles = "Tumor", germlinefiles = "Germline")

#### exlpanation

Plots SNP array data

Input:
  -ASCAT object 
  -tumorfiles: start of filename for tumor data plots (no plotting if NULL)
  -germlinefiles: start of filename for germline data plots (no plotting if NULL)

### ascat.GCcorrect

#### Usage

ascat.GCcorrect = function(ASCATobj, GCcontentfile = NULL)

#### exlpanation

Correct the signal for GC biais

Input:
  - ASCAT object
  - GC model file

Output:
  - ASCAT object

Note that probes not present in the GCcontentfile will be lost from the results


### ascat.aspcf

#### Usage

ascat.aspcf = function(ASCATobj, selectsamples = 1:length(ASCATobj$samples), ascat.gg = NULL, penalty = 25)

#### exlpanation

run the ASPCF segmentation

Input:
  - ASCAT object  
  - selectsamples: a vector containing the sample(number)s for which segmentation will be apply
  - germline genotypes (NULL if germline data is available)
  - penalty: penalty of introducing an additional ASPCF breakpoint (expert parameter, don't adapt unless you know what you're doing)

Output:
  - ASCAT object


### ascat.plotSegmentedData

#### Usage

ascat.plotSegmentedData = function(ASCATobj, filenames = "ASPCF")

#### exlpanation

Plots SNP array data

Input:
  - ASCAT object (e.g. from ASCAT.ASPCF) and plots the SNP array data before and after segmentation
  - filenames: start of the names of the output files
 

### ascat.runAscat

#### Usage

ascat.runAscat = function(ASCATobj, gamma = 0.55, sunrisefiles = "sunrise", profilefiles = "ASCATprofile", rawprofilefiles = "rawprofile", aberrationreliabilityfiles = "aberrationreliability")

#### exlpanation

The ASCAT main function, calculating the allele-specific copy numbers

Input:
  - ASCAT object
  - gamma: technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 % aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
  - sunrisefiles: filenames of 'sunrise plots' (plots with distance matrices)
  - profilefiles: filenames of ASCAT profiles (plots with final allele-specific copy number)
  - rawprofilefiles: filenames of raw copy number profiles

Output:
  - ASCAT output object containing:
    1. nA: copy number of the A allele
    2. nB: copy number of the B allele
    3. aberrantcellfraction: the aberrant cell fraction of all arrays
    4. ploidy: the ploidy of all arrays
    5. failedarrays: arrays on which ASCAT analysis failed
    6. nonaberrantarrays = arrays on which ASCAT analysis indicates that they so virtually no aberrations
    7. segments: an array containing the copy number segments of each sample (not including failed arrays)
    8. segments_raw: an array containing the copy number segments of each sample without any rounding applied
    
Note: for copy number only probes, nA contains the copy number value and nB = 0





-------------------------------------



## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI and Louis Letourneau from MUGQIC. This tutorial use materails form the [BioDicovery.com website](http://www.biodiscovery.com/video-library/), the [Alkan, Coe & Eichler's review article](http://www.nature.com/nrg/journal/v12/n5/full/nrg2958.html) and the [Zhao _et al._ artcile](http://www.biomedcentral.com/1471-2105/14/S11/S1). I would like to thank and acknowledge all of them.


