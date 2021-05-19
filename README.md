# Introduction to DNA-Seq processing for cancer data - CNVs
***By Mathieu Bourgey, Ph.D***  
*https://bitbucket.org/mugqic/mugqic_pipelines*

================================

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

====================================

In this workshop, we will present the main steps that are commonly used to process and to analyze cancer sequencing data. We will focus on whole genome data and SNParray data and we will provide command lines that allow detecting Copy Number Variations (CNV). 


## Introduction

Goals of this session are:

 1. to perform a copy number variation analysis (CNV) on a normal/tumour pair of alignment files (BAMs) produced by the mapping of Illumina short read sequencing data.
 2. to perform a copy number variation analysis (CNV) on a normal/tumour pair of SNParray data files (BAF and LRR) produced by the processing of Illumina chip.

**What are the advantages and limitations of using each type of technology ?**
[solution](solutions/__dataDiff.md)

---------------------------------
# NGS data analysis
In this session  we will produce the main steps of the analysis of NGS data in order to detect CNV and also to estimate cellularity, ploidy of the tumor sample. Sequenza is the tool we will use to perform this analysis. Sequenza is able to perform the CNV analysis in both WGS and WES data. 

It consists of a two Python pre-processing steps followed by a third step in R to infer the samples estimates and to plot the results for interpretation.

In a second time we will also been using IGV to visualise and manually inspect the copy number variation, we inferred in the first part, for validation purposes.

## Data Source
We will be working on a CageKid sample pair, patient C0053.
The CageKid project is part of ICGC and is focused on renal cancer in many of it's forms.
The raw data can be found on EGA and calls, RNA and DNA, can be found on the ICGC portal. 
For more details about [CageKid](http://www.cng.fr/cagekid/)

To ensure reasonable analysis times, we will perform the analysis on a heavily subset pair of BAM files. These files contain just 60Mb of chromosome 2 


# Prepare the Environment
We will use a dataset derived from whole genome sequencing of a clear-cell renal carcinoma patient (Kidney cancer)

```{.bash}
## set environement

#launch docker
docker run --privileged -v /tmp:/tmp --network host -it -w $PWD -v $HOME:$HOME \
--user $UID:$GROUPS -v /etc/group:/etc/group  -v /etc/passwd:/etc/passwd \
-e DISPLAY=$DISPLAY -v /etc/fonts/:/etc/fonts/  c3genomics/genpipes:0.8


export REF=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/

cd $HOME/ebicancerworkshop2021/CNV/NGS

```
### Software requirements
These are all already installed, but here are the original links.

  * [R](https://www.r-project.org/)
  * [sequenza R package](https://cran.r-project.org/web/packages/sequenza/index.html)
  * [samtools](http://sourceforge.net/projects/samtools/)
  * [python 2.7](https://www.python.org/)
  * [IGV](http://www.broadinstitute.org/software/igv/download)

We should load the corresponding modules 

```{.bash}
module load mugqic/samtools/1.4.1 mugqic/R_Bioconductor/3.6.0_3.9 mugqic/python/2.7.14

```

## Original Setup

The initial structure of your folders should look like this:
```
<ROOT>
|-- C0053/                         # fastqs from the center (down sampled)
    `-- normal                     # The blood sample directory
        `-- normal_chr2_60Mb.bam   # Alignment file as  generated through the SNV session
        `-- normal_chr2_60Mb.bai   # Index of the alignment file
    `-- tumor                      # The tumor sample directory
        `-- tumor_chr2_60Mb.bam    # Alignment file as  generated through the SNV session
        `-- tumor_chr2_60Mb.bai    # Index of the alignment file
|-- saved_results                  # Pre-computed files
`-- scripts                        # cheat sheet folder
```
## Generate a seqz file

More info on Sequenza see [User Guide](https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html)

A seqz file contains genotype information, alleles and mutation frequency, and other features. This file is used as input for the R-based part of Sequenza.

The seqz file could be generated from a pileUp file (as shown in SNV session) or directly from the bam file. 

**Should we start from the mpileUp of the SNV session ?** [solution](solutions/__pileUp.md)

### transforming bam files in seqz file
As we haven’t already generated the pileup files, and we are not interested in storing the pileup for further use, we can use the function
`bam2seqz` which converting on the fly to pileup using samtools without storing the pileup file.

**What the impact of converting data on the fly ?** [solution](solutions/__seqz1.md)

```{.bash}
## sequenza preprocessing step 1 - bam 2 seqz format
###Already preprocessed
mkdir -p sequenza

# sequenza-utils bam2seqz \
# -n C0053/normal/normal_chr2_60Mb.bam \
# -t C0053/tumor/tumor_chr2_60Mb.bam \
# --fasta ${REF}/genome/Homo_sapiens.GRCh37.fa  \
# -gc ${REF}/annotations/Homo_sapiens.GRCh37.gc50Base.txt.gz \
# -q 20 \
# -N 20 \
# -C 2:106000000-166000000 | gzip > \
# sequenza/C0053.seqz.gz 

```

To reduce the size of the seqz file, we'll use of a binning function provided in sequenza-utils.py. This binning decreases the memory requirement to load the data into R, and it also speeds up the processing of the sample.

```{.bash}
## sequenza preprocessing step 2 - seqz binning 500bp
sequenza-utils seqz_binning \
 -w 500 \
 -s sequenza/C0053.seqz.gz | gzip > \
 sequenza/C0053.seqz.bin500.gz

```

We can look at the first few lines of the output in the file `sequenza/C0053.seqz.gz` with:

```{.bash}
zless -S sequenza/C0053.seqz.gz

```
This output has one line for each position in the BAMs and includes information on the position, depths, allele frequencies, zygosity, GC in the location.

Note that since many projects might already have been processed with VarScan2, it can be convenient to be able to import such results. For this purpose a simple function is provided within the R package, to convert the output of the somatic and copynumber programs of the VarScan2 suite into the seqz format.

## Exploring the seqz file and depth ratio normalization details

After the aligned sequence data have been pre-processed, the sequenza R package handles all the normalization and analysis steps. Thus, the remainder of this vignette will take place in R.

Let's launch R

```{.bash}
R

```

You should now see the R prompt identified with `>`.


Load the sequenza package

```{.R}
library("sequenza")

```


### Read the seqz file
The seqz file can be read all at once, but processing one chromosome at a time is less demanding on computational resources, especially while processing whole genome data, and might be preferable in case of limited computational resources.

**Why could we process each chromosome individulally ?** [solution](solutions/__seqz2.md)

The function `sequenza.extract` is designed to efficiently access the seqz file and take care of normalization steps. The arguments enable customization of a set of actions listed below:  

 - binning depth ratio and B allele frequency in a desired window size (allowing a desired number of overlapping windows);
 - performing a fast, allele specific segmentation using the copynumber package;
 - filter mutations by frequency and noise.

```{.R}
data.file = "sequenza/C0053.seqz.bin500.gz"
seqzdata = sequenza.extract(data.file)

```

After the raw data is processed, the size of the data is considerably reduced. Typically, the R object resulting from `sequenza.extract` can be stored as a file of a few megabytes, even for whole genome sequencing data.

### Inference of cellularity and ploidy
After the raw data is processed, imported into R, and normalized, we can apply the parameter inference implemented in the package. The function `sequenza.fit` performs the inference using the calculated B allele frequency and depth ratio of the obtained segments. 

```{.R}
CP.example = sequenza.fit(seqzdata)

```

###  Results of model fitting
The last part of the workflow is to apply the estimated parameters. There is an all-in-one function that plots and saves the results, giving control on file names and output directory


```{.R}
sequenza.results(sequenza.extract = seqzdata, 
 cp.table = CP.example, 
 sample.id = "C0053", 
 out.dir="sequenza/results")

```

We can now quit R and explore the generated results

```{.R}
q("yes")

```

## Quit the docker environment


```{.bash}
exit

```
## Sequenza Analysis Results and Visualisation
One of the first and most important estimates that Sequenza provides is the tumour cellularity (the estimated percentage of tumour cells in the tumour genome). This estimate is based on the B allele frequency and depth ratio through the genome and is an important metric to know for interpretation of Sequenza results and for other analyses.


Lets look at the cellularity estimate for our analysis by opening model fit.pdf with the command:

```{.bash}
evince CNV/NGS/sequenza/results/C0053_model_fit.pdf &

```

**What is the graph telling us ?** [solution](solutions/__results1.md)




Let’s now look at the CNV inferences through our genomic block. Open the  genome copy number visualisation file with:

```{.bash}
evince CNV/NGS/sequenza/results/C0053_genome_view.pdf &

```
This file contains three “pages” of copy number events through the entire genomic block. The first page shows copy numbers of the A (red) and B (blue) alleles, the second page shows overall copy number changes and the third page shows the B allele frequency and depth ratio through genomic block.

**What are these graphs telling us ?** [solution](solutions/__results2.md)


We can see how this is a very easy to read output and it lets us immediately see the frequency and severity of copy number events through the genome.



## CNV Visualisation/Confirmation in IGV
Let’s see if we can visualise the CNV events. We will now open IGV and see if we can observe the predicted increase in copy number alterations within our genomic region.

```{.bash}
igv &

```

IGV will take 30 seconds or so to open so just be patient.

For a events of this size (several Mb), we should not be able to easily observe it just by looking at the raw read alignments. In order to see coverage at large scale I rpre-generate the tdf file of each bam files. This means that we can aggregate the average read depth over relatively large chunks of the genome and compare these values between the normal and tumour genomes.

Once IGV is open just load the normal et tumor bam files and zoom on the region `2:100000000-170000000`

**What IGV profiles are telling us ?** [solution](solutions/__visu1.md)

**How can you explain the peak and drop observed in both nromal and tumor ?** [solution](solutions/__visu2.md)

**Are IGV profiles in concordence with sequenza results ?** [solution](solutions/__visu3.md)

If we look at the tumor profiles we can see :   

  1 - the 3 copies state correspond to a mean coverage of 60x
  2 - the 2 copies state correspond to a mean coverage of 50x 
  3 - the 1 copy state correspond to a mean coverage of 40x. 

**How could you explain these values ?** [solution](solutions/__visu4.md)


## Whole genome data

Let’s compare the small genomic block we ran with the entire genome which will be processed through the SNParray technology.



----------------------------------

# SNParray data analysis

SNParray analysis are very similar to NGS data analysis and it is always good to use 2 different technologies to confirm your findings.

start from one LRR signal file and one BAF signal file for each of the germline and matched tumor samples from an individul.

Many software are avaiable for doing CNV call from SNParray. Here is a non-exhaustive list of  softaware that could be used:

#### Proprietary softwares
  1. [GenomeStudio/CNVpartition](http://support.illumina.com/array/array_software/genomestudio/downloads.ilmn) - Illumina
  2. [Genotyping Console/Birdsuite](http://www.broadinstitute.org/science/programs/medical-and-population-genetics/birdsuite/birdsuite-0) - Affymetrix
#### Affymetrix oriented softwares
  1. [Genome Alteration Detection Algorithm (GADA)](https://r-forge.r-project.org/R/?group_id=915)
  2. [Cokgen](http://mendel.gene.cwru.edu/laframboiselab/software.php)
#### Commercial softwares
  1. [Partek Genomics Suite](http://www.partek.com/)
  2. [Golden Helix SNP](http://www.goldenhelix.com/SNP_Variation/CNV_Analysis_Package/index.html)
#### Freely available general software
  1. [PennCNV](http://www.openbioinformatics.org/penncnv/)
  2. [QuantiSNP](https://sites.google.com/site/quantisnp/)
#### Freely available cancer oriented software
  1. [**Allele-Specific Copy number Analysis of Tumors (ASCAT)**](http://heim.ifi.uio.no/bioinf/Projects/ASCAT/)
  2. [OncoSNP](https://sites.google.com/site/oncosnp/)


**What are the major cancer factors that could biais a CNV analysis  ?** [solution](solutions/__challenges.md)

**What are the steps to proceed this analysis ?** [solution](solutions/__SnpAnalysisSummary.md)

## Prepare the Environment
We will use a dataset derived from whole genome sequencing of a clear-cell renal carcinoma patient (Kidney cancer)

```{.bash}
## set environement

#launch docker
docker run --privileged -v /tmp:/tmp --network host \
  -it -w $PWD -v $HOME:$HOME --user $UID:$GROUPS \
  -v /etc/group:/etc/group  -v /etc/passwd:/etc/passwd \
  -e DISPLAY=$DISPLAY -v /etc/fonts/:/etc/fonts/  c3genomics/genpipes:0.8


cd $HOME/ebicancerworkshop2021/CNV/SNParray

```
### Software requirements
These are all already installed, but here are the original links.

  * [R](https://www.r-project.org/)
  * [ASCAT R package*](https://github.com/Crick-CancerGenomics/ascat)

We should load the corresponding modules 

```{.bash}
module load mugqic/R_Bioconductor/3.6.0_3.9

```


## Original Setup

The initial structure of your folders should look like this:
```
<ROOT>
|-- C0053/                         # fastqs from the center (down sampled)
    `-- normal                     # The blood sample directory
        `-- normal2.BAF.tsv         # Beta Allele frequency file
        `-- normal2.LRR.tsv         # Log R Ratio file 
    `-- tumor                      # The tumor sample directory
        `-- tumor2.BAF.tsv          # Beta Allele frequency file
        `-- tumor2.LRR.tsv          # Log R Ratio file
|-- saved_results                  # Pre-computed files
`-- scripts                        # cheat sheet 
```

## SNP data analysis for CNV detection 
In our case, the data are in LRR and BAF format so we'll skip the first processing steps 

###  Plot probe LRR and BAF QC
This steps aim to load the data and plot the BAF and LRR signal to ensure high quality signals.

First let's launch R:

```{.bash}
R

```
Load the ASCAT in R from the folder scripts


```{.R}
library(ASCAT)

```

Load the data into an ASCAT object

```{.R}
ascat.bc = ascat.loadData(
  "C0053/tumor/tumor2.LRR.tsv", 
  "C0053/tumor/tumor2.BAF.tsv", 
  "C0053/normal/normal2.LRR.tsv",
  "C0053/normal/normal2.BAF.tsv")

```

Plot the raw data 

```{.R}
ascat.plotRawData(ascat.bc)

```

look at the graphs `tumor2.germline.png` and `tumor2.tumor.png`

**what stands out from these graphs ?** [solution](solutions/__rawPlot.md)

### Filtering and segmentation of LRR and BAF signal
The next step filter out SNPs which are found to be homozygous for both tumor and normal. it also allows to perform the segmentation of both LRR and BAF signal. The main points of this segmentation is estimate a models of segmentation that should fit between the 2 signals

```{.R}
ascat.seg = ascat.aspcf(ascat.bc)

```

Plot the result off the segmentation

```{.R}
ascat.plotSegmentedData(ascat.seg)

```

look at the `tumor2.ASPCF.png'. We can see after fitering out the homozygous probes the signal is very clean and confirm the presence of several CNA.

### Estimation of the model paremters
This function will use the computed segmentation model and estimate the following sample paramters: 
  1. aberrant cell fraction (cellularity)
  2. tumor ploidy
  3. absolute allele-specific copy number calls (for each allelic probes of the SNP)


```{.R}
ascat.output = ascat.runAscat(ascat.seg)

```
Next save these estimates into a file


```{.R}
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

```

Look at the 'sample.Param_estimate.tsv' text files.  


The estimated aberrant_cell_fraction is 0.53 which means approximately 50% of the cell in the  tumor sample come from the normal.  

The etimated polidy is 2.32 which means approximately a third of the genome is triploid.


### CNV calling from segments
The last step will determine the copy number by simply counting the total number of allele reported to the sample general ploidy.


```{.R}
CNA=rep(".",dim(ascat.output$segments)[1])
CNA[rowSums(ascat.output$segments[,5:6]) > round(ascat.output$ploidy)]="DUP"
CNA[rowSums(ascat.output$segments[,5:6]) < round(ascat.output$ploidy)]="DEL"
output.table=data.frame(ascat.output$segments,CNA=CNA)
```

And we finally save the CNV results into a file

```{.R}
write.table(output.table[output.table$CNA != ".",],"sample_CNVcalls.tsv",quote=F,sep="\t",col.names=T,row.names=F)
q(save="yes")
```


Now you can quite the docker environment

```{.bash}
exit
```


Explore the result files

**What are these results telling us ?** [solution](solutions/__ascat1.md)

The ASCAT analysis have been done on the same sample than the sequenza analysis

**Are the two analyses concordent ?** [solution](solutions/__ascat2.md)

-------------------------------------

## Aknowledgments
I would like to thank and acknowledge Louis Letourneau and Dr. Velimir Gayevskiy for this help and for sharing his material. The format of the tutorial has been inspired from Mar Gonzalez Porta. I also want to acknowledge Joel Fillon, Louis Letrouneau (again), Robert Eveleigh, Edouard Henrion, Francois Lefebvre, Maxime Caron and Guillaume Bourque for the help in building these pipelines and working with all the various datasets.
