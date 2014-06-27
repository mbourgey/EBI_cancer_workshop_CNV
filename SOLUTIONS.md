2014 EBI Cancer Genomics workshop - CNV analysis - Solutions
==============
Mathieu Bourgey, Ph.D

MCGill University and Genome Quebec Innovation Center

mbourgey@genomequebec.com 

--------------------------------------------------------------

## Technology differences

### SNP array
#### Pros
SNP array data offers many advantages in terms of cost, throughput,  turn around time and processing. But the key advantage of SNP arrays is the use of SNP BAF in addtion to LRR which really increase the CNV detection sensitivity and a better identification of the event type.


#### Cons
SNP array are limited to detecting copy-number differences of sequences present in the reference assembly used to design the probes, provide no information on the location of duplicated copies. A major limitation of SNP array is the size range of events detected. In order to provide good CNV calls a minimum nbumber of probes is required which allows only the detection of large events and disable the resolution of breakpoints at the single-base-pair level.


### NGS
#### Pros
The most important benefits of NGS technologies are a genome-wide analysis without a priori information, the specificity and linear dynamic range response of NGS data offer many advantages for estimation of copy number. Additionally NGS allows fine detection of small event and several methods of analysis coming from the CGH arrays domain are available.

#### Cons
NGS data are limited  by the cost, the turn around time and the processing complexity. Moreover the major limitations of NGS approach for CNV are the use of uniform reads distribution assumption (False) and the inherente noise in the data introduced by the quality of the reference sequence during the alignment step.

### NGS (CGH) vs. SNP-array data representation
![NGS (CGH) vs. SNP-array data representation](images/CGH_SNParray.jpg "CGH_SNParray")


### Conclusion
When studying CNV, using a combined approach of SNParray and NGS is actually the best alternative.



----------------------------------




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

  [step 4 - analysis](Bin count)

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

[step 6 to 9  - analysis](DNAcopy)



-----------------------------------



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


--------------------------------------


## DNACopy explanation

Here is the  minimal list of functions in the DNAcopy package that are used to run a basic CNV calling algorithm

| Function | Explanation | 
|:---------|-----------:|
| CNA | Creates a _Copy Number Array_ data object | 
| DNAcopy | R object resulting of the segmentation step |
| smooth.CNA | Smooth the signal to reduce outliers points |
| segment | Find segments harboring similar signal using the CBS algorithm |
| _plot.DNAcopy_ | optional: plot the result of the segmentation |
| _segments.summary_ | optional: provides statistics from segements |




-------------------------------------------


## DNAcopy analysis code
#### 1. Open R
 
```
R
```

#### 2. Load the DNAcopy package in R 
```
library("DNAcopy")
```

#### 3. Load the bin count data in R 
```
data=read.table("sample_binCount_1kb.tsv",header=T)
head(data)
```
```
  chr start  end sample_raw ref_raw sample_normalized ref_normalized
1   1     0  199          0       0                 0              0
2   1   200  399          0       0                 0              0
3   1   400  599          0       0                 0              0
4   1   600  799          0       0                 0              0
5   1   800  999          0       0                 0              0
6   1  1000 1199          0       0                 0              0
  ln.sample.ref.
1            NaN
2            NaN
3            NaN
4            NaN
5            NaN
6            NaN

```


#### 4. Clean count data (remove uncovered region) 
```
dataClean=data[data[,4] > 0 | data[,5] > 0,1:5]
head(dataClean)
```

```
         chr    start      end sample_raw ref_raw
13549734  19 50500200 50500399         70     120
13549735  19 50500400 50500599        185     328
13549736  19 50500600 50500799        188     349
13549737  19 50500800 50500999        174     363
13549738  19 50501000 50501199        146     312
13549739  19 50501200 50501399        198     341
```

#### 5. Normalize count 
```
dataNorm=cbind(dataClean,dataClean[,4]/sum(dataClean[,4]),dataClean[,5]/sum(dataClean[,5]))
head(dataNorm)
```

```
         chr    start      end sample_raw ref_raw
13549734  19 50500200 50500399         70     120
13549735  19 50500400 50500599        185     328
13549736  19 50500600 50500799        188     349
13549737  19 50500800 50500999        174     363
13549738  19 50501000 50501199        146     312
13549739  19 50501200 50501399        198     341
         dataClean[, 4]/sum(dataClean[, 4]) dataClean[, 5]/sum(dataClean[, 5])
13549734                       4.080481e-05                       3.686617e-05
13549735                       1.078413e-04                       1.007675e-04
13549736                       1.095901e-04                       1.072191e-04
13549737                       1.014291e-04                       1.115202e-04
13549738                       8.510718e-05                       9.585203e-05
13549739                       1.154193e-04                       1.047614e-04
```


#### 6. Compute LogRatio 
```
Chr=dataNorm[,1]
Pos=dataNorm[,2]
logR=log2(dataNorm[,6]/dataNorm[,7])
logR[1:6]
```

```
[1]  0.14275958  0.09694868  0.03126659 -0.13555779 -0.16964886  0.13851792
```




#### 7. Create a copy number alteration object from the data 
```
CNA.object=CNA(logR,Chr,Pos, data.type="logratio",sampleid="TNratio")
CNA.object
```

```
Number of Samples 1 
Number of Probes  9818 
Data Type         logratio
```



#### 8. Smooth the LRR signal 
```
smoothed.CNA.object=smooth.CNA(CNA.object)
smoothed.CNA.object
```

```
Number of Samples 1 
Number of Probes  9818 
Data Type         logratio 
```



#### 9. Generate segments 
```
segment1=segment(smoothed.CNA.object, verbose=1)
head(segment1$output)
```

```
       ID chrom loc.start  loc.end num.mark seg.mean
1 TNratio    19  50500200 50595200      476   0.0082
2 TNratio    19  50596400 50598600        7  -0.7104
3 TNratio    19  50598800 50609200        9   0.4972
4 TNratio    19  50609400 50636000       15  -0.7169
5 TNratio    19  50636600 50637600        6  -0.0457
6 TNratio    19  50637800 50638800        6   0.5814
```



#### 10. Call CNA segments 
```
duplicationTh=log2(2)
deletionTh=log2(0.5)
minBin=10
CNAcall=segment1$output[segment1$output[,5] >= minBin & (segment1$output[,6] <= deletionTh | segment1$output[,6] >= duplicationTh) ,]
head(CNAcall)
```

```
        ID chrom loc.start  loc.end num.mark seg.mean
81 TNratio    19  52133600 52149400       80  -5.6947

```



#### 11. Save the CNA 
```
CNAtype=rep(".",dim(CNAcall)[1])
names(CNAtype)="CNA_type
CNAtype[CNAcall[,6] >= duplicationTh]="DUP"
CNAtype[CNAcall[,6] <= deletionTh]="DEL"
CNAcallfinal=cbind(CNAcall,CNAtype)
write.table(CNAcallfinal,"sampleCNAcall.tsv",sep="\t",quote=F,col.names=T,row.names=F)
q(save="yes")
```

you can then look at the call file using this command:

```
head sampleCNAcall.tsv
```

```
ID	chrom	loc.start	loc.end	num.mark	seg.mean	CNAtype
TNratio	19	52133600	52149400	80  	-5.6947	DEL
```



--------------------------------------------


## DNAcopy results examination

This file shows the presence of large deleton event in the chr19 region (19  52133600 52149400).


If you reduce the threshold value you will start to see many region that will pop-up all along the candidat region.
In fact one of the major issue of most of the CNV caller from NGS data: due to not so uniform coverage, local variation occurs in the read depth which can produce noise in the resulls:

1. scattered calls
2. false positive calls
3. false negative calls




----------------------------------------------


## Cancer challenges

### Tumor ploidy

It is well known that many types of tumors frequently have genomic aberrations involving gain or loss of whole or large parts of chromosomes. Thus, the average ploidy or total genomic content of tumor cells cannot be assumed to be 2N.

Conventional microarray copy number analysis is based on comparing the probe intensities to those of a set of diploid reference samples. This works well for detecting aberrations in diploid non-cancer samples as the normalized intensity of copy number two should coincide for query and reference data. However, many individual tumors have such extensive genomic aberrations that the assumption that the query cells have a genomic content of 2N on average is severely violated. 

### Tumor heterogeneity

Tumor samples are a mix of cancer cells and genetically normal cells due to samping step.

The proportion of tumor cells can vary considerably, complicating the analysis since the measured signal from any locus will be a combined signal from both tumor and non-tumor cells. If the proportion of tumor cells is too low, aberrations will remain undetected. 

### Tumor cell heterogeneity

Copy number aberrations in tumor cells may arise several times throughout tumor development, and may give rise to different subclones. The extent of the proliferative advantage and the time between occurrence of the aberration and tumor collection influence the proportion of tumor cells with each aberration [15].

The average copy numbers of heterogenic genomic regions may be non-integer. Long non-integer regions may severely disturb a model-based copy number analysis as they do not fit the pre-determined relationship between signal and copy number.


**We choose ASCAT because it handles samples with aneuploidy and the presence of normal cells and facilitates detection of tumor cell heterogeneity**


--------------------------------------------------------


## SNP data analysis for CNV detection 

Starting from a CEL files (1-9) or from LRR and BAF files (3-9):

#### 1. Normalize probe intensity

   Normalise sample probe intensity against some reference HapMap pooled samples in order to generate a LRR values 

#### 2. Generate genotyping from CEL files

   Contrast the normalized signal intensity for each allelic probe of a same SNP in order to generate BAF values

#### 3. Probe filtering

   Remove homozygote SNPs

#### 4. simultaneous segmenetation of LRR and BAF signal

   Models of segmentation should fit between the 2 signal


#### 5. Estimation of the model paremters
  1. aberrant cell fraction
  2. tumor ploidy
  3. absolute allele-specific copy number calls (for each allelic probes of the SNP)
  

#### 6. CNV calling from segments

   determine the copy number by simply counting the total number of allele reported to the sample ploidy 

[step 3 to 6  - analysis](ASCAT)





---------------------------------------------------------

## ASCAT analysis code

### 1. Go to the working directory  and Open R 
```
cd /home/training/EBI_CNV_workshop2014/SNP_CNV/
R
```


### 2. Load the ASCAT in R from the folder ../src 
```
source("../src/ascat-2.2.R")
```



### 3. Load the data into an ASCAT object 
```
ascat.raw = ascat.loadData("../SNParray/tumor2.LRR.tsv","../SNParray/tumor2.BAF.tsv", "../SNParray/normal2.LRR.tsv","../SNParray/normal2.BAF.tsv")
```



### 4. Plot the raw data 
```
ascat.plotRawData(ascat.raw)
```



### 5. Apply the segmentation 
```
ascat.seg = ascat.aspcf(ascat.raw)
```




### 6. Plot the segments 
```
ascat.plotSegmentedData(ascat.seg)
```



### 7. Compute allele specific copy number 
```
ascat.output = ascat.runAscat(ascat.seg)
```



### 8. Save aberant cell fraction and ploidy in a file 
```
params.estimate=data.frame(Sample=names(ascat.output$aberrantcellfraction),Aberrant_cell_fraction=round(ascat.output$aberrantcellfraction,2),Ploidy=round(ascat.output$ploidy,2))
write.table(params.estimate,"sample.Param_estimate.tsv",sep="\t",quote=F,col.names=T,row.names=F)
```



### 9. Call CNA segments 
```
CNA=rep(".",dim(ascat.output$segments)[1])
CNA[rowSums(output.table[,5:6]) > round(ascat.output$ploidy)]="DUP"
CNA[rowSums(output.table[,5:6]) < round(ascat.output$ploidy)]="DEL"
output.table=data.frame(ascat.output$segments,CNA=CNA)
```


### 10. Save the CNA 
```
write.table(output.table[output.table$CNA != ".",],"sample_CNVcalls.tsv",quote=F,sep="\t",col.names=T,row.names=F)
q(save="yes")
```



-------------------------------------------------------


## ASCAT results examination

These files shows the presence of large deletion and duplication evenet all along the genome of this indivduals.

Particularly in the chromosome 12 a very impressive duplication can be observed.

The only limitation of this approach is the size of event that could be detected. Few kb events or less will be missed.


----------------------------------------------------------




## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI and Louis Letourneau from MUGQIC. This tutorial use materails form the [BioDicovery.com website](http://www.biodiscovery.com/video-library/), the [Alkan, Coe & Eichler's review article](http://www.nature.com/nrg/journal/v12/n5/full/nrg2958.html) and the [Zhao _et al._ artcile](http://www.biomedcentral.com/1471-2105/14/S11/S1). I would like to thank and acknowledge all of them.


