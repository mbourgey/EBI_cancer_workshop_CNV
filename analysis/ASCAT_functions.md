# ASCAT version 2.2, 11/12/2012

## ascat.loadData

### Usage

ascat.loadData = function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, Germline_BAF_file = NULL, chrs = c(1:22,"X","Y"), gender = NULL, sexchromosomes = c("X","Y"))

### exlpanation

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

## ascat.plotRawData

### Usage

ascat.plotRawData = function(ASCATobj, tumorfiles = "Tumor", germlinefiles = "Germline")

### exlpanation

Plots SNP array data

Input:
  -ASCAT object 
  -tumorfiles: start of filename for tumor data plots (no plotting if NULL)
  -germlinefiles: start of filename for germline data plots (no plotting if NULL)

## ascat.GCcorrect

### Usage

ascat.GCcorrect = function(ASCATobj, GCcontentfile = NULL)

### exlpanation

Correct the signal for GC biais

Input:
  - ASCAT object
  - GC model file

Output:
  - ASCAT object

Note that probes not present in the GCcontentfile will be lost from the results


## ascat.aspcf

### Usage

ascat.aspcf = function(ASCATobj, selectsamples = 1:length(ASCATobj$samples), ascat.gg = NULL, penalty = 25)

### exlpanation

run the ASPCF segmentation

Input:
  - ASCAT object  
  - selectsamples: a vector containing the sample(number)s for which segmentation will be apply
  - germline genotypes (NULL if germline data is available)
  - penalty: penalty of introducing an additional ASPCF breakpoint (expert parameter, don't adapt unless you know what you're doing)

Output:
  - ASCAT object


## ascat.plotSegmentedData

### Usage

ascat.plotSegmentedData = function(ASCATobj, filenames = "ASPCF")

### exlpanation

Plots SNP array data

Input:
  - ASCAT object (e.g. from ASCAT.ASPCF) and plots the SNP array data before and after segmentation
  - filenames: start of the names of the output files
  
## ascat.runAscat

### Usage

ascat.runAscat = function(ASCATobj, gamma = 0.55, sunrisefiles = "sunrise", profilefiles = "ASCATprofile", rawprofilefiles = "rawprofile", aberrationreliabilityfiles = "aberrationreliability")

### exlpanation

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


