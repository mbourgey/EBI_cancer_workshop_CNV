#### Ploidy

It is well known that many types of tumors frequently have genomic aberrations involving gain or loss of whole or large parts of chromosomes. Thus, the average ploidy or total genomic content of tumor cells cannot be assumed to be 2N.

Conventional microarray copy number analysis is based on comparing the probe intensities to those of a set of diploid reference samples. This works well for detecting aberrations in diploid non-cancer samples as the normalized intensity of copy number two should coincide for query and reference data. However, many individual tumors have such extensive genomic aberrations that the assumption that the query cells have a genomic content of 2N on average is severely violated. 

#### Cellularity (purity)

Tumor samples are a mix of cancer cells and genetically normal cells due to samlping step.

The proportion of tumor cells can vary considerably, complicating the analysis since the measured signal from any locus will be a combined signal from both tumor and normal cells. If the proportion of tumor cells is too low, aberrations will remain undetected. 

#### Clonality

Copy number aberrations in tumor cells may arise several times throughout the tumor development, and may give rise to different subclones. The extent of the proliferative advantage and the time between occurrence of the aberration and tumor collection influence the proportion of subclones in the tumor cells.

The average copy numbers of heterogenic regions may be non-integer. Long non-integer regions may severely disturb a model-based copy number analysis as they do not fit the pre-determined relationship between signal and copy number.



**We choose ASCAT because it handles samples with aneuploidy and the presence of normal cells and facilitates detection of tumor cell heterogeneity**


