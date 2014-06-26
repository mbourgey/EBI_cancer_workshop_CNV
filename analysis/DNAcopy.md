# CNV analysis
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

[solution](../solutions/4.DnaCopy_explanation.md)


### Analysis

If necessary take a look at these basic R commands:

[Surviving in R](http://personality-project.org/r/r.commands.html)


To call CNV, we will apply a smoothing step, followed by a circulary binary segmentation and retain the segement that shows a sufficient ratio deviation (sugested: _[0.5;2]_ less than a 1/2 deviation in both direction) for a sufficient number of consecutive bins (sugested: _10 bins_ )

Try to do the analysis on your own:

1. Open R [solution](../solutions/4.1DnaCopy_Ropen.md)
2. Load the DNAcopy package in R [solution](../solutions/4.2DnaCopy_loadDC.md)
3. Load the bin count data in R [solution](../solutions/4.3DnaCopy_loadData.md)
4. Clean count data (remove uncovered region) [solution](../solutions/4.4DnaCopy_cleanData.md)
5. Normalize count [solution](../solutions/4.5DnaCopy_normData.md)
6. Compute LogRatio [solution](../solutions/4.6DnaCopy_LRR.md)
7. Create a copy number alteration object from the data [solution](../solutions/4.7DnaCopy_CNA.md)
8. Smooth the LRR signal [solution](../solutions/4.8DnaCopy_smooth.md)
9. Generate segments [solution](../solutions/4.9DnaCopy_segments.md)
10. Call CNA segments [solution](../solutions/4.10DnaCopy_calls.md)
11. Save the CNA [solution](../solutions/4.11DnaCopy_output.md)

###### What can we see from the result file ?
[solution](../solutions/4.12DnaCopy_exam.md)

To confirm these result, it is recommended to proceed to addtional validation steps:

1. Visual inspection of calls
2. Lab validation


