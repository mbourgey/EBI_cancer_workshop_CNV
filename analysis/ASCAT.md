# CNV analysis

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
[ASCAT functions](./ASCAT_functions.md)


### Analysis

If necessary take a look at these basic R commands:

[Surviving in R](http://personality-project.org/r/r.commands.html)


To call CNV, we will apply the segmentation and the allele specific copy number estimation step and then ectract the results and plot the data

Try to do the analysis on your own:

1. Go to the working directory  and Open R [solution](solutions/7.1ascat_Ropen.md)
2. Load the ASCAT in R from the folder ../src [solution](solutions/7.2ascat_loadM.md)
3. Load the data into an ASCAT object [solution](solutions/7.3ascat_loadD.md)
4. Plot the raw data [solution](solutions/7.4ascat_plotRaw.md)
5. Apply the segmentation [solution](solutions/7.5ascat_segments.md)
6. Plot the segments [solution](solutions/7.6ascat_plotSeg.md)
7. Compute allele specific copy number [solution](solutions/7.7ascat_ASCP.md)
8. Save aberant cell fraction and ploidy in a file [solution](solutions/7.8ascat_savePlo.md)
9. Call CNA segments [solution](solutions/7.9ascat_CNA.md)
10. Save the CNA [solution](solutions/7.10ascat_saveCNA.md)

###### What can we see from result files ?
[solution](solutions/7.11ascat_exam.md)



