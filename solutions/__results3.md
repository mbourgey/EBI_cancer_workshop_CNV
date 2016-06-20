Let's first look at the model fit

```{.bash}
evince saved_results/preComputed/C0053_full_model_fit.pdf

```

We can see that the parameter estimates based on the whole genome (cellularity 0.47, ploidy 2.5) correspond to thoise estiamtate only on a small part of geneome. which show the efficiency of the method.  


Now let's look at the genome view  

```{.bash}
evince saved_results/preComputed/C0053_full_genome_view.pdf

```

We can see from the graph that major part of the genome is diploid (chr 1, 4, 5-partial, 6, 8, 9-partial, 10, 14, 15, 17, 18, 21 and 22). Some few deletion can be observed (chr 2 and 3) and many large dupications (chr 2, 3, 5, 7, 11, 13, 16, 19, 20). This large amount of duplicated region explain why the ploidy estimate is 2.5.  

We can focus on the chromosome 2

```{.bash}
evince saved_results/preComputed/C0053_full_chromosome_view.pdf

```
Go to the second page of the pdf and look at the bottom panel. We can see that the overall chromosme is triploid. At 130Mb we can see there is deletion of 2 of the 3 copies. This deletion stand until the position 137Mb where the deletion continue but for only 1 of the 3 copies. This deletion stands up to the position 170Mb where the copy number came back to the the basal 3 copies of the chromosome. This fully correspond to what we observed while working only on 60Mb of the chromosome.   

  