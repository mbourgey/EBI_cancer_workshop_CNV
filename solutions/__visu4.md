To explain why we are only seeing such a small coverage increase, we need to turn to our good friend mathematics!

 - If a 3 copies state correspond to 60x of coverage. 
 - we expect 2 copies at 40x 
 - 1 copy at 20x 
 - and 0 copy at 0X
 
But here the cellularity vairaible is missing. We estimated a cellularity of ~ 50% which mean that half reads observed in tumor come from the normal cells (30x). So if focus only on reads that really come from tumor cells we would expect:

 - If a 3 copies state correspond to 30x of coverage. 
 - we expect 2 copies at 20x 
 - 1 copy at 10x 
 - and 0 copy at 0X
 
And now if we focus on all the read present in the tumor sample:

 - the 3 copies state corresponds to 30x of read from the tumor cells + 30x from normal cells = 60x. 
 - the 2 copies state  corresponds to 20x of read from the tumor cells + 30x from normal cells = 50x. 
 - the 1 copy state  corresponds to 10x of read from the tumor cells + 30x from normal cells = 40x.

__Mathematics is the key !__