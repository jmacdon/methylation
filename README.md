methylation
===========

This is a simple package for use with the Bioconductor minfi/bumphunter packages, 
to plot regions of methylation using the Gviz package.

Using minfi as an example, one might read in a bunch of Illumina methylation arrays,
normalize, etc, and then use bumphunter to find regions that appear to be differentially
methylated. Once you have a set of regions that might be differentially methylated, 
wouldn't it be sweet to be able to make a plot that shows the genomic region, any genes/transcripts
in that region, and plot the methylation data to show where the differential methylation
is occurring? That's what this package does.
