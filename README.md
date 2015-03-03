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

In addition, this package uses the Bioconductor ReportingTools package to create HTML tables that
present the results for a given analysis in an HTML table, with clickable links to the above plots. 
This tends to be a nice format to give to collaborators, so they can easily peruse the results, and 
also look at graphical representations of significant genomic regions.

If you also have gene expression data, you can look for genes in CIS (within say 1 Mb of a region that shows differential methylation) that have changes in expression that correlate with changes in methylation. This package allows you to generate HTML tables for each methylation region, showing the genes in CIS, statistics that measure correlation between methylation and expression, as well as plots showing the relationship.
