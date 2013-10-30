\name{methByRegion}
\alias{methByRegion}
\title{Create plots showing differential methylation and correlation to expression data.}
\usage{
  methByRegion(bmpsObj, eset, samps, contname, longname,
    txdb, gene.data = NULL, chip.db = NULL, fitcol,
    cut = 0.001, dontuse = "", orgpkg)
}
\arguments{
  \item{bmpsObj}{The output from \code{bumphunter}}

  \item{eset}{Usually a GenomicRatioSet, created by a call
  to \code{preprocessQuantile} from the minfi package}

  \item{samps}{A data.frame that maps samples to phenotype.
  This data.frame MUST contain columns named Category and
  Gender!}

  \item{contname}{A contrast name, used to name the
  directory where these data will be stored. Usually of the
  form 'this_vs_that'.}

  \item{longname}{A long form of the contrast name, usually
  of the form 'This versus that'}

  \item{txdb}{A transcript database (e.g.,
  TxDb.Hsapiens.UCSC.hg19.knownGene)}

  \item{gene.data}{Default is NULL. Otherwise, a data.frame
  or matrix containing gene expression data. The columns of
  this data.frame must correspond exactly to columns of the
  methylation data. Row names must be array IDs or Entrez
  Gene IDs.}

  \item{chip.db}{The chip-level array data package
  corresponding to the gene expression data. If NULL, the
  assumption will be made that the row.names are ENTREZ
  GENE IDs.}

  \item{fitcol}{Which column of the MArrayLM object
  corresponds to the coefficient tested by
  \code{bumphunter}?}

  \item{cut}{The p-value cutoff used to select significant
  'bumps'.}

  \item{dontuse}{Which Categories from the samps data.frame
  should we NOT use? If only two Category levels, use "".}

  \item{orgpkg}{The organism-level annotation package
  (e.g., org.Hs.eg.db)}

  \item{fit}{An MArrayLM object, created by fitting the
  same model as used by \code{bumphunter}, but probe-wise
  using the limma package}
}
\value{
  This returns an HTMLReportRef that can be used to create
  an index.html page.
}
\description{
  A function to create an HTML page with links to plots
  showing methylation status over each differentially
  methylated region, a dotplot showing sex-stratified mean
  methylation, and an HTML table showing correlation
  between methylation status and gene expression for all
  genes within 1 Mb of the methylation site. This is the
  main function for this package, and is likely the only
  one that an end user should use.
}
\details{
  This function is intended to create plots corresponding
  to a region of the genome that is considered to be
  differentially methylated, usually after running
  \code{bumphunter} to detect differentially methylated
  regions. The methylation region plot will consist of the
  probe-wise methylation, with a smoothed line to indicate
  the portion that is differentially methylated. The
  dotplot will show sex-stratified differential methylation
  (as it is usually safer to do differential methylation
  for each sex separately). If there are expression data
  available, there will be links to tables showing the
  correlation between methylation and gene expression.
}
\author{
  James W. MacDonald (\email{jmacdon@u.washington.edu})
}

