\name{makeMethPlot}
\alias{makeMethPlot}
\title{Create methylation region plot}
\usage{
  makeMethPlot(bumpsObj, eset, row, txdb, samps,
    dontuse = "", sexfirst = "Male")
}
\arguments{
  \item{bumpsObj}{The output from bumphunter}

  \item{eset}{Usually a \code{GenomicRatioSet} created by
  the minfi package}

  \item{row}{Which row of the bumpsObj are we plotting?}

  \item{txdb}{A transcript database object (e.g.,
  TxDb.Hsapiens.UCSC.hg19.knownGene)}

  \item{samps}{A data.frame containing sample information.
  This data.frame must contain two columns; one called
  Gender, listing the gender of each sex, and one called
  Category, listing the category for each sample}

  \item{dontuse}{A character vector of Category levels that
  won't be plotted. If using all levels, this argument is
  "".}

  \item{sexfirst}{Which sex should be plotted first? In
  general this is the sex for which the bump being shown
  was significant.}

  \item{fitobj}{An MArrayLM object, after fitting the same
  model as bumphunter. This is the source of the
  methylation data in the plot}
}
\value{
  This function doesn't return anything. It is only called
  for the side effect of creating a plot.
}
\description{
  Create a plot for a region showing differential
  methylation.
}
\details{
  This function is used to plot a genomic region where
  there is differential methylation that includes the
  methylation data, a smoothed line, and any genes in the
  region. This function is in general not called by the end
  user. Instead use \code{methByRegion}.
}
\author{
  James W. MacDonald (\email{jmacdon@u.washington.edu})
}

