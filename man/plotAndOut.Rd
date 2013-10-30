\name{plotAndOut}
\alias{plotAndOut}
\title{Create tables and plots showing correlation between methylation and gene expression data.}
\usage{
  plotAndOut(lstitm, eset, prb, samps, file, contname,
    orgpkg)
}
\arguments{
  \item{lstitm}{A data.frame containing gene expression
  data for genes that are within the 1 Mb region centered
  on the methylation 'bump'.}

  \item{eset}{Usually a GenomicRatioSet, created by a call
  to \code{preprocessQuantile} from the mimfi package}

  \item{prb}{A vector of probe IDs that correspond to the
  differently methylated (Illumina 450K) probe IDs.}

  \item{samps}{A data.frame that maps samples to phenotype.
  This data.frame MUST contain columns named Category and
  Gender!}

  \item{file}{A filename. In general this is the genomic
  region, separated by underscores (e.g.,
  chr1_12345_23456)}

  \item{contname}{A contrast name. Usually this is
  lowercase and separated by underscores (e.g.,
  this_vs_that)}

  \item{orgpkg}{An organism level annotation package (e.g.,
  org.Hs.eg.db)}
}
\value{
  This returns an HTMLReportRef.
}
\description{
  A function to create tables showing correlation between
  methylation status and gene expression.
}
\details{
  This function is not intended to be called by the end
  user. Instead it is intended to create the HTML table
  showing correlation between gene gene expression and
  methylation status, as well as xyplots that show the
  correlation graphically.
}
\author{
  James W. MacDonald (\email{jmacdon@u.washington.edu})
}

