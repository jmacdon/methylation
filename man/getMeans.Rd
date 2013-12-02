\name{getMeans}
\alias{getMeans}
\title{Get mean methylation data}
\usage{
  getMeans(bmptab, eset)
}
\arguments{
  \item{bmptab}{The table item from a bumps object created
  from running \code{bumphunter} on methylation data}

  \item{eset}{A GenomicRatioSet or MethylSet.}
}
\value{
  A data.frame containing the mean expression for all
  probes within all differentially methylated regions.
}
\description{
  Compute means from methylation data
}
\details{
  This function is an internal function and not intended
  for use by end users. The purpose is to compute mean
  methylation data that will be used for creating dotplots
}
\author{
  James W. MacDonald (\email{jmacdon@u.washington.edu})
}

