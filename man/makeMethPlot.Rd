\name{makeMethPlot}
\alias{makeMethPlot}
\title{Create methylation region plot}
\usage{
  makeMethPlot(bumpsObj, eset, fitobj, row, fitcol, mart)
}
\arguments{
  \item{bumpsObj}{The table data.frame from the output of
  bumphunter (e.g., bumphunter()$table)}

  \item{eset}{Usually a \code{GenomicRatioSet} created by
  the minfi package}

  \item{fitobj}{An MArrayLM object, after fitting the same
  model as bumphunter. This is the source of the
  methylation data in the plot}

  \item{row}{Which row of the bumpsObj are we plotting?}

  \item{fitcol}{Which column of the fitobj corresponds to
  the coefficient tested by bumphunter?}

  \item{mart}{A connection to the correct biomaRt
  database.}
}
\value{
  This returns a plot of the region in question.
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

