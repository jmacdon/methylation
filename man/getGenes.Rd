\name{getGenes}
\alias{getGenes}
\title{Get gene symbols}
\usage{
  getGenes(plus, minus, mart, peakloc)
}
\arguments{
  \item{plus}{Ensembl IDs on the plus strand}

  \item{minus}{Ensembl IDs on the minus strand}

  \item{mart}{A connection to the correct biomaRt
  database.}

  \item{peakloc}{The location of the peak itself}
}
\value{
  This creates a TextOverlay object for inclusion to the
  methylation region plot
}
\description{
  Get gene symbols given an Ensembl ID
}
\details{
  This function is used by the \code{makeMethPlot} function
  in order to annotate genes. This is an internal function
  and is not intended to be called by the end user.
}
\author{
  James W. MacDonald (\email{jmacdon@u.washington.edu})
}

