\name{bwplotfun}
\alias{bwplotfun}
\title{Sex stratified mean methylation}
\usage{
  bwplotfun(samps, bumpavg, dontuse = "AGA")
}
\arguments{
  \item{samps}{A data.frame that maps samples to phenotype.
  There must be both a Category and Gender column.}

  \item{bumpavg}{Mean methylation data, usually from a call
  to \code{getMeans}}

  \item{dontuse}{Category levels that are not to be used in
  the dotplot. If there are only two levels, use "".}
}
\value{
  Nothing is returned. Only called for the side effect of
  creating a dotplot.
}
\description{
  Generate sex-stratified dotplot of methylation for a
  given region.
}
\details{
  This is an internal function and not intended for direct
  use. This is intended to create a dotplot stratified by
  sex, to show the mean methylation for a given region of
  the genome, presumably because there is differential
  methylation for at least on sex.
}
\author{
  James W. MacDonald (\email{jmacdon@u.washington.edu})
}

