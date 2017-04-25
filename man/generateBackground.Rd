\name{generateBackground}
\alias{generateBackground}
\title{
  Generate median background read depth
}
\description{
  Normalizes a set of columns representing read counts from different
  samples by their mean.  Then calculates a statistic across the rows of
  normalized counts.
}
\usage{
  generateBackground(sample.names, gr, fn=median)
}
\arguments{
  \item{sample.names}{
    A vector of metadata column names in gr to be used as background samples
  }
  \item{gr}{
    A GRanges object containing the read counts
  }
  \item{fn}{
    The statistic to be applied across the rows of normalized counts. 
    Defaults to median, but the standard deviation can also be
    calculated in this way.
    }
}
\value{
  The value of fn applied across the rows of normalized read counts. 
}
\examples{
  data(exomecounts)
  sample.names <- grep("HG.+",colnames(mcols(exomecounts)),value=TRUE)
  exomecounts$bg <- generateBackground(sample.names, exomecounts, median) 
  exomecounts$bg.sd <- generateBackground(sample.names, exomecounts, sd) 
}
