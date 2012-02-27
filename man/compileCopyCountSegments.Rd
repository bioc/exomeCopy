\name{compileCopyCountSegments}
\alias{compileCopyCountSegments}
\title{
  Compile segments across samples
}
\description{
  A short function which extracts the segments of constant copy count using the
  \code{\link{copyCountSegments}} function on a named list of named lists containing
  fitted ExomeCopy objects.  See vignette for a full example of multiple
  samples and chromosomes.
}
\usage{
  compileCopyCountSegments(fit.list)
}
\arguments{
  \item{fit.list}{
    A named list of named lists.  The outer list indexes patients, while
    the inner list indexes sequences/chromosomes.
  }
}
\value{
  A RangedData object of all segments across samples and chromosomes.
  }
\examples{
  example(exomeCopy)
  fit.list <- list(sample1 = list(chr1 = fit))
  CNV.segments <- compileCopyCountSegments(fit.list)
}