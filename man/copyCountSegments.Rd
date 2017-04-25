\name{copyCountSegments}
\alias{copyCountSegments}
\title{
  Segments of identical copy count from exomeCopy
}
\description{
  Unpacks an ExomeCopy object and returns a GRanges object with
  segments of identical predicted copy count in genomic coordinates.
}
\usage{
  copyCountSegments(object)
}
\arguments{
  \item{object}{
    ExomeCopy object
  }
}
\details{
  The log odds column is calculated by summing the log ratios over the
  contained ranges.  The log ratios at each range is the log of the
  emission probability for the given read count for the predicted state
  divided by the emission probability for the normal state.  The higher
  the value, the more likely that the read counts in this range could
  not have been generated from the normal state.
}
\value{
  Returns a GRanges object with the predicted copy count, the log
  odds of predicted copy count over normal copy count, a combined
  p-value, the number of genomic ranges spanned by the segment, the
  number of targeted basepairs in the segment, and the sample name. 
}
\seealso{
  \code{\link{exomeCopy}}
  \code{\link{ExomeCopy-class}}
}
\examples{

example(exomeCopy)
copyCountSegments(fit)

}
