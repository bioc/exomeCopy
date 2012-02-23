\name{copyCountSegments}
\alias{copyCountSegments}
\title{
  Segments of identical copy count from exomeCopy
}
\description{
  Unpacks an ExomeCopy object and returns a RangedData object with
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
\value{
  Returns a RangedData object with the predicted copy count and
  the number of genomic ranges spanned by the segment.
}
\seealso{
  \code{\link{exomeCopy}}
  \code{\link{ExomeCopy-class}}
}
\examples{

example(exomeCopy)
copyCountSegments(fit)

}