\name{exomeCopy-package}
\alias{exomeCopy-package}
\alias{exomeCopy-package}
\docType{package}
\title{
  Detection of CNV in exome/targeted sequencing data
}
\description{
  A hidden Markov model for the detection of copy number variants (CNV) in exome/targeted sequencing read depth data.  The package uses positional covariates, such as background read depth and GC-content, to simultaneously normalize and segment the samples into regions of constant copy count.
}
\details{

  \tabular{ll}{
    Package: \tab exomeCopy\cr
    Type: \tab Package\cr
    Version: \tab 0.99.3\cr
    Date: \tab 2011-10-27\cr
    License: \tab GPL (>= 2)\cr
    LazyLoad: \tab yes\cr
    Depends: \tab methods, graphics, IRanges, GenomicRanges, Rsamtools
    (>= 1.4.3) \cr
    Suggests: \tab Biostrings\cr
  }

  exomeCopy fits a hidden Markov model to observed read counts using
  covariates.  It returns the Viterbi path, the most likely path of
  hidden states, which is the predicted copy count at each window. 
}
\author{
  Michael Love <love@molgen.mpg.de>
}
\references{
  "Modeling read counts for CNV detection in exome sequencing data",
  Michael I. Love, Alena Mysickova, Ruping Sun, Vera Kalscheuer, Martin
  Vingron, Stefan A. Haas.  Under review. 
}
\keyword{ package }
\seealso{
  \code{\link{exomeCopy}}
}