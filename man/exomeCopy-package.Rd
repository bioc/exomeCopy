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
    Version: \tab 1.1.2\cr
    Date: \tab 2011-11-23\cr
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
  Love, Michael I.; Mysickova, Alena; Sun, Ruping; Kalscheuer, Vera;
  Vingron, Martin; and Haas, Stefan A. (2011) "Modeling Read Counts for
  CNV Detection in Exome Sequencing Data," Statistical Applications in
  Genetics and Molecular Biology: Vol. 10 : Iss. 1, Article 52. DOI: 10.2202/1544-6115.1732
  \url{http://www.bepress.com/sagmb/vol10/iss1/art52}.
}
\keyword{ package }
\seealso{
  \code{\link{exomeCopy}}
}