\name{ExomeCopy-class}
\Rdversion{1.1}
\docType{class}
\alias{ExomeCopy-class}
\alias{plot,ExomeCopy,missing-method}
\alias{show,ExomeCopy-method}

\title{Class "ExomeCopy"}
\description{
  Object returned by exomeCopy
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ExomeCopy")}.
}
\section{Slots}{
  \describe{
    \item{\code{sample.name}:}{Object of class \code{"character"}: an
      identifier for the sample }
    \item{\code{type}:}{Object of class \code{"character"}: the type of
      model used, either
      "exomeCopy" or "exomeCopyVar" }
    \item{\code{path}:}{Object of class \code{"Rle"}: the index of
      the predicted state for each window }
    \item{\code{ranges}:}{Object of class \code{"IRangesList"}:
      the corresponding ranges for the observed counts and covariates }   
    \item{\code{O.norm}:}{Object of class \code{"numeric"}: the input
      vector of counts divided by X * beta }
    \item{\code{log.odds}:}{Object of class \code{"numeric"}: the log
      ratio of the emission probabilities along the predicted path over
      the emission probabilities of the normal state }
    \item{\code{fx.par}:}{Object of class \code{"list"}: a list of the
      settings S, d, normal.state and fit.var}   
  \item{\code{init.par}:}{Object of class \code{"list"}: a list of the initial parameters goto.cnv, goto.normal, beta.hat and phi.hat }
    \item{\code{final.par}:}{Object of class \code{"list"}: a list of the final parameters goto.cnv, goto.normal, beta (and gamma for exomeCopyVar) }
    \item{\code{counts}:}{Object of class \code{"numeric"}: the number
      of evaluations of the log likelihood performed by \code{\link{optim}} }
    \item{\code{convergence}:}{Object of class \code{"numeric"}: the
      integer for convergence of \code{\link{optim}}, 0 for convergence }
    \item{\code{nll}:}{Object of class \code{"numeric"}: the final value
    of the negative log likelihood }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "ExomeCopy", y = "missing")}: ... }
    \item{show}{\code{signature(object = "ExomeCopy")}: ... }
  }
}
\seealso{
  \code{\link{exomeCopy}}
}
\examples{
showClass("ExomeCopy")
}
\keyword{classes}