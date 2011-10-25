\name{plot.ExomeCopy}
\alias{plot.ExomeCopy}
\title{
  Plot function for exomeCopy
}
\description{
  Plots the predicted copy count segments of an ExomeCopy object
}
\usage{
plot.ExomeCopy(x, points = TRUE, cols = NULL, show.legend = TRUE,
main = "exomeCopy predicted segments", xlab = "genomic position",
ylab = "normalized read count", xlim = NULL, ylim = NULL, cex = 0.5, lwd = 4, ...)
}
\arguments{
  \item{x}{
    The ExomeCopy object.
  }
  \item{points}{
    Logical, whether normalized read counts should be drawn.
  }
  \item{cols}{
    A vector of the same length as b, specifying a color for each of the
    states of the HMM.
  }
  \item{show.legend}{
    Logical, whether a default legend should be shown.
  }
  \item{main}{
    main title
  }
  \item{xlab}{
    x axis label
  }
  \item{ylab}{
   y axis label
  }
  \item{xlim}{
    x limits
  }
  \item{ylim}{
   y limits
  }
  \item{cex}{
    size of the points (if plotted)
  }
  \item{lwd}{
    line width
  }
  \item{...}{
    Other arguments passed to plot()
  }
}
\seealso{
  \code{\link{exomeCopy}}
  \code{\link{ExomeCopy-class}}
  \code{\link{copyCountSegments}}
}
\examples{

example(exomeCopy)
plot(fit)

}