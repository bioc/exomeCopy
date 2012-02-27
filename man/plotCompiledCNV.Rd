\name{plotCompiledCNV}
\alias{plotCompiledCNV}
\title{
  Plot compiled CNV segments for one sequence/chromosome
}
\description{
  This function takes a RangedData object as produced by
  \code{\link{compileCopyCountSegments}} and plots the CNV segments for
  one sequence/chromosomes across the samples with CNV segments.

  The segments in the normal state should be removed as shown below in
  the example to produce a cleaned RangedData object.  See the vignette
  for a more complete example.
}
\usage{
  plotCompiledCNV(CNV.segments, seq.name, xlim=NULL, col=NULL,
  copy.counts=0:6, normal.state = 2)
}
\arguments{
  \item{CNV.segments}{
    A RangedData object as produced by
  \code{\link{compileCopyCountSegments}} and with normal state removed.
  }
  \item{seq.name}{
    The name of the sequence to plot
  }
  \item{xlim}{
    The genomic coordinates for the x axis.  If not included, the
  plotting window will cover the range of the CNVs in \code{CNV.segments}
  }
  \item{col}{
    The colors to use for the different copy count states
  }
  \item{copy.counts}{
    The corresponding copy counts for the colors
  }
  \item{normal.state}{
    The copy count of the normal state
  }
}
\value{
  Produces a plot.
  }
\examples{
  example(compileCopyCountSegments)
  CNV.clean <- CNV.segments[CNV.segments$copy.count != 2,]
  chr.start <- start(range(fit@ranges)[["chr1"]])
  chr.end <- end(range(fit@ranges)[["chr1"]])
  plotCompiledCNV(CNV.clean, "chr1", xlim=c(chr.start,chr.end))
}