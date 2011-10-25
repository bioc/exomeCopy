\name{subdivideGRanges}
\alias{subdivideGRanges}
\title{
  Subdivide ranges of a GRanges object into nearly equal width ranges
}
\description{
  Takes an input GRanges object and, splits each range into multiple
  ranges of nearly equal width.  For an input range of width \code{w} and
  subdividing size \code{s}, it will subdivide the range into
  \code{max(1,floor(w/s))} nearly equal width ranges.   The output is
  then a new GRanges object.  This function can be used to split the targeted region (such as exons in exome enrichment experiments) into nearly equal width ranges.
}
\usage{
  subdivideGRanges(x,subsize=100)
}
\arguments{
  \item{x}{
    An object of type GRanges.
  }
  \item{subsize}{
    The desired width for the ranges in the output GRanges object.
  }
}
\value{
  A GRanges object with ranges from the input GRanges object subdivided
  to nearly \code{subsize}.
  }
\seealso{
  \code{\link{GRanges}}
}
\examples{
  ## read in target region BED file
  target.file <- system.file("extdata", "targets.bed", package="exomeCopy")
  target.df <- read.delim(target.file, header=FALSE,
col.names=c("seqname","start","end")) 

  ## create GRanges object with 5 ranges over 2 sequences
  target <- GRanges(seqname=target.df$seqname,
               IRanges(start=target.df$start,end=target.df$end))

  ## subdivide into 7 smaller genomic ranges
  target.sub <- subdivideGRanges(target)
}