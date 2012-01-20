\name{countBamInGRanges}
\alias{countBamInGRanges}
\title{
  Count reads from BAM file in genomic ranges
}
\description{
  Counts the number of reads with a specified minimum mapping quality
  from BAM files in genomic ranges specified by a GRanges object.  This
  is a convenience function for counting the reads in ranges covering
  the targeted regions, such as the exons in exome enrichment
  experiments, from each sample.  These read counts are used by
  \code{\link{exomeCopy}} in predicting CNVs in samples.  

  With the default setting (\code{read.width}=1), only the read starts
  are used for counting purposes (the leftmost position regardless of
  the strandedness of the read). 

  With the accurate read width, or with \code{get.width} = TRUE, then
  the function returns the number of overlapping reads, as returned by 
  \code{\link{countOverlaps}} in the GenomicRanges package. 

  The function \code{\link{subdivideGRanges}} can be used first to
  subdivide ranges of different size into ranges of nearly equal width.

  The BAM file requires a associated index file (see the man page for
  \code{\link{indexBam}} in the Rsamtools package). 
}
\usage{
  countBamInGRanges(bam.file,granges,min.mapq=1,read.width=1,get.width=FALSE,remove.dup=FALSE)
}
\arguments{
  \item{bam.file}{
    The path of the BAM file for the sample to be counted.
  }
  \item{granges}{
    An object of type GRanges with the ranges in which to count reads.
  }
  \item{min.mapq}{
    The minimum mapping quality to count a read.  Defaults to 1.  Set to
  0 for counting all reads.
  }  
  \item{read.width}{
    The width of a read, used in counting overlaps of mapped reads with
    the genomic ranges.  The default is 1, resulting in the counting of
    only read starts in genomic ranges.  If the length of fixed width
    reads is used, e.g. 100 for 100bp reads, then the function will
    return the count of all overlapping reads with the genomic ranges.
    However, counting all overlapping reads introduces dependency
    between the counts in adjacent windows. 
  }
  \item{get.width}{
    Whether the function should retrieve the read width from the CIGAR
    encoding rather than assign the value from read.width.  If TRUE, then
    read.width is ignored.
  }
  \item{remove.dup}{
    If true, the function will count only one read for each unique
    combination of position, strand and read width. 
  }
}
\value{
  An integer vector giving the number of reads over the input GRanges
  }
\seealso{
  \code{\link{Rsamtools}}
  \code{\link{GRanges}}
  \code{\link{subdivideGRanges}}
}
\examples{
  ## get subdivided genomic ranges covering targeted region
  ## using subdivideGRanges()
  example(subdivideGRanges)

  ## BAM file included in Rsamtools package
  bam.file <- system.file("extdata", "mapping.bam", package="exomeCopy")

  ## create RangedData object to store read counts
  rdata <- RangedData(space=seqnames(target.sub),ranges=ranges(target.sub))

  ## extract read counts from the BAM file in these genomic ranges
  rdata[["sample"]] <- countBamInGRanges(bam.file,target.sub)
}