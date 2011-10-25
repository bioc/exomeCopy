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

  With the default setting (\code{read.width}=1), only the read starts are used
  for counting purposes (the leftmost position regardless of the
  strandedness of the read).

  The function \code{\link{subdivideGRanges}} can be used first to
  subdivide ranges of different size into ranges of nearly equal width.
}
\usage{
  countBamInGRanges(bam.file,granges,min.mapq=1,read.width=1)
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
}
\value{
  A vector giving the number of reads over the input GRanges
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

  ## finally, we can also calculate and store GC-content for exomeCopy
  reference.file <- system.file("extdata", "reference.fa", package="exomeCopy")
  target.dnastringset <- scanFa(reference.file,target.sub)
  getGCcontent <- function(x) {
    GC.count <- letterFrequency(x,"GC")
    all.count <- letterFrequency(x,"ATGC")
    as.vector(ifelse(all.count==0,NA,GC.count/all.count))
  }
  rdata[["GC"]] <- getGCcontent(target.dnastringset)
}