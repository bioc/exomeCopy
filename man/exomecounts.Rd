\name{exomecounts}
\docType{data}
\alias{exomecounts}
\title{Sample counts from 16 exome sequencing samples from 1000 Genomes Project}
\description{
  This data set gives sample read counts in 1000 genomic ranges for 16
  exome sequencing samples from the PUR population of the 1000 Genomes
  Project, along with the GC-content in the ranges.  For instructions on
  how to prepare read count and covariate data, please see the example
  code in the man pages for \code{\link{subdivideGRanges}} and
  \code{\link{countBamInGRanges}}. 

  The genomic ranges are generated from small portion of the CCDS regions of
  chromosome 1 (hg19).  The CCDS regions are subdivided evenly into
  ranges around 100bp using the \code{\link{subdivideGRanges}} function
  with default settings. Only ranges with positive counts across samples
  are retained.  These regions were downloaded as a BED file from the
  UCSC Genome Browser
  (\url{http://genome.ucsc.edu/cgi-bin/hgGateway}). The mapping files 
  for the exome sequencing data and descriptions of the experiments are
  available at the 1000 Genomes Project website
  (\url{http://www.1000genomes.org/data}).  The directories used are
  listed in the file \code{1000Genomes_files.txt} in the \code{extdata} directory.

  The column names are the sample names from the 1000 Genomes Project.
  Library format is paired-end reads and sample counts reflect both
  sequenced reads counted in their respective genomic ranges.  
} 
\usage{data(exomecounts)}
\format{A RangedData object.}
\source{1000 Genomes Project and Consensus Coding Sequence Project}
\references{
  1000 Genomes Project Consortium. A map of human genome variation from
  population-scale sequencing. Nature 467, 1061-1073 (2010). 
  \url{http://dx.doi.org/10.1038/nature09534}.

  1000 Genomes Project: Release of phase 1 exome alignments
  \url{http://www.1000genomes.org/announcements/release-phase-1-exome-alignments-2011-07-19}

  Pruitt, K. D. et al. The consensus coding sequence (CCDS) project:
  Identifying a common protein-coding gene set for the human and mouse
  genomes. Genome research 19, 1316-1323 (2009). 
  \url{http://dx.doi.org/10.1101/gr.080531.108}.
}
\keyword{datasets}
\seealso{
  \code{\link{RangedData}}
}