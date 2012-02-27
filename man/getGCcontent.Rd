\name{getGCcontent}
\alias{getGCcontent}
\title{
  Get the GC content of target ranges from a reference FASTA file
}
\description{
  A short function using scanFa from the Rsamtools package on a target GRanges and a reference FASTA file
}
\usage{
  getGCcontent(target, reference.file)
}
\arguments{
  \item{target}{
    GRanges object
  }
  \item{reference.file}{
    the path to the reference FASTA file
  }
}
\value{
  Returns a vector of the ratio of G and C basepairs to total basepairs (not
  counting N's).
}
\examples{
  target.file <- system.file("extdata","targets.bed",package="exomeCopy")
  target.df <- read.delim(target.file,header=FALSE,col.names=c("seqname","start","end")) 
  target <- GRanges(seqname=target.df$seqname,IRanges(start=target.df$start+1,end=target.df$end))
  reference.file <- system.file("extdata","reference.fa",package="exomeCopy")
  GCcontent <- getGCcontent(target, reference.file)
}