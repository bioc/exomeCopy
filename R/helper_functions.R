logistic <- function(x) {
  exp(x) / (1 + exp(x))
}

logit <- function(x) {
  log(x) - log(1-x)
}

copyCountSegments <- function(object) {
  x <- object@path
  changes <- which(x[-length(x)] != x[-1]) + 1
  range.start <- c(1,changes)
  range.end <- c(changes-1,length(x))
  RangedData(
             IRanges(start=start(object@ranges)[[1]][range.start],
             end=end(object@ranges)[[1]][range.end]),
             space=space(object@ranges)[1],
             universe=universe(object@ranges),
             copy.count=object@fx.par$S[x[range.start]],
             nranges=(range.end - range.start + 1))
}

plot.ExomeCopy <- function (x,points=TRUE,cols=NULL,show.legend=TRUE,main="exomeCopy predicted segments",xlab="genomic position",ylab="normalized read count",xlim=NULL,ylim=NULL,cex=1,lwd=4,...) {
  if (is.null(cols)) {
    cols <- rep("black",length(x@fx.par$S))
    cols[x@fx.par$S < x@fx.par$d] <- "red"
    cols[x@fx.par$S > x@fx.par$d] <- "blue"
  }
  if (is.null(ylim)) ylim <- c(0,5*mean(x@O.norm))
  if (is.null(xlim)) xlim <- c(as.numeric(start(range(x@ranges))),as.numeric(end(range(x@ranges))))
  plot(0,0,type="n",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=main,...)
  if (points) {
    points(mid(x@ranges[[1]]),x@O.norm*x@fx.par$d,col=cols[x@path],cex=cex)
  }
  ccs <- copyCountSegments(x)
  segments(as.numeric(start(ranges(ccs))),ccs$copy.count,
           as.numeric(end(ranges(ccs))),ccs$copy.count,lwd=lwd,col=cols[match(ccs$copy.count,x@fx.par$S)])
  if (show.legend) {
    legend("topright",col=rev(cols),pch=1,legend=rev(x@fx.par$S),title="copy count",bg="white",cex=.75)
  }
}

countBamInGRanges <- function(bam.file,granges,min.mapq=1,read.width=1) {
  rds.counts <- numeric(length(granges))
  seq.names <- as.character(unique(seqnames(granges)))
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)
  for (seq.name in seq.names) {
    if (seq.name %in% seq.names.in.bam) {
      granges.subset <- granges[seqnames(granges)==seq.name]
      strand(granges.subset) <- "*"
      rds <- scanBam(bam.file,param=ScanBamParam(what=c("pos","mapq"),which=range(granges.subset)))
      mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
      if (sum(mapq.test) > 0) {
        rds.ranges <- GRanges(seq.name,IRanges(start=rds[[1]]$pos[mapq.test],width=read.width))
        rds.counts.seq.name <- countOverlaps(granges.subset,rds.ranges)
        rds.counts[as.logical(seqnames(granges)==seq.name)] <- rds.counts.seq.name
      } else {
        rds.counts[as.logical(seqnames(granges)==seq.name)] <- 0
      }
    } else {
      rds.counts[as.logical(seqnames(granges)==seq.name)] <- 0
    }
  }
  rds.counts
}

subdivide.range <- function(start.pos,end.pos,subsize=100) {
  width <- end.pos - start.pos + 1
  if (width < 2*subsize) {
    stop("Width is less than 2 times subsize")
  }
  nchunks <- floor(width/subsize)
  relative.starts <- round(0:(nchunks-1)*width/nchunks)
  relative.ends <- round(1:nchunks*width/nchunks)-1
  return(list(start.pos=start.pos+relative.starts,end.pos=start.pos+relative.ends,length=nchunks))
}

subdivideIRanges <- function(x,subsize=100) {
  if (length(x) == 0) {
    return(x)
  }
  start.pos <- start(x)
  end.pos <- end(x)
  widths <- width(x)
  nsubranges <- pmax(1,floor(widths/100))
  out.start.pos <- numeric(sum(nsubranges))
  out.end.pos <- numeric(sum(nsubranges))
  out.idx <- 1
  for (i in 1:length(x)) {
    if (widths[i] < 2*subsize) {
      out.start.pos[out.idx] <- start.pos[i]
      out.end.pos[out.idx] <- end.pos[i]
      out.idx <- out.idx + 1
    } else {
      sr <- subdivide.range(start.pos[i],end.pos[i],subsize)
      out.start.pos[out.idx:(out.idx+sr$length-1)] <- sr$start.pos
      out.end.pos[out.idx:(out.idx+sr$length-1)] <- sr$end.pos
      out.idx <- out.idx + sr$length
    }
  }
  IRanges(start=out.start.pos,end=out.end.pos)
}

subdivideGRanges <- function(x,subsize=100) {
  if (length(x) == 0) {
    return(x)
  }
  gr_list <- lapply(levels(seqnames(x)),
                    function(seqlvl) {
                      rg <- ranges(x[seqnames(x) == seqlvl])
                      GRanges(seqnames=seqlvl,
                              ranges=subdivideIRanges(rg,subsize),seqlengths=seqlengths(x))
                    })
  do.call(c, gr_list)
}