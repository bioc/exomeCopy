logistic <- function(x) {
  exp(x) / (1 + exp(x))
}

logit <- function(x) {
  log(x) - log(1-x)
}

getGCcontent <- function(target, reference.file) {
  x <- scanFa(reference.file, target)
  GC.count <- letterFrequency(x,"GC")
  all.count <- letterFrequency(x,"ATGC")
  as.vector(ifelse(all.count==0,NA,GC.count/all.count))
}

generateBackground <- function(sample.names, rdata, fn=median) {
  C <- as.data.frame(unlist(values(rdata)[,colnames(rdata) %in% sample.names]))
  C.norm <- sweep(C,2,colMeans(C),"/")
  apply(C.norm,1,fn)
}

copyCountSegments <- function(object) {
  x <- object@path
  changes <- which(x[-length(x)] != x[-1]) + 1
  range.start <- c(1,changes)
  range.end <- c(changes-1,length(x))
  targeted.bp <- sapply(1:length(range.start),function(i) sum(width(object@ranges)[[1]][range.start[i]:range.end[i]]))
  RangedData(
             IRanges(start=start(object@ranges)[[1]][range.start],
             end=end(object@ranges)[[1]][range.end]),
             space=space(object@ranges)[1],
             universe=universe(object@ranges),
             copy.count=object@fx.par$S[as.numeric(x[range.start])],
             nranges=(range.end - range.start + 1),
             targeted.bp=targeted.bp,
             sample.name=object@sample.name)
}

compileCopyCountSegments <- function(fit.list) {
  sample.res <- lapply(fit.list,function(seq.list) {
    seq.res <- lapply(seq.list, copyCountSegments)
    do.call(rbind,seq.res)
  })
  do.call(rbind,sample.res)
}

plot.ExomeCopy <- function (x,points=TRUE,cols=NULL,show.legend=TRUE,main="exomeCopy predicted segments",xlab="genomic position",ylab="normalized read count",xlim=NULL,ylim=NULL,cex=1,lwd=4,...) {
  if (length(cols) != length(x@fx.par$S)) {
    warning("Supplied colors are not the same length as S (copy number states)")
  }
  if (is.null(cols) | (length(cols) != length(x@fx.par$S))) {
    cols <- rep("black",length(x@fx.par$S))
    cols[x@fx.par$S < x@fx.par$d] <- "red"
    cols[x@fx.par$S > x@fx.par$d] <- "blue"
  }
  if (is.null(ylim)) ylim <- c(0,5*mean(x@O.norm))
  if (is.null(xlim)) xlim <- c(start(range(x@ranges))[[1]],end(range(x@ranges))[[1]])
  plot(0,0,type="n",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=main,...)
  if (points) {
    points(mid(x@ranges[[1]]),x@O.norm*x@fx.par$d,col=cols[as.numeric(x@path)],cex=cex)
  }
  ccs <- copyCountSegments(x)
  segments(start(ranges(ccs))[[1]],ccs$copy.count,
           end(ranges(ccs))[[1]],ccs$copy.count,lwd=lwd,
           col=cols[match(ccs$copy.count,x@fx.par$S)])
  if (show.legend) {
    legend("topright",col=rev(cols),pch=1,legend=rev(x@fx.par$S),title="copy count",bg="white",cex=.75)
  }
}

plotCompiledCNV <- function(CNV.segments, seq.name, xlim=NULL, col=NULL, copy.counts=0:6, normal.state = 2) {
  CNV.subset <- CNV.segments[CNV.segments$space == seq.name,]
  seq.range <- range(CNV.subset)[[seq.name]]
  if (is.null(xlim)) {
    xlim <- c(start(seq.range),end(seq.range))
  }
  if (is.null(col)) {
    col <- ifelse(copy.counts < normal.state,"red",ifelse(copy.counts > normal.state,"blue","black"))
  }
  samples <- unique(CNV.subset$sample.name)
  plot(0,0,type="n",xlim=xlim,ylim=c(0,length(samples)+1),xlab="genomic position",ylab="",main=seq.name,yaxt="n")
  axis(2,at=1:length(samples),labels=samples,las=2,cex.axis=.6)
  for (i in 1:length(samples)) {
    sample.name <- samples[i]
    CNV.sample <- CNV.subset[CNV.subset$sample.name == sample.name,]
    copy.count.idx <- match(CNV.sample$copy.count,copy.counts)  
    segments(start(CNV.sample),rep(i,length(CNV.sample)),end(CNV.sample),rep(i,length(CNV.sample)),lwd=3,col=col[copy.count.idx])
  }
}

countBamInGRanges <- function (bam.file, granges, min.mapq = 1, read.width = 1, stranded.start = FALSE, get.width = FALSE, remove.dup = FALSE) {
  rds.counts <- integer(length(granges))
  seq.names <- unique(as.character(seqnames(granges)))
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)
  for (seq.name in seq.names) {
    if (seq.name %in% seq.names.in.bam) {
      granges.subset <- granges[seqnames(granges) == seq.name]
      strand(granges.subset) <- "*"
      scan.what <- c("pos","mapq")
      if (stranded.start | get.width) {
        scan.what <- c(scan.what, "qwidth")
      }
      if (stranded.start | remove.dup) {
        scan.what <- c(scan.what, "strand")
      }
      rds <- scanBam(bam.file, param = ScanBamParam(what = scan.what, which = range(granges.subset)))
      if (min.mapq > 0) {
        mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
      } else {
        mapq.test <- rep(TRUE,length(rds[[1]]$mapq))
      }
      if (remove.dup) {
        if (get.width) {
          # this check is fast and accurate, assuming that read widths
          # are less than 1e4 bp and positions are less than 1e12 bp
          # the precision for doubles is 16 significant digits
          mapq.test <- mapq.test & !duplicated(rds[[1]]$pos + as.numeric(rds[[1]]$strand)/10 + rds[[1]]$qwidth/1e5)
        } else {
          mapq.test <- mapq.test & !duplicated(rds[[1]]$pos + as.numeric(rds[[1]]$strand)/10)
        }
      }
      if (sum(mapq.test) > 0) {
        if (stranded.start) {
          rds.ranges <- GRanges(seq.name, IRanges(start = ifelse(rds[[1]]$strand[mapq.test] %in% c("+","*"), rds[[1]]$pos[mapq.test] , rds[[1]]$pos[mapq.test] + rds[[1]]$qwidth[mapq.test] - read.width), end = ifelse(rds[[1]]$strand[mapq.test] %in% c("+","*"), rds[[1]]$pos[mapq.test] + read.width - 1, rds[[1]]$pos[mapq.test] + rds[[1]]$qwidth[mapq.test] - 1)))
        }
        else if (get.width) {
          rds.ranges <- GRanges(seq.name, IRanges(start = rds[[1]]$pos[mapq.test], width = rds[[1]]$qwidth[mapq.test]))
        }
        else {
          rds.ranges <- GRanges(seq.name, IRanges(start = rds[[1]]$pos[mapq.test], width = read.width))
        }
        rds.counts.seq.name <- countOverlaps(granges.subset, rds.ranges)
        rds.counts[as.logical(seqnames(granges) == seq.name)] <- rds.counts.seq.name
      }
      else {
        rds.counts[as.logical(seqnames(granges) == seq.name)] <- 0
      }
    }
    else {
      rds.counts[as.logical(seqnames(granges) == seq.name)] <- 0
    }
  }
  if (sum(rds.counts) == 0) {
    warning("No reds found with minimum mapping quality")
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
  nsubranges <- pmax(1,floor(widths/subsize))
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

subdivideGRanges <- function (x, subsize=100) {
  if (length(x) == 0) {
    return(x)
  }
  gr_list <- lapply(levels(seqnames(x)), function(seqlvl) {
    if (!any(seqnames(x) == seqlvl)) {
      return(GRanges())
    }
    rg <- ranges(x[seqnames(x) == seqlvl])
    GRanges(seqnames = seqlvl, ranges = subdivideIRanges(rg,subsize), seqlengths = seqlengths(x))
  })
  do.call(c, gr_list)
}

