setMethod("plot", signature(x="ExomeCopy",y="missing"), plot.ExomeCopy)

setMethod("show","ExomeCopy",function(object) {
  cat("\nExomeCopy object\n")
  cat(paste("sample name: ",object@sample.name,"\n",sep=""))
  cat(paste("type: ",object@type,"\n",sep=""))
  normal.state <- object@fx.par$normal.state
  normal.percent <- round(100*sum(object@path==normal.state)/length(object@path),2)
  cat(paste("percent normal state: ",normal.percent,"%\n",sep=""))
  })
