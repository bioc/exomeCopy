setClass("ExomeCopy",
         representation =
         representation(sample.name="character",
                        type="character",
                        path="Rle",
                        ranges="GRanges",
                        O.norm="numeric",
                        log.odds="numeric",                       
                        fx.par="list",
                        init.par="list",
                        final.par="list",
                        counts="numeric",
                        convergence="numeric",
                        nll="numeric"))


