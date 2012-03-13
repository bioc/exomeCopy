\name{exomeCopy}
\alias{exomeCopy}
\title{
  Fit the exomeCopy or exomeCopyVar model to the observed counts.
}
\description{
  Fits a hidden Markov model to observed read counts using positional
  covariates. It returns an object containing the fitted parameters and
  the Viterbi path, the most likely path of hidden states, which is the
  predicted copy count at each window.  

  \code{exomeCopy} is designed to run on read counts from consecutive
  genomic ranges on a single chromosome, as it tries to identify higher
  or lower read depth relative to a baseline.  Please see the vignette
  for an example of how to prepare input data for \code{exomeCopy}, how
  to loop the function over multiple chromosomes and samples, and how to
  extract the resulting predicted CNVs.  

  \code{exomeCopy} requires as input a \code{\link{RangedData}} object
  containing read counts in genomic ranges along with the covariates.
  Some convenience functions are provided for preparing input for
  \code{exomeCopy}: 

  \enumerate{
  \item \code{\link{subdivideGRanges}}, to subdivide a \code{\link{GRanges}}
  object containing the genomic ranges of the targeted region into
  genomic ranges of nearly equal width,
  \item \code{\link{countBamInGRanges}}, to count the number of read starts
  from a BAM read mapping file in a \code{\link{GRanges}} object,
  \item \code{\link{getGCcontent}}, to get the GC-content in the ranges
  given a FASTA file of the reference sequence,
  \item \code{\link{generateBackground}}, to calculate median normalized
  read depth over a set of control samples, and also any statistic over
  normalized read depth.
}
}
\usage{
  exomeCopy(rdata, sample.name, X.names, Y.names, fit.var=FALSE, reltol
  = 0.0001, S = 0:4, d = 2, goto.cnv = 1e-4, goto.normal = 1/20,
init.phi="norm")    
}
\arguments{
  \item{rdata}{
    A RangedData object with the sample counts and positional covariates
    over the genomic ranges.
  }
  \item{sample.name}{
    The name of the value column of rdata with the sample read counts.
  }
  \item{X.names}{
    The names of the value columns of rdata with covariates for
  estimating mu. 
  }
  \item{Y.names}{
    (optional) the names of the value columns of rdata with covariates for
  estimating phi, only required if fit.var = TRUE. 
  }
  \item{fit.var}{
    A logical, whether the model should fit the overdispersion parameter
  phi with a linear combination of covariates (exomeCopyVar) or with a
  scalar (exomeCopy).  Defaults to FALSE (exomeCopy). 
  }
  \item{reltol}{
    The relative tolerance for convergence used in the \code{\link{optim}} function
  for optimizing the parameter settings.  From testing, the default
  value was sufficient for fitting parameters, but lower relative
  tolerances can be used.
  }
  \item{S}{
    A vector of possible copy numbers for the different states.
  }
  \item{d}{
    The expected copy number for the normal state.  This should be set
    to 2 for autosomes and 1 for haploid data.
  }
  \item{goto.cnv}{
    The initial setting for probability to transfer to a CNV state.
  }
  \item{goto.normal}{
    The initial setting for probability to transfer to a normal state.
  }
  \item{init.phi}{
    Either "norm" or "counts": initialize phi with the moment estimate
    using residuals from a linear model of read counts on covariates or
    with the raw counts.
  }
}
\details{
  \code{exomeCopy} fits transitional and emission parameters of an HMM to best
  explain the observed counts of a sample from exome or targeted
  sequencing.  The set of underlying copy number states, S, in the
  sample must be provided before running the algorithm.

  The emission probabilities are given as a negative binomial
  distribution using positional covariates, such as background read
  depth, quadratic terms for GC-content, and range width, which are
  stored in a matrix \eqn{X}.  Optionally, for fitting the variance of the
  distribution, the standard deviation and/or variance of the background
  set can be included in a matrix \eqn{Y}.  All covariates are normalized
  within \code{exomeCopy} for improved optimization.
  
  For the observed count at range \eqn{t}, \eqn{O_t}, the emission
  probability is given by:

  \deqn{f \sim  \textrm{NB}(O_t, \mu_{ti}, \phi)}{f = NB(O_t, mu_ti, phi)}

  The mean parameter \eqn{\mu_{ti}}{mu_ti} is given by:

  \deqn{ \mu_{ti} = \frac{S_i}{d} ( x_{t*} \vec{\beta} ) }{ mu_ti
  = (S_i / d)(x_t* beta)}

  Here \eqn{S_i} is the i-th possible copy number state, \eqn{d} is the
  expected background copy number (d = 2 for diploid sequence), and
  \eqn{\vec{\beta}}{beta} is a vector of coefficients fitted by the
  model.  \eqn{x_{t*}}{x_t*} is the t-th row of the matrix \eqn{X}.

  \eqn{\mu}{mu} must be positive, so it is replaced with a small positive
  number if the value is less than zero.

  For exomeCopyVar, which also fits the variance, the emission
  probability is given by:

  \deqn{f \sim  \textrm{NB}(O_t, \mu_{ti}, \phi_{t})}{f = NB(O_t, mu_ti,
  phi_t)}

  where

  \deqn{ \phi_{t} = y_{t*} \vec{\gamma} }{ phi_t = y_t* gamma }

  or a small positive number if this is less than zero.
  
  Two transition probabilities are fitted in the model: the
  probabilities of transitioning to a normal state and to a CNV state.

  \code{exomeCopy} calls \code{\link{negLogLike}} to evaluate the likelihood of
  the HMM.  The parameters are fit using Nelder-Mead optimization with
  the \code{\link{optim}} function on the negative likelihood. The
  Viterbi path is calculated by calling \code{\link{viterbiPath}}.  
}
\value{
  Returns an \code{\link{ExomeCopy-class}} object.  See this page for the
  slot descriptions.  Also see the vignette and
  \code{\link{copyCountSegments}} on how to extract segments. 
}
\references{
  Love, Michael I.; Mysickova, Alena; Sun, Ruping; Kalscheuer, Vera;
  Vingron, Martin; and Haas, Stefan A. (2011) "Modeling Read Counts for
  CNV Detection in Exome Sequencing Data," Statistical Applications in
  Genetics and Molecular Biology: Vol. 10 : Iss. 1, Article 52. DOI: 10.2202/1544-6115.1732
  \url{http://www.degruyter.com/view/j/sagmb.2011.10.issue-1/1544-6115.1732/1544-6115.1732.xml}.

  References for HMM algorithms and use of HMM for segmentation of
  genomic data by copy number:

  Rabiner, L. R. (1989): "A tutorial on hidden Markov models and
  selected applications in speech recognition," Proceedings of the
  IEEE, 77, 257, 286, \url{http://dx.doi.org/10.1109/5.18626}.

  Fridlyand, J., A. M. Snijders, D. Pinkel, D. G. Albertson, and Jain
  (2004): "Hidden Markov models approach to the analysis of array CGH
  data," Journal of Multivariate Analysis, 90, 132, 153, 
  \url{http://dx.doi.org/10.1016/j.jmva.2004.02.008}.

  Marioni, J. C., N. P. Thorne, and S. Tavare (2006):"BioHMM: a
  heterogeneous hidden Markov model for segmenting array CGH data."
  Bioinformatics, 22, 1144, 1146, 
  \url{http://view.ncbi.nlm.nih.gov/pubmed/16533818}. 
}
\author{
  Michael Love
}
\seealso{
  \code{\link{ExomeCopy-class}}
  \code{\link{subdivideGRanges}}
  \code{\link{countBamInGRanges}}
  \code{\link{copyCountSegments}}
  \code{\link{plot.ExomeCopy}}
  \code{\link{negLogLike}}
  \code{\link{IRanges}}
  \code{\link{RangedData}}
}
\examples{

## The following is an example of running exomeCopy on simulated
## read counts using the model parameters defined above.  For an example
## using real exome sequencing read counts (with simulated CNV) please
## see the vignette.

## create RangedData for storing genomic ranges and covariate data
## (background, background stdev, GC-content)

m <- 5000
rdata <- RangedData(IRanges(start=0:(m-1)*100+1,width=100),
space=rep("chr1",m), universe="hg19", bg=rexp(m,1), bg.sd=rexp(m,1),
gc=rnorm(m,50,10))  

## create read depth distributional parameters mu and phi
rdata$gc.sq <- rdata$gc^2
X <- cbind(bg=rdata$bg,gc=rdata$gc,gc.sq=rdata$gc.sq)
Y <- cbind(bg.sd=rdata$bg.sd)
beta <- c(20,10,2,-.01)
gamma <- c(.1,.05)
rdata$mu <- beta[1] + scale(X) \%*\% beta[2:4]
rdata$mu[rdata$mu<1e-8] <- 1e-8
rdata$phi <- gamma[1] + scale(Y) \%*\% gamma[2]
rdata$phi[rdata$phi<1e-8] <- 1e-8

## create observed counts with simulated heterozygous duplication
cnv.nranges <- 200
bounds <- (round(m/2)+1):(round(m/2)+cnv.nranges)
O <- rnbinom(nrow(rdata),mu=rdata$mu,size=1/rdata$phi)
O[bounds] <- O[bounds] + rbinom(cnv.nranges,prob=0.5,size=O[bounds])
rdata[["sample1"]] <- O

## run exomeCopy() and list segments
fit <- exomeCopy(rdata,"sample1",X.names=c("bg","gc","gc.sq"))

## see man page for copyCountSegments() for summary of
## the predicted segments of constant copy count, and
## for plot.ExomeCopy() for plotting fitted objects

}
