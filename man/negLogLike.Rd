\name{negLogLike}
\alias{negLogLike}
\alias{viterbiPath}
\title{
  Generalized negative log likelihood and Viterbi algorithms
}
\description{
  negLogLike: Returns the negative log likelihood calculated with the
  forward equations.
  
  viterbiPath: Calculates the most likely sequence of hidden states for the
  Markov model given the current parameters. 
}
\usage{
  negLogLike(par,fx.par,data,nstates,stFn,trFn,emFn)
  viterbiPath(par,fx.par,data,nstates,stFn,trFn,emFn)
}
\arguments{
  \item{par}{
    A list of parameters, over which the likelihood will be optimized.
  }
  \item{fx.par}{
    A list of fixed parameters.
  }
  \item{data}{
    A list of data objects, which must contain a vector O, which
    represents the observed sequence of the HMM.
  }
  \item{nstates}{
    The number of states of the HMM.
  }
  \item{stFn}{
    A function which takes arguments par, fx.par, data, and nstates, and
    returns a vector of length nstates of starting probabilities.
  }
  \item{trFn}{
    A function which takes arguments par, fx.par, data, and nstates, and
    returns a matrix of dimension (nstates,nstates) of the transition
    probabilities. 
  }
  \item{emFn}{
    A function which takes arguments par, fx.par, data, and nstates, and
    returns a matrix of dimension (nstates,length(O)) of the emission
    probabilities. 
  }
}
\value{
  negLogLike: The negative log likelihood of the HMM.  The likelihood is
  slightly modified to account for ranges with read counts which have
  zero probability of originating from any of the states.  In this case
  the likelihood is lowered and the range is skipped.  
  
  viterbiPath: The Viterbi path through the states given the parameters. 
}
\references{
  On the forward equations and the Viterbi algorithm:
  
  Rabiner, L. R. (1989): "A tutorial on hidden Markov models and
  selected applications in speech recognition," Proceedings of the
  IEEE, 77, 257, 286, \url{http://dx.doi.org/10.1109/5.18626}.

}
\examples{
## functions for starting, transition, and emission probabilities
stFn <- function(par,fx.par,data,nstates) rep(1/nstates,nstates)
trFn <- function(par,fx.par,data,nstates) {
  A <- matrix(1/(nstates*10),ncol=nstates,nrow=nstates)
  diag(A) <- 1 - rowSums(A)
  A
}
emFn <- function(par,fx.par,data,nstates) {
  t(sapply(1:nstates,function(j) dnorm(data$O,par$means[j],fx.par$sdev)))
}

## simulate some observations from two states
Q <- c(rep(1,100),rep(2,100),rep(1,100),rep(2,100))
T <- length(Q)
means <- c(-0.5,0.5)
sdev <- 1
O <- rnorm(T,means[Q],sdev)

## use viterbiPath() to recover the state chain using parameters
viterbi.path <- viterbiPath(par=list(means=means),fx.par=list(sdev=sdev),data=list(O=O),nstates=2,stFn,trFn,emFn)
plot(O,pch=Q,col=c("darkgreen","orange")[viterbi.path])
}