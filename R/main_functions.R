stFn <- function(par,fx.par,data,nstates) {
  # this is a quick fix to cap transition probabilities near 0 and 1
  par[1:2] <- pmin(pmax(par[1:2],-100),100)


  goto.cnv <- logistic(par[1])
  goto.normal <- logistic(par[2])
  normal.state <- fx.par$normal.state
  min.normal.stay <- .5
  start.probs <- numeric(nstates)
  if (goto.cnv < 1/(nstates-1)) {
        start.probs[-normal.state] <- goto.cnv
      } else {
        start.probs[-normal.state] <- (1 - min.normal.stay)/(nstates-1)
      }
  start.probs[normal.state] <- 1 - sum(start.probs[-normal.state])
  return(start.probs)
}

trFn <- function(par,fx.par,data,nstates) {
  # this is a quick fix to cap transition probabilities near 0 and 1
  par[1:2] <- pmin(pmax(par[1:2],-100),100)

  goto.cnv <- logistic(par[1])
  goto.normal <- logistic(par[2])
  normal.state <- fx.par$normal.state
  min.normal.stay <- .5
  A <- matrix(1e-8,ncol=nstates,nrow=nstates)
  for (i in 1:nstates) {
    if (i != normal.state) {
      if (i+1 <= nstates) {
        A[i,i+1] <- goto.cnv
      }
      if (i-1 >= 1) {
        A[i,i-1] <- goto.cnv
      }
      A[i,normal.state] <- goto.normal
    }
    if (i == normal.state) {
      if (goto.cnv < 1/(nstates-1)) {
        A[i,-i] <- goto.cnv
      } else {
        A[i,-i] <- (1 - min.normal.stay)/(nstates-1)
      }
    }
    A[i,i] <- 1 - sum(A[i,-i])
  }
  return(A)
}


emFn <- function(par,fx.par,data,nstates) {
  S <- fx.par$S
  d <- fx.par$d
  fit.var <- fx.par$fit.var
  beta <- par[3:(2+ncol(data$X))]
  mu <- data$X %*% beta
  mu[mu < 1] <- 1
  if (!fit.var) {
    phi <- exp(par[(3+ncol(data$X))])
    phi <- max(1e-6,phi)
    phi <- min(phi,1e4)
  } else {
    gamma <- par[(3+ncol(data$X)):length(par)]
    phi <- data$Y %*% gamma
    phi[phi < 1e-6] <- 1e-6
    phi[phi > 1e4] <- 1e4
  }
  emit.probs <- t(sapply(1:nstates,function(j) dnbinom(data$O,mu=(mu*S[j]/d + ifelse(S[j]==0,1,0)),size=1/phi)))
  return(emit.probs)
}

negLogLike <- function(par,fx.par,data,nstates,stFn,trFn,emFn) {
  start.probs <- stFn(par,fx.par,data,nstates)
  A <- trFn(par,fx.par,data,nstates)
  emit.probs <- emFn(par,fx.par,data,nstates)
  T <- length(data$O)
  if (length(start.probs) != nstates) {
    stop("vector of starting probabilities not equal to number of states")
  }
  if (all(dim(A) != nstates)) {
    stop("transition matrix must have the same number of rows and columns as number of states")
  }
  if ((nrow(emit.probs) != nstates) | (ncol(emit.probs) != T)) {
    stop("emission probabilities matrix must have a row for each state and a column for each position")
  }
  negloglike.call <- .C("negloglike",tmax=as.integer(T),nstates=as.integer(nstates),start.probs=as.double(start.probs),A=as.double(A),emit.probs=as.double(emit.probs),alpha=double(nstates),alpha.new=double(nstates),nll=as.double(1))
  return(negloglike.call$nll)
}

viterbiPath <- function(par,fx.par,data,nstates,stFn,trFn,emFn) {
  start.probs <- stFn(par,fx.par,data,nstates)
  A <- trFn(par,fx.par,data,nstates)
  emit.probs <- emFn(par,fx.par,data,nstates)
  T <- length(data$O)
  if (length(start.probs) != nstates) {
    stop("vector of starting probabilities not equal to number of states")
  }
  if (all(dim(A) != nstates)) {
    stop("transition matrix must have the same number of rows and columns as number of states")
  }
  if ((nrow(emit.probs) != nstates) | (ncol(emit.probs) != T)) {
    stop("emission probabilities matrix must have a row for each state and a column for each position")
  }
  V <- matrix(0,nrow=nstates,ncol=T)
  V.path <- matrix(0,nrow=nstates,ncol=T)
  viterbi.call <- .C("viterbi",tmax=as.integer(T),nstates=as.integer(nstates),start.probs=as.double(start.probs),A=as.double(A),emit.probs=as.double(emit.probs),V=as.double(V),V.path=as.integer(V.path),path=as.integer(numeric(T)),trans.prob=as.double(numeric(nstates^2)),trans.prob.max=as.double(numeric(nstates)),trans.prob.whichmax=as.integer(numeric(nstates)))
  return(viterbi.call$path + 1)
}

exomeCopy <- function(rdata, sample.name, X.names, Y.names, fit.var=FALSE, reltol=1e-4, S=0:4, d=2, goto.cnv=1e-4, goto.normal=1/20, init.phi="norm") {
  if (!sample.name %in% colnames(rdata)) {
    stop("sample.name is not a column in rdata.")
  }
  O <- rdata[[sample.name]]
  if (any(O != round(O) | O < 0)) {
    stop("Sample counts must be non-negative integers")
  }
  if (mean(O == 0) > 0.9) {
    warning("More than 90% of sample counts are zero, exomeCopy will return NULL")
    return(NULL)
  }
  if (any(S!=round(S)|S<0) | any(d!=round(d)|d<0)) {
    stop("S and d must be non-negative integers")
  }
  if (!all(d %in% S)) {
    stop("The normal state, d, must be one of the possible copy states in S")
  } 
  if (!(all(X.names %in% colnames(rdata)))) {
    stop("all X.names must be columns in rdata")
  }
  if (length(rdata) != 1) {
    stop("The current rdata argument has ranges over the following chromosomes/spaces: ",paste(names(rdata),collapse=", "),".\n The rdata argument should contain ranged data over a single space.\n You can pass exomeCopy a ranged data from a single space with this syntax: rdata['chr1'].\n See vignette for example of running on multiple chromosomes.")
  }
  if (nrow(rdata) < 100) {
    warning("exomeCopy was tested for thousands of ranges covering a targeted region on a single chromosome.  The results might not be reliable for less than a hundred ranges.")
  }
  if (is.unsorted(start(rdata))) {
    stop("Genomic ranges for ranged data must be sorted.")
  }
  normal.state = which(S==d)
  X <- as.matrix(as.data.frame(values(rdata))[,X.names])
  colnames(X) <- X.names
  if (fit.var) {
    if (!(all(Y.names %in% colnames(rdata)))) {
      stop("Y.names must be variable names in rdata")
    }
    Y <- as.matrix(as.data.frame(values(rdata))[,Y.names])
    colnames(Y) <- Y.names
  }
  controls <- list(reltol=reltol,maxit=10000)
  X.full <- cbind(intercept=rep(1,nrow(X)),scale(X))
  lmfit <- lm(O ~ X.full + 0)
  beta.hat <- lmfit$coefficients
  names(beta.hat) <- colnames(X.full)
  if (init.phi == "norm") {
    phi.hat <- (var(lmfit$residuals) - mean(O))/mean(O)^2
  } else if (init.phi == "counts") { 
    phi.hat <- (var(O)-mean(O))/mean(O)^2
  }
  phi.hat <- ifelse(phi.hat < 1e-6,1e-6,phi.hat)
  init.par <- list(goto.cnv=goto.cnv,goto.normal=goto.normal,beta.hat=beta.hat,phi.hat=phi.hat)
  fx.par <- list(S=S,d=d,normal.state=normal.state,fit.var=fit.var)
  nstates <- length(S)
  if (!fit.var) {
    data <- list(O=O,X=X.full)
    nm.fit <- optim(c(logit(goto.cnv),logit(goto.normal),beta.hat,log(phi.hat)),function(par) negLogLike(par,fx.par,data,nstates,stFn,trFn,emFn),method="Nelder-Mead",control=controls)
  } else {
    Y.full <- cbind(intercept=rep(1,nrow(Y)),scale(Y))
    data <- list(O=O,X=X.full,Y=Y.full)
    gamma.hat <- c(phi.hat,rep(0,ncol(Y)))
    nm.fit <- optim(c(logit(goto.cnv),logit(goto.normal),beta.hat,gamma.hat),function(par) negLogLike(par,fx.par,data,nstates,stFn,trFn,emFn),method="Nelder-Mead",control=controls)
  }
  goto.cnv.hat <- logistic(nm.fit$par[1])
  goto.normal.hat <- logistic(nm.fit$par[2])
  A <- trFn(nm.fit$par[1:2],fx.par,data,nstates) 
  beta.hat <- nm.fit$par[3:(2+ncol(X.full))]
  names(beta.hat) <- colnames(X.full)
  mu.hat <- as.numeric(X.full %*% beta.hat)
  mu.hat[mu.hat < 1] <- 1
  if (!fit.var) {
    phi.hat <- exp(nm.fit$par[(3+ncol(X.full))])
    gamma.hat <- NULL
    type="exomeCopy"
  } else {
    gamma.hat <- nm.fit$par[(3+ncol(X.full)):length(nm.fit$par)]
    names(gamma.hat) <- colnames(Y.full)
    phi.hat <- as.numeric(Y.full %*% gamma.hat)
    phi.hat[phi.hat < 1e-6] <- 1e-6
    type="exomeCopyVar"
  }
  path <- viterbiPath(nm.fit$par,fx.par,data,nstates,stFn,trFn,emFn)
  emit.probs <- emFn(nm.fit$par,fx.par,data,nstates)
  log.odds <- log(emit.probs[cbind(path,seq(path))]+1e-6) - log(emit.probs[normal.state,]+1e-6)
  final.par <- list(goto.cnv=goto.cnv.hat,goto.normal=goto.normal.hat,beta=beta.hat,gamma=gamma.hat,phi=ifelse(fit.var,NULL,phi.hat))
  fit <- new("ExomeCopy",sample.name=sample.name,type=type,path=Rle(path),ranges=ranges(rdata),O.norm=as.numeric(O/mu.hat),log.odds=log.odds,fx.par=fx.par,init.par=init.par,final.par=final.par,counts=nm.fit$counts,convergence=nm.fit$convergence,nll=nm.fit$value) 
  return(fit)
}
