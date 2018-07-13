
\name{bacisThetaPosterior}
\alias{bacisThetaPosterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the posterior distribution of \eqn{\theta} in the classification model.
}
\description{
The classification model is conducted based on the BaCIS method and the posterior distribution of \eqn{\theta} is returned for further analyses.
}
\usage{

bacisThetaPosterior(numGroup, tau1, tau2, phi1, phi2,
                    AdaptiveCluster, MCNum, nDat, xDat)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{numGroup}{Number of subgroups in the trial.}
  \item{tau1}{The precision parameter of subgroups clustering for the classification model.}
  \item{tau2}{The precision prior for the latent variable for the classification.}
  \item{phi1}{Center for the low response rate cluster.}
  \item{phi2}{Center for the high response rate cluster.}
  \item{AdaptiveCluster}{The parameter to determine whether to use the adaptive clustering method.}
  \item{MCNum}{The number of MCMC sampling iterations.}
  \item{nDat}{The vector of total sample sizes of all subgroups.}
  \item{xDat}{The vector of the response numbers of all subgroups.}

}

\value{
The classification model is conducted using the input parameter values and subgroup outcomes. The posterior distribution of \eqn{\theta} is returned. The returned value is an matrix in which each column corresponds the data of one subgroup.

}

\author{
%%  ~~who you are~~
Nan Chen and J. Jack Lee /
Department of Biostatistics
UT MD Anderson Cancer Center
}

\examples{

## Conduct subgroup classification and
## compute the posterior distribution of \eqn{\theta}.

library(bacistool)
result<-bacisThetaPosterior(numGroup=5,
                      tau1=NA,
                      tau2=.001,
                      phi1=0.1, phi2=0.3,
                      AdaptiveCluster=TRUE,
                      MCNum=500,
                      nDat=c(25,25,25,25,25),
                      xDat=c(3,4,3,8,7)
)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line



