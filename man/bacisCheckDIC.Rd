
\name{bacisCheckDIC}
\alias{bacisCheckDIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the DIC value for the classification model.
}
\description{
In this function, the classification model is applied using the input parameter values and the DIC value is calculated.
}
\usage{

bacisCheckDIC(numGroup, tau1, tau2, phi1, phi2,
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
The classification model is applied using the input parameter values and the DIC value is returned.

}

\author{
%%  ~~who you are~~
Nan Chen and J. Jack Lee /
Department of Biostatistics
UT MD Anderson Cancer Center
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
## An example to compute the DIC value.
library(bacistool)
result<-bacisCheckDIC(numGroup=5,
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



