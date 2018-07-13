
\name{bacisClassification}
\alias{bacisClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conduct classification for subgroups.
}
\description{
The classification model is conducted based on the BaCIS method and the subgroupos are classified into two clusters: high respone rate cluster and low response rate cluster.
}
\usage{

bacisClassification(numGroup, tau1, tau2, phi1, phi2,
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
The classification model is applied using the input parameter values and subgroup outcomes. The classification results are returned. The return list includes highResponseGroup and lowResponseGroup index vlaues.

}

\author{
%%  ~~who you are~~
Nan Chen and J. Jack Lee /
Department of Biostatistics
UT MD Anderson Cancer Center
}

\examples{

## An example to conduct subgroup classification.

library(bacistool)
result<-bacisClassification(numGroup=5,
                      tau1=NA,
                      tau2=.001,
                      phi1=0.1, phi2=0.3,
                      AdaptiveCluster=TRUE,
                      MCNum=1000,
                      nDat=c(25,25,25,25,25),
                      xDat=c(3,4,3,8,7))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


