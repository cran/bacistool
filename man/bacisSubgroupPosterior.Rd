
\name{bacisSubgroupPosterior}
\alias{bacisSubgroupPosterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the posterior distribution of response rates of subgroups using the BaCIS method.
}
\description{
In this function, a trial computation is conducted based on the BaCIS model. It calls the OpenBUGS for the Bayesian MCMC sampling for the subgroup classification and hierarchical model information borrowing. The response rate posterior distributions of subgroups are returned from this function.
}
\usage{

bacisSubgroupPosterior(numGroup, tau1, tau2, phi1, phi2, tau4, alpha, beta,
                      clusterCutoff, MCNum, nDat, xDat, seed)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numGroup}{Number of subgroups in the trial.}
  \item{tau1}{The precision parameter of subgroups clustering for the classification model.}
  \item{tau2}{The precision prior for the latent variable for the classification.}
  \item{phi1}{Center for the low response rate cluster.}
  \item{phi2}{Center for the high response rate cluster.}
  \item{tau4}{The precision prior for the center of the cluster in the information borrowing model.}
  \item{alpha}{Hyperprior parameters alpha to control the magnitude of information borrowing model.}
  \item{beta}{Hyperprior parameters beta to control the magnitude of the information borrowing model.}
  \item{clusterCutoff}{The cutoff value of the cluster classification. If its value is NA, adaptive classification is applied.}
  \item{MCNum}{The number of MCMC sampling iterations.}
  \item{nDat}{The vector of total sample sizes of all subgroups.}
  \item{xDat}{The vector of the response numbers of all subgroups.}
  \item{seed}{Random seed value. If its value is NA, a time dependent random seed is generated and applied.}

}

\value{
The MCMC sampling data of the response rate posterior distributions of all subgroups is returned as an matrix format. Each column of the return matrix corresponds to the response rate distribution of one subgroup.

}

\author{
%%  ~~who you are~~
Nan Chen and J. Jack Lee /
Department of Biostatistics
UT MD Anderson Cancer Center
}



\examples{

## Compute the response rate posterior distributioni
## of each subgroup using the BaCIS method
\dontrun{
library(bacistool)
result<-bacisSubgroupPosterior(numGroup=5,
                               tau1=NA,
                               tau2=.001,
                               phi1=0.1, phi2=0.3,
                               tau4=0.1,
                               alpha=50,
                               beta=20,
                               clusterCutoff=NA,
                               MCNum=5000,
                               nDat=c(25,25,25,25,25),
                               xDat=c(3,4,3,8,7))


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.




\keyword{bacistool}% use one of  RShowDoc("KEYWORDS")
\keyword{bacistool}% __ONLY ONE__ keyword per line
