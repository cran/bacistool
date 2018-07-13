\name{bacisOneTrial}
\alias{bacisOneTrial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Running one trial computation based on the BaCIS model.
}
\description{
The bacisOneTrial function takes data and parameter values as input. It conducts a trial computation based on the BaCIS model.  It calls the OpenBUGS for the Bayesian MCMC sampling for the subgroup classification and hierarchical model information borrowing. It illustrates plots of the classficaiton results and the posterior response distributions of subgroups, and returns the inference results.
}
\usage{

bacisOneTrial(numGroup, tau1, tau2, phi1, phi2, tau4, alpha, beta,
            AdaptiveCluster, cutOff, MCNum, nDat,xDat, cols, clusterCols, yLim)

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
  \item{AdaptiveCluster}{The parameter to determine whether to use the adaptive clustering method.}
  \item{cutOff}{The posterior cutoff value of the final inference for each subgroup.}
  \item{MCNum}{The number of MCMC sampling iterations.}
  \item{nDat}{The vector of total sample sizes of all subgroups.}
  \item{xDat}{The vector of the response numbers of all subgroups.}
  \item{cols}{The color vector of all subgroups in the illustration.}
  \item{clusterCols}{The color vector of all clusters in the illustration.}
  \item{yLim}{The maximum Y-axis value in the illustration.}




}

\value{
The trial simulation illustrates the plot of posterior distribution of classificaiton, posterior response rates of all subgroups, and the posterior response distributions of two clusters.

It also return a matrix including the following information of all subgroups:


  \item{Prob(p_i>phi_1)}{Posterior probability of response probability being greater than phi_1.}
  \item{Prob(p_i>phi_2)}{Posterior probability of response probability being greater than phi_2.}
  \item{theta>0}{Posterior probability of latent variable being greater than 0.}
  \item{Classified to high response cluster}{0: Classified into the lower response cluster, 1: classified into the high response cluster.}
  \item{The treatment is effective }{0: The subgroup is not effective, 1: the subgroup is effective.}
  \item{Posterior Resp.}{Posterior response rates of subgroups.}
  \item{Observed Resp.}{Observe response rates of subgroups.}
  \item{Number of response}{Number of responses of subgroups.}
  \item{Total sample size}{Total sample sizes of subgroups.}


}

\author{
%%  ~~who you are~~
Nan Chen and J. Jack Lee /
Department of Biostatistics
UT MD Anderson Cancer Center
}

\examples{
\dontrun{
## An example for running a simulation trial using the BaCIS method
library(bacistool)
bacisOneTrial(
  numGroup=5,
  tau1=NA,
  tau2=.001,
  phi1=0.1, phi2=0.3,
  tau4=0.1,
  alpha=50,
  beta=2,
  AdaptiveCluster=FALSE,
  cutOff=0.92,
  MCNum=2000,
  nDat=c(25,25,25,25,25),
  xDat=c(2,3,7,6,10),
  cols=c("brown","red","orange","blue","green"),
  clusterCols=c(6,4),
  yLim=22
)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


