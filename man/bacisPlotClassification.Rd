
\name{bacisPlotClassification}
\alias{bacisPlotClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the posterior density of \eqn{\theta} in the classification model.
}
\description{
The classification model is conducted based on the BaCIS method and the posterior density of \eqn{\theta} is plotted.
}
\usage{

bacisPlotClassification(numGroup, tau1, tau2, phi1, phi2,
                        clusterCutoff, MCNum, nDat, xDat, cols, seed)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{numGroup}{Number of subgroups in the trial.}
  \item{tau1}{The precision parameter of subgroups clustering for the classification model.}
  \item{tau2}{The precision prior for the latent variable for the classification.}
  \item{phi1}{Center for the low response rate cluster.}
  \item{phi2}{Center for the high response rate cluster.}
  \item{clusterCutoff}{The cutoff value of the cluster classification. If its value is NA, adaptive classification is applied.}
  \item{MCNum}{The number of MCMC sampling iterations.}
  \item{nDat}{The vector of total sample sizes of all subgroups.}
  \item{xDat}{The vector of the response numbers of all subgroups.}
  \item{cols}{The color vector of all subgroups in the illustration.}
  \item{seed}{Random seed value. If its value is NA, a time dependent random seed is generated and applied.}
}

\value{
The classification model is conducted using the input parameter values and subgroup outcomes. The posterior density of \eqn{\theta} is plotted.

}

\author{
%%  ~~who you are~~
Nan Chen and J. Jack Lee /
Department of Biostatistics
UT MD Anderson Cancer Center
}



\examples{

## Compute the posterior distribution of \eqn{\theta}.
library(bacistool)
bacisPlotClassification(numGroup=5,
                            tau1=NA,
                            tau2=.001,
                            phi1=0.1, phi2=0.3,
                            clusterCutoff=NA,
                            MCNum=5000,
                            nDat=c(25,25,25,25,25),
                            xDat=c(3,4,3,8,7),
                            cols = c("brown", "red", "orange", "blue", "green")

)

}




\keyword{bacistool}% use one of  RShowDoc("KEYWORDS")
\keyword{bacistool}% __ONLY ONE__ keyword per line
