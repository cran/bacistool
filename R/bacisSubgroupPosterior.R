#' @export
#' @importFrom rjags jags.model coda.samples dic.samples
#' @importFrom graphics legend lines par plot points
#' @importFrom stats density rbeta

bacisSubgroupPosterior<- function(numGroup = 5,
                                  tau1 = NA,
                                  tau2 = .001,
                                  phi1 = 0.1,
                                  phi2 = 0.3,
                                  tau4 = 0.1,
                                  alpha = 50,
                                  beta = 2,
                                  AdaptiveCluster = FALSE,
                                  cutOff = 0.92,
                                  MCNum = 50000,
                                  nDat = c(25, 25, 25, 25, 25),
                                  xDat = c(2, 3, 7, 6, 10)
)

{
  if (is.na(tau1))
  {
    sd <- (logit(phi2) - logit(phi1)) / 6
    tau1 <- 1 / sd / sd
  }
  #print(tau1)

  xLim <- max(xDat / nDat) + 0.3
  #print(xDat)
  t <- SubgroupPost(
    xDat = xDat,
    nDat = nDat,
    numGroup = numGroup,
    pp1 = phi1,
    pp2 = phi2,
    alpha = alpha,
    beta = beta,
    tau1 = tau1,
    tau2 = tau2,
    tau4 = tau4,
    AdaptiveCluster = AdaptiveCluster,
    MCNum = MCNum
  )
  return(t)
}
