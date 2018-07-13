#' @export
#' @importFrom rjags jags.model coda.samples dic.samples
#' @importFrom graphics legend lines par plot points
#' @importFrom stats density rbeta

# Posterior distribution of Theta
bacisThetaPosterior<- function(numGroup = 5,
                               tau1 = NA,
                               tau2 = .001,
                               phi1 = 0.1,
                               phi2 = 0.3,
                               AdaptiveCluster = FALSE,
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
  t <- ModelOne(
    xDat = xDat,
    nDat = nDat,
    numGroup = numGroup,
    pp1 = phi1,
    pp2 = phi2,
    tau1 = tau1,
    tau2 = tau2,
    AdaptiveCluster = AdaptiveCluster,
    MCNum = MCNum
  )

  return(t$allTheta)
}
