#' @export
#' @importFrom rjags jags.model coda.samples dic.samples
#' @importFrom graphics legend lines par plot points
#' @importFrom stats density rbeta

bacisOneTrial <- function(numGroup = 5,
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
                          xDat = c(2, 3, 7, 6, 10),
                          cols = c("brown", "red", "orange", "blue", "green"),
                          clusterCols = c(6, 4),
                          yLim = 22)

{
  if (is.na(tau1))
  {
    sd <- (logit(phi2) - logit(phi1)) / 6
    tau1 <- 1 / sd / sd
  }
  #print(tau1)

  xLim <- max(xDat / nDat) + 0.3
  #print(xDat)
  t <- OneTrial(
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
    MCNum = MCNum,
    xLim = xLim,
    yLim = yLim,
    cols = cols,
    clusterCols = clusterCols
  )
  cutOffSel <- t$cutOffSel
  reject <- t[[2]]
  allProb2 <- t$allProb2
  select <- t[[1]]
  cluster <- select > cutOffSel
  decision <- reject > cutOff
  result <-
    rbind(
      reject,
      allProb2,
      select,
      cluster,
      decision,
      round(t$allResp, 3),
      round(xDat / nDat, 3),
      xDat,
      nDat
    )
  rownames(result) <-
    c(
      "Prob(p_i>phi_1)",
      "Prob(p_i>phi_2)",
      "Prob(theta>0)",
      "Classified to high response cluster",
      "The treatment is effective",
      "Posterior Resp.",
      "Observed Resp.",
      "Number of response",
      "Total sample size"
    )
  colnames(result) <- c(1:numGroup)
  result <- round(result, 3)
  #print(result)
  #return(list(DIC = t$dic, result = result))
  return(result)
}
