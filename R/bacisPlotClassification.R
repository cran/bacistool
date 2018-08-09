

bacisPlotClassification<- function(numGroup = 5,
                                   tau1 = NA,
                                   tau2 = .001,
                                   phi1 = 0.1,
                                   phi2 = 0.3,
                                   clusterCutoff = NA,
                                   MCNum = 50000,
                                   nDat = c(25, 25, 25, 25, 25),
                                   xDat = c(2, 3, 7, 7, 10),
                                   cols = c("brown", "red", "orange", "blue", "green"),
                                   seed = NA
)

{
  if (is.na(seed))
  {
    set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
  }
  else{
    set.seed(seed)
  }
  if (is.na(tau1))
  {
    sd <- (logit(phi2) - logit(phi1)) / 6
    tau1 <- 1 / sd / sd
    cat("The value of tau1 is set at:", tau1, "\n")
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
    clusterCutoff = clusterCutoff,
    MCNum = MCNum
  )

  allD<-t$allTheta

  x <- c(-100, 100)
  y <- c(0, 0.04)
  par(mfrow = c(1, 1))
  plot(
    x,
    y,
    col = "white",
    main = "Distributions for theta's",
    xlab = expression(theta),
    ylab = "Density"
  )
  legendStr <- c()
  y <- rep(0, numGroup)
  x <- xDat / nDat

  # Get sampling data and draw the plot
  for (i in 1:numGroup) {

    theta.Sampled <- allD[, i]

    lines(density(theta.Sampled),
          col = cols[i],
          lwd = 3)
    legendStr <- c(legendStr, paste("Arm", i))
  }
  legend("topright", legendStr, lwd = 2, col = cols)


}
