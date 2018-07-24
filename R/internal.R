#' @importFrom rjags jags.model coda.samples dic.samples
#' @importFrom graphics legend lines par plot points
#' @importFrom stats density rbeta


model1Str<-function()
{
  mod1 <- "model
  {
  for (i in 1:numGroups)
  {
  y[i] ~ dbin(p[i],n[i])
  logit(p[i]) <- eta[i]
  eta[i] ~ dnorm(gamma[ind[i]], tau1)
  ind[i] <- step(theta[i]) + 1
  }


  for (i in 1:numGroups)
  {
  theta[i] ~ dnorm(0, tau2)
  }

  #mu1<-logit(phi1)
  #mu2<-logit(phi2)


  gamma[1] <-logit(phi1)
  gamma[2] <-logit(phi2)
}"
  return(mod1)
}

model2Str<-function()
{
  mod2 <- "model
  {
  for (i in 1:numGroups)
  {
  y[i] ~ dbin(p[i],n[i]);
  logit(p[i]) <- eta[i];
  eta[i] ~ dnorm(mu,tau3)
  pg[i] <- step(p[i] - targetResp)
  # Probability that the response rate for each group is greater than targetResp	which is xed response probability *
  }
  #Priors
  mu ~ dnorm(mu0, tau4)
  tau3 ~ dgamma(alpha, beta)
  }
  "
  return(mod2)
}



sigmoid <- function(x)
{
  1.0 / (1 + exp(-x))
}

logit <- function(x)
{
  log(x / (1 - x))
}


compESS <- function(Prec, y, p.obs){
  a <- 1
  b <- 7
  c <- 16 - Prec*(1+y)
  d <- 12 - Prec*(1+y)*(1-y)
  C <- rbind( c(0,0,-d), c(1,0,-c), c(0,1,-b) )
  roots <- eigen(C, symmetric=FALSE, only.values=TRUE)$values
  roots[roots<0] <- 0.0001
  test <- abs( (y/roots)-p.obs )
  o <- order( test, decreasing=FALSE )
  return(roots[o[1]])
}


# Simulate one trial using the BaCIS model
OneTrial<-
  function(xDat,
           nDat,
           numGroup,
           pp1,
           pp2,
           alpha,
           beta,
           tau1,
           tau2,
           tau4,
           clusterCutoff = NA,
           MCNum = 20000,
           xLim = 1.0,
           yLim = 22,
           cols = c("brown", "red", "orange", "blue", "green"),
           clusterCols = c(6, 4))
  {
    targetInd <- 1.5

    mydata <-
      list(
        y = xDat,
        n = nDat,
        numGroups = numGroup,
        phi1 = pp1,
        phi2 = pp2,
        tau1 = tau1,
        tau2 = tau2
      )


    #browser()
    num <- length(xDat)
    inits <- function() {
      list(
        theta = rep(0, num),
        ind = rep(1, num),
        p = rep((pp1 + pp2) / 2, num)
      )
    }
    parameters <- c("theta", "p", "ind")

    jSeed <- floor(runif(1, 1,10000))

    mText1 <- model1Str()
    modelSpec1 <-textConnection(mText1)
    jags <- jags.model(modelSpec1,
                       data = mydata,
                       n.chains = 1,
                       n.adapt = MCNum / 5,quiet = TRUE,
                       inits=list(.RNG.name= "base::Wichmann-Hill",
                                  .RNG.seed= jSeed)
                       )

    chain<-coda.samples(jags,
                         parameters,
                         n.iter=MCNum*5,
                         thin=5
    )
    close(modelSpec1)

    jSeed <- floor(runif(1, 1,10000))
    mText1 <- model1Str()
    modelSpec1 <-textConnection(mText1)
    jags <- jags.model(modelSpec1,
                       data = mydata,
                       n.chains = 5,
                       n.adapt = MCNum / 5,quiet = TRUE,
                       inits=list(.RNG.name= "base::Wichmann-Hill",
                                  .RNG.seed= jSeed)
                       )

    dic  <- dic.samples(jags,
                         n.iter=MCNum
    )
    close(modelSpec1)

    value<-sum(dic[[1]])+sum(dic[[2]])
    cat("DIC: ", value, "\n")


    #fit.result<-fit.r$Stats
    rejectProb <- c()
    #cat("DIC", fit.r$DIC[1,3], "\n")
    #dic <- 9 ##MCSim$DIC
    probTheta <- c()
    allInd <- c()
    allResp <- c()

    # Draw the plot
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
    chain<-chain[[1]]
    # Get sampling data and draw the plot
    for (i in 1:numGroup) {
      #p.sampled<-samplesSample(paste("p[",i,"]",sep=""))
      p.sampled <- chain[, i + numGroup]

      #theta.Sampled<-samplesSample(paste("theta[",i,"]",sep=""))
      theta.Sampled <- chain[, i + 2 * numGroup]

      probTheta <-
        c(probTheta,
          length(theta.Sampled[theta.Sampled > 0]) / length(theta.Sampled))

      #ind.sampled<-samplesSample(paste("ind[",i,"]",sep=""))
      ind.sampled <- chain[ , i]
      indProb <- sum(ind.sampled > targetInd) / length(ind.sampled)
      allInd <- c(allInd, indProb)
      r1 <- p.sampled[ind.sampled == 1]
      lines(density(theta.Sampled),
            col = cols[i],
            lwd = 3)
      legendStr <- c(legendStr, paste("Arm", i))
    }
    legend("topright", legendStr, lwd = 2, col = cols)
    Sys.sleep(3)


    #cat("Adaptive Clustering: ", AdaptiveCluster, "\n")

    # Group Index
    indThreshold <- clusterCutoff
    if (is.na(clusterCutoff))
    {
      meanRate <- mean(xDat / nDat)
      delta <- meanRate - (pp1 + pp2) / 2
      indThreshold <- 1 - 1 / (1 + exp(-2 / (pp2 - pp1) * delta))
      cat("The value of the threshold avlue for classification is: ", indThreshold,"\n")
    }
    #cat("Threshold value for classification: ", indThreshold, "\n")
    highGroup <- which(allInd > indThreshold)
    lowGroup <- which(allInd <= indThreshold)

    sampledP <- matrix(0, MCNum, numGroup)
    priorA <- 1
    priorB <- 1
    pp1S <- pp2S <- c()
    pp1T <- pp2T <- c()
    #High response
    if (length(highGroup) == 1)
    {
      index <- highGroup[1]
      sampledP[, index] <-
        rbeta(MCNum, priorA + xDat[index], priorB + nDat[index] - xDat[index])
      pp2T <- sampledP[, index]
    }

    # The second model
    if (length(highGroup) > 1)
    {
      mydata <-
        list(
          y = xDat[highGroup],
          n = nDat[highGroup],
          numGroups = length(highGroup),
          targetResp = pp1,
          mu0 = logit(pp2),
          tau4 = tau4,
          alpha = alpha,
          beta = beta
        )

      num <- length(highGroup)
      inits <- function() {
        list(
          p = rep(pp2, num),
          mu = logit(pp2),
          tau3 = alpha / beta
        )
      }
      parameters <- c("p", "mu")
      jSeed <- floor(runif(1, 1,10000))
      mText2 <- model2Str()
      modelSpec2 <-textConnection(mText2)
      jags <- jags.model(modelSpec2,
                         data = mydata,
                         n.chains = 1,
                         n.adapt = MCNum / 5, quiet = TRUE,
                         inits=list(.RNG.name= "base::Wichmann-Hill",
                                    .RNG.seed= jSeed)
                         )

      chain<-coda.samples(jags,
                          parameters,
                          n.iter=MCNum*5,
                          thin=5
      )
      close(modelSpec2)
      chain<-chain[[1]]

      #print(chain)
      #fit.result<- BRugsFit(modelFile="model2.txt",data=mydata,numChains=4,
      #                      para=c("p","mu"), nBurnin=MCNum/5,nIter=MCNum,nThin=4,
      #                      DIC=FALSE,BRugsVerbose=FALSE)$Stats

      for (i in 1:length(highGroup)) {
        #p<-samplesSample(paste("p[",i,"]",sep=""))
        p <- chain[, i+1]
        #print(length(p))
        sampledP[, highGroup[i]] <- p
      }
      #pp2S<-samplesSample("mu")
      pp2S <- chain[,1]
      pp2T <- sapply(pp2S, sigmoid)
    }


    #Low response
    if (length(lowGroup) == 1)
    {
      index <- lowGroup[1]
      sampledP[, index] <-
        rbeta(MCNum, priorA + xDat[index], priorB + nDat[index] - xDat[index])
      pp1T <- sampledP[, index]
    }
    if (length(lowGroup) > 1)
    {
      mydata <-
        list(
          y = xDat[lowGroup],
          n = nDat[lowGroup],
          numGroups = length(lowGroup),
          targetResp = pp1,
          mu0 = logit(pp1),
          tau4 = tau4,
          alpha = alpha,
          beta = beta
        )

      num <- length(highGroup)
      inits <- function() {
        list(
          p = rep(pp1, num),
          mu = logit(pp1),
          tau3 = alpha / beta
        )
      }
      parameters <- c("p", "mu")
      jSeed <- floor(runif(1, 1,10000))
      mText2 <- model2Str()
      modelSpec2 <-textConnection(mText2)
      jags <- jags.model(modelSpec2,
                         data = mydata,
                         n.chains = 1,
                         n.adapt = MCNum / 5, quiet = TRUE,
                         inits=list(.RNG.name= "base::Wichmann-Hill",
                                    .RNG.seed= jSeed)
                         )

      chain<-coda.samples(jags,
                          parameters,
                          n.iter=MCNum*5,
                          thin=5
      )

      close(modelSpec2)

      chain<-chain[[1]]
      for (i in 1:length(lowGroup)) {
        #p<-samplesSample(paste("p[",i,"]",sep=""))
        p <- chain[, i+1]
        sampledP[, lowGroup[i]] <- p
      }
      pp1S <- chain[,1]
      pp1T <- sapply(pp1S, sigmoid)
    }

    # Draw the plot
    par(mfrow = c(1, 2))
    y <- c(0, yLim)
    x <- c(0, xLim)
    plot(
      x,
      y,
      col = "white",
      xlab = "Response Rate",
      ylab = "Density",
      main = "Response Rates of All Groups"
    )
    legendStr <- c()
    y <- rep(0, numGroup)
    xObs <- xDat / nDat
    rejectProb <- c()
    allProb2 <- c()
    allESS <- c()
    for (i in 1:numGroup) {
      p.sampled <- sampledP[, i]
      allResp <- c(allResp, mean(p.sampled))
      prob <- sum(p.sampled > pp1) / length(p.sampled)
      prob2 <- sum(p.sampled > pp2) / length(p.sampled)
      size <- compESS(1/var(p.sampled),xDat[i], xObs[i])
      allESS <- c(allESS, size)

      rejectProb <- c(rejectProb, prob)
      allProb2 <- c(allProb2, prob2)
      lines(density(p.sampled), col = cols[i], lwd = 3)
      yy <- -yLim / 80
      if (i > 1)
      {
        for (j in 1:(i - 1))
        {
          if (xObs[i] == xObs[j])
          {
            yy <- yy + 0.8
          }
        }
      }
      #print(allInd)
      points(
        xObs[i],
        yy,
        col = cols[i],
        pch = 17 + 2 * (allInd[i] > indThreshold),
        cex = 2
      )
      # points(xObs[i], -yLim/80+1.8,col=cols[(allInd[i]>indThreshold)+1], pch=19, cex=2)
      legendStr <- c(legendStr, paste("Arm", i))
    }
    legend("topright", legendStr, lwd = 3, col = cols)

    x <- c(0, xLim)
    y <- c(0, yLim)
    plot(
      x,
      y,
      col = "white",
      xlab = "Response Rate",
      ylab = "Density",
      main = "Response Rates of Two Clusters"
    )
    if (length(pp1T) > 0)
      lines(density(pp1T), col = clusterCols[1], lwd = 3)
    if (length(pp2T) > 0)
      lines(density(pp2T), col = clusterCols[2], lwd = 3)

    return (
      list(
        fitR = allInd,
        rejectProb = rejectProb,
        probTheta = probTheta,
        allResp = allResp,
        dic = dic,
        cutOffSel = indThreshold,
        allProb2 = allProb2,
        allESS = allESS
      )
    )
  }




# Subgroup posterior distribution
ModelOne <-
  function(xDat,
           nDat,
           numGroup,
           pp1,
           pp2,
           tau1,
           tau2,
           clusterCutoff,
           MCNum = 20000)
  {

    targetInd <- 1.5

    mydata <-
      list(
        y = xDat,
        n = nDat,
        numGroups = numGroup,
        phi1 = pp1,
        phi2 = pp2,
        tau1 = tau1,
        tau2 = tau2
      )


    #browser()
    num <- length(xDat)
    inits <- function() {
      list(
        theta = rep(0, num),
        ind = rep(1, num),
        p = rep((pp1 + pp2) / 2, num)
      )
    }
    parameters <- c("theta", "p", "ind")

    #jags <- jags.model(paste(tempdir(),'\\model1.txt',sep=""),
    jSeed <- floor(runif(1, 1,10000))
    mText1 <- model1Str()
    modelSpec1 <-textConnection(mText1)
    jags <- jags.model(modelSpec1,
                       data = mydata,
                       n.chains = 1,
                       n.adapt = MCNum / 5,
                       quiet = TRUE,
                       inits=list(.RNG.name= "base::Wichmann-Hill",
                                  .RNG.seed= jSeed)
                       )

    chain<-coda.samples(jags,
                        parameters,
                        n.iter=MCNum*5,
                        thin=5
    )
    close(modelSpec1)



    rejectProb <- c()
    jSeed <- floor(runif(1, 1,10000))
    mText1 <- model1Str()
    modelSpec1 <-textConnection(mText1)
    jags <- jags.model(modelSpec1,
                       data = mydata,
                       n.chains = 5,
                       n.adapt = MCNum / 5,quiet = TRUE,
                       inits=list(.RNG.name= "base::Wichmann-Hill",
                                  .RNG.seed= jSeed)
                       )
    dic  <- dic.samples(jags,
                        n.iter=MCNum
    )
    close(modelSpec1)
    probTheta <- c()
    allInd <- c()
    allResp <- c()

    y <- rep(0, numGroup)
    x <- xDat / nDat
    allTheta<-matrix(0,MCNum,numGroup)
    chain<-chain[[1]]


    # Get sampling data and draw the plot
    for (i in 1:numGroup) {
      #p.sampled<-samplesSample(paste("p[",i,"]",sep=""))
      p.sampled <- chain[, i + numGroup]
      #theta.Sampled<-samplesSample(paste("theta[",i,"]",sep=""))

      theta.Sampled <- chain[, i + 2 * numGroup]
      allTheta[,i]<-theta.Sampled
      probTheta <-
        c(probTheta,
          length(theta.Sampled[theta.Sampled > 0]) / length(theta.Sampled))

      #ind.sampled<-samplesSample(paste("ind[",i,"]",sep=""))
      ind.sampled <- chain[ , i]
      indProb <- sum(ind.sampled > targetInd) / length(ind.sampled)

      allInd <- c(allInd, indProb)
      r1 <- p.sampled[ind.sampled == 1]

    }

    #cat("Adaptive Clustering: ", AdaptiveCluster, "\n")

    # Group Index
    indThreshold <- clusterCutoff
    if (is.na(clusterCutoff))
    {
      meanRate <- mean(xDat / nDat)
      delta <- meanRate - (pp1 + pp2) / 2
      indThreshold <- 1 - 1 / (1 + exp(-2 / (pp2 - pp1) * delta))
      cat("The value of the threshold avlue for classification is: ", indThreshold,"\n")
    }
    #cat("Threshold value for classification: ", indThreshold, "\n")
    highGroup <- which(allInd > indThreshold)
    lowGroup <- which(allInd <= indThreshold)
    return(list(DIC = dic, highResponseGroup=highGroup, lowResponseGroup=lowGroup, allTheta=allTheta))

  }



# Subgroup posterior distribution
SubgroupPost <-
  function(xDat,
           nDat,
           numGroup,
           pp1,
           pp2,
           alpha,
           beta,
           tau1,
           tau2,
           tau4,
           clusterCutoff,
           MCNum = 20000)
  {

    targetInd <- 1.5

    mydata <-
      list(
        y = xDat,
        n = nDat,
        numGroups = numGroup,
        phi1 = pp1,
        phi2 = pp2,
        tau1 = tau1,
        tau2 = tau2
      )


    #browser()
    num <- length(xDat)
    inits <- function() {
      list(
        theta = rep(0, num),
        ind = rep(1, num),
        p = rep((pp1 + pp2) / 2, num)
      )
    }
    parameters <- c("theta", "p", "ind")

    jSeed <- floor(runif(1, 1,10000))
    mText1 <- model1Str()
    modelSpec1 <-textConnection(mText1)
    jags <- jags.model(modelSpec1,
                       data = mydata,
                       n.chains = 1,
                       n.adapt = MCNum / 5,quiet = TRUE,
                       inits=list(.RNG.name= "base::Wichmann-Hill",
                                  .RNG.seed= jSeed)
                       )

    chain<-coda.samples(jags,
                        parameters,
                        n.iter=MCNum*5,
                        thin=5
    )
    close(modelSpec1)



    rejectProb <- c()

    probTheta <- c()
    allInd <- c()
    allResp <- c()

    y <- rep(0, numGroup)
    x <- xDat / nDat

    chain<-chain[[1]]
    # Get sampling data and draw the plot
    for (i in 1:numGroup) {
      #p.sampled<-samplesSample(paste("p[",i,"]",sep=""))
      p.sampled <- chain[, i + numGroup]

      #theta.Sampled<-samplesSample(paste("theta[",i,"]",sep=""))
      theta.Sampled <- chain[, i + 2 * numGroup]

      probTheta <-
        c(probTheta,
          length(theta.Sampled[theta.Sampled > 0]) / length(theta.Sampled))

      #ind.sampled<-samplesSample(paste("ind[",i,"]",sep=""))
      ind.sampled <- chain[ , i]
      indProb <- sum(ind.sampled > targetInd) / length(ind.sampled)
      allInd <- c(allInd, indProb)
      r1 <- p.sampled[ind.sampled == 1]

    }


    #cat("Adaptive Clustering: ", AdaptiveCluster, "\n")

    # Group Index
    indThreshold <- clusterCutoff
    if (is.na(clusterCutoff))
    {
      meanRate <- mean(xDat / nDat)
      delta <- meanRate - (pp1 + pp2) / 2
      indThreshold <- 1 - 1 / (1 + exp(-2 / (pp2 - pp1) * delta))
      cat("The value of the threshold avlue for classification is: ", indThreshold,"\n")
    }
    #cat("Threshold value for classification: ", indThreshold, "\n")
    highGroup <- which(allInd > indThreshold)
    lowGroup <- which(allInd <= indThreshold)

    sampledP <- matrix(0, MCNum, numGroup)
    priorA <- 1
    priorB <- 1
    pp1S <- pp2S <- c()
    pp1T <- pp2T <- c()
    #High response
    if (length(highGroup) == 1)
    {
      index <- highGroup[1]
      sampledP[, index] <-
        rbeta(MCNum, priorA + xDat[index], priorB + nDat[index] - xDat[index])
      pp2T <- sampledP[, index]
    }

    # The second model
    if (length(highGroup) > 1)
    {
      mydata <-
        list(
          y = xDat[highGroup],
          n = nDat[highGroup],
          numGroups = length(highGroup),
          targetResp = pp1,
          mu0 = logit(pp2),
          tau4 = tau4,
          alpha = alpha,
          beta = beta
        )

      num <- length(highGroup)
      inits <- function() {
        list(p = rep(pp2, num),
             mu = logit(pp2),
             tau3 = alpha / beta)
      }
      parameters <- c("p", "mu")
      jSeed <- floor(runif(1, 1,10000))
      mText2 <- model2Str()
      modelSpec2 <-textConnection(mText2)
      jags <- jags.model(modelSpec2,
                         data = mydata,
                         n.chains = 1,
                         n.adapt = MCNum / 5,quiet = TRUE,
                         inits=list(.RNG.name= "base::Wichmann-Hill",
                                    .RNG.seed= jSeed)
                         )

      chain<-coda.samples(jags,
                          parameters,
                          n.iter=MCNum*5,
                          thin=5
      )
      close(modelSpec2)



      #fit.result<- BRugsFit(modelFile="model2.txt",data=mydata,numChains=4,
      #                      para=c("p","mu"), nBurnin=MCNum/5,nIter=MCNum,nThin=4,
      #                      DIC=FALSE,BRugsVerbose=FALSE)$Stats
      chain<-chain[[1]]
      for (i in 1:length(highGroup)) {
        #p<-samplesSample(paste("p[",i,"]",sep=""))
        p <- chain[, i+1]
        #print(length(p))
        sampledP[, highGroup[i]] <- p
      }

    }


    #Low response
    if (length(lowGroup) == 1)
    {
      index <- lowGroup[1]
      sampledP[, index] <-
        rbeta(MCNum, priorA + xDat[index], priorB + nDat[index] - xDat[index])
      pp1T <- sampledP[, index]
    }
    if (length(lowGroup) > 1)
    {
      mydata <-
        list(
          y = xDat[lowGroup],
          n = nDat[lowGroup],
          numGroups = length(lowGroup),
          targetResp = pp1,
          mu0 = logit(pp1),
          tau4 = tau4,
          alpha = alpha,
          beta = beta
        )

      num <- length(highGroup)
      inits <- function() {
        list(p = rep(pp1, num),
             mu = logit(pp1),
             tau3 = alpha / beta)
      }
      parameters <- c("p", "mu")
      jSeed <- floor(runif(1, 1,10000))
      mText2 <- model2Str()
      modelSpec2 <-textConnection(mText2)
      jags <- jags.model(modelSpec2,
                         data = mydata,
                         n.chains = 1,
                         n.adapt = MCNum / 5,quiet = TRUE,
                         inits=list(.RNG.name= "base::Wichmann-Hill",
                                    .RNG.seed= jSeed)
                         )

      chain<-coda.samples(jags,
                          parameters,
                          n.iter=MCNum*5,
                          thin=5
      )

      close(modelSpec2)
      chain<-chain[[1]]
      for (i in 1:length(lowGroup)) {
        #p<-samplesSample(paste("p[",i,"]",sep=""))
        p <- chain[, i+1]
        sampledP[, lowGroup[i]] <- p
      }
    }

    return(sampledP)
  }
