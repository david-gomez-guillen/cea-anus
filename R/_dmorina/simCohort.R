simCohort <- function(probs1, stopYear = 85, stepTime = 1, M = 1e5, Nsim = 500, 
                      screenCoverage = NULL, screenSensi = NULL, screenProbs = NULL,
                      screenPrice.md, screenPrice.nmd, screenPrice.i, 
                      costCoefs.md, costCoefs.nmd, ### Direct medical; Direct non medical 
                      costCoefs.i, utilityCoefs, screenPeriod = 1, seed = 1234, 
                      dnaScSensi = NULL, dnaScCost.md = NULL, dnaScCost.nmd = NULL, dnaScCost.i = NULL, papScSensi = NULL, dnaScAgeGroups = NULL, 
                      newScreenPrice.md, newScreenPrice.nmd, newScreenPrice.i, newScreenCoverage, ScreenType = 0) 
  ### ScreenType = 0 -> No screening
  ### ScreenType = 1 -> Conventional screening
  ### ScreenType = 2 -> Biomarker based screening
{
  ll <- as.list(match.call())[-1]     ## 
  myfor <- formals(simCohort)         ## formals with default arguments
  for (v in names(myfor)){
    if (!(v %in% names(ll)))
      ll <- append(ll, myfor[v])      ## if arg is missing I add it
  }
  if(!is.null(seed)) 
  {
    set.seed(seed)
  }else{
    set.seed(1234)
  }
  if (!ScreenType %in% c(0, 1, 2)) stop("Wrong screen type")
  
  browser()
  totNumSteps <-(stopYear-as.numeric(substr(probs1[, 1], 1, 2))[1])/stepTime
  ageGroups <- unique(data.frame(Start=as.numeric(substr(probs1[, 1], 1, 2)), End=as.numeric(substr(probs1[, 1], 4, 5))))
  probs1 <- split(probs1[, 2:dim(probs1)[2]], probs1[, 1])
  
  K <- length(probs1)    ### Age groups
  P <- ncol(probs1[[1]]) ### Health states
  res <- list()
  func <- function(sim){
    probs <- probs1
    cumProbs <- lapply(probs, function(x)
      t(apply(x, 1, cumsum)))
    
    stageNi <- rep(0, P)
    names(stageNi) <- names(probs[[1]])
    stageN <- rep(list(stageNi), totNumSteps)
    attr(stageN[[1]], "newCases")         <- 0 # No new cases in the first step!
    attr(stageN[[1]], "newASCUS")         <- 0
    attr(stageN[[1]], "newHGAIN")         <- 0
    attr(stageN[[1]], "newSurvivalCases") <- 0
    attr(stageN[[1]], "ACdeaths")         <- 0
    attr(stageN[[1]], "Otherdeaths")      <- 0
    attr(stageN[[1]], "utility")          <- M
    attr(stageN[[1]], "probs_k")          <- rep(0, 6)
    attr(stageN[[1]], "probs_u")          <- rep(0, 6)
    attr(stageN[[1]], "nCito")            <- M
    attr(stageN[[1]], "nVPH")             <- M
    stageN[[1]][1]                        <- attr(stageN[[1]], "probs_u")[1] <- M
    
    attr(stageN[[1]], "md_cost")          <- 0
    attr(stageN[[1]], "nmd_cost")         <- 0
    attr(stageN[[1]], "i_cost")           <- 0
    
    currentAge <- ageGroups[1, 1]
    currentAgeGroup <- 1 # Age group
    screenCount <- 1
    for(i in 1:(totNumSteps-1)){ # Step loop
      if(currentAge > ageGroups$End[currentAgeGroup])
        currentAgeGroup <- currentAgeGroup + 1 # Change the age group if pertinent
      J <- which(stageN[[i]] > 0)
      stageN[[i+1]] <- stageN[[i]] # Initialize next step
      attr(stageN[[i+1]], "nCito") <- 0
      attr(stageN[[i+1]], "nVPH")  <- 0
      attr(stageN[[i+1]], "probs_k") <- attr(stageN[[i]], "probs_k")
      attr(stageN[[i+1]], "probs_u") <- attr(stageN[[i]], "probs_u")
      newSurvivalCases <- 0
      cost.md  <- 0
      cost.nmd <- 0
      cost.i   <- 0
      for(j in J){ # Stage loop
        if(j %in% 1:5){ # ASCUS/LGAIN - HGAIN
          nU <- attr(stageN[[i]], "probs_u")[j]
          nK <- attr(stageN[[i]], "probs_k")[j]
          pU <- runif(nU)
          pK <- runif(nK)
          qU <- findInterval(pU, unlist(cumProbs[[currentAgeGroup]][j,])) + 1
          qK <- findInterval(pK, unlist(cumProbs[[currentAgeGroup]][j,])) + 1
          outU <- sum(qU != j)
          outK <- sum(qK != j)
          attr(stageN[[i+1]], "probs_u")[j] <-
            attr(stageN[[i+1]], "probs_u")[j] - outU
          attr(stageN[[i+1]], "probs_k")[j] <-
            attr(stageN[[i+1]], "probs_k")[j] - outK
          if(j == 2)
            attr(stageN[[i+1]], "newASCUS") <-
            sum(qU == 3) + sum(qK == 3)
          if(j == 3)
            attr(stageN[[i+1]], "newHGAIN") <-
            sum(qU == 4) + sum(qK == 4)
          if(j == 4) # Store the number of subjects arriving to AC from HGAIN (new cases)
            attr(stageN[[i+1]], "newCases") <- sum(qU == 5) + sum(qK == 5)
          for(l in union(unique(qU), unique(qK))){
            if(l != j){
              if(l %in% 6:8)
                stageN[[i+1]][l] <- stageN[[i+1]][l] + sum(qU==l) + sum(qK==l)
              else{
                attr(stageN[[i+1]], "probs_u")[l] <- 
                  attr(stageN[[i+1]], "probs_u")[l] + sum(qU==l)
                attr(stageN[[i+1]], "probs_k")[l] <- 
                  attr(stageN[[i+1]], "probs_k")[l] + sum(qK==l)
              }
            }
          }
          
          # if((screening)&&(screenCount%%screenPeriod == 0) && ScreenType == 1){## CONVENTIONAL SCREENING
          #   pU <- runif(attr(stageN[[i+1]], "probs_u")[j])
          #   pK <- runif(attr(stageN[[i+1]], "probs_k")[j])
          #   nScreenedU <- sum(pU < screenCoverage[currentAgeGroup])
          #   nScreenedK <- sum(pK < screenCoverage[currentAgeGroup])
          #   if(dnaScreening && currentAgeGroup%in%dnaScAgeGroups){
          #     cost.md <- cost.md + nScreenedU*dnaScCost.md
          #     cost.nmd <- cost.nmd + nScreenedU*dnaScCost.nmd
          #     cost.i <- cost.i + nScreenedU*dnaScCost.i
          #     p2 <- runif(nScreenedU)
          #     attr(stageN[[i+1]], "nVPH") <- ifelse(!is.null(attr(stageN[[i+1]], "nVPH")), attr(stageN[[i+1]], "nVPH")+nScreenedU+nScreenedK,
          #                                           nScreenedU+nScreenedK)
          #     nDetected <- sum(p2 < dnaScSensi[j])
          #     attr(stageN[[i+1]], "nCito") <- ifelse(!is.null(attr(stageN[[i+1]], "nCito")), attr(stageN[[i+1]], "nCito")+nDetected,
          #                                            nDetected)
          #     p22 <- runif(nDetected)
          #     nPapDetected <- sum(p22 < papScSensi[j])
          #     cost.md <- cost.md + nDetected*screenPrice.md
          #     cost.nmd <- cost.nmd + nDetected*screenPrice.nmd
          #     cost.i <- cost.i + nDetected*screenPrice.i
          #     cost.md <- cost.md + nPapDetected*costCoefs.md[j]
          #     cost.nmd <- cost.nmd + nPapDetected*costCoefs.nmd[j]
          #     cost.i <- cost.i + nPapDetected*costCoefs.i[j]
          #     p3 <- runif(nPapDetected)
          #     nRecovered <- sum(p3 < screenProbs[j])
          #     attr(stageN[[i+1]], "probs_u")[j] <- 
          #       attr(stageN[[i+1]], "probs_u")[j] - nPapDetected
          #     attr(stageN[[i+1]], "probs_k")[j] <- 
          #       attr(stageN[[i+1]], "probs_k")[j] + nPapDetected - nRecovered
          #     if(j <= 5)
          #       attr(stageN[[i+1]], "probs_u")[1] <-
          #       attr(stageN[[i+1]], "probs_u")[1] + nRecovered
          #     else{
          #       stageN[[i+1]][6] <- stageN[[i+1]][6] + nRecovered
          #       newSurvivalCases <- newSurvivalCases + nRecovered
          #     }
          #   }
          #   else{
          #     cost.md <- cost.md + nScreenedU*screenPrice.md
          #     cost.nmd <- cost.nmd + nScreenedU*screenPrice.nmd
          #     cost.i <- cost.i + nScreenedU*screenPrice.i
          #     p2 <- runif(nScreenedU)
          #     attr(stageN[[i+1]], "nVPH") <- 0
          #     attr(stageN[[i+1]], "nCito") <- ifelse(!is.null(attr(stageN[[i+1]], "nCito")), attr(stageN[[i+1]], "nCito")+nScreenedU+nScreenedK,
          #                                            nScreenedU+nScreenedK)
          #     nDetected <- sum(p2 < screenSensi[j])
          #     cost.md <- cost.md + nDetected*costCoefs.md[j]
          #     cost.nmd <- cost.nmd + nDetected*costCoefs.nmd[j]
          #     cost.i <- cost.i + nDetected*costCoefs.i[j]
          #     p3 <- runif(nDetected)
          #     nRecovered <- sum(p3 < screenProbs[j])
          #     attr(stageN[[i+1]], "probs_u")[j] <- 
          #       attr(stageN[[i+1]], "probs_u")[j] - nDetected
          #     attr(stageN[[i+1]], "probs_k")[j] <- 
          #       attr(stageN[[i+1]], "probs_k")[j] + nDetected - nRecovered
          #     if(j <= 5)
          #       attr(stageN[[i+1]], "probs_u")[1] <-
          #       attr(stageN[[i+1]], "probs_u")[1] + nRecovered
          #     else{
          #       stageN[[i+1]][6] <- stageN[[i+1]][6] + nRecovered
          #       newSurvivalCases <- newSurvivalCases + nRecovered
          #     }
          #   }
          #   cost.md <- cost.md + nScreenedK*newScreenPrice.md
          #   cost.nmd <- cost.nmd + nScreenedK*newScreenPrice.nmd
          #   cost.i <- cost.i + nScreenedK*newScreenPrice.i
          # }
          
          ## GENERATE SYMPTOMS (COST COMPUTATIONS)
          if(j == 5){
            cost.md <- cost.md + costCoefs.md[j]*attr(stageN[[i+1]], "probs_k")[j]
            cost.nmd <- cost.nmd + costCoefs.nmd[j]*attr(stageN[[i+1]], "probs_k")[j]
            cost.i <- cost.i + costCoefs.i[j]*attr(stageN[[i+1]], "probs_k")[j]
            p <- runif(attr(stageN[[i+1]], "probs_u")[j])
            nSymp <- sum(p < 1) ## All AC cases present symptoms
            cost.md <- cost.md + nSymp*costCoefs.md[j]
            cost.nmd <- cost.nmd + nSymp*costCoefs.nmd[j]
            cost.i <- cost.i + nSymp*costCoefs.i[j]
            p <- runif(nSymp)
            nCured <- sum(p < screenProbs[j])
            attr(stageN[[i+1]], "probs_u")[j] <-
              attr(stageN[[i+1]], "probs_u")[j] - nSymp
            attr(stageN[[i+1]], "probs_k")[j] <-
              attr(stageN[[i+1]], "probs_k")[j] + nSymp - nCured
            stageN[[i+1]][6] <- stageN[[i+1]][6] + nCured
            newSurvivalCases <- newSurvivalCases + nCured
          }
        }
        else{
          p <- runif(stageN[[i]][j]) # Generate Uniform(0,1) random numbers for each individual at the stage
          q <- findInterval(p, unlist(cumProbs[[currentAgeGroup]][j,])) + 1 # Assign the generated numbers to the cumulative probabilities
          out <- sum(q != j) # Number of individuals leaving the stage
          stageN[[i+1]][j] <- stageN[[i+1]][j] - out # Subtract the number of individuals leaving the stage
          for(l in unique(q)){ # Loop for assigning the leaving individuals to their new stage
            if(l != j){
              if(l %in% 3:5)
                attr(stageN[[i+1]], "probs_u")[l-2] <-
                  attr(stageN[[i+1]], "probs_u")[l-2] + sum(q == l)
              else
                stageN[[i+1]][l] <- stageN[[i+1]][l] + sum(q == l)
            }
          }
          if(j != 6)
            newSurvivalCases <- newSurvivalCases + sum(q == 6)
          if((screenCount%%screenPeriod == 0) && (j == 6) && ScreenType==1){
            p <- runif(stageN[[i+1]][j])
            nScreened <- sum(p < screenCoverage[currentAgeGroup])
            if(currentAgeGroup%in%dnaScAgeGroups){
              cost.md <- cost.md + nScreened*dnaScCost.md
              cost.nmd <- cost.nmd + nScreened*dnaScCost.nmd
              cost.i <- cost.i + nScreened*dnaScCost.i
              cost.md <- cost.md + nScreened*screenPrice.md
              cost.nmd <- cost.nmd + nScreened*screenPrice.nmd
              cost.i <- cost.i + nScreened*screenPrice.i
              attr(stageN[[i+1]], "ACdeaths") <- stageN[[i+1]][7] - stageN[[i]][7]
              names(attr(stageN[[i+1]], "ACdeaths")) <- NULL
              attr(stageN[[i+1]], "Otherdeaths") <- stageN[[i+1]][8] - stageN[[i]][8]
              names(attr(stageN[[i+1]], "Otherdeaths")) <- NULL
            }
          }
        }
      }
      for(j in 1:5)
        stageN[[i+1]][j] <- attr(stageN[[i+1]], "probs_u")[j] + 
        attr(stageN[[i+1]], "probs_k")[j]
      
      ### first HR-DNA-HPV + CYTOLOGY
      nScreenedU <- sum(stageN[[i+1]][1:6])
      cost.md    <- cost.md  + nScreenedU*dnaScCost.md
      cost.nmd   <- cost.nmd + nScreenedU*dnaScCost.nmd
      cost.i     <- cost.i   + nScreenedU*dnaScCost.i
      p2         <- runif(nScreenedU)
      attr(stageN[[i+1]], "nVPH") <- nScreenedU
      nDetected <- sum(p2 < dnaScSensi[j]) ## TODO: La prova determinant sembla ser la CITO, no DNA, confirmar amb Mireia!
      attr(stageN[[i+1]], "nCito") <- nScreenedU
      if(nDetected != 0)
      {
        p22          <- runif(nDetected)
        nPapDetected <- sum(p22 < papScSensi[j])
        cost.md      <- cost.md  + nDetected*screenPrice.md
        cost.nmd     <- cost.nmd + nDetected*screenPrice.nmd
        cost.i       <- cost.i   + nDetected*screenPrice.i
        cost.md      <- cost.md  + nPapDetected*costCoefs.md[j]
        cost.nmd     <- cost.nmd + nPapDetected*costCoefs.nmd[j]
        cost.i       <- cost.i   + nPapDetected*costCoefs.i[j]
        p3           <- runif(nPapDetected)
        nRecovered   <- sum(p3 < screenProbs[j])
        attr(stageN[[i+1]], "probs_u")[j] <- 
          attr(stageN[[i+1]], "probs_u")[j] - nPapDetected
        attr(stageN[[i+1]], "probs_k")[j] <- 
          attr(stageN[[i+1]], "probs_k")[j] + nPapDetected - nRecovered
        if(j <= 4)
          attr(stageN[[i+1]], "probs_u")[1] <-
          attr(stageN[[i+1]], "probs_u")[1] + nRecovered
        else{
          stageN[[i+1]][6] <- stageN[[i+1]][6] + nRecovered
          newSurvivalCases <- newSurvivalCases + nRecovered
        }
        
        cost.md  <- cost.md  + nScreenedK*newScreenPrice.md
        cost.nmd <- cost.nmd + nScreenedK*newScreenPrice.nmd
        cost.i   <- cost.i   + nScreenedK*newScreenPrice.i
      }
      
      attr(stageN[[i+1]], "newSurvivalCases")   <- newSurvivalCases
      attr(stageN[[i+1]], "md_cost")            <- cost.md
      attr(stageN[[i+1]], "nmd_cost")           <- cost.nmd
      attr(stageN[[i+1]], "i_cost")             <- cost.i
      names(stageN)[[i]]                        <- currentAge
      currentAge                                <- currentAge + stepTime
      attr(stageN[[i+1]], "utility")            <- utilityCoefs%*%stageN[[i+1]]
      attr(stageN[[i+1]], "ACdeaths")           <- stageN[[i+1]][7] - stageN[[i]][7]
      names(attr(stageN[[i+1]], "ACdeaths"))    <- NULL
      attr(stageN[[i+1]], "Otherdeaths")        <- stageN[[i+1]][8] - stageN[[i]][8]
      names(attr(stageN[[i+1]], "Otherdeaths")) <- NULL
      screenCount                               <- screenCount + 1
    }
    names(stageN)[[length(stageN)]] <- currentAge
    return(stageN)
  }
  
  res <- lapply(seq(1,Nsim,1), func)
  listDat <- lapply(res, function(x){
    y <- as.data.frame(t(as.matrix(as.data.frame(x))))
    aux <- unlist(lapply(x, function(x) attr(x, "newCases")))
    attr(y, "newCases") <- aux
    aux <- unlist(lapply(x, function(x) attr(x, "newASCUS")))
    attr(y, "newASCUS") <- aux
    aux <- unlist(lapply(x, function(x) attr(x, "newHGAIN")))
    attr(y, "newHGAIN") <- aux
    aux2 <- unlist(lapply(x, function(x) attr(x, "newSurvivalCases")))
    attr(y, "newSurvivalCases") <- aux2
    aux <- unlist(lapply(x, function(x) attr(x, "md_cost")))
    attr(y, "md_cost") <- aux
    aux <- unlist(lapply(x, function(x) attr(x, "nmd_cost")))
    attr(y, "nmd_cost") <- aux
    aux <- unlist(lapply(x, function(x) attr(x, "i_cost")))
    attr(y, "i_cost") <- aux
    aux <- unlist(lapply(x, function(x) attr(x, "utility")))
    attr(y, "utility") <- aux
    aux <- unlist(lapply(x, function(x) attr(x, "ACdeaths")))
    attr(y, "ACdeaths") <- aux
    aux <- unlist(lapply(x, function(x) attr(x, "Otherdeaths")))
    attr(y, "Otherdeaths") <- aux
    aux <- unlist(lapply(x, function(x) attr(x, "nCito")))
    attr(y, "nCito") <- aux
    aux <- unlist(lapply(x, function(x) attr(x, "nVPH")))
    attr(y, "nVPH") <- aux
    rownames(y) <- NULL
    y
  })
  fullDat <- data.frame()
  attr(fullDat, "newCases")         <- data.frame()
  attr(fullDat, "newASCUS")         <- data.frame()
  attr(fullDat, "newHGAIN")         <- data.frame()
  attr(fullDat, "newSurvivalCases") <- data.frame()
  attr(fullDat, "md_cost")          <- data.frame()
  attr(fullDat, "nmd_cost")         <- data.frame()
  attr(fullDat, "i_cost")           <- data.frame()
  attr(fullDat, "utility")          <- data.frame()
  attr(fullDat, "ACdeaths")         <- data.frame()
  attr(fullDat, "Otherdeaths")      <- data.frame()
  attr(fullDat, "nCito")            <- data.frame()
  attr(fullDat, "nVPH")             <- data.frame()
  for(i in 1:length(listDat)){
    fullDat <- rbind(fullDat, cbind(listDat[[i]], sim=i))
    attr(fullDat, "newCases") <- rbind(attr(fullDat, "newCases"),
                                       attr(listDat[[i]], "newCases"))
    attr(fullDat, "newASCUS") <- rbind(attr(fullDat, "newASCUS"),
                                       attr(listDat[[i]],"newASCUS"))
    attr(fullDat, "newHGAIN") <- rbind(attr(fullDat, "newHGAIN"),
                                       attr(listDat[[i]],"newHGAIN"))
    attr(fullDat, "newSurvivalCases") <- rbind(attr(fullDat, "newSurvivalCases"),
                                               attr(listDat[[i]], "newSurvivalCases"))
    attr(fullDat, "md_cost") <- rbind(attr(fullDat, "md_cost"),
                                      attr(listDat[[i]], "md_cost"))
    attr(fullDat, "nmd_cost") <- rbind(attr(fullDat, "nmd_cost"),
                                       attr(listDat[[i]], "nmd_cost"))
    attr(fullDat, "i_cost") <- rbind(attr(fullDat, "i_cost"),
                                     attr(listDat[[i]], "i_cost"))
    attr(fullDat, "utility") <- rbind(attr(fullDat, "utility"),
                                      attr(listDat[[i]], "utility"))
    attr(fullDat, "ACdeaths") <- rbind(attr(fullDat, "ACdeaths"),
                                       attr(listDat[[i]], "ACdeaths"))
    attr(fullDat, "Otherdeaths") <- rbind(attr(fullDat, "Otherdeaths"),
                                          attr(listDat[[i]], "Otherdeaths"))
    attr(fullDat, "nCito") <- rbind(attr(fullDat, "nCito"),
                                    attr(listDat[[i]], "nCito")) 
    attr(fullDat, "nVPH") <- rbind(attr(fullDat, "nVPH"),
                                   attr(listDat[[i]], "nVPH")) 
  }
  
  fullDat$age <- ageGroups[1, 1]:(stopYear-1)
  rownames(attr(fullDat, "newCases")) <- paste0("sim", 1:Nsim)
  newCases1 <- attr(fullDat, "newCases")
  rownames(attr(fullDat, "newASCUS")) <- paste0("sim", 1:Nsim)
  newASCUS <- attr(fullDat, "newASCUS")
  rownames(attr(fullDat, "newHGAIN")) <- paste0("sim", 1:Nsim)
  newHGAIN <- attr(fullDat, "newHGAIN")
  rownames(attr(fullDat, "newSurvivalCases")) <- paste0("sim", 1:Nsim)
  newSurvivalCases1 <- attr(fullDat, "newSurvivalCases")
  rownames(attr(fullDat, "md_cost")) <- paste0("sim", 1:Nsim)
  mdcost1 <- attr(fullDat, "md_cost")
  rownames(attr(fullDat, "nmd_cost")) <- paste0("sim", 1:Nsim)
  nmdcost1 <- attr(fullDat, "nmd_cost")
  rownames(attr(fullDat, "i_cost")) <- paste0("sim", 1:Nsim)
  icost1 <- attr(fullDat, "i_cost")
  rownames(attr(fullDat, "utility")) <- paste0("sim", 1:Nsim)
  utility1 <- attr(fullDat, "utility")
  rownames(attr(fullDat, "ACdeaths")) <- paste0("sim", 1:Nsim)
  ACdeaths1 <- attr(fullDat, "ACdeaths")
  rownames(attr(fullDat, "Otherdeaths")) <- paste0("sim", 1:Nsim)
  Otherdeaths1 <- attr(fullDat, "Otherdeaths")
  rownames(attr(fullDat, "nCito")) <- paste0("sim", 1:Nsim)
  nCito1 <- attr(fullDat, "nCito")
  rownames(attr(fullDat, "nVPH")) <- paste0("sim", 1:Nsim)
  nVPH1 <- attr(fullDat, "nVPH")
  attr(fullDat, "Call") <- ll
  return(fullDat)
}