#### Analysis
library(gdata)
library(WriteXLS)

source("R/byAgeGroup.R")
source("R/byAgeGroup2.R")
source("R/simCohort.R")
source("R/qalys.R")
source("R/le.R")
source("R/costs.R")
source("R/cer.R")
source("R/icer.R")
source("R/plotPrevalence.R")
source("R/plotIncidence.R")
source("R/plotMortality.R")
source("R/plotASCUSIncidence.R")
source("R/plotHGAINIncidence.R")

prev.obs <- read.xls("Data/HNObsPrev.xls") ### HPV prevalence
inc.obs  <- read.xls("Data/HNObsInc.xls")  ### AC incidence
probs    <- read.xls("Data/probs.xls")     ### Transition probabilities matrix

## SCENARIO 1: CONVENTIONAL SCREENING
res1 <- simCohort(probs1=probs, Nsim = 10,
                  screenPrice.md = 0, screenPrice.nmd = 0, screenPrice.i = 0,
                  costCoefs.md = c(0, 0, 248.2, 1461.0, 5417.3, 0, 0, 0),
                  costCoefs.nmd = c(0, 0, 79.5, 189.6, 189.6, 0, 0, 0),
                  costCoefs.i = c(0, 0, 0, 0, 0, 0, 0, 0), 
                  utilityCoefs = c(1, 1, 0.987, 0.76, 0.67, 0.938, 0, 0), 
                  screenProbs = c(0,0,1,0.7064,0.3986,0,0,0), ScreenType = 0)

qalys1_1 <- qalys(res1); qalys1_2 <- qalys(res1, disc.rate = 3)
costs1_1 <- costs(res1); costs1_2 <- costs(res1, disc.rate = 3)
le1_1 <- le(res1); le1_2 <- le(res1, disc.rate = 3)
res1_1 <- byAgeGroup2(aggregate(cbind(apply(attr(res1, "nCito"), 2, mean)), by=list(rownames(cbind(apply(attr(res1, "nCito"), 2, mean)))), mean))[, 2]
res1_2 <- byAgeGroup2(aggregate(cbind(apply(attr(res1, "nVPH"), 2, mean)), by=list(rownames(cbind(apply(attr(res1, "nVPH"), 2, mean)))), mean))[, 2]
res1_3 <- byAgeGroup(plotPrevalence(res1))[, 2]
res1_4 <- byAgeGroup2(plotASCUSIncidence(res1)[[1]])[, 2]
res1_5 <- byAgeGroup(plotASCUSIncidence(res1)[[2]])[, 2]
res1_6 <- byAgeGroup2(plotHGAINIncidence(res1)[[1]])[, 2]
res1_7 <- byAgeGroup(plotHGAINIncidence(res1)[[2]])[, 2]
res1_8 <- byAgeGroup2(plotIncidence(res1)[[1]])[, 2]
res1_9 <- byAgeGroup(plotIncidence(res1)[[2]])[, 2]
res1_10 <- byAgeGroup2(plotMortality(res1)[[1]])[, 2]
res1_11 <- byAgeGroup(plotMortality(res1)[[2]])[, 2]
result <- c(qalys1_1[1], qalys1_2[1], qalys1_1[2], qalys1_2[2], le1_1, le1_2, NA, NA,costs1_1[1], costs1_2[1], costs1_1[2], costs1_2[2], costs1_1[3], costs1_2[3],
            costs1_1[4], costs1_2[4], costs1_1[5], costs1_2[5], costs1_1[6], costs1_2[6], 
            res1_1, res1_2, res1_3, res1_4, res1_5, res1_6, res1_7, res1_8, res1_9, res1_10,
            res1_11, sum(res1_1), sum(res1_2), mean(res1_3), sum(res1_4), mean(res1_5),
            sum(res1_6), mean(res1_7), sum(res1_8), mean(res1_9), sum(res1_10), mean(res1_11))
WriteXLS(as.data.frame(result), paste0(getwd(), '/res1.xls'))
