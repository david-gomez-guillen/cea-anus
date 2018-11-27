library(Hmisc)

simCohort <- function(M, period, nsim, seed=NULL, disc=3)
{
  if(!is.null(seed)) 
  {
    set.seed(seed)
  }else{
    set.seed(1234)
  }
  costs1  <- vector()
  qalys1  <- vector()
  le1     <- vector()
  costsd1 <- vector()
  qalysd1 <- vector()
  led1    <- vector()
  costs2  <- vector()
  qalys2  <- vector()
  le2     <- vector()
  costsd2 <- vector()
  qalysd2 <- vector()
  led2    <- vector()

  ### conventional screening vs biomarker-based screening
  for (i in 1:nsim)
  {
    costsd1[i] <- 0
    qalysd1[i] <- 0
    led1[i]    <- 0
    costsd2[i] <- 0
    qalysd2[i] <- 0
    led2[i]    <- 0
    ndna       <- M
    ncyto      <- M
    rcyto      <- rbinom(M, 1, 0.6)                              ### TODO: benign cyto (0) with probability 0.4
    nanos      <- table(rcyto)[2]
    ranos      <- rbinom(nanos, 1, 0.4)                          ### TODO: normal anoscopy (0) with probability 0.6
    nbiops     <- table(ranos)[2]
    rbiops     <- rMultinom(rbind(c(0.3, 0.3, 0.4)), nbiops)     ### TODO: LSIL with prob 0.3, HSIL with prob 0.3, INV CA with prob 0.4
    
    nbiom      <- M
    rcyto2     <- rMultinom(rbind(c(0.7, 0.1, 0.1, 0.1)), nbiom) ### TODO: Benign cyto or ASCUS/LSIL + NEG BIOM with prob 0.7, 
                                                                 ### cyto HSIL + NEG BIOM with prob 0.1, Benign cyto + POS BIOM with prob 0.1,
                                                                 ### abnormal cyto (LSIL/HSIL) + POS BIOM wiht prob 0.1
    nanos2     <- sum(table(rcyto2)[2:4])
    ranos2     <- rbinom(nanos2, 1, 0.6)                         ### TODO: normal anoscopy (0) with probability 0.4
    nbiops2    <- table(ranos2)[2]
    rbiops2    <- rMultinom(rbind(c(0.3, 0.3, 0.4)), nbiops2)    ### TODO: LSIL with prob 0.3, HSIL with prob 0.3, INV CA with prob 0.4
    
    ### UNDISCOUNTED QALYS and COSTS (TODO: All the costs are set to 1€)
    ### conventional screening
    costs1[i] <- ncyto*1 + ndna*1 + table(rcyto)[1]*(period-1)*1*1 + nanos*period*1 + nbiops*period*1 + table(ranos)[1]*(period-1)*1*1 + 
      table(rbiops)[1]*2*(period-1)*1*1 + ### LSIL CONTROL 6Months: Cost of Cytology + HR-HPV-DNA
      table(rbiops)[2]*period*1 +     ### HSIL: Local treatment
      table(rbiops)[3]*period*1       ### Surgery
    qalys1[i] <- table(rcyto)[1]*period*1 + table(ranos)[1]*period*0.97 + table(rbiops)[1]*period*0.87 + 
      table(rbiops)[2]*period*0.76 + table(rbiops)[3]*period*0.67
    le1[i] <- table(rcyto)[1]*period + table(ranos)[1]*period + table(rbiops)[1]*period + 
      table(rbiops)[2]*period + table(rbiops)[3]*period
    ### biomarker-based screening
    costs2[i] <- ncyto*1 + nbiom*1 + table(rcyto2)[1]*(period-1)*1*1 + nanos2*period*1 + nbiops2*period*1 + table(ranos2)[1]*(period-1)*1*1 + 
      table(rbiops2)[1]*2*(period-1)*1*1 + ### LSIL CONTROL 6Months: Cost of Cytology + Biomarker
      table(rbiops2)[2]*period*1 +     ### HSIL: Local treatment
      table(rbiops2)[3]*period*1       ### Surgery
    qalys2[i] <- table(rcyto2)[1]*period*1 + table(ranos2)[1]*period*0.97 + 
      table(rbiops2)[1]*period*0.87 + table(rbiops2)[2]*period*0.76 + table(rbiops2)[3]*period*0.67
    le2[i] <- table(rcyto2)[1]*period + table(ranos2)[1]*period + 
      table(rbiops2)[1]*period + table(rbiops2)[2]*period + table(rbiops2)[3]*period
    ### DISCOUNTED QALYS and COSTS (TODO: All the costs are set to 1€)
    for (j in 1:period)
    {
      costsd1[i] <- costsd1[i] + (costs1[i]/period)/((1+disc/100)^(j-1)) 
      qalysd1[i] <- qalysd1[i] + (qalys1[i]/period)/((1+disc/100)^(j-1))
      led1[i]    <- led1[i] + (le1[i]/period)/((1+disc/100)^(j-1))
      costsd2[i] <- costsd2[i] + (costs2[i]/period)/((1+disc/100)^(j-1)) 
      qalysd2[i] <- qalysd2[i] + (qalys2[i]/period)/((1+disc/100)^(j-1))
      led2[i]    <- led2[i] + (le2[i]/period)/((1+disc/100)^(j-1))
    }
  }
  convent <- data.frame(sim=c(1:nsim), Costs=as.vector(costs1), QALYS=as.vector(qalys1), LE=as.vector(le1), Costs_Disc=as.vector(costsd1),
                        QALYS_Disc=as.vector(qalysd1), LE_Disc=as.vector(led1))
  attr(convent, "ncytos") <- ncyto + (table(rcyto)[1] + table(ranos)[1] + 2*table(rbiops)[1])*(period-1)
  attr(convent, "ndnas")  <- attr(convent, "ncytos")
  attr(convent, "nanos")  <- nanos*period 
  attr(convent, "nHSIL")  <- table(rbiops)[2]
  attr(convent, "nInvCan")<- table(rbiops)[3]
  biomark <- data.frame(sim=c(1:nsim), Costs=as.vector(costs2), QALYS=as.vector(qalys2), LE=as.vector(le2), Costs_Disc=as.vector(costsd2),
                        QALYS_Disc=as.vector(qalysd2), LE_Disc=as.vector(led2))
  attr(biomark, "ncytos") <- ncyto + (table(rcyto2)[1] + table(ranos2)[1] + 2*table(rbiops2)[1])*(period-1)
  attr(biomark, "nbiom")  <- attr(biomark, "ncytos")
  attr(biomark, "nanos")  <- nanos2*period
  attr(biomark, "nHSIL")  <- table(rbiops2)[2]
  attr(biomark, "nInvCan")<- table(rbiops2)[3]
  return(list(convent, biomark))
}

res <- simCohort(100000, 20, 5, 3) ### size, period, nsim, disc (%)
res[[1]] ### Conventional screening
res[[2]] ### Biomarker-based screening

### Number of cytologies
attr(res[[1]], "ncytos") ### Conventional screening
attr(res[[2]], "ncytos") ### Biomarker-based screening

### Number of anoscopies
attr(res[[1]], "nanos") ### Conventional screening
attr(res[[2]], "nanos") ### Biomarker-based screening

### Number of HSIL
attr(res[[1]], "nHSIL")
attr(res[[2]], "nHSIL")

### Number of invasive cancer
attr(res[[1]], "nInvCan")
attr(res[[2]], "nInvCan")

### summary
res_summ1 <- colMeans(res[[1]])
res_summ2 <- colMeans(res[[2]])
res_summ  <- rbind(res_summ1, res_summ2)
res_summ[, 1] <- c(1, 2) # 1=Conventional; 2=Biomarker
res_summ