plotHGAINIncidence <- function(..., current=NULL, labels=NULL)
{
  if (dev.cur() != 1) dev.off()
  dots <- list(...)
  if (length(dots) < 1) stop("At least one scenario should be defined")
  
  CIN2Incidence1 <- function(i, x, aux){
    num <- aux[i]
    den <- sum(x[i, 1:4])
    num/den
  }
  newCin2Cases1     <- attr(dots[[1]], "newHGAIN")
  newCin2CasesMean1 <- apply(newCin2Cases1, 2, mean)
  
  tmp <- lapply(as.list(unique(dots[[1]]$age)), function(i, dat){
    aux <- dots[[1]][dots[[1]]$age == i, ]
    apply(aux[,1:(ncol(aux)-2)], 2, mean)
  }, dat = dots[[1]])
  meanDat1 <- as.data.frame(t(as.matrix(as.data.frame(tmp))))
  names(newCin2CasesMean1) <- 1:length(newCin2CasesMean1)
  totCIN2Inc1 <- unlist(lapply(as.list(1:nrow(meanDat1)), CIN2Incidence1,
                           x = meanDat1, aux = newCin2CasesMean1))
  age <- min(dots[[1]]$age):max(dots[[1]]$age)
  par(xpd=TRUE)
  colorets <- "#000000"
  plot(age, totCIN2Inc1*100000, pch = 20, xlab = "Age", ylab = "HGAIN incidence",
       main = "HGAIN Incidence by age", axes = FALSE, type = "l", col=colorets)
  axis(2)
  axis(1)
  if (length(dots)>1)
  {
    colorets <- c(colorets, rainbow(length(dots)-1))
    for (j in 2:length(dots))
    {
      newCin2Cases1     <- attr(dots[[j]], "newHGAIN")
      newCin2CasesMean1 <- apply(newCin2Cases1, 2, mean)
      
      tmp <- lapply(as.list(unique(dots[[j]]$age)), function(i, dat){
        aux <- dots[[j]][dots[[j]]$age == i, ]
        apply(aux[,1:(ncol(aux)-2)], 2, mean)
      }, dat = dots[[j]])
      meanDat1 <- as.data.frame(t(as.matrix(as.data.frame(tmp))))
      names(newCin2CasesMean1) <- 1:length(newCin2CasesMean1)
      totCIN2Inc1 <- unlist(lapply(as.list(1:nrow(meanDat1)), CIN2Incidence1,
                                   x = meanDat1, aux = newCin2CasesMean1))
      age <- min(dots[[j]]$age):max(dots[[j]]$age)
      lines(age, totCIN2Inc1*100000, col=colorets[j])
    }
  }
  if (!is.null(current)) 
  {
    colorets  <- c(colorets, "#00ff00")
    h         <- as.numeric(substr(current[1, 1], 4, 5))[1]-as.numeric(substr(current[1, 1], 1, 2))[1]
    age       <- as.numeric(substr(current[1, 1], 1, 2)):as.numeric(substr(current[dim(current)[1], 1], 4, 5))
    CIN2incid <- rep(current[, 2], each=h+1)
    lines(age, CIN2incid, col=colorets[length(colorets)])
  }
  if (is.null(labels)) labels <- ""
  if (!is.null(current)) labels <- c(labels, "Current HGAIN inc.")
  legend("topright", labels, lty=1, col=colorets)
  return(list(data.frame(age=seq(20:84), val=cbind(as.numeric(apply(newCin2Cases1, 2, mean)))), 
              data.frame(age=seq(20:84), val=cbind(as.numeric(totCIN2Inc1*100000)))))
}