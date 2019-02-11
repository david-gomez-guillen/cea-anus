qalys <- function(scenario, disc.rate=0)
{
  res <- as.data.frame(apply(attr(scenario, "utility"), 2, mean))
  colnames(res) <- "QALY"
  res$QALY <- res$QALY/(1+disc.rate/100)^(as.numeric(rownames(res))-as.numeric(rownames(res))[1])
  cat("Total QALYs per strategy\n")
  cat("=======================\n")
  cat(mean(sum(res)))
  cat("\n\n")
  cat("Total QALYs per person\n")
  cat("=======================\n")
  cat(mean(sum(res))/scenario[1, 1])
  cat("\n\n")
  return(c(mean(sum(res)), mean(sum(res))/scenario[1, 1]))
}