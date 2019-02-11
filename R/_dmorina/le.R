le <- function(scenario, disc.rate=0)
{
  utCoefs <- c(1, 1, 1, 1, 1, 1, 0, 1)
  attr(scenario, "Call")$utilityCoefs <- utCoefs
  res2 <- as.data.frame(apply(attr(scenario, "utility"), 2, mean))
  res1 <- do.call("simCohort", attr(scenario, "Call"))
  res <- mean(apply(attr(res1, "utility"), 1, sum))/attr(res1, "Call")$M
  res <- res/(1+disc.rate/100)^(as.numeric(rownames(res2))[1]+0.5)
  return(res)
}