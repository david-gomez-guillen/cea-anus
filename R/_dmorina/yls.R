yls <- function(scenario1, scenario2, disc.rate = 0)
{
  if (scenario1[1, 1] != scenario2[1, 1]) stop("Different cohort sizes")
  res1 <- le(scenario1, disc.rate)
  res2 <- le(scenario2, disc.rate)
  
  res <- res2 - res1
  cat("Total YLS per strategy\n")
  cat("======================\n")
  cat(res*scenario1[1, 1])
  cat("\n\n")
  cat("Total YLS per person\n")
  cat("====================\n")
  cat(res)
  cat("\n\n")
}