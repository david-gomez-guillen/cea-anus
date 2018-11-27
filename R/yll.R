yll <- function(scenario)
{
  CCdeathsMean1 <- apply(attr(scenario, "CCdeaths"), 2, mean)
  CCyears1 <- (max(scenario$age)+1-as.integer(names(CCdeathsMean1)))%*%CCdeathsMean1
  return(CCyears1)
}