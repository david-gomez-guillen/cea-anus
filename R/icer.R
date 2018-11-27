icer <- function(data1, data2)
{
  c1   <- apply(attr(data1, "md_cost"), 1, sum) + apply(attr(data1, "nmd_cost"), 1, sum) + apply(attr(data1, "i_cost"), 1, sum)
  c2   <- apply(attr(data2, "md_cost"), 1, sum) + apply(attr(data2, "nmd_cost"), 1, sum) + apply(attr(data2, "i_cost"), 1, sum)
  u1   <- apply(attr(data1, "utility"), 1, sum)
  u2   <- apply(attr(data2, "utility"), 1, sum)
  icer <- abs((c1-c2)/(u2-u1))
  res  <- c(est = mean(icer), max = max(icer), min = min(icer))
  
  return(res)
}