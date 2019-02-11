cer <- function(data1)
{
  c1  <- apply(attr(data1, "md_cost"), 1, sum) + apply(attr(data1, "nmd_cost"), 1, sum) + apply(attr(data1, "i_cost"), 1, sum)
  u1  <- apply(attr(data1, "utility"), 1, sum)
  cer <- c1/u1
  res <- c(est = mean(cer), max = max(cer), min = min(cer))
  return(res)
}