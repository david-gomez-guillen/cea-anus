source('decisionTree.R')

N <- 100

trees <- list()
for(f in list.files('./trees')) {
  if (!startsWith(f, '_')) {
    name <- substr(f,1,(nchar(f) - 5))
    print(name)
    t <- parseYAML(paste0('trees/',f))
    assign(name, t)
    trees[name] <- t
  }
}

generateContext <- function(N) {
  context <- list(
    p_cyto_benign=runif(N,.7,.9),
    p_cyto_hpv16la_benign=runif(N, .7, .9),
    p_no_hsil=runif(N,.8,.9),
    p_hsil=runif(N,.1,.2),
    p_regression=runif(N,.4,.6),
    p_cyto_benign_hpv16_neg=runif(N,.7,.8),
    p_cyto_benign_hpv1618_neg=runif(N,.7,.75), # Beware probability consistency!
    p_cyto_benign_hpvhr_neg=runif(N,.75,.8),
    p_cyto_benign_ascus_lsil_arnm_neg=runif(N,.7,.8),
    
    dist_hra = rdirichlet(N, c(60, 3, 1)),
    
    c_cyto = runif(N, 200, 400),
    c_cyto_hpv16la = runif(N, 200, 400),
    c_followup = runif(N, 1000, 2000),
    c_hra = runif(N, 20000, 30000),
    c_surgery = runif(N, 10000, 20000),
    
    u_cyto_benign = runif(N, .9, 1),
    u_no_hsil = runif(N, .8, .9),
    u_hsil = runif(N, .7, .8),
    u_surgery = runif(N, .3, .6)
  )
  return(context)
}

# conventional$compareStrategies(alternatives=list(
#                                 hpv16la=conventional_hpv16la
#                                ),
#                                context=generateContext(N))

# profvis({
# cat(conventional$toString())
# conventional$getNetHealthBenefit(wtp=30000, context)
# conventional$getUsedVariables()
# conventional_hpv16la$getExpectedCost(context)
# conventional$getExpectedEffectiveness(context)
# conventional$getCostEffectiveness(context)
# })

# profvis({
  # context <- generateContext(N)
  # results <- conventional$runIterations(context=context, 'name')
  # N_surgery <- length(results[results %in% c('surgery', 'followup_surgery')])
  # 
  # context <- generateContext(N-N_surgery)
  # context[['p_cyto_benign']] <- runif(N-N_surgery,.8,.95)
  # results <- conventional$runIterations(context=context)
  # N_surgery <- N_surgery + length(results[results %in% c('surgery', 'followup_surgery')])
  # 
  # 
  # context <- generateContext(N-N_surgery)
  # context[['p_cyto_benign']] <- runif(N-N_surgery,.9,1)
  # results <- conventional$runIterations(context=context)
  # N_surgery <- N_surgery + length(results[results %in% c('surgery', 'followup_surgery')])
# })

