source('decisionTree.R')

trees <- list()
for(f in list.files('./trees')) {
  if (!startsWith(f, '_') && endsWith(f, '.yaml')) {
    name <- substr(f,1,(nchar(f) - 5))
    print(name)
    t <- parseYAML(paste0('trees/',f))
    assign(name, t)
    trees[name] <- t
  }
}

generateContext <- function(N, year=1) {
  context <- list(
    p_cyto_benign=switch(year,
                         '1'=runif(N,.7,.9),
                         '2'=runif(N,.8,.95),
                         '3'=runif(N,.9,1)),
    p_cyto_hpv16la_benign=runif(N, .7, .9),
    p_cyto_benign_hpv16_neg=runif(N,.7,.8),
    p_cyto_benign_hpv16la_neg=runif(N,.7,.8),
    p_cyto_benign_hpv1618_neg=runif(N,.7,.75), # Beware probability consistency!
    p_cyto_benign_hpv1618la_neg=runif(N,.7,.75), # Beware probability consistency!
    p_cyto_benign_hpvhr_neg=runif(N,.75,.8),
    p_cyto_benign_hpvhrla_neg=runif(N,.75,.8),
    p_cyto_benign_hpvhrhc_neg=runif(N,.75,.8),
    p_cyto_benign_ascus_lsil_arnm_neg=runif(N,.7,.8),
    p_no_hsil=runif(N,.8,.9),
    p_hsil=runif(N,.1,.2),
    p_regression=runif(N,.4,.6),
    
    dist_hra = rdirichlet(N, c(60, 3, 1)),
    
    c_cyto = runif(N, 200, 400),
    c_cyto_hpvla = runif(N, 200, 400),
    c_cyto_hpvhc = runif(N, 100, 300),
    c_cyto_arnm = runif(N, 300, 500),
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

N <- 1000

# CEAsummary<- conventional$compareStrategies(alternatives=list(
#                                 hpv16=conventional_hpv16la,
#                                 hpv1618=conventional_hpv1618la,
#                                 hpvhrla=conventional_hpvhrla,
#                                 hpvhrhc=conventional_hpvhrhc,
#                                 arnm=arnme6e7
#                                ),
#                                context=generateContext(N))

print(CEAsummary)

# profvis(
  for (tree_name in names(trees)) {
      t <- trees[[tree_name]]
      context <- generateContext(N, year=1)
      results <- conventional$runIterations(context=context, 'name')
      N_hra <- length(results[results %in% c('hra_annual_followup', 'regression_annual_followup', 'followup_surgery', 'surgery')])

      context <- generateContext(N-N_hra, year=2)
      results <- conventional$runIterations(context=context)
      N_hra <- N_hra + length(results[results %in% c('hra_annual_followup', 'regression_annual_followup', 'followup_surgery', 'surgery')])

      context <- generateContext(N-N_hra, year=3)
      results <- conventional$runIterations(context=context)
      N_hra <- N_hra + length(results[results %in% c('hra_annual_followup', 'regression_annual_followup', 'followup_surgery', 'surgery')])
      cat(paste0(tree_name, ': ', N_hra, ' out of ', N, '\n'))
  }
# )

