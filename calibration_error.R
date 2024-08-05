library(dplyr)
library(ggplot2)

setwd('~/Documents/models_ce/anus')

source('load_models.R')
source('markov.R')

set.seed(123)

age.groups <- as.numeric(substr(names(strat.ctx), 2, 3))
CALIB.AGE.GROUPS <- names(strat.ctx)[age.groups >= START.AGE & age.groups < MAX.AGE]
CALIB.PARAMS <- c('p_cancer___hsil_annual', 'p_hsil_regression_annual', 'p_hsil_annual', 'survival_5year')

# PLACEHOLDER
TARGET.INC.HSIL <- c(0.002489277, 0.002469715, 0.002304056, 0.002192828, 0.001754127, 0.001508943, 0.001389184, 0.001325221)

calibration.error <- function(pars) {
  agg.output <- calculate.outputs(pars)
  if (!is.null(agg.output)) {
    sim.inc.hsil <- agg.output$incidence_hsil
    sim.prev.hsil <- agg.output$prevalence_hsil
    
    error <- sum((sim.inc.hsil - TARGET.INC.HSIL)^2)
  } else {
    error <- 1e10
  }
  print(error)
  return(error)
}

calculate.outputs <- function(pars) {
  initial.state <- sapply(markov$nodes,
                          function(n) if (n$name=='hiv_positive') 1 else 0)
  
  calib.strat.ctx <- calib.vec.to.ctx(pars, strat.ctx)
  
  suppressWarnings(sink())
  sink('/dev/null')
  markov.result <- tryCatch(simulate('hiv_msm',
                                     trees['conventional'],
                                     markov,
                                     calib.strat.ctx,
                                     initial.state,
                                     discount.rate=.03),
                            error=function(e) NULL)
  sink()
  
  if (!is.null(markov.result)) {
    output <- markov.result$info[[1]]$additional.info[c('n_healthy', 'n_hsils', 'n_cancers', 'n_new_cancers', 'n_new_deaths_cancer', 'n_new_detected_false_hsils', 'n_new_detected_true_hsils')]
    output <- output[c(-1),]
    output$incidence_hsil <- (output$n_new_detected_true_hsils + output$n_new_detected_false_hsils) / sum(output[1:3])
    output$prevalence_hsil <- output$n_hsils / sum(output[1:3])
    output$incidence_cancer <- output$n_new_cancers / sum(output[1:3])
    output$prevalence_cancer <- output$n_cancers / sum(output[1:3])
    output$mortality_cancer <- output$n_new_deaths_cancer / output$n_cancers
    
    # Aggregate by 5-year age groups (10 semesters)
    output$age.group <- (seq(0,nrow(output)-1) %/% 10)
    agg.output <- output %>% 
      group_by(age.group) %>% 
      summarise_at(c('incidence_hsil', 'prevalence_hsil'), mean)
  } else {
    agg.output <- NULL
  }
  return(agg.output)
}

calib.vec.to.ctx <- function(pars, strat.ctx) {
  calib.strat.ctx <- strat.ctx

  if (length(pars) != length(CALIB.AGE.GROUPS) * length(CALIB.PARAMS)) stop('Parameter vector with wrong length for calibration')
  
  ix <- 1
  for(ag in CALIB.AGE.GROUPS) {
    for(p in CALIB.PARAMS) {
      calib.strat.ctx[[ag]][[p]] <- pars[ix]
      ix <- ix + 1
    }
  }
  
  calib.strat.ctx <- refresh.context(CALIB.PARAMS, calib.strat.ctx, excel.strata.df)
  return(calib.strat.ctx)
}

ctx.to.calib.vec <- function(strat.ctx) {
  ix <- 1
  pars <- numeric(0)
  for(ag in CALIB.AGE.GROUPS) {
    for(p in CALIB.PARAMS) {
      pars <- c(pars, strat.ctx[[ag]][[p]])
      ix <- ix + 1
    }
  }
  
  return(pars)
}

get.initial.guess <- function() {
  return(ctx.to.calib.vec(strat.ctx))
}

### Calibration tests

# res <- optim(ctx.to.calib.vec(strat.ctx), calibration.error)
# 
# initial.output <- calculate.outputs(ctx.to.calib.vec(strat.ctx))$incidence_hsil
# calibrated.output <- calculate.outputs(res$par)$incidence_hsil
# 
# df <- data.frame(x=rep(1:8,3),
#                  error=c(TARGET.INC.HSIL, initial.output, calibrated.output),
#                  type=c(rep('target', 8), rep('initial', 8), rep('calibrated', 8)))
# ggplot(df, aes(x=x, y=error, color=type)) + 
#   geom_line() +
#   scale_color_manual(name='HSIL incidence',
#                      breaks=c('target', 'initial', 'calibrated'),
#                      values=c('black', 'red', 'green')) +
#   theme_minimal()
