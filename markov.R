library(CEAModel)
library(ggplot2)
library(reshape2)
library(plotly)
library(dplyr)

EPSILON <- 1e-7  # Maximum allowed error due to precision errors
EXCEL.START.AGE <- 25  # Start age for the first sheet of the excel context file
# Do linear interpolation of contiguous strata when current simulation period falls between two of them
INTERPOLATE.STRATA <- TRUE

DISCOUNT.RATE <- .03

DEFAULT.START.AGE <- 40
DEFAULT.MAX.AGE <- 70

if (INTERPOLATE.STRATA) {
  get.context.stratum <- function(strat.context, year, period) {
    current.stratum <- ((year-EXCEL.START.AGE) %/% 5) + 1
    prev.stratum <- ((year-period-EXCEL.START.AGE) %/% 5) + 1
    if (current.stratum == prev.stratum) {
      ctx <- strat.context[[min(prev.stratum, length(strat.context))]]
    } else {
      current.w <- year %% 5
      prev.w <- period - year %% 5
      current.ctx <- strat.context[[min(current.stratum, length(strat.context))]]
      prev.ctx <- strat.context[[min(prev.stratum, length(strat.context))]]
      ctx <- lapply(names(current.ctx), function(n) {
        return((current.w * current.ctx[[n]] + prev.w * prev.ctx[[n]]) / period)
      })
      names(ctx) <- names(current.ctx)
    }
    return(ctx)
  }
  
  get.matrix.stratum <- function(strat.matrices, year, period) {
    current.stratum <- ((year-EXCEL.START.AGE) %/% 5) + 1
    prev.stratum <- ((year-period-EXCEL.START.AGE) %/% 5) + 1
    if (current.stratum == prev.stratum) {
      mtx <- strat.matrices[[min(prev.stratum, length(strat.matrices))]]
    } else {
      current.w <- year %% 5
      prev.w <- period - year %% 5
      current.mtx <- strat.matrices[[min(current.stratum, length(strat.matrices))]]
      prev.mtx <- strat.matrices[[min(prev.stratum, length(strat.matrices))]]
      mtx <- lapply(names(current.mtx), function(n) {
        return((current.w * current.mtx[[n]] + prev.w * prev.mtx[[n]]) / period)
      })
      names(mtx) <- names(current.mtx)
    }
    return(mtx)
  }
} else {
  get.context.stratum <- function(strat.context, year, period) {
    year <- floor(year-period+EPSILON)  #  TODO: Numerical instability, find more robust alternative
    ctx <- strat.context[[min(((year-EXCEL.START.AGE) %/% 5) + 1, length(strat.context))]]
    return(ctx)
  }
  
  get.matrix.stratum <- function(strat.matrices, year, period) {
    year <- floor(year-period+EPSILON)  #  TODO: Numerical instability, find more robust alternative
    mtx <- strat.matrices[[min(((year-EXCEL.START.AGE) %/% 5) + 1, length(strat.matrices))]]
    return(mtx)
  }
}

# get.context.stratum <- function(strat.context, year) {
#   ctx <- strat.context[[min(((year-25) %/% 5) + 1, length(strat.context))]]
#   return(ctx)
# }
# 
# get.matrix.stratum <- function(strat.matrices, year) {
#   mtx <- strat.matrices[[min(((year-25) %/% 5) + 1, length(strat.matrices))]]
#   return(mtx)
# }

setup.markov <- function(trees, strat.ctx, costs, utilities) {
  tpMatrices <- list()
  extended.strat.ctx <- list()
  for(stratum in names(strat.ctx)) {
    . <- function(x) {
      if (is.na(x) || length(x) == 0)
        return(0)
      else return(x)
    }
    
    context.full <- strat.ctx[[stratum]]
    context <- lapply(context.full, function(e) e[1]) # Base values
    # Assign markov probabilities according to tree outcomes
    
    # Outcomes for non-HSIL
    outcomes <- trees$hiv_msm$summarize(context, prevalence=0)  
    row.names(outcomes) <- outcomes$name
    
    context$p_semestral_followup___no_hsil <- .(sum(outcomes[c('semestral_followup_no_hsil'), 'prob']))
    context$p_semestral_followup_irc___no_hsil <- .(sum(outcomes[c('semestral_followup_no_hsil_irc'), 'prob']))
    context$p_surgery_no_cancer <- .(sum(outcomes[c('surgery_no_cancer'), 'prob']))
    
    # Outcomes for HSIL
    outcomes.undetected <- trees$hiv_msm$summarize(context, prevalence=1) 
    row.names(outcomes.undetected) <- outcomes.undetected$name
    
    context$p_semestral_followup___hsil <- .(sum(outcomes.undetected[c('semestral_followup_hsil'), 'prob']))
    context$p_semestral_followup_irc___hsil <- .(sum(outcomes.undetected[c('semestral_followup_hsil_irc'), 'prob']))
    context$p_cancer___undetected_hsil <- .(sum(outcomes.undetected[c('surgery_cancer'), 'prob']))
    context$p_surgery_no_cancer___undetected_hsil <- .(sum(outcomes.undetected[c('surgery_no_cancer'), 'prob']))
    
    
    # Outcomes for semestral followup (no HSIL)
    outcomes.sem <- trees$semestral_followup$summarize(context, prevalence=0) 
    row.names(outcomes.sem) <- outcomes.sem$name
    
    context$p_hsil_detected___semestral_followup_no_hsil <- .(sum(outcomes.sem[c('semestral_followup_hsil'), 'prob']))
    context$p_hsil_detected___semestral_followup_no_hsil_irc <- .(sum(outcomes.sem[c('semestral_followup_hsil_irc'), 'prob']))
    context$p_surgery_no_cancer___semestral_followup_no_hsil <- .(sum(outcomes.sem[c('surgery_no_cancer'), 'prob']))
    # Assuming IRC does not affect cancer probabilities
    context$p_surgery_no_cancer___semestral_followup_irc_no_hsil <- .(sum(outcomes.sem[c('surgery_no_cancer'), 'prob']))
    
    # Outcomes for semestral followup (HSIL)
    outcomes.sem.hsil <- trees$semestral_followup$summarize(context, prevalence=1) 
    row.names(outcomes.sem.hsil) <- outcomes.sem.hsil$name
    
    context$p_hsil_detected___semestral_followup_hsil <- .(sum(outcomes.sem.hsil[c('semestral_followup_hsil'), 'prob']))
    context$p_hsil_detected___semestral_followup_hsil_irc <- .(sum(outcomes.sem.hsil[c('semestral_followup_hsil_irc'), 'prob']))
    context$p_cancer___semestral_followup_hsil <- .(sum(outcomes.sem.hsil[c('surgery_cancer'), 'prob']))
    
    context$p_cancer___semestral_followup_irc_hsil <- .(sum(outcomes.sem.hsil[c('surgery_cancer'), 'prob'])) / context$hr_cancer_irc
    context$p_surgery_no_cancer___semestral_followup_hsil <- .(sum(outcomes.sem.hsil[c('surgery_no_cancer'), 'prob']))
    context$p_surgery_no_cancer___semestral_followup_irc_hsil <- .(sum(outcomes.sem.hsil[c('surgery_no_cancer'), 'prob']))
    
    # Replace parameters that change with treatment
    if (endsWith(trees$hiv_msm$name, '_t')) {
      context$p_hsil___semestral_followup_no_hsil <- context$p_hsil___semestral_followup_no_hsil_treatment
      context$p_hsil___semestral_followup_no_hsil_irc <- context$p_hsil___semestral_followup_no_hsil_irc_treatment
      context$p_undetected_hsil <- context$p_undetected_hsil_treatment
    }
    
    hiv_msm.cost <- weighted.mean(outcomes$cost, outcomes$prob)
    
    extended.strat.ctx[[stratum]] <- context
    tpMatrices[[stratum]] <- markov$evaluateTpMatrix(context)
  # print(tpMatrices[[stratum]])
  # print(stratum)
    costs[[stratum]]['hiv_positive'] <- hiv_msm.cost
    # browser()
    for(matrix.name in names(tpMatrices)) {
      missing.names <- mapply(function(e){if (is.na(as.numeric(e)) && e != '#') e else NA}, tpMatrices[[matrix.name]][[stratum]])
      missing.names <- missing.names[!is.na(missing.names)]
      missing.names <- missing.names[!duplicated(missing.names)]
      if (length(missing.names) > 0)
        stop(paste0('Variables are missing in matrix: ', paste0(missing.names, collapse=', ')))
    }
  }
  
  return(list(tpMatrices=tpMatrices,
              strat.ctx=extended.strat.ctx,
              costs=costs,
              utilities=utilities))
}

calculate.iteration.measures <- function(additional.info, year, iter, current.state, tpMatrix, cost, eff, ctx) {
  ac_incidence <- current.state %*% tpMatrix$strategy[, 'cancer']
  ac_mortality <- current.state %*% tpMatrix$other[, 'death_cancer']
  
  additional.info <- rbind(additional.info,
                           list(
                             year=year,
                             iter=iter,
                             cost=cost,
                             eff=eff,
                             ac_incidence=ac_incidence,
                             ac_mortality=ac_mortality))
  return(additional.info)
}

simulate.markov <- function(trees, 
                            markov,
                            initial.state, 
                            strat.ctx, 
                            start.age=DEFAULT.START.AGE, 
                            max.age=DEFAULT.MAX.AGE,
                            iters.per.year=2,
                            discount.rate=DISCOUNT.RATE) {
  additional.info <- data.frame(year=start.age-1, 
                                iter=1,
                                cost=0,
                                eff=0,
                                ac_incidence=0,
                                ac_mortality=0)
  
  period <- 1/iters.per.year
  # ASSUMPTION: all node info is equal for all strata
  ctx <- get.context.stratum(strat.ctx, start.age+period, period)
  ctx <- lapply(ctx, function(e) e[1]) # Base values
  evaluated.markov <- markov$evaluate(ctx)
  costs <- sapply(names(strat.ctx), function(stratum) sapply(evaluated.markov$nodes, function(n)n$info$cost), simplify=FALSE, USE.NAMES=TRUE)
  # costs['lynch'] <- lynch.cost
  # costs['postmenopausal_asymptomatic'] <- asymptomatic.cost
  # costs['postmenopausal_bleeding'] <- bleeding.cost
  utilities <- sapply(names(strat.ctx), function(stratum) sapply(evaluated.markov$nodes, function(n)n$info$outcome), simplify=FALSE, USE.NAMES=TRUE)
  
  setup.result <- setup.markov(trees, strat.ctx, costs, utilities)
  tpMatrices <- setup.result$tpMatrices
  strat.ctx <- setup.result$strat.ctx
  costs <- setup.result$costs
  utilities <- setup.result$utilities
  
  current.state <- t(initial.state)
  overall.cost <- 0
  overall.eff <- 0
  states <- data.frame()
  # for(year in seq(start.age, max.age)) {
  #   for(iter in seq(iters.per.year)) {
  
  n.periods <- ceiling((max.age - start.age) * iters.per.year)
  if (((max.age - start.age) * iters.per.year) %% 1 != 0) {
    warning('Simulation time is not multiple of monthly periodicity ',
            12/iters.per.year, 
            ', rounding to highest multiple: ',
            (max.age - start.age), ' -> ', formatC(n.periods/iters.per.year, digits = 2, format='f'))
  }  
  
  for(iter in seq(n.periods)) {
    real.year <- iter / iters.per.year + start.age
    year <- floor(real.year)
    states <- rbind(states, current.state)
    tpMatrix <- get.matrix.stratum(tpMatrices, real.year, period)
    ctx <- get.context.stratum(strat.ctx, real.year, period)
    # tpMatrix <- get.matrix.stratum(tpMatrices, year)
    # ctx <- get.context.stratum(strat.ctx, year)
    
    # year.costs <- get.context.stratum(costs, year)
    # year.utilities <- get.context.stratum(utilities, year)
    iter.costs <- unlist(get.context.stratum(costs, real.year, period))
    iter.utilities <- unlist(get.context.stratum(utilities, real.year, period))
    # Discount used to be ...^(year-start.age), adapted to variable simulation steps
    current.cost <- sum(current.state * iter.costs) * (1-discount.rate)^(iter/iters.per.year)  
    overall.cost <- overall.cost + current.cost
    current.eff <- sum(current.state * iter.utilities) * (1-discount.rate)^(iter/iters.per.year)
    overall.eff <- overall.eff + current.eff
    
    additional.info <- calculate.iteration.measures(additional.info, year, iter, current.state, tpMatrix, current.cost, current.eff, ctx)
    
    next.state <- current.state %*% tpMatrix$strategy %*% tpMatrix$other
    
    if (any(next.state < -EPSILON)) {
      print(trees$hiv_msm$name)
      print(tpMatrix$strategy)
      stop('States with negative populations, probabilities might have errors.')
    }
    else if (any(next.state < 0)) {
      # ASSUMPTION: Small rounding errors, renormalize
      # TODO: Reconsider if other alternative might be better
      next.state[next.state < 0] <- 0
      next.state <- next.state / sum(next.state)
    }
    current.state <- next.state
  }
  states$age <- as.numeric(row.names(states)) + start.age - 1
  melted.states <- reshape2::melt(states, id.vars='age')
  
  strategy.name <- paste0(sapply(trees, function(t) t$name), collapse = '-')
  p <- ggplot(melted.states, aes(x=age, y=value, color=variable)) + 
    geom_line() + 
    ylim(0,1) + 
    xlab('Age') + 
    ylab('Cohort (%)') +
    ggtitle(strategy.name)
  p <- ggplotly(p)
  
  results.df <- data.frame(strategy=strategy.name,
                           C=overall.cost,
                           E=overall.eff / iters.per.year, 
                           stringsAsFactors = FALSE)
  return(list(
    plot=p,
    states=states,
    additional.info=additional.info,
    summary=results.df
  )
  )
}

simulate <- function(type, 
                     strategies, 
                     markov, 
                     strat.ctx, 
                     initial.state, 
                     start.age, 
                     max.age, 
                     discount.rate=DISCOUNT.RATE) {
  results.df <- data.frame()
  markov.outputs <- list()
  additional.info <- list()
  for(tree in strategies) {
    used.trees <- list()
    used.trees[[type]] <- tree
    
    if (endsWith(tree$name, '_t')) {
      used.trees$semestral_followup <- trees[['semestral_followup_t']]
    } else {
      used.trees$semestral_followup <- trees[['semestral_followup']]
    }
    
    periodicity.months <- tree$root$info$periodicity_months
    if (is.null(periodicity.months)) {
      iters.per.year <- 2
      strat.ctx.period <- strat.ctx
    } else {
      iters.per.year <- 12 %/% periodicity.months
      strat.ctx.period <- lapply(strat.ctx, function(ctx) {
        ctx$periodicity_months <- periodicity.months
        ctx
      })
      strat.ctx.period <- refresh.context('periodicity_months', strat.ctx.period, excel.strata.df)
    }
    
    result <- simulate.markov(used.trees, markov, initial.state, strat.ctx,
                              start.age=start.age, max.age=max.age,
                              iters.per.year=iters.per.year,
                              discount.rate=discount.rate)
    markov.outputs[[tree$name]] <- result
    additional.info[[tree$name]] <- result$additional.info
    results.df <- rbind(results.df, result$summary)
  }
  
  r <- CEAModel::analyzeCE(results.df, cost.label='€', eff.label='QALY', plot=TRUE)
  r$info <- markov.outputs
  return(r)
}
