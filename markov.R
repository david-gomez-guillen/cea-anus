library(CEAModel)
library(ggplot2)
library(reshape2)
library(plotly)
library(dplyr)

EPSILON <- 1e-7  # Maximum allowed error due to precision errors
DISCOUNT.RATE <- .03

DEFAULT.START.AGE <- 40
DEFAULT.MAX.AGE <- 70

get.context.stratum <- function(strat.context, year) {
  ctx <- strat.context[[min(((year-25) %/% 5) + 1, length(strat.context))]]
  return(ctx)
}

get.matrix.stratum <- function(strat.matrices, year) {
  mtx <- strat.matrices[[min(((year-25) %/% 5) + 1, length(strat.matrices))]]
  return(mtx)
}

setup.markov <- function(trees, strat.ctx, costs, utilities) {
  tpMatrices <- list()
  extended.strat.ctx <- list()
  for(stratum in names(strat.ctx)) {
    . <- function(x) {
      if (length(x) == 0)
        return(0)
      else return(x)
    }
    
    context.full <- strat.ctx[[stratum]]
    context <- lapply(context.full, function(e) e[1]) # Base values
    # Assign markov probabilities according to tree outcomes
    prevalence <- strat.ctx[[stratum]]$p_cancer
    outcomes <- trees$hiv_msm$summarize(context, prevalence=prevalence)
    row.names(outcomes) <- outcomes$name
    outcomes.undetected <- trees$hiv_msm$summarize(context, prevalence=1) # Outcomes for the previously undetected cancers (so prevalence=100%)
    row.names(outcomes.undetected) <- outcomes.undetected$name
    
    context$p_undetected_cancer <- .(sum(outcomes['followup_cancer', 'prob']))
    context$p_semestral_followup_undetected <- .(sum(outcomes[c('semestral_followup_cancer'), 'prob']))
    context$p_semestral_followup <- .(sum(outcomes[c('semestral_followup_no_cancer'), 'prob']))
    context$p_detected_cancer <- .(sum(outcomes['surgery_cancer', 'prob']))
    # TODO: Handle surgery_no_cancer -> Assuming like followup_no_cancer
    
    context$p_semestral_followup___undetected <- .(sum(outcomes.undetected[c('semestral_followup_cancer'), 'prob']))
    context$p_detected_cancer___undetected <- .(sum(outcomes['surgery_cancer', 'prob']))
    
    hiv_msm.cost <- weighted.mean(outcomes$cost, outcomes$prob)
    
    extended.strat.ctx[[stratum]] <- context
    tpMatrices[[stratum]] <- markov$evaluateTpMatrix(context)
    costs[[stratum]]['hiv_positive'] <- hiv_msm.cost
    
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
  if (iter == 1) {
    incidence <- current.state[, 'hiv_positive'] * tpMatrix$strategy['hiv_positive', 'cancer'] +
      current.state[, 'undetected_cancer'] * tpMatrix$strategy['undetected_cancer', 'cancer'] +
      current.state[, 'semestral_followup1'] * tpMatrix$strategy['semestral_followup1', 'cancer'] +
      current.state[, 'semestral_followup2'] * tpMatrix$strategy['semestral_followup2', 'cancer'] +
      current.state[, 'semestral_followup1_undetected_cancer'] * tpMatrix$strategy['semestral_followup1_undetected_cancer', 'cancer'] +
      current.state[, 'semestral_followup2_undetected_cancer'] * tpMatrix$strategy['semestral_followup2_undetected_cancer', 'cancer']
  } else {
    incidence <- 0
  }
  
  mortality <- current.state[, 'cancer'] * tpMatrix$other['cancer', 'death_cancer'] +
    current.state[, 'semestral_followup1_undetected_cancer'] * tpMatrix$other['semestral_followup1_undetected_cancer', 'death_cancer'] +
    current.state[, 'semestral_followup2_undetected_cancer'] * tpMatrix$other['semestral_followup2_undetected_cancer', 'death_cancer'] +
    current.state[, 'undetected_cancer'] * tpMatrix$other['undetected_cancer', 'death_cancer']
  
  additional.info <- rbind(additional.info,
                           list(
                             year=year,
                             iter=iter,
                             cost=cost,
                             eff=eff,
                             incidence=incidence,
                             mortality=mortality))
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
                                incidence=0,
                                mortality=0)
  
  # ASSUMPTION: all node info is equal for all strata
  ctx <- get.context.stratum(strat.ctx, start.age)
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
  for(year in seq(start.age, max.age)) {
    for(iter in seq(iters.per.year)) {
      states <- rbind(states, current.state)
      tpMatrix <- get.matrix.stratum(tpMatrices, year)
      ctx <- get.context.stratum(strat.ctx, year)
      
      year.costs <- get.context.stratum(costs, year)
      year.utilities <- get.context.stratum(utilities, year)
      
      current.cost <- sum(current.state * year.costs) * (1-discount.rate)^(year-start.age)
      overall.cost <- overall.cost + current.cost
      current.eff <- sum(current.state * year.utilities) * (1-discount.rate)^(year-start.age)
      overall.eff <- overall.eff + current.eff
      
      additional.info <- calculate.iteration.measures(additional.info, year, iter, current.state, tpMatrix, current.cost, current.eff, ctx)
      
      next.state <- current.state
      if (iter == 1) {  
        # Yearly strategies, only performed in the first iteration of the year
        next.state <- next.state %*% tpMatrix$strategy
      }
      next.state <- next.state %*% tpMatrix$other
  
      if (any(next.state < -EPSILON)) {
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
    trees <- list()
    trees[[type]] <- tree
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
    
    result <- simulate.markov(trees, markov, initial.state, strat.ctx,
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
