options(java.parameters = "-Xss256m")

#setwd(model.wd)
library(here)

source('markov.R')
source('load_models.R')

get.strategies <- function() {
  #return(names(strategies$hiv_msm))
  strats <- names(strategies$hiv_msm)
  return(list(
    ARN=strats[startsWith(strats, 'arn')],
    `ASCUS/LSIL diff`=strats[startsWith(strats, 'ascus')],
    Conventional=strats[startsWith(strats, 'conventional')],
    `No intervention`='no_intervention'
    )
  )
}

get.strata <- function() {
  return(names(strat.ctx))
}

get.parameters <- function() {
  indep.pars <- getIndependentParameters(paste0(here(), '/params/context.xlsx'))
  indep.pars <- indep.pars[!indep.pars %in% c('periodicity_times_in_year')]
  is.age.specific <- sapply(indep.pars, function(p) {
    par.values <- sapply(strat.ctx, function(ctx) ctx[[p]])
    return(!all(par.values==par.values[1]))
  })

  par.list <- list()
  for(p in indep.pars) {
    if (!is.age.specific[p]) {
      par.list <- append(par.list, list(
        list(
          name=p,
          base.value=strat.ctx[[1]][[p]]
        )
      ))
    } else {
      for(stratum in get.strata()) {
        par.list <- append(par.list, list(
          list(
            name=p,
            base.value=strat.ctx[[stratum]][[p]],
            stratum=stratum
          )
        ))
      }
    }
  }
  par.list <- append(par.list, list(
    list(
      name='discount',
      base.value=0.03
    )
    # list(
    #   name='starting.age',
    #   base.value=40
    # ),
    # list(
    #   name='max.age',
    #   base.value=80
    # )
  ))

  par.list <- lapply(par.list, function(li) {
    if (any(startsWith(li$name, c('p_', '.p_', 'hr_', '.rate_', '.sensitivity_', '.specificity_')))) {
      li$min.value <- 0
      li$max.value <- 1
      li$class <- 'Probabilities'
    } else if (any(startsWith(li$name, c('c_', '.c_')))) {
      li$min.value <- 0
      if (is.null(li$max.value)) li$max.value <- li$base.value * 4
      li$class <- 'Costs'
    } else if (startsWith(li$name, 'u_')) {
      li$min.value <- 0
      li$max.value <- 1
      li$class <- 'Utilities'
    } else if (li$name == 'discount') {
      li$min.value <- 0
      li$max.value <- 0.05
      li$class <- 'Other Parameters'
    } else {
      li$class <- 'Other Parameters'
    }
    return(li)
  })
  class.order <- c('Costs', 'Probabilities', 'Utilities', 'Other Parameters')
  par.list <- par.list[order(sapply(par.list, function(p) match(p$class, class.order)), sapply(par.list, function(p) p$name))]

  return(par.list)
}

run.simulation <- function(strategies, pars, start.age=40, max.age=80) {
  # discount <- pars[[which(sapply(pars, function(p) p$name=='discount'))]]$base.value
  # starting.age <- pars[[which(sapply(pars, function(p) p$name=='starting.age'))]]$base.value
  # max.age <- pars[[which(sapply(pars, function(p) p$name=='max.age'))]]$base.value

  for(pn in names(pars)) {
    # if (pn %in% names(strat.ctx[[1]])) {
      par <- pars[[pn]]
      if (is.list(par)) {
        # strat.ctx <- lapply(par, function(stratum) {
        #   strat.ctx[[stratum]][[pn]] <- par[[stratum]]
        # })
        for (stratum in names(par)) {
          strat.ctx[[stratum]][[pn]] <- par[[stratum]]
        }
        # strat.ctx[[par$stratum]][[par$name]] <- par$base.value
      } else {
        strat.ctx <- lapply(strat.ctx, function(ctx) {
          ctx[[pn]] <- par
          return(ctx)
        })
      }
      # par <- pars[[which(sapply(pars, function(x) x$name==pn))]]


      # if (is.null(par$stratum)) {
      #   strat.ctx <- lapply(strat.ctx, function(ctx) {
      #     ctx[[pn]] <- par$base.value
      #     return(ctx)
      #   })
      # } else {
      #   strat.ctx[[par$stratum]][[par$name]] <- par$base.value
      # }
  }
  strat.ctx <- refresh.context(independent.pars, strat.ctx, excel.strata.df, context.setup)

  sim.strategies <- trees[names(trees) %in% strategies]
  initial.state <- sapply(markov$nodes,
                          function(n) if (n$name=='hiv_positive') 1 else 0)
  results <- simulate('hiv_msm',
                      sim.strategies,
                      markov,
                      strat.ctx,
                      initial.state,
                      start.age=start.age, #pars$starting.age,
                      max.age=max.age, #pars$max.age,
                      discount.rate=pars$discount)

  results$summary$hsil_incidence <- sapply(results$info, function(strat) {
    return(mean(strat$additional.info$incidence_hsil) * 2)
  })
  results$summary$n_cyto <- sapply(results$info, function(strat) {
    return(mean(strat$additional.info$n_cyto))
  })
  results$summary$n_hpv <- sapply(results$info, function(strat) {
    return(mean(strat$additional.info$n_hpv))
  })
  
  results$summary$strategy <- str_extract(results$summary$strategy, '(.*)?-.*', group=1)
  results$summary$domination <- NULL
  return(results)
}

get.calibration.schemes <- function() {
  return(list(
    standard=list(
      description='Standard calibration',
      parameters='p_hsil_regression_annual',
      target=list(
        `Cancer incidence`=c(0.0006868571, 0.0009250000, 0.0009764286, 0.0010300000, 0.0010844000, 0.0011364000, 0.0011884000, 0.0012404000)
      ),
      strata=get.strata()[2:11],
      error_function=calibration.error,
      latent_space_training_set=generate.training.dataset,
      latent_space_training_epochs=50,
      latent_space_latent_dim=7,
      plots=function(){}
    )
  ))
}

calibration.error <- function(pars, target) {
  calibration.strategy <- 'no_intervention'
  target.inc <- target$`Cancer incidence` 
  results <- run.simulation(calibration.strategy, pars, start.age=35)
  ac.inc <- results$info[[calibration.strategy]]$additional.info$incidence_cancer
  ac.inc <- ac.inc[2:length(ac.inc)]
  ac.inc <- colMeans(matrix(ac.inc, nrow=10))
  ac.inc <- ac.inc[2:length(ac.inc)]
  names(ac.inc) <- get.strata()[4:11]
  error <- sum((ac.inc-target.inc)^2)
  result <- list(
    error=error,
    output=list(ac.incidence=ac.inc)
  )
  return(result)
}

generate.training.dataset <- function(initial_guess, n=1000, ...) {
  f.pars <- list(...)
  variation <- f.pars$variation

  n_params <- length(initial_guess)

  dataset <- matrix(NA, nrow=n, ncol=n_params)

  for(i in 1:n) {
	  factors <- runif(n_params, min=1-.25, max=1+.25)
	  dataset[i,] <- pmin(1, initial_guess * factors)
  }

  dataset <- dataset[sample(nrow(dataset)),]

  return(dataset)
}

# params <- lapply(get.parameters(), function(p) p$base.value)
# names(params) <- sapply(get.parameters(), function(p) p$name)
# run.simulation(c('conventional'), params)
