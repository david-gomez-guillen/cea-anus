library(dplyr)
library(ggplot2)
library(cmaes)
library(gridExtra)

sink()

# setwd('~/Documents/models_ce/anus')

source('load_models.R')
source('markov.R')

set.seed(123)

age.groups <- as.numeric(substr(names(strat.ctx), 2, 3))
CALIB.START.AGE <- DEFAULT.START.AGE$hiv_msm - 5
CALIB.MAX.AGE <- DEFAULT.MAX.AGE$hiv_msm
CALIB.AGE.GROUPS <- names(strat.ctx)[age.groups >= CALIB.START.AGE & age.groups < CALIB.MAX.AGE]
K <- length(CALIB.AGE.GROUPS)
CALIB.PARAMS <- c(
  # 'p_cancer___hsil_annual_treatment'
  # ,'p_cancer___hsil_annual_no_treatment'
  'p_cancer___hsil_annual'
  #,'p_hsil_annual___cd4_g500'
  #,'p_hsil_regression_annual'
  #,'p_hsil_annual'
  # ,'survival_5year'
)
# CALIB.PARAMS <- c('p_cancer___hsil_annual')


# First calibration HSIL incidence multiplied by 4.5 so the mean incidence equal to Deshmuk (0.102827)
# TARGET.INC.HSIL <- c(0.03266429, 0.06450535, 0.07406059, 0.09361854, 0.12650886, 0.14022268, 0.14507399, 0.14596170)
TARGET.INC.HSIL <- 0.054200565 #0.102827


### Forcing hsil incidence
# for(i in seq_along(TARGET.INC.HSIL)) {
#   strat.ctx[CALIB.AGE.GROUPS][[i]]$p_undetected_hsil <- TARGET.INC.HSIL[i]
# }

 # Source: Clifford 2021 - "A meta‐analysis of anal cancer incidence by risk group: Toward a unified anal cancer risk scale"
 # Ages: 30-44, 45-59, >=60
TARGET.INC.AC <- c(66.2, 99.7, 107.5) / 100000
# target.df <- data.frame(x=c(37, 52, 70), y=c(66.2, 99.7, 107.5)/100000)
target.df <- data.frame(x=c(30, 37, 44, 45, 52, 59, 79),
                        y=c(60.7,66.2,72,92.5,99.7,107.4,128.2)/100000)

lin.interp <- function(x, target.df) {
  ret.values <- sapply(x, function(xi) {
    if (xi < target.df$x[1]) {
      ind <- 1
    } else if (xi >= target.df$x[nrow(target.df)]) {
      ind <- nrow(target.df) - 1
    } else {
      ind <- which(target.df$x > xi)[1] - 1
    }

    slope <- (target.df$y[ind+1]-target.df$y[ind]) / (target.df$x[ind+1] - target.df$x[ind])
    intercept <- target.df$y[ind] - slope * target.df$x[ind]
    return(intercept + slope*xi)
  })
  return(ret.values)
}
# x <- (DEFAULT.START.AGE$hiv_msm):(DEFAULT.MAX.AGE$hiv_msm)
# df <- data.frame(x=x, y=lin.interp(x, target.df), group=c(sapply(1:length(x), function(i) rep(i, 5))))
# df <- df %>% group_by(group) %>% summarise_at(vars(y), list(inc.avg=mean))
# TARGET.INC.AC <- df$inc.avg
# plot(TARGET.INC.AC)
TARGET.INC.AC <- lin.interp(seq((CALIB.START.AGE),(CALIB.MAX.AGE-1),5), target.df)
# plot(x=c(seq((DEFAULT.START.AGE$hiv_msm),(DEFAULT.MAX.AGE$hiv_msm-1),5)), y=c(TARGET.INC.AC), xlim=c(30,80), ylim=c(0.0005, 0.0013))
# plot(x=c(seq((DEFAULT.START.AGE$hiv_msm),(DEFAULT.MAX.AGE$hiv_msm-1),5), target.df$x), y=c(TARGET.INC.AC, target.df$y), ylim=c(0.0005, 0.0013))
# TARGET.INC.AC <- c(0.013300098, 0.016618926, 0.019937754, 0.023256581, 0.03231475, 0.04137292, 0.050431089, 0.059489258)

eval.count <- 0

calibration.error <- function(pars) {
  result <- run.calib.simulation(pars)
  eval.count <<- eval.count + 1
  # agg.output <- calculate.outputs(result)
  
  if (!is.null(result)) {
    # sim.inc.hsil <- agg.output$incidence_hsil
    # sim.prev.hsil <- agg.output$prevalence_hsil
    #
    # error <- sum((sim.inc.hsil - TARGET.INC.HSIL)^2)
    INCREASING.CONSTRAINT.PENALTY <- 1e7

    ac.inc <- calculate.ac.incidence(result)
    hsil.inc <- calculate.hsil.incidence(result)
    error <- sum((((ac.inc - TARGET.INC.AC)/TARGET.INC.AC)^2)) +
      # sum((((hsil.inc - TARGET.INC.HSIL)/TARGET.INC.HSIL)^2))
      #((mean(hsil.inc) - TARGET.INC.HSIL)/TARGET.INC.HSIL)^2 +
      # max(-sum(diff(ac.inc)) * INCREASING.CONSTRAINT.PENALTY, 0)
      # sum((pmax(diff(ac.inc), 0))^2) * INCREASING.CONSTRAINT.PENALTY
      0
  } else {
    error <- 1e10
  }
  print(error)
  if (error < best.solution.val) {
    best.solution <<- pars
    best.solution.val <<- error
    print('BEST: ')
    print(best.solution)

    plot.calibration.state(best.solution)
  }
  return(error)
}

calculate.ac.incidence <- function(par) {
  if(is.list(par)) {
    markov.result <- par
  } else {
    markov.result <- run.calib.simulation(par)
  }
  
  output <- markov.result$info[[1]]$additional.info[c('year', 'incidence_cancer')]

  output$age.group <- cut(output$year, c(seq(CALIB.START.AGE,CALIB.MAX.AGE,5), 200), right=FALSE, labels=FALSE)

  output <- output[3:nrow(output),]
  
  agg.output <- output %>% group_by(age.group) %>% summarise_at(c('incidence_cancer'), mean) * 2 # Data is semestral, annual mean is 2 * semestral mean
  return(agg.output$incidence_cancer)
}


calculate.hsil.incidence <- function(par) {
  if(is.list(par)) {
    markov.result <- par
  } else {
    markov.result <- run.calib.simulation(par)
  }
  output <- markov.result$info[[1]]$additional.info[c('year', 'incidence_hsil')]

  output$age.group <- cut(output$year, c(seq(CALIB.START.AGE,CALIB.MAX.AGE,5), 200), right=FALSE, labels=FALSE)

  output <- output[3:nrow(output),]
  
  agg.output <- output %>%
    group_by(age.group) %>%
    summarise_at(c('incidence_hsil'), mean) * 2 # Data is semestral, annual mean is 2 * semestral mean
  return(agg.output$incidence_hsil)
}

calculate.hsil.prevalence <- function(par) {
  if(is.list(par)) {
    markov.result <- par
  } else {
    markov.result <- run.calib.simulation(par)
  }
  output <- markov.result$info[[1]]$additional.info[c('year', 'n_hsils')]
  # output$year.period <- c(1, rep(1:40, each=2))
  agg.output <- output %>%
    group_by(year) %>%
    summarise_at(c('n_hsils'), sum)

  agg.output$age.group <- cut(agg.output$year, c(seq(40,79,5), 200), right=FALSE, labels=FALSE)
  agg.output <- agg.output %>%
    group_by(age.group) %>%
    summarise_at(c('n_hsils'), mean)
  return(agg.output$n_hsils)
}

run.calib.simulation <- function(pars) {
  initial.state <- sapply(markov$nodes,
                          function(n) if (n$name=='hiv_positive') 1 else 0)
  
  calib.strat.ctx <- calib.vec.to.ctx(pars, strat.ctx)

  suppressWarnings(sink())
  sink('/dev/null')
  markov.result <- tryCatch(simulate('hiv_msm',
                                     trees['no_intervention'],
                                     markov,
                                     calib.strat.ctx,
                                     initial.state,
                                     discount.rate=.03,
                                     start.age = CALIB.START.AGE,
                                     max.age = CALIB.MAX.AGE),
                            error=function(e) {print(e)})
  sink()
  return(markov.result)
}

calculate.outputs <- function(par) {
  if(is.list(par)) {
    markov.result <- par
  } else {
    markov.result <- run.calib.simulation(par)
  }
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

plot.calibration.state <- function(pars) {
  n.params <- length(CALIB.PARAMS)
  other.pars <- pars[-seq(1,K*n.params,n.params)]
  ac.pars <- pars[seq(1,K*n.params,n.params)]
  
  # df <- data.frame(age=rep(1:8, each=2), val=other.pars, probs=c(rep(c('hsil_regression', 'hsil'), 8)))
  # plt <- ggplot(df, aes(x=age, y=val, color=probs)) + 
  #   geom_line() + 
  #   scale_x_continuous(breaks=1:8, labels=paste0(seq(40,75,5), '-', seq(44,79,5))) +
  #   theme_minimal() +
  #   theme(panel.grid.minor=element_blank())
  
  dff <- data.frame(age=rep(1:K), val=ac.pars, probs=c(rep(c('cancer'), K)))
  pltt <- ggplot(dff, aes(x=age, y=val, color=probs)) + 
    geom_line() +
    scale_x_continuous(breaks=1:K, labels=paste0(seq(CALIB.START.AGE,CALIB.MAX.AGE,5), '-', seq(CALIB.START.AGE+4,CALIB.MAX.AGE+4,5))) +
    theme_minimal() +
    theme(panel.grid.minor=element_blank())
  
  # browser()
  hsil.inc <- calculate.hsil.incidence(pars)
  df2 <- data.frame(age=rep(1:K, 2),
                    val=c(hsil.inc, rep(TARGET.INC.HSIL, K)),
                    type=c(rep('simulation', K), rep('target', K)))
  plt2 <- ggplot(df2[df2$age!=1,], aes(x=age, y=val, linetype=type)) + geom_line() +
    ylab('HSIL incidence') +
    ggtitle(paste0('Mean HSIL incidence: ', mean(hsil.inc))) +
    scale_x_continuous(breaks=1:K, labels=paste0(seq(DEFAULT.START.AGE$hiv_msm-5,DEFAULT.MAX.AGE$hiv_msm,5), '-', seq(DEFAULT.START.AGE$hiv_msm-1,DEFAULT.MAX.AGE$hiv_msm+4,5))) +
    theme_minimal() +
    # ylim(0, max(mean(hsil.inc), TARGET.INC.HSIL)) +
    theme(panel.grid.minor=element_blank())
  plt2
  
  ac.inc <- calculate.ac.incidence(pars)
  df3 <- data.frame(age=rep(1:K, 2),
                    val=c(ac.inc, TARGET.INC.AC),
                    linetype=c(rep('simulation', K), rep('target', K)))
  plt3 <- ggplot(df3[df3$age!=1,], aes(x=age, y=val, linetype=linetype)) + 
    geom_line() + ylab('Cancer incidence') +
    scale_x_continuous(breaks=1:K, labels=paste0(seq(CALIB.START.AGE,CALIB.MAX.AGE,5), '-', seq(CALIB.START.AGE+4,CALIB.MAX.AGE+4,5))) +
    theme_minimal() +
    theme(panel.grid.minor=element_blank())
  # grid.arrange(plt, pltt, plt2, plt3, ncol=1)
  grid.arrange(pltt, plt2, plt3, ncol=1)
}

# plot.calibration.state(c(0.007963232, 0.004910898, 0.006778633, 0.006961322, 0.007391491, 0.007806981,
#                          0.008209872, 0.008587149, 0.008995657))

# plot.calibration.state(c(0.007273147, 0.004666855, 0.004906142, 0.006724213, 0.007014540, 0.007467946,
#                          0.007822012, 0.008166392, 0.008589980, 0.008943130))

plot.calibration.state(sapply(strat.ctx, function(ctx) ctx$p_cancer___hsil_annual)[2:11])
browser()
initial.guess <- ctx.to.calib.vec(strat.ctx)

# initial.guess[seq(1,24,3)] <- rev(initial.guess[seq(1,24,3)])

### Calibration tests

eval.count <- 0
start.time <- Sys.time()
# res <- optim(initial.guess,
#              calibration.error,
#              method='L-BFGS-B',
#              lower=0,
#              upper=pmin(1, initial.guess*2),
#              control=list(parscale=rep(1e1, length(initial.guess))))

best.solution <- NULL
best.solution.val <- 1e10

plot.calibration.state(initial.guess)

# res <- cma_es(initial.guess,
#              calibration.error,
#              lower=initial.guess*.5,
#              upper=pmin(1, initial.guess*3),
#              control=list(parscale=rep(1e1, length(initial.guess))))


cat('Elapsed time: ', as.numeric(difftime(Sys.time(), start.time, units='min')), ' min')
cat('Evaluations: ', eval.count)

# CALIBRATED VALUES FOR NON-CONSTANT HSIL INCIDENCE
# calibrated.params <- c(
# 0.0033, 0.2259, 0.1757, 0.0026, 0.3416, 0.1799, 0.0032, 0.1599, 0.0632, 0.0029, 0.0762, 0.0915, 0.0025, 0.0806, 0.0549, 0.0031, 0.0840, 0.0283, 0.0033, 0.0278, 0.0362, 0.0022, 0.0167, 0.1626
# )

# calibrated.params <- res$par

uncalibrated.params <- c(
  sapply(strat.ctx[CALIB.AGE.GROUPS], function(ctx) ctx$p_cancer___hsil_annual)
  # rep(1/131, 8)
  # ,
  # sapply(strat.ctx[CALIB.AGE.GROUPS], function(ctx) ctx$p_hsil_regression_annual)
  # ,
  # sapply(strat.ctx[CALIB.AGE.GROUPS], function(ctx) ctx$p_hsil_annual)
)


# initial.output <- calculate.ac.incidence(ctx.to.calib.vec(strat.ctx))
# calibrated.output <- calculate.ac.incidence(calibrated.params)
#
# initial.output.hsil <- calculate.hsil.incidence(ctx.to.calib.vec(strat.ctx))
# calibrated.output.hsil <- calculate.hsil.incidence(calibrated.params)

initial.output <- calculate.ac.incidence(uncalibrated.params)
calibrated.output <- calculate.ac.incidence(ctx.to.calib.vec(strat.ctx))

initial.output.hsil <- calculate.hsil.incidence(uncalibrated.params)
calibrated.output.hsil <- calculate.hsil.incidence(ctx.to.calib.vec(strat.ctx))

# N <- length(initial.output)
# df <- data.frame(x=c(
#                      1:N,
#                      1:N,
#                      rep(1:N, 2),
#                      rep(1:N, 2),
#                      NULL
#                      ),
#                  error=c(
#                          TARGET.INC.AC,
#                          TARGET.INC.HSIL,
#                          initial.output, calibrated.output,
#                          initial.output.hsil, calibrated.output.hsil,
#                          NULL
#                          ) * 1e5,
#                  type=c(
#                         rep('target', N),
#                         rep('target', N),
#                         rep('initial', N), rep('calibrated', N),
#                         rep('initial', N), rep('calibrated', N),
#                         NULL
#                         ),
#                  measure=c(
#                            rep('ac', N),
#                            rep('hsil', N),
#                            rep('ac', N), rep('ac', N),
#                            rep('hsil', N), rep('hsil', N),
#                            NULL
#                            ))
# plt <- ggplot(df, aes(x=x, y=error, color=type, linetype=measure)) +
#   geom_line() +
#   scale_x_continuous(breaks=1:8, labels=paste0(40+((1:8)-1)*5, '-', 40+((2:9)-1)*5)) +
#   # coord_cartesian(ylim=c(0, 0.0012)) +
#   scale_color_manual(name='',
#                      breaks=c('target', 'initial', 'calibrated'),
#                      values=c('black', 'red', 'blue'),
#                      labels=c('Target (AC incidence)', 'Initial', 'Calibrated')) +
#   scale_linetype_manual(name='',
#                         breaks=c('ac', 'hsil'),
#                         values=c('solid', 'dashed'),
#                         labels=c('AC incidence', 'HSIL incidence')) +
#   xlab('Age group') +
#   ylab('Incidence (per 100,000)') +
#   theme_minimal() +
#   theme(panel.grid = element_blank())
# print(plt)
# ggplotly(plt)
# 
# x <- ctx.to.calib.vec(strat.ctx)
# print(paste0('Mean HSIL incidence: ', mean(calculate.hsil.incidence(x))))




# initial.output <- calculate.ac.incidence(ctx.to.calib.vec(strat.ctx))
# # print(initial.output)
# # initial.output.hsils <- calculate.hsil.prevalence(ctx.to.calib.vec(strat.ctx))
# N <- length(initial.output)
# df <- data.frame(x=rep(1:N,2),
#                  incidence=100000*c(TARGET.INC.AC, initial.output),
#                  type=c(rep('target', N), rep('initial', N)))
# plt <- ggplot(df, aes(x=x, y=incidence, color=type)) +
#   geom_line() +
#   scale_x_continuous(breaks=1:8, labels=paste0(40+((1:8)-1)*5, '-', 40+((2:9)-1)*5-1)) +
#   # coord_cartesian(ylim=c(0, 0.0012)) +
#   scale_color_manual(name='AC incidence',
#                      breaks=c('target', 'initial'),
#                      values=c('black', 'red')) +
#   theme_minimal()
# ggplotly(plt)
