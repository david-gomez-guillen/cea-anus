library(dplyr)
library(ggplot2)

# setwd('~/Documents/models_ce/anus')

source('load_models.R')
source('markov.R')

set.seed(123)

age.groups <- as.numeric(substr(names(strat.ctx), 2, 3))
CALIB.AGE.GROUPS <- names(strat.ctx)[age.groups >= DEFAULT.START.AGE$hiv_msm & age.groups < DEFAULT.MAX.AGE$hiv_msm]
# CALIB.PARAMS <- c('p_cancer___hsil_annual', 'p_hsil_regression_annual', 'p_hsil_annual', 'survival_5year')
CALIB.PARAMS <- c('p_cancer___hsil_annual')

# PLACEHOLDER
# TARGET.INC.HSIL <- c(0.002489277, 0.002469715, 0.002304056, 0.002192828, 0.001754127, 0.001508943, 0.001389184, 0.001325221)
# TARGET.INC.HSIL <- c(0.013300098, 0.016618926, 0.019937754, 0.023256581, 0.03533414, 0.047411699, 0.059489258, 0.071566816)


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
x <- 40:79
df <- data.frame(x=x, y=lin.interp(x, target.df), group=c(sapply(1:8, function(i) rep(i, 5))))
df <- df %>% group_by(group) %>% summarise_at(vars(y), list(inc.avg=mean))
TARGET.INC.AC <- df$inc.avg
TARGET.INC.AC <- lin.interp(seq(40,79,5), target.df)
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
    
    ac.inc <- calculate.ac.incidence(result)
    error <- sum((((ac.inc - TARGET.INC.AC)/TARGET.INC.AC)^2))
  } else {
    error <- 1e10
  }
  print(error)
  return(error)
}

calculate.ac.incidence <- function(par) {
  if(is.list(par)) {
    markov.result <- par
  } else {
    markov.result <- run.calib.simulation(par)
  }
  # browser()
  output <- markov.result$info[[1]]$additional.info[c('year', 'incidence_cancer')]
  
  output$age.group <- cut(output$year, c(seq(40,79,5), 200), right=FALSE, labels=FALSE)
  
  agg.output <- output %>% 
    group_by(age.group) %>% 
    summarise_at(c('incidence_cancer'), mean) * 2 # Data is semestral, annual mean is 2 * semestral mean
  return(agg.output$incidence_cancer)
}


calculate.hsil.incidence <- function(par) {
  if(is.list(par)) {
    markov.result <- par
  } else {
    markov.result <- run.calib.simulation(par)
  }
  output <- markov.result$info[[1]]$additional.info[c('year', 'incidence_hsil')]
  
  output$age.group <- cut(output$year, c(seq(40,79,5), 200), right=FALSE, labels=FALSE)
  
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
                                     discount.rate=.03),
                            error=function(e) NULL)
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

initial.guess <- ctx.to.calib.vec(strat.ctx)

### Calibration tests

# eval.count <- 0
# start.time <- Sys.time()
# res <- optim(initial.guess,
#              calibration.error,
#              method='L-BFGS-B',
#              lower=0,
#              upper=pmin(1, initial.guess*2),
#              control=list(parscale=rep(1e0, length(initial.guess))))
# cat('Elapsed time: ', as.numeric(difftime(Sys.time(), start.time, units='min')), ' min')
# cat('Evaluations: ', eval.count)

# calibrated.params <- c(-0.0139209272984589, 0.21232154057876, 0.364071566588375, 0.717170166876237, -0.0622384853973795, 0.314116846940157, 0.186051428946503, 0.740628449313684, -0.135785937725696, 0.14925722025157, 0.243267527370599, 0.685908114352199, -0.0893746736930389, 0.0358020274651913, 0.0126455789263969, 0.650912765904558, -0.0330301253200985, 0.0647795838910365, 0.0927301615440421, 0.674989033250804, 0.168174095163349, 0.0586796353790068, 0.0909229577061081, 0.624857773495392, -0.149165565705991, 0.06126115901712, 0.0809641281482692, 0.662236155916717, 0.0400284915095845, 0.0555599660575189, 0.106234844080886, 0.582454650757503)
# calibrated.params <- res$par

calibrated.params <- CALIB.VALUES

initial.output <- calculate.ac.incidence(ctx.to.calib.vec(strat.ctx))
calibrated.output <- calculate.ac.incidence(calibrated.params)

initial.output.hsil <- calculate.hsil.incidence(ctx.to.calib.vec(strat.ctx))
calibrated.output.hsil <- calculate.hsil.incidence(calibrated.params)

N <- length(initial.output)
df <- data.frame(x=c(
                     # 1:N,
                     # rep(1:N, 2),
                     rep(1:N, 2),
                     NULL
                     ),
                 error=c(
                         # TARGET.INC.AC,
                         # initial.output, calibrated.output,
                         initial.output.hsil, calibrated.output.hsil,
                         NULL
                         ) * 1e5,
                 type=c(
                        # rep('target', N),
                        # rep('initial', N), rep('calibrated', N),
                        rep('initial', N), rep('calibrated', N),
                        NULL
                        ),
                 measure=c(
                           # rep('ac', N),
                           # rep('ac', N), rep('ac', N),
                           rep('hsil', N), rep('hsil', N),
                           NULL
                           ))
plt <- ggplot(df, aes(x=x, y=error, color=type, linetype=measure)) +
  geom_line() +
  scale_x_continuous(breaks=1:8, labels=paste0(40+((1:8)-1)*5, '-', 40+((2:9)-1)*5)) +
  # coord_cartesian(ylim=c(0, 0.0012)) +
  scale_color_manual(name='',
                     breaks=c('target', 'initial', 'calibrated'),
                     values=c('black', 'red', 'blue'),
                     labels=c('Target (AC incidence)', 'Initial', 'Calibrated')) +
  scale_linetype_manual(name='',
                        breaks=c('ac', 'hsil'),
                        values=c('solid', 'dashed'),
                        labels=c('AC incidence', 'HSIL incidence')) +
  xlab('Age group') +
  ylab('Incidence (per 100,000)') +
  theme_minimal()
plt
ggplotly(plt)




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
