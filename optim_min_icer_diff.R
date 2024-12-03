library(dplyr)

source('load_models.R')
source('markov.R')

VARIATION <- .2

strategy <- 'arnme6e7_hpvhr_t_tca'
reference <- 'conventional_t_tca'

initial.state <- sapply(markov$nodes,
                        function(n) if (n$name=='hiv_positive') 1 else 0)

sim.strats <- strategies$hiv_msm
sim.strats <- sim.strats[names(sim.strats) %in% c(strategy, reference)]

icer.diff <- function(par.val, par.name) {
  strat.ctx.i <- lapply(strat.ctx, function(ctx) {
    ctx[par.name] <- par.val
    ctx
  })
  strat.ctx.i <- refresh.context(par, strat.ctx.i, excel.strata.df)

  x <- simulate('hiv_msm',
                sim.strats,
                markov,
                strat.ctx.i,
                initial.state,
                discount.rate = 0.03)$summary

  st <- x[startsWith(x$strategy,strategy),]
  ref <- x[startsWith(x$strategy,reference),]
  nhb.diff <- ((st$E-ref$E) - (st$C-ref$C)/25000)^2
  # icer.diff <- ((st$C-ref$C)/(st$E-ref$E) - 25000)^2
  # cat(par.val, nhb.diff, icer.diff, '\n')
  # print(x)
  return(nhb.diff)
}

par.list <- c('u_hiv_p', 'p_cyto_b___hsil', 'p_arnmhr_p___hsil', 'p_undetected_hsil_treatment_whole_followup',
              'u_cancer_delayed', 'c_hra_treatment', 'p_arnmhr_p___no_hsil', 'p_cyto_b___no_hsil',
              'c_arn_kit', 'p_no_hsil___hsil_tca', 'survival_5year', 'p_hsil_regression_annual', 'c_surgery_delayed',
              'p_hra_hsil___cyto_hsil__hsil')

df <- data.frame()

for(par in par.list) {
  cat('Finding', par, 'value with NHB~0...\n')
  initial.guess <- strat.ctx$y40_44[[par]]

  suppressWarnings(sink())
  sink('/dev/null')
  lower.limit <- initial.guess*(1-VARIATION)
  upper.limit <- ifelse(startsWith(par, 'c_'),
                         initial.guess*(1+VARIATION),
                         min(initial.guess*(1+VARIATION), 1))

  res <- optimize(icer.diff, lower=lower.limit, upper=upper.limit, par.name=par)
  sink()

  df <- rbind(df, data.frame(par=par, sq.nhb.diff=res$objective, base.value=initial.guess, value=paste0(ifelse(initial.guess-res$minimum>0, '<=', '>='), res$minimum)))
}

print(df)
