library(dplyr)
library(gtools)

source('load_models.R')
source('markov.R')

VARIATION <- 3

# strategy <- 'arnme6e7_hpvhr_t_irc'
# reference <- 'conventional_t_irc'
strategy <- 'arnme6e7_hpvhr_t_tca'
# reference <- 'conventional_t_tca'
reference <- 'arnme6e7_hpvhr_t_irc'

initial.state <- sapply(markov$nodes,
                        function(n) if (n$name=='hiv_positive') 1 else 0)

sim.strats <- strategies$hiv_msm
sim.strats <- sim.strats[names(sim.strats) %in% c(strategy, reference)]

calculate.output <- function(par.val, par.name, age.group='-') {
  # strat.ctx.i <- lapply(strat.ctx, function(ctx) {
  #   if (is.null(age.group) || age.group == )
  #   ctx[par.name] <- par.val
  #   ctx
  # })
  strat.ctx.i <- lapply(names(strat.ctx), function(ctx.name) {
    ctx.i <- strat.ctx[[ctx.name]]
    if (age.group == '-' || age.group == ctx.name)
      ctx.i[par.name] <- par.val
    ctx.i
  })
  names(strat.ctx.i) <- names(strat.ctx)
  strat.ctx.i <- refresh.context(par, strat.ctx.i, excel.strata.df)

  x <- simulate('hiv_msm',
                sim.strats,
                markov,
                strat.ctx.i,
                initial.state,
                discount.rate = 0.03)$summary
  return(x)
}

calculate.nhb <- function(x) {
  st <- x[startsWith(x$strategy,strategy),]
  ref <- x[startsWith(x$strategy,reference),]
  nhb <- (st$E-ref$E) - (st$C-ref$C)/25000
  return(nhb)
}

compute.icer <- function(x) {
  st <- x[startsWith(x$strategy,strategy),]
  ref <- x[startsWith(x$strategy,reference),]
  icer <- (st$C-ref$C)/(st$E-ref$E)
  return(icer)
}

calculate.icer <- function(par.val, par.name, age.group='-') {
  sink('/dev/null')
  icer <- tryCatch({
    x <- calculate.output(par.val, par.name, age.group)
    if (is.na(x$IE[2])) icer <- 1e10
    else icer <- compute.icer(x)
    icer
  }, error=function(e) {1e10})
  sink()
  print(paste0(par.val, ': ', icer))
  return(icer)
}

icer.diff <- function(par.val, par.name, age.group='-') {
  sink('/dev/null')
  x <- calculate.output(par.val, par.name, age.group)
  sink()
  # st <- x[startsWith(x$strategy,strategy),]
  # ref <- x[startsWith(x$strategy,reference),]
  
  
  error <- (calculate.nhb(x))^2
  print(paste0(par.val, ': ', compute.icer(x)))
  # error <- ((st$C-ref$C)/(st$E-ref$E) - 25000)^2
  # print((st$C-ref$C)/(st$E-ref$E))
  return(error)
}

# par.list <- c('c_arn_kit')
# par.list <- pars[startsWith(pars, 'u_')]
# par.list <- c('u_hiv_p',
#               'p_cyto_b___hsil',
#               'p_arnmhr_p___hsil',
#               'p_undetected_hsil_treatment_whole_followup',
#               'u_cancer_delayed',
#               'c_hra_treatment',
#               'p_arnmhr_p___no_hsil',
#               'p_cyto_b___no_hsil',
#               'c_arn_kit',
#               'p_no_hsil___hsil_tca',
#               'survival_5year',
#               'p_hsil_regression_annual',
#               'c_surgery_delayed',
#               'p_hra_hsil___cyto_hsil__hsil')

# par.list <- c('u_hsil___irc', 'p_cyto_b___hsil', 'p_undetected_hsil_treatment_whole_followup', 'p_arnmhr_p___hsil')
# par.list <- c(
#   'c_cyto',
#   'c_arn_kit',
#   'c_hra_treatment',
#   'c_irc',
#   'c_surgery',
#   'c_surgery_delayed',
#   'p_cyto_b___hsil',
#   'p_cyto_b___no_hsil',
#   'p_cyto_hsil___no_hsil',
#   'p_death_other_annual',
#   'p_arnmhr_p___hsil',
#   'p_arnmhr_p___no_hsil',
#   'p_hra_hsil___cyto_hsil__hsil',
#   'p_hsil_regression_annual',
#   'p_undetected_hsil_treatment_whole_followup',
#   'survival_5year',
#   'p_no_hsil___hsil_irc',
#   'u_hiv_p',
#   'u_hsil___irc',
#   'u_cancer',
#   'u_cancer_delayed'
# )

# Age-specific
# par.list <- c('p_hsil_regression_annual')


# conventional vs arn-hr (TCA)
# par.list <- c('c_arn_kit', 'u_hiv_p_cd4_g500', 'p_cd4_g500', 'p_cyto_b___hsil',
#               'p_cancer___hsil_annual', 'p_arnmhr_p___hsil', 'p_undetected_hsil_treatment_whole_followup')

# arn-hr (IRC) vs arn-hr (TCA)
par.list <- c('p_no_hsil___hsil_tca', 'p_no_hsil___hsil_irc')

# ARN-HPV-HR
# par.list <- c('u_hsil___irc', 'u_hsil___tca', 'u_hiv_p___irc', 'u_hiv_p___tca',
#               'p_no_hsil___hsil_irc', 'p_no_hsil___hsil_tca')


df <- data.frame()

for(par in par.list) {
  if (length(unique(sapply(strat.ctx, function(ctx)ctx[[par]]))) == 1) {
    age.groups <- '-'
  } else {
    age.groups <- names(strat.ctx)[age.groups >= DEFAULT.START.AGE$hiv_msm & age.groups < DEFAULT.MAX.AGE$hiv_msm]
  }
  for(age.group in age.groups) {
    # if (age.group == '-') age.group <- NULL
    if (age.group == '-') {
      initial.guess <- strat.ctx$y40_44[[par]]
      cat('Finding', par, 'value with NHB~0...\n')
    } else {
      initial.guess <- strat.ctx[[age.group]][[par]]
      cat('Finding', par, 'value for', age.group, 'with NHB~0...\n')
    }
    
    suppressWarnings(sink())
    # sink('/dev/null')
    lower.limit <- max(0, initial.guess*(1-VARIATION))
    upper.limit <- ifelse(!any(startsWith(par, c('p', '.p', 'u'))),
                           initial.guess*(1+VARIATION),
                           min(initial.guess*(1+VARIATION), 1))
    
    # opt.val <- tryCatch({
    #   res <- binsearch(calculate.icer, c(lower.limit, upper.limit), target=25000, par.name=par, age.group=age.group)
    #   opt.val <- mean(res$where)
    #   opt.val
    # }, error=function(e) {
      res <- optimize(icer.diff,
                      lower=lower.limit,
                      upper=upper.limit,
                      par.name=par,
                      age.group=age.group,
                      tol=1e-2)
      opt.val <- res$minimum
      opt.val
    # })
    
    sink('/dev/null')
    optim.output <- calculate.output(opt.val, par, age.group=age.group)
    sink()
    
    df <- rbind(df, data.frame(par=par,
                               age.group=age.group,
                               # sq.nhb.diff=res$objective,
                               # nhb=calculate.nhb(optim.output),
                               icer=compute.icer(optim.output),
                               base.value=initial.guess,
                               critical.value=opt.val))
    print(df)
  }
}

print(df)
# View(df)
# df <- df[,c('par', 'icer', 'base.value', 'value')]
write.xlsx(df, 'output/results/critical_values.xlsx')
