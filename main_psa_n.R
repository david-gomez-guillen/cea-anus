library(officer)
library(rvg)
library(magrittr)

# setwd('~/Documents/models_ce/anus')

source('load_models.R')
source('markov_psa.R')

PSA.SEED <- 12345
N.ITERS <- 1000
N.CORES <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(N.CORES)) N.CORES <- 8
DISCOUNT.RATE <- .03
DEBUG <- T  # 1 core if TRUE

JITTER.X <- .05
JITTER.Y <- .05

EXCLUDED.PARAMS <- names(strat.ctx$y25_29[strat.ctx$y25_29 %in% c(0,1)])
EXCLUDED.PARAMS <- c(EXCLUDED.PARAMS, 
                     'periodicity_months', 'periodicity_times_in_year', 'p_vaccination', 'p_coverage'
)

pars <- independent.pars
pars <- pars[!pars %in% EXCLUDED.PARAMS]

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  n.cores <- N.CORES
} else {
  n.cores <- as.numeric(args[1])
}

GET.SD.FUNC <- function(divisor) {
  function(p.name, value) value/divisor
}

SD.ESTIMATE.FUNCTIONS <- list(
  # sd_10=GET.SD.FUNC(10)
  # ,
  sd_5=GET.SD.FUNC(5)
  # ,
  # sd_2=GET.SD.FUNC(2)
  # ,
  # sd_1=GET.SD.FUNC(1)
  # ,
  # medium=function(par.name, value) {
  #   # These parameters are considered less trustworthy than the rest
  #   if (par.name %in% c('.pipelle_success',
  #                  '.sensitivity_hysteroscopy_bleeding',
  #                  '.sensitivity_molecular',
  #                  '.sensitivity_pipelle')) {
  #     return(value/10)
  #   } else {
  #     return(value/20)
  #   }
  # }
)

SIMULATION.OPTIONS <- list(
  cito_vs_cotest_arn=list(
   population='hiv_msm',
   reference='conventional_t_tca',
   strategy='arnme6e7_hpvhr_t_tca',
   reference.name='Conventional (TCA)',
   strategy.name='Co-test mRNA (TCA)'
  )
  #,
  # cito_vs_cotest_pcr=list(
  #  population='hiv_msm',
  #  reference='conventional_t_tca',
  #  strategy='conventional_hpvhrla_t_tca',
  #  reference.name='Conventional (TCA)',
  #  strategy.name='Co-test PCR (TCA)'
  # )
  # ,
  #  cotest_pcr_irc_vs_tca=list(
  #    population='hiv_msm',
  #    reference='conventional_hpvhrla_t_irc',
  #    strategy='conventional_hpvhrla_t_tca',
  #    reference.name='Co-test PCR (IRC)',
  #    strategy.name='Co-test PCR (TCA)'
  #  )
  # ,
  #  cotest_arn_irc_vs_tca=list(
  #   population='hiv_msm',
  #   reference='arnme6e7_hpvhr_t_irc',
  #   strategy='arnme6e7_hpvhr_t_tca',
  #   reference.name='ARNmE6/E7 HPV-HR (IRC)',
  #   strategy.name='ARNmE6/E7 HPV-HR (TCA)'
  #  )
)

psa.pars <- list(
  all=pars
  # ,
  # all_no_costs=pars[!startsWith(pars, 'c')]
  # ,
  # costs=pars[startsWith(pars, 'c_')]
  # ,
  # utilities=pars[startsWith(pars, 'u_')]
  # ,
  # probs=pars[startsWith(pars, 'p_') |
  #              startsWith(pars, 'survival_')]
  # ,
  # cd4=list(list('p_cd4_l200', 'p_cd4_200_500', 'p_cd4_g500'))
  # sensitivities=pars[(startsWith(pars, 'p_hpv') | startsWith(pars, 'p_arn')) & endsWith(pars, '___hsil')]
  # ,
  # specificities=pars[(startsWith(pars, 'p_hpv') | startsWith(pars, 'p_arn')) & endsWith(pars, '___no_hsil')]
  # ,
  # sensitivities_specificities=pars[(startsWith(pars, 'p_hpv') | startsWith(pars, 'p_arn')) & (endsWith(pars, '___no_hsil') | endsWith(pars, '___hsil'))]
  # ,
  # irc_tca=c('c_irc', 'c_tca_single', 'n_tca', 'p_arnmhr_p___hsil', 'p_arnmhr_p___no_hsil', 
  #                 'p_no_hsil___hsil_irc', 'p_no_hsil___hsil_tca', 'u_hsil___irc', 'u_hsil___tca')
  # ,
  # tornados=c('u_hiv_p_cd4_g500', 'p_cd4_g500', 'p_cyto_b___hsil',
  #            'p_arnmhr_p___hsil', 'p_undetected_hsil_treatment_whole_followup')
  # arn_sens_spec=c('p_arnmhr_p___hsil', 'p_arnmhr_p___no_hsil')
  # critical=c('p_cancer___hsil_annual', 'p_arnmhr_p___hsil', 'u_hiv_p_cd4_g500')
  # critical=c('p_cancer___hsil_annual', 'p_hpvhrla_p___hsil', 'p_undetected_hsil_treatment_whole_followup', 'u_hiv_p_cd4_g500')
)
# for(p in pars[!pars %in% psa.pars$tornados]) {
#   psa.pars[[paste0('tornado_plus_', p)]] <- c(psa.pars$tornados, p)
# }
# psa.pars$tornados <- NULL


# print(strat.ctx$y40_44$p_cancer___hsil_annual)
# for(stratum in names(strat.ctx)) {
#   strat.ctx[[stratum]]$p_cancer___hsil_annual <- strat.ctx[[stratum]]$p_cancer___hsil_annual * 2
# }
# print(strat.ctx$y40_44$p_cancer___hsil_annual)


#### MULTIVARIATE PSAS

RESUME.PAR <- ''

trees.all <- trees

if (DEBUG) {
  n.cores <- 1
  cl <- NULL
} else {
  out.file <- ifelse(Sys.getenv('SLURM_JOB_ID') == '', 'workers.out', paste0('workers-', Sys.getenv('SLURM_JOB_ID'), '.out'))
  cl <- makeCluster(n.cores, outfile=out.file)
  on.exit({
    stopCluster(cl)
  }, add=TRUE)
  registerDoParallel(cl)
env <- foreach:::.foreachGlobals
rm(list=ls(name=env), pos=env)  
  varlist <- ls(pos=1)
  # varlist <- varlist[!varlist %in% c('results', 'trees')]
  clusterExport(cl=cl,varlist=varlist[(varlist %in% c('conventional', 'ascus_lsil_diff_arnme6e7_16')) | (!startsWith(varlist, 'conventional_') & !startsWith(varlist, 'ascus_lsil_diff') & !startsWith(varlist, 'arnme6e7') & !varlist %in% c('trees', 'strategies'))],
                envir=environment())
  clusterEvalQ(cl=cl,{
    library(data.table)
    library(stringr)
    library(ggplot2)
    library(plotly)
    library(pbapply)
    library(gtools)
  })
}

for(param.set.name in names(psa.pars)) {
  param.set <- psa.pars[[param.set.name]]
  if (RESUME.PAR != '')
    param.set <- param.set[match(RESUME.PAR, param.set):length(param.set)]
  
  for(option.name in names(SIMULATION.OPTIONS)) {
    trees <- trees.all
    options <- SIMULATION.OPTIONS[[option.name]]
    output.dir <- paste0(getwd(), '/output/psa_multivariate/', options$population)
    
    if (!is.null(cl)) {
      # We prune the 'tree' list to avoid copying it to each process
      trees <- trees[names(trees) %in% c(options$strategy, options$reference) | startsWith(names(trees), 'semestral_')]
      clusterExport(cl=cl,varlist='trees',
                    envir=environment())
    }
    
    for(sd.estimate.name in names(SD.ESTIMATE.FUNCTIONS)) {
      sd.estimate.func <- SD.ESTIMATE.FUNCTIONS[[sd.estimate.name]]
      csv.data <- data.frame()
      filename <- paste0(options$strategy, '__sd_', sd.estimate.name, '__par_', param.set.name)
      
      if (!is.null(cl))
        clusterExport(cl=cl,varlist=c('param.set', 'options', 'sd.estimate.func'),
                      envir=environment())
      
      results <- psa.n(param.set, 
                       strat.ctx,
                       options$population, 
                       options$strategy,
                       options$reference,
                       markov,
                       excel.file = 'params/context.xlsx',
                       sd.estimate.func=sd.estimate.func,
                       context.setup.func=context.setup,
                       n.cores = n.cores,
                       cluster = cl,
                       n.iters=N.ITERS,
                       seed=PSA.SEED,
                       discount.rate=DISCOUNT.RATE,
                       jitter.x=JITTER.X,
                       jitter.y=JITTER.Y)
      
      # Remove constant parameters
      # for(col in names(results$summary)) {n
      #   if (length(unique(results$summary[[col]])) == 1)
      #     results$summary[[col]] <- NULL
      # }
      # Remove all strata but one for non-stratified parameters
      # par.cols <- names(results$summary)
      # for(p in param.set) {
      #   p.cols <- par.cols[endsWith(par.cols, paste0('_',p))]
      #   p.df <- results$summary[p.cols]
      #   p.is.stratified <- !all(sapply(1:ncol(p.df), function(c) all(p.df[1]==p.df[c])))
      #   if (!p.is.stratified)
      #     results$summary[p.cols[2:length(p.cols)]] <- NULL
      # }
# browser()
      # filename <- paste0(output.dir, '/', filename, '.csv')
      # write.csv(results$summary, filename)
      # results <- NULL
      # 
      # df.results <- data.frame()
      # for(f in list.files(output.dir, pattern = paste0(filename, '\\d*.csv'))) {
      #   df.results <- rbind(df.results, read.csv(paste0(output.dir, '/', f)))
      #   unlink(paste0(output.dir, '/', f))
      # }
      # results <- generate.psa.summary(df.results, options$strategy, options$population, options$reference, strat.ctx, param.set, jitter.x=JITTER.X, jitter.y=JITTER.Y)

      #'output/psa_', psa.type, '/', sim.options$population, '/', strategy, '__', option.name, '__sd_', sd.estimate.name, '__par_', suffix, '.pptx'
      filename <- paste0(options$strategy, '__', option.name, '__sd_', sd.estimate.name, '__par_', param.set.name)
      store.results.psa(results, 'multivariate', options$population, options$strategy.name, filename, discount=DISCOUNT.RATE)

      build.plots('multivariate', options$population, options$strategy, options$reference, strat.ctx, option.name, options, sd.estimate.name, suffix=param.set.name, param.set.name=param.set.name)
    }
  }
}


