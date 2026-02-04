library(officer)
library(rvg)
library(magrittr)

# setwd('~/Documents/models_ce/anus')

source('load_models.R')
source('markov_psa.R')

PSA.SEED <- 12345
N.ITERS <- 1000
N.CORES <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(N.CORES)) N.CORES <- 2
DISCOUNT.RATE <- .03
DEBUG <- F  # 1 core if TRUE

JITTER.X <- .05
JITTER.Y <- .05

EXCLUDED.PARAMS <- names(strat.ctx$y25_29[strat.ctx$y25_29 %in% c(0,1)])
EXCLUDED.PARAMS <- c(EXCLUDED.PARAMS, 
                     'periodicity_months', 'periodicity_times_in_year', 'p_vaccination', 'p_coverage'
)

pars <- independent.pars
pars <- pars[!pars %in% EXCLUDED.PARAMS]


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
  # followup=list(
  #   population='hiv_msm',
  #   reference='conventional_t_tca',
  #   strategy='arnme6e7_hpvhr_t_tca',
  #   reference.name='Conventional (TCA)',
  #   strategy.name='ARNmE6/E7 HPV-HR (TCA)'
  # )
  # ,
  # treatment=list(
  #   population='hiv_msm',
  #   reference='conventional_t_tca',
  #   strategy='ascus_lsil_diff_arnme6e7_16_t_tca',
  #   reference.name='Conventional (TCA)',
  #   strategy.name='ASCUS/LSIL diff - ARNME6E7 16 (TCA)'
  # )
  # ,
  cito_vs_cotest_pcr=list(
    population='hiv_msm',
    reference='conventional_t_tca',
    strategy='conventional_hpvhrla_t_tca',
    reference.name='Cytology + TCA',
    strategy.name='Co-testing DNA-PCR + TCA'
  )
)

psa.pars <- list(
  # all=pars
  # ,
  # costs=pars[startsWith(pars, 'c_')]
  # ,
  # utilities=pars[startsWith(pars, 'u_')]
  # ,
  # probs=pars[startsWith(pars, 'p_') |
  #              startsWith(pars, 'survival_')]
  # ,
  test=c('rate_death_other_annual')
)



#### UNIVARIATE PSAS

RESUME.PAR <- ''

trees.all <- trees

if (DEBUG) {
  N.CORES <- 1
  cl <- NULL
} else {
  out.file <- ifelse(Sys.getenv('SLURM_JOB_ID') == '', 'workers.out', paste0('workers-', Sys.getenv('SLURM_JOB_ID'), '.out'))
  cl <- makeCluster(N.CORES, outfile=out.file)
  on.exit({
    stopCluster(cl)
  }, add=TRUE)
  registerDoParallel(cl)
  
  varlist <- ls(pos=1)
  varlist <- varlist[!varlist %in% c('trees', 'strategies')]
  
  clusterExport(cl=cl,varlist=varlist,
                envir=environment())
  clusterEvalQ(cl=cl,{
    library(data.table)
    library(stringr)
    library(ggplot2)
    library(plotly)
    library(pbapply)
  })
}

for(param.set.name in names(psa.pars)) {
  param.set <- psa.pars[[param.set.name]]
  if (RESUME.PAR != '')
	param.set <- param.set[match(RESUME.PAR, param.set):length(param.set)]

  trees <- trees.all
  for(option.name in names(SIMULATION.OPTIONS)) {
    options <- SIMULATION.OPTIONS[[option.name]]
    output.dir <- paste0(getwd(), '/output/psa_univariate/', options$population)
    
    if (!is.null(cl)) {
      # We prune the 'tree' list to avoid copying it to each process
      trees <- trees[names(trees) %in% c(options$strategy, options$reference) | startsWith(names(trees), 'semestral_')]
      clusterExport(cl=cl,varlist='trees',
                    envir=environment())
    }
    
    for(sd.estimate.name in names(SD.ESTIMATE.FUNCTIONS)) {
      sd.estimate.func <- SD.ESTIMATE.FUNCTIONS[[sd.estimate.name]]
      for(par in param.set) {
        filename.preffix <- paste0(options$strategy, '__sd_', sd.estimate.name, '__par_', par, '__')
        results <- psa.n(par,
                         strat.ctx,
                         options$population,
                         options$strategy,
                         options$reference,
                         markov,
                         excel.file = 'params/context.xlsx',
                         sd.estimate.func=sd.estimate.func,
                         context.setup.func=context.setup,
                         n.cores = N.CORES,
                         cluster = cl,
                         # n.cores = 1,
                         # cluster = NULL,
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
        # filename <- paste0(output.dir, '/', filename.preffix, '.csv')
        # write.csv(results$summary, filename)
        # results <- NULL

        # df.results <- data.frame()
        # for(f in list.files(output.dir, pattern = paste0(filename.preffix, '\\d*.csv'))) {
        #   df.results <- rbind(df.results, read.csv(paste0(output.dir, '/', f)))
        #   unlink(paste0(output.dir, '/', f))
        # }
        # results <- generate.psa.summary(df.results, options$strategy, options$population, options$reference, strat.ctx, param.set, jitter.x=JITTER.X, jitter.y=JITTER.Y)

        filename <- paste0(options$strategy, '__', option.name, '__sd_', sd.estimate.name, '__par_', par)
        store.results.psa(results, 'univariate', options$population, options$strategy.name, filename)
      }
      
      build.plots('univariate', options$population, options$strategy, options$reference, strat.ctx, option.name, options, sd.estimate.name, suffix=param.set.name, pars=param.set)
    }
  }
}





