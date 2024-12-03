library(CEAModel)
library(xlsx)
library(tools)

source('distributions.R')
source('excel_params.R')

### CALIBRATION SETTINGS (see below)
UPDATE.CALIBRATED.PARAMS <- T
###

load.all.trees <- function() {
  trees <- list()
  for(f in list.files('models/')) {
    if (!startsWith(f, '_') && endsWith(f, '.yaml')) {
      if (endsWith(f, '_t.yaml')) {
        # Two versions of each treatment strategy (TCA + IRC)
        name <- paste0(substr(f,1,(nchar(f) - 5)), '_tca')
        t <- loadDecisionTree(paste0('models/', f))
        t$name <- name
        # assign(name, t, envir = .GlobalEnv)
        trees[[name]] <- t

        name <- paste0(substr(f,1,(nchar(f) - 5)), '_irc')
        t <- loadDecisionTree(paste0('models/', f))
        t$name <- name
        # assign(name, t, envir = .GlobalEnv)
        trees[[name]] <- t
      } else {
        # One version (followup)
        name <- substr(f,1,(nchar(f) - 5))
        t <- loadDecisionTree(paste0('models/', f))
        # assign(name, t, envir = .GlobalEnv)
        trees[[name]] <- t
      }
    }
  }
  return(trees)
}

cat('Loading trees...\n')
trees <- load.all.trees()

strategies <- list()
strategies$hiv_msm <- trees[
                          names(trees) == 'no_intervention' |
                          startsWith(names(trees), 'conventional') |
                          startsWith(names(trees), 'arnm') |
                          startsWith(names(trees), 'ascus_lsil')
                      ]

# strategies$hiv_msm <- trees[c('no_intervention', 'conventional')]
# strategies$hiv_msm <- trees[1:5]

cat('Loading context(s)...\n')
# Context loading
context.setup <- function(strat.ctx) {
  strat.ctx <- lapply(strat.ctx, function(ctx) {
    return(ctx)
  })
  return(strat.ctx)
}

strat.ctx.path <- 'params/context.xlsx'
strat.ctx.hash.file <- paste0(strat.ctx.path, '.hash')
if (file.exists(strat.ctx.hash.file) && readLines(strat.ctx.hash.file) == md5sum(strat.ctx.path)) {
  full.strat.ctx <- readRDS('params/_full.strat.ctx.RData')
  full.strat.metadata <- readRDS('params/_full.strat.metadata.RData')
  excel.strata.df <- readRDS('params/_excel.strata.df.RData')
  excel.strata.df.full <- readRDS('params/_excel.strata.df.full.RData')
  strat.ctx <- readRDS('params/_strat.ctx.RData')
} else {
  full.strat.data <- loadStratifiedContextFile(strat.ctx.path)
  full.strat.ctx <- full.strat.data$ctx
  full.strat.metadata <- full.strat.data$metadata
  excel.strata.df <- list()
  excel.strata.df.full <- list()
  .strata <- names(xlsx::getSheets(xlsx::loadWorkbook(strat.ctx.path)))
  for(stratum in .strata) {
    excel.strata.df[[stratum]] <- read.xlsx(strat.ctx.path, sheetName = stratum, keepFormulas = T)[,c(1,2)]
    excel.strata.df.full[[stratum]] <- read.xlsx(strat.ctx.path, sheetName = stratum, keepFormulas = T)
  }
  strat.ctx <- refresh.context(c(), full.strat.ctx, excel.strata.df, context.setup)

  # Save data for faster loading next time if unchanged
  saveRDS(full.strat.ctx, file='params/_full.strat.ctx.RData')
  saveRDS(full.strat.metadata, file='params/_full.strat.metadata.RData')
  saveRDS(excel.strata.df, file='params/_excel.strata.df.RData')
  saveRDS(excel.strata.df.full, file='params/_excel.strata.df.full.RData')
  saveRDS(strat.ctx, file='params/_strat.ctx.RData')
  cat(paste0(md5sum(strat.ctx.path), '\n'), file=strat.ctx.hash.file)
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

get.calibrated.params <- function() {
  return(c('p_cancer___hsil_annual'))
}

update.calibrated.params <- function(strat.ctx) {
  DEFAULT.CALIBRATED.START.AGE <- 40
  DEFAULT.CALIBRATED.MAX.AGE <- 80
  age.groups <- as.numeric(substr(names(strat.ctx), 2, 3))
  CALIB.AGE.GROUPS <- names(strat.ctx)[age.groups >= DEFAULT.CALIBRATED.START.AGE & age.groups < DEFAULT.CALIBRATED.MAX.AGE]
  CALIB.PARAMS <- c('p_cancer___hsil_annual')
  CALIB.VALUES <- c(0.0058691355195081, 0.00344980252811424, 0.00324775583133558, 0.00258766335262408, 0.00193687668776452, 0.0018456703895995, 0.00184691109687812, 0.0018916537948187)

  strat.ctx <- calib.vec.to.ctx(CALIB.VALUES, strat.ctx)
  refresh.context(CALIB.PARAMS, strat.ctx, excel.strata.df)
  return(strat.ctx)
}

if (UPDATE.CALIBRATED.PARAMS) {
  cat('Calibrated parameters enabled, updating context...\n')
  strat.ctx <- update.calibrated.params(strat.ctx)
}

base.ctx <- lapply(strat.ctx[[1]], function(l)l[1])

independent.pars <- getIndependentParameters(strat.ctx.path)

cat('Loading markov model...\n')
markov <- loadMarkovModels(paste0('models/markov.xlsx'))


# for(s in strategies) {
#   plt <- s$show(strat.ctx$y40_44, prevalence=1, nodeInfoFields = c('path.prob.100k', 'cost'))
#   visNetwork::visSave(plt, paste0('output/hsil_', s$name, '.html'), selfcontained = TRUE)
#   plt <- s$show(strat.ctx$y40_44, prevalence=0, nodeInfoFields = c('path.prob.100k', 'cost'))
#   visNetwork::visSave(plt, paste0('output/no_hsil_', s$name, '.html'), selfcontained = TRUE)
# }
