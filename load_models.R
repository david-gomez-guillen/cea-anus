library(CEAModel)
library(xlsx)
library(tools)

source('distributions.R')
source('excel_params.R')

load.all.trees <- function() {
  trees <- list()
  for(f in list.files('models/')) {
    if (!startsWith(f, '_') && endsWith(f, '.yaml')) {
      if (endsWith(f, '_t.yaml')) {
        # Two versions of each treatment strategy (TCA + IRC)
        name <- paste0(substr(f,1,(nchar(f) - 5)), '_tca')
        t <- loadDecisionTree(paste0('models/', f))
        t$name <- name
        assign(name, t, envir = .GlobalEnv)
        trees[[name]] <- t
        
        name <- paste0(substr(f,1,(nchar(f) - 5)), '_irc')
        t <- loadDecisionTree(paste0('models/', f))
        t$name <- name
        assign(name, t, envir = .GlobalEnv)
        trees[[name]] <- t
      } else {
        # One version (followup)
        name <- substr(f,1,(nchar(f) - 5))
        t <- loadDecisionTree(paste0('models/', f))
        assign(name, t, envir = .GlobalEnv)
        trees[[name]] <- t
      }
    }
  }
  return(trees)
}

cat('Loading trees...\n')
trees <- load.all.trees()

strategies <- trees[startsWith(names(trees), 'conventional') |
                      startsWith(names(trees), 'arnm') |
                      startsWith(names(trees), 'ascus_lsil')]


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
  excel.strata.df <- readRDS('params/_excel.strata.df.RData')
  excel.strata.df.full <- readRDS('params/_excel.strata.df.full.RData')
  strat.ctx <- readRDS('params/_strat.ctx.RData')
} else {
  full.strat.ctx <- loadStratifiedContextFile(strat.ctx.path)
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
  saveRDS(excel.strata.df, file='params/_excel.strata.df.RData')
  saveRDS(excel.strata.df.full, file='params/_excel.strata.df.full.RData')
  saveRDS(strat.ctx, file='params/_strat.ctx.RData')
  cat(paste0(md5sum(strat.ctx.path), '\n'), file=strat.ctx.hash.file)
}

calib.strat.ctx.path <- 'params/context_calib.xlsx'
calib.strat.ctx.hash.file <- paste0(calib.strat.ctx.path, '.hash')
if (file.exists(calib.strat.ctx.path)) {
  if (file.exists(calib.strat.ctx.hash.file) && readLines(calib.strat.ctx.hash.file) == md5sum(calib.strat.ctx.path)) {
    calib.full.strat.ctx <- readRDS('params/_calib.full.strat.ctx.RData')
    calib.excel.strata.df <- readRDS('params/_calib.excel.strata.df.RData')
    calib.excel.strata.df.full <- readRDS('params/_calib.excel.strata.df.full.RData')
    calib.strat.ctx <- readRDS('params/_calib.strat.ctx.RData')
  } else {
    calib.full.strat.ctx <- loadStratifiedContextFile(calib.strat.ctx.path)
    calib.excel.strata.df <- list()
    calib.excel.strata.df.full <- list()
    .strata <- names(xlsx::getSheets(xlsx::loadWorkbook(calib.strat.ctx.path)))
    for(stratum in .strata) {
      calib.excel.strata.df[[stratum]] <- read.xlsx(calib.strat.ctx.path, sheetName = stratum, keepFormulas = T)[,c(1,2)]
      calib.excel.strata.df.full[[stratum]] <- read.xlsx(calib.strat.ctx.path, sheetName = stratum, keepFormulas = T)
    }
    calib.strat.ctx <- refresh.context(c(), calib.full.strat.ctx, calib.excel.strata.df, context.setup)
    
    # Save data for faster loading next time if unchanged
    saveRDS(calib.full.strat.ctx, file='params/_calib.full.strat.ctx.RData')
    saveRDS(calib.excel.strata.df, file='params/_calib.excel.strata.df.RData')
    saveRDS(calib.excel.strata.df.full, file='params/_calib.excel.strata.df.full.RData')
    saveRDS(calib.strat.ctx, file='params/_calib.strat.ctx.RData')
    cat(paste0(md5sum(calib.strat.ctx.path), '\n'), file=calib.strat.ctx.hash.file)
  }
} else {
  suppressWarnings({
    rm('calib.full.strat.ctx')
    rm('calib.excel.strata.df')
    rm('calib.excel.strata.df.full')
    rm('calib.strat.ctx')
  })
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
