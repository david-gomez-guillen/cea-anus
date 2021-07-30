library(CEAModel)
library(ggplot2)
library(reshape2)
library(plotly)
library(dplyr)

EPSILON <- 1e-7  # Maximum allowed error due to precision errors
DISCOUNT.RATE <- .03
PROPHYLACTIC.SURGERY.AGE <- 40  # 35-40
MENOPAUSE.AGE <- 50

LYNCH.START.AGE <- 25
LYNCH.MAX.AGE <- 70
ASYMPTOMATIC.START.AGE <- 50
ASYMPTOMATIC.MAX.AGE <- 85
BLEEDING.START.AGE <- 50
BLEEDING.MAX.AGE <- 85
COMBINED.START.AGE <- 25
COMBINED.MAX.AGE <- 80

DEFAULT.START.AGE <- list(
  asymptomatic=50,
  bleeding=50,
  lynch=25
)
DEFAULT.MAX.AGE <- list(
  asymptomatic=85,
  bleeding=85,
  lynch=60
)

POPULATION.REFERENCES <- list(
  asymptomatic='tree_asymptomatic_current',
  bleeding='tree_bleeding_current2',
  lynch='tree_lynch_current2'
)

POPULATION.INITIAL.STATES <- list(
  asymptomatic='postmenopausal_asymptomatic',
  bleeding='postmenopausal_bleeding',
  lynch='lynch'
)

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
    context.full <- strat.ctx[[stratum]]
    context <- lapply(context.full, function(e) e[1]) # Base values
    # Assign markov probabilities according to tree outcomes
    if ('lynch' %in% names(trees)) {
      prevalence <- strat.ctx[[stratum]]$.p_cancer___lynch
      outcomes <- trees$lynch$summarize(context, prevalence=prevalence)
      context$._p_hysterectomy_unnecessary_lynch <- outcomes[outcomes$name == 'no_cancer_hysterectomy', 'prob']
      context$._p_hysterectomy_necessary_lynch <- outcomes[outcomes$name == 'cancer_removed', 'prob']
      context$p_hysterectomy_lynch <- context$._p_hysterectomy_unnecessary_lynch + context$._p_hysterectomy_necessary_lynch
      
      cancer.states <- c('cancer', 'undetected_cancer_discharge', 'undetected_cancer_evaluation')
      context$p_cancer_lynch <- sum(outcomes[outcomes$name %in% cancer.states, 'prob'])
      context$p_cancer_lynch_s1 <- sum(outcomes[outcomes$name %in% 'cancer_removed', 'prob'])
      context$p_cancer_lynch_s2 <- sum(outcomes[outcomes$name %in% 'cancer', 'prob']) * context$.p_cancer_stage_2 / p_cancer_s2_4
      context$p_cancer_lynch_s3 <- sum(outcomes[outcomes$name %in% 'cancer', 'prob']) * context$.p_cancer_stage_3 / p_cancer_s2_4
      context$p_cancer_lynch_s4 <- sum(outcomes[outcomes$name %in% 'cancer', 'prob']) * (1 - (context$.p_cancer_stage_2 +
                                                                                                context$.p_cancer_stage_3) / p_cancer_s2_4)
      
      context$p_discharge <- 0 # Only discharged for three years in asymptomatic
      context$p_death_cancer_lynch <- 0
      lynch.cost <- weighted.mean(outcomes$cost, outcomes$prob)
    } else {
      context$p_hysterectomy_lynch <- 0
      context$p_cancer_lynch <- 0
      context$p_cancer_lynch_s1 <- 0
      context$p_cancer_lynch_s2 <- 0
      context$p_cancer_lynch_s3 <- 0
      context$p_cancer_lynch_s4 <- 0
      context$p_death_cancer_lynch <- 0
      context$p_evaluation_bleeding_endo_thin <- 0
      context$p_cancer_bleeding_endo_thin <- 0
      context$p_hysterectomy_bleeding_endo_thin <- 0
      lynch.cost <- 0
    } 
    if ('asymptomatic' %in% names(trees)) {
      prevalence <- strat.ctx[[stratum]]$.p_cancer___asymptomatic
      outcomes <- trees$asymptomatic$summarize(context, prevalence=prevalence)
      strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_asymptomatic <- outcomes[outcomes$name == 'no_cancer_hysterectomy', 'prob']
      strat.ctx[[stratum]]$._p_hysterectomy_necessary_asymptomatic <- outcomes[outcomes$name == 'cancer_removed', 'prob']
      context$p_hysterectomy_asymptomatic <- strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_asymptomatic + strat.ctx[[stratum]]$._p_hysterectomy_necessary_asymptomatic
      
      cancer.states <- c('cancer', 'undetected_cancer_discharge', 'undetected_cancer_evaluation')
      context$p_cancer_asymptomatic <- sum(outcomes[outcomes$name %in% cancer.states, 'prob'])
      context$p_cancer_asymptomatic_s1 <- sum(outcomes[outcomes$name %in% 'cancer_removed', 'prob'])
      context$p_cancer_asymptomatic_s2 <- sum(outcomes[outcomes$name %in% 'cancer', 'prob']) * context$.p_cancer_stage_2 / p_cancer_s2_4
      context$p_cancer_asymptomatic_s3 <- sum(outcomes[outcomes$name %in% 'cancer', 'prob']) * context$.p_cancer_stage_3 / p_cancer_s2_4
      context$p_cancer_asymptomatic_s4 <- sum(outcomes[outcomes$name %in% 'cancer', 'prob']) * (1 - (context$.p_cancer_stage_2 +
                                                                                                       context$.p_cancer_stage_3) / p_cancer_s2_4)
      
      discharge.states <- c('no_cancer_discharge')
      context$p_discharge <- sum(outcomes[outcomes$name %in% discharge.states, 'prob'])
      context$p_death_cancer_asymptomatic <- 0
      asymptomatic.cost <- weighted.mean(outcomes$cost, outcomes$prob)
    } else {
      context$p_hysterectomy_asymptomatic <- 0
      context$p_cancer_asymptomatic <- 0
      context$p_cancer_asymptomatic_s1 <- 0
      context$p_cancer_asymptomatic_s2 <- 0
      context$p_cancer_asymptomatic_s3 <- 0
      context$p_cancer_asymptomatic_s4 <- 0
      context$p_death_cancer_asymptomatic <- 0
      context$p_evaluation_bleeding_endo_thin <- 0
      context$p_cancer_bleeding_endo_thin <- 0
      context$p_hysterectomy_bleeding_endo_thin <- 0
      asymptomatic.cost <- 0
    }
    if ('bleeding' %in% names(trees)) {
      prevalence <- strat.ctx[[stratum]]$.p_cancer___bleeding
      outcomes <- trees$bleeding$summarize(context, prevalence=prevalence)
      outcomes.undetected <- trees$bleeding$summarize(context, prevalence=1) # Outcomes for the previously undetected cancers (so prevalence=100%)
      # TODO: Replicate for other populations
      strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_bleeding <- outcomes[outcomes$name == 'no_cancer_hysterectomy', 'prob']
      if (length(strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_bleeding) == 0)
        strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_bleeding <- 0
      strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding <- sum(outcomes[outcomes$name %in% c('cancer_removed', 'cancer_s1', 'cancer_s2', 'cancer_s3', 'cancer_s4'), 'prob'])
      strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% c('cancer_removed', 'cancer_s1', 'cancer_s2', 'cancer_s3', 'cancer_s4'), 'prob'])
      if (length(strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding) == 0)
        strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding <- 0
      context$p_hysterectomy_bleeding <- strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_bleeding
      
      cancer.states <- c('cancer', 
                         'cancer_s1',
                         'cancer_s2',
                         'cancer_s3',
                         'cancer_s4'
                         # 'undetected_cancer_discharge', 
                         # 'undetected_cancer_evaluation', 
                         # 'undetected_cancer_evaluation_bleeding_endo_thin'
                         )
      context$p_cancer_bleeding <- sum(outcomes[outcomes$name %in% cancer.states, 'prob'])
      if (length(context$p_cancer_bleeding) == 0) context$p_cancer_bleeding <- 0
      context$p_cancer_bleeding_s1 <- sum(outcomes[outcomes$name %in% 'cancer_s1', 'prob'])
      if (length(context$p_cancer_bleeding_s1) == 0) context$p_cancer_bleeding_s1 <- 0
      context$p_cancer_bleeding_s2 <- sum(outcomes[outcomes$name %in% 'cancer_s2', 'prob'])
      if (length(context$p_cancer_bleeding_s2) == 0) context$p_cancer_bleeding_s2 <- 0
      context$p_cancer_bleeding_s3 <- sum(outcomes[outcomes$name %in% 'cancer_s3', 'prob'])
      if (length(context$p_cancer_bleeding_s3) == 0) context$p_cancer_bleeding_s3 <- 0
      context$p_cancer_bleeding_s4 <- sum(outcomes[outcomes$name %in% 'cancer_s4', 'prob'])
      if (length(context$p_cancer_bleeding_s4) == 0) context$p_cancer_bleeding_s4 <- 0
      
      # Calculate probabilities for cancer detection in previously undetected cancers
      context$p_cancer_bleeding_s1_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s1', 'prob'])
      if (length(context$p_cancer_bleeding_s1_undetected) == 0) context$p_cancer_bleeding_s1_undetected <- 0
      context$p_cancer_bleeding_s2_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s2', 'prob'])
      if (length(context$p_cancer_bleeding_s2_undetected) == 0) context$p_cancer_bleeding_s2_undetected <- 0
      context$p_cancer_bleeding_s3_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s3', 'prob'])
      if (length(context$p_cancer_bleeding_s3_undetected) == 0) context$p_cancer_bleeding_s3_undetected <- 0
      context$p_cancer_bleeding_s4_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s4', 'prob'])
      if (length(context$p_cancer_bleeding_s4_undetected) == 0) context$p_cancer_bleeding_s4_undetected <- 0
      context$p_cancer_bleeding_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% cancer.states, 'prob'])
      if (length(context$p_cancer_bleeding_undetected) == 0) context$p_cancer_bleeding_undetected <- 0
      
      context$p_undetected_cancer_bleeding <- sum(outcomes[outcomes$name %in% c('undetected_cancer_discharge', 
                                                                                'undetected_cancer_evaluation', 
                                                                                'undetected_cancer_evaluation_bleeding_endo_thin'), 'prob'])
      if (length(context$p_undetected_cancer_bleeding) == 0) context$p_undetected_cancer_bleeding <- 0
      
      context$p_discharge <- 0 # Only discharged for three years in asymptomatic
      context$p_evaluation_bleeding_endo_thin <- outcomes[outcomes$name=='no_cancer_evaluation_bleeding_endo_thin', 'prob']
      if (length(context$p_evaluation_bleeding_endo_thin) == 0) context$p_evaluation_bleeding_endo_thin <- 0
      context$p_evaluation_bleeding_endo_thin_undetected <- outcomes.undetected[outcomes.undetected$name=='no_cancer_evaluation_bleeding_endo_thin', 'prob']
      if (length(context$p_evaluation_bleeding_endo_thin_undetected) == 0) context$p_evaluation_bleeding_endo_thin_undetected <- 0
      context$p_death_cancer_bleeding <- 0
      bleeding.cost <- weighted.mean(outcomes$cost, outcomes$prob)
      
      if (trees$bleeding$name == 'tree_bleeding_current') {
        # TODO: Generalize
        tree_endo_thin <- tree_bleeding_current_endo_thin
        tn <- outcomes[outcomes$name == 'no_cancer_evaluation_bleeding_endo_thin',]
        fp <- outcomes[outcomes$name == 'undetected_cancer_evaluation_bleeding_endo_thin',]
        prevalence <- fp / (fp + tn)
        outcomes <- tree_endo_thin$summarize(context, prevalence=prevalence)
        outcomes.undetected <- tree_endo_thin$summarize(context, prevalence=1) # Outcomes for the previously undetected cancers (so prevalence=100%)
        strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_bleeding_endo_thin <- outcomes[outcomes$name == 'no_cancer_hysterectomy', 'prob']
        strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding_endo_thin <- sum(outcomes[outcomes$name %in% c('cancer_removed', 'cancer_s1', 'cancer_s2', 'cancer_s3', 'cancer_s4'), 'prob'])
        strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding_endo_thin_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% c('cancer_removed', 'cancer_s1', 'cancer_s2', 'cancer_s3', 'cancer_s4'), 'prob'])
        context$p_hysterectomy_bleeding_endo_thin <- strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_bleeding_endo_thin #+ strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding_endo_thin
        
        # cancer.states <- c('cancer', 'undetected_cancer_discharge', 'undetected_cancer_evaluation')
        context$p_cancer_bleeding_endo_thin <- sum(outcomes[outcomes$name %in% cancer.states, 'prob'])
        context$p_cancer_bleeding_endo_thin_s1 <- sum(outcomes[outcomes$name %in% 'cancer_s1', 'prob'])
        context$p_cancer_bleeding_endo_thin_s2 <- sum(outcomes[outcomes$name %in% 'cancer_s2', 'prob'])
        context$p_cancer_bleeding_endo_thin_s3 <- sum(outcomes[outcomes$name %in% 'cancer_s3', 'prob'])
        context$p_cancer_bleeding_endo_thin_s4 <- sum(outcomes[outcomes$name %in% 'cancer_s4', 'prob'])
        context$p_undetected_cancer_bleeding_endo_thin <- sum(outcomes[outcomes$name %in% c('undetected_cancer_discharge', 
                                                                                            'undetected_cancer_evaluation', 
                                                                                            'undetected_cancer_evaluation_bleeding_endo_thin'), 'prob'])
        context$p_cancer_bleeding_endo_thin_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% cancer.states, 'prob'])
        context$p_cancer_bleeding_endo_thin_s1_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s1', 'prob'])
        context$p_cancer_bleeding_endo_thin_s2_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s2', 'prob'])
        context$p_cancer_bleeding_endo_thin_s3_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s3', 'prob'])
        context$p_cancer_bleeding_endo_thin_s4_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s4', 'prob'])
        
        context$p_death_cancer_bleeding <- 0
        bleeding.cost <- bleeding.cost + weighted.mean(outcomes$cost, outcomes$prob)
      } else if (trees$bleeding$name == 'tree_bleeding_current2') {
        # TODO: Generalize
        tree_endo_thin2 <- tree_bleeding_current_endo_thin2
        tn <- outcomes[outcomes$name == 'no_cancer_evaluation_bleeding_endo_thin','prob']
        fp <- outcomes[outcomes$name == 'undetected_cancer_evaluation_bleeding_endo_thin','prob']
        prevalence <- fp / (fp + tn)
        outcomes <- tree_endo_thin2$summarize(context, prevalence=prevalence)
        outcomes.undetected <- tree_endo_thin2$summarize(context, prevalence=1) # Outcomes for the previously undetected cancers (so prevalence=100%)
        strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_bleeding_endo_thin <- outcomes[outcomes$name == 'no_cancer_hysterectomy', 'prob']
        strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding_endo_thin <- sum(outcomes[outcomes$name %in% c('cancer_removed', 'cancer_s1', 'cancer_s2', 'cancer_s3', 'cancer_s4'), 'prob'])
        strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding_endo_thin_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% c('cancer_removed', 'cancer_s1', 'cancer_s2', 'cancer_s3', 'cancer_s4'), 'prob'])
        context$p_hysterectomy_bleeding_endo_thin <- strat.ctx[[stratum]]$._p_hysterectomy_unnecessary_bleeding_endo_thin #+ strat.ctx[[stratum]]$._p_hysterectomy_necessary_bleeding_endo_thin
        
        # cancer.states <- c('cancer', 'undetected_cancer_discharge', 'undetected_cancer_evaluation')
        context$p_cancer_bleeding_endo_thin <- sum(outcomes[outcomes$name %in% cancer.states, 'prob'])
        context$p_cancer_bleeding_endo_thin_s1 <- sum(outcomes[outcomes$name %in% 'cancer_s1', 'prob'])
        context$p_cancer_bleeding_endo_thin_s2 <- sum(outcomes[outcomes$name %in% 'cancer_s2', 'prob'])
        context$p_cancer_bleeding_endo_thin_s3 <- sum(outcomes[outcomes$name %in% 'cancer_s3', 'prob'])
        context$p_cancer_bleeding_endo_thin_s4 <- sum(outcomes[outcomes$name %in% 'cancer_s4', 'prob'])
        context$p_undetected_cancer_bleeding_endo_thin <- sum(outcomes[outcomes$name %in% c('undetected_cancer_discharge', 
                                                                                            'undetected_cancer_evaluation', 
                                                                                            'undetected_cancer_evaluation_bleeding_endo_thin'), 'prob'])
        context$p_cancer_bleeding_endo_thin_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% cancer.states, 'prob'])
        context$p_cancer_bleeding_endo_thin_s1_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s1', 'prob'])
        context$p_cancer_bleeding_endo_thin_s2_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s2', 'prob'])
        context$p_cancer_bleeding_endo_thin_s3_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s3', 'prob'])
        context$p_cancer_bleeding_endo_thin_s4_undetected <- sum(outcomes.undetected[outcomes.undetected$name %in% 'cancer_s4', 'prob'])
        
        context$p_death_cancer_bleeding <- 0
        bleeding.cost <- bleeding.cost + weighted.mean(outcomes$cost, outcomes$prob)
      } else{
        context$p_cancer_bleeding_endo_thin <- 0
        context$p_cancer_bleeding_endo_thin_s1 <- 0
        context$p_cancer_bleeding_endo_thin_s2 <- 0
        context$p_cancer_bleeding_endo_thin_s3 <- 0
        context$p_cancer_bleeding_endo_thin_s4 <- 0
        context$p_cancer_bleeding_endo_thin_undetected <- 0
        context$p_cancer_bleeding_endo_thin_s1_undetected <- 0
        context$p_cancer_bleeding_endo_thin_s2_undetected <- 0
        context$p_cancer_bleeding_endo_thin_s3_undetected <- 0
        context$p_cancer_bleeding_endo_thin_s4_undetected <- 0
        context$p_undetected_cancer_bleeding_endo_thin <- 0
        context$p_hysterectomy_bleeding_endo_thin <- 0
        context$._p_hysterectomy_unnecessary_bleeding_endo_thin <- 0
        context$._p_hysterectomy_necessary_bleeding_endo_thin <- 0
        context$._p_hysterectomy_necessary_bleeding_endo_thin_undetected <- 0
      }
    } else {
      context$p_hysterectomy_bleeding <- 0
      context$._p_hysterectomy_necessary_bleeding_undetected <- 0
      context$p_hysterectomy_bleeding_endo_thin <- 0
      context$._p_hysterectomy_necessary_bleeding_endo_thin_undetected <- 0
      context$p_cancer_bleeding <- 0
      context$p_cancer_bleeding_s1 <- 0
      context$p_cancer_bleeding_s2 <- 0
      context$p_cancer_bleeding_s3 <- 0
      context$p_cancer_bleeding_s4 <- 0
      context$p_cancer_bleeding_undetected <- 0
      context$p_cancer_bleeding_s1_undetected <- 0
      context$p_cancer_bleeding_s2_undetected <- 0
      context$p_cancer_bleeding_s3_undetected <- 0
      context$p_cancer_bleeding_s4_undetected <- 0
      context$p_cancer_bleeding_endo_thin <- 0
      context$p_cancer_bleeding_endo_thin_s1 <- 0
      context$p_cancer_bleeding_endo_thin_s2 <- 0
      context$p_cancer_bleeding_endo_thin_s3 <- 0
      context$p_cancer_bleeding_endo_thin_s4 <- 0
      context$p_cancer_bleeding_endo_thin_undetected <- 0
      context$p_cancer_bleeding_endo_thin_s1_undetected <- 0
      context$p_cancer_bleeding_endo_thin_s2_undetected <- 0
      context$p_cancer_bleeding_endo_thin_s3_undetected <- 0
      context$p_cancer_bleeding_endo_thin_s4_undetected <- 0
      context$p_undetected_cancer_bleeding_endo_thin <- 0
      context$p_death_cancer_bleeding <- 0
      context$p_cancer_bleeding_endo_thin <- 0
      bleeding.cost <- 0
    }
    
    ####### TEMPORARY
    context$p_bleeding <- 0
    context$p_stop_bleeding <- 0
    #######
    extended.strat.ctx[[stratum]] <- context
    tpMatrices[[stratum]] <- markov$evaluateTpMatrix(context)
    costs[[stratum]]['postmenopausal_asymptomatic'] <- asymptomatic.cost
    costs[[stratum]]['postmenopausal_bleeding'] <- bleeding.cost
    costs[[stratum]]['lynch'] <- lynch.cost
    
    # if (trees$bleeding$name == 'tree_bleeding_hysterectomy') {
    #   row <- tpMatrices[[stratum]]$strategies['postmenopausal_bleeding',]
    #   excess.prob <- sum(row[!names(row) %in% 'postmenopausal_bleeding']) - 1
    #   tpMatrices[[stratum]]$strategies['postmenopausal_bleeding','postmenopausal_bleeding'] <- 0
    #   tpMatrices[[stratum]]$strategies['postmenopausal_bleeding','postmenopausal_hysterectomy'] <- tpMatrices[[stratum]]$strategies['postmenopausal_bleeding','postmenopausal_hysterectomy'] - excess.prob
    # }
    
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

simulate.markov <- function(trees, 
                            markov,
                            initial.state, 
                            strat.ctx, 
                            start.age=25, 
                            max.age=60,
                            prophylactic.surgery.age=PROPHYLACTIC.SURGERY.AGE,
                            menopause.age=MENOPAUSE.AGE,
                            discount.rate=DISCOUNT.RATE) {
  additional.info <- data.frame(year=start.age-1, 
                                cost=0,
                                eff=0,
                                unnecessary_hysterectomies=0, 
                                necessary_hysterectomies=0, 
                                total_hysterectomies=0, 
                                incidence=0,
                                incidence_s1=0,
                                incidence_s2=0,
                                incidence_s3=0,
                                incidence_s4=0,
                                mortality=0,
                                mortality_s1=0,
                                mortality_s2=0,
                                mortality_s3=0,
                                mortality_s4=0)
  
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
    states <- rbind(states, current.state)
    tpMatrix <- get.matrix.stratum(tpMatrices, year)
    ctx <- get.context.stratum(strat.ctx, year)
    
    year.costs <- get.context.stratum(costs, year)
    year.utilities <- get.context.stratum(utilities, year)
    if (year == start.age) {
      year.costs['postmenopausal_bleeding'] <- year.costs['postmenopausal_bleeding'] - ctx$.c_visit + ctx$.c_first_visit
    }
    
    if (year == menopause.age) {
      tpMatrix$other['premenopausal_asymptomatic', 'postmenopausal_asymptomatic'] <- tpMatrix$other['premenopausal_asymptomatic', 'premenopausal_asymptomatic']
      tpMatrix$other['premenopausal_asymptomatic', 'premenopausal_asymptomatic'] <- 0
    }
    current.cost <- sum(current.state * year.costs) * (1-discount.rate)^(year-start.age)
    overall.cost <- overall.cost + current.cost
    current.eff <- sum(current.state * year.utilities) * (1-discount.rate)^(year-start.age)
    overall.eff <- overall.eff + current.eff
    
# TODO: Generalize
custom.data.extractor <- function(additional.info, current.state, tpMatrix, cost, eff, ctx) {
  unnecessary_hysterectomies <- current.state[,'postmenopausal_bleeding'][[1]] * ctx$p_hysterectomy_bleeding +
             current.state[,'postmenopausal_bleeding_endo_thin'][[1]] * ctx$p_hysterectomy_bleeding_endo_thin
  necessary_hysterectomies <- current.state[,'postmenopausal_bleeding'][[1]] * ctx$p_cancer_bleeding +
      current.state[,'postmenopausal_bleeding_endo_thin'][[1]] * ctx$p_cancer_bleeding_endo_thin +
      current.state[,'undetected_cancer_bleeding'][[1]] * ctx$p_cancer_bleeding_undetected +
      current.state[,'undetected_cancer_bleeding_endo_thin'][[1]] * ctx$p_cancer_bleeding_endo_thin_undetected
  total_hysterectomies <- unnecessary_hysterectomies + necessary_hysterectomies
  
  if ('endometrial_cancer' %in% colnames(tpMatrix$strategies)) {
    incidence <- current.state[, 'postmenopausal_bleeding'] * tpMatrix$strategies['postmenopausal_bleeding', 'endometrial_cancer'] +
      current.state[, 'postmenopausal_bleeding_endo_thin'] * tpMatrix$strategies['postmenopausal_bleeding_endo_thin', 'endometrial_cancer']
    incidence_s1 <- NULL
    incidence_s2 <- NULL
    incidence_s3 <- NULL
    incidence_s4 <- NULL
    
    mortality <- current.state[, 'endometrial_cancer'] * tpMatrix$strategies['endometrial_cancer', 'death_cancer']
    mortality_s1 <- NULL
    mortality_s2 <- NULL
    mortality_s3 <- NULL
    mortality_s4 <- NULL
  } else {
    incidence_s1 <- current.state[, 'postmenopausal_bleeding'][[1]] * sum(tpMatrix$strategies['postmenopausal_bleeding', 'endometrial_cancer_s1']) +
      current.state[, 'postmenopausal_bleeding_endo_thin'][[1]] * sum(tpMatrix$strategies['postmenopausal_bleeding_endo_thin', 'endometrial_cancer_s1'])
    incidence_s2 <- current.state[, 'postmenopausal_bleeding'][[1]] * sum(tpMatrix$strategies['postmenopausal_bleeding', 'endometrial_cancer_s2']) +
      current.state[, 'postmenopausal_bleeding_endo_thin'][[1]] * sum(tpMatrix$strategies['postmenopausal_bleeding_endo_thin', 'endometrial_cancer_s2'])
    incidence_s3 <- current.state[, 'postmenopausal_bleeding'][[1]] * sum(tpMatrix$strategies['postmenopausal_bleeding', 'endometrial_cancer_s3']) +
      current.state[, 'postmenopausal_bleeding_endo_thin'][[1]] * sum(tpMatrix$strategies['postmenopausal_bleeding_endo_thin', 'endometrial_cancer_s3'])
    incidence_s4 <- current.state[, 'postmenopausal_bleeding'][[1]] * sum(tpMatrix$strategies['postmenopausal_bleeding', 'endometrial_cancer_s4']) +
      current.state[, 'postmenopausal_bleeding_endo_thin'][[1]] * sum(tpMatrix$strategies['postmenopausal_bleeding_endo_thin', 'endometrial_cancer_s4'])
    incidence <- incidence_s1 + incidence_s2 + incidence_s3 + incidence_s4
    
    mortality_s1 <- current.state[, 'endometrial_cancer_s1'][[1]] * tpMatrix$other['endometrial_cancer_s1', 'death_cancer']
    mortality_s2 <- current.state[, 'endometrial_cancer_s2'][[1]] * tpMatrix$other['endometrial_cancer_s2', 'death_cancer']
    mortality_s3 <- current.state[, 'endometrial_cancer_s3'][[1]] * tpMatrix$other['endometrial_cancer_s3', 'death_cancer']
    mortality_s4 <- current.state[, 'endometrial_cancer_s4'][[1]] * tpMatrix$other['endometrial_cancer_s4', 'death_cancer']
    mortality <- mortality_s1 + mortality_s2 + mortality_s3 + mortality_s4
  }
  additional.info <- rbind(additional.info,
                           list(
                             year=year,
                             cost=cost,
                             eff=eff,
                             unnecessary_hysterectomies=unnecessary_hysterectomies,
                             necessary_hysterectomies=necessary_hysterectomies,
                             total_hysterectomies=total_hysterectomies,
                             incidence=incidence,
                             incidence_s1=incidence_s1,
                             incidence_s2=incidence_s2,
                             incidence_s3=incidence_s3,
                             incidence_s4=incidence_s4,
                             mortality=mortality,
                             mortality_s1=mortality_s1,
                             mortality_s2=mortality_s2,
                             mortality_s3=mortality_s3,
                             mortality_s4=mortality_s4))
  return(additional.info)
}

additional.info <- custom.data.extractor(additional.info, current.state, tpMatrix, current.cost, current.eff, ctx)
    next.state <- current.state %*% tpMatrix$time_pass %*% tpMatrix$other %*% tpMatrix$strategies

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
  states$age <- as.numeric(row.names(states)) + start.age - 1
  states$postmenopausal_asymptomatic_healthy <- states$postmenopausal_asymptomatic + states$postmenopausal_asymptomatic_1y + states$postmenopausal_asymptomatic_2y
  states$postmenopausal_asymptomatic <- NULL
  states$postmenopausal_asymptomatic_1y <- NULL
  states$postmenopausal_asymptomatic_2y <- NULL
  states$survive_all <- states$survive_s1 + states$survive_s2 + states$survive_s3 + states$survive_s4
  states$endometrial_cancer_all <- states$endometrial_cancer_s1 + states$endometrial_cancer_s2 + states$endometrial_cancer_s3 + states$endometrial_cancer_s4
  states$death_all <- states$death_cancer + states$death_other
  states$postmenopausal_bleeding_all <- states$postmenopausal_bleeding + states$postmenopausal_bleeding_endo_thin
  states$postmenopausal_bleeding <- NULL
  states$postmenopausal_bleeding_endo_thin <- NULL
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
                           E=overall.eff, stringsAsFactors = FALSE)
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
    result <- simulate.markov(trees, markov, initial.state, strat.ctx,
                              start.age=start.age, max.age=max.age,
                              discount.rate=discount.rate)
    markov.outputs[[tree$name]] <- result
    additional.info[[tree$name]] <- result$additional.info
    results.df <- rbind(results.df, result$summary)
  }
  
  r <- CEAModel::analyzeCE(results.df, cost.label='€', eff.label='QALY')
  r$info <- markov.outputs
  return(r)
}
