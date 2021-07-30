source('load_models.R')
source('markov.R')

cat('Running markov model...\n')
markov.outputs <- list()
# Lynch simulation

# initial.state <- sapply(markov$nodes,
#                         function(n) if (n$name=='lynch') 1 else 0)
# markov.result.lynch <- simulate('lynch',
#                                 strategies.lynch,
#                                 markov,
#                                 strat.ctx,
#                                 initial.state,
#                                 start.age=LYNCH.START.AGE,
#                                 max.age=LYNCH.MAX.AGE)
# print(markov.result.lynch$plot)

# PMB simulation

# strat.ctx <- lapply(strat.ctx, function(ctx) {ctx$.p_bleeding <- 1; ctx;})
# strat.ctx <- refresh.context('.p_bleeding', strat.ctx, excel.strata.df, context.setup)

# calib.strat.ctx <- loadStratifiedContextFile('params/context_by_age_adjusted_ly_dogc_calib_mol_spec95.xlsx')
# calib.excel.strata.df <- list()
# .strata <- names(xlsx::getSheets(xlsx::loadWorkbook(strat.ctx.path)))
# for(stratum in .strata) {
#   calib.excel.strata.df[[stratum]] <- read.xlsx(strat.ctx.path, sheetName = stratum, keepFormulas = T)[,c(1,2)]
# }

initial.state <- sapply(markov$nodes,
                        function(n) if (n$name=='postmenopausal_bleeding') 1 else 0)
markov.result.bleeding <- simulate('bleeding',
                                   strategies.bleeding,
                                   markov,
                                   strat.ctx,
                                   initial.state,
                                   start.age=BLEEDING.START.AGE,
                                   max.age=BLEEDING.MAX.AGE,
                                   discount.rate = .03)
print(markov.result.bleeding$plot)

# Asymptomatic simulation

# initial.state <- sapply(markov$nodes, 
#                         function(n) if (n$name=='postmenopausal_asymptomatic') 1 else 0)
# markov.result.asymptomatic <- simulate('asymptomatic',
#                                        strategies.asymptomatic,
#                                        markov,
#                                        strat.ctx,
#                                        initial.state,
#                                        start.age=ASYMPTOMATIC.START.AGE,
#                                        max.age=ASYMPTOMATIC.MAX.AGE)
# print(markov.result.asymptomatic$plot)

cat('Done.\n')