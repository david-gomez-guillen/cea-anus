source('load_models.R')
source('markov.R')

cat('Running markov model...\n')
markov.outputs <- list()

initial.state <- sapply(markov$nodes,
                        function(n) if (n$name=='hiv_positive') 1 else 0)
markov.result <- simulate('hiv_msm',
                           strategies,
                           markov,
                           strat.ctx,
                           initial.state,
                           start.age=40,
                           max.age=80,
                           discount.rate = .03)
print(markov.result$plot)

cat('Done.\n')