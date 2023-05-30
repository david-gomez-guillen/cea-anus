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
print(markov.result$summary)

# x <- markov.result$info$conventional$additional.info
# xt <- markov.result$info$conventional_t$additional.info
# df <- data.frame(iter=x$iter,
#                  strategy='conventional',
#                  ac_incidence=x$ac_incidence,
#                  ac_mortality=x$ac_mortality)
# df <- rbind(df, data.frame(iter=xt$iter,
#                            strategy='conventional_t',
#                            ac_incidence=xt$ac_incidence,
#                            ac_mortality=xt$ac_mortality))
# df <- melt(df, id.vars = c('iter', 'strategy'))
# plt <- ggplot(df, aes(x=iter, y=value, color=strategy, linetype=variable)) + 
#   geom_line()
# 
# print(plt)

cat('Done.\n')