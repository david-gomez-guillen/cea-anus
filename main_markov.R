library(dplyr)
#setwd('~/Documents/models_ce/anus')
source('load_models.R')
source('markov.R')

cat('Running markov model...\n')
markov.outputs <- list()

initial.state <- sapply(markov$nodes,
                        function(n) if (n$name=='hiv_positive') 1 else 0)
markov.result <- simulate('hiv_msm',
                           strategies$hiv_msm,
                           markov,
                           strat.ctx,
                           initial.state,
                           discount.rate=.03)

# Remove redundant suffix for strategy names
markov.result$summary$strategy <- gsub('(.*?)-(.*)', '\\1', markov.result$summary$strategy)

# print(markov.result$plot)
# print(markov.result$summary)
# print(ggplotly(markov.result$plot + theme(legend.position = 'none')))

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


cat('Disc0\n')
markov.result.0 <- simulate('hiv_msm',
                          strategies$hiv_msm,
                          markov,
                          strat.ctx,
                          initial.state,
                          discount.rate = .0)

strat.ctx <- lapply(strat.ctx, function(ctx) {
  ctx[c('c_arnm16','c_arnm161845','c_arnmhr')] <- 6
  ctx
})

markov.result.3.6e <- simulate('hiv_msm',
                          strategies$hiv_msm,
                          markov,
                          strat.ctx,
                          initial.state,
                          discount.rate = .03)

markov.result.0.6e <- simulate('hiv_msm',
                               strategies$hiv_msm,
                               markov,
                               strat.ctx,
                               initial.state,
                               discount.rate = .0)
source('load_models.R')

summary.df <- lapply(names(markov.result$info), function(strat) {
  info <- markov.result$info[[strat]]$additional.info
  info.sum <- apply(info, 2, sum)
  info.sum$strategy <- strat
  info.sum$eff <- info.sum$eff / 2

  info.0 <- markov.result.0$info[[strat]]$additional.info
  info.0.6e <- markov.result.0.6e$info[[strat]]$additional.info
  info.3.6e <- markov.result.3.6e$info[[strat]]$additional.info

  info.sum$eff.0 <- sum(info.0$eff)/2
  info.sum$cost.0.25e <- sum(info.0$cost)
  info.sum$cost.0.6e <- sum(info.0.6e$cost)
  info.sum$cost.3.6e <- sum(info.3.6e$cost)
  info.sum
}) %>% bind_rows()

summary.df <- summary.df[,c(19, 3, 21, 23, 22, 4, 20, 5:18)]
summary.df[,c(8:21)] <- summary.df[,c(8:21)] * 100000
names(summary.df)[2:7] <- c('cost.disc3.test25', 'cost.disc0.test25', 'cost.disc3.test6', 'cost.disc0.test6', 'qalys.disc3', 'qalys.disc0')
# write.xlsx(summary.df, 'output/summary_anus.xlsx', row.names = F)
