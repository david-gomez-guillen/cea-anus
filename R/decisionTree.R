library(CEAModel)
library(dplyr)
library(plotly)

trees <- list()
for(f in list.files('R/models')) {
  if (!startsWith(f, '_') && endsWith(f, '.yaml')) {
    name <- substr(f,1,(nchar(f) - 5))
    # print(name)
    t <- loadDecisionTree(paste0('R/models/',f))
    assign(name, t)
    trees[[name]] <- t
    # t$show(showProbs=F)
  }
}


# template.ctx <- do.call(CEAModel::generateEmptyContextFromModels, unlist(trees))
# CEAModel::saveContextFile(template.ctx, 'R/models/context.xlsx')

strat.ctx <- loadStratifiedContextFile('R/models/context.xlsx')
ctx <- strat.ctx[[1]]

tree.outputs <- list()

results <- lapply(trees, function(tree) tree$analyze(ctx)) %>% bind_rows()
tree.result <- compareStrategies(results)
ggplotly(tree.result$plot)



### Molecular test cost sweep

WTP <- 25000
cohort <- 'lynch'
target <- 'molecular'
reference <- 'current'

sweep.molecular.cost <- lapply(seq(200, 2000, 50), function(cost) {
  ctx.sweep <- ctx
  ctx.sweep$c_molecular_test <- cost
  output <- get(paste0('tree_', cohort, '_', target))$analyze(ctx.sweep)
  output$strategy <- paste0(output$strategy, ' (cost=', cost, ')')
  results.df <- get(paste0('results.', cohort))
  reference.result <- paste0('tree_', cohort, '_', reference)
  icer <- compareStrategies(rbind(results.df[results.df$strategy==reference.result,], output))$summary[2,'ICER']
  data.frame(cost=cost, icer=icer)
}) %>% bind_rows()

p <- ggplot(sweep.molecular.cost, aes(x=cost, y=icer)) + 
  geom_hline(yintercept = WTP, color='red') + 
  geom_line() + 
  geom_point() +
  # scale_y_continuous(breaks=seq(20000, 30000, 2000)) +
  xlab('Molecular test cost ($)') +
  ylab('ICER ($/QALY)') +
  ggtitle(paste0('Molecular test cost sweep on ', cohort, ' ', target, ' (tree)'))

ggplotly(p)


tree_lynch_molecular$show(ctx, nodeInfoFields = 'outcome') %>% visNetwork::visSave(file = 'test')


