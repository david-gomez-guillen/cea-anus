library(CEAModel)
library(ggplot2)
library(reshape2)
library(plotly)
# library(dirmult)

trees <- list()
for(f in list.files('./R/models')) {
  if (!startsWith(f, '_') && endsWith(f, '.yaml')) {
    name <- substr(f,1,(nchar(f) - 5))
    print(name)
    t <- loadDecisionTree(paste0('R/models/',f))
    assign(name, t)
    trees[[name]] <- t
    # t$show(showProbs=F)
  }
}

markov <- loadMarkovModels('R/models/markov.xlsx')

ctx <- loadContextFile('R/models/context.xlsx')



N.YEARS <- 40

results.df <- data.frame()
for(tree in trees) {
  outcomes <- tree$calculateOutcomes(ctx)
  ctx$p_cancer <- outcomes[outcomes$name=='surgery','prob'] * ctx$p_surgery_cancer
  ctx$c_healthy <- weighted.mean(outcomes$cost, outcomes$prob)
  costs <- c(ctx$c_healthy, ctx$c_cancer, ctx$c_survive, 0, 0)
  utilities <- c(ctx$u_healthy, ctx$u_cancer, ctx$u_survive, 0, 0)
  tpMatrix <- markov$evaluateTpMatrix(ctx)
  
  overall.cost <- 0
  overall.qalys <- 0
  state <- c(1, 0, 0, 0, 0)
  names(state) <- sapply(markov$nodes, function(n) n$name)
  states <- data.frame(t(state))
  for(year in seq(N.YEARS)) {
    state <- state %*% tpMatrix
    overall.cost <- overall.cost + sum(state * costs)
    overall.qalys <- overall.qalys + sum(state * utilities)
    
    states <- rbind(states, state)
  }
  # print(states)
  # print(overall.cost)
  # print(overall.qalys)
  states$iteration <- row.names(states)
  melted.states <- melt(states, id.vars='iteration')
  
  p <- ggplot(melted.states, aes(x=as.numeric(iteration), y=value, color=variable)) + 
    geom_line() + 
    ylim(0,1) + 
    xlab('Year') + 
    ylab('Cohort (%)')
  print(ggplotly(p, main='a'))
  
  results.df <- rbind(results.df, data.frame(strategy=tree$name,
                                       C=overall.cost,
                                       E=overall.qalys, stringsAsFactors = FALSE))
}

r <- CEAModel::compareStrategies(results.df)

print(ggplotly(r$plot))
# CEAsummary<- CEAModel::compareStrategies(
#                                 conventional=conventional,
#                                 hpv16=conventional_hpv16la,
#                                 hpv1618=conventional_hpv1618la,
#                                 hpvhrla=conventional_hpvhrla,
#                                 hpvhrhc=conventional_hpvhrhc,
#                                 arnm_16=arnme6e7_hpv16,
#                                 arnm_1618=arnme6e7_hpv1618,
#                                 arnm_hr=arnme6e7_hpvhr,
#                                context=generateContext(N))
# 
# print(CEAsummary$plot)
