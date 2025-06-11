library(dplyr)
library(Cairo)
library(ggplot2)
library(gridExtra)
library(scales)
library(data.table)
library(officer)

SUFFIX <- ''
N.CORES.DEFAULT <- 1
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  N.CORES <- N.CORES.DEFAULT
} else {
  N.CORES <- as.numeric(args[1])
}
sd.estimate.func <- function(p.name, value) value/5

source('load_models.R')

# source('main_markov.R')

# Molecular cost PSAs

source('markov_psa.R')
library(GenBinomApps)
cost.df <- data.frame()
THRESHOLD.WTP <- 25000
N <- 1000

evaluated.trees <- c(
  'arnme6e7_hpvhr_t_tca'
)
cost.vals <- seq(0, 26, 2)

for(tree in evaluated.trees) {
  for(i in cost.vals) {
    sweep.strat.ctx <- lapply(strat.ctx, function(ctx) {
      ctx$c_arn_kit <- i
      return(ctx)
    })
    results <- psa.1('c_arn_kit', 
                     sweep.strat.ctx,
                     excel.file = 'params/context.xlsx',
                     population='hiv_msm', 
                     strategy=tree,
                     markov=markov,
                     reference='conventional_t_tca',
                     sd.estimate.func=sd.estimate.func,
                     context.setup.func=context.setup,
                     n.cores=N.CORES,
                     n.iters=N,
                     seed=1234,
                     discount.rate=.03)
    
    df <- results$c_arn_kit$summary
    icers <- df$ICER
    xs <- seq(0, 50000, 500)
    # ys <- sapply(xs, function(wtp) mean(icers<wtp))
    ys <- sapply(xs, function(wtp) mean(df$IE-df$IC/wtp>0))
    ys.ci <- lapply(xs, function(wtp)
      clopper.pearson.ci(sum(df$IE-df$IC/wtp>0), N, alpha = .05, CI='two.sided')
      )
    ys.l <- sapply(ys.ci, function(ci) ci$Lower.limit)
    ys.h <- sapply(ys.ci, function(ci) ci$Upper.limit)
    # pct.wtp <- mean(icers < THRESHOLD.WTP)
    pct.wtp <- mean(df$IE - df$IC/THRESHOLD.WTP > 0)
    cost.df <- rbind(cost.df, data.frame(strategy=tree,
                                         cost=i,
                                         wtp=xs,
                                         pct=ys,
                                         pct.l=ys.l,
                                         pct.h=ys.h))
    results <- NULL
  }
}
write.csv(cost.df, 'output/arn_costs.csv')

all.cost.df <- read.csv('output/arn_costs.csv')
diff.costs <- cost.vals[2] - cost.vals[1]
slides <- read_pptx('Plantilla_HPVinformationCentre.pptx')

for(tree in evaluated.trees) {
  cost.df <- all.cost.df[all.cost.df$strategy==tree,]
  # MINIMIZE.ICER = TRUE if the goal is to minimize target (i.e. desirable ICER as low as possible)
  # MINIMIZE.ICER = FALSE if the goal is to maximize target, or minimize the reference (i.e. ICER as high as possible)
  MINIMIZE.ICER <- T
  # if (!MINIMIZE.ICER) cost.df$pct <- 1 - cost.df$pct
  
  # lower.pcts <- cost.df[cost.df$wtp==22000,]
  upper.pcts <- cost.df[cost.df$wtp==25000,]
  labels <- sapply(cost.vals, function(cost) {
    # return(paste0(cost, ' € [', 
    #               formatC(lower.pcts[lower.pcts$cost==cost, 'pct'] * 100, digits=0, format = 'd'),
    #               '% - ',
    #               formatC(upper.pcts[upper.pcts$cost==cost, 'pct'] * 100, digits=0, format = 'd'),
    #               '%]'))
    return(paste0(cost, ' € [', 
                  # formatC(lower.pcts[lower.pcts$cost==cost, 'pct'] * 100, digits=0, format = 'd'),
                  # '% - ',
                  formatC(upper.pcts[upper.pcts$cost==cost, 'pct'] * 100, digits=0, format = 'd'),
                  '%]'))
  })
  
  brks <- c(seq(0,50000,10000), 25000)
  brks <- brks[order(brks)]
  
  p <- ggplot(cost.df, aes(x=wtp, y=pct, color=as.factor(cost))) + 
    geom_line() +
    # geom_ribbon(aes(ymin=pct.l, ymax=pct.h), alpha=.1) +
    # geom_vline(xintercept = 22000, linetype=2) +
    geom_vline(xintercept = 25000, linetype=2) +
    # annotate("rect", xmin = 22000, xmax = 25000, ymin = 0, ymax = 1, alpha = .3) +
    scale_x_continuous(breaks=brks, labels=formatC(brks, big.mark = ',', format='d')) +
    scale_y_continuous(breaks=seq(0,1,.1), labels=paste0(seq(0,1,.1)*100, '%'), limits = c(0,1), expand=c(0,0)) +
    scale_color_hue(name='ARN kit cost',
                         labels=labels) +
    xlab('WTP [€/QALY]') + 
    ylab('Percentage of cost-effective simulations for current strategy') +
    ggtitle(paste0(tree, ' vs conventional_t_tca')) +
    theme_minimal()
    
    # coord_cartesian(xlim=c(0,50000))
  p
  
  ggsave(paste0(getwd(), '/output/cost_sweep_psa_arnhr_', tree, '.pdf'), 
         p,
         width = 195,
         height = 130,
         units = 'mm',
         device = cairo_pdf
  )
  ggplotly(p)
  
  color <- '#3388ff'
  wtp.cost.df <- cost.df[cost.df$wtp==THRESHOLD.WTP,]
  if (MINIMIZE.ICER) {
    upper.range <- wtp.cost.df[wtp.cost.df$pct < .5,'cost']
  } else {
    upper.range <- wtp.cost.df[wtp.cost.df$pct > .5,'cost']
  }
  cost.upper <- upper.range[1]
  if (cost.upper != cost.vals[1]) {
    cost.lower <- cost.upper - diff.costs
    pct.upper <- wtp.cost.df[wtp.cost.df$cost==cost.upper, 'pct']
    pct.lower <- wtp.cost.df[wtp.cost.df$cost==cost.lower, 'pct']
    slope <- (pct.upper - pct.lower) / ( cost.upper - cost.lower)
    intercept <- pct.upper - slope * cost.upper
    p50.cost <- (.5 - intercept) / slope
  } else {
    p50.cost <- numeric(0)
  }
  
  base.cost <- strat.ctx[[1]]$c_arn_kit
  cost.upper <- (base.cost %/% diff.costs + 1) * diff.costs
  cost.lower <- base.cost %/% diff.costs * diff.costs
  pct.upper <- wtp.cost.df[wtp.cost.df$cost==cost.upper, 'pct']
  pct.lower <- wtp.cost.df[wtp.cost.df$cost==cost.lower, 'pct']
  slope <- (pct.upper - pct.lower) / ( cost.upper - cost.lower)
  intercept <- pct.upper - slope * cost.upper
  base.pct <- slope * base.cost + intercept
  
  base.cost2 <- 25
  cost.upper <- (base.cost2 %/% diff.costs + 1) * diff.costs
  cost.lower <- base.cost2 %/% diff.costs * diff.costs
  pct.upper <- wtp.cost.df[wtp.cost.df$cost==cost.upper, 'pct']
  pct.lower <- wtp.cost.df[wtp.cost.df$cost==cost.lower, 'pct']
  slope <- (pct.upper - pct.lower) / ( cost.upper - cost.lower)
  intercept <- pct.upper - slope * cost.upper
  base.pct2 <- slope * base.cost2 + intercept
  
  p2 <- ggplot(wtp.cost.df, aes(x=cost, y=pct)) + 
    geom_point(color=color) +
    geom_line(color=color) +
    geom_ribbon(aes(ymin=pct.l, ymax=pct.h), alpha=.25, fill=color) +
    geom_hline(yintercept = .5, linetype=3)
  if (length(p50.cost) > 0) {
    p2 <- p2 +
      geom_segment(x=p50.cost, y=.5, xend=p50.cost, yend=0, linetype=3) +
      geom_label(x=p50.cost,
                y=0,
                hjust=0,
                vjust=0,
                label=paste0('Threshold cost: ', formatC(p50.cost, digits = 2, format='f'), ' €'))
  }
  p2 <- p2 +
    geom_segment(x=base.cost, y=base.pct, xend=base.cost, yend=0, linetype=2) +
    geom_label(x=base.cost,
               y=0,
               hjust=1,
               vjust=0,
               label=paste0('Base cost: ', formatC(base.cost, digits = 0, format='d'), ' €')) +
    geom_label(x=base.cost-.2,
               y=base.pct,
               hjust=1,
               vjust=1,
               label=paste0(formatC(base.pct*100, digits = 1, format='f'), '%')) +
    geom_segment(x=base.cost2, y=base.pct2, xend=base.cost2, yend=0, linetype=2) +
    geom_label(x=base.cost2,
               y=0,
               hjust=1,
               vjust=0,
               label=paste0('Base cost: ', formatC(base.cost2, digits = 0, format='d'), ' €')) +
    geom_label(x=base.cost2-.2,
               y=base.pct2,
               hjust=1,
               vjust=1,
               label=paste0(formatC(base.pct2*100, digits = 1, format='f'), '%')) +
    scale_x_continuous(breaks=cost.vals, expand=c(.01, .01)) +
    scale_y_continuous(breaks=seq(0,1,.1), labels=paste0(seq(0,1,.1)*100, '%'), limits = c(0,1), expand=c(0,0)) +
    xlab('ARN kit cost [€]') +
    ylab('Percentage of cost-effective simulations') +
    ggtitle(paste0(tree, ' vs conventional_t_tca')) +
    theme_minimal()
  p2
  
  ggsave(paste0(getwd(), '/output/cost_sweep_psa_2_', tree, '.pdf'), 
         p2,
         width = 220,
         height = 160,
         units = 'mm',
         device = cairo_pdf
  )
  ggplotly(p2)
  
    slides <- slides %>%
    add_slide(layout='Titulo y objetos', master='Plantilla ICO') %>%
    ph_with(rvg::dml(ggobj=p), ph_location_type('body', id=1))
  
  slides <- slides %>%
    add_slide(layout='Titulo y objetos', master='Plantilla ICO') %>%
    ph_with(rvg::dml(ggobj=p2), ph_location_type('body', id=1))
  
}
slides %>% print(paste0('output/arn_cost_sweep.pptx'))
  