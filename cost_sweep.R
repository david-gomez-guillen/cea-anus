library(dplyr)
# library(Cairo)
library(ggplot2)
library(gridExtra)
library(scales)
library(data.table)
library(officer)
library(openxlsx)

SUFFIX <- ''
N.CORES.DEFAULT <- 5
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
N <- 20
PAR <- 'c_arn_kit'
COST.SWEEP <- seq(0, 30, 2)

# evaluated.trees <- c(
#   'arnme6e7_hpvhr_t_tca'
# )
evaluated.trees <- list(
  #  list(reference='conventional_t_tca',
  #     strategy='conventional_hpvhrla_t_tca')
    # list(reference='conventional_t_tca',
    #    strategy='arnme6e7_hpvhr_t_tca')
  #  list(reference='conventional_hpvhrhc_t_tca',
  #     strategy='conventional_hpvhrla_t_tca')
    #  list(reference='conventional_hpvhrhc_t_tca',
    #   strategy='arnme6e7_hpvhr_t_tca')
      list(reference='conventional_hpvhrla_t_tca',
       strategy='arnme6e7_hpvhr_t_tca')
)

for(tree.pair in evaluated.trees) {
  for(i in COST.SWEEP) {
    sweep.strat.ctx <- lapply(strat.ctx, function(ctx) {
      ctx[[PAR]] <- i
      return(ctx)
    })
    cat('+++++ PSA for', PAR, '= ', i, '\n')
    results <- psa.1(PAR, 
                     sweep.strat.ctx,
                     excel.file = 'params/context.xlsx',
                     population='hiv_msm', 
                     strategy=tree.pair$strategy,
                     markov=markov,
                     reference=tree.pair$reference,
                     sd.estimate.func=sd.estimate.func,
                     context.setup.func=context.setup,
                     n.cores=N.CORES,
                     n.iters=N,
                     seed=1234,
                     discount.rate=.03)
    
    df <- results[[PAR]]$summary
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
    cost.df <- rbind(cost.df, data.frame(strategy=tree.pair$strategy,
                                         cost=i,
                                         wtp=xs,
                                         pct=ys,
                                         pct.l=ys.l,
                                         pct.h=ys.h))
    results <- NULL
  }
}

unlist.strat.ctx <- unlist(strat.ctx)
date.suffix <- format(Sys.Date(), '%Y%m%d')
sheet.data <- list(
  'Results'=cost.df,
  'Base parameters'=data.frame(
    parameter=c(names(unlist.strat.ctx), 'discount'), 
    base.value=c(unlist.strat.ctx, 0.03)
  )
)
openxlsx::write.xlsx(sheet.data, 
                    paste0('output/costs__', date.suffix, '.xlsx'), 
                    rowNames = F, 
                    colWidths='auto')


all.cost.df <- openxlsx::read.xlsx(paste0('output/costs__', date.suffix, '.xlsx'), sheet='Results')
# all.cost.df <- openxlsx::read.xlsx('output/costs_c6_hpvhrla_vs_mrna__20251014.xlsx', sheet='Results')
diff.costs <- COST.SWEEP[2] - COST.SWEEP[1]
slides <- read_pptx('Plantilla_HPVinformationCentre.pptx')

for(tree.pair in evaluated.trees) {
  cost.df <- all.cost.df[all.cost.df$strategy==tree.pair$strategy,]
  # MINIMIZE.ICER = TRUE if the goal is to minimize target (i.e. desirable ICER as low as possible)
  # MINIMIZE.ICER = FALSE if the goal is to maximize target, or minimize the reference (i.e. ICER as high as possible)
  MINIMIZE.ICER <- T
  # if (!MINIMIZE.ICER) cost.df$pct <- 1 - cost.df$pct
  
  # lower.pcts <- cost.df[cost.df$wtp==22000,]
  upper.pcts <- cost.df[cost.df$wtp==25000,]
  labels <- sapply(COST.SWEEP, function(cost) {
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
    scale_color_hue(name=PAR,
                         labels=labels) +
    xlab('WTP [€/QALY]') + 
    ylab('Percentage of cost-effective simulations for current strategy') +
    ggtitle(paste0(tree.pair$strategy, ' vs ', tree.pair$reference)) +
    theme_minimal()
    
    # coord_cartesian(xlim=c(0,50000))
  p
  
  # ggsave(paste0(getwd(), '/output/cost_sweep_psa_arnhr_', tree, '.pdf'), 
  #        p,
  #        width = 195,
  #        height = 130,
  #        units = 'mm',
  #        device = cairo_pdf
  # )
  
  # grDevices::cairo_pdf(paste0(getwd(), '/output/cost_sweep_mrna_', tree.pair$strategy, '_vs_', tree.pair$reference, '__', date.suffix, '.pdf'),
  #                       width = 9,
  #                       height= 9)
  # print(p)
  # dev.off()
  # ggplotly(p)
  
  color <- '#3388ff'
  wtp.cost.df <- cost.df[cost.df$wtp==THRESHOLD.WTP,]
  if (MINIMIZE.ICER) {
    upper.range <- wtp.cost.df[wtp.cost.df$pct < .5,'cost']
  } else {
    upper.range <- wtp.cost.df[wtp.cost.df$pct > .5,'cost']
  }
  cost.upper <- upper.range[1]
  if (!is.na(cost.upper) && cost.upper != COST.SWEEP[1]) {
    cost.lower <- cost.upper - diff.costs
    pct.upper <- wtp.cost.df[wtp.cost.df$cost==cost.upper, 'pct']
    pct.lower <- wtp.cost.df[wtp.cost.df$cost==cost.lower, 'pct']
    slope <- (pct.upper - pct.lower) / ( cost.upper - cost.lower)
    intercept <- pct.upper - slope * cost.upper
    p50.cost <- (.5 - intercept) / slope
  } else {
    p50.cost <- numeric(0)
  }
  
  base.cost <- strat.ctx[[1]][[PAR]]
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
    scale_x_continuous(breaks=COST.SWEEP, expand=c(.01, .01)) +
    scale_y_continuous(breaks=seq(0,1,.1), labels=paste0(seq(0,1,.1)*100, '%'), limits = c(0,1), expand=c(0,0)) +
    xlab(paste0(PAR, ' [€]')) +
    ylab('Percentage of cost-effective simulations') +
    ggtitle(paste0(tree.pair$strategy, ' vs ', tree.pair$reference)) +
    theme_minimal()
  p2
  
  # ggsave(paste0(getwd(), '/output/cost_sweep_psa_2_', tree, '.pdf'), 
  #        p2,
  #        width = 220,
  #        height = 160,
  #        units = 'mm',
  #        device = cairo_pdf
  # )
  
  # grDevices::cairo_pdf(paste0(getwd(), '/output/cost_sweep2_mrna_', tree.pair$strategy, '_vs_', tree.pair$reference, '__', date.suffix, '.pdf'),
  #                       width = 9,
  #                       height= 9)
  # print(p2)
  # dev.off()
  # ggplotly(p2)
  
    slides <- slides %>%
    add_slide(layout='Titulo y objetos', master='Plantilla ICO') %>%
    ph_with(rvg::dml(ggobj=p), ph_location_type('body', type_idx=1))
  
  slides <- slides %>%
    add_slide(layout='Titulo y objetos', master='Plantilla ICO') %>%
    ph_with(rvg::dml(ggobj=p2), ph_location_type('body', type_idx=1))
  
}
slides %>% print(paste0('output/cost_sweep__', date.suffix, '.pptx'))
  