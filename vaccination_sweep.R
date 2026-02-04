library(dplyr)
library(officer)
#setwd('~/Documents/models_ce/anus')
source('load_models.R')
source('markov.R')

cat('Running markov model...\n')
markov.outputs <- list()

unlist.strat.ctx <- unlist(strat.ctx)

# strat.ctx <- lapply(strat.ctx, function(ctx) {ctx$p_vaccination <- .2; ctx})
# strat.ctx <- refresh.context('p_vaccination', strat.ctx, excel.strata.df)

# EXPERIMENTS <- list(
#   list(
#     reference='conventional_t_tca',
#     strategy='arnme6e7_hpvhr_t_tca',
#     file.suffix='cito_vs_cotest_arn'
#   ),
#   list(
#     reference='conventional_t_tca',
#     strategy='conventional_hpvhrla_t_tca',
#     file.suffix='cito_vs_cotest_pcr'
#   ),
#   list(
#     reference='conventional_t_tca',
#     strategy='conventional_hpvhrhc_t_tca',
#     file.suffix='cito_vs_cotest_hyb'
#   ),
#   list(
#     reference='conventional_t_tca',
#     strategy='ascus_lsil_diff_arnme6e7_hpvhr_t_tca',
#     file.suffix='cito_vs_reflex_arn'
#   ),
#   list(
#     reference='conventional_t_tca',
#     strategy='ascus_lsil_diff_hpvhrla_t_tca',
#     file.suffix='cito_vs_reflex_pcr'
#   ),
#   list(
#     reference='conventional_t_tca',
#     strategy='ascus_lsil_diff_hpvhrhc_t_tca',
#     file.suffix='cito_vs_reflex_hyb'
#   )
# )

EXPERIMENTS <- list(
  list(
    reference='conventional_t_tca',
    strategy='arnme6e7_hpvhr_t_tca',
    file.suffix='c6_scenario1_cito_vs_arn'
  ),
  list(
    reference='conventional_t_tca',
    strategy='conventional_hpvhrla_t_tca',
    file.suffix='c6_scenario1_cito_vs_pcr'
  ),
  list(
    reference='conventional_t_tca',
    strategy='conventional_hpvhrhc_t_tca',
    file.suffix='c6_scenario1_cito_vs_hyb'
  )
)

coverage.values <- seq(0, 1, .2)
eff.values <- seq(0, 1, .1)
eff.values <- eff.values[order(eff.values)]

for(exp in EXPERIMENTS) {
  reference <- exp$reference
  strategy <- exp$strategy
  file.suffix <- exp$file.suffix
  
  cat(paste0('Experiment: ', reference, ' vs ', strategy, '\n'))

  vac.df <- data.frame()
  for(coverage in coverage.values) {
    cat(paste0('Coverage: ', coverage, '\n'))
    for(vac.eff in eff.values) {
      cat(paste0('Efficacy: ', vac.eff, '\n'))
      for(ctx.name in names(strat.ctx)) {
        strat.ctx[[ctx.name]]$p_vaccination <- vac.eff * coverage
      }
      strat.ctx <- refresh.context('p_vaccination', strat.ctx, excel.strata.df)
      initial.state <- sapply(markov$nodes,
                              function(n) if (n$name=='hiv_positive') 1 else 0)
      sink('/dev/null')
      markov.result <- simulate('hiv_msm',
                                # strategies$hiv_msm,
                                # strategies$hiv_msm['no_intervention'],
                                # strategies$hiv_msm[c('conventional_hpvhrla_t_tca', 'arnme6e7_hpvhr_t_tca')],

                                # strategies$hiv_msm[c('conventional_t_irc', 'conventional_t_tca')],
                                strategies$hiv_msm[c(reference, strategy)],
                                # strategies$hiv_msm[c('ascus_lsil_diff_hpvhrla_t_tca', 'conventional_hpvhrhc_t_tca')],
                                # strategies$hiv_msm[c('conventional_hpvhrhc_t_tca', 'conventional_hpvhrla_t_tca')],

                                # strategies$hiv_msm[c('arnme6e7_hpvhr_t_tca', 'arnme6e7_hpvhr_t_irc')],
                                # strategies$hiv_msm[c('arnme6e7_hpvhr_t_tca', 'conventional_t_tca')],
                                markov,
                                strat.ctx,
                                initial.state,
                                discount.rate=.03)
      sink()

      markov.result$summary$strategy <- gsub('(.*?)-(.*)', '\\1', markov.result$summary$strategy)
      row.names(markov.result$summary) <- markov.result$summary$strategy
      
      ic <- markov.result$summary[markov.result$summary$strategy == strategy, 'C'] - markov.result$summary[markov.result$summary$strategy == reference, 'C']
      ie <- markov.result$summary[markov.result$summary$strategy == strategy, 'E'] - markov.result$summary[markov.result$summary$strategy == reference, 'E']
      icer <- ic / ie
      nhb <- ie - ic/25000

      vac.df <- rbind(vac.df, data.frame(coverage=coverage, vac.eff=vac.eff, ic=ic, ie=ie, icer=icer, nhb=nhb))
      # Remove redundant suffix for strategy names
      # markov.result$summary$strategy <- gsub('(.*?)-(.*)', '\\1', markov.result$summary$strategy)
      
      # ggplotly(markov.result$plot)
    }
  }

  date.suffix <- format(Sys.Date(), '%Y%m%d')
  sheet.data <- list(
    'Results'=vac.df,
    'Base parameters'=data.frame(
      parameter=c(names(unlist.strat.ctx), 'discount'), 
      base.value=c(unlist.strat.ctx, 0.03)
    )
  )
  openxlsx::write.xlsx(sheet.data, 
                      paste0('output/vaccination_', file.suffix, '__', date.suffix, '.xlsx'), 
                      rowNames = F, 
                      colWidths='auto')


  slides <- read_pptx('Plantilla_HPVinformationCentre.pptx')

  vac.plt <- ggplot(vac.df, aes(x=coverage, y=icer, color=as.factor(vac.eff))) + 
    geom_line() + 
    geom_hline(yintercept = 25000, linetype='dashed') +
    xlab('Vaccine coverage') +
    ylab('ICER [€/QALY]') +
    theme_minimal() + 
    coord_cartesian(xlim=c(0,1), ylim=c(0, 100000)) +
    scale_x_continuous(labels=function(x) paste0(as.numeric(x)*100, '%')) +
    scale_y_continuous(labels=function(y) formatC(y, format='d', big.mark = ',')) +
    scale_color_discrete(name='Vaccine effectiveness', labels=function(x) paste0(as.numeric(x)*100, '%'))

  slides <- slides %>%
    add_slide(layout='Titulo y objetos', master='Plantilla ICO') %>%
    ph_with(rvg::dml(ggobj=vac.plt), ph_location_type('body', id=1))

  vac.plt2 <- ggplot(vac.df, aes(x=coverage, y=nhb, color=as.factor(vac.eff))) + 
    geom_line() + 
    geom_hline(yintercept = 0, linetype='dashed') +
    xlab('Vaccine coverage') +
    ylab('NHB [QALY]') +
    theme_minimal() + 
    coord_cartesian(xlim=c(0,1), ylim=c(min(min(vac.df$nhb),0),min(max(vac.df$nhb),1))) +
    scale_x_continuous(labels=function(x) paste0(as.numeric(x)*100, '%')) +
    #scale_y_continuous(labels=function(y) formatC(y, format='d', big.mark = ',')) +
    scale_color_discrete(name='Vaccine effectiveness', labels=function(x) paste0(as.numeric(x)*100, '%'))

  slides <- slides %>%
    add_slide(layout='Titulo y objetos', master='Plantilla ICO') %>%
    ph_with(rvg::dml(ggobj=vac.plt2), ph_location_type('body', id=1))

  slides %>% print(paste0('output/vaccination_', file.suffix, '__', date.suffix, '.pptx'))
}