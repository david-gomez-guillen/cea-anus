library(dplyr)
library(officer)
#setwd('~/Documents/models_ce/anus')
source('load_models.R')
source('markov.R')

cat('Running markov model...\n')
markov.outputs <- list()

strat.ctx <- lapply(strat.ctx, function(ctx) {ctx$p_vaccination <- .2; ctx})
strat.ctx <- refresh.context('p_vaccination', strat.ctx, excel.strata.df)

vac.df <- data.frame()
for(coverage in seq(0, 1, .2)) {
  cat(paste0('Coverage: ', coverage, '\n'))
  for(vac.eff in seq(0, 1, .05)) {
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
                              strategies$hiv_msm[c('conventional_t_tca', 'arnme6e7_hpvhr_t_tca')],
                              # strategies$hiv_msm[c('arnme6e7_hpvhr_t_tca', 'arnme6e7_hpvhr_t_irc')],
                              # strategies$hiv_msm[c('arnme6e7_hpvhr_t_tca', 'conventional_t_tca')],
                              markov,
                              strat.ctx,
                              initial.state,
                              discount.rate=.03)
    sink()
    
    vac.df <- rbind(vac.df, data.frame(coverage=coverage, vac.eff=vac.eff, icer=markov.result$summary[2, 'ICER']))
    # Remove redundant suffix for strategy names
    # markov.result$summary$strategy <- gsub('(.*?)-(.*)', '\\1', markov.result$summary$strategy)
    
    # ggplotly(markov.result$plot)
  }
}

write.csv(vac.df, 'output/vac_df1.csv', row.names = FALSE)

vac.df <- read.csv('output/vac_df1.csv')

slides <- read_pptx('Plantilla_HPVinformationCentre.pptx')

vac.plt <- ggplot(vac.df, aes(x=vac.eff, y=icer, color=as.factor(coverage))) + 
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

slides %>% print(paste0('output/vaccination.pptx'))
