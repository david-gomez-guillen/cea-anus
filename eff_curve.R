library(ggplot2)

cost.label <- '€'
eff.label <- 'QALY'
nudge.x.icer <- .01
nudge.y.icer <- .005

x.label <- ifelse(is.null(cost.label), "Cost", paste0("Cost [€]"))
y.label <- ifelse(is.null(eff.label), "Effectiveness", 
                  paste0("Effectiveness [QALY]"))
ce.label <- ifelse(is.null(cost.label) && is.null(eff.label), 
                   "", paste0(" ", cost.label, "/", eff.label))
fullSummary <- markov.result$summary

# fullSummary <- fullSummary[grepl('la', fullSummary$strategy, perl=TRUE) | grepl('diff_hpv16', fullSummary$strategy, perl=TRUE),]
# fullSummary <- fullSummary[grepl('arnm', fullSummary$strategy, perl=TRUE) | 
#                              grepl('^conventional-', fullSummary$strategy, perl=TRUE) |
#                              grepl('^conventional_t-', fullSummary$strategy, perl=TRUE),]
# fullSummary <- fullSummary[grepl('_t$', fullSummary$strategy, perl=TRUE),]
fullSummary <- analyzeCE(fullSummary, plot=TRUE)$summary

plot.df <- fullSummary
undominated.df <- fullSummary[fullSummary$domination == 
                                "undominated", ]
if (nrow(undominated.df) > 0) {
  undominated.df$label <- paste0("ICER=", formatC(round(undominated.df$ICER, 
                                                        digits = 2), big.mark = ",", format = "d"), ce.label)
  undominated.df[1, "label"] <- ""
} else {
  undominated.df$label <- numeric(0)
}
y.range <- diff(range(plot.df$E))
x.range <- diff(range(plot.df$C))
plt <- ggplot2::ggplot(plot.df, ggplot2::aes(x = C, 
                                             y = E)) + ggplot2::geom_line(data = undominated.df) + 
  ggplot2::geom_text(data = undominated.df, mapping = ggplot2::aes(label = label), 
                     nudge_x = -nudge.x.icer * x.range, nudge_y = nudge.y.icer * 
                       y.range)

strat.names <- c(
  'conventional-semestral_followup',
  'conventional_hpv16la-semestral_followup',
  'conventional_hpv1618la-semestral_followup',
  'conventional_hpvhrla-semestral_followup',
  'conventional_hpvhrhc-semestral_followup',
  'arnme6e7_hpv16-semestral_followup',
  'arnme6e7_hpv161845-semestral_followup',
  'arnme6e7_hpvhr-semestral_followup',
  'ascus_lsil_diff_hpv16la-semestral_followup',
  'ascus_lsil_diff_hpv1618la-semestral_followup',
  'ascus_lsil_diff_hpvhrla-semestral_followup',
  'ascus_lsil_diff_hpvhrhc-semestral_followup',
  'ascus_lsil_diff_arnme6e7_16-semestral_followup',
  'ascus_lsil_diff_arnme6e7_161845-semestral_followup',
  'ascus_lsil_diff_arnme6e7_hr-semestral_followup',
  'conventional_t_irc-semestral_followup_t_irc',
  'conventional_hpv16la_t_irc-semestral_followup_t_irc',
  'conventional_hpv1618la_t_irc-semestral_followup_t_irc',
  'conventional_hpvhrla_t_irc-semestral_followup_t_irc',
  'conventional_hpvhrhc_t_irc-semestral_followup_t_irc',
  'arnme6e7_hpv16_t_irc-semestral_followup_t_irc',
  'arnme6e7_hpv161845_t_irc-semestral_followup_t_irc',
  'arnme6e7_hpvhr_t_irc-semestral_followup_t_irc',
  'ascus_lsil_diff_hpv16la_t_irc-semestral_followup_t_irc',
  'ascus_lsil_diff_hpv1618la_t_irc-semestral_followup_t_irc',
  'ascus_lsil_diff_hpvhrla_t_irc-semestral_followup_t_irc',
  'ascus_lsil_diff_hpvhrhc_t_irc-semestral_followup_t_irc',
  'ascus_lsil_diff_arnme6e7_16_t_irc-semestral_followup_t_irc',
  'ascus_lsil_diff_arnme6e7_161845_t_irc-semestral_followup_t_irc',
  'ascus_lsil_diff_arnme6e7_hr_t_irc-semestral_followup_t_irc',
  'conventional_t_tca-semestral_followup_t_tca',
  'conventional_hpv16la_t_tca-semestral_followup_t_tca',
  'conventional_hpv1618la_t_tca-semestral_followup_t_tca',
  'conventional_hpvhrla_t_tca-semestral_followup_t_tca',
  'conventional_hpvhrhc_t_tca-semestral_followup_t_tca',
  'arnme6e7_hpv16_t_tca-semestral_followup_t_tca',
  'arnme6e7_hpv161845_t_tca-semestral_followup_t_tca',
  'arnme6e7_hpvhr_t_tca-semestral_followup_t_tca',
  'ascus_lsil_diff_hpv16la_t_tca-semestral_followup_t_tca',
  'ascus_lsil_diff_hpv1618la_t_tca-semestral_followup_t_tca',
  'ascus_lsil_diff_hpvhrla_t_tca-semestral_followup_t_tca',
  'ascus_lsil_diff_hpvhrhc_t_tca-semestral_followup_t_tca',
  'ascus_lsil_diff_arnme6e7_16_t_tca-semestral_followup_t_tca',
  'ascus_lsil_diff_arnme6e7_161845_t_tca-semestral_followup_t_tca',
  'ascus_lsil_diff_arnme6e7_hr_t_tca-semestral_followup_t_tca'
)

strat.labels <- c(
  'Conventional',
  'HPV-16 LA',
  'HPV-16/18 LA',
  'HPV-HR LA',
  'HPV HC',
  'ARNmE6/E7 (HPV-16)',
  'ARNmE6/E7 (HPV-16/18/45)',
  'ARNmE6/E7 (HPV-HR)',
  'ASCUS/LSIL diff (HPV-16 LA)',
  'ASCUS/LSIL diff (HPV-16/18 LA)',
  'ASCUS/LSIL diff (HPV-HR LA)',
  'ASCUS/LSIL diff (HPV HC)',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16)',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16/18/45)',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-HR)',
  'Conventional [T=IRC]',
  'HPV-16 LA [T=IRC]',
  'HPV-16/18 LA [T=IRC]',
  'HPV-HR LA [T=IRC]',
  'HPV HC [T=IRC]',
  'ARNmE6/E7 (HPV-16) [T=IRC]',
  'ARNmE6/E7 (HPV-16/18/45) [T=IRC]',
  'ARNmE6/E7 (HPV-HR) [T=IRC]',
  'ASCUS/LSIL diff (HPV-16 LA) [T=IRC]',
  'ASCUS/LSIL diff (HPV-16/18 LA) [T=IRC]',
  'ASCUS/LSIL diff (HPV-HR LA) [T=IRC]',
  'ASCUS/LSIL diff (HPV HC) [T=IRC]',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16) [T=IRC]',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16/18/45) [T=IRC]',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-HR) [T=IRC]',
  'Conventional [T=TCA]',
  'HPV-16 LA [T=TCA]',
  'HPV-16/18 LA [T=TCA]',
  'HPV-HR LA [T=TCA]',
  'HPV HC [T=TCA]',
  'ARNmE6/E7 (HPV-16) [T=TCA]',
  'ARNmE6/E7 (HPV-16/18/45) [T=TCA]',
  'ARNmE6/E7 (HPV-HR) [T=TCA]',
  'ASCUS/LSIL diff (HPV-16 LA) [T=TCA]',
  'ASCUS/LSIL diff (HPV-16/18 LA) [T=TCA]',
  'ASCUS/LSIL diff (HPV-HR LA) [T=TCA]',
  'ASCUS/LSIL diff (HPV HC) [T=TCA]',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16) [T=TCA]',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16/18/45) [T=TCA]',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-HR) [T=TCA]'
)

names(strat.labels) <- strat.names

# strat.shapes <- c(22, 25, 25, 24, 24, 23, 23, 23, 23, 23, 23)
strat.shapes <- c(rep(22, 5), rep(23, 3), rep(24, 7), rep(22, 5), rep(23, 3), rep(24, 7), rep(22, 5), rep(23, 3), rep(24, 7))
names(strat.shapes) <- strat.names

# strat.fill <- c('black', '#00ccff', '#0000ff', '#00ff00', '#007700', '#6f0e18', '#961b1e', '#b91f26', '#ee3f3f', '#f48e8a', '#f9c7c2')
strat.fill <- rep(c(
  '#ffffcc',
  '#ffeda0',
  '#fed976',
  '#feb24c',
  '#fd8d3c',
  '#abe098', # ARN
  '#57c84d',
  '#6ef62c',
  '#eff3ff', # ASCUSS diff
  '#c6dbef',
  '#9ecae1',
  '#6baed6',
  '#4292c6',
  '#2171b5',
  '#084594'
  ), 3
)
names(strat.fill) <- strat.names

strat.colors <- c(rep('black', 15), rep('red', 15), rep('blue', 15))
names(strat.colors) <- strat.names

strat.indices <- strat.names %in% fullSummary$strategy
strat.names <- strat.names[strat.indices]
strat.colors <- strat.colors[strat.indices]
strat.fill <- strat.fill[strat.indices]
strat.shapes <- strat.shapes[strat.indices]

fullSummary$strategy <- factor(fullSummary$strategy, levels=strat.names, ordered=TRUE)

plt <- plt + 
  geom_point(data=fullSummary, size=5, aes(fill=strategy, color=strategy, shape=strategy), stroke=1) + 
  scale_shape_manual(values=strat.shapes,
                     labels=strat.labels,
                     name='Strategy') +
  scale_color_manual(values=strat.colors, 
                     labels=strat.labels,
                     name='Strategy') +
  scale_fill_manual(values=strat.fill,
                    labels=strat.labels,
                    name='Strategy') +
  theme(legend.position = 'bottom')
  # coord_cartesian(ylim=c(16.94, 16.96))
print(plt)
ggplotly(plt)
