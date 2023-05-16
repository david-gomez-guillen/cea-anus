library(ggplot2)

cost.label <- '€'
eff.label <- 'QALY'
nudge.x.icer <- .01
nudge.y.icer <- .01

x.label <- ifelse(is.null(cost.label), "Cost", paste0("Cost [€]"))
y.label <- ifelse(is.null(eff.label), "Effectiveness", 
                  paste0("Effectiveness [QALY]"))
ce.label <- ifelse(is.null(cost.label) && is.null(eff.label), 
                   "", paste0(" ", cost.label, "/", eff.label))
fullSummary <- markov.result$summary
plot.df <- fullSummary
undominated.df <- fullSummary[fullSummary$domination == 
                                "undominated", ]
undominated.df$label <- paste0("ICER=", formatC(round(undominated.df$ICER, 
                                                      digits = 2), big.mark = ",", format = "d"), ce.label)
undominated.df[1, "label"] <- ""
y.range <- diff(range(plot.df$E))
x.range <- diff(range(plot.df$C))
plt <- ggplot2::ggplot(plot.df, ggplot2::aes(x = C, 
                                             y = E)) + ggplot2::geom_line(data = undominated.df) + 
  ggplot2::geom_text(data = undominated.df, mapping = ggplot2::aes(label = label), 
                     nudge_x = -nudge.x.icer * x.range, nudge_y = nudge.y.icer * 
                       y.range)

strat.names <- c(
  'conventional',
  'conventional_hpv16la',
  'conventional_hpv1618la',
  'conventional_hpvhrla',
  'conventional_hpvhrhc',
  'arnme6e7_hpv16',
  'arnme6e7_hpv161845',
  'arnme6e7_hpvhr',
  'ascus_lsil_diff_hpv16',
  'ascus_lsil_diff_hpv1618',
  'ascus_lsil_diff_hpvhrla',
  'ascus_lsil_diff_hpvhrhc',
  'ascus_lsil_diff_arn16',
  'ascus_lsil_diff_arn161845',
  'ascus_lsil_diff_arnhr',
  'conventional_t',
  'conventional_hpv16la_t',
  'conventional_hpv1618la_t',
  'conventional_hpvhrla_t',
  'conventional_hpvhrhc_t',
  'arnme6e7_hpv16_t',
  'arnme6e7_hpv161845_t',
  'arnme6e7_hpvhr_t',
  'ascus_lsil_diff_hpv16_t',
  'ascus_lsil_diff_hpv1618_t',
  'ascus_lsil_diff_hpvhrla_t',
  'ascus_lsil_diff_hpvhrhc_t',
  'ascus_lsil_diff_arn16_t',
  'ascus_lsil_diff_arn161845_t',
  'ascus_lsil_diff_arnhr_t'
)

strat.labels <- c(
  'Conventional',
  'Conventional (HPV-16)',
  'Conventional (HPV-16/18)',
  'Conventional (HPV-HR LA)',
  'Conventional (HPV-HR HC)',
  'ARNmE6/E7 (HPV-16)',
  'ARNmE6/E7 (HPV-16/18/45)',
  'ARNmE6/E7 (HPV-HR)',
  'ASCUS/LSIL diff (HPV-16)',
  'ASCUS/LSIL diff (HPV-16/18)',
  'ASCUS/LSIL diff (HPV-HR LA)',
  'ASCUS/LSIL diff (HPV-HR HC)',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16)',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16/18/45)',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-HR)',
  'Conventional [T]',
  'Conventional (HPV-16) [T]',
  'Conventional (HPV-16/18) [T]',
  'Conventional (HPV-HR LA) [T]',
  'Conventional (HPV-HR HC) [T]',
  'ARNmE6/E7 (HPV-16) [T]',
  'ARNmE6/E7 (HPV-16/18/45) [T]',
  'ARNmE6/E7 (HPV-HR) [T]',
  'ASCUS/LSIL diff (HPV-16) [T]',
  'ASCUS/LSIL diff (HPV-16/18) [T]',
  'ASCUS/LSIL diff (HPV-HR LA) [T]',
  'ASCUS/LSIL diff (HPV-HR HC) [T]',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16) [T]',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-16/18/45) [T]',
  'ASCUS/LSIL diff (ARNmE6/E7 HPV-HR) [T]'
)

names(strat.labels) <- strat.names

# strat.shapes <- c(22, 25, 25, 24, 24, 23, 23, 23, 23, 23, 23)
strat.shapes <- c(rep(22, 5), rep(23, 3), rep(24, 7), rep(22, 5), rep(23, 3), rep(24, 7))
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
  ), 2
)
names(strat.fill) <- strat.names

strat.colors <- c(rep('black', 15), rep('red', 15))
names(strat.colors) <- strat.names

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
                    name='Strategy')
print(plt)
