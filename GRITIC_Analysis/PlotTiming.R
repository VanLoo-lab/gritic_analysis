PlotTiming <- function(s, gp_adj = gp_adjust, linewidth=0.1){
  ggplot(s) + 
    # draw the segments
    geom_segment(aes(x=`Segment_Start`+`Adj`, xend=`Segment_End`+`Adj`, y=`Gain_Timing`, yend=`Gain_Timing`), colour="black", linewidth=linewidth) + 
    scale_x_continuous(breaks = as.numeric(gp_adj$Adj), labels = gp_adj$Chromosome, expand = c(0, 0), limits = c(0, 3088269832)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_line(linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
    theme(axis.line = element_line(color = 'black'))
}
