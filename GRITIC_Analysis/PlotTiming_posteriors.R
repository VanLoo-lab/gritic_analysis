# This plotting function is updated; now the input object is a full dataset to be filtered
# colnames needed: Anom_Sample_ID, Segment_ID, Segment_Start, Segment_End, Posterior_Sample_Index, Route, Gain_Timing, Adj
# Object needed: gp_adjust, timing_object

PlotTiming_posteriors <- function(timing_obj=timing_keep, sampleid="SAMPLE_ID"){
  s <- timing_obj %>% filter(Anom_Sample_ID %in% sampleid) %>% droplevels()
  sample_obj <- s %>% group_by(Segment_ID, Posterior_Sample_Index, Route, Anom_Sample_ID) %>%
    mutate(groupInd = with_order(order_by = Gain_Timing, fun = row_number, x = Gain_Timing))
  sample_obj$groupInd <- as.factor(sample_obj$groupInd)
  
  ggplot(sample_obj) + 
    # draw the segments
    geom_segment(aes(x=`Segment_Start`+`Adj`, xend=`Segment_End`+`Adj`, y=`Gain_Timing`, yend=`Gain_Timing`, color=groupInd), 
                 linewidth=0.1) + 
    scale_x_continuous(breaks = as.numeric(gp_adjust$Adj), 
                       labels = gp_adjust$Chromosome, 
                       expand = c(0, 0), 
                       limits = c(0, 3088269832)) + 
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1)) +
    guides(color = guide_legend(title="Gain\nNumber", override.aes = list(linewidth=1)))+
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_line(linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
    labs(title=paste(sampleid, sample_obj$Cancer_Type[[1]], sep = ", "),
         x ="Genome Position", y = "Gain Timing") + 
    theme(axis.line = element_line(color = 'black'), plot.title = element_text(hjust = 0.5)) +
    geom_abline(aes(slope = 0, intercept = mean(sample_obj$WGD_Timing), linetype="WGD"), color="red") + 
    scale_linetype_manual(values = c("WGD" = "dashed"), name = NULL) + 
    theme(text = element_text(family = "Times New Roman", face = "bold"),
          axis.text.x = element_text(size = 8))
}
