library(tidyverse)
library(magrittr)
library(data.table)
library(ggpubr)
library(patchwork)

library(extrafont)
extrafont::loadfonts()

setwd("/Users/slai4/Library/CloudStorage/OneDrive-InsideMDAnderson/VanLooLab/SyncGain/version_04092024")

source("../Scripts/perm_plots.R")
source("../Scripts/PlotTiming_posteriors.R")
source("../Scripts/PlotTiming.R")

# Format data --------
# run this part only once and save results
temp_gsize <- read.csv("../resources/genome_size.csv")
gp_adjust <- data.frame(Chromosome = temp_gsize$Chromosome, Adj = temp_gsize$Position-temp_gsize$Length)
rm(temp_gsize)

non_parsimony_penalty <- 'True'
timing <- read_tsv(paste0("../in_data/posterior_table_filtered_penalty_", non_parsimony_penalty, ".tsv"))
timing$Chromosome <- as.character(timing$Chromosome)
timing <- left_join(timing, gp_adjust, by = c("Chromosome"="Chromosome"))
timing$Center <- round(timing$Segment_Start+timing$Segment_End)/2 + timing$Adj
colnames(timing) # make sure the IDs are compatible with the script
colnames(timing)[16] <- "Anom_Sample_ID"

# keep cancer types with at leat 15 samples
cancer_id <- timing %>% group_by(Cancer_Type, Anom_Sample_ID) %>% summarise(n()) %>% 
  {table(.$Cancer_Type)} %>% as.data.frame() %>% filter(Freq>=15) %>% pull(Var1) %>% as.character()
timing_freq15 <- timing %>% filter(Cancer_Type %in% cancer_id)

sample_class <- data.frame(Anom_Sample_ID=unique(timing_freq15$Anom_Sample_ID))
WGD <- timing_freq15 %>% group_by(Anom_Sample_ID) %>% slice(1) %>% select(Anom_Sample_ID, WGD_Status)

# keep informative samples with at least 5 chrs with gains
sample_keep <- timing_freq15 %>% group_by(Anom_Sample_ID) %>% summarise(length(unique(Chromosome))) %>% filter(`length(unique(Chromosome))`>=5) %>% pull(Anom_Sample_ID) %>% as.character()
sinfo <- data.frame(Anom_Sample_ID=sample_keep, Info="Informative")

timing_keep <- timing_freq15 %>% filter(Anom_Sample_ID %in% sample_keep) %>% droplevels()
colnames(timing_keep)
colnames(timing_keep)[2] <- "Chr"
colnames(timing_keep)
cols <- c("Segment_ID", "Posterior_Sample_Index", "Chr", "Anom_Sample_ID", "Route")
timing_keep %<>% mutate_each_(funs(factor(.)),cols)
saveRDS(timing_keep, file = paste0("timing_keep_non_parsimony_penalty_",non_parsimony_penalty,"_04092024.rds"))
rm(timing, timing_keep, timing_freq15)

timing <- readRDS(paste0("timing_keep_non_parsimony_penalty_",non_parsimony_penalty,"_04092024.rds"))
sample_class <- left_join(sample_class, sinfo)
sample_class[is.na(sample_class$Info),]$Info <- "Uninformative"
sample_class <- left_join(sample_class, WGD)

# No need to do this anymore; the updated data only contains good samples
# generate subset with only good samples ----
# bad <- read_csv("../version_0606/lesky_posterior_export/bad_samples.tsv")
# timing_good <- filter(timing, !Anom_Sample_ID %in% bad$Sample_ID) %>% droplevels()
# length(levels(timing$Anom_Sample_ID))
# length(levels(timing_good$Anom_Sample_ID))
# saveRDS(timing_good, file = "timing_keep_non_parsimony_penalty_true_good_04092024.rds")
# # only run this if the plots are for good samples only
# sinfo <- filter(sinfo, !Anom_Sample_ID %in% bad$Sample_ID)
# WGD <- filter(WGD, !Anom_Sample_ID %in% bad$Sample_ID)
# sample_class <- filter(sample_class, !Anom_Sample_ID %in% bad$Sample_ID)

# Run weighted distance and permutation on server ------
min_dist_real <- read_delim(paste0("../in_data/min_dist_median_non_parsimony_penalty_",non_parsimony_penalty,".txt"),delim = "\t", col_names = FALSE)
colnames(min_dist_real) <- c("Anom_Sample_ID", "Min_dist", "Timing_Sync")
sample_id <- min_dist_real$Anom_Sample_ID

# notice that the function was updated to print out ID, sync_level and timing; be careful when generating perm_non_parsimony_penalty_*.txt
perm_WD <- read_delim(paste0("../in_data/perm_non_parsimony_penalty_",non_parsimony_penalty,".txt"), col_names = FALSE, delim = " ")[,1:1000] %>% t() %>% as.data.frame()
colnames(perm_WD) <- sample_id
pvalue <- rep(NA,length(sample_id))
for (i in 1:length(pvalue)){
  pvalue[i] <- sum(perm_WD[,i]<min_dist_real$Min_dist[i], na.rm = TRUE)/1000
}

issync <- ifelse(pvalue<0.05, TRUE, FALSE)
table(issync)
synctable <- data.frame(sample_id,Sync_Gain=issync)
sinfo <- left_join(sinfo, synctable, by=c("Anom_Sample_ID"="sample_id"))
sample_class <- left_join(sample_class, sinfo[, c(1,3)])
sample_class$Label[is.na(sample_class$Sync_Gain)] <- "Uninformative"
sample_class$Label[sample_class$Sync_Gain==TRUE] <- "Synchronous"
sample_class$Label[sample_class$Sync_Gain==FALSE] <- "Asynchronous"
table(sample_class$WGD_Status, sample_class$Label)
write.csv(sample_class, "../in_data/sample_class_posterior_true.csv", row.names = FALSE)

# Draw final plots figure 1a ------
sample_class <- read_csv(paste0("../in_data/sample_class_posterior_",non_parsimony_penalty,".csv"))
df <- table(sample_class$WGD_Status, sample_class$Label, useNA = "ifany") %>% data.frame()
colnames(df) <- c("WGD_Status", "Sync_Status", "Value")
df$Freq <- NULL
df$Freq[df$WGD_Status == "TRUE"] <- df$Value[df$WGD_Status == "TRUE"]/sum(df$Value[df$WGD_Status == "TRUE"])
df$Freq[df$WGD_Status == "FALSE"] <- df$Value[df$WGD_Status == "FALSE"]/sum(df$Value[df$WGD_Status == "FALSE"])
table(sample_class$WGD_Status)
wgdtrue <- table(sample_class$WGD_Status)[2] %>% as.integer()
wgdfalse <- table(sample_class$WGD_Status)[1] %>% as.integer()
Fig.1a <- ggplot(df, aes(x=WGD_Status, y=Freq, fill = Sync_Status)) +
  geom_bar(position = "fill", stat = "identity", color=NA, width=0.9, linewidth = 0.5) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .02))) +
  labs(y= "Tumors")+
  scale_x_discrete(labels = c(paste0("Near-Diploid\n(n = ", wgdfalse, ")"), paste0("WGD\n(n = ", wgdtrue,")")), expand = expansion(add = .55)) +
  scale_fill_manual(name = NULL, labels = c("Non-punctuated", "Punctuated", "Uninformative"), values=c("#EFC000FF", "#0073C2FF", "#868686FF"))+
  geom_text(aes(label = paste0(round(Freq*100, 2),"%")), 
            # position = position_stack(vjust = 0.5), size = 7, fontface = "bold") +
            position = position_stack(vjust = 0.5), size = 7) +
  theme(panel.background = element_blank(), text = element_text(size = 24)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(vjust = +2),
        axis.text.y = element_blank(), axis.text.x = element_text(size = 20, vjust = -1, color = "black")) + 
  theme(axis.line = element_line(color="black", linewidth = 1), axis.ticks.y = element_blank(), axis.ticks = element_line(colour = "black", size = 1.2), axis.ticks.length.x = unit(5, "mm"))+
  theme(legend.position="right") + theme(legend.key.size = unit(1.2, "cm")) + 
  guides(fill=guide_legend(nrow=3, byrow=TRUE)) + 
  theme(text=element_text(family="HelveticaNeue")) +
  # theme(text = element_text(face = "bold"))+
  theme(plot.margin = margin(5,5,10,10))

# theme(legend.position = c(0.25, 0.1), legend.background = element_rect(fill=alpha("white", 0.4)))
Fig.1a
ggsave(paste0("Fig.1a_non_parsimony_penalty_",non_parsimony_penalty,".eps"), plot = Fig.1a, width = 7, height = 7)


# Fig 1a without uninformative samples ------
df1 <- table(sample_class$WGD_Status, sample_class$Label, useNA = "ifany") %>% data.frame()
df1 <- df1[1:4,]
colnames(df1) <- c("WGD_Status", "Sync_Status", "Value")
df1$Freq <- NULL
df1$Freq[df1$WGD_Status == "TRUE"] <- df1$Value[df1$WGD_Status == "TRUE"]/sum(df1$Value[df1$WGD_Status == "TRUE"])
df1$Freq[df1$WGD_Status == "FALSE"] <- df1$Value[df1$WGD_Status == "FALSE"]/sum(df1$Value[df1$WGD_Status == "FALSE"])
Fig.1a1 <- ggplot(df1, aes(x=WGD_Status, y=Freq, fill = Sync_Status)) +
  geom_bar(position = "fill", stat = "identity", color=NA, width=0.9, linewidth = 0.5) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .02))) +
  labs(y= "Tumors")+
  scale_x_discrete(labels = c(paste0("Near-Diploid\n(n = ", wgdfalse, ")"), paste0("WGD\n(n = ", wgdtrue,")")), expand = expansion(add = .55)) +
  scale_fill_manual(name = NULL, labels = c("Non-punctuated", "Punctuated", "Uninformative"), values=c("#EFC000FF", "#0073C2FF", "#868686FF"))+
  geom_text(aes(label = paste0(round(Freq*100, 2),"%")), 
            # position = position_stack(vjust = 0.5), size = 7, fontface = "bold") +
            position = position_stack(vjust = 0.5), size = 7) +
  theme(panel.background = element_blank(), text = element_text(size = 24)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(vjust = +2),
        axis.text.y = element_blank(), axis.text.x = element_text(size = 20, vjust = -1, color = "black")) + 
  theme(axis.line = element_line(color="black", linewidth = 1), axis.ticks.y = element_blank(), axis.ticks = element_line(colour = "black", size = 1.2), axis.ticks.length.x = unit(5, "mm"))+
  theme(legend.position="right") + theme(legend.key.size = unit(1.2, "cm")) + 
  guides(fill=guide_legend(nrow=3, byrow=TRUE)) + 
  theme(text=element_text(family="HelveticaNeue")) +
  # theme(text = element_text(face = "bold"))+
  theme(plot.margin = margin(5,5,10,10))

# theme(legend.position = c(0.25, 0.1), legend.background = element_rect(fill=alpha("white", 0.4)))
Fig.1a1
# extrafont::loadfonts()
ggsave(paste0("Fig.1a1_non_parsimony_penalty_",non_parsimony_penalty,".eps"), plot = Fig.1a1, width = 7, height = 7)


# Draw final plots figure 1c ------
WGD_Sync <- sample_class %>% filter(WGD_Status==TRUE, Label == "Synchronous")
WGD_Sync <- left_join(WGD_Sync, min_dist_real)
timing_WGD_sync <-timing %>% filter(WGD_Status == TRUE) %>% group_by(Anom_Sample_ID) %>% 
  mutate(WGD_Timing_avg=mean(WGD_Timing)) %>% slice(1) %>% select(Anom_Sample_ID, WGD_Timing_avg)
WGD_Sync <- left_join(WGD_Sync, timing_WGD_sync)
WGD_Sync$Relative_Sync <- ifelse(WGD_Sync$Timing_Sync>WGD_Sync$WGD_Timing_avg, "Post_WGD", "Pre_WGD")
wgdsynctable <- table(WGD_Sync$Relative_Sync) %>% as.data.frame()
wgdsynctable$frac <- wgdsynctable$Freq/sum(wgdsynctable$Freq)
wgdsynctable$Label <- "Sync"
write_csv(wgdsynctable, "pre_post_label.csv")

wgdsynctable <- read_csv("pre_post_label.csv")
Fig.1c <- ggplot(wgdsynctable, aes(x=Label, y=frac, fill = Var1)) +
  geom_bar(position = "fill", stat = "identity", color=NA, width=0.9, linewidth = 0.5) +
  geom_text(aes(label = paste0(round(frac*100, 2),"%  n=", Freq)), 
            # position = position_stack(vjust = 0.5), size = 7, fontface = "bold")+
            position = position_stack(vjust = 0.5), size = 7)+
  labs(y = "Synchronous Gains in WGD Tumors")+
  scale_x_discrete(labels = element_blank(), expand = expansion(add = .55)) +
  scale_fill_manual(name = NULL, labels = c("Post_WGD", "Pre_WGD"), values=c("#8d5ea1", "#0b8437"))+
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .02))) +
  theme(panel.background = element_blank(), text = element_text(size = 24)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(vjust = +2), axis.text = element_blank()) + 
  theme(axis.line = element_line(color="black", linewidth = 1), axis.ticks = element_blank()) + 
  theme(legend.position="right") + theme(legend.key.size = unit(1.2, "cm")) + 
  # theme(text = element_text(face = "bold"))+
  theme(text=element_text(family="HelveticaNeue")) +
  theme(plot.margin = margin(10,30,20,20))
Fig.1c
ggsave(paste0("Fig.1c_non_parsimony_penalty_",non_parsimony_penalty,".eps"), plot = Fig.1c, width = 6.5, height = 7)

# not used this time -----
(Fig.1b <- PlotTiming_posteriors(timing_obj = timing, sampleid = "DO50842"))
ggsave("Fig.1b.pdf", plot = Fig.1b, width = 7, height = 3.5)

(Fig.1d <- PlotTiming_posteriors(timing_obj = timing, sampleid = "HMF001677A"))
ggsave("Fig.1d.pdf", plot = Fig.1d, width = 7, height = 3.5)

for (i in 1:nrow(WGD_Sync)) {
  if(WGD_Sync$Relative_Sync[i]=="Pre_WGD"){
    Figtmp <- PlotTiming_posteriors(timing_obj = timing, sampleid = WGD_Sync$Anom_Sample_ID[i])
    ggsave(paste0("sync plots pre/",WGD_Sync$Anom_Sample_ID[i], ".pdf"), plot = Figtmp, width = 8.5, height = 3.5)
  }
}




# assign gains to the sync event by threshold -----
sample_class <- read_csv(paste0("../in_data/sample_class_posterior_",non_parsimony_penalty,".csv"))
timing <- readRDS(paste0("timing_keep_non_parsimony_penalty_",non_parsimony_penalty,"_04092024.rds"))
min_dist_real <- read_delim(paste0("../in_data/min_dist_median_non_parsimony_penalty_",non_parsimony_penalty,".txt"),delim = "\t", col_names = FALSE)
colnames(min_dist_real) <- c("Anom_Sample_ID", "Min_dist", "Timing_Sync")
sample_id <- min_dist_real$Anom_Sample_ID
min_dist_real <- left_join(min_dist_real, sample_class)
# p1 <- hist(min_dist_real$Min_dist[min_dist_real$Sync_Gain==FALSE], breaks = 100)
# p2 <- hist(min_dist_real$Min_dist[min_dist_real$Sync_Gain==TRUE], breaks = 100)
# plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,0.5))
# plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,0.5), add=T) 
# 
# sampleid = "CPCT02130085T"
# d0 <- min_dist_real$Min_dist[min_dist_real$Anom_Sample_ID==sampleid]
# t0 <- min_dist_real$Timing_Sync[min_dist_real$Anom_Sample_ID==sampleid]
# p1 <- PlotTiming_posteriors(timing_obj = timing, sampleid = sampleid) + geom_abline(slope = 0, intercept = t0)
# ts <- timing %>% filter(Anom_Sample_ID %in% sampleid) %>% droplevels() %>% group_by(Segment_ID, Posterior_Sample_Index, Route) %>%
#   mutate(groupInd = with_order(order_by = Gain_Timing, fun = row_number, x = Gain_Timing)) %>% 
#   ungroup() %>% group_by(Segment_ID, groupInd) %>% mutate(timing_median = median(Gain_Timing))
# ts$groupInd <- as.factor(ts$groupInd)
# # unique(ts %>% select(Segment_ID, timing_median))
# ts <- ts %>% mutate(belong_to_sync = (timing_median<=(t0+0.1) & timing_median>=(t0-0.1)))
# p2 <- PlotTiming_posteriors(timing_obj = ts[ts$belong_to_sync==TRUE,], sampleid = sampleid) + geom_abline(slope = 0, intercept = t0)
# p1+p2+plot_layout(ncol = 1)



annotate_gains <- function(timing_obj=timing, sampleid = "CPCT02130085T"){
  d0 <- min_dist_real$Min_dist[min_dist_real$Anom_Sample_ID==sampleid]
  t0 <- min_dist_real$Timing_Sync[min_dist_real$Anom_Sample_ID==sampleid]
  ts <- timing %>% filter(Anom_Sample_ID %in% sampleid) %>% droplevels() %>% group_by(Segment_ID, Posterior_Sample_Index, Route) %>%
    mutate(groupInd = with_order(order_by = Gain_Timing, fun = row_number, x = Gain_Timing)) %>% 
    ungroup() %>% group_by(Segment_ID, groupInd) %>% mutate(timing_median = median(Gain_Timing))
  ts$groupInd <- as.factor(ts$groupInd)
  ts <- ts %>% mutate(belong_to_sync = (timing_median<=(t0+d0) & timing_median>=(t0-d0)))
  ts
}
sync_sample_id <- min_dist_real$Anom_Sample_ID[min_dist_real$Sync_Gain]

library(data.table)
list_of_dfs <- vector("list", length(sync_sample_id))
for (i in seq_along(sync_sample_id)) {
  list_of_dfs[[i]] <- annotate_gains(timing_obj = timing, sampleid = sync_sample_id[i])
}
timing_annotated <- rbindlist(list_of_dfs)
saveRDS(timing_annotated, file=paste0("../in_data/timing_sync_non_parsimony_penalty_",non_parsimony_penalty,"_annotated_04092024.rds"))
