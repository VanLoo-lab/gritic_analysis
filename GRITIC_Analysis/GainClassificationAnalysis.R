
GenomeSizes <- read.table("../resources/GenomeSizes.txt", header = T, sep = "\t")
source("GainClassifier.R")

## Packages
library(dplyr)
library(ggplot2)

## Gain Table w/ Synchronicity Annotation (Need to add Gain_Index)
#change to prior false to run without non-parsimony penalty
SyncAnnotatedGains <- readRDS("in_data/timing_sync_prior_true_annotated_04092024.rds")
SyncAnnotatedGains <- SyncAnnotatedGains %>% group_by(Anom_Sample_ID,Segment_ID,Posterior_Sample_Index) %>% mutate(Gain_Index=rank(Gain_Timing))

### TESTING
SVDirectory = c("PCAWG_PATH","TCGA_PATH")
SVSuffix = c(".pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz", ".pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz")
df <- EventClassification(SyncAnnotatedGains, SegmentMergeRule = "none", SpaceAdjacencyThreshold = 60000, TimeAdjacencyThreshold = 0.025, SVDirectory = SVDirectory, SVSuffix = SVSuffix)

#MergePlot(df, "0009b464-b376-4fbc-8a56-da538269a02f")
#MergePlot(df, unique(df[[2]]$Anom_Sample_ID)[2])

ClassifiedSegments <- df[[4]]

GenomeSizes$Chr <- factor(GenomeSizes$Chr, levels = c(1:22, "X", "Y"))
ClassifiedSegments$Chr <- factor(ClassifiedSegments$Chr, levels = c(1:22, "X", "Y"))

ClassifiedSegments %>% 
  filter(Anom_Sample_ID == "10bb1a92-901e-4a14-80f4-5e88f997754b") %>%
  ggplot() + 
  geom_segment(data = GenomeSizes, aes(x = Start, xend = End, y = 0, yend = 0), color = "transparent") + 
  geom_rect(data = GenomeSizes, aes(xmin = Centromere_Start, xmax = Centromere_End, ymin = -Inf, ymax = Inf), fill = "grey80", color = "transparent") + 
  geom_rect(aes(xmin = Segment_Start_After_Merge, xmax = Segment_End_After_Merge, ymin = Timing_After_Merge - 0.005, ymax = Timing_After_Merge + 0.005, fill = factor(Topology_Class)), linewidth = 0.5) + 
  ggh4x::facet_grid2(Anom_Sample_ID~Chr, space = "free_x", scale = "free_x", switch = "x") + 
  coord_cartesian(clip = "off") + 
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) + 
  scale_x_continuous(expand = c(0,0), breaks = c(0,1)) + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(
    y = "Gain Timing",
    x = "Genome Position",
    fill = "Gain"
  ) +
  theme_minimal() + 
  theme(
    panel.spacing.x = unit(0,units = "lines"),
    panel.grid = element_blank(),
    panel.background = element_rect(linewidth = 0.25, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.25),
    plot.title = element_text(face = "bold"),
    strip.placement = "inside"
  )

df[[4]]$Length_Class <- factor(ClassifiedSegments$Length_Class, levels = c("<100KB", "100KB-500KB", "500KB-1MB", "1MB-5MB", "5MB-10MB", "10MB-50MB", ">50MB"))
ggplot(ClassifiedSegments) + geom_bar(aes(x = Length_Class, y = ..count.., fill = Form_Class))

#ggsave("10bb1a92-901e-4a14-80f4-5e88f997754b_Classified_SV.pdf", width = 16, height = 8)

#write.table(df[[1]], file = "20240303-ARL_SegmentationAfterMerge.txt", sep = "\t", col.names = T, row.names = F, quote = F)

df[[4]] %>%
  select(Anom_Sample_ID, SVFileExists) %>%
  unique() %>%
  ggplot() + 
  geom_bar(aes(x=1, fill = SVFileExists)) + 
  coord_cartesian(expand = F, clip = "off") + 
  theme_classic() + 
  theme(
    axis.text.x = element_blank()
  ) + 
  labs(
    y = "Has SV File",
    x = ""
  )

#ggsave("SampleWithSVFile.pdf", width = 2.2, height = 4.5)



df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Majority_Synchronous) %>%
  mutate(Total = n()) %>%
  group_by(Majority_Synchronous, Topology_Class) %>%
  summarize(
    Perc = n() / Total
  ) %>% unique() -> p1.dat
p1.dat

df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  ggplot() + 
  geom_bar(aes(x=Majority_Synchronous, y = (..count..)/sum(..count..), fill = Topology_Class)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  #coord_cartesian(expand = F, clip = "off") + 
  theme_classic() + 
  labs(
    y = "Fraction of Gains",
    x = "Synchronous"
  )

#ggsave("Sync_Topo_Type.pdf", width = 3, height = 5)


df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Majority_Synchronous) %>%
  mutate(Total = n()) %>%
  group_by(Majority_Synchronous, Form_Class) %>%
  summarize(
    Perc = n() / Total
  ) %>% unique() -> p2.dat
p2.dat

df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  ggplot() + 
  geom_bar(aes(x=Majority_Synchronous, y = (..count..)/sum(..count..), fill = Form_Class)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_fill_brewer(palette = "Dark2") + 
  #coord_cartesian(expand = F, clip = "off") + 
  theme_classic() + 
  labs(
    y = "Fraction of Gains",
    x = "Synchronous"
  )

#ggsave("Sync_Seg_Type.pdf", width = 4, height = 5)

# Contingency Tables
## Synchronicity - Length
df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Majority_Synchronous, Length_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Majority_Synchronous, Length_Class, fill = list(Count = 0)) %>%
  mutate(
    Length_Class = factor(Length_Class, levels = c("<100KB", "100KB-500KB", "500KB-1MB", "1MB-5MB", "5MB-10MB", "10MB-50MB", ">50MB")),
    Synchronicity = ifelse(Majority_Synchronous == TRUE, "Punctuated", "Not Punctuated")
  ) %>%
  tidyr::pivot_wider(names_from = Length_Class, values_from = Count) %>% as.data.frame() -> p1.mat

rownames(p1.mat) <- p1.mat$Synchronicity
p1.mat$Synchronicity <- NULL
p1.mat$Majority_Synchronous <- NULL
p1.mat <- as.matrix(p1.mat)

fisher.test(p1.mat, simulate.p.value = T)
x2.1 <- chisq.test(p1.mat)

df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Majority_Synchronous, Length_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Majority_Synchronous, Length_Class, fill = list(Count = 0)) %>%
  mutate(
    Length_Class = factor(Length_Class, levels = c("<100KB", "100KB-500KB", "500KB-1MB", "1MB-5MB", "5MB-10MB", "10MB-50MB", ">50MB")),
    Synchronicity = ifelse(Majority_Synchronous == TRUE, "Punctuated", "Not Punctuated")
  ) %>%
  ggplot() + 
  ggtitle(paste("X2 =", round(x2.1$statistic, 4), ifelse(x2.1$p.value < 0.0001, ", p < 0.0001", paste(", p =", round(x2.1$p.value, 4))))) +
  geom_tile(aes(x = Synchronicity, y = Length_Class, fill = log10(Count + 0.1)), color = "black") +
  geom_label(aes(x = Synchronicity, y = Length_Class, label = Count), color = "black", label.r = unit(0, "lines"), label.size = unit(0, "lines"), size = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  coord_cartesian(expand = F, clip = "off") + 
  theme_classic() + 
  theme(
    plot.title = element_text(size = 7),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = "transparent"),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "none"
  ) -> p1

## Synchronicity - Form
df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Majority_Synchronous, Form_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Majority_Synchronous, Form_Class, fill = list(Count = 0)) %>%
  mutate(
    Form_Class = factor(Form_Class, levels = rev(c("Whole Chromosome", "Supra-Centromeric Arm", "Whole Arm", "Sub-Arm Segment (Telomere-Bounded)", "Sub-Arm Segment (Centromere-Bounded)", "Sub-Arm Segment (Interstitial)"))),
    Synchronicity = ifelse(Majority_Synchronous == TRUE, "Punctuated", "Not Punctuated")
  ) %>%
  tidyr::pivot_wider(names_from = Form_Class, values_from = Count) %>% as.data.frame() -> p2.mat

rownames(p2.mat) <- p2.mat$Synchronicity
p2.mat$Synchronicity <- NULL
p2.mat$Majority_Synchronous <- NULL
p2.mat <- as.matrix(p2.mat)

fisher.test(p2.mat, simulate.p.value = T)
x2.2 <- chisq.test(p2.mat)

df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Majority_Synchronous, Form_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Majority_Synchronous, Form_Class, fill = list(Count = 0)) %>%
  mutate(
    Form_Class = factor(Form_Class, levels = rev(c("Whole Chromosome", "Supra-Centromeric Arm", "Whole Arm", "Sub-Arm Segment (Telomere-Bounded)", "Sub-Arm Segment (Centromere-Bounded)", "Sub-Arm Segment (Interstitial)"))),
    Synchronicity = ifelse(Majority_Synchronous == TRUE, "Punctuated", "Not Punctuated")
  ) %>%
  ggplot() + 
  ggtitle(paste("X2 =", round(x2.2$statistic, 4), ifelse(x2.2$p.value < 0.0001, ", p < 0.0001", paste(", p =", round(x2.2$p.value, 4))))) +
  geom_tile(aes(x = Synchronicity, y = Form_Class, fill = log10(Count + 0.1)), color = "black") +
  geom_label(aes(x = Synchronicity, y = Form_Class, label = Count), color = "black", label.r = unit(0, "lines"), label.size = unit(0, "lines"), size = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  coord_cartesian(expand = F, clip = "off") + 
  theme_classic() + 
  theme(
    plot.title = element_text(size = 7),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = "transparent"),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "none"
  ) -> p2

## Topology - Form
df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Topology_Class, Form_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Topology_Class, Form_Class, fill = list(Count = 0)) %>%
  mutate(
    Topology_Class = factor(Topology_Class, levels = c("Intra-Chromosomal", "Inter-Chromosomal", "None")),
    Form_Class = factor(Form_Class, levels = rev(c("Whole Chromosome", "Supra-Centromeric Arm", "Whole Arm", "Sub-Arm Segment (Telomere-Bounded)", "Sub-Arm Segment (Centromere-Bounded)", "Sub-Arm Segment (Interstitial)")))
  ) %>%
  tidyr::pivot_wider(names_from = Form_Class, values_from = Count) %>% as.data.frame() -> p3.mat

rownames(p3.mat) <- p3.mat$Topology_Class
p3.mat$Topology_Class <- NULL
p3.mat <- as.matrix(p3.mat)

fisher.test(p3.mat, simulate.p.value = T)
x2.3 <- chisq.test(p3.mat)

df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Topology_Class, Form_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Topology_Class, Form_Class, fill = list(Count = 0)) %>%
  mutate(
    Topology_Class = factor(Topology_Class, levels = c("Intra-Chromosomal", "Inter-Chromosomal", "None")),
    Form_Class = factor(Form_Class, levels = rev(c("Whole Chromosome", "Supra-Centromeric Arm", "Whole Arm", "Sub-Arm Segment (Telomere-Bounded)", "Sub-Arm Segment (Centromere-Bounded)", "Sub-Arm Segment (Interstitial)")))
  ) %>%
  ggplot() + 
  ggtitle(paste("X2 =", round(x2.3$statistic, 4), ifelse(x2.3$p.value < 0.0001, ", p < 0.0001", paste(", p =", round(x2.3$p.value, 4))))) +
  geom_tile(aes(x = Topology_Class, y = Form_Class, fill = log10(Count + 0.1)), color = "black") +
  geom_label(aes(x = Topology_Class, y = Form_Class, label = Count), color = "black", label.r = unit(0, "lines"), label.size = unit(0, "lines"), size = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  coord_cartesian(expand = F, clip = "off") + 
  theme_classic() + 
  theme(
    plot.title = element_text(size = 7),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = "transparent"),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none"
  ) -> p3

#ggsave(plot = p3, "TopologyCrosstab.pdf", width = 4, height = 5)

## Length - Form (! important statistical modification)
df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Length_Class, Form_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Length_Class, Form_Class, fill = list(Count = 0)) %>%
  mutate(
    Length_Class = factor(Length_Class, levels = c("<100KB", "100KB-500KB", "500KB-1MB", "1MB-5MB", "5MB-10MB", "10MB-50MB", ">50MB")),
    Form_Class = factor(Form_Class, levels = rev(c("Whole Chromosome", "Supra-Centromeric Arm", "Whole Arm", "Sub-Arm Segment (Telomere-Bounded)", "Sub-Arm Segment (Centromere-Bounded)", "Sub-Arm Segment (Interstitial)")))
  ) %>%
  tidyr::pivot_wider(names_from = Form_Class, values_from = Count) %>% as.data.frame() -> p4.mat

rownames(p4.mat) <- p4.mat$Length_Class
p4.mat$Length_Class <- NULL
p4.mat[,c("Supra-Centromeric Arm", "Whole Arm", "Whole Chromosome")] <- NULL # Removed as these cause dependence in the data set. This is denoted in the figure as a (!)
p4.mat <- as.matrix(p4.mat)

fisher.test(p4.mat, simulate.p.value = T)
x2.4 <- chisq.test(p4.mat, simulate.p.value = T)

df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Length_Class, Form_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Length_Class, Form_Class, fill = list(Count = 0)) %>%
  mutate(
    Length_Class = factor(Length_Class, levels = c("<100KB", "100KB-500KB", "500KB-1MB", "1MB-5MB", "5MB-10MB", "10MB-50MB", ">50MB")),
    Form_Class = factor(Form_Class, levels = rev(c("Whole Chromosome", "Supra-Centromeric Arm", "Whole Arm", "Sub-Arm Segment (Telomere-Bounded)", "Sub-Arm Segment (Centromere-Bounded)", "Sub-Arm Segment (Interstitial)")))
  ) %>%
  ggplot() + 
  ggtitle(paste("X2 =", round(x2.4$statistic, 4), ifelse(x2.4$p.value < 0.0001, ", p < 0.0001", paste(", p =", round(x2.4$p.value, 4))), "(!)")) +
  geom_tile(aes(x = Length_Class, y = Form_Class, fill = log10(Count + 0.1)), color = "black") +
  geom_label(aes(x = Length_Class, y = Form_Class, label = Count), color = "black", label.r = unit(0, "lines"), label.size = unit(0, "lines"), size = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  coord_cartesian(expand = F, clip = "off") + 
  theme_classic() + 
  theme(
    plot.title = element_text(size = 7),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = "transparent"),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none"
  ) -> p4

## Topology - Synchronicity
df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Topology_Class, Majority_Synchronous) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Topology_Class, Majority_Synchronous, fill = list(Count = 0)) %>%
  mutate(
    Topology_Class = factor(Topology_Class, levels = c("Intra-Chromosomal", "Inter-Chromosomal", "None")),
    Synchronicity = ifelse(Majority_Synchronous == TRUE, "Punctuated", "Not Punctuated")
  ) %>%
  tidyr::pivot_wider(names_from = Topology_Class, values_from = Count) %>% as.data.frame() -> p5.mat

rownames(p5.mat) <- p5.mat$Synchronicity
p5.mat$Synchronicity <- NULL
p5.mat$Majority_Synchronous <- NULL
p5.mat <- as.matrix(p5.mat)


fisher.test(p5.mat, simulate.p.value = T)
x2.5 <- chisq.test(p5.mat, simulate.p.value = T)

df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Majority_Synchronous, Topology_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Majority_Synchronous, Topology_Class, fill = list(Count = 0)) %>%
  mutate(
    Topology_Class = factor(Topology_Class, levels = c("Intra-Chromosomal", "Inter-Chromosomal", "None")),
    Synchronicity = ifelse(Majority_Synchronous == TRUE, "Punctuated", "Not Punctuated")
  ) %>%
  ggplot() + 
  ggtitle(paste("X2 =", round(x2.5$statistic, 4), ifelse(x2.5$p.value < 0.0001, ", p < 0.0001", paste(", p =", round(x2.5$p.value, 4))))) +
  geom_tile(aes(x = Synchronicity, y = Topology_Class, fill = log10(Count + 0.1)), color = "black") +
  geom_label(aes(x = Synchronicity, y = Topology_Class, label = Count), color = "black", label.r = unit(0, "lines"), label.size = unit(0, "lines"), size = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  coord_cartesian(expand = F, clip = "off") + 
  theme_classic() + 
  theme(
    plot.title = element_text(size = 7),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = "transparent"),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none"
  ) -> p5

## Length - Topology
df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Length_Class, Topology_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Length_Class, Topology_Class, fill = list(Count = 0)) %>%
  mutate(
    Length_Class = factor(Length_Class, levels = c("<100KB", "100KB-500KB", "500KB-1MB", "1MB-5MB", "5MB-10MB", "10MB-50MB", ">50MB")),
    Topology_Class = factor(Topology_Class, levels = c("Intra-Chromosomal", "Inter-Chromosomal", "None")),
  ) %>%
  tidyr::pivot_wider(names_from = Length_Class, values_from = Count) %>% as.data.frame() -> p6.mat

rownames(p6.mat) <- p6.mat$Topology_Class
p6.mat$Topology_Class <- NULL
p6.mat <- as.matrix(p6.mat)


fisher.test(p6.mat, simulate.p.value = T)
x2.6 <- chisq.test(p6.mat, simulate.p.value = T)

df[[4]] %>%
  filter(SVFileExists == TRUE) %>%
  group_by(Length_Class, Topology_Class) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  tidyr::complete(Length_Class, Topology_Class, fill = list(Count = 0)) %>%
  mutate(
    Length_Class = factor(Length_Class, levels = c("<100KB", "100KB-500KB", "500KB-1MB", "1MB-5MB", "5MB-10MB", "10MB-50MB", ">50MB")),
    Topology_Class = factor(Topology_Class, levels = c("Intra-Chromosomal", "Inter-Chromosomal", "None")),
  ) %>%
  ggplot() + 
  ggtitle(paste("X2 =", round(x2.5$statistic, 4), ifelse(x2.5$p.value < 0.0001, ", p < 0.0001", paste(", p =", round(x2.5$p.value, 4))))) +
  geom_tile(aes(x = Topology_Class, y = Length_Class, fill = log10(Count+0.1)), color = "black") +
  geom_label(aes(x = Topology_Class, y = Length_Class, label = Count), color = "black", label.r = unit(0, "lines"), label.size = unit(0, "lines"), size = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  coord_cartesian(expand = F, clip = "off") + 
  theme_classic() + 
  theme(
    plot.title = element_text(size = 7),
    panel.border = element_rect(color = "black", linewidth = 0.25, fill = "transparent"),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none"
  ) -> p6

#ggsave(plot = p6, "LengthCrosstab.pdf", width = 4, height = 5)

pdf("Crosstabs_PriorTrue_NoMerge.pdf", width = 6.5*1.5, height = 9*1.25)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, align = "hv", labels = "AUTO", nrow = 3)
dev.off()

