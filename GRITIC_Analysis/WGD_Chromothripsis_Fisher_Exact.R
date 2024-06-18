# Load Packages
library(dplyr)
library(stats)
library(ggmosaic)

# Load data

#raw PCAWG data
PCAWG_IDs <- read.table("WGD_Chromothripsis_Fisher/PCAWG_ID.txt")$V1
PCAWG_Chromothripsis <- read.table("WGD_Chromothripsis_Fisher/PCAWG_chromothripsisOverlap.txt", sep = "\t", header = T)

for(non_parsimony_penalty in c("True", "False")){

  Synchronous_Gain_Calls <- read.csv(paste0("../in_data/sample_class_posterior_",non_parsimony_penalty,".csv"), sep = ",", header = T)
  
  # Pull informative samples in PCAWG from synchronous gain calls
  Synchronous_Gain_Calls_Filtered <- Synchronous_Gain_Calls %>% filter(Anom_Sample_ID %in% PCAWG_IDs & Label != "Uninformative")
  Synchronous_Gain_Calls_Filtered_Boolean <- Synchronous_Gain_Calls_Filtered %>% group_by(Anom_Sample_ID) %>% summarize(WGD = WGD_Status, Synchronous = ifelse(Label == "Synchronous", TRUE, FALSE))
  
  # Use informative PCAWG synchronous gain calls to filter matching PCAWG chromothripsis samples
  PCAWG_Chromothripsis_Filtered <- PCAWG_Chromothripsis %>% filter(samplename %in% Synchronous_Gain_Calls_Filtered$Anom_Sample_ID)
  PCAWG_Chromothripsis_Filtered_Boolean <- PCAWG_Chromothripsis_Filtered %>% group_by(samplename) %>% summarize(Chromothripsis = any(FinalCalls == "Chromothripsis"))
  
  # Match data
  ## Omit synchronous gain data not present in chromothripsis calls
  # All_Calls_Boolean <- cbind(
  #   Synchronous_Gain_Calls_Filtered_Boolean, 
  #   PCAWG_Chromothripsis_Filtered_Boolean[match(Synchronous_Gain_Calls_Filtered_Boolean$Anom_Sample_ID, PCAWG_Chromothripsis_Filtered_Boolean$samplename),"Chromothripsis"]
  # ) %>% na.omit()
  
  ## Include synchronous gain data not present in chromothripsis calls and set chromothripsis status to FALSE
  All_Calls_Boolean <- cbind(
    Synchronous_Gain_Calls_Filtered_Boolean,
    PCAWG_Chromothripsis_Filtered_Boolean[match(Synchronous_Gain_Calls_Filtered_Boolean$Anom_Sample_ID, PCAWG_Chromothripsis_Filtered_Boolean$samplename),"Chromothripsis"]
  )
  All_Calls_Boolean[is.na(All_Calls_Boolean$Chromothripsis), 4] <- FALSE
  
  # Create matrices
  WGD_CHROMO <- matrix(c(0,0,0,0), ncol = 2, nrow = 2, dimnames = list(c("WGD_FALSE", "WGD_TRUE"), c("CHROMO_FALSE", "CHROMO_TRUE")))
  WGD_SYNC <- matrix(c(0,0,0,0), ncol = 2, nrow = 2, dimnames = list(c("WGD_FALSE", "WGD_TRUE"), c("SYNC_FALSE", "SYNC_TRUE")))
  CHROMO_SYNC <- matrix(c(0,0,0,0), ncol = 2, nrow = 2, dimnames = list(c("CHROMO_FALSE", "CHROMO_TRUE"), c("SYNC_FALSE", "SYNC_TRUE")))
  
  # Fill matrices
  for(i in 1:nrow(All_Calls_Boolean)){
    Row <- All_Calls_Boolean[i,]
    WGD <- Row$WGD
    SYNC <- Row$Synchronous
    CHROMO <- Row$Chromothripsis
    
    WGD_CHROMO[WGD+1,CHROMO+1] <- WGD_CHROMO[WGD+1,CHROMO+1] + 1
    WGD_SYNC[WGD+1,SYNC+1] <- WGD_SYNC[WGD+1,SYNC+1] + 1
    CHROMO_SYNC[CHROMO+1,SYNC+1] <- CHROMO_SYNC[CHROMO+1,SYNC+1] + 1
  }
  
  # 2x2 Fisher Exact Tests w/ Plots + Tables
  png(file = paste0("WGD_Chromothripsis_Fisher/WGD_Chromothripsis_Fisher_",non_parsimony_penalty,".png"))
  WGD_CHROMO_RES <- fisher.test(WGD_CHROMO)
  WGD_CHROMO_TBL <- data.frame(p = signif(WGD_CHROMO_RES$p, 2), OR = signif(WGD_CHROMO_RES$estimate, 2), CI.lower = signif(WGD_CHROMO_RES$conf.int[1], 2), CI.upper = signif(WGD_CHROMO_RES$conf.int[2], 2))
  rownames(WGD_CHROMO_TBL) <- NULL
  par(mfrow = c(2,1), mar = c(2,2,2,2))
  mosaicplot(WGD_CHROMO, color = T)
  gridExtra::grid.table(WGD_CHROMO_TBL)
  dev.off()
  ###
  png(file = paste0("WGD_Chromothripsis_Fisher/WGD_SynchronousGains_Fisher_",non_parsimony_penalty,".png"))
  WGD_SYNC_RES <- fisher.test(WGD_SYNC)
  WGD_SYNC_TBL <- data.frame(p = signif(WGD_SYNC_RES$p, 2), OR = signif(WGD_SYNC_RES$estimate, 2), CI.lower = signif(WGD_SYNC_RES$conf.int[1], 2), CI.upper = signif(WGD_SYNC_RES$conf.int[2], 2))
  rownames(WGD_SYNC_TBL) <- NULL
  par(mfrow = c(2,1), mar = c(2,2,2,2))
  mosaicplot(WGD_SYNC, color = T)
  gridExtra::grid.table(WGD_SYNC_TBL)
  dev.off()
  ###
  png(file = paste0("WGD_Chromothripsis_Fisher/Chromothripsis_SynchronousGains_Fisher_",non_parsimony_penalty,".png"))
  CHROMO_SYNC_RES <- fisher.test(CHROMO_SYNC)
  CHROMO_SYNC_TBL <- data.frame(p = signif(CHROMO_SYNC_RES$p, 2), OR = signif(CHROMO_SYNC_RES$estimate, 2), CI.lower = signif(CHROMO_SYNC_RES$conf.int[1], 2), CI.upper = signif(CHROMO_SYNC_RES$conf.int[2], 2))
  rownames(CHROMO_SYNC_TBL) <- NULL
  par(mfrow = c(2,1), mar = c(2,2,2,2))
  mosaicplot(CHROMO_SYNC, color = T)
  gridExtra::grid.table(CHROMO_SYNC_TBL)
  dev.off()
  
  
  ggplot(data = All_Calls_Boolean) +
    ggtitle(paste("non_parsimony_penalty:", non_parsimony_penalty)) +
    geom_mosaic(aes(x = product(Chromothripsis, Synchronous), fill = Synchronous), alpha = 0.6) + 
    geom_mosaic_text(aes(x = product(Chromothripsis, Synchronous), fill = Synchronous, label = after_stat(.wt))) +
    scale_fill_brewer(palette = "Dark2") + 
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 12),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0),
      axis.text.x = element_text(hjust = 0.5, vjust = 1),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "Punctuated")
  
  ggsave(paste0("WGD_Chromothripsis_Fisher/Chromothripsis_SynchronousGains_Mosaic_", non_parsimony_penalty, ".pdf"))
  
  ggplot(data = All_Calls_Boolean) +
    ggtitle(paste("non_parsimony_penalty:", non_parsimony_penalty)) +
    geom_mosaic(aes(x = product(WGD, Synchronous), fill = Synchronous), alpha = 0.6) + 
    geom_mosaic_text(aes(x = product(WGD, Synchronous), fill = Synchronous, label = after_stat(.wt))) +
    scale_fill_brewer(palette = "Dark2") + 
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 12),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0),
      axis.text.x = element_text(hjust = 0.5, vjust = 1),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "Punctuated")
  
  ggsave(paste0("WGD_Chromothripsis_Fisher/WGD_SynchronousGains_Mosaic_", non_parsimony_penalty, ".pdf"))
  
  ggplot(data = All_Calls_Boolean) +
    ggtitle(paste("non_parsimony_penalty:", non_parsimony_penalty)) +
    geom_mosaic(aes(x = product(Chromothripsis, WGD), fill = WGD), alpha = 0.6) + 
    geom_mosaic_text(aes(x = product(Chromothripsis, WGD), fill = WGD, label = after_stat(.wt))) +
    scale_fill_brewer(palette = "Dark2") + 
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 12),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0),
      axis.text.x = element_text(hjust = 0.5, vjust = 1),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "WGD")
  
  ggsave(paste0("WGD_Chromothripsis_Fisher/Chromothripsis_WGD_Mosaic_", non_parsimony_penalty, ".pdf"))
  
}
