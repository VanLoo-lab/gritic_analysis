# Genome Size/Length Information
#GenomeSizes <- read.table("GenomeSizes.txt", header = T, sep = "\t")

# Define Function
EventClassification <- function(GainSet, SegmentMergeRule = "adj_idx_sync", SpaceAdjacencyThreshold = 100, TimeAdjacencyThreshold = 0.05, SVDirectory = NULL, SVSuffix = NULL){
  # Pre-processing data set
  print("Reducing gain space")
  GainSetReduced <- GainSet %>% 
    dplyr::group_by(Anom_Sample_ID,Segment_ID,Posterior_Sample_Index) %>% 
    mutate(Gain_Index=rank(Gain_Timing)) %>%
    ungroup() %>%
    dplyr::select(Anom_Sample_ID, Segment_ID, Chr, Segment_Start, Segment_End, Gain_Index, timing_median, belong_to_sync) %>%
    unique() %>%
    arrange(Anom_Sample_ID, Chr, Segment_Start, Gain_Index)
  
  # Splitting and merging by sample, chromosome, and gain index to merge adjacent segments which may be split by segmental gains occurring later.
  print(paste("Grouping segments by merge rule:", SegmentMergeRule))
  
  if(SegmentMergeRule == "none"){
    print("Skipping Merge")
    Grouped <- GainSetReduced %>%
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index) %>%
      dplyr::mutate(
        Diff_Seg_Left = Segment_Start - lag(Segment_End, default = Inf),
        Test = Diff_Seg_Left > SpaceAdjacencyThreshold,
        Majority_Synchronous = sum(belong_to_sync) / n() >= 0.5,
        Segment_Start_After_Merge = Segment_Start,
        Segment_End_After_Merge = Segment_End,
        Segment_ID_After_Merge = Segment_ID,
        Timing_After_Merge = timing_median,
        Grouped_Timing_Before_Merge = mean(timing_median)
      )
    
    Merged <- Grouped %>% arrange(Anom_Sample_ID, Chr, Segment_Start_After_Merge, Gain_Index)
  }
  
  if(SegmentMergeRule == "adj_idx_sync"){
    print("Merging Segments")
    Grouped <- GainSetReduced %>%
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index) %>%
      dplyr::mutate(
        Diff_Seg_Left = Segment_Start - lag(Segment_End, default = Inf),
        Test = Diff_Seg_Left > SpaceAdjacencyThreshold,
        Merge_Group = cumsum(Diff_Seg_Left > SpaceAdjacencyThreshold)
      ) %>%
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index, Merge_Group) %>%
      mutate(
        Majority_Synchronous = sum(belong_to_sync) / n() >= 0.5,
        Segment_Start_After_Merge = min(Segment_Start),
        Segment_End_After_Merge = max(Segment_End),
        Segment_ID_After_Merge = paste(Chr, Segment_Start_After_Merge, Segment_End_After_Merge, sep = "-"),
        Grouped_Timing_Before_Merge = mean(timing_median)
      )
    
    Merged_Segments <- Grouped %>% 
      filter(Majority_Synchronous == TRUE) %>% 
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index, Merge_Group) %>%
      dplyr::mutate(
        Segment_Start_After_Merge = min(Segment_Start),
        Segment_End_After_Merge = max(Segment_End),
        Segment_ID_After_Merge = paste(Chr, Segment_Start_After_Merge, Segment_End_After_Merge, sep = "-"),
        Timing_After_Merge = mean(timing_median)
      ) %>% 
      dplyr::ungroup() %>%
      dplyr::select(Anom_Sample_ID, Segment_ID_After_Merge, Chr, Gain_Index, Segment_Start_After_Merge, Segment_End_After_Merge, Timing_After_Merge, Majority_Synchronous) %>%
      unique()
    
    NonMerged_Segments <- Grouped %>%
      filter(Majority_Synchronous == FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Segment_Start_After_Merge = Segment_Start,
        Segment_End_After_Merge = Segment_End,
        Segment_ID_After_Merge = Segment_ID,
        Timing_After_Merge = timing_median
      ) %>%
      dplyr::select(Anom_Sample_ID, Segment_ID_After_Merge, Chr, Gain_Index, Segment_Start_After_Merge, Segment_End_After_Merge, Timing_After_Merge, Majority_Synchronous)
    
    Merged <- rbind(Merged_Segments, NonMerged_Segments) %>% arrange(Anom_Sample_ID, Chr, Segment_Start_After_Merge, Gain_Index)
  }
  
  if(SegmentMergeRule == "adj_idx_time"){
    print("Merging Segments")
    Grouped <- GainSetReduced %>%
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index) %>%
      dplyr::mutate(
        Diff_Seg_Left = Segment_Start - lag(Segment_End, default = Inf),
        Diff_Time_Left = abs(timing_median - lag(timing_median, default = Inf)),
        Test = Diff_Seg_Left > SpaceAdjacencyThreshold | (Diff_Seg_Left <= SpaceAdjacencyThreshold & Diff_Time_Left > TimeAdjacencyThreshold),
        Merge_Group = cumsum(Diff_Seg_Left > SpaceAdjacencyThreshold | (Diff_Seg_Left <= SpaceAdjacencyThreshold & Diff_Time_Left > TimeAdjacencyThreshold))
      ) %>%
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index, Merge_Group) %>%
      mutate(
        Majority_Synchronous = sum(belong_to_sync) / n() >= 0.5,
        Segment_Start_After_Merge = min(Segment_Start),
        Segment_End_After_Merge = max(Segment_End),
        Segment_ID_After_Merge = paste(Chr, Segment_Start_After_Merge, Segment_End_After_Merge, sep = "-"),
        Grouped_Timing_Before_Merge = mean(timing_median)
      )
    
    Merged <- Grouped %>% 
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index, Merge_Group) %>%
      dplyr::mutate(
        Segment_Start_After_Merge = min(Segment_Start),
        Segment_End_After_Merge = max(Segment_End),
        Segment_ID_After_Merge = paste(Chr, Segment_Start_After_Merge, Segment_End_After_Merge, sep = "-"),
        Timing_After_Merge = mean(timing_median)
      ) %>% 
      dplyr::ungroup() %>%
      dplyr::select(Anom_Sample_ID, Segment_ID_After_Merge, Chr, Gain_Index, Segment_Start_After_Merge, Segment_End_After_Merge, Timing_After_Merge, Majority_Synchronous) %>%
      unique()
    
  }
  
  if(SegmentMergeRule == "adj_copy_time"){
    print("Merging Segments")
    Grouped <- GainSetReduced %>%
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index) %>%
      dplyr::mutate(
        Diff_Seg_Left = Segment_Start - lag(Segment_End, default = Inf),
        Diff_Time_Left = abs(timing_median - lag(timing_median, default = Inf)),
        Test = Diff_Seg_Left > SpaceAdjacencyThreshold | (Diff_Seg_Left <= SpaceAdjacencyThreshold & Diff_Time_Left > TimeAdjacencyThreshold),
        Merge_Group = cumsum(Diff_Seg_Left > SpaceAdjacencyThreshold | (Diff_Seg_Left <= SpaceAdjacencyThreshold & Diff_Time_Left > TimeAdjacencyThreshold))
      ) %>%
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index, Merge_Group) %>%
      mutate(
        Majority_Synchronous = sum(belong_to_sync) / n() >= 0.5,
        Segment_Start_After_Merge = min(Segment_Start),
        Segment_End_After_Merge = max(Segment_End),
        Segment_ID_After_Merge = paste(Chr, Segment_Start_After_Merge, Segment_End_After_Merge, sep = "-"),
        Grouped_Timing_Before_Merge = mean(timing_median)
      )
    
    Merged <- Grouped %>% 
      dplyr::group_by(Anom_Sample_ID, Chr, Gain_Index, Merge_Group) %>%
      dplyr::mutate(
        Segment_Start_After_Merge = min(Segment_Start),
        Segment_End_After_Merge = max(Segment_End),
        Segment_ID_After_Merge = paste(Chr, Segment_Start_After_Merge, Segment_End_After_Merge, sep = "-"),
        Timing_After_Merge = mean(timing_median)
      ) %>% 
      dplyr::ungroup() %>%
      dplyr::select(Anom_Sample_ID, Segment_ID_After_Merge, Chr, Gain_Index, Segment_Start_After_Merge, Segment_End_After_Merge, Timing_After_Merge, Majority_Synchronous) %>%
      unique()
    
  }
    
    ### Structural Variants
    print("Checking structural variants")
    
    SVFiles <- vector()
    SVNamePlates <- vector()
    
    if(!is.null(SVDirectory)){
      for(i in 1:length(SVDirectory)){
        files <- list.files(SVDirectory[i], pattern = SVSuffix[i], full.names = F)[!grepl("md5", list.files(SVDirectory[i], pattern = SVSuffix[i], full.names = T))]
        paths <- list.files(SVDirectory[i], pattern = SVSuffix[i], full.names = T)[!grepl("md5", list.files(SVDirectory[i], pattern = SVSuffix[i], full.names = T))]
        SVFiles <- append(SVFiles, paths)
        SVNamePlates <- append(SVNamePlates, gsub(SVSuffix[i], "", files))
      }
      
      AllSVs <- read.table(file = gzfile(SVFiles[1]), header = T) %>% mutate(Sample_ID = SVNamePlates[1])
      
      for(i in 2:length(SVFiles)){
        AllSVs <- rbind(AllSVs, read.table(file = gzfile(SVFiles[i]), header = T) %>% mutate(Sample_ID = SVNamePlates[i]))
      }
    } else{
      AllSVs <- data.frame(Sample_ID = vector(), chrom1 = vector(), chrom2 = vector())
    }
    
    Merged$SVFileExists <- Merged$Anom_Sample_ID %in% SVNamePlates
    Merged$SVFile <- SVFiles[match(Merged$Anom_Sample_ID, SVNamePlates)]
    
    CNA_SVs <- lapply(unique(Merged$Anom_Sample_ID), FUN = function(id){
      cna <- Merged[Merged$Anom_Sample_ID == id,]
      sv <- AllSVs[AllSVs$Sample_ID == id,]
      df <- data.frame(sample_id = vector(), segment_id = vector(), sv_class = vector())
      for(seg_id in unique(cna$Segment_ID_After_Merge)){
        seg <- cna[cna$Segment_ID_After_Merge == seg_id,]
        seg_chr <- unique(seg$Chr)
        seg_start <- unique(seg$Segment_Start_After_Merge)
        seg_end <- unique(seg$Segment_End_After_Merge)
        sv_chr <- sv[sv$chrom1 == seg_chr | sv$chrom2 == seg_chr,]
        
        cna_svs_dat <- sv_chr[(sv_chr$chrom1 == seg_chr & (abs(seg_end - sv_chr$start1) < 10 | abs(seg_start - sv_chr$end1) < 10)) | (sv_chr$chrom2 == seg_chr & (abs(seg_end - sv_chr$start2) < 10 | abs(seg_start - sv_chr$end2) < 10)), ]
        df <- rbind(df,data.frame(sample_id = id, segment_id = seg_id, sv_class = ifelse(nrow(cna_svs_dat >= 1), paste(cna_svs_dat$svclass, collapse = " "), "None")))
      }
      return(df)
    }
    ) %>% data.table::rbindlist()
    
    Merged$SV_Class<- CNA_SVs[match(paste(Merged$Anom_Sample_ID, Merged$Segment_ID_After_Merge), paste(CNA_SVs$sample_id, CNA_SVs$segment_id)),]$sv_class
    
    print("Classifying segment type, length, and topology")
    apply(Merged, 1, function(seg){
      # Set up variables
      sample_id = seg["Anom_Sample_ID"]
      seg_id = seg["Segment_ID_After_Merge"]
      seg_chr = seg["Chr"]
      seg_start = as.numeric(seg["Segment_Start_After_Merge"])
      seg_end = as.numeric(seg["Segment_End_After_Merge"])
      seg_length = seg_end - seg_start
      
      centromere_start = GenomeSizes[GenomeSizes$Chr == seg_chr,]$Centromere_Start
      centromere_end = GenomeSizes[GenomeSizes$Chr == seg_chr,]$Centromere_End
      
      chr_buffer = GenomeSizes[GenomeSizes$Chr == seg_chr,]$Length * 0.01 # 1% buffer on all calls
      chr_length = GenomeSizes[GenomeSizes$Chr == seg_chr,]$Length
      
      # Identify segment features
      WholeChromosome <- ifelse(seg_length >= chr_length - chr_buffer, TRUE, FALSE)
      SpansCentromere <- ifelse( (seg_start <= centromere_start - chr_buffer) & (seg_end >= centromere_end + chr_buffer) , TRUE, FALSE)
      TelomereBounded <- ifelse( (seg_start <= chr_buffer) | (seg_end >= chr_length - chr_buffer) , TRUE, FALSE)
      CentromereBounded <- ifelse( (seg_start >= centromere_start - chr_buffer & seg_start <= centromere_end + chr_buffer) | (seg_end >= centromere_start - chr_buffer & seg_end <= centromere_end + chr_buffer) , TRUE, FALSE)
      
      sv_file_exists <- seg["SVFileExists"]
      sv_class <- seg["SV_Class"]
      
      # Classify segment type
      ## Whole chromosome if WholeChromosome is TRUE
      ## Supra-Centromeric if TelomereBounded & SpansCentromere are TRUE & WholeChromosome is FALSE
      ## Arm if TelomereBounded & CentromereBounded are TRUE
      ## Telomere-Bounded if TelomereBounded is TRUE & CentromereBounded is FALSE & SpansCentromere is FALSE
      ## Centromere-Bounded if CentromereBounded is TRUE & TelomereBounded is FALSE
      ## Interstitial if Telomere-Bounded is FALSE & CentromereBounded is FALSE
      
      Form_Class <- ifelse(WholeChromosome == TRUE, "Whole Chromosome",
                              ifelse(WholeChromosome == FALSE & TelomereBounded == TRUE & SpansCentromere == TRUE, "Supra-Centromeric Arm",
                                     ifelse(TelomereBounded == TRUE & CentromereBounded == TRUE, "Whole Arm", 
                                            ifelse(TelomereBounded == TRUE & CentromereBounded == FALSE & SpansCentromere == FALSE, "Sub-Arm Segment (Telomere-Bounded)", 
                                                   ifelse(TelomereBounded == FALSE & CentromereBounded == TRUE, "Sub-Arm Segment (Centromere-Bounded)", 
                                                          ifelse(TelomereBounded == FALSE & CentromereBounded == FALSE, "Sub-Arm Segment (Interstitial)", "Other"))))))
      # Classify lengths into bins
      Length_Class <- ifelse(seg_length < 100000, "<100KB",
                             ifelse(seg_length >= 100000 & seg_length < 500000, "100KB-500KB",
                                    ifelse(seg_length >= 500000 & seg_length < 1000000, "500KB-1MB",
                                           ifelse(seg_length >= 1000000 & seg_length < 5000000, "1MB-5MB",
                                                  ifelse(seg_length >= 5000000 & seg_length < 10000000, "5MB-10MB",
                                                         ifelse(seg_length >= 10000000 & seg_length < 50000000,  "10MB-50MB", 
                                                                ifelse(seg_length >= 50000000,  ">50MB", "Other")))))))
      
      
      # Classify topology
      ## NA if SVFileExists is FALSE
      ## Undetermined if SVFileExists is TRUE & SV_Class is None
      ## Inter-Chromosomal if SVFileExists is TRUE & grepl("TRA", SV_Class) is TRUE
      ## Intra-Chromosomal if SVFileExists is TRUE & SV_Class is not None & grepl("TRA", SV_Class) is FALSE
      ## Other
      Topology_Class <- ifelse(sv_file_exists == FALSE, NA, 
                               ifelse(sv_file_exists == TRUE & sv_class == "None", "None",
                                      ifelse(sv_file_exists == TRUE & grepl("TRA", sv_class) == TRUE, "Inter-Chromosomal",
                                             ifelse(sv_file_exists == TRUE & sv_class != "None" & grepl("TRA", sv_class) == FALSE, "Intra-Chromosomal", "Other"))))
      # Return
      return(
        data.frame(
          WholeChromosome, 
          SpansCentromere, 
          TelomereBounded, 
          CentromereBounded, 
          Length_Class,
          Form_Class,
          Topology_Class)
      )
      
      
    }) %>% data.table::rbindlist() -> ClassifiedSegments
    
    ClassifiedSegments <- cbind(Merged, ClassifiedSegments)
    
    print("Complete")
    return(list(Merged, Grouped, SegmentMergeRule, ClassifiedSegments))
}


MergePlot <- function(MergedGainSet = NULL, Sample = NULL, Chrs = c(1:22,"X","Y")){
  if(df[[3]] == "adj_idx_time"){
    p <- ggplot() + 
      geom_segment(data = GenomeSizes %>% filter(Chr %in% Chrs), aes(x = Start, xend = End, y = 0, yend = 0), color = "transparent") + 
      geom_rect(data = GenomeSizes %>% filter(Chr %in% Chrs), aes(xmin = Centromere_Start, xmax = Centromere_End, ymin = -Inf, ymax = Inf), fill = "grey80", color = "transparent") + 
      geom_rect(data = MergedGainSet[[2]] %>% filter(Anom_Sample_ID == Sample & Chr %in% Chrs), aes(xmin = Segment_Start, xmax = Segment_End, ymin = timing_median - 0.005, ymax = timing_median + 0.005, fill = factor(Gain_Index)), color = "transparent", linewidth = 0.25, alpha = 0.33) +
      geom_rect(data = MergedGainSet[[1]] %>% filter(Anom_Sample_ID == Sample & Chr %in% Chrs), aes(xmin = Segment_Start_After_Merge, xmax = Segment_End_After_Merge, ymin = Timing_After_Merge - 0.005, ymax = Timing_After_Merge + 0.005, fill = factor(Gain_Index)), linewidth = 0.5) + 
      geom_segment(data = MergedGainSet[[2]] %>% filter(Anom_Sample_ID == Sample & Chr %in% Chrs & timing_median != Grouped_Timing_Before_Merge), aes(x = Segment_Start + (Segment_End-Segment_Start)/2, xend = Segment_Start + (Segment_End-Segment_Start)/2, y = timing_median + (0.005*sign(Grouped_Timing_Before_Merge-timing_median)), yend = Grouped_Timing_Before_Merge - (0.005 * sign(Grouped_Timing_Before_Merge-timing_median))), linewidth = 0.25, arrow = grid::arrow(length = unit(1, "mm"))) + 
      ggh4x::facet_grid2(Anom_Sample_ID~Chr, space = "free_x", scale = "free_x", switch = "x") + 
      coord_cartesian(clip = "off") + 
      scale_y_continuous(limits = c(0,1), breaks = c(0,1)) + 
      scale_x_continuous(expand = c(0,0), breaks = c(0,1)) + 
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
  }
  
  if(df[[3]] == "adj_idx_sync"){
    p <- ggplot() + 
      geom_segment(data = GenomeSizes %>% filter(Chr %in% Chrs), aes(x = Start, xend = End, y = 0, yend = 0), color = "transparent") + 
      geom_rect(data = GenomeSizes %>% filter(Chr %in% Chrs), aes(xmin = Centromere_Start, xmax = Centromere_End, ymin = -Inf, ymax = Inf), fill = "grey80", color = "transparent") + 
      geom_rect(data = MergedGainSet[[2]] %>% filter(Anom_Sample_ID == Sample & Chr %in% Chrs), aes(xmin = Segment_Start, xmax = Segment_End, ymin = timing_median - 0.005, ymax = timing_median + 0.005, fill = factor(Gain_Index)), color = "transparent", linewidth = 0.25, alpha = 0.33) +
      geom_rect(data = MergedGainSet[[1]] %>% filter(Anom_Sample_ID == Sample & Chr %in% Chrs), aes(xmin = Segment_Start_After_Merge, xmax = Segment_End_After_Merge, ymin = Timing_After_Merge - 0.005, ymax = Timing_After_Merge + 0.005, fill = factor(Gain_Index)), linewidth = 0.5) + 
      geom_segment(data = MergedGainSet[[2]] %>% filter(Anom_Sample_ID == Sample & Chr %in% Chrs & Majority_Synchronous == TRUE & timing_median != Grouped_Timing_Before_Merge), aes(x = Segment_Start + (Segment_End-Segment_Start)/2, xend = Segment_Start + (Segment_End-Segment_Start)/2, y = timing_median + (0.005*sign(Grouped_Timing_Before_Merge-timing_median)), yend = Grouped_Timing_Before_Merge - (0.005 * sign(Grouped_Timing_Before_Merge-timing_median))), linewidth = 0.25, arrow = grid::arrow(length = unit(1, "mm"))) + 
      ggh4x::facet_grid2(.~Chr, space = "free_x", scale = "free_x", switch = "x") + 
      coord_cartesian(clip = "off") + 
      scale_y_continuous(limits = c(0,1), breaks = c(0,1)) + 
      scale_x_continuous(expand = c(0,0)) + 
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
  }
  
  print(p)
}