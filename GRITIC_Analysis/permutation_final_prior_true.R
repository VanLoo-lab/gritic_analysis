args = commandArgs(trailingOnly=TRUE)
jobid=as.integer(args[1])
set.seed(jobid)

library(tidyverse)
library(magrittr)
library(pbapply)
library(data.table)

#set the working directory
PATH_WD = ''
setwd(PATH_WD)
#obtained from https://www.nature.com/articles/ncomms5114
source("Curveball.R")


#change to non_parsimony_penalty false to run without non-parsimony penalty
# Step 1: read the formatted object and do permutation ------
non_parsimony_penalty <- 'True'
timing <- read_rds(paste0("../in_data/timing_keep_non_parsimony_penalty_",non_parsimony_penalty,"_04092024.rds"))
timing_p <- data.frame()
cancer_id <- unique(timing$Cancer_Type)
for (CT in cancer_id){
  timing_CT <- dplyr::filter(timing, Cancer_Type == CT) %>% droplevels()
  meta <- dplyr::select(timing_CT, Chr, Anom_Sample_ID)
  
  mat <- as.data.frame.matrix(table(timing_CT$Chr, timing_CT$Anom_Sample_ID))
  mat[mat != 0] <- 1 # binary matrix
  mat_p <- as.data.frame(curve_ball(mat)) # permutation process
  rownames(mat_p) <- rownames(mat)
  colnames(mat_p) <- colnames(mat)
  
  hp=list()
  for (col in 1:ncol(mat)) {hp[[col]]=(which(mat[, col]==1))} # get which chromosomes each sample has
  hp_p=list()
  for (col in 1:ncol(mat_p)) {hp_p[[col]]=(which(mat_p[, col]==1))}
  
  meta$New_ID <- NA
  for (sample in 1:ncol(mat_p)){
    chr_s <- hp_p[[sample]]
    for (c in chr_s) {
      pool <- which(mat[c, ] == 1) # all possible samples with chromosome c for swapping
      pool_select <- if_else(length(pool) == 1, pool[1], sample(pool,1)) # fixed a bug here! when pool is of length 1, should return pool instead of choosing one from 1:pool
      meta$New_ID[meta$Chr == rownames(mat)[c] & meta$Anom_Sample_ID == colnames(mat)[pool_select]] <- colnames(mat_p)[sample]
      mat[c,pool_select] <- 0 # used chromosome should be removed from the pool; no replacement
    }
  }
  
  timing_CT$Anom_Sample_ID <- meta$New_ID # this cancer type is randomized
  # notice that after permutation, the column WGD_Status, WGD_Timing etc. is not useful anymore
  timing_p <- rbind(timing_p, timing_CT)
}
message("Permutation Finished.")
message(Sys.time())
# now the object timing_p is a permutation from timing

# Step 2: get 60 posterior samples and get max_ratio distribution ------
sample_id=levels(timing$Anom_Sample_ID) # be aware that here timing_p doesn't have levels in Anom_Sample_ID
source("min_dist_nearest_timing.R")

res <- pblapply(sample_id, function(x) {
  data.table(sample = x, min_dist = min_dist_nearest_timing(timing_obj = timing_p, sampleid = x, n_posterior = 60))
}, cl = 3)
res <- rbindlist(res)
res <- separate_wider_delim(res, cols = min_dist, names = c("min_dist", "timing"), delim = "|")
res$min_dist <- as.numeric(res$min_dist)
res$timing <- as.numeric(res$timing)/101

# save as a txt file; each sample has a median from posterior draws and the timing for the most sync time.
write.table(res, file=paste0(PATH_WD,"permutation_final_non_parsimony_penalty_",non_parsimony_penalty,"/min_dist_median_p", jobid, ".txt"), quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
message("Min_dist_median saved.")
message(Sys.time())

