library(tidyverse)
library(magrittr)
library(pbapply)
library(data.table)

PATH_WD=""
setwd(PATH_WD)

# Step 1: read the formatted object ------
non_parsimony_penalty <- 'True'
timing <- read_rds(paste0("../in_data/timing_keep_non_parsimony_penalty_",non_parsimony_penalty,"_04092024.rds"))

# Step 2: get 100 posterior samples and get min_dist distribution ------
sample_id=levels(timing$Anom_Sample_ID)
source("min_dist_nearest_timing.R")

res <- pblapply(sample_id, function(x) {
  data.table(sample = x, min_dist = min_dist_nearest_timing(timing_obj = timing, sampleid = x, n_posterior = 100))
}, cl = 40)
res <- rbindlist(res)
res <- separate_wider_delim(res, cols = min_dist, names = c("min_dist", "timing"), delim = "|")
res$min_dist <- as.numeric(res$min_dist)
res$timing <- as.numeric(res$timing)/101
# save as a txt file; each sample has a median from 100 posterior draws
write.table(res, file=paste0(PATH_WD, paste0("min_dist_median_non_parsimony_penalty_",non_parsimony_penalty,".txt")), quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
message("Min_dist_median and timing saved.")
