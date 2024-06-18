# Function: 
# select N times from posteriors;
# run Weighted distance method and only keep the nearest gain to the line;
# return the median from a vector of N min_dist
min_dist_nearest_timing <- function(timing_obj=timing_keep, sampleid="SAMPLE1", n_posterior=50){
  s <- dplyr::filter(timing_obj, Anom_Sample_ID %in% sampleid) %>% droplevels()
  s <- dplyr::group_by(s, Segment_ID)
  setDT(s)
  min_dist <- rep(NA,n_posterior)
  synctiming <- rep(NA,n_posterior)
  for (i in 1:n_posterior){
    # for a given sample, choose a Sample_Index randomly for each unique Segment_ID to make up a posterior-sample
    s0 <- s[s[, .I[Posterior_Sample_Index == sample(levels(Posterior_Sample_Index), 1)], by = Segment_ID]$V1]
    # Then calculate weighted sum minimum but only keep nearest gain per each segment
    sq <- (0:100)/100
    dist.w <- rep(NA,length(sq))
    length_total <- sum(s0$Segment_Width[s0[ , .I[1], by = Segment_ID]$V1])
    for (j in 1:length(sq)){
      s00 <- s0
      s00[, dist := abs(Gain_Timing-sq[j])]
      s00 <- s00[s00[ , .I[which.min(dist)], by = Segment_ID]$V1]
      dist.w[j] <-sum(s00$dist*s00$Segment_Width)/length_total
    }
    min_dist[i] <- min(dist.w)
    synctiming[i] <- which.min(dist.w)
  }
  return(paste(median(min_dist), median(synctiming), sep = "|"))
}
