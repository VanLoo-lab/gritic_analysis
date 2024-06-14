library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(data.table)
library(stringr)


setwd("/nemo/project/proj-vanloo-secure/working/bakert/GRITIC/Analysis_Pipeline/code")


hartwig_metadata <- fread("/camp/project/proj-vanloo-secure/Hartwig/DR-163-update2/metadata/metadata.tsv")

all_hartwig_copy_numbers <- list()
for (metadata_index in 1:nrow(hartwig_metadata)) {
  set_name <- as.character(hartwig_metadata[metadata_index,'setName'])
  sample_id <- as.character(hartwig_metadata[metadata_index,'sampleId'])
  cn_path <- paste0("/camp/project/proj-vanloo-secure/Hartwig/DR-163-update2/somatics/",sample_id,"/purple/",sample_id,".purple.cnv.somatic.tsv")
  if(file.exists(cn_path)){
  cn_table <- fread(cn_path) %>% dplyr::rename(Chromosome=chromosome,Segment_Start=start,Segment_End=end)%>% dplyr::mutate(Major_CN=as.integer(round(majorAlleleCopyNumber)),Minor_CN=as.integer(round(minorAlleleCopyNumber)))%>% dplyr::select(Chromosome,Segment_Start,Segment_End,Major_CN,Minor_CN)
  cn_table <- cn_table %>% mutate(Sample_ID=sample_id)
  all_hartwig_copy_numbers[[metadata_index]] <- cn_table
  }
  print(metadata_index)

}

all_pcawg_copy_numbers <- list()
pcawg_metadata <- fread("/camp/project/proj-vanloo-secure/PCAWG/ICGC_annotations/summary_table_combined_annotations_v2.txt")
sample_id_store <- character()
wgd_status_store <- logical()
for (metadata_index in 1:nrow(pcawg_metadata)) {
  sample_id <- as.character(pcawg_metadata[metadata_index,'samplename'])  
  donor_id <- as.character(pcawg_metadata[metadata_index,'icgc_donor_id'])  
  project_code <- as.character(pcawg_metadata[metadata_index,'projectcode']) 
  is_tcga <- !is.na(pcawg_metadata$tcga_donor_uuid[metadata_index])
cn_dir <- "/camp/project/proj-vanloo-secure/PCAWG/Hartwig_Pipeline/"
    if(is_tcga){
cn_path <- paste0(cn_dir,'TCGA/',project_code,'/',donor_id,'/purple/',donor_id,'T.purple.cnv.somatic.tsv')
wgd_path <- paste0(cn_dir,'TCGA/',project_code,'/',donor_id,'/purple/',donor_id,'T.purple.purity.tsv')
  }else{
cn_path <- paste0(cn_dir,'ICGC/',project_code,'/',donor_id,'/purple/',donor_id,'T.purple.cnv.somatic.tsv')
wgd_path <- paste0(cn_dir,'ICGC/',project_code,'/',donor_id,'/purple/',donor_id,'T.purple.purity.tsv')
  }
print(metadata_index)
if (!file.exists(cn_path)){
  next
}
  cn_table <- fread(cn_path) %>% dplyr::rename(Chromosome=chromosome,Segment_Start=start,Segment_End=end)%>% dplyr::mutate(Major_CN=as.integer(round(majorAlleleCopyNumber)),Minor_CN=as.integer(round(minorAlleleCopyNumber)))%>% dplyr::select(Chromosome,Segment_Start,Segment_End,Major_CN,Minor_CN)
  cn_table <- cn_table %>% mutate(Sample_ID=sample_id)
  all_pcawg_copy_numbers[[metadata_index]] <- cn_table
  
  sample_wgd_status <- tolower(fread(wgd_path)$wholeGenomeDuplication[1]) =='true'
  sample_id_store <- c(sample_id_store,sample_id)
  wgd_status_store <- c(wgd_status_store,sample_wgd_status)

}
hartwig_copy_number <- do.call("rbind", all_hartwig_copy_numbers)
pcawg_copy_number <- do.call("rbind", all_pcawg_copy_numbers)

pcawg_wgd_status_table <- data.frame(Sample_ID=sample_id_store,WGD_Status=wgd_status_store)


hartwig_copy_number <- hartwig_copy_number %>% dplyr::mutate(Tumor_Type="Metastatic")
pcawg_copy_number <- pcawg_copy_number %>% dplyr::mutate(Tumor_Type="Primary")


cn_table <- rbind(hartwig_copy_number,pcawg_copy_number)
cn_table <- cn_table %>% dplyr::mutate(Segment_Width=Segment_End-Segment_Start,CN_State=paste0(Major_CN,"+",Minor_CN),Total_CN=Major_CN+Minor_CN)%>%filter(Major_CN>=0 & Minor_CN>=0 & Major_CN >= Minor_CN)


hartwig_wgd_status_table <- read_tsv('/camp/project/proj-vanloo-secure/Hartwig/DR-163-update2/metadata/wgd_status.tsv')
wgd_status_table <- rbind(pcawg_wgd_status_table,hartwig_wgd_status_table)
cn_table_status <- inner_join(cn_table,wgd_status_table)
fwrite(cn_table_status,'../output/cn_table_status.tsv',sep="\t")


