library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(cowplot)
library(ggpubr)



out_dir <- '../plots/pre_post_correlations_corrected/'
load_counts = function(event_type,min_samples_per_cancer_type=50,apply_penalty='False'){
    if(event_type=='Gain'){
        event_path <- paste0('../output/arm_pre_post_clean/gain_timing_',apply_penalty,'_bin_counts.tsv')
    }
    else if (event_type=='Loss'){
        event_path <- '../output/arm_pre_post_clean/loss_timing_bin_counts.tsv'
    }
    event_table <- fread(event_path)
    
    count_table <- event_table %>% dplyr::select(Cancer_Type,Sample_ID,WGD_Status)%>% distinct() %>% group_by(Cancer_Type,WGD_Status)%>% tally()
    count_table <- count_table %>% pivot_wider(names_from=WGD_Status,values_from = n,names_prefix = 'WGD_') %>% dplyr::mutate(WGD_TRUE=replace_na(WGD_TRUE,0),WGD_FALSE=replace_na(WGD_FALSE,0))
    count_table <- count_table %>% dplyr::filter(WGD_TRUE[1]>=min_samples_per_cancer_type&WGD_FALSE[1]>=min_samples_per_cancer_type)
    event_table <- event_table %>% dplyr::filter(Cancer_Type %in% count_table$Cancer_Type)
    event_table <- inner_join(event_table,count_table)
    acrocentric_chr <- c('13p','14p','15p','21p','22p')
    event_table <- event_table %>% dplyr::filter(!Bin_ID %in% acrocentric_chr)
    return(event_table)
}



setwd("/nemo/project/proj-vanloo-secure/working/bakert/GRITIC/Analysis_Pipeline/code")

for (apply_penalty in c('True','False')){
gain_counts <- load_counts('Gain',apply_penalty=apply_penalty) 

loss_counts <- load_counts('Loss')



gain_proportions <- gain_counts %>% group_by(Cancer_Type,Bin_ID)%>% dplyr::summarise(Prop_Pre_Gain=sum(Pre_WGD_Gain>0.5)/WGD_TRUE[1],Prop_Post_Gain=sum(Post_WGD_Gain>0.5)/WGD_TRUE[1],Prop_Non_Gain=sum(Non_WGD_Gain>0.5)/WGD_FALSE[1])
loss_proportions <- loss_counts %>% group_by(Cancer_Type,Bin_ID)%>% dplyr::summarise(Prop_Pre_Loss=sum(Pre_WGD_Loss>0.5)/WGD_TRUE[1],Prop_Post_Loss=sum(Post_WGD_Loss>0.5)/WGD_TRUE[1],Prop_Non_Loss=sum(Non_WGD_Loss>0.5)/WGD_FALSE[1])
loss_proportions <- loss_proportions %>% dplyr::mutate(Corrected_Prop_Post_WGD_Loss=Prop_Post_Loss/(1-Prop_Pre_Loss))
gain_proportions<- gain_proportions %>% group_by(Cancer_Type)%>%mutate(Prop_Pre_Gain_Norm=Prop_Pre_Gain/mean(Prop_Pre_Gain),Prop_Post_Gain_Norm=Prop_Post_Gain/mean(Prop_Post_Gain))
loss_proportions <- loss_proportions %>%group_by(Cancer_Type)%>%mutate(Prop_Pre_Loss_Norm=Prop_Pre_Loss/mean(Prop_Pre_Loss),Prop_Post_Loss_Norm=Prop_Post_Loss/mean(Prop_Post_Loss),Corrected_Prop_Post_WGD_Loss_Norm=Corrected_Prop_Post_WGD_Loss/mean(Corrected_Prop_Post_WGD_Loss))

wgd_counts_gain_combined <- gain_counts %>% dplyr::select(Cancer_Type,WGD_TRUE,WGD_FALSE)%>% dplyr::distinct() %>% summarise(WGD_TRUE=sum(WGD_TRUE),WGD_FALSE=sum(WGD_FALSE))
wgd_counts_loss_combined <- loss_counts %>% dplyr::select(Cancer_Type,WGD_TRUE,WGD_FALSE)%>% dplyr::distinct() %>% summarise(WGD_TRUE=sum(WGD_TRUE),WGD_FALSE=sum(WGD_FALSE))

gain_proportions_combined <- gain_counts %>% dplyr::select(-WGD_TRUE,-WGD_FALSE) %>% dplyr::mutate(WGD_TRUE=wgd_counts_gain_combined$WGD_TRUE[1],WGD_FALSE=wgd_counts_gain_combined$WGD_FALSE[1]) %>% group_by(Bin_ID)%>%dplyr::summarise(Prop_Pre_Gain=sum(Pre_WGD_Gain>0.5)/WGD_TRUE[1],Prop_Post_Gain=sum(Post_WGD_Gain>0.5)/WGD_TRUE[1],Prop_Non_Gain=sum(Non_WGD_Gain>0.5)/WGD_FALSE[1])
loss_proportions_combined <- loss_counts %>% dplyr::select(-WGD_TRUE,-WGD_FALSE) %>% dplyr::mutate(WGD_TRUE=wgd_counts_loss_combined$WGD_TRUE[1],WGD_FALSE=wgd_counts_loss_combined$WGD_FALSE[1]) %>% group_by(Bin_ID)%>%dplyr::summarise(Prop_Pre_Loss=sum(Pre_WGD_Loss>0.5)/WGD_TRUE[1],Prop_Post_Loss=sum(Post_WGD_Loss>0.5)/WGD_TRUE[1],Prop_Non_Loss=sum(Non_WGD_Loss>0.5)/WGD_FALSE[1])
loss_proportions_combined <- loss_proportions_combined %>% dplyr::mutate(Corrected_Prop_Post_WGD_Loss=Prop_Post_Loss/(1-Prop_Pre_Loss))

loss_proportions <- loss_proportions%>% dplyr::mutate(Arm_Color=ifelse(Bin_ID%in%c('17p','9p'),Bin_ID,'Other'))
gain_proportions %>% ggplot(aes(x=Prop_Non_Gain,y=Prop_Pre_Gain+Prop_Post_Gain))+geom_point()+geom_abline(slope=1,intercept=0)


base_theme <-theme_classic()
color_scheme <- scale_color_manual(values=c('#85C2FF','#E675A2','#878787'))

pre_non_loss_plot <- loss_proportions %>% ggplot(aes(x=Prop_Pre_Loss,y=Prop_Non_Loss))+geom_point(aes(color=Arm_Color))+base_theme+labs(x='Pre-WGD Loss Proportion',y='Non-WGD Loss Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed') 
pre_post_loss_plot <- loss_proportions %>% ggplot(aes(x=Prop_Pre_Loss,y=Corrected_Prop_Post_WGD_Loss))+geom_point(aes(color=Arm_Color))+base_theme +labs(x='Pre-WGD Loss Proportion',y='Corrected Post-WGD Loss Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed') +ylim(0,1)
non_post_loss_plot <- loss_proportions %>% ggplot(aes(x=Prop_Non_Loss,y=Corrected_Prop_Post_WGD_Loss))+geom_point(aes(color=Arm_Color)) +base_theme+labs(x='Non-WGD Loss Proportion',y='Corrected Post-WGD Loss Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed') 

pre_uncorrected_post_loss_plot <- loss_proportions %>% ggplot(aes(x=Prop_Pre_Loss,y=Prop_Post_Loss))+geom_point(aes(color=Arm_Color))+base_theme +labs(x='Pre-WGD Loss Proportion',y='Uncorrected Post-WGD Loss Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed') +ylim(0,1)
non_uncorrrected_post_loss_plot <- loss_proportions %>% ggplot(aes(x=Prop_Non_Loss,y=Prop_Post_Loss))+geom_point(aes(color=Arm_Color)) +base_theme+labs(x='Non-WGD Loss Proportion',y='Uncorrected Post-WGD Loss Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed') 

loss_proportions %>% ggplot(aes(x=Prop_Pre_Loss,y=Corrected_Prop_Post_WGD_Loss))+geom_point(aes(color=Arm_Color))+base_theme +labs(x='Pre-WGD Loss Proportion',y='Corrected Post-WGD Loss Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed') +ylim(0,1)+facet_wrap(~Cancer_Type)

ggsave(paste0(out_dir,'pre_post_loss_cancer_type.pdf'),width=12,height=6)
gain_proportions %>% ggplot(aes(x=Prop_Pre_Gain,y=Prop_Post_Gain))+geom_point(color='#878787')+base_theme +labs(x='Pre-WGD Gain Proportion',y='Post-WGD Gain Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho') +geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed')+facet_wrap(~Cancer_Type)
ggsave(paste0(out_dir,'pre_post_gain_cancer_type_apply_penalty',apply_penalty,'.pdf'),width=12,height=6)

pre_non_gain_plot <- gain_proportions %>% ggplot(aes(x=Prop_Pre_Gain,y=Prop_Non_Gain))+geom_point(color='#878787')+base_theme+labs(x='Pre-WGD Gain Proportion',y='Non-WGD Gain Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed') 
pre_post_gain_plot <- gain_proportions %>% ggplot(aes(x=Prop_Pre_Gain,y=Prop_Post_Gain))+geom_point(color='#878787')+base_theme +labs(x='Pre-WGD Gain Proportion',y='Post-WGD Gain Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho') +geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed')
non_post_gain_plot <- gain_proportions %>% ggplot(aes(x=Prop_Non_Gain,y=Prop_Post_Gain))+geom_point(color='#878787') +base_theme+labs(x='Non-WGD Gain Proportion',y='Post-WGD Gain Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed')

pre_close_gain_plot <- gain_proportions %>% ggplot(aes(x=Prop_Pre_Gain,y=Prop_Close_Gain))+geom_point(color='#878787')+base_theme+labs(x='Pre-WGD Gain Proportion',y='Close-WGD Gain Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed')
post_close_gain_plot <- gain_proportions %>% ggplot(aes(x=Prop_Post_Gain,y=Prop_Close_Gain))+geom_point(color='#878787')+base_theme+labs(x='Post-WGD Gain Proportion',y='Close-WGD Gain Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed')
non_close_gain_plot <- gain_proportions %>% ggplot(aes(x=Prop_Non_Gain,y=Prop_Close_Gain))+geom_point(color='#878787')+base_theme+labs(x='Non-WGD Gain Proportion',y='Close-WGD Gain Proportion',color='Arm')+color_scheme+stat_cor(method = "spearman",cor.coef.name='rho')+geom_abline(slope=1,intercept=0,color='red',size=1,linetype='dashed')


plot_grid(pre_post_gain_plot, pre_post_loss_plot,ncol=2)
ggsave(paste0(out_dir,'pre_post_gain_loss_main_apply_penalty_',apply_penalty,'.pdf'),width=12,height=6)
plot_grid(pre_post_loss_plot,pre_non_loss_plot,non_post_loss_plot,pre_uncorrected_post_loss_plot,non_uncorrrected_post_loss_plot,ncol=3)
ggsave(paste0(out_dir,'all_loss_correlations.pdf'),width=16,height=8)
plot_grid(pre_post_gain_plot,pre_non_gain_plot,non_post_gain_plot,ncol=3)
ggsave(paste0(out_dir,'all_gain_correlations_apply_penalty_',apply_penalty,'.pdf'),width=10,height=5)
plot_grid(pre_uncorrected_post_loss_plot,ncol=1)
ggsave(paste0(out_dir,'pre_post_uncorrected.pdf'),width=8,height=5)
plot_grid(pre_post_loss_plot,ncol=1)
ggsave(paste0(out_dir,'pre_post_corrected.pdf'),width=8,height=4)
plot_grid(pre_post_loss_plot,pre_non_loss_plot,non_post_loss_plot,ncol=3)
ggsave(paste0(out_dir,'all_loss_correlations_corrected.pdf'),width=11,height=5)

davoli_arm_scores <- fread('../resources/davoli_arm_scores.txt')

davoli_arm_scores <- davoli_arm_scores %>% dplyr::select(Arm,Density_TSG_in_the_arm,Density_OG_in_the_arm,Density_Essenatil_in_the_arm,N_TSGs_in_the_arm,N_Essential_in_the_arm,N_OGs_in_the_arm,Total_N_Genes_analyzed) %>%dplyr::rename(Bin_ID=Arm,Essential_Density=Density_Essenatil_in_the_arm,TSG_Density=Density_TSG_in_the_arm,OG_Density=Density_OG_in_the_arm,N_TSGs=N_TSGs_in_the_arm,N_Essential=N_Essential_in_the_arm,N_OGs=N_OGs_in_the_arm,N_Genes=Total_N_Genes_analyzed)
davoli_arm_scores <- davoli_arm_scores %>% dplyr::mutate(Density_TSG_in_the_arm=N_TSGs/N_Genes,Density_TSG_in_the_arm=N_OGs/N_Genes,Density_Essential_in_the_arm=N_Essential/N_Genes)

gain_proportions_davoli <- inner_join(gain_proportions_combined,davoli_arm_scores)
loss_proportions_davoli <- inner_join(loss_proportions_combined,davoli_arm_scores)

pre_gain_tsg_plot <- gain_proportions_davoli %>% ggplot(aes(x=TSG_Density,y=Prop_Pre_Gain)) + geom_point(color='#878787') +stat_cor(method="spearman",cor.coef.name='rho',size=6) + base_theme +labs(x='TSG Density',y='Pre-WGD Gain Proportion') + theme(text = element_text(size=20))
pre_gain_og_plot <- gain_proportions_davoli %>% ggplot(aes(x=OG_Density,y=Prop_Pre_Gain)) + geom_point(color='#878787') +stat_cor(method="spearman",cor.coef.name='rho',size=6)  + base_theme +labs(x='OG Density',y='Pre-WGD Gain Proportion') + theme(text = element_text(size=20))
post_gain_tsg_plot <- gain_proportions_davoli %>% ggplot(aes(x=TSG_Density,y=Prop_Post_Gain)) + geom_point(color='#878787') +stat_cor(method="spearman",cor.coef.name='rho',size=6)  + base_theme +labs(x='TSG Density',y='Post-WGD Gain Proportion') + theme(text = element_text(size=20))
post_gain_og_plot <- gain_proportions_davoli %>% ggplot(aes(x=OG_Density,y=Prop_Post_Gain)) + geom_point(color='#878787') +stat_cor(method="spearman",cor.coef.name='rho',size=6) + base_theme +labs(x='OG Density',y='Post-WGD Gain Proportion') + theme(text = element_text(size=20))

non_gain_tsg_plot <- gain_proportions_davoli %>% ggplot(aes(x=TSG_Density,y=Prop_Non_Gain)) + geom_point(color='#878787') +stat_cor(method="spearman",cor.coef.name='rho',size=6) + base_theme +labs(x='TSG Density',y='Non-WGD Gain Proportion') + theme(text = element_text(size=20))
non_gain_og_plot <- gain_proportions_davoli %>% ggplot(aes(x=OG_Density,y=Prop_Non_Gain)) + geom_point(color='#878787') +stat_cor(method="spearman",cor.coef.name='rho',size=6) + base_theme +labs(x='OG Density',y='Non-WGD Gain Proportion') + theme(text = element_text(size=20))

plot_grid(pre_gain_tsg_plot,pre_gain_og_plot,post_gain_tsg_plot,post_gain_og_plot,non_gain_tsg_plot,non_gain_og_plot,ncol=2)
ggsave(paste0(out_dir,'gain_davoli_',apply_penalty,'.pdf'),width=12,height=16)



pre_loss_tsg_plot <- loss_proportions_davoli %>% ggplot(aes(x=TSG_Density,y=Prop_Pre_Loss)) + geom_point(color='#878787') +stat_cor(size=6,method="spearman",cor.coef.name='rho') + base_theme +labs(x='TSG Density',y='Pre-WGD Loss Proportion')+ theme(text = element_text(size=20))
pre_loss_og_plot <- loss_proportions_davoli %>% ggplot(aes(x=OG_Density,y=Prop_Pre_Loss)) + geom_point(color='#878787') +stat_cor(size=6,method="spearman",cor.coef.name='rho')  + base_theme +labs(x='OG Density',y='Pre-WGD Loss Proportion')+ theme(text = element_text(size=20))
post_loss_tsg_plot <- loss_proportions_davoli %>% ggplot(aes(x=TSG_Density,y=Prop_Post_Loss)) + geom_point(color='#878787') +stat_cor(size=6,method="spearman",cor.coef.name='rho')  + base_theme +labs(x='TSG Density',y='Post-WGD Loss Proportion')+ theme(text = element_text(size=20))
post_loss_og_plot <- loss_proportions_davoli %>% ggplot(aes(x=OG_Density,y=Prop_Post_Loss)) + geom_point(color='#878787') +stat_cor(size=6,method="spearman",cor.coef.name='rho') + base_theme +labs(x='OG Density',y='Post-WGD Loss Proportion')+ theme(text = element_text(size=20))
corrected_post_loss_tsg_plot <- loss_proportions_davoli %>% ggplot(aes(x=TSG_Density,y=Corrected_Prop_Post_WGD_Loss)) + geom_point(color='#878787') +stat_cor(size=6)  + base_theme +labs(x='TSG Density',y='Corrected Post-WGD Loss Proportion')+ theme(text = element_text(size=20))
corrected_post_loss_og_plot <- loss_proportions_davoli %>% ggplot(aes(x=OG_Density,y=Corrected_Prop_Post_WGD_Loss)) + geom_point(color='#878787') +stat_cor(size=6) + base_theme +labs(x='OG Density',y='Corrected Post-WGD Loss Proportion')+ theme(text = element_text(size=20))

non_loss_tsg_plot <- loss_proportions_davoli %>% ggplot(aes(x=TSG_Density,y=Prop_Non_Loss)) + geom_point(color='#878787') +stat_cor(method="spearman",cor.coef.name='rho',size=6)  + base_theme +labs(x='TSG Density',y='Non-WGD Loss Proportion') + theme(text = element_text(size=20))
non_loss_og_plot <- loss_proportions_davoli %>% ggplot(aes(x=OG_Density,y=Prop_Non_Loss)) + geom_point(color='#878787') +stat_cor(method="spearman",cor.coef.name='rho',size=6) + base_theme +labs(x='OG Density',y='Non-WGD Loss Proportion') + theme(text = element_text(size=20))


plot_grid(pre_loss_tsg_plot,pre_loss_og_plot,post_loss_tsg_plot,post_loss_og_plot,corrected_post_loss_tsg_plot,corrected_post_loss_og_plot,non_loss_tsg_plot,non_loss_og_plot,ncol=2)
ggsave(paste0(out_dir,'loss_davoli.pdf'),width=12,height=20)

}
