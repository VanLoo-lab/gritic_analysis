import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

#font size 18, pdf font type 42, nimbus sans font
plt.rcParams.update({'font.size': 18, 'pdf.fonttype': 42, 'font.family': 'Nimbus Sans'})

import DataTools

def get_pre_post_by_chromosome(pre_post_table):

    pre_post_chr = pre_post_table.groupby(['Chromosome','Major_CN']).apply(lambda x: np.average(x['Pre_Probability'], weights=x['Segment_Width'])).reset_index()
    #pre_post_chr = pre_post.groupby(['Chromosome','Major_CN']).apply(lambda x: np.average(x['Pre_Probability'])).reset_index()
    pre_post_chr.columns = ['Chromosome','Major_CN','Pre_Probability']
    #pivot wider on Major CN
    pre_post_chr = pre_post_chr.pivot(index='Chromosome', columns='Major_CN', values='Pre_Probability').reset_index()
    pre_post_chr.columns = ['Chromosome','Pre_Probability_3','Pre_Probability_4']

    chromosome_order = list(map(str, range(1,23))) 
    pre_post_chr['Chromosome'] = pd.Categorical(pre_post_chr['Chromosome'], categories=chromosome_order, ordered=True)
    pre_post_chr = pre_post_chr.sort_values('Chromosome')
    return pre_post_chr

def load_pre_post_table(apply_penalty):   
    posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty)
    posterior_table = posterior_table[posterior_table['Major_CN'].isin([3,4]) & posterior_table['WGD_Status']]
    seg_info = posterior_table[['Sample_ID','Segment_ID','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN','N_Mutations']].drop_duplicates()
    
    pre_post = pd.read_csv(f'../output/pre_post_non/pre_post_non_gains_apply_penalty_{apply_penalty}_major_cn_All_mode_absolute.tsv',sep="\t")
    
    pre_post = seg_info.merge(pre_post, how='inner')
    pre_post['Segment_Width'] = pre_post['Segment_End'] - pre_post['Segment_Start']
    return pre_post

def plot_scatter(pre_post_chr_all,apply_penalty):
    cor,p_value = pearsonr(pre_post_chr_all['Pre_Probability_3'], pre_post_chr_all['Pre_Probability_4'])
    fig,ax = plt.subplots(figsize=(10,10))
    ax.scatter(pre_post_chr_all['Pre_Probability_3'], pre_post_chr_all['Pre_Probability_4'], s=100, color='#e87641')
    ax.set_xlabel('Major CN 3 Pre-WGD Gain Probability')
    ax.set_ylabel('Major CN 4 Pre-WGD Gain Probability')
    #add labels for chromosomes
    for i, txt in enumerate(pre_post_chr_all['Chromosome']):
        ax.annotate(txt, (pre_post_chr_all['Pre_Probability_3'][i], pre_post_chr_all['Pre_Probability_4'][i]), fontsize=15)
    #annoate with correlation
    ax.annotate(f"Pearson's correlation: {cor:.2f}\n P= {p_value:.2e}", (0.05,0.85), xycoords='axes fraction', fontsize=20)
    plt.savefig(f'../plots/pre_post_major_cn/chromosome_scatter_{apply_penalty}.pdf', bbox_inches='tight')

def plot_bar(pre_post_chr_all,apply_penalty):
    fig,axs = plt.subplots(2,1,figsize=(20,10),sharex=True)
    #barplot with bootstrap error bars
    axs[0].bar(pre_post_chr_all['Chromosome'], pre_post_chr_all['Pre_Probability_3'], color='#4392F1', label='Pre-probability of 3 CN', yerr=[pre_post_chr_all['Pre_Probability_3'] - pre_post_chr_all['Pre_Probability_3_Low_CI'], pre_post_chr_all['Pre_Probability_3_High_CI'] - pre_post_chr_all['Pre_Probability_3']], capsize=10)
    axs[0].set_ylabel('Pre-WGD Gain Probability')
    axs[0].set_title('Major Copy Number = 3')
    axs[1].bar(pre_post_chr_all['Chromosome'], pre_post_chr_all['Pre_Probability_4'], color='#E06052', label='Pre-probability of 4 CN', yerr=[pre_post_chr_all['Pre_Probability_4'] - pre_post_chr_all['Pre_Probability_4_Low_CI'], pre_post_chr_all['Pre_Probability_4_High_CI'] - pre_post_chr_all['Pre_Probability_4']], capsize=10)
    axs[1].set_ylabel('Pre-WGD Gain Probability')
    axs[1].set_title('Major Copy Number = 4')
    axs[1].set_xlabel('Chromosome')
    plt.savefig(f'../plots/pre_post_major_cn/chromosome_bars_{apply_penalty}.pdf', bbox_inches='tight')

def get_bootstrap_summary(pre_post_sample_dict,n_bootstraps=250):
    pre_post_chr_bootstraps_store = []
    for i in range(n_bootstraps):
        bootstrap_table = DataTools.get_bootstrap_table(pre_post_sample_dict)

        pre_post_chr_bootstrap = get_pre_post_by_chromosome(bootstrap_table)
        pre_post_chr_bootstrap['Bootstrap'] = i
        pre_post_chr_bootstraps_store.append(pre_post_chr_bootstrap)
        print(i)
    pre_post_chr_bootstraps_store = pd.concat(pre_post_chr_bootstraps_store)
    #pre_post_chr_bootstraps_store = pd.concat(pre_post_chr_bootstraps_store)
    pre_post_chr_bootstrap_summary = pre_post_chr_bootstraps_store.groupby(['Chromosome']).agg({'Pre_Probability_3': [lambda x: np.percentile(x, 5/2), lambda x: np.percentile(x, 100-5/2)], 'Pre_Probability_4': [lambda x: np.percentile(x, 11/2), lambda x: np.percentile(x, 100-11/2)]}).reset_index()
    #flatten multiindex
    pre_post_chr_bootstrap_summary.columns = ['Chromosome','Pre_Probability_3_Low_CI','Pre_Probability_3_High_CI','Pre_Probability_4_Low_CI','Pre_Probability_4_High_CI']
    return pre_post_chr_bootstrap_summary

if __name__ == '__main__':

    for apply_penalty in [False,True]:
        pre_post_table = load_pre_post_table(apply_penalty)
        pre_post_chr = get_pre_post_by_chromosome(pre_post_table)
        pre_post_sample_dict = DataTools.get_sample_dict(pre_post_table)
        
        pre_post_chr_bootstrap_summary = get_bootstrap_summary(pre_post_sample_dict)
        pre_post_chr_all = pre_post_chr.merge(pre_post_chr_bootstrap_summary, how='inner')
        print(pre_post_chr_all)
        plot_scatter(pre_post_chr_all,apply_penalty)

