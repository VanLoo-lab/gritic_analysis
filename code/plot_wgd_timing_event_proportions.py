import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#font size 22, pdf font type 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams.update({'font.size': 18})
plt.rcParams['font.sans-serif'] = 'Nimbus Sans'
plt.rcParams['pdf.fonttype'] = 42
import DataTools
from event_proportions_relative_to_wgd_probabilistic import load_gain_timing_table,load_loss_timing_table,load_chromosome_table
from scipy.stats import spearmanr,pearsonr

def format_p_value(p):
    if p<0.01:

        p_str = f'{p:.1e}'
        p_str = f'{p_str.split("e")[0]}$\\times$10$^{{{int(p_str.split("e")[1])}}}$'

        return p_str

    else:
        #2sf
        return f'{p:.2f}'

def event_past_tense(event_type):
    if event_type == 'Gain':
        return 'Gained'
    elif event_type == 'Loss':
        return 'Lost'
    else:
        raise ValueError(f'Event {event_type} not recognized')

def event_plural(event_type):
    if event_type == 'Gain':
        return 'Gains'
    elif event_type == 'Loss':
        return 'Losses'
    else:
        raise ValueError(f'Event {event_type} not recognized')
def run_scatter(proportions,event_type,event_timing,norm,ax,annotate=True,title=''):
    
    event_col = f'{event_timing}_WGD_{event_type}'
    if norm == 'Width' or norm == 'Width_Timing':
        event_col = f'{event_col}_Norm'
    if norm == 'Seg' or norm == 'Seg_Timing':
        event_col = f'{event_col}_Seg_Norm'
    
    proportion_val = proportions[event_col].values
    wgd_timing = proportions['WGD_Timing'].values
    if 'Timing' in norm:
        
        proportion_val = proportion_val[np.logical_and(wgd_timing>=0.1,wgd_timing<=0.9)]
        wgd_timing = wgd_timing[np.logical_and(wgd_timing>=0.1,wgd_timing<=0.9)]
        if 'Pre' in event_timing:
            proportion_val = proportion_val/wgd_timing
        if 'Post' in event_timing:
            proportion_val = proportion_val/(1-wgd_timing)
        proportion_val = np.log2(proportion_val)

    ax.scatter(wgd_timing,proportion_val,color='#070dba',alpha=0.5,s=5,rasterized=True)
    rho,p = spearmanr(wgd_timing,proportion_val)
    if annotate:
        ax.annotate(f'$\\rho$={rho:.2f}, p={format_p_value(p)}',xy=(0.5,0.5),xytext=(0.5,0.9),xycoords='axes fraction',textcoords='axes fraction',ha='center',va='center',fontsize=20)
    ax.set_xlabel('WGD Timing')

    if norm == 'None':
        ax.set_ylabel(f'Proportion of Genome\n{event_past_tense(event_type)} {event_timing}-WGD')
    if norm == 'Seg':
        ax.set_ylabel(f'Proportion of {event_type} Events {event_timing}-WGD')
    if norm == 'Width':
        ax.set_ylabel(f'Proportion of {event_plural(event_type)} {event_timing}-WGD')
    if norm == 'Seg_Timing':
        ax.set_ylabel(f'{event_type} Events {event_timing}-WGD\nLog2 Normalized Proportion')
    if norm == 'Width_Timing':
        ax.set_ylabel(f'{event_plural(event_type)} {event_timing}-WGD \n Log2 Normalized proportion')
    if title != '':
        ax.set_title(title,pad=28)
def plot_proportions(proportions,event_type,norm,out_path):
    fig,axs = plt.subplots(2,1,figsize=(12,10))
    run_scatter(proportions,event_type,'Pre',norm,axs[0])
    run_scatter(proportions,event_type,'Post',norm,axs[1])
    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path,dpi=300)
    plt.close(fig)

def plot_proportions_post_only(proportions,event_type,norm,out_path):
    fig,axs = plt.subplots(1,1,figsize=(12,5))
    run_scatter(proportions,event_type,'Post',norm,axs)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path,dpi=300)
    plt.close(fig)

def plot_proportions_cancer_type(proportions,event_type,event_timing,norm,out_path,min_samples=20):
    cancer_type_counts = proportions.groupby('Cancer_Type').size()
    cancer_types = cancer_type_counts[cancer_type_counts>=min_samples].index
    cancer_types = sorted(list(cancer_types))
    proportions = proportions[proportions['Cancer_Type'].isin(cancer_types)]
    
    fig,axs = plt.subplots(5,5,figsize=(33,22))
    axs = axs.flatten()
    processed_count = 0
    for i,cancer_type in enumerate(cancer_types):
        cancer_proportions = proportions[proportions['Cancer_Type']==cancer_type]
        n_samples = cancer_proportions.shape[0]
        run_scatter(cancer_proportions,event_type,event_timing,norm,axs[i],annotate=True,title=f'{cancer_type} (n={n_samples})')
        processed_count+=1
    #hide extra axes
    for i in range(processed_count,len(axs)):
        axs[i].axis('off')
    fig.tight_layout()
    fig.savefig(out_path,dpi=300)
    plt.close(fig)

def load_wgd_timing():

    wgd_timing = DataTools.get_wgd_calling_info()
    wgd_timing = wgd_timing[wgd_timing['WGD_Status']]
    wgd_timing = wgd_timing[['Sample_ID','Timing']].rename(columns={'Timing':'WGD_Timing'})
    good_samples = DataTools.get_good_samples()
    wgd_timing = wgd_timing[wgd_timing['Sample_ID'].isin(good_samples)]

    classifications = DataTools.load_classifications()[['Sample_ID','Cancer_Type']]
    wgd_timing = wgd_timing.merge(classifications,on='Sample_ID')
    return wgd_timing

def get_sample_proportions(sample_table,event,genome_size,min_threshold=0.05):
    pre_event = f'Pre_WGD_{event}'
    post_event = f'Post_WGD_{event}'
    sample_sum = sample_table.groupby('Sample_ID').apply(lambda x: np.sum(np.multiply(x[pre_event]+x[post_event],x['Segment_End']-x['Segment_Start']))).reset_index().rename(columns={0:'Total_Event_Size'})
    sample_sum['Proportion'] = sample_sum['Total_Event_Size']/genome_size
    sample_sum = sample_sum[sample_sum['Proportion']>min_threshold]
    sample_table = sample_table.merge(sample_sum,on='Sample_ID',how='inner')


    sample_pre_proportions = sample_table.groupby('Sample_ID').apply(lambda x: np.sum(np.multiply(x['Segment_End']-x['Segment_Start'],x[pre_event]))/genome_size).reset_index().rename(columns={0:pre_event})
    sample_post_proportions = sample_table.groupby('Sample_ID').apply(lambda x: np.sum(np.multiply(x['Segment_End']-x['Segment_Start'],x[post_event]))/genome_size).reset_index().rename(columns={0:post_event})
    sample_pre_proportions_norm = sample_table.groupby('Sample_ID').apply(lambda x: (np.sum(np.multiply(x['Segment_End']-x['Segment_Start'],x[pre_event]))+10e6)/(np.sum(np.multiply(x['Segment_End']-x['Segment_Start'],x[post_event]+x[pre_event]))+10e6)).reset_index().rename(columns={0:f'{pre_event}_Norm'})
    sample_post_proportions_norm = sample_table.groupby('Sample_ID').apply(lambda x: (np.sum(np.multiply(x['Segment_End']-x['Segment_Start'],x[post_event]))+10e6)/(np.sum(np.multiply(x['Segment_End']-x['Segment_Start'],x[post_event]+x[pre_event]))+10e6)).reset_index().rename(columns={0:f'{post_event}_Norm'})

    sample_pre_proportions_seg_norm = sample_table.groupby('Sample_ID').apply(lambda x: (np.sum(x[pre_event])+1)/(np.sum(x[post_event]+x[pre_event])+2)).reset_index().rename(columns={0:f'{pre_event}_Seg_Norm'})
    sample_post_proportions_seg_norm = sample_table.groupby('Sample_ID').apply(lambda x: (np.sum(x[post_event])+1)/(np.sum(x[post_event]+x[pre_event])+2)).reset_index().rename(columns={0:f'{post_event}_Seg_Norm'})
    sample_proportions = sample_pre_proportions.merge(sample_post_proportions,on='Sample_ID').merge(sample_pre_proportions_norm,on='Sample_ID').merge(sample_post_proportions_norm,on='Sample_ID')
    sample_proportions = sample_proportions.merge(sample_pre_proportions_seg_norm,on='Sample_ID').merge(sample_post_proportions_seg_norm,on='Sample_ID')
   
    return sample_proportions

def plot_pre_wgd_gain_comparison(gain_proportions):
    for mode in ['absolute','sum']:
        
        penalty_true_table = gain_proportions[True,mode][['Sample_ID','Pre_WGD_Gain','Post_WGD_Gain']].rename(columns={'Pre_WGD_Gain':'Pre_WGD_Gain_True','Post_WGD_Gain':'Post_WGD_Gain_True'})
        penalty_false_table = gain_proportions[False,mode][['Sample_ID','Pre_WGD_Gain','Post_WGD_Gain']].rename(columns={'Pre_WGD_Gain':'Pre_WGD_Gain_False','Post_WGD_Gain':'Post_WGD_Gain_False'})
        comparison_table = penalty_true_table.merge(penalty_false_table,on='Sample_ID')
        for event_timing in ['Pre_WGD','Post_WGD']:
            out_path = f'../plots/wgd_timing_event_proportion/Gain/comparison_{event_timing}_{mode}.pdf'
            fig,ax = plt.subplots(1,1,figsize=(10,10))
            ax.scatter(comparison_table[f'{event_timing}_Gain_False'],comparison_table[f'{event_timing}_Gain_True'],color='#070dba',alpha=0.5,s=5,rasterized=True)
            correlation,p_value = pearsonr(comparison_table[f'{event_timing}_Gain_False'],comparison_table[f'{event_timing}_Gain_True'])
            max_proportion = max(comparison_table[f'{event_timing}_Gain_False'].max(),comparison_table[f'{event_timing}_Gain_True'].max())
            ax.plot([0,max_proportion],[0,max_proportion],color='red',linestyle='--',linewidth=3)

            ax.set_xlabel(f"Proportion of Genome Gained {event_timing.replace('_',' ')} Without Non-Parsimony Penalty")
            ax.set_ylabel(f"Proportion of Genome Gained {event_timing.replace('_',' ')} With Non-Parsimony Penalty")
            ax.set_title(f'Pearson Correlation: {correlation:.2f}, p={format_p_value(p_value)}')
            fig.tight_layout()
            fig.savefig(out_path,dpi=300)
            plt.close(fig)
if __name__ == '__main__':
    output_dir = '../plots/wgd_timing_event_proportion'
    for event in ['Gain','Loss']:
        os.makedirs(f'{output_dir}/{event}',exist_ok=True)
    chromosome_table = load_chromosome_table()
    genome_size = (chromosome_table['Chrom_End']-chromosome_table['Chrom_Start']).sum()

    loss_timing_data = load_loss_timing_table()
    loss_timing_proportions = get_sample_proportions(loss_timing_data,'Loss',genome_size)
    
    gain_timing_proportions_store = {} 
    for mode in ['absolute','sum']:
        for apply_penalty in [True,False]:
            gain_timing_table = load_gain_timing_table(apply_penalty,'All',mode)
            gain_timing_proportions_store[apply_penalty,mode] = get_sample_proportions(gain_timing_table,'Gain',genome_size)
    
    plot_pre_wgd_gain_comparison(gain_timing_proportions_store)
    wgd_timing_store = {True:load_wgd_timing(),False:load_wgd_timing()}
    
    

    loss_plotting_proportions = loss_timing_proportions.merge(wgd_timing_store[False],how='inner')
    
    for norm in ['None','Seg','Width','Seg_Timing','Width_Timing']:
        plot_proportions(loss_plotting_proportions,'Loss',norm,f'{output_dir}/Loss/proportions_wgd_timing_{norm}.pdf')
        plot_proportions_post_only(loss_plotting_proportions,'Loss',norm,f'{output_dir}/Loss/proportions_pre_wgd_timing_{norm}.pdf')
        plot_proportions_cancer_type(loss_plotting_proportions,'Loss','Pre',norm,f'{output_dir}/Loss/proportions_post_wgd_timing_cancer_type_norm_{norm}.pdf')

        plot_proportions_cancer_type(loss_plotting_proportions,'Loss','Post',norm,f'{output_dir}/Loss/proportions_post_wgd_timing_cancer_type_norm_{norm}.pdf')
        for apply_penalty in  [True,False]:
            for mode in ['absolute','sum']:
                plotting_proportions = gain_timing_proportions_store[apply_penalty,mode].merge(wgd_timing_store[apply_penalty],how='inner')
                
                plot_proportions(plotting_proportions,'Gain',norm,f'{output_dir}/Gain/proportions_wgd_timing_penalty_{apply_penalty}_norm_{norm}_mode_{mode}.pdf')
                plot_proportions_post_only(plotting_proportions,'Gain',norm,f'{output_dir}/Gain/proportions_post_wgd_timing_penalty_{apply_penalty}_norm_{norm}_mode_{mode}.pdf')
                plot_proportions_cancer_type(plotting_proportions,'Gain','Pre',norm,f'{output_dir}/Gain/proportions_pre_wgd_timing_penalty_{apply_penalty}_cancer_type_norm_{norm}_mode_{mode}.pdf')
                plot_proportions_cancer_type(plotting_proportions,'Gain','Post',norm,f'{output_dir}/Gain/proportions_post_wgd_timing_penalty_{apply_penalty}_cancer_type_norm_{norm}_mode_{mode}.pdf')


        

