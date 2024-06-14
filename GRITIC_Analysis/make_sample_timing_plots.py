import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
#font size 18, pdf font type 42, nimbus sans font
plt.rcParams.update({'font.size': 18, 'pdf.fonttype': 42, 'font.family': 'Nimbus Sans'})

import DataTools

def get_chrom_abs_positions():
    chrom_arm_positions = pd.read_csv('../resources/chrom_arm_positions.tsv',sep="\t")
    chrom_positions = chrom_arm_positions.groupby('Chromosome').agg({'Arm_End':'max'}).reset_index().rename(columns={'Arm_End':'Chrom_Size'})
    chrom_order = list(map(str,range(1,23))) + ['X','Y']
    chrom_positions['Chromosome'] = pd.Categorical(chrom_positions['Chromosome'],chrom_order)
    chrom_positions = chrom_positions.sort_values(by='Chromosome')
    abs_start = np.cumsum(chrom_positions['Chrom_Size'].values)
    abs_start = np.insert(abs_start,0,0)

    chrom_positions['Chrom_Abs_Start'] = abs_start[:-1]
    chrom_positions['Chrom_Abs_End'] = abs_start[1:]
    return chrom_positions

def add_gain_index(sample_table):
    
    return sample_table
def plot_sample_table(sample_table,apply_penalty):
    gain_index_colors = {'1st':'#ffeda0','2nd':'#feb24c','3+':'#f03b20'}
    sample_table = sample_table.merge(chrom_abs_positions,on='Chromosome')
    sample_table['Segment_Start_Abs'] = sample_table['Chrom_Abs_Start'] + sample_table['Segment_Start']
    sample_table['Segment_End_Abs'] = sample_table['Chrom_Abs_Start'] + sample_table['Segment_End']
    
    fig,ax = plt.subplots(figsize=(10,5))
    for index,row in sample_table.iterrows():
        ax.plot([row['Segment_Start_Abs'],row['Segment_End_Abs']],[row['Gain_Timing'],row['Gain_Timing']],color=gain_index_colors[row['Gain_Index_Str']],lw=3,alpha=0.5)
    
    wgd_timing = sample_table['WGD_Timing'].values[0]
    if not pd.isnull(wgd_timing):
        ax.plot([0,chrom_abs_positions['Chrom_Abs_End'].values[-1]],[wgd_timing,wgd_timing],color='red',lw=3,linestyle='--')
    ax.set_xticks(chrom_plotting_positions['Chrom_Abs_Start'].values)
    ax.set_xticklabels(chrom_plotting_positions['Chromosome'].values)
    #make legend
    legend_handles = []
    legend_handle_names = []
    for gain_index in ['1st','2nd','3+']:
        if gain_index in sample_table['Gain_Index_Str'].values:
            legend_handles.append(plt.Line2D([0],[0],color=gain_index_colors[gain_index],lw=3))
            legend_handle_names.append(gain_index)
    if not pd.isnull(wgd_timing):
        legend_handles.append(plt.Line2D([0],[0],color='red',lw=3,linestyle='--'))
        legend_handle_names.append('WGD')
    ax.legend(legend_handles,legend_handle_names,loc='upper left',prop={'size': 12})
    ax.set_ylim(-0.1,1.1)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Gain Timing')
    tumor_type = 'Primary' if sample_table['Dataset'].values[0] =='PCAWG' else 'Metastatic'
    plt_title = f'{sample_table["Cancer_Type"].values[0]} - {tumor_type}'
    plt.title(plt_title)
    out_path = f'../plots/sample_timing_plots/prior_{apply_penalty}/{sample_table["Sample_ID"].values[0]}_sample_timing.pdf'
    print(out_path)
    fig.tight_layout()
    plt.savefig(out_path,bbox_inches='tight')
    plt.savefig(out_path.replace('.pdf','.png'),bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    chrom_abs_positions = get_chrom_abs_positions()
    chrom_plotting_positions = chrom_abs_positions[~chrom_abs_positions['Chromosome'].isin(['18','20','22'])].copy()
    sample_index = int(sys.argv[1])
    for apply_penalty in [True,False]:
        posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty)
        posterior_table['Gain_Index_Str'] = posterior_table['Gain_Index'].map(lambda x: {0:'1st',1:'2nd',2:'3+'}[min(2,x)])

        sample_ids = sorted(posterior_table['Sample_ID'].unique())
    
        sample_table = posterior_table[posterior_table['Sample_ID']==sample_ids[sample_index]].copy()
        sample_table = add_gain_index(sample_table)
        plot_sample_table(sample_table,apply_penalty)
        print(f'Finished {sample_ids[sample_index]}, {apply_penalty}')



