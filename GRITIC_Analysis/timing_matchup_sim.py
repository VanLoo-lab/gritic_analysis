import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Sans'
import matplotlib.pyplot as plt
#font size 18
matplotlib.rcParams.update({'font.size': 18})
import argparse
import os
from pathlib import Path
from scipy.stats import pearsonr
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_segment_shuffle(metadata_table):
    sample_table = metadata_table[['Sample_ID','Major_CN','Minor_CN','Segment_ID']].copy().drop_duplicates()
    shuffle_store = []
    for cn_state,cn_state_table in sample_table.groupby(['Major_CN','Minor_CN']):

        cn_state_table_shuffle = cn_state_table.sample(frac=1,replace=False).copy().reset_index(drop=True)
        cn_state_table['New_Segment_ID'] = cn_state_table_shuffle['Segment_ID'].values
        cn_state_table['New_Sample_ID'] = cn_state_table_shuffle['Sample_ID'].values
 
        shuffle_store.append(cn_state_table)
    shuffle_table = pd.concat(shuffle_store,axis=0)
    return shuffle_table



def load_metadata_table(input_dir,dataset,shuffle=False):
    metadata_table = pd.read_csv(input_dir /f'complete_metadata_table_{dataset}.tsv',sep="\t")
    metadata_table['Segment_ID'] = metadata_table['Segment_ID'].str.replace('Segment_','')
    metadata_table['CN_State'] = metadata_table['Major_CN'].astype(str) + '_' + metadata_table['Minor_CN'].astype(str)
    metadata_table['N_Gains'] = np.maximum(0,metadata_table['Major_CN']-1) + np.maximum(0,metadata_table['Minor_CN']-1)
    #extract real gain timing from metadata table and convert to array
    metadata_table['Real_Gain_Timing'] = metadata_table['Gain_Timing'].apply(lambda x: np.array(list(eval(x).values())))
    metadata_table['Real_Gain_Timing'] = metadata_table.apply(lambda x: np.append(x['Real_Gain_Timing'],np.repeat(x['WGD_Timing'],x['N_Gains']-x['Real_Gain_Timing'].size)),axis=1)
    metadata_table['Real_Gain_Timing'] = metadata_table['Real_Gain_Timing'].apply(lambda x: np.sort(x))
    metadata_table = metadata_table[['Sample_ID','CN_State','Major_CN','Minor_CN','Segment_ID','Real_Gain_Timing','WGD_Timing']].rename(columns={'WGD_Timing':'Real_WGD_Timing'})
    if shuffle:
        segment_shuffle = get_segment_shuffle(metadata_table)
        metadata_table = metadata_table.merge(segment_shuffle,how='inner',on=['Sample_ID','Major_CN','Minor_CN','Segment_ID'])
        metadata_table =metadata_table.drop(columns=['Segment_ID','Sample_ID']).rename(columns={'New_Sample_ID':'Sample_ID','New_Segment_ID':'Segment_ID'})

    return metadata_table

def load_posterior_table(input_dir,dataset,apply_penalty):
    posterior_table = pd.read_csv(input_dir /f'complete_posterior_table_penalty_{apply_penalty}_{dataset}.tsv',sep="\t",dtype={'Chromosome':str})
    posterior_table = posterior_table[~posterior_table['WGD_Timing'].isnull()]
    posterior_table = posterior_table[posterior_table['N_Mutations']>=20]
    posterior_table['N_Gains'] = np.maximum(0,posterior_table['Major_CN']-1) + np.maximum(0,posterior_table['Minor_CN']-1)
    posterior_table_measured_timing = posterior_table.groupby(['Sample_ID','Segment_ID','Posterior_Sample_Index']).apply(lambda x: np.append(x['Gain_Timing'].values,np.repeat(x['WGD_Timing'].values[0],x['N_Gains'].values[0]-x['Gain_Timing'].values.size)))
    posterior_table_measured_timing = posterior_table_measured_timing.apply(lambda x: np.sort(x))
    posterior_table_measured_timing = posterior_table_measured_timing.reset_index().rename(columns={0:'Measured_Gain_Timing'})
    return posterior_table_measured_timing


my_parser = argparse.ArgumentParser(allow_abbrev=False)
my_parser.add_argument('--dataset', action='store', type=str, required=True)
my_parser.add_argument('--apply_penalty', action='store', type=str2bool, required=True)
args = my_parser.parse_args()
dataset = args.dataset
apply_penalty = args.apply_penalty

input_dir = Path(f'../in_data/{dataset}/')

plot_dir = Path(f'../plots/state_matchup_histograms/{dataset}_penalty_{apply_penalty}/')

os.makedirs(plot_dir,exist_ok=True)
metadata_table = load_metadata_table(input_dir,dataset,shuffle=False)


posterior_table = load_posterior_table(input_dir,dataset,apply_penalty)

timing_table = metadata_table.merge(posterior_table,how='inner')



timing_table['Measured_Gain_Timing_WGD_Filtered'] = timing_table.apply(lambda x: x['Measured_Gain_Timing'][~np.isclose(x['Real_WGD_Timing'],x['Real_Gain_Timing'])],axis=1)
timing_table['Real_Gain_Timing_WGD_Filtered'] = timing_table.apply(lambda x: x['Real_Gain_Timing'][~np.isclose(x['Real_WGD_Timing'],x['Real_Gain_Timing'])],axis=1)


combined_hist_store = []
fig,axs = plt.subplots(6,5,figsize=(20,20))
axs = axs.flatten()
plt_count = 0
for cn_state,cn_state_table in timing_table.groupby('CN_State'):
    true_array = np.concatenate(cn_state_table['Real_Gain_Timing_WGD_Filtered'].values)
    measured_array = np.concatenate(cn_state_table['Measured_Gain_Timing_WGD_Filtered'].values)
    r,p = pearsonr(true_array,measured_array)
    bins = np.linspace(0,1,101)
    hist2d, xedges, yedges = np.histogram2d(true_array, measured_array, bins=bins)
    hist2d = hist2d.T
    combined_hist_store.append(hist2d)
    hist2d_norm = hist2d/ np.max(hist2d,axis=0,keepdims=True)

    axs[plt_count].imshow(hist2d_norm,origin='lower',extent=[0,1,0,1])
    #axs[plt_count].set_title(cn_state.replace('_','+'))
    axs[plt_count].set_title(f'{cn_state.replace("_","+")}')
    plt_count += 1
#hide empty plots
for i in range(plt_count,len(axs)):
    axs[i].axis('off')
fig.supxlabel('True Timing')
fig.supylabel('Measured Timing')
plt.tight_layout()
plt.savefig(plot_dir /f'match_up_cn_states_separate.png',dpi=300)
plt.savefig(plot_dir /f'match_up_cn_states_separate.pdf',dpi=300)
plt.close()

combined_hist = np.sum(combined_hist_store,axis=0)
combined_hist_norm = combined_hist/ np.max(combined_hist,axis=0,keepdims=True)

fig,ax = plt.subplots(1,1,figsize=(10,10))
im = ax.imshow(combined_hist_norm,origin='lower',extent=[0,1,0,1])
ax.set_title('All CN States')
ax.set_xlabel('True Timing')
ax.set_ylabel('Measured Timing')
#colorbar
cbar = plt.colorbar(im,ax=ax)
cbar.set_label('Normalized Posterior Probability')
plt.savefig(plot_dir /f'match_up_cn_states_combined.png',dpi=300)
plt.savefig(plot_dir /f'match_up_cn_states_combined.pdf',dpi=300)