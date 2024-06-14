import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#font size 18,pdf font type 42, Nimbus Sans Font
plt.rcParams.update({'font.size': 18, 'pdf.fonttype': 42, 'font.family': 'Nimbus Sans'})

import DataTools

from scipy.stats import pearsonr

RNG = np.random.default_rng()
def load_posterior_table(apply_penalty,proportion_cut_off = 0.5):
    posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty)
    
    

    route_counts = posterior_table[['Sample_ID','Segment_ID','Posterior_Sample_Index','Route']].drop_duplicates()
    route_counts = route_counts.groupby(['Sample_ID','Segment_ID'])['Route'].value_counts(normalize=True).reset_index(name='Proportion')
    max_counts = route_counts.groupby(['Sample_ID','Segment_ID'])['Proportion'].max().reset_index()
    max_counts = max_counts[max_counts['Proportion'] < proportion_cut_off][['Sample_ID','Segment_ID']]
    #max_counts = max_counts.sample(frac=0.1)
    posterior_table = posterior_table.merge(max_counts, on=['Sample_ID','Segment_ID'], how='inner')
    return posterior_table

def generate_posterior_route_samples(route_map):
    if len(set(route_map.values())) == 1:
        raise ValueError('Only one route found')
    while True:
        posterior_index_one,posterior_index_two = RNG.choice(list(route_map.keys()),2)
        if route_map[posterior_index_one] != route_map[posterior_index_two]:
            return posterior_index_one,posterior_index_two

def get_all_gains(sample_table):
    wgd_status = ~np.isnan(sample_table['WGD_Timing'].iloc[0])
    n_gains = (sample_table['Major_CN'].iloc[0]-1) + max(0,sample_table['Minor_CN'].iloc[0]-1)
    independent_gains = sample_table['Gain_Timing'].values
    if not wgd_status:
        assert n_gains == len(independent_gains)
        return np.sort(independent_gains),np.ones(n_gains).astype(bool)
    
    if n_gains == len(independent_gains):
        return np.sort(independent_gains),np.ones(n_gains).astype(bool)
    wgd_gains = [sample_table['WGD_Timing'].iloc[0]]*(n_gains - len(independent_gains))
    gain_timing = np.concatenate([independent_gains,wgd_gains])
    is_independent_gain = np.concatenate([np.ones(len(independent_gains)),np.zeros(n_gains - len(independent_gains))])
    sort_indices = np.argsort(gain_timing)
    gain_timing = gain_timing[sort_indices]
    is_independent_gain = is_independent_gain[sort_indices].astype(bool)
    return gain_timing,is_independent_gain
def generate_segment_comparison(segment_table,mode):
    assert mode in ['first','last','all']
    route_map = segment_table[['Posterior_Sample_Index','Route']].drop_duplicates()
    route_map  = {posterior_sample_index:route for posterior_sample_index,route in zip(route_map['Posterior_Sample_Index'],route_map['Route'])}
    posterior_index_one,posterior_index_two = generate_posterior_route_samples(route_map)

    sample_one = segment_table[segment_table['Posterior_Sample_Index'] == posterior_index_one]
    sample_two = segment_table[segment_table['Posterior_Sample_Index'] == posterior_index_two]

    if mode == 'first':
        return [sample_one['Gain_Timing'].min()], [sample_two['Gain_Timing'].min()]
    if mode == 'last':
        return [sample_one['Gain_Timing'].max()], [sample_two['Gain_Timing'].max()]
    sample_one_gains,sample_one_independent = get_all_gains(sample_one)
    sample_two_gains,sample_two_independent = get_all_gains(sample_two)
    either_independent = sample_one_independent | sample_two_independent
    sample_one_gains = sample_one_gains[either_independent]
    sample_two_gains = sample_two_gains[either_independent]
    return sample_one_gains,sample_two_gains
    

def generate_comparison_table(posterior_table,mode):
    comparison_data = {'Sample_ID':[],'Segment_ID':[],'Sample_One_Gain_Timing':[],'Sample_Two_Gain_Timing':[]}
    for (sample_id,segment_id),segment_table in posterior_table.groupby(['Sample_ID','Segment_ID']):
        
        segment_comparison_one,segment_comparison_two = generate_segment_comparison(segment_table,mode=mode)
        assert len(segment_comparison_one) == len(segment_comparison_two) 
        comparison_data['Sample_ID'].extend([sample_id]*len(segment_comparison_one))
        comparison_data['Segment_ID'].extend([segment_id]*len(segment_comparison_one))
        comparison_data['Sample_One_Gain_Timing'].extend(segment_comparison_one)
        comparison_data['Sample_Two_Gain_Timing'].extend(segment_comparison_two)
    return pd.DataFrame(comparison_data)

out_dir = '../plots/route_difference'

for apply_penalty in [True,False]:
    posterior_table = load_posterior_table(apply_penalty)

    mode_title_mapping = {'all':'All Gains','first':'First Gain','last':'Last Gain'}
    for min_n_mutations in [0,100]:
        posterior_table_mutation = posterior_table[posterior_table['N_Mutations'] > min_n_mutations].copy()
        fig,ax = plt.subplots(1,3,figsize=(20,10))

        for i,mode in enumerate(['all','first','last']):
            comparison_table = generate_comparison_table(posterior_table_mutation,mode=mode)
            ax[i].scatter(comparison_table['Sample_One_Gain_Timing'],comparison_table['Sample_Two_Gain_Timing'],s=30,alpha=0.05,color='blue',rasterized=True)
            ax[i].set_xlabel('Route A Gain Timing')
            ax[i].set_ylabel('Route B Gain Timing')
            ax[i].set_title(mode_title_mapping[mode])
            ax[i].plot([0,1],[0,1],'k--')
            r,p = pearsonr(comparison_table['Sample_One_Gain_Timing'],comparison_table['Sample_Two_Gain_Timing'])
            ax[i].text(0.1,0.9,f'r={r:.2f}, p={p:.2e}',transform=ax[i].transAxes)
        plt.tight_layout()
        plt.savefig(f'{out_dir}/route_difference_{min_n_mutations}_mutations_apply_penalty_{apply_penalty}.pdf')
        plt.savefig(f'{out_dir}/route_difference_{min_n_mutations}_mutations_apply_penalty_{apply_penalty}.png',dpi=300)

