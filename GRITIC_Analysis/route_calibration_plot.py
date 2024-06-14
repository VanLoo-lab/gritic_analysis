import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Sans'

import matplotlib.pyplot as plt
#font size 22
plt.rcParams.update({'font.size': 22})
import matplotlib.patches as patches
import DataTools
import itertools
RNG = np.random.default_rng(42)
def load_route_data(parsimony_only,parsimony_penalty,parsimony_routes):
    print(parsimony_only)
    if parsimony_only:
        dataset = 'state_validation_parsimony_only'
    else:
        dataset = 'state_validation_all_routes'

    route_table_path = f'../../output/{dataset}/complete_route_table_{dataset}.tsv'  
    metadata_path = f'../../output/{dataset}/complete_metadata_table_{dataset}.tsv' 
    route_table = pd.read_csv(route_table_path, sep='\t',dtype={'Chromosome':str})
    route_table = route_table[['Sample_ID','Segment_ID','WGD_Status','Major_CN','Minor_CN','Route','N_Mutations','Probability','Average_N_Events','Average_Pre_WGD_Losses','Average_Post_WGD_Losses']].drop_duplicates()
    route_table = route_table[route_table['WGD_Status']]
    route_table = route_table[route_table['N_Mutations']>=20]
 
    metadata = pd.read_csv(metadata_path, sep='\t')
    metadata = metadata[['Sample_ID','Segment_ID','Route']].drop_duplicates()
    metadata['Route'] = metadata['Route'].str.slice(0,9)
    metadata['Segment_ID'] = metadata['Segment_ID'].str.replace('Segment_','')
    metadata = metadata.rename(columns={'Route':'True_Route'})
  
    route_table = route_table.merge(metadata, on=['Sample_ID','Segment_ID'], how='inner')
    
    route_table['Correct'] = route_table['Route'] == route_table['True_Route']
    
    if parsimony_penalty:
        route_table.loc[:,'Probability'] = DataTools.apply_prior_quick(route_table,2.7)
    if parsimony_routes:
        wgd_routes = pd.read_csv('../../Analysis/output/all_routes_n_events.tsv',sep="\t")
        wgd_routes.loc[:,'Route'] = wgd_routes.loc[:,'Route'].str.slice(0,9)
        parsimony_routes = wgd_routes[wgd_routes['Parsimonious'] & (wgd_routes['WGD_Status'])]
        route_table = route_table[route_table['Route'].isin(parsimony_routes['Route'])]
    return route_table

def get_bin_data(route_table, bins):
    bin_data = {'Major_CN':[],'Minor_CN':[],'Bin_Low':[], 'Bin_High':[], 'N_Segments':[],'N_Segments_Max':[], 'Proportion_Correct':[]}
    #filter for max probability per Sample_ID, Segment_ID
    route_table = route_table.copy()
    route_table.loc[:,'Probability_Bin'] = pd.cut(np.minimum(np.maximum(route_table['Probability'],1e-6),1-1e-6), bins=bins, labels=False)
   
    for (major_cn,minor_cn),cn_table in route_table.groupby(['Major_CN','Minor_CN']):
        max_probability_table = cn_table.groupby(['Sample_ID','Segment_ID'])['Probability'].max().reset_index()
        for i,bin_table in cn_table.groupby('Probability_Bin'):
            max_probability_bin_table = max_probability_table[(max_probability_table['Probability'] >= bins[i]) & (max_probability_table['Probability'] < bins[i+1])]
            
            bin_data['Major_CN'].append(major_cn)
            bin_data['Minor_CN'].append(minor_cn)
            bin_data['Bin_Low'].append(bins[i])
            bin_data['Bin_High'].append(bins[i+1])
            bin_data['N_Segments'].append(len(bin_table.index))
            bin_data['N_Segments_Max'].append(len(max_probability_bin_table.index))
            
            bin_data['Proportion_Correct'].append(np.mean(bin_table['Correct']))
    bin_data = pd.DataFrame(bin_data)
    bin_data.loc[:,'Bin_ID'] = bin_data['Bin_Low'].round(2).astype(str) + '-' + bin_data['Bin_High'].round(2).astype(str)
    return pd.DataFrame(bin_data)

def plot_calibrations_axis(bin_data,bins,ax,plot_error_bars=True,legend=True):
    
    bin_data = bin_data.copy()
    bin_data = bin_data[bin_data['N_Segments']>=100]
    bin_data.loc[:,'Low_CI_Width'] = bin_data['Proportion_Correct'] - bin_data['Low_CI']
    bin_data.loc[:,'High_CI_Width'] = bin_data['High_CI'] - bin_data['Proportion_Correct']
    print(bin_data)
    ax.plot([0,1],[0,1],color='red',linestyle='--',label='Perfect Calibration',lw=5)
    ax.scatter((bin_data['Bin_Low']+bin_data['Bin_High'])/2, bin_data['Proportion_Correct'],s=120,color='#64a4ed')
    if plot_error_bars:
        #add error bars
        ax.errorbar((bin_data['Bin_Low']+bin_data['Bin_High'])/2, bin_data['Proportion_Correct'], yerr=[bin_data['Low_CI_Width'],bin_data['High_CI_Width']], fmt='none', ecolor='black', elinewidth=2, capsize=14)

    if legend:
        ax.legend()
        ax.set_xlabel('Measured Probability')
        ax.set_ylabel('True Probability')

def plot_calibrations_cn(bin_data,bins,parsimony_only,parsimony_penalty,parsimony_routes,plot_error_bars=True):
    bin_data = bin_data.copy()
    bin_data.loc[:,'CN_State'] = bin_data['Major_CN'].astype(str) + '+' + bin_data['Minor_CN'].astype(str)
    
    #two rows, one column shared x axis
    fig, axs = plt.subplots(6,5, sharex=True, figsize=(32,26))
    axs = axs.flatten()
    processed_count =0
    for cn_state,cn_state_table in bin_data.groupby('CN_State'):
        ax = axs[processed_count]
        plot_calibrations_axis(cn_state_table,bins,ax,plot_error_bars=plot_error_bars,legend=False)
        ax.set_title(cn_state)
        processed_count +=1

    plt.tight_layout()
    #hide last four axes
    for i in range(len(axs)-4,len(axs)):
        axs[i].axis('off')
    output_path = f'../plots/route_calibrations/calibration_curve_cn_parsimony_only_{parsimony_only}_parsimony_penalty_{parsimony_penalty}_parsimony_routes_{parsimony_routes}.png'
    plt.savefig(output_path)
    plt.savefig(output_path.replace('.png','.pdf'))
def plot_calibrations_combined(bin_data,bins,parsimony_only,parsimony_penalty,parsimony_routes,plot_error_bars=True):

    #two rows, one column shared x axis
    fig, axs = plt.subplots(2,1, figsize=(10,10))
    plot_calibrations_axis(bin_data,bins,axs[0],plot_error_bars=plot_error_bars)
    axs[1].bar(bin_data['Bin_ID'],bin_data['N_Segments_Max'],color='#64a4ed')
    axs[1].set_xticklabels(bin_data['Bin_ID'],rotation=45)
    axs[1].set_ylabel('Number of Segments')
    axs[1].set_yscale('log')
    plt.tight_layout()
    output_path = f'../plots/route_calibrations/calibration_curve_combined_parsimony_only_{parsimony_only}_parsimony_penalty_{parsimony_penalty}_parsimony_routes_{parsimony_routes}.png'
    plt.savefig(output_path)
    plt.savefig(output_path.replace('.png','.pdf'))
    print(output_path)


def generate_bootstrapped_bin_data(route_table, bins,n_bootstraps):
    sample_dict = DataTools.get_sample_dict(route_table)
    bootstrap_bin_data = []
    for i in range(n_bootstraps):
        print(f'Bootstrap {i}')
        bootstrap_sample_table = DataTools.get_bootstrap_table(sample_dict)
        bootstrap_sample = get_bin_data(bootstrap_sample_table, bins)
        bootstrap_sample.loc[:,'Bootstrap'] = i
        bootstrap_bin_data.append(bootstrap_sample)
    bootstrap_bin_data = pd.concat(bootstrap_bin_data)
    return bootstrap_bin_data
def get_bootstrapped_bin_data(route_table, bins,parsimony_only,parsimony_penalty, parsimony_routes,n_bootstraps,overwrite=True):
    output_path = f'../output/route_calibrations/bootstrap_bin_data_{n_bootstraps}_parsimony_only_{parsimony_only}_parsimony_penalty_{parsimony_penalty}_parsimony_routes_{parsimony_routes}.tsv'
    if os.path.exists(output_path) and not overwrite:
        bootstrap_bin_data = pd.read_csv(output_path, sep='\t')
    else:
        bootstrap_bin_data = generate_bootstrapped_bin_data(route_table, bins, n_bootstraps)
        bootstrap_bin_data.to_csv(output_path, sep='\t', index=False)
    return bootstrap_bin_data

def summarise_bootstraps(bootstrap_bin_data,use_cn=False):
    if use_cn:
        group_cols = ['Major_CN','Minor_CN','Bin_ID']
    else:
        group_cols = ['Bin_ID']
    low_ci_columns = [*group_cols,'Low_CI']
    high_ci_columns = [*group_cols,'High_CI']


    low_ci = bootstrap_bin_data.groupby(group_cols)['Proportion_Correct'].quantile(0.05/2).reset_index()
    low_ci.columns = low_ci_columns
    high_ci = bootstrap_bin_data.groupby(group_cols)['Proportion_Correct'].quantile(1-0.05/2).reset_index()
    high_ci.columns = high_ci_columns
    
    
    ci = pd.merge(low_ci, high_ci)
    return ci


def combine_major_cn(bin_data,use_bootstraps):
    
    bin_data_summary = bin_data[['Bin_ID','Bin_Low','Bin_High']].drop_duplicates()
    if not use_bootstraps:
        bin_data_combined = bin_data.groupby(['Bin_ID']).apply(lambda x: np.average(x['Proportion_Correct'], weights=x['N_Segments'])).reset_index()
        n_segments = bin_data.groupby(['Bin_ID'])[['N_Segments','N_Segments_Max']].sum().reset_index()
        bin_data_combined = bin_data_combined.merge(n_segments, on='Bin_ID')
        bin_data_combined.columns = ['Bin_ID','Proportion_Correct','N_Segments','N_Segments_Max']
    else:
        bin_data_combined = bin_data.groupby(['Bin_ID','Bootstrap']).apply(lambda x: np.average(x['Proportion_Correct'], weights=x['N_Segments'])).reset_index()
        n_segments = bin_data.groupby(['Bin_ID','Bootstrap'])[['N_Segments','N_Segments_Max']].sum().reset_index()
        bin_data_combined = bin_data_combined.merge(n_segments, on=['Bin_ID','Bootstrap'])
        bin_data_combined.columns = ['Bin_ID','Bootstrap','Proportion_Correct','N_Segments','N_Segments_Max']
    return bin_data_summary.merge(bin_data_combined, on='Bin_ID')

#all possible combinations of 3xTrue False
run_params = list(itertools.product([False, True], repeat=3))

run_number = int(sys.argv[1])

parsimony_only, parsimony_penalty,parsimony_routes = run_params[run_number]
print('PARAMS')
print('Parsimony only: ', parsimony_only)
print('Parsimony penalty: ', parsimony_penalty)
print('Parsimony routes: ', parsimony_routes)

route_table = load_route_data(parsimony_only=parsimony_only, parsimony_penalty=parsimony_penalty,parsimony_routes=parsimony_routes)
route_table =route_table[route_table['Major_CN']>2]
route_table['CN_State'] = route_table['Major_CN'].astype(str) + route_table['Minor_CN'].astype(str)

bins = np.linspace(0,1,21)
small_bins = np.linspace(0,1,11)

bin_data = get_bin_data(route_table, bins)

bin_data_major_cn_combined = combine_major_cn(bin_data,use_bootstraps=False)
bin_data_major_cn_combined.to_csv(f'../output/route_calibrations/bin_data_major_cn_combined_parsimony_only_{parsimony_only}_parsimony_penalty_{parsimony_penalty}_parsimony_routes_{parsimony_routes}.tsv',sep='\t',index=False)

bootstrapped_bin_data = get_bootstrapped_bin_data(route_table, bins,parsimony_only,parsimony_penalty,parsimony_routes,n_bootstraps=250,overwrite=True)

bootstrapped_bin_data_major_cn_combined = combine_major_cn(bootstrapped_bin_data,use_bootstraps=True)
bootstrapped_bin_data = summarise_bootstraps(bootstrapped_bin_data,use_cn=True)

bootstrapped_bin_data_major_cn_combined = summarise_bootstraps(bootstrapped_bin_data_major_cn_combined,use_cn=False)



bin_data = bin_data.merge(bootstrapped_bin_data)

bin_data_major_cn_combined = bin_data_major_cn_combined.merge(bootstrapped_bin_data_major_cn_combined, on=['Bin_ID'])


plot_calibrations_combined(bin_data_major_cn_combined,bins,parsimony_only,parsimony_penalty,parsimony_routes,plot_error_bars=True)
plot_calibrations_cn(bin_data,bins,parsimony_only,parsimony_penalty,parsimony_routes,plot_error_bars=True)


