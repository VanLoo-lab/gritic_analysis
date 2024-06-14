import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
#font size 18,pdf font type 42, font nimbus sans
plt.rcParams.update({'font.size': 12, 'pdf.fonttype': 42, 'font.family': 'Nimbus Sans'})

import DataTools

import itertools

#import fisher exact test
from scipy.stats import fisher_exact
RNG = np.random.default_rng()
def load_pre_post_table(apply_penalty):
    #wgd_timing_table = DataTools.load_route_table(apply_penalty=False)[['Sample_ID','WGD_Timing']].drop_duplicates()
    #early_wgd = wgd_timing_table[wgd_timing_table['WGD_Timing']>0.7]
    pre_post_table_path =  f'../output/pre_post_non/pre_post_non_gains_apply_penalty_{apply_penalty}_major_cn_All_mode_absolute.tsv'
    pre_post_table = pd.read_csv(pre_post_table_path, sep='\t')
    pre_post_table['Chromosome'] = pre_post_table['Segment_ID'].str.split('-').str[0]
    pre_post_table['Segment_Start'] = pre_post_table['Segment_ID'].str.split('-').str[1].astype(int)
    pre_post_table['Segment_End'] = pre_post_table['Segment_ID'].str.split('-').str[2].astype(int)
    
    pre_post_table['Segment_Width'] = pre_post_table['Segment_End']-pre_post_table['Segment_Start']
    pre_post_table = pre_post_table[pre_post_table['Non_Probability']<1e-6]
    pre_post_table = pre_post_table[pre_post_table['Major_CN'].isin([3,4])]
    
    
    pre_post_table_pre_wgd = pre_post_table.groupby(['Sample_ID','Chromosome','Major_CN']).apply(lambda x: np.average(x['Pre_Probability'],weights=x['Segment_Width'])).reset_index().rename(columns={0:'Pre_Probability'})
    pre_post_table_post_wgd = pre_post_table.groupby(['Sample_ID','Chromosome','Major_CN']).apply(lambda x: np.average(x['Post_Probability'],weights=x['Segment_Width'])).reset_index().rename(columns={0:'Post_Probability'})
    pre_post_table = pre_post_table_pre_wgd.merge(pre_post_table_post_wgd,on=['Sample_ID','Chromosome','Major_CN'])
    
    
    pre_post_table['Assignment'] = np.where(pre_post_table['Pre_Probability']>0.5,'Pre','Post')
    

    return pre_post_table[['Sample_ID','Chromosome','Major_CN','Assignment']].copy()

def get_assignment_store(pre_post_table):
    aggremeent_matrix = np.zeros((2,2))
    assignment_store = {('Pre','Pre'):0,('Pre','Post'):0,('Post','Pre'):0,('Post','Post'):0}
    sample_n_chromosomes = pre_post_table.groupby('Sample_ID')['Chromosome'].nunique().to_dict()
    for (sample_id,chromosome),sample_chromosome_table in pre_post_table.groupby(['Sample_ID','Chromosome']):
        if len(sample_chromosome_table['Major_CN'].unique())<2:
            continue
        major_3_statuses = sample_chromosome_table[sample_chromosome_table['Major_CN']==3]['Assignment'].tolist()
        major_4_statuses = sample_chromosome_table[sample_chromosome_table['Major_CN']==4]['Assignment'].tolist()
        n_pairs = len(major_3_statuses)*len(major_4_statuses)
        for pair in itertools.product(major_3_statuses,major_4_statuses):
            assignment_store[pair] += 1
    return assignment_store
def get_permuted_table(pre_post_table):
    permuted_table_store = []
    for (sample_id,major_cn),sample_table in pre_post_table.groupby(['Sample_ID','Major_CN']):
        permuted_table = sample_table.copy()
        
        permuted_table['Assignment'] = RNG.permutation(permuted_table['Assignment'].tolist())
        
        permuted_table_store.append(permuted_table)
    permuted_table = pd.concat(permuted_table_store)
    return permuted_table

def get_permuted_assignment_store(pre_post_table,n_samples,permute_over_sample=True):
    permuted_assignment_store =  {('Pre','Pre'):[],('Pre','Post'):[],('Post','Pre'):[],('Post','Post'):[]}
    for i in range(n_samples):
        if permute_over_sample:
            permuted_table = get_permuted_table(pre_post_table)
        else:
            permuted_table = pre_post_table.copy()
            permuted_major_3 = permuted_table[permuted_table['Major_CN']==3].copy()
            permuted_major_4 = permuted_table[permuted_table['Major_CN']==4].copy()

            permuted_major_3['Assignment'] = RNG.permutation(permuted_major_3['Assignment'].tolist())
            permuted_major_4['Assignment'] = RNG.permutation(permuted_major_4['Assignment'].tolist())
            permuted_table = pd.concat([permuted_major_3,permuted_major_4])

        permuted_assignments = get_assignment_store(permuted_table)
        for key in permuted_assignment_store.keys():
            permuted_assignment_store[key].append(permuted_assignments[key])
    return permuted_assignment_store
def plot_bars(true_assigments,permuted_assignment_store,out_path):


    key_order = [('Pre','Pre'),('Pre','Post'),('Post','Pre'),('Post','Post')]
    key_label = ['Major 3 Pre WGD\nMajor 4 Pre WGD','Major 3 Pre WGD\nMajor 4 Post WGD','Major 3 Post WGD\nMajor 4 Pre WGD','Major 3 Post WGD\nMajor 4 Post WGD']

    true_data = [true_assigments[key] for key in key_order]
    permuted_data = [np.median(permuted_assignment_store[key]) for key in key_order]
    permuted_lower = [np.quantile(permuted_assignment_store[key],0.025) for key in key_order]
    permuted_upper = [np.quantile(permuted_assignment_store[key],0.975) for key in key_order]

    fig,ax = plt.subplots(figsize=(7,14/3))
    ax.bar(np.arange(len(key_order)),true_data,color='#EC4E20',label='True',width=0.25,align='edge')
    ax.bar(np.arange(len(key_order))+0.25,permuted_data,color='#5F5F58',label='Permuted',width=0.25,align='edge',yerr=[np.array(permuted_data)-np.array(permuted_lower),np.array(permuted_upper)-np.array(permuted_data)],capsize=5)
    ax.set_xticks(np.arange(len(key_order))+0.25/2)
    ax.set_xticklabels(key_label)
    ax.set_ylabel('Chromosome Timing Pairs')
    #remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
    fig.tight_layout()
    plt.savefig(out_path)

if __name__ == '__main__':
    
    n_permutations = 250
    for apply_penalty in [True,False]:
        
        pre_post_table = load_pre_post_table(apply_penalty)
        true_assigments = get_assignment_store(pre_post_table)
        contigency_table= np.array([[np.round(true_assigments[('Pre','Pre')]).astype(int),np.round(true_assigments[('Pre','Post')]).astype(int)],[np.round(true_assigments[('Post','Pre')]).astype(int),np.round(true_assigments[('Post','Post')]).astype(int)]])
    
    
        sample_dict = DataTools.get_sample_dict(pre_post_table)

        permuted_assignment_store_no_sample = get_permuted_assignment_store(pre_post_table,n_permutations,permute_over_sample=False)
        plot_bars(true_assigments,permuted_assignment_store_no_sample,f'../plots/major_3_4_agreement/permuted_no_sample_apply_penalty_{apply_penalty}.pdf')

        permuted_assignment_store_sample = get_permuted_assignment_store(pre_post_table,n_permutations,permute_over_sample=True)
        plot_bars(true_assigments,permuted_assignment_store_sample,f'../plots/major_3_4_agreement/permuted_sample_apply_penalty_{apply_penalty}.pdf')







