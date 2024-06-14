import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import DataTools
import multiprocessing as mp
import sys
RNG = np.random.default_rng(int(sys.argv[1])*42)
def get_cancer_type_counts(posterior_table,filter=True):

    cancer_type_counts = posterior_table[['Cancer_Type','WGD_Status','Sample_ID']].drop_duplicates().groupby(['Cancer_Type','WGD_Status']).count().reset_index()
    #pivot wider
    cancer_type_counts = cancer_type_counts.pivot(index='Cancer_Type', columns='WGD_Status', values='Sample_ID').reset_index(drop=False)
    cancer_type_counts = cancer_type_counts.rename(columns={True:'WGD',False:'No_WGD'}).fillna(0)
    #convert to integer
    cancer_type_counts.loc[:,'WGD'] = cancer_type_counts.loc[:,'WGD'].astype(int)
    cancer_type_counts.loc[:,'No_WGD'] = cancer_type_counts.loc[:,'No_WGD'].astype(int)
    #filter for cancer types with at least 10 samples and not XXX
    cancer_type_counts = cancer_type_counts[(cancer_type_counts['WGD'] >= 10) & (cancer_type_counts['No_WGD'] >= 10)]
    cancer_type_counts = cancer_type_counts[cancer_type_counts['Cancer_Type'] != 'XXX']
    return cancer_type_counts
   
def get_pre_post_sum_table(args):
    posterior_table,close_threshold = args

    sample_table_wgd_summary = posterior_table['WGD_Timing'].median()
    sample_table_wgd_ci_low = posterior_table['WGD_Timing'].quantile(0.11/2)
    sample_table_wgd_ci_high = posterior_table['WGD_Timing'].quantile(1-0.11/2)
    pre_probabilities = posterior_table.groupby(['Sample_ID','Segment_ID','Posterior_Sample_Index']).apply(lambda x: (x['Gain_Timing']<x['WGD_Timing']-close_threshold).any()).reset_index()
    pre_probabilities = pre_probabilities.groupby(['Sample_ID','Segment_ID'])[0].mean().reset_index().rename(columns={0:'Pre_Probability'})
    post_probabilities = posterior_table.groupby(['Sample_ID','Segment_ID','Posterior_Sample_Index']).apply(lambda x: (x['Gain_Timing']>x['WGD_Timing']+close_threshold).any()).reset_index()
    post_probabilities = post_probabilities.groupby(['Sample_ID','Segment_ID'])[0].mean().reset_index().rename(columns={0:'Post_Probability'})

    close_probabilities = posterior_table.groupby(['Sample_ID','Segment_ID','Posterior_Sample_Index']).apply(lambda x: ((x['Gain_Timing']>x['WGD_Timing']-close_threshold) & (x['Gain_Timing']<x['WGD_Timing']+close_threshold)).any()).reset_index()
    close_probabilities = close_probabilities.groupby(['Sample_ID','Segment_ID'])[0].mean().reset_index().rename(columns={0:'Close_Probability'})
    #merge pre and post and close probabilities
    pre_post_sums = pre_probabilities.merge(post_probabilities, how='outer', on=['Sample_ID','Segment_ID'])
    pre_post_sums = pre_post_sums.merge(close_probabilities, how='outer', on=['Sample_ID','Segment_ID'])
    pre_post_sums = pre_post_sums.fillna(0)
    posterior_table_info = posterior_table[['Sample_ID','Segment_ID','Major_CN','Minor_CN','N_Mutations']].drop_duplicates()
    pre_post_sums = pre_post_sums.merge(posterior_table_info, how='inner', on=['Sample_ID','Segment_ID'])
    pre_post_sums.loc[:,'WGD_Timing_Median'] = sample_table_wgd_summary
    pre_post_sums.loc[:,'WGD_Timing_CI_Low'] = sample_table_wgd_ci_low
    pre_post_sums.loc[:,'WGD_Timing_CI_High'] = sample_table_wgd_ci_high
    return pre_post_sums

def get_wgd_dist_store(posterior_table_wgd):
    wgd_dist_store ={cancer_type:[] for cancer_type in posterior_table_wgd['Cancer_Type'].unique()}
    for sample_id,sample_table in posterior_table_wgd.groupby('Sample_ID'):
        cancer_type = sample_table['Cancer_Type'].iloc[0]
        wgd_dist_store[cancer_type].append(sample_table['WGD_Timing'].to_numpy())
    return wgd_dist_store

def get_permuted_posterior_table_no_wgd(posterior_table_no_wgd,wgd_dist_store,sample_index):
    sample_tables = []
    processed_count = 0
    for sample_id,sample_table in posterior_table_no_wgd.groupby('Sample_ID'):
        sample_table = sample_table.copy()
        
        cancer_type = sample_table['Cancer_Type'].iloc[0]
        wgd_dist = wgd_dist_store[cancer_type][RNG.integers(len(wgd_dist_store[cancer_type]))]
        #to do this should probably be changed so that each segment has a different wgd timing dstribution
        wgd_table = pd.DataFrame({'WGD_Timing':RNG.choice(wgd_dist,size=100,replace=True),'Posterior_Sample_Index':np.arange(100)})
        sample_table = sample_table.drop(columns=['WGD_Timing']).merge(wgd_table, how='inner', on='Posterior_Sample_Index')
        sample_table.loc[:,'Sample_ID'] = f'{sample_id}_permuted_{sample_index}'
        sample_tables.append(sample_table)
        processed_count += 1
    return pd.concat(sample_tables)
def get_posterior_permutations(posterior_table_no_wgd,posterior_table_wgd,n_samples=30):
    posterior_permutations = []
    wgd_dist_store = get_wgd_dist_store(posterior_table_wgd)
    
    for i in range(n_samples):
        permuted_posterior_table_no_wgd =get_permuted_posterior_table_no_wgd(posterior_table_no_wgd,wgd_dist_store,i)
        posterior_permutations.append(permuted_posterior_table_no_wgd)
    posterior_permutations = pd.concat(posterior_permutations)
    return  posterior_permutations


def get_run_types():
    run_types = []
    for apply_penalty in [True,False]:
        for permute in [True,False]:
            run_types.append((apply_penalty,permute))
    return run_types

if __name__ == '__main__':
    run_types = get_run_types()
    run_index = int(sys.argv[1])

    apply_penalty,permute = run_types[run_index]

    output_path = f'../output/pre_post_close/{apply_penalty}_with_permutations_{permute}_pre_post.tsv'
    temp_output_path = output_path.replace('.tsv','.temp.tsv')

    posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty)
    posterior_table_wgd = posterior_table[posterior_table['WGD_Status']]
    
    posterior_table_no_wgd = posterior_table[~posterior_table['WGD_Status']]

    
    if permute:
        cancer_type_counts = get_cancer_type_counts(posterior_table)
        posterior_table_wgd = posterior_table_wgd[posterior_table_wgd['Cancer_Type'].isin(cancer_type_counts['Cancer_Type'])]
        posterior_table_no_wgd = posterior_table_no_wgd[posterior_table_no_wgd['Cancer_Type'].isin(cancer_type_counts['Cancer_Type'])]
        posterior_table_no_wgd_permuted = get_posterior_permutations(posterior_table_no_wgd,posterior_table_wgd,n_samples=5)
        posterior_table = pd.concat([posterior_table_wgd,posterior_table_no_wgd_permuted])
    else:
        posterior_table = posterior_table_wgd
    

    n_samples = len(posterior_table['Sample_ID'].unique())
    processsed_count = 0
    sample_dict = DataTools.get_sample_dict(posterior_table)
    del posterior_table
    store = []
    sample_ids = list(sample_dict.keys())
    args = [(sample_dict[sample_id],0) for sample_id in sample_ids]
    print('Starting multiprocessing')
    with mp.Pool(mp.cpu_count()) as pool:
        for result in pool.imap_unordered(get_pre_post_sum_table,args):
            store.append(result)
    pre_post_table = pd.concat(store)
    pre_post_table.to_csv(output_path,sep='\t',index=False)


