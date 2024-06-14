import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import DataTools
import multiprocessing as mp
import sys



def get_pre_post_sum_table_wgd_absolute(sample_posterior_table):
    pre_probabilities = sample_posterior_table.groupby(['Sample_ID','Segment_ID','Posterior_Sample_Index']).apply(lambda x: (x['Gain_Timing']<x['WGD_Timing']).any()).reset_index()
    pre_probabilities = pre_probabilities.groupby(['Sample_ID','Segment_ID'])[0].mean().reset_index().rename(columns={0:'Pre_Probability'})
    
    
    post_probabilities = sample_posterior_table.groupby(['Sample_ID','Segment_ID','Posterior_Sample_Index']).apply(lambda x: (x['Gain_Timing']>x['WGD_Timing']).any()).reset_index()
    post_probabilities = post_probabilities.groupby(['Sample_ID','Segment_ID'])[0].mean().reset_index().rename(columns={0:'Post_Probability'})
    #merge pre and post and close probabilities
    pre_post_sums = pre_probabilities.merge(post_probabilities, how='outer', on=['Sample_ID','Segment_ID'])

    pre_post_sums = pre_post_sums.fillna(0)
    sample_posterior_table_info = sample_posterior_table[['Sample_ID','Segment_ID','Major_CN','Minor_CN','N_Mutations']].drop_duplicates()
    pre_post_sums = pre_post_sums.merge(sample_posterior_table_info, how='inner', on=['Sample_ID','Segment_ID'])
    pre_post_sums.loc[:,'Non_Probability'] = 0
    return pre_post_sums

def get_pre_post_sum_table_wgd_summed(sample_posterior_table):
    pre_probabilities = sample_posterior_table.groupby(['Sample_ID','Segment_ID','Posterior_Sample_Index']).apply(lambda x: (x['Gain_Timing']<x['WGD_Timing']).sum()).reset_index()
    pre_probabilities = pre_probabilities.groupby(['Sample_ID','Segment_ID'])[0].mean().reset_index().rename(columns={0:'Pre_Probability'})
    
    
    post_probabilities = sample_posterior_table.groupby(['Sample_ID','Segment_ID','Posterior_Sample_Index']).apply(lambda x: (x['Gain_Timing']>x['WGD_Timing']).sum()).reset_index()
    post_probabilities = post_probabilities.groupby(['Sample_ID','Segment_ID'])[0].mean().reset_index().rename(columns={0:'Post_Probability'})
  
    #merge pre and post and close probabilities
    pre_post_sums = pre_probabilities.merge(post_probabilities, how='outer', on=['Sample_ID','Segment_ID'])
    
    pre_post_sums = pre_post_sums.fillna(0)
    sample_posterior_table_info = sample_posterior_table[['Sample_ID','Segment_ID','Major_CN','Minor_CN','N_Mutations']].drop_duplicates()
    pre_post_sums = pre_post_sums.merge(sample_posterior_table_info, how='inner', on=['Sample_ID','Segment_ID'])
    pre_post_sums.loc[:,'Non_Probability'] = 0
    
    return pre_post_sums

def get_pre_post_sum_table_no_wgd(sample_posterior_table,mode):
    if mode == 'absolute':
        segment_table = sample_posterior_table[['Sample_ID','Segment_ID','Major_CN','Minor_CN','N_Mutations']].drop_duplicates().reset_index(drop=True)
        segment_table['Pre_Probability'] = 0
        segment_table['Post_Probability'] = 0
        segment_table['Non_Probability'] = 1
    else:
        segment_table = sample_posterior_table.groupby(['Sample_ID','Segment_ID','Major_CN','Minor_CN','N_Mutations']).size().reset_index().rename(columns={0:'Non_Probability'})
        segment_table['Non_Probability'] = segment_table['Non_Probability']/100
        segment_table['Pre_Probability'] = 0
        segment_table['Post_Probability'] = 0
    return segment_table
def get_pre_post_sum_table(sample_posterior_table,mode):
    wgd_status = sample_posterior_table['WGD_Status'].iloc[0]
    if wgd_status:
        if mode == 'absolute':
            return get_pre_post_sum_table_wgd_absolute(sample_posterior_table)
        else:
            return get_pre_post_sum_table_wgd_summed(sample_posterior_table)
    else:
        return get_pre_post_sum_table_no_wgd(sample_posterior_table,mode)

def get_cohort_table(sample_dict,mode):
    store = []
    sample_ids = list(sample_dict.keys())
    args = [(sample_dict[sample_id],mode) for sample_id in sample_ids]
    with mp.Pool(mp.cpu_count()) as pool:
        for result in pool.starmap(get_pre_post_sum_table,args):
            store.append(result)
    pre_post_table = pd.concat(store)
    pre_post_table['Mode'] = mode
    return pre_post_table

  


if __name__ == '__main__':

    for apply_penalty in [True,False]:
        for major_cn in ['All',3,4,5,6,7]:
            
            posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty)
            if major_cn != 'All':
                posterior_table = posterior_table[posterior_table['Major_CN']==major_cn]
            
            sample_dict = DataTools.get_sample_dict(posterior_table)
            for mode in ['absolute','sum']:
                output_path = f'../output/pre_post_non/pre_post_non_gains_apply_penalty_{apply_penalty}_major_cn_{major_cn}_mode_{mode}.tsv'
                pre_post_table = get_cohort_table(sample_dict,mode)
                pre_post_table.to_csv(output_path,index=False,sep='\t')
                print(pre_post_table)
                