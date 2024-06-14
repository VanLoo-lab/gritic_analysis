import pandas as pd
import numpy as np
import DataTools
import multiprocessing as mp
import itertools
from event_proportions_relative_to_wgd_probabilistic import load_gain_timing_table,load_loss_timing_table


def get_bin_counts(sample_table,genome_bins,event_type,sample_id):
    
    pre_event_type,post_event_type,non_event_type = f'Pre_WGD_{event_type}',f'Post_WGD_{event_type}',f'Non_WGD_{event_type}'
    
    bin_counts = sample_table.merge(genome_bins,on=['Chromosome'],how='inner')
    if event_type == 'Loss':
        autosomes = list(map(str,range(1,23)))
        bin_counts = bin_counts[bin_counts['Chromosome'].isin(autosomes)]
    #filter so that Bin intersects with Segment
    bin_counts = bin_counts[(bin_counts['Bin_Start']<=bin_counts['Segment_End'])&(bin_counts['Bin_End']>=bin_counts['Segment_Start'])]
  
    #get overlap width between segment and bin for each row
    bin_counts['Overlap_Width'] = bin_counts[['Segment_Start','Segment_End','Bin_Start','Bin_End']].apply(lambda x: min(x['Segment_End'],x['Bin_End'])-max(x['Segment_Start'],x['Bin_Start']),axis=1)
    #get overlap fraction between segment and bin for each row
    bin_counts['Overlap_Fraction'] = bin_counts['Overlap_Width']/(bin_counts['Bin_End']-bin_counts['Bin_Start'])

    #filter to bins with Pre_WGD_Gain or Post_WGD_Gain > 0.5 and round to 1
    bin_counts[pre_event_type] = bin_counts[pre_event_type].apply(lambda x: 1 if x>=0.5 else 0)
    bin_counts[post_event_type] = bin_counts[post_event_type].apply(lambda x: 1 if x>=0.5 else 0)
    bin_counts[non_event_type] = bin_counts[non_event_type].apply(lambda x: 1 if x>=0.5 else 0)
   

    pre_wgd_bins = bin_counts.groupby(['Bin_ID']).apply(lambda x: np.sum(x[pre_event_type]*x['Overlap_Fraction']))
    post_wgd_bins = bin_counts.groupby(['Bin_ID']).apply(lambda x: np.sum(x[post_event_type]*x['Overlap_Fraction']))
    non_wgd_bins = bin_counts.groupby(['Bin_ID']).apply(lambda x: np.sum(x[non_event_type]*x['Overlap_Fraction']))

    bin_tables = [pre_wgd_bins,post_wgd_bins,non_wgd_bins]

    bin_counts = pd.concat(bin_tables,axis=1).reset_index().rename(columns={0:pre_event_type,1:post_event_type,2:non_event_type})

    #filter to bins with Pre_WGD_Gain or Post_WGD_Gain > 0.5 and round to 1
    bin_counts[pre_event_type] = bin_counts[pre_event_type].apply(lambda x: 1 if x>=0.5 else 0)
    bin_counts[post_event_type] = bin_counts[post_event_type].apply(lambda x: 1 if x>=0.5 else 0)
    bin_counts[non_event_type] = bin_counts[non_event_type].apply(lambda x: 1 if x>=0.5 else 0)
    bin_counts['Sample_ID'] = sample_id
    return bin_counts

def get_sample_metadata():
    good_samples =DataTools.get_good_samples()
    wgd_status_calls = DataTools.get_wgd_calling_info()
    wgd_status_calls = wgd_status_calls[wgd_status_calls['Sample_ID'].isin(good_samples)]
    wgd_status_calls = wgd_status_calls[~((wgd_status_calls['Major_CN_Mode']==2) & (~wgd_status_calls['WGD_Status']))]
    wgd_status_calls = wgd_status_calls[['Sample_ID','WGD_Status']]
    classifications = DataTools.load_classifications()[['Sample_ID','Cancer_Type']]
    sample_metadata = wgd_status_calls.merge(classifications,how='inner')
    return sample_metadata

def add_missing_events(bin_counts,sample_metadata,genome_arm_bins):
    all_samples = list(sample_metadata['Sample_ID'].unique())
    all_bins = list(genome_arm_bins['Bin_ID'].unique())
    pairs= itertools.product(all_samples,all_bins)
    pairs = list(pairs)
    
    full_sample_table = pd.DataFrame({'Sample_ID':[pair[0] for pair in pairs],'Bin_ID':[pair[1] for pair in pairs]})
    bin_counts = full_sample_table.merge(bin_counts,on=['Sample_ID','Bin_ID'],how='left').fillna(0)
    return bin_counts

if __name__ == '__main__':
    out_dir = '../output/arm_pre_post_clean'
    sample_metadata = get_sample_metadata()
    
    genome_bins = DataTools.load_genome_arm_bins()
    loss_timing_data = load_loss_timing_table()
     
    loss_timing_dict = DataTools.get_sample_dict(loss_timing_data)
    gain_timing_store = {True:load_gain_timing_table(True,'All','absolute'),False:load_gain_timing_table(False,'All','absolute')}
    gain_timing_dict = {True:DataTools.get_sample_dict(gain_timing_store[True]),False:DataTools.get_sample_dict(gain_timing_store[False])}

    sample_id_test = list(loss_timing_dict.keys())[0]
    loss_bin_timing_test = get_bin_counts(loss_timing_dict[sample_id_test],genome_bins,'Loss',sample_id_test)
    gain_bin_timing_test = get_bin_counts(gain_timing_dict[True][sample_id_test],genome_bins,'Gain',sample_id_test)
    gain_bin_timing_test_false = get_bin_counts(gain_timing_dict[False][sample_id_test],genome_bins,'Gain',sample_id_test)
    
    
    #multiprocessing version
    pool = mp.Pool(processes=mp.cpu_count())
    loss_timing_bin_counts = pool.starmap(get_bin_counts,[(loss_timing_dict[sample_id],genome_bins,'Loss',sample_id) for sample_id in loss_timing_dict.keys()])
    pool.close()
    pool.join()
    loss_timing_bin_counts = pd.concat(loss_timing_bin_counts,axis=0)
  
    loss_timing_bin_counts = add_missing_events(loss_timing_bin_counts,sample_metadata,genome_bins)
    loss_timing_bin_counts = loss_timing_bin_counts[~loss_timing_bin_counts['Bin_ID'].str.contains('X')]
    loss_timing_bin_counts = loss_timing_bin_counts.merge(sample_metadata,how='inner')
    loss_timing_bin_counts.to_csv(f'{out_dir}/loss_timing_bin_counts.tsv',sep='\t',index=False)
    print('Loss timing bin counts complete')
    
    
    pool = mp.Pool(processes=mp.cpu_count())
    gain_timing_true_bin_counts = pool.starmap(get_bin_counts,[(gain_timing_dict[True][sample_id],genome_bins,'Gain',sample_id) for sample_id in gain_timing_dict[True].keys()])
    pool.close()
    pool.join()
    gain_timing_true_bin_counts = pd.concat(gain_timing_true_bin_counts,axis=0)
    gain_timing_true_bin_counts = add_missing_events(gain_timing_true_bin_counts,sample_metadata,genome_bins)
    gain_timing_true_bin_counts = gain_timing_true_bin_counts.merge(sample_metadata,how='inner')
    gain_timing_true_bin_counts.to_csv(f'{out_dir}/gain_timing_True_bin_counts.tsv',sep='\t',index=False)

    pool = mp.Pool(processes=mp.cpu_count())
    gain_timing_false_bin_counts = pool.starmap(get_bin_counts,[(gain_timing_dict[False][sample_id],genome_bins,'Gain',sample_id) for sample_id in gain_timing_dict[False].keys()])
    pool.close()
    pool.join()
    gain_timing_false_bin_counts = pd.concat(gain_timing_false_bin_counts,axis=0)
    gain_timing_false_bin_counts = add_missing_events(gain_timing_false_bin_counts,sample_metadata,genome_bins)
    gain_timing_false_bin_counts = gain_timing_false_bin_counts.merge(sample_metadata,how='inner')
    gain_timing_false_bin_counts.to_csv(f'{out_dir}/gain_timing_False_bin_counts.tsv',sep='\t',index=False)

    
    