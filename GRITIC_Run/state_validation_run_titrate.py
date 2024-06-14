import numpy as np
import pandas as pd
import string
from gritic import gritictimer,sampletools,dataloader
import os
import sys
sys.path.insert(0,'../code')
import DataTools

import PCAWGDataLoader
import HartwigDataLoader
import SampleSimulator

import argparse

import ast

from collections import Counter

from scipy.stats import poisson
import matplotlib.pyplot as plt
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
RNG = np.random.default_rng()
class SampleLoader:
    def __init__(self) -> None:
        self.good_samples = DataTools.get_good_samples()
        self.pcawg_ids = self.load_pcawg_ids()
        
        self.hartwig_ids,self.hartwig_set_names = self.load_hartwig_ids()
        
        self.pcawg_weighting = len(self.pcawg_ids)/(len(self.pcawg_ids)+len(self.hartwig_ids))
    def load_pcawg_ids(self):
        pcawg_metadata = pd.read_csv('/PCAWG_DIR/consensus.20170217.purity.ploidy.txt.gz',sep="\t")
        wgd_calling_info = pd.read_csv('../in_data/wgd_calling_info.tsv',sep="\t")
        pcawg_metadata = pcawg_metadata[pcawg_metadata['samplename'].isin(wgd_calling_info['Sample_ID'])]
        pcawg_metadata = pcawg_metadata[pcawg_metadata['samplename'].isin(self.good_samples)]
        return list(pcawg_metadata['samplename'])
    
    def load_hartwig_ids(self):
        hartwig_metadata_path = '/HARTWIG_DIR/metadata.tsv'
        hartwig_metadata_table = pd.read_csv(hartwig_metadata_path,encoding='iso-8859-1',sep="\t")
        wgd_calling_info = pd.read_csv('../in_data/wgd_calling_info.tsv',sep="\t")
        hartwig_metadata_table = hartwig_metadata_table[hartwig_metadata_table['sampleId'].isin(wgd_calling_info['Sample_ID'])]
        hartwig_metadata_table = hartwig_metadata_table[hartwig_metadata_table['sampleId'].isin(self.good_samples)]
        return list(hartwig_metadata_table['sampleId']),list(hartwig_metadata_table['setName'])
    
    def get_random_sample(self):
        if RNG.random()<self.pcawg_weighting:
            sample_id,dataset = RNG.choice(self.pcawg_ids),'PCAWG'
            sample_data = PCAWGDataLoader.load_sample_data(sample_id)
        else:
            sample_id_index = RNG.integers(0,len(self.hartwig_ids))
            sample_id,set_name,dataset = self.hartwig_ids[sample_id_index],self.hartwig_set_names[sample_id_index],'Hartwig'
            sample_data = HartwigDataLoader.load_sample_data(sample_id,set_name,'PURPLE')
        return sample_id,dataset,sample_data



class Random_Gain_Timing:
    def __init__(self) -> None:
        pass
    def get_random_gain_timing(self,n_timings):
        return RNG.uniform(0.05,0.95,size=n_timings)



def run_sample(sample,output_dir,wgd_override=None):

    gritictimer.process_sample(sample,output_dir,plot_trees=False,wgd_override=wgd_override)


def generate_complex_segment_frame(nrpcc,max_routes_per_state=500):

    mutations_range = [5,10,20,30,50,80,100,250,1000]

    estimated_missing_rate = 1-poisson.cdf(2,nrpcc)


    
    wgd_routes = pd.read_csv('in_data/cn_states_n_routes.tsv',sep="\t").rename(columns={'Route_ID':'Route'})
    wgd_routes = wgd_routes[wgd_routes['WGD_Status']].copy()
    
    wgd_routes = wgd_routes[(wgd_routes['Major_CN']>=3)].copy()
    n_routes = wgd_routes.groupby(['Major_CN','Minor_CN']).size().reset_index().rename(columns={0:'N_Routes'})
    cn_states = n_routes[n_routes['N_Routes']<=max_routes_per_state][['Major_CN','Minor_CN']].drop_duplicates()
    
    chosen_states = cn_states.sample(n=1,random_state=RNG).iloc[0]

    #sample one row
    wgd_routes =wgd_routes[(wgd_routes['Major_CN']==chosen_states['Major_CN']) & (wgd_routes['Minor_CN']==chosen_states['Minor_CN'])].copy()

    
    wgd_routes =wgd_routes.groupby(['Major_CN','Minor_CN']).apply(lambda x: x.sample(n=1,random_state=RNG)).reset_index(drop=True)
    wgd_routes =wgd_routes[['Major_CN','Minor_CN','Route']].copy()
    
    wgd_routes_store = []
    for n_mutations in mutations_range:
        n_mutations_target = n_mutations/estimated_missing_rate

        n_mutations_sim =int(n_mutations_target*3)
        
        wgd_routes_n_mutations = wgd_routes.copy()
        wgd_routes_n_mutations['N_Mutations'] = n_mutations_sim
        wgd_routes_n_mutations['N_Mutations_Target'] = n_mutations
        wgd_routes_store.append(wgd_routes_n_mutations)
    wgd_routes = pd.concat(wgd_routes_store,ignore_index=True)
  
    
    wgd_routes['Segment_Width'] = wgd_routes['N_Mutations']+5
    #adding 50Mb here stops it overlapping with any real segments
    wgd_routes['Segment_Start'] = np.cumsum(wgd_routes['Segment_Width'])-wgd_routes['Segment_Width']+int(50e6)
    wgd_routes['Segment_End'] = wgd_routes['Segment_Start']+wgd_routes['Segment_Width']-1
    wgd_routes['Chromosome'] = '21'
    wgd_routes['Segment_ID'] = wgd_routes['Chromosome'] + '-' + wgd_routes['Segment_Start'].astype(str) + '-' + wgd_routes['Segment_End'].astype(str)
    
    return wgd_routes


def make_sim_template(mutation_table,purity,nrpcc):

    cn_table_summary = mutation_table.groupby(['Segment_ID','Major_CN','Minor_CN','Chromosome','Segment_Start','Segment_End']).agg({'Position':'count'}).reset_index()
    cn_table_summary = cn_table_summary.rename(columns={'Position':'N_Mutations'})
    
    
    wgd_segments = cn_table_summary[(cn_table_summary['N_Mutations']>=10) & (cn_table_summary['Major_CN']==2) ].copy()
    if len(wgd_segments.index)==0:
        return pd.DataFrame()
    complex_segments = generate_complex_segment_frame(nrpcc)

    
    if len(wgd_segments.index)>0:
        wgd_segments['Route'] = np.nan
    
    sim_template = pd.concat([wgd_segments,complex_segments],ignore_index=True)
    sim_template['Total_CN'] = sim_template['Major_CN']+sim_template['Minor_CN']
    sim_template['Adjusted_Coverage'] = nrpcc*(sim_template['Total_CN']*purity + 2*(1-purity))/purity
    sim_template['Segment_Width'] = sim_template['Segment_End']-sim_template['Segment_Start']
    sim_template['Gain_Timing_Target'] = pd.Series([np.array([])]*len(sim_template.index),index=sim_template.index)
    
    return sim_template




def get_sample_data(max_sample_attempts=1000):
    for _ in range(max_sample_attempts):
        try:
            dataset_sample_id,dataset,sample_data = sample_loader.get_random_sample()
        except ValueError:
            print('couldnt get sample data for',sample_id)
            continue
        
        major_cn_mode = sampletools.get_major_cn_mode_from_cn_table(sample_data['CN_Table'])
        if major_cn_mode != 2:
            continue
        return dataset_sample_id,dataset,sample_data
    raise ValueError('Couldnt get sample data,max attempts reached')

def get_n_gains(segment_row):
    #all samples are true wgd
    baseline_cn = 1
    n_gains = max(segment_row['Major_CN'] - baseline_cn,0) + max(segment_row['Minor_CN'] - baseline_cn,0)
    return n_gains

def get_sample_run_posterior_table(posterior_table_path,purity,wgd_status,phasing,apply_reads_correction):
    
    
    posterior_table = pd.read_csv(posterior_table_path, sep='\t',dtype={'Chromosome':str})
    posterior_table['Purity'] = purity
    posterior_table['WGD_Status'] = wgd_status
    posterior_table['Phasing'] = phasing
    posterior_table['Apply_Reads_Correction'] = apply_reads_correction
    posterior_table = posterior_table[posterior_table['Major_CN']>2]
    posterior_subtable = posterior_table[['Segment_ID','Major_CN','Minor_CN','Posterior_Sample_Index','Gain_Timing','WGD_Timing','N_Mutations']].copy()
    posterior_subtable = posterior_subtable.rename(columns={'N_Mutations':'True_N_Mutations'})
    new_subtable_store = []
    for (segment_id,posterior_index),segment_subtable in posterior_subtable.groupby(['Segment_ID','Posterior_Sample_Index']):
        gain_timing = list(segment_subtable['Gain_Timing'].values)
        wgd_timing = segment_subtable['WGD_Timing'].values[0]
        n_gains = get_n_gains(segment_subtable.iloc[0])
        gain_timing.extend([wgd_timing]*max(n_gains - len(gain_timing),0))
        gain_timing = sorted(gain_timing)
        new_subtable_data = {'Segment_ID':[segment_id]*n_gains,'Gain_Index':list(range(n_gains)),'Gain_Timing':gain_timing,'Posterior_Sample_Index':[posterior_index]*n_gains,'WGD_Timing':[wgd_timing]*n_gains}
        
        
        new_subtable = pd.DataFrame(new_subtable_data)
        new_subtable_store.append(new_subtable)
    new_subtable_store = pd.concat(new_subtable_store)

    posterior_table = posterior_table.drop(columns=['Node','Gain_Timing','WGD_Timing','Gain_Index','Node_Phasing']).drop_duplicates()

    posterior_table = posterior_table.merge(new_subtable_store, on=['Segment_ID','Posterior_Sample_Index'],how='inner')
   
    return posterior_table

def get_processed_posterior_table(posterior_table,posterior_table_path):
  
    posterior_table = posterior_table.drop_duplicates(subset=['Segment_ID', 'Posterior_Sample_Index','Gain_Index'],keep='first').copy()
    posterior_table = posterior_table[posterior_table['True_Independent_Gain']]
    
    posterior_table['Error'] =  posterior_table['Gain_Timing'] - posterior_table['True_Gain_Timing']
    posterior_table['CN_State'] = posterior_table['Major_CN'].astype(str) + '+' + posterior_table['Minor_CN'].astype(str)
    posterior_table = posterior_table.rename(columns={'N_Mutations':'N_Mutations_Observed'})
    
    group_cols = ['Segment_ID','Chromosome', 'Segment_Start', 'Segment_End', 'Major_CN', 'Minor_CN', 'WGD_Status', 'Purity', 'Phasing', 'Sample_ID', 'N_Mutations_Observed', 'N_Mutations_Full','N_Mutations_Target','Coverage','True_Route','Apply_Reads_Correction','NRPCC']
    group_cols_with_posterior = group_cols + ['Posterior_Sample_Index']
    
    posterior_table_error= posterior_table.groupby(group_cols)['Error'].apply(lambda x: np.mean(x)).reset_index().rename(columns={'Error':'Mean_Error'})
    posterior_table_absolute_error = posterior_table.groupby(group_cols)['Error'].apply(lambda x: np.mean(np.abs(x))).reset_index().rename(columns={'Error':'Mean_Absolute_Error'})
    
    posterior_table_true_pre_wgd_gains = posterior_table.groupby(group_cols_with_posterior).apply(lambda x: np.sum((~np.isclose(x['True_Gain_Timing'], x['True_WGD_Timing'])) & (x['True_Gain_Timing'] < x['True_WGD_Timing']))).reset_index().rename(columns={0:'True_Pre_WGD_Gains'}).groupby(group_cols).agg({'True_Pre_WGD_Gains':np.mean}).reset_index()
    posterior_table_true_post_wgd_gains = posterior_table.groupby(group_cols_with_posterior).apply(lambda x: np.sum((~np.isclose(x['True_Gain_Timing'], x['True_WGD_Timing'])) & (x['True_Gain_Timing'] > x['True_WGD_Timing']))).reset_index().rename(columns={0:'True_Post_WGD_Gains'}).groupby(group_cols).agg({'True_Post_WGD_Gains':np.mean}).reset_index()
    posterior_table_measured_pre_wgd_gains = posterior_table.groupby(group_cols_with_posterior).apply(lambda x: np.sum((~np.isclose(x['Gain_Timing'], x['True_WGD_Timing'])) & (x['Gain_Timing'] < x['True_WGD_Timing']))).reset_index().rename(columns={0:'Measured_Pre_WGD_Gains'}).groupby(group_cols).agg({'Measured_Pre_WGD_Gains':np.mean}).reset_index()
    posterior_table_measured_post_wgd_gains = posterior_table.groupby(group_cols_with_posterior).apply(lambda x: np.sum((~np.isclose(x['Gain_Timing'], x['True_WGD_Timing'])) & (x['Gain_Timing'] > x['True_WGD_Timing']))).reset_index().rename(columns={0:'Measured_Post_WGD_Gains'}).groupby(group_cols).agg({'Measured_Post_WGD_Gains':np.mean}).reset_index()
    
    posterior_table_route_correct = posterior_table.drop_duplicates(subset=['Segment_ID', 'Posterior_Sample_Index'],keep='first').groupby(group_cols_with_posterior).apply(lambda x: x['True_Route'].str.slice(0,9) == x['Route']).reset_index().rename(columns={0:'Route_Correct'}).groupby(group_cols).agg({'Route_Correct':np.mean}).reset_index()
    posterior_table_route_counts= posterior_table.drop_duplicates(subset=['Segment_ID', 'Posterior_Sample_Index'],keep='first').groupby(group_cols)['Route'].apply(lambda x: str(dict(Counter(x)))).reset_index().rename(columns={'Route':'Route_Counts'})
    
    posterior_table_processed = posterior_table_error.merge(posterior_table_absolute_error, how='inner')
    posterior_table_processed = posterior_table_processed.merge(posterior_table_true_pre_wgd_gains, how='inner')
    posterior_table_processed = posterior_table_processed.merge(posterior_table_true_post_wgd_gains, how='inner')
    posterior_table_processed = posterior_table_processed.merge(posterior_table_measured_pre_wgd_gains, how='inner')
    posterior_table_processed = posterior_table_processed.merge(posterior_table_measured_post_wgd_gains, how='inner')
    posterior_table_processed = posterior_table_processed.merge(posterior_table_route_correct, how='inner')
    posterior_table_processed = posterior_table_processed.merge(posterior_table_route_counts, how='inner')

    subclone_table_path = '/'.join(posterior_table_path.split('/')[:-1]) + '/test_subclone_table.tsv'
    if os.path.exists(subclone_table_path):
        subclone_table = pd.read_csv(subclone_table_path, sep='\t')
        n_subclones = len(subclone_table)
        min_subclone_ccf = subclone_table['Subclone_CCF'].min()
    else:
        n_subclones = 0
        min_subclone_ccf = pd.NA

    posterior_table_processed['N_Subclones'] = n_subclones
    posterior_table_processed['Min_Subclone_CCF'] = min_subclone_ccf
    
    for col in ['True_Pre_WGD_Gains', 'True_Post_WGD_Gains']:
        posterior_table_processed[col] = posterior_table_processed[col].astype(int)
    posterior_table_processed['True_Route'] = posterior_table_processed['True_Route'].str.slice(0,9)
    posterior_table_processed = posterior_table_processed.rename(columns={'Error':'Mean_Error'})
    return posterior_table_processed

def get_processed_metadata(metadata):
    metadata = metadata.copy()
    
    metadata_store = {'Sample_ID':sample_id,'Segment_ID':[],'Gain_Index':[],'True_Gain_Timing':[],'N_Mutations':[],'N_Mutations_Target':[],'Adjusted_Coverage':[],'True_Route':[],'True_WGD_Timing':[],'True_Independent_Gain':[],'Major_CN':[],'Minor_CN':[]}
    for segment_index,segment_row in metadata.iterrows():
        gain_dict = ast.literal_eval(segment_row['Gain_Timing'])
        n_independent_gains = len(gain_dict)
        n_gains = get_n_gains(segment_row)

        
        gain_timing = list(gain_dict.values())
        gain_timing.extend([segment_row['WGD_Timing']]*max(n_gains - n_independent_gains,0))
        independent_gain_store = [1]*n_independent_gains + [0]*max(n_gains - n_independent_gains,0)
        independent_gain_store = np.array(independent_gain_store)
        independent_gain_store = independent_gain_store[np.argsort(gain_timing)].astype(bool)
        gain_timing = sorted(gain_timing)

        metadata_store['Gain_Index'].extend(list(range(n_gains)))
        metadata_store['True_Gain_Timing'].extend(gain_timing)
        metadata_store['True_Independent_Gain'].extend(independent_gain_store)
        for col in ['Segment_ID','N_Mutations','N_Mutations_Target','Adjusted_Coverage','True_Route','True_WGD_Timing','Major_CN','Minor_CN']:
            if col == 'True_Route':
                metadata_store[col].extend([segment_row['Route']]*n_gains)
            elif col == 'True_WGD_Timing':
                metadata_store[col].extend([segment_row['WGD_Timing']]*n_gains)
            else:
                metadata_store[col].extend([segment_row[col]]*n_gains)

    
    metadata_df = pd.DataFrame(metadata_store)
    metadata_df = metadata_df.rename(columns={'Adjusted_Coverage':'Coverage','N_Mutations':'N_Mutations_Full'})
    metadata_df['Segment_ID'] = metadata_df['Segment_ID'].str.replace('Segment_','')
    return metadata_df



def n_mutations_selection(base_mutation_table,segment_metadata):
    new_mutation_table_store = [base_mutation_table[base_mutation_table['Major_CN']==2].copy()]
    base_mutation_table = base_mutation_table.merge(segment_metadata[['Segment_ID','N_Mutations_Target']],on='Segment_ID',how='inner')
    
    for (segment_id,n_mutations_target),segment_table in base_mutation_table.groupby(['Segment_ID','N_Mutations_Target']):
        if segment_table['Major_CN'].iloc[0] == 2:
            continue
        n_mutations = len(segment_table)
        if n_mutations > n_mutations_target:
            segment_table = segment_table.sample(n=int(n_mutations_target),random_state=RNG,replace=False).copy()
        if len(segment_table)<n_mutations_target:
            continue
        new_mutation_table_store.append(segment_table)
        
    new_mutation_table = pd.concat(new_mutation_table_store)
    return new_mutation_table


if __name__ == '__main__':

    processed_count = 0
    my_parser = argparse.ArgumentParser(allow_abbrev=False)
    my_parser.add_argument('--base_output_dir', action='store', type=str, required=True)

    args = my_parser.parse_args()

    random_gain_generator =  Random_Gain_Timing()
    sample_loader = SampleLoader()    

    sample_id = SampleSimulator.generate_random_sample_id(length=12)

    while True:
        dataset_sample_id,dataset,sample_data = get_sample_data()
        if dataloader.get_major_cn_mode_from_cn_table(sample_data['CN_Table']) ==2:
            break
    
    wgd_status = True
    wgd_timing = RNG.uniform(0.01,0.99)
    purity = np.round(RNG.uniform(0.1,0.9),2)
 
    nrpcc = RNG.choice([1,3,5,10,15,20,25])

    base_sim_template = make_sim_template(sample_data['Mutation_Table'],purity,nrpcc)
    if len(base_sim_template.index)==0:
        raise ValueError('No segments in sim template')
    

    base_sim_template['Sample_ID']= sample_id
    base_sim_template['WGD_Timing'] = wgd_timing
    
    output_dir = f"{args.base_output_dir}/{sample_id}"
    os.makedirs(output_dir,exist_ok=True)
    
    base_sim_template.to_csv(f'{output_dir}/test_sim_template.tsv',sep="\t",index=False)
    if sample_data['Subclone_Table'] is not None:
        sample_data['Subclone_Table'].to_csv(f'{output_dir}/test_subclone_table.tsv',sep="\t",index=False)
    
    sim_template = base_sim_template.copy()
    sim_template['Purity'] = purity

  
   
    
    base_mutation_table,segment_metadata = SampleSimulator.simulate_sample(sim_template,sample_data['Subclone_Table'],wgd_status=wgd_status,sex='XX')
    base_mutation_table = base_mutation_table[base_mutation_table['Tumor_Alt_Count']>=3].copy()
    base_mutation_table = base_mutation_table[(base_mutation_table['Tumor_Ref_Count']+ base_mutation_table['Tumor_Alt_Count'])>=10].copy()

    base_mutation_table = n_mutations_selection(base_mutation_table,segment_metadata)
    

    
    base_mutation_table['Sample_ID'] = sample_id

    base_cn_table = base_mutation_table[['Segment_ID','Chromosome','Major_CN','Minor_CN']].drop_duplicates()   
    base_cn_table['Segment_Start'] = base_cn_table['Segment_ID'].apply(lambda x: int(x.split('-')[1]))
    base_cn_table['Segment_End'] = base_cn_table['Segment_ID'].apply(lambda x: int(x.split('-')[2]))
    
    base_mutation_table.to_csv(f'{output_dir}/mutation_table.tsv',sep="\t",index=False)
    segment_metadata.to_csv(f'{output_dir}/segment_metadata.tsv',sep="\t",index=False)
    processed_metadata = get_processed_metadata(segment_metadata)
    processed_metadata['NRPCC'] = nrpcc
    if len(segment_metadata.index)==0:
        raise ValueError('No segments in sim template')
    
   
   
    purity_sample_id = f'{sample_id}_purity_{purity:.2f}'

    
    segment_metadata['Dataset_ID'] = dataset_sample_id
    segment_metadata['Dataset'] = dataset
    
    wgd_override_store = RNG.permutation([True])
    apply_reads_correction_store = RNG.permutation([True])
    phasing = False
    apply_reads_correction = True
    wgd_override = True
   
                

    run_mutation_table = base_mutation_table.drop(columns=['Segment_ID','Major_CN','Minor_CN','Total_CN','True_Multiplicity','True_VAF'])
    run_output_dir = f'{output_dir}/purity_{purity:.2f}_wgd_{wgd_override}_phasing_{phasing}_correction_{apply_reads_correction}'
    
    run_sample_id = f'{purity_sample_id}_wgd_{wgd_override}_phasing_{phasing}_correction_{apply_reads_correction}'
    sample = sampletools.Sample(run_mutation_table,base_cn_table,sample_data['Subclone_Table'],run_sample_id,purity,merge_cn=False,apply_reads_correction=apply_reads_correction)
    
    run_sample(sample,run_output_dir,wgd_override=wgd_override)
    processed_count+=1
    for apply_penalty in [True,False]:
        posterior_table_path = f'{run_output_dir}/{run_sample_id}_posterior_timing_table_penalty_{apply_penalty}.tsv'
        posterior_table = get_sample_run_posterior_table(posterior_table_path,purity,wgd_override,phasing,apply_reads_correction)
        
        #print(posterior_table.head(50))
        posterior_table = posterior_table.merge(processed_metadata,on=['Segment_ID','Gain_Index','Major_CN','Minor_CN'],how='inner')
        
        
        
        posterior_table_path = f'{output_dir}/{run_sample_id}_posterior_timing_table_penalty_{apply_penalty}.tsv'
        posterior_table.to_csv(posterior_table_path,sep="\t",index=False)
        processed_posterior_table = get_processed_posterior_table(posterior_table,posterior_table_path)
        processed_posterior_table_path = f'{output_dir}/{run_sample_id}_processed_posterior_timing_table_penalty_{apply_penalty}.tsv'
        processed_posterior_table.to_csv(processed_posterior_table_path,sep="\t",index=False)

    #remove run output dir
    os.system(f'rm -r {run_output_dir}')

