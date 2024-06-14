import numpy as np
import pandas as pd
import string
from gritic import gritictimer,sampletools
import os
import sys
sys.path.insert(0,'../code')
import DataTools

import PCAWGDataLoader
import HartwigDataLoader
import SampleSimulator

import argparse

from scipy.stats import kstest
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



    def __init__(self) -> None:
        pass
    def get_random_gain_timing(self,n_timings):
        return RNG.uniform(0.05,0.95,size=n_timings)



def run_sample(sample,output_dir,wgd_override=None):

    gritictimer.process_sample(sample,output_dir,plot_trees=False,wgd_override=wgd_override)

def get_valid_complex_states(max_routes_per_state=500):
    wgd_routes = pd.read_csv('../in_data/cn_states_n_routes.tsv',sep="\t").rename(columns={'Route_ID':'Route'})
    
    wgd_routes = wgd_routes[wgd_routes['WGD_Status']].copy()
    
    wgd_routes = wgd_routes[(wgd_routes['Major_CN']>=3)].copy()
    #wgd_routes = wgd_routes[(wgd_routes['Major_CN']<=3) & (wgd_routes['Minor_CN']<=0)].copy()
    n_routes = wgd_routes.groupby(['Major_CN','Minor_CN']).size().reset_index().rename(columns={0:'N_Routes'})
    cn_states = n_routes[n_routes['N_Routes']<=max_routes_per_state][['Major_CN','Minor_CN']].drop_duplicates()
    return cn_states
def make_sim_template_wgd(mutation_table,parsimony_only=False,major_gains_wgd=False,uniform_sampling=False):
    if parsimony_only and major_gains_wgd:
        raise ValueError('Cannot have both parsimony only and major gains wgd')
    wgd_routes = pd.read_csv('../in_data/all_routes_n_events.tsv',sep="\t")
    wgd_routes = wgd_routes[wgd_routes['WGD_Status']]

    if parsimony_only:
        wgd_routes = wgd_routes[wgd_routes['Parsimonious']]
    if major_gains_wgd:
        wgd_routes = wgd_routes[wgd_routes['Major_Gains_WGD']]

    mutation_table.loc[:,'Total_Coverage'] = mutation_table['Tumor_Ref_Count']+mutation_table['Tumor_Alt_Count']
 
    cn_table_summary = mutation_table.groupby(['Segment_ID','Major_CN','Minor_CN','Chromosome','Segment_Start','Segment_End']).agg({'Position':'count','Total_Coverage':'mean'}).reset_index()
    cn_table_summary.rename(columns={'Position':'N_Mutations','Total_Coverage':'Adjusted_Coverage'},inplace=True)
    cn_table_summary.loc[:,'Adjusted_Coverage'] = cn_table_summary['Adjusted_Coverage'].round().astype(int)
    
    complex_segments = cn_table_summary[(cn_table_summary['N_Mutations']>=10) & (cn_table_summary['Major_CN']>2) ]

    if len(complex_segments.index) == 0:
        return complex_segments
    valid_complex_states = get_valid_complex_states()
  
    complex_segments = complex_segments.merge(valid_complex_states,how='inner',on=['Major_CN','Minor_CN'])
    complex_segments = complex_segments.merge(wgd_routes,how='inner').drop(columns=['Average_N_Events'])
    if len(complex_segments.index) == 0:
        return complex_segments
    complex_segments['Segment_ID'] = complex_segments['Segment_ID'] + '|' + pd.Series(np.arange(len(complex_segments.index)).astype(str))
    if parsimony_only:
        
        complex_segments = complex_segments[complex_segments['Parsimonious']]

    if not uniform_sampling:
        complex_segments['Segment_ID'] = complex_segments['Segment_ID'].str.split('|').str[0]
        complex_segments = complex_segments.groupby('Segment_ID').apply(lambda x: x.sample(n=1,random_state=RNG)).reset_index(drop=True)

    #complex_segments = complex_segments.drop(columns=['Parsimonious'])
    
    
    wgd_segments = cn_table_summary[(cn_table_summary['N_Mutations']>=10) & (cn_table_summary['Major_CN']==2) ].copy()
    if len(wgd_segments) == 0:
        return wgd_segments
    if len(wgd_segments.index)>0:
        wgd_segments.loc[:,'Route'] = np.nan
    
    sim_template = pd.concat([wgd_segments,complex_segments],ignore_index=True)

    sim_template.loc[:,'Gain_Timing_Target'] = pd.Series([np.array([])]*len(sim_template.index),index=sim_template.index)
    return sim_template

def load_wgd_dist(dataset,sample_id):
    if dataset =='PCAWG':
        wgd_timing_path= '/GRITIC_PCAWG_Output/complete_route_table_PCAWGe.tsv'
    if dataset =='Hartwig':
        wgd_timing_path = '/GRITIC_Hartwig_Output/complete_route_table_Hartwi.tsv'
    wgd_timing_table = pd.read_csv(wgd_timing_path,sep="\t",usecols=['Sample_ID','WGD_Timing','WGD_Timing_CI_Low','WGD_Timing_CI_High']).drop_duplicates()
    wgd_timing_table = wgd_timing_table.dropna()
    if sample_id =='Random':
        sample_id = RNG.choice(wgd_timing_table['Sample_ID'].unique())
    if not sample_id in wgd_timing_table['Sample_ID'].values:
        raise ValueError("Sample ID doesn't have WGD timing")
    wgd_timing_row = wgd_timing_table[wgd_timing_table['Sample_ID']==sample_id].iloc[0]
    std = np.maximum(wgd_timing_row['WGD_Timing_CI_High']-wgd_timing_row['WGD_Timing'],wgd_timing_row['WGD_Timing']-wgd_timing_row['WGD_Timing_CI_Low'])/1.6
    return wgd_timing_row['WGD_Timing'],std

def resample_metadata(segment_metadata,n_tests=1000):
    segment_metadata = segment_metadata.copy()
    wgd_timing = segment_metadata['WGD_Timing'].iloc[0]
    segment_metadata['Base_Segment_ID'] = segment_metadata['Segment_ID'].str.split('|').str[:-1].str.join('_')
    segment_metadata['Gain_Timing_Array'] = segment_metadata['Gain_Timing'].apply(lambda x: np.array(list(eval(x).values())))
    segment_metadata['Segment_Width'] = segment_metadata['Segment_End']-segment_metadata['Segment_End']
    seg_store = [[] for _ in range(n_tests)]

    wgd_timing = segment_metadata['WGD_Timing'].iloc[0]
    for base_segment_id,segment_group in segment_metadata.groupby('Base_Segment_ID'):
        segment_group = segment_group.copy()
        segment_group['Weight'] = np.where(segment_group['Gain_Timing_Array'].apply(lambda x: np.any(x<wgd_timing)),wgd_timing,1-wgd_timing)
        segment_group['Weight'] = np.sqrt(segment_group['Weight'])
        for i in range(n_tests):
            sampled_segment_group = segment_group.sample(n=1,random_state=RNG,weights=segment_group['Weight'])
            seg_store[i].append(sampled_segment_group)
   
    scores = []
    for i in range(n_tests):
        seg_store[i] = pd.concat(seg_store[i])
        gain_timing = np.concatenate(seg_store[i]['Gain_Timing_Array'].to_numpy())
        
        weights = np.concatenate([[wgd_timing]*x.size for x in seg_store[i]['Gain_Timing_Array'].values])
        
        prop_before_wgd = np.average(gain_timing<wgd_timing,weights=weights)

        if np.abs(prop_before_wgd-wgd_timing)>0.1:
            score = 0
        else:
            score = kstest(gain_timing,'uniform',args=(0,1))[1]
            

        scores.append(score)
    scores = np.array(scores)
    scores[scores <0.05] = 0
    if np.isclose(np.sum(scores),0):
        print('all scores are zero')
        raise ValueError('All scores are zero')
    #chosen_index =RNG.choice(np.arange(n_tests),p=scores/np.sum(scores))
    chosen_index = np.argmax(scores)
    segment_metadata = seg_store[chosen_index]
    
    segment_metadata = segment_metadata.drop(columns=['Base_Segment_ID','Gain_Timing_Array'])
    
    return segment_metadata.copy()
class Random_Gain_Timing:
    def __init__(self) -> None:
        pass
    def get_random_gain_timing(self,n_timings):
        return RNG.uniform(0.05,0.95,size=n_timings)


if __name__ == '__main__':
    processed_count = 0
    my_parser = argparse.ArgumentParser(allow_abbrev=False)
    my_parser.add_argument('--base_output_dir', action='store', type=str, required=True)
    my_parser.add_argument('--parsimony_only', action='store', type=str2bool, required=False, default=False)
    my_parser.add_argument('--remove_subclones', action='store', type=str2bool, required=False, default=False)
    my_parser.add_argument('--sync_gains_test', action='store', type=str2bool, required=False, default=False)
    my_parser.add_argument('--uniform_sampling', action='store', type=str2bool, required=False, default=False)
    my_parser.add_argument('--apply_reads_correction', action='store', type=str2bool, required=False, default=True)
    my_parser.add_argument('--test_reads_correction', action='store', type=str2bool, required=False, default=False)
    my_parser.add_argument('--test_wgd_constraint', action='store', type=str2bool, required=False, default=False)
    my_parser.add_argument('--wgd_mode', action='store', type=str, required=False, default='WGD_Only')
    my_parser.add_argument('--max_n_samples', action='store', type=int, required=False, default=5000)
    

    args = my_parser.parse_args()
    os.makedirs(args.base_output_dir,exist_ok=True)
    random_gain_generator =  Random_Gain_Timing()
    sample_loader = SampleLoader()    
    while processed_count<50:
        if len(os.listdir(args.base_output_dir))>args.max_n_samples:
            print('reached max samples')
            break
        sample_id = SampleSimulator.generate_random_sample_id(length=12)

        try:
            dataset_sample_id,dataset,sample_data = sample_loader.get_random_sample()
        except ValueError:
            print('couldnt get sample data for',sample_id)
            continue
        if sample_data['Mutation_Table'] is None:
            continue
        
        major_cn_mode = sampletools.get_major_cn_mode_from_cn_table(sample_data['CN_Table'])
        if major_cn_mode != 2:
            continue
        
        if args.wgd_mode == 'WGD_Only':
            wgd_status = True
            if args.uniform_sampling:
                try:
                    wgd_mean,wgd_std = load_wgd_dist(dataset,dataset_sample_id)
                except ValueError:
                    continue
                
                wgd_timing = np.clip(RNG.normal(loc=wgd_mean,scale=wgd_std),0.01,0.99)
            else:
                wgd_timing = RNG.uniform(0.01,0.99)
                #wgd_timing = np.clip(RNG.normal(loc=wgd_mean,scale=wgd_std),0.01,0.99)
   
            
        
        if args.wgd_mode =='WGD_Call_Test':
            
            wgd_status = RNG.choice([True,False])
            if wgd_status:
                wgd_mean,wgd_std = load_wgd_dist(dataset,sample_id='Random')
                wgd_timing = np.clip(RNG.normal(loc=wgd_mean,scale=wgd_std),0,1)
                sample_id = f'{sample_id}_WGD'
            else:
                wgd_timing = None
                sample_id = f'{sample_id}_No_WGD'
        
        if major_cn_mode ==1 and args.wgd_mode =='WGD_Call_Test':
            continue
        if major_cn_mode >2 and wgd_status:
            major_gains_wgd = True
        else:
            major_gains_wgd = False
        
        if major_gains_wgd and args.parsimony_only:
            print('cannot have major gains WGD with parsimony only, skipping')
            continue
       
        
        sim_template = make_sim_template_wgd(sample_data['Mutation_Table'],parsimony_only=args.parsimony_only,major_gains_wgd=major_gains_wgd,uniform_sampling=args.uniform_sampling)
        
        if len(sim_template.index)==0:
            continue
        output_dir = f"{args.base_output_dir}/{sample_id}"
        os.makedirs(output_dir,exist_ok=True)
        sim_template.loc[:,'Purity'] = sample_data['Sample_Purity']
        sim_template.loc[:,'Sample_ID']= sample_id
        
        if wgd_status:
            sim_template.loc[:,'WGD_Timing'] = wgd_timing
        elif not wgd_status and args.wgd_mode == 'WGD_Call_Test':
            sim_template.loc[:,'WGD_Timing'] = random_gain_generator.get_random_gain_timing(n_timings=len(sim_template.index))
        else:
            sim_template.loc[:,'WGD_Timing'] = None
        
        if args.remove_subclones and not sample_data['Subclone_Table'] is None:
        
            subclone_proportion = sample_data['Subclone_Table']['Subclone_Fraction'].sum()
            sim_template.loc[:,'N_Mutations'] = np.round(sim_template['N_Mutations']*(1-subclone_proportion)).astype(int)
            sample_data['Subclone_Table'] = None
        
        if args.sync_gains_test:
            timing_samples = random_gain_generator.get_random_gain_timing(n_timings=2000)
            if RNG.random()<0.5:
                synchronous_sample = False
                sync_gain_timing = np.nan
                sync_array= [np.array([RNG.choice(timing_samples)]) for _ in range(len(sim_template.index))]
            else:
                sync_gain_timing = timing_samples[0]
                timing_samples = timing_samples[1:]
                synchronous_sample = True
                
                sync_array = []
                for _,row in sim_template.iterrows():
                    if row['Major_CN']==2:
                        sync_array.append(np.array([]))
                    else:
                        if RNG.random()<0.9:
                            sync_array.append(np.array([sync_gain_timing]))
                        else:
                            sync_array.append(np.array([RNG.choice(timing_samples)]))
            sim_template.loc[:,'Gain_Timing_Target'] = pd.Series(sync_array, index=sim_template.index)
        sim_template = sim_template[sim_template['Chromosome'].isin(list(map(str,range(1,23))) + ['X'])]
        sim_template.to_csv(f'{output_dir}/test_sim_template.tsv',sep="\t",index=False)
        
        if sample_data['Subclone_Table'] is not None:
            sample_data['Subclone_Table'].to_csv(f'{output_dir}/test_subclone_table.tsv',sep="\t",index=False)
        
        n_subclones = 0 if sample_data['Subclone_Table'] is None else len(sample_data['Subclone_Table'])
        
        mutation_table,segment_metadata = SampleSimulator.simulate_sample(sim_template,sample_data['Subclone_Table'],wgd_status=True,sex=sample_data['Sex'])
        segment_metadata['N_Subclones'] = n_subclones
        
        if args.uniform_sampling:
            try:
                segment_metadata = resample_metadata(segment_metadata)
                mutation_table = mutation_table[(mutation_table['Major_CN']<=2) | (mutation_table['Segment_ID'].isin(segment_metadata['Segment_ID']))]

                segment_metadata['Segment_ID'] = segment_metadata['Segment_ID'].astype(str).str.split('|').str[0]
        
                mutation_table['Segment_ID'] = mutation_table['Segment_ID'].str.split('|').str[0]
            except ValueError:
                continue
        else:
            pass
        
        mutation_table = mutation_table[(mutation_table['Major_CN']<=2) | (mutation_table['Segment_ID'].isin(segment_metadata['Segment_ID']))]
        

        
        mutation_table.to_csv(f'{output_dir}/test_mutation_table.tsv',sep="\t",index=False)
        segment_metadata.to_csv(f'{output_dir}/test_segment_metadata.tsv',sep="\t",index=False)
        
        if args.wgd_mode == 'WGD_Call_Test':
            mutation_table = mutation_table[mutation_table['Major_CN']==2]
        if len(segment_metadata.index)==0:
            continue

        segment_metadata.loc[:,'Dataset_ID'] = dataset_sample_id
        segment_metadata.loc[:,'Dataset'] = dataset
        if args.sync_gains_test:
            segment_metadata.loc[:,'Sync_Gain_Timing'] = sync_gain_timing
            segment_metadata.loc[:,'Synchronous_Sample'] = synchronous_sample
        if args.sync_gains_test and len(segment_metadata['Chromosome'].unique())<3:
            continue

        
        metadata_path = f'{output_dir}/{sample_id}_metadata.tsv'
        observed_purity = sample_data['Sample_Purity']
        segment_metadata.loc[:,'Observed_Purity'] = observed_purity

        
        if not args.test_wgd_constraint:
            output_dir = f"{args.base_output_dir}/{sample_id}"
            unphased_sample_id = f'{sample_id}'
            unphased_mutation_table = mutation_table.copy()
            unphased_mutation_table['Phasing'] = pd.NA
            unphased_sample = sampletools.Sample(unphased_mutation_table,sample_data['CN_Table'],sample_data['Subclone_Table'],sample_id,observed_purity)
            run_sample(unphased_sample,output_dir)
            
            segment_metadata.to_csv(metadata_path,sep="\t",index=False)
            if sample_data['Subclone_Table'] is not None:
                sample_data['Subclone_Table'].to_csv(f'{output_dir}/test_subclone_table.tsv',sep="\t",index=False)
            os.system(f'rm -r {output_dir}/{sample_id}_timing_dicts')

            if args.test_reads_correction:
                output_dir_no_correction = f"{args.base_output_dir}/{sample_id}_no_correction"
                no_correction_sample_id = f'{sample_id}_no_correction'
                no_correction_mutation_table =unphased_mutation_table.copy()
                no_correction_sample = sampletools.Sample(no_correction_mutation_table,sample_data['CN_Table'],sample_data['Subclone_Table'],no_correction_sample_id,observed_purity,apply_reads_correction=False)
                run_sample(no_correction_sample,output_dir_no_correction)
                os.system(f'rm -r {output_dir_no_correction}/{no_correction_sample_id}_timing_dicts')
        else:
            for wgd_override in [True,False]:
                output_dir = f"{args.base_output_dir}/{sample_id}_WGD_{wgd_override}"
                unphased_sample_id = f'{sample_id}_WGD_{wgd_override}'
                unphased_mutation_table = mutation_table.copy()
                unphased_mutation_table['Phasing'] = pd.NA
                unphased_sample = sampletools.Sample(unphased_mutation_table,sample_data['CN_Table'],sample_data['Subclone_Table'],unphased_sample_id,observed_purity)
                run_sample(unphased_sample,output_dir,wgd_override=wgd_override)
                segment_metadata.to_csv(f'{output_dir}/{unphased_sample_id}_metadata.tsv',sep="\t",index=False)
                if sample_data['Subclone_Table'] is not None:
                    sample_data['Subclone_Table'].to_csv(f'{output_dir}/test_subclone_table.tsv',sep="\t",index=False)
                os.system(f'rm -r {output_dir}/{unphased_sample_id}_timing_dicts')
       
                
