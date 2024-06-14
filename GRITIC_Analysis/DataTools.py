import pandas as pd
import numpy as np
import platform
from pathlib import Path
import hashlib
import sys
import bz2
import pickle

from collections import Counter

def apply_penalty_quick(route_table,l):

    penalty_table = route_table[['Sample_ID','Segment_ID','Route','Average_N_Events','Probability']].copy().drop_duplicates()
    penalty_table['Probability_Corrected'] = penalty_table.groupby(['Sample_ID','Segment_ID'], group_keys=False).apply(lambda g: np.multiply(g.Probability,np.exp(-l*g.Average_N_Events))/np.sum(np.multiply(g.Probability,np.exp(-l*g.Average_N_Events))))
    return penalty_table['Probability_Corrected'].to_numpy()

def load_timing_from_dict(segment_path):
    input_file = bz2.BZ2File(segment_path,'rb')
    timing_dict = pickle.load(input_file)
    input_file.close()
    return timing_dict

'''def load_route_table(wgd_status=True,nrows=None,cancer_types=None,apply_penalty=True,include_pcawg=True,include_hartwig=True,min_mutations=20):

    if include_hartwig:

        dataset = 'Hartwig'
        hartwig_route_table = pd.read_csv(f'../in_data/{dataset}/complete_route_table_{dataset}.tsv',sep="\t",nrows=nrows,dtype={'Chromosome': str,'Node_Phasing':str})

            
        hartwig_route_table.loc[:,'Dataset'] = 'Hartwig'
    else:
        hartwig_route_table = pd.DataFrame()
    
    if include_pcawg:
       
        dataset = 'PCAWG'
        pcawg_route_table =  pd.read_csv(f'../in_data/{dataset}/complete_route_table_{dataset}.tsv',sep="\t",nrows=nrows, dtype={'Chromosome': str,'Node_Phasing':str})
        pcawg_route_table.loc[:,'Dataset'] = 'PCAWG'
    else:
        pcawg_route_table = pd.DataFrame()
    
    
    complete_route_table = pd.concat([hartwig_route_table,pcawg_route_table])
    
    route_cols = ['Sample_ID','Segment_ID','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN','Route','N_Mutations' ,'Probability', 'Average_N_Events', 'Average_Pre_WGD_Losses', 'Average_Post_WGD_Losses','Dataset','WGD_Status','WGD_Timing','WGD_Timing_CI_Low','WGD_Timing_CI_High']
    complete_route_table = complete_route_table[route_cols].drop_duplicates()
    baseline_cn = np.where(complete_route_table['WGD_Status'],2.1,1.1)
    complete_route_table = complete_route_table[complete_route_table['Major_CN']>baseline_cn]
    complete_route_table = complete_route_table[complete_route_table['N_Mutations']>=min_mutations]
    
    if wgd_status:
        complete_route_table = complete_route_table[complete_route_table['WGD_Status']]
    #convert cols to int
    int_cols = ['Segment_Start','Segment_End','Major_CN','Minor_CN']
    complete_route_table[int_cols] = complete_route_table[int_cols].astype(int)
    
    complete_route_table.loc[:,'Segment_Width'] = complete_route_table['Segment_End'] - complete_route_table['Segment_Start']
    
    classifications = load_classifications()

    #merge classifications
    complete_route_table = complete_route_table.merge(classifications[['Sample_ID','Cancer_Type']], on=['Sample_ID'],how='inner')
    
    
    if cancer_types is not None:
        complete_route_table = complete_route_table[complete_route_table['Cancer_Type'].isin(cancer_types)]
    
    
    if apply_penalty:
        penalty_corrected_probabilities = apply_penalty_quick(complete_route_table,2.7)
        complete_route_table['Probability'] = penalty_corrected_probabilities

    good_samples = get_good_samples()
    complete_route_table = complete_route_table[complete_route_table['Sample_ID'].isin(good_samples)]

    chromosomes = list(map(str,range(1,23))) + ['X']
    complete_route_table = complete_route_table[complete_route_table['Chromosome'].isin(chromosomes)]

    wgd_status_calls = get_wgd_calling_info()
    wgd_status_calls = wgd_status_calls[wgd_status_calls['Sample_ID'].isin(good_samples)]
    bad_samples = wgd_status_calls[((wgd_status_calls['Major_CN_Mode']==2) & (~wgd_status_calls['WGD_Status']))]['Sample_ID'].values
    
    complete_route_table = complete_route_table[~complete_route_table['Sample_ID'].isin(bad_samples)]
    return complete_route_table'''

def load_route_table(apply_penalty,cancer_types=None):
    route_table = pd.read_csv(f'../in_data/route_table_filtered_penalty_{apply_penalty}.tsv',sep="\t",dtype={'Chromosome': str,'Node_Phasing':str})
    if cancer_types is not None:
        route_table = route_table[route_table['Cancer_Type'].isin(cancer_types)]
    return route_table

#function with full data
'''def load_posterior_table(include_pcawg=True,include_hartwig=True,apply_penalty=False,min_mutations=20):
    if include_pcawg:
        dataset = 'PCAWG'
        pcawg_posterior_table = pd.read_csv(f'../in_data/{dataset}/complete_posterior_table_penalty_{apply_penalty}_{dataset}.tsv',sep="\t",dtype={'Chromosome': str})
        pcawg_posterior_table.loc[:,'Dataset'] = 'PCAWG'
    else:
        pcawg_posterior_table = pd.DataFrame()
    if include_hartwig:

        dataset ='Hartwig'
        hartwig_posterior_table = pd.read_csv(f'../in_data/{dataset}/complete_posterior_table_penalty_{apply_penalty}_{dataset}.tsv',sep="\t",dtype={'Chromosome': str})
        hartwig_posterior_table.loc[:,'Dataset'] = 'Hartwig'
    else:
        hartwig_posterior_table = pd.DataFrame()

    complete_posterior_table = pd.concat([pcawg_posterior_table,hartwig_posterior_table])

    complete_posterior_table = complete_posterior_table[complete_posterior_table['N_Mutations']>=min_mutations]
    
    
    int_cols = ['Segment_Start','Segment_End','Major_CN','Minor_CN']
    complete_posterior_table[int_cols] = complete_posterior_table[int_cols].astype(int)
    complete_posterior_table.loc[:,'Segment_Width'] = complete_posterior_table['Segment_End'] - complete_posterior_table['Segment_Start']
    complete_posterior_table.loc[:,'WGD_Status'] = ~complete_posterior_table['WGD_Timing'].isnull()


    classifications = load_classifications()

    #merge classifications
    complete_posterior_table = complete_posterior_table.merge(classifications[['Sample_ID','Cancer_Type']], on=['Sample_ID'],how='inner')
    good_samples = get_good_samples()
    complete_posterior_table = complete_posterior_table[complete_posterior_table['Sample_ID'].isin(good_samples)]
    chromosomes = list(map(str,range(1,23))) + ['X']
    complete_posterior_table = complete_posterior_table[complete_posterior_table['Chromosome'].isin(chromosomes)]

    wgd_status_calls = get_wgd_calling_info()
    wgd_status_calls = wgd_status_calls[wgd_status_calls['Sample_ID'].isin(good_samples)]
    bad_samples = wgd_status_calls[((wgd_status_calls['Major_CN_Mode']==2) & (~wgd_status_calls['WGD_Status']))]['Sample_ID'].values
    complete_posterior_table = complete_posterior_table[~complete_posterior_table['Sample_ID'].isin(bad_samples)]

    return complete_posterior_table'''

def load_posterior_table(apply_penalty):
    posterior_table = pd.read_csv(f'../in_data/posterior_table_filtered_penalty_{apply_penalty}.tsv',sep="\t",dtype={'Chromosome': str})
    return posterior_table


def load_genome_bins(bin_size):
    genome = pd.read_csv('../resources/chrom_arm_positions.tsv',sep="\t")
    genome = genome.groupby('Chromosome').agg({'Arm_Start':'min','Arm_End':'max'}).reset_index().rename(columns={'Arm_Start':'Chromosome_Start','Arm_End':'Chromosome_End'})
    #filter X,Y,MT
    genome = genome[~genome['Chromosome'].isin(['X','Y','MT'])]
    bin_store = []
    for index,row in genome.iterrows():
        bins = np.arange(row['Chromosome_Start'],row['Chromosome_End'],bin_size)
        bin_store.append(pd.DataFrame({'Chromosome':row['Chromosome'],'Bin_Start':bins[:-1],'Bin_End':bins[1:]}))
    
    bins = pd.concat(bin_store)
    #chromosome is string
    bins.loc[:,'Chromosome'] = bins['Chromosome'].astype(str)
    bins.loc[:,'Bin_ID'] = bins['Chromosome'] + '_' + bins['Bin_Start'].astype(str) + '_' + bins['Bin_End'].astype(str)
    return bins
def load_genome_arm_bins():
    genome = pd.read_csv('../resources/chrom_arm_positions.tsv',sep="\t")
    genome = genome[genome['Arm'].isin(['p','q'])]
    genome.loc[:,'Bin_ID'] = genome['Chromosome'] + genome['Arm']
    #filter X,Y,MT
    genome = genome[~genome['Chromosome'].isin(['Y','MT'])]
    genome.loc[:,'Chromosome'] = genome['Chromosome'].astype(str)
    #rename Arm_Start, Arm_End to Bin_Start, Bin_End
    genome = genome.rename(columns={'Arm_Start':'Bin_Start','Arm_End':'Bin_End'})
    return genome

def load_chromosome_sizes():
    arm_bins = load_genome_arm_bins()
    chromosome_bins = arm_bins.groupby('Chromosome').agg({'Bin_Start':'min','Bin_End':'max'}).reset_index()
    chromosome_bins.loc[:,'Chromosome'] = chromosome_bins['Chromosome'].astype(str)
    #rename Bin_Start, Bin_End to Chromosome_Start, Chromosome_End
    chromosome_bins = chromosome_bins.rename(columns={'Bin_Start':'Chromosome_Start','Bin_End':'Chromosome_End'})
    chromosome_bins.loc[:,'Chromosome_Width'] = chromosome_bins['Chromosome_End'] - chromosome_bins['Chromosome_Start']
    return chromosome_bins


def produce_joint_mapping(present_samples,group_cols,min_threshold=0.8):
    cancer_type_mapping ={}
    for group_vars,group in present_samples.groupby(group_cols):
        joint_cancer_type = [tuple(x) for x in group[['Cancer_Type_Joint','Cancer_Type_Code_Joint']].values]
        joint_cancer_type_counts = Counter(joint_cancer_type)
        
        joint_cancer_type_counts_norm = {k:v/len(group) for k,v in joint_cancer_type_counts.items()}
        max_group,max_proportion = max(joint_cancer_type_counts_norm.items(), key=lambda x: x[1])
        if max_proportion >= min_threshold:
            cancer_type_mapping[group_vars] = max_group
    return cancer_type_mapping
def get_missing_mapping(classifications,full_samples,name):
    classifications = classifications.copy().rename(columns={'Cancer_Type':'Cancer_Type_Joint','Cancer_Type_Code':'Cancer_Type_Code_Joint'}).drop(columns=['Donor_ID'])
    
    #get set of columns that can group by cancer type according to the original cohort metadata
    group_cols = list(full_samples.columns)
    for col in ['Sample_ID','Donor_ID']:
        if col in group_cols:
            group_cols.remove(col)
    for col in group_cols:
        full_samples[col] = full_samples[col].astype(str).str.strip()
    missing_samples = full_samples[~full_samples['Sample_ID'].isin(classifications['Sample_ID'])]

    sample_mapping = pd.read_csv(f'../resources/cancer_type_unified_{name}.txt',sep="\t",encoding='latin-1').fillna('')

    sample_mapping = {tuple([row[x].strip() for x in group_cols]):row['Cancer_Type_Overall'] for index,row in sample_mapping.iterrows()}

    remapped_samples = []

    for group,group_samples in missing_samples.groupby(group_cols):
        
        if group in sample_mapping:
            
            cancer_type = sample_mapping[group]
            if cancer_type == '':
                continue
            group_samples['Cancer_Type_Joint'] = cancer_type
            remapped_samples.append(group_samples.copy())
  
    
    code_mapping = classifications[['Cancer_Type_Joint','Cancer_Type_Code_Joint']].drop_duplicates()

    remapped_samples = pd.concat(remapped_samples).merge(code_mapping,on=['Cancer_Type_Joint'],how='inner')
    
    remapped_samples = remapped_samples.drop(columns=group_cols).rename(columns={'Cancer_Type_Joint':'Cancer_Type'})
    
    still_missing = missing_samples[~missing_samples['Sample_ID'].isin(remapped_samples['Sample_ID'])]
    print('mapped',len(remapped_samples),'still missing',len(still_missing),'for',name)
    return remapped_samples

#loading cancer types with full data
'''def load_cancer_type_classifications():
    classifications = pd.read_csv('/camp/home/bakert/secure-working/GRITIC/Analysis/resources/pcawg_hartwig_unified_classifications.tsv',sep="\t")
    classifications['Sample_ID'] = np.where(classifications['cohort']=='PCAWG',classifications['icgc_aliquot_id'],classifications['sample_id'])
    
    classifications = classifications[['Sample_ID','patient_id','cancer_type','cancer_type_code']].rename(columns={'patient_id':'Donor_ID','cancer_type':'Cancer_Type','cancer_type_code':'Cancer_Type_Code'})

    pcawg_table = pd.read_csv('/PCAWG_PATH/ICGC_annotations/summary_table_combined_annotations_v2.txt',sep="\t")
    
    pcawg_table = pcawg_table[['samplename','icgc_donor_id','cancer_type','histology_abbreviation']].rename(columns={'samplename':'Sample_ID','icgc_donor_id':'Donor_ID','cancer_type':'Cancer_Type','histology_abbreviation':'Histology'})
    hartwig_metadata = pd.read_csv('/HARTWIG_PATH/metadata/metadata.tsv',sep="\t").rename(columns={'sampleId':'Sample_ID','hmfPatientId':'Donor_ID'})[['Sample_ID','Donor_ID','primaryTumorLocation','primaryTumorSubLocation','primaryTumorType','primaryTumorSubType']].fillna('')
    missing_pcawg = get_missing_mapping(classifications,pcawg_table,'PCAWG')
    
    missing_hartwig = get_missing_mapping(classifications,hartwig_metadata,'Hartwig')
    
    classifications = pd.concat([classifications,missing_pcawg,missing_hartwig])
    return classifications
def load_classifications(wgd_only=False):

    cancer_type_classifications = load_cancer_type_classifications()

    wgd_status_table = load_wgd_status()
    classifications_with_wgd_status = cancer_type_classifications.merge(wgd_status_table,on=['Sample_ID'],how='inner')
    
    
    if wgd_only:
        raise NotImplementedError('wgd_only not implemented')
       
    good_samples = get_good_samples()
    classifications_with_wgd_status = classifications_with_wgd_status[classifications_with_wgd_status['Sample_ID'].isin(good_samples)]
    
    wgd_status_calls = get_wgd_calling_info()
    wgd_status_calls = wgd_status_calls[wgd_status_calls['Sample_ID'].isin(good_samples)]
    bad_samples = wgd_status_calls[((wgd_status_calls['Major_CN_Mode']==2) & (~wgd_status_calls['WGD_Status']))]['Sample_ID'].values
    classifications_with_wgd_status = classifications_with_wgd_status[~classifications_with_wgd_status['Sample_ID'].isin(bad_samples)]
    
    return classifications_with_wgd_status'''

def load_classifications():
    classifications = pd.read_csv('../in_data/cancer_type_classifications.tsv',sep="\t")
    return classifications
def load_wgd_status():
    
    wgd_table_hartwig_path = '../in_data/Hartwig/complete_wgd_calling_info_Hartwig.tsv'
    wgd_table_hartwig = pd.read_csv(wgd_table_hartwig_path,sep="\t")
    wgd_table_pcawg_path = '../in_data/PCAWG/complete_wgd_calling_info_PCAWG.tsv'
    wgd_table_pcawg = pd.read_csv(wgd_table_pcawg_path,sep="\t")
    wgd_table = pd.concat([wgd_table_hartwig,wgd_table_pcawg])
    #wgd_status_table = wgd_table[['Sample_ID','WGD_Status']]
    return wgd_table

def process_pcawg_cn_table(cn_table):
    filter_cols = ['chromosome','start','end','major_cn','minor_cn']
    cn_table = cn_table[filter_cols]
    cn_table = cn_table.rename(columns={'chromosome':'Chromosome', 'start':'Segment_Start', 'end':'Segment_End', 'major_cn':'Major_CN', 'minor_cn':'Minor_CN'})
    cn_table.loc[:,'Chromosome'] = cn_table['Chromosome'].astype(str)
    cn_table = cn_table.dropna().reset_index(drop=True)
    #Major_CN, Minor_CN are floats, convert to int
    cn_table.loc[:,'Major_CN'] = cn_table['Major_CN'].astype(int)
    cn_table.loc[:,'Minor_CN'] = cn_table['Minor_CN'].astype(int)
    return cn_table
def load_pcawg_copy_number():
    input_dir = Path('/PCAWG_PATH/ICGC_consensus_copynumber/20170119_release')
    cn_tables = []
    for path in input_dir.glob('*.txt'):
        cn_table = pd.read_csv(path,sep="\t")
        sample_id = path.stem.split('.')[0]
        cn_table = process_pcawg_cn_table(cn_table)
        cn_table.loc[:,'Sample_ID'] = sample_id
        cn_tables.append(cn_table)
    copy_number = pd.concat(cn_tables).reset_index(drop=True)
    return copy_number

def load_sample_copy_number(include_hartwig=True,include_pcawg=True):
    cn_tables =[]
    if include_hartwig:
        hartwig_table = pd.read_csv('/HARTWIG_PATH/hartwig_copy_number.tsv',sep="\t").drop(columns=['Tumor_Type'])
        cn_tables.append(hartwig_table.reset_index(drop=True))
    if include_pcawg:
        pcawg_table = load_pcawg_copy_number().reset_index(drop=True)
        cn_tables.append(pcawg_table)
    if len(cn_tables) == 0:
        raise ValueError("No copy number tables selected")
    copy_number = pd.concat(cn_tables)
    #segment start and end are floats, convert to int
    copy_number.loc[:,'Segment_Start'] = copy_number['Segment_Start'].astype(int)
    copy_number.loc[:,'Segment_End'] = copy_number['Segment_End'].astype(int)
    #set width
    copy_number.loc[:,'Segment_Width'] = copy_number['Segment_End'] - copy_number['Segment_Start']
    return copy_number
#loading wgd calling info with full data
'''def get_wgd_calling_info():
    hartwig_wgd_calling_info = pd.read_csv('../in_data/Hartwig/complete_wgd_calling_info_Hartwig.tsv',sep="\t")
    pcawg_wgd_calling_info = pd.read_csv('../in_data/PCAWG/complete_wgd_calling_info_PCAWG.tsv',sep="\t")
    wgd_calling_info = pd.concat([hartwig_wgd_calling_info,pcawg_wgd_calling_info])
    return wgd_calling_info'''

def get_wgd_calling_info():
    wgd_calling_info = pd.read_csv('../in_data/wgd_calling_info.tsv',sep="\t")
    return wgd_calling_info

def get_bootstrap_table(sample_dict):
    RNG = np.random.default_rng()
    sample_ids = sorted(list(sample_dict.keys()))
    bootstrap_ids = RNG.choice(sample_ids,size=len(sample_ids),replace=True)
    samples = []
    for i,sample_id in enumerate(bootstrap_ids):
        sample_table = sample_dict[sample_id].copy()
        sample_table.loc[:,'Sample_ID'] = f'{sample_id}_{i}'
        samples.append(sample_table)
    return pd.concat(samples)
def get_sample_dict(sample_table):
    sample_dict = {}
    for sample_id,sample_table in sample_table.groupby('Sample_ID'):
        sample_dict[sample_id] = sample_table
    return sample_dict
#getting non-filtered samples with full data
'''def get_good_samples(nrpcc_cutoff=5):
    hartwig_metadata = pd.read_csv('/HARTWIG_PATH/metadata/metadata.tsv',sep="\t").rename(columns={'sampleId':'Sample_ID','hmfPatientId':'Donor_ID'})[['Sample_ID','Donor_ID']]
    pcawg_metadata = pd.read_csv('/PCAWG_PATH/ICGC_annotations/summary_table_combined_annotations_v2.txt',sep="\t").rename(columns={'samplename':'Sample_ID','icgc_donor_id':'Donor_ID'})
    n_samples = len(hartwig_metadata) + len(pcawg_metadata)
    
    pcawg_nrpcc_table = pd.read_csv('../in_data/PCAWG/complete_nrpcc_table_PCAWG.tsv',sep="\t")
    hartwig_nrpcc_table = pd.read_csv('../in_data/Hartwig/complete_nrpcc_table_Hartwig.tsv',sep="\t")
    combined_nrpcc_table = pd.concat([pcawg_nrpcc_table,hartwig_nrpcc_table])
    
    #assert set(pcawg_metadata['Sample_ID']).symmetric_difference(set(pcawg_nrpcc_table['Sample_ID'])) == set()
    assert set(hartwig_metadata['Sample_ID']).symmetric_difference(set(hartwig_nrpcc_table['Sample_ID'])) == set()
    #filter out non-preferred multisample PCAWG cases
    pcawg_metadata = pcawg_metadata.merge(pcawg_nrpcc_table,on=['Sample_ID'],how='inner')
    is_preferred_pcawg = pcawg_metadata.groupby('Donor_ID')['is_preferred'].any()
    assert is_preferred_pcawg.all()

    pcawg_metadata = pcawg_metadata[pcawg_metadata['is_preferred']][['Sample_ID','Donor_ID']]

   
   
    
    hartwig_metadata = hartwig_metadata.merge(hartwig_nrpcc_table,on=['Sample_ID'],how='inner')
    pcawg_metadata = pcawg_metadata.merge(pcawg_nrpcc_table,on=['Sample_ID'],how='inner')
    
    metadata = pd.concat([hartwig_metadata,pcawg_metadata])
   
    assert len(metadata['Sample_ID'].unique()) == len(metadata)
    n_pre_filtered_samples = len(metadata)
    
    metadata = metadata[metadata['NRPCC']>=nrpcc_cutoff].reset_index()
    x = len(metadata)

    metadata = metadata.sort_values(by=['Donor_ID','NRPCC'],ascending=[True,False]).drop_duplicates(subset=['Donor_ID'],keep='first').reset_index(drop=True)

    samples_filtered = n_samples - len(metadata)
   
    good_samples = list(metadata['Sample_ID'].unique())
    return good_samples'''

def get_good_samples():
    with open('../in_data/good_samples.txt','r') as f:
        good_samples = f.read().splitlines()
    return good_samples
    
