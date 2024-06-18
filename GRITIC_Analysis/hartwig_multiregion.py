import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#font size 18, pdf font type 42, nimbus sans font
plt.rcParams.update({'font.size': 18, 'pdf.fonttype': 42, 'font.family': 'Nimbus Sans'})

import DataTools

from collections import Counter

def load_nrpcc_table():
    hartwig_nrpcc_table = pd.read_csv('../in_data/Hartwig/complete_nrpcc_table_Hartwig.tsv',sep="\t")
    return hartwig_nrpcc_table
def load_classifications():
    classifications = pd.read_csv('../in_data/cancer_type_classifications_all_samples.tsv',sep="\t")
    return classifications[['Donor_ID','Sample_ID','Cancer_Type']]

def load_chr_sizes():
    chr_arm_positions_path = '../resources/chrom_arm_positions.tsv'
    chr_arm_positions = pd.read_csv(chr_arm_positions_path,sep='\t')
    chr_sizes = chr_arm_positions.groupby(['Chromosome'])['Arm_End'].max().reset_index()
    chr_sizes['Chromosome'] = chr_sizes['Chromosome'].astype(str)
    return chr_sizes
def load_timing_table(apply_penalty,min_coverage=0.3):
    #read in the data
    posterior_path = f'../in_data/Hartwig/complete_posterior_table_penalty_{apply_penalty}.tsv'
    posterior_table = pd.read_csv(posterior_path,sep='\t',dtype={'Chromosome':str})
    min_cn = np.where(posterior_table['WGD_Timing'].isnull(),1,2)
    posterior_table = posterior_table[posterior_table['Major_CN']>=min_cn]
    posterior_table = posterior_table[posterior_table['N_Mutations']>=20]
    posterior_table = posterior_table.groupby(['Sample_ID','Segment_ID','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN','Posterior_Sample_Index'])['Gain_Timing'].min().reset_index()
    posterior_table['Segment_Width'] = posterior_table['Segment_End']-posterior_table['Segment_Start']

    chr_sizes = load_chr_sizes()
    posterior_chromosome_coverage = posterior_table[['Sample_ID','Chromosome','Segment_ID','Segment_Width']].drop_duplicates().groupby(['Sample_ID','Chromosome']).apply(lambda x : np.sum(x['Segment_Width'])).reset_index().rename(columns={0:'Chromosome_Coverage'})
    posterior_chromosome_coverage = posterior_chromosome_coverage.merge(chr_sizes,on='Chromosome')
    posterior_chromosome_coverage['Chromosome_Coverage'] = posterior_chromosome_coverage['Chromosome_Coverage']/posterior_chromosome_coverage['Arm_End']
    posterior_chromosome_coverage = posterior_chromosome_coverage[posterior_chromosome_coverage['Chromosome_Coverage']>=min_coverage]
    posterior_table = posterior_table.merge(posterior_chromosome_coverage[['Sample_ID','Chromosome']],on=['Sample_ID','Chromosome'],how='inner')
    posterior_table = posterior_table.groupby(['Sample_ID','Chromosome']).apply(lambda x : np.average(x['Gain_Timing'],weights=x['Segment_Width'])).reset_index().rename(columns={0:'Gain_Timing'})
    return posterior_table

def run_background_test(sample_timing_table,test_chrs):
    timing_test_chr = sample_timing_table[sample_timing_table['Chromosome'].isin(test_chrs)]['Gain_Timing'].values
    timing_background_chr = sample_timing_table[~sample_timing_table['Chromosome'].isin(test_chrs)]['Gain_Timing']
    return np.mean(timing_test_chr-np.mean(timing_background_chr))

def run_comparison(timing_table,sample_ids):
    comparison_store = []
    sample_timing_tables = {}
    valid_samples = []
    for sample_id in sample_ids:
        sample_timing_tables[sample_id] = timing_table[timing_table['Sample_ID']==sample_id]
        if len(sample_timing_tables[sample_id])>0:
            valid_samples.append(sample_id)
    if len(valid_samples)<2:
        return comparison_store
    sample_chrs = {sample_id:set(sample_timing_tables[sample_id]['Chromosome']) for sample_id in sample_ids}
    all_chrs = [chr for sample_id in sample_ids for chr in sample_chrs[sample_id]]
    shared_chrs = set([chr for chr,count in Counter(all_chrs).items() if count>1])

    for sample_id in sample_ids:
        test_chrs = sample_chrs[sample_id] - shared_chrs
        background_chrs = sample_chrs[sample_id] - test_chrs
        if len(test_chrs)==0 or len(background_chrs)==0:
            continue
        comparison_store.extend([run_background_test(sample_timing_tables[sample_id],test_chrs)])
    return comparison_store
def plot_comparison_store(comparison_store,apply_penalty):
    bins = np.linspace(-1,1,21)
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax.hist(comparison_store,bins=bins,color='#42398f',edgecolor='black')
    ax.set_xlabel('Timing of gains unique to sample - timing of gains shared with other samples')
    ax.set_ylabel('Frequency')
    fig.tight_layout()

    plt.savefig(f'../plots/multiregion_hartwig/gain_difference_apply_penalty_{apply_penalty}.pdf',bbox_inches='tight')

if __name__ == '__main__':

    nrpcc_table = load_nrpcc_table()
    nrpcc_table = nrpcc_table[nrpcc_table['NRPCC']>=5]
    wgd_status_calls = DataTools.get_wgd_calling_info()
    bad_samples = wgd_status_calls[((wgd_status_calls['Major_CN_Mode']==2) & (~wgd_status_calls['WGD_Status']))]['Sample_ID'].values

    for apply_penalty in [True,False]:
        timing_table = load_timing_table(apply_penalty)
        timing_table = timing_table.merge(nrpcc_table,on='Sample_ID')
        timing_table = timing_table[~timing_table['Sample_ID'].isin(bad_samples)]

        comparison_store = []
        classifications = load_classifications()
        for donor_id,donor_table in classifications.groupby('Donor_ID'):
            sample_ids = list(donor_table['Sample_ID'])
            
            if len(sample_ids)<2:
                continue

            chr_difference_timing = run_comparison(timing_table,sample_ids)
            comparison_store.extend(chr_difference_timing)
        print(apply_penalty,np.mean(np.array(comparison_store)>0),np.mean(comparison_store),np.sum(np.array(comparison_store)>0),len(comparison_store))
        plot_comparison_store(comparison_store,apply_penalty)

        

    
