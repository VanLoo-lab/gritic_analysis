import os

import pandas as pd
import numpy as np

def get_combined_segments(all_sample_cn):
    combined_segment_data = {'Chromosome':[],'Segment_Start':[],'Segment_End':[]}
    for chromosome,chromosome_table in all_sample_cn.groupby('Chromosome'):
        all_cn_breakpoints = np.concatenate([chromosome_table['Segment_Start'],chromosome_table['Segment_End']])
        all_cn_breakpoints = np.sort(np.unique(all_cn_breakpoints))
        segment_starts = all_cn_breakpoints[:-1]
        segment_ends = all_cn_breakpoints[1:]
        combined_segment_data['Chromosome'].extend([chromosome]*segment_starts.size)
        combined_segment_data['Segment_Start'].extend(segment_starts)
        combined_segment_data['Segment_End'].extend(segment_ends)
    combined_segment_data = pd.DataFrame(combined_segment_data)
    combined_segment_data.loc[:,'Segment_Width'] = combined_segment_data['Segment_End']-combined_segment_data['Segment_Start']
    return combined_segment_data

def load_allele_specific_cn():
    all_sample_cn = []
    picnic_dir = '../input_data/PICNIC_data/'
    for filename in os.listdir(picnic_dir):
        if not filename.endswith('ascatformat.txt'):
            continue
        filepath = os.path.join(picnic_dir,filename)
        
        sample_cn = pd.read_csv(filepath,sep=" ",dtype={'chr':str})
        
        sample_cn = sample_cn.rename(columns={'name':'Anom_Sample_ID','chr':'Chromosome','start':'Segment_Start','end':'Segment_End','minorA':'Minor_CN','majorB':'Major_CN'})
        sample_cn = sample_cn.drop(columns='num.probes')
        sample_cn.loc[:,'Total_CN']= sample_cn['Major_CN']+sample_cn['Minor_CN']
        
        all_sample_cn.append(sample_cn)
    
    #order of Sample IDs taken from 
    #https://github.com/ucbtsl1/lopez_etal_2019_wgd-cancer/blob/master/Figure3_plot.R (22/07/2022)
    defined_ids = ["T1","T4_DC14","T50_TC13","T50_TC16","T50_TC17","T50_TC17b","T50_TC35","T4_DC25","T4_TC13","T4_TC16","T4_TC17","T4_TC17b","T4_TC35","T50_DC14","T50_DC25"]
    anom_ids = ["D_HCT001_BS","D_HCT001_SU_T1-R1","D_HCT001_SU_T10-R1","D_HCT001_SU_T11-R1","D_HCT001_SU_T12-R1","D_HCT001_SU_T13-R1","D_HCT001_SU_T14-R1","D_HCT001_SU_T2-R1","D_HCT001_SU_T3-R1","D_HCT001_SU_T4-R1","D_HCT001_SU_T5-R1","D_HCT001_SU_T6-R1","D_HCT001_SU_T7-R1","D_HCT001_SU_T8-R1","D_HCT001_SU_T9-R1"]
    id_mapping = pd.DataFrame({'Anom_Sample_ID':anom_ids,'Sample_ID':defined_ids})
    combined_cn =pd.concat(all_sample_cn).reset_index(drop=True)

    combined_cn = combined_cn.merge(id_mapping,how='inner')
    
    combined_cn = combined_cn.drop(columns=['Anom_Sample_ID'])
    combined_cn.loc[:,'Segment_Width'] = combined_cn['Segment_End']-combined_cn['Segment_Start']
    combined_cn = combined_cn[['Sample_ID','Chromosome', 'Segment_Start', 'Segment_End','Segment_Width', 'Major_CN','Minor_CN','Total_CN']]
    combined_cn.loc[:,'Tetraploid'] = combined_cn['Sample_ID'].str.contains('TC')
    return combined_cn

def load_allele_specific_cn():
    all_sample_cn = []
    picnic_dir = '../input_data/PICNIC_data/'
    for filename in os.listdir(picnic_dir):
        if not filename.endswith('ascatformat.txt'):
            continue
        filepath = os.path.join(picnic_dir,filename)
        
        sample_cn = pd.read_csv(filepath,sep=" ",dtype={'chr':str})
        
        sample_cn = sample_cn.rename(columns={'name':'Anom_Sample_ID','chr':'Chromosome','start':'Segment_Start','end':'Segment_End','minorA':'Minor_CN','majorB':'Major_CN'})
        sample_cn = sample_cn.drop(columns='num.probes')
        sample_cn.loc[:,'Total_CN']= sample_cn['Major_CN']+sample_cn['Minor_CN']
        
        all_sample_cn.append(sample_cn)
    
    #order of Sample IDs taken from 
    #https://github.com/ucbtsl1/lopez_etal_2019_wgd-cancer/blob/master/Figure3_plot.R (22/07/2022)
    defined_ids = ["T1","T4_DC14","T50_TC13","T50_TC16","T50_TC17","T50_TC17b","T50_TC35","T4_DC25","T4_TC13","T4_TC16","T4_TC17","T4_TC17b","T4_TC35","T50_DC14","T50_DC25"]
    anom_ids = ["D_HCT001_BS","D_HCT001_SU_T1-R1","D_HCT001_SU_T10-R1","D_HCT001_SU_T11-R1","D_HCT001_SU_T12-R1","D_HCT001_SU_T13-R1","D_HCT001_SU_T14-R1","D_HCT001_SU_T2-R1","D_HCT001_SU_T3-R1","D_HCT001_SU_T4-R1","D_HCT001_SU_T5-R1","D_HCT001_SU_T6-R1","D_HCT001_SU_T7-R1","D_HCT001_SU_T8-R1","D_HCT001_SU_T9-R1"]
    id_mapping = pd.DataFrame({'Anom_Sample_ID':anom_ids,'Sample_ID':defined_ids})
    combined_cn =pd.concat(all_sample_cn).reset_index(drop=True)

    combined_cn = combined_cn.merge(id_mapping,how='inner')
    
    combined_cn = combined_cn.drop(columns=['Anom_Sample_ID'])
    combined_cn.loc[:,'Segment_Width'] = combined_cn['Segment_End']-combined_cn['Segment_Start']
    combined_cn = combined_cn[['Sample_ID','Chromosome', 'Segment_Start', 'Segment_End','Segment_Width', 'Major_CN','Minor_CN','Total_CN']]
    combined_cn.loc[:,'Tetraploid'] = combined_cn['Sample_ID'].str.contains('TC')
    return combined_cn
def load_mutation_data():
    all_sample_mutations = []
    picnic_dir = '../input_data/PICNIC_data/'
    for filename in os.listdir(picnic_dir):
        if not filename.endswith('_completetable.txt'):
            continue
        filepath = os.path.join(picnic_dir,filename)
        sample_id = filename.replace("-",".").replace('_completetable.txt','')
        sample_mutations = pd.read_csv(filepath,sep=" ",dtype={'chr':str})
        sample_mutations = sample_mutations[sample_mutations['is_SNV']]
        
        sample_mutations = sample_mutations.rename(columns={'chr':'Chromosome','start':'Position',f'{sample_id}.ref_count':'Tumor_Ref_Count',f'{sample_id}.var_count':'Tumor_Alt_Count'})
        sample_mutations = sample_mutations[['Chromosome','Position','Tumor_Ref_Count','Tumor_Alt_Count']]
        sample_mutations.loc[:,'Chromosome'] = sample_mutations['Chromosome'].str.replace('chr','')
        sample_mutations.loc[:,'Anom_Sample_ID'] =sample_id
        all_sample_mutations.append(sample_mutations)
    
    #order of Sample IDs taken from 
    #https://github.com/ucbtsl1/lopez_etal_2019_wgd-cancer/blob/master/Figure3_plot.R (22/07/2022)
    defined_ids = ["T1","T4_DC14","T50_TC13","T50_TC16","T50_TC17","T50_TC17b","T50_TC35","T4_DC25","T4_TC13","T4_TC16","T4_TC17","T4_TC17b","T4_TC35","T50_DC14","T50_DC25"]
    anom_ids = ["D_HCT001_BS","D_HCT001_SU_T1-R1","D_HCT001_SU_T10-R1","D_HCT001_SU_T11-R1","D_HCT001_SU_T12-R1","D_HCT001_SU_T13-R1","D_HCT001_SU_T14-R1","D_HCT001_SU_T2-R1","D_HCT001_SU_T3-R1","D_HCT001_SU_T4-R1","D_HCT001_SU_T5-R1","D_HCT001_SU_T6-R1","D_HCT001_SU_T7-R1","D_HCT001_SU_T8-R1","D_HCT001_SU_T9-R1"]
    anom_ids = [anom_id.replace('-','.') for anom_id in anom_ids]
    id_mapping = pd.DataFrame({'Anom_Sample_ID':anom_ids,'Sample_ID':defined_ids})
    combined_mutations =pd.concat(all_sample_mutations).reset_index(drop=True)
    print(combined_mutations)
    combined_mutations = combined_mutations.merge(id_mapping,how='inner')
    
    combined_mutations = combined_mutations.drop(columns=['Anom_Sample_ID'])
    combined_mutations = combined_mutations[['Sample_ID','Chromosome','Position','Tumor_Ref_Count','Tumor_Alt_Count']].dropna()
    
    combined_mutations = combined_mutations[(combined_mutations['Tumor_Alt_Count']+combined_mutations['Tumor_Ref_Count']>=20)]
    #print(combined_mutations[combined_mutations['Tumor_Alt_Count']>=3])
    
    combined_mutations = combined_mutations[(combined_mutations['Tumor_Alt_Count']>=3)]
    print(combined_mutations)
    combined_mutations.loc[:,'Tetraploid'] = combined_mutations['Sample_ID'].str.contains('TC')
    return combined_mutations
def assess_loh(sample_cn,combined_segments):
    bad_count = 0
    total_count = 0
    for seg_index,segment_row in combined_segments.iterrows():
        overlapping_rows = sample_cn[(sample_cn['Chromosome']==segment_row['Chromosome'])&(sample_cn['Segment_Start']<=segment_row['Segment_Start'])&(sample_cn['Segment_End']>=segment_row['Segment_End'])]
        assert len(overlapping_rows.index) == len(overlapping_rows['Sample_ID'].unique())
        n_tetraploid = overlapping_rows['Tetraploid'].sum()
        n_diploid = len(overlapping_rows.index) -n_tetraploid

        tetraploid_loh = (overlapping_rows[overlapping_rows['Tetraploid']]['Minor_CN']==0).sum()
        diploid_loh = (overlapping_rows[~overlapping_rows['Tetraploid']]['Minor_CN']==0).sum()
        if tetraploid_loh>0:
            total_count +=segment_row['Segment_Width']
            if diploid_loh ==0:
                bad_count +=segment_row['Segment_Width']
        


def assess_pre_gain(sample_cn,combined_segments):
    bad_count = 0
    total_count = 0
    for seg_index,segment_row in combined_segments.iterrows():
        overlapping_rows = sample_cn[(sample_cn['Chromosome']==segment_row['Chromosome'])&(sample_cn['Segment_Start']<=segment_row['Segment_Start'])&(sample_cn['Segment_End']>=segment_row['Segment_End'])]
        assert len(overlapping_rows.index) == len(overlapping_rows['Sample_ID'].unique())
        n_tetraploid = overlapping_rows['Tetraploid'].sum()
        n_diploid = len(overlapping_rows.index) -n_tetraploid

        tetraploid_pre_gain = (overlapping_rows[overlapping_rows['Tetraploid']]['Major_CN']==4).sum()
        diploid_gain = (overlapping_rows[~overlapping_rows['Tetraploid']]['Major_CN']>1).sum()
        if tetraploid_pre_gain>0:
            total_count +=segment_row['Segment_Width']
            if diploid_gain ==0:
                bad_count +=segment_row['Segment_Width']

def track_cn_changes(sample_cn,combined_segments):
    passage_pairs = [("T4_DC14","T50_DC14"),("T4_TC13","T50_TC13"),("T4_TC16","T50_TC16"),("T4_TC17","T50_TC17"),("T4_TC17b","T50_TC17b"),("T4_TC35","T50_TC35"),("T4_DC25","T50_DC25")]

    prev_cn_table_data = {'Clone_Name':[],'Passage_50_Major':[],'Passage_50_Minor':[],'Passage_4_Major':[],'Passage_4_Minor':[],'Width':[]}
    for passage_pair in passage_pairs:

        if not 'TC' in passage_pair[0]:
                    continue
        clone_number = passage_pair[1].replace("T50_TC","")
        clone_name = f'Clone {clone_number}'

        prev_cn_store = {}
        for seg_index,segment_row in combined_segments.iterrows():
            overlapping_rows = sample_cn[(sample_cn['Chromosome']==segment_row['Chromosome'])&(sample_cn['Segment_Start']<=segment_row['Segment_Start'])&(sample_cn['Segment_End']>=segment_row['Segment_End'])]
            
            try:
                passage_4 = overlapping_rows[(overlapping_rows['Sample_ID']==passage_pair[0])].iloc[0]
                passage_50 = overlapping_rows[(overlapping_rows['Sample_ID']==passage_pair[1])].iloc[0]
            except IndexError:
                continue
            passage_50_cn = passage_50['Major_CN'],passage_50['Minor_CN']
            passage_4_cn = passage_4['Major_CN'],passage_4['Minor_CN']
            if passage_50_cn not in prev_cn_store.keys():
                prev_cn_store[passage_50_cn] = {}
            if passage_4_cn not in prev_cn_store[passage_50_cn].keys():
                prev_cn_store[passage_50_cn][passage_4_cn] = 0
            prev_cn_store[passage_50_cn][passage_4_cn]+= segment_row['Segment_Width']
            
        for passage_50_cn in prev_cn_store.keys():
            for passage_4_cn,total_width in prev_cn_store[passage_50_cn].items():
                prev_cn_table_data['Clone_Name'].append(clone_name)
                prev_cn_table_data['Passage_50_Major'].append(passage_50_cn[0])
                prev_cn_table_data['Passage_50_Minor'].append(passage_50_cn[1])
                prev_cn_table_data['Passage_4_Major'].append(passage_4_cn[0])
                prev_cn_table_data['Passage_4_Minor'].append(passage_4_cn[1])
                prev_cn_table_data['Width'].append(total_width)
    return pd.DataFrame(prev_cn_table_data)
 
if __name__=='main':            
    sample_cn = load_allele_specific_cn()
    combined_segments = get_combined_segments(sample_cn)

    cn_changes = track_cn_changes(sample_cn,combined_segments)
    cn_changes.to_csv('../in_data/tetraploid_passages_cn_changes.tsv',sep="\t",index=False)
