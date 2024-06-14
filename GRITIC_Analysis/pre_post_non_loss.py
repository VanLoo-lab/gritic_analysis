import pandas as pd
import DataTools
import numpy as np
cn_table = pd.read_csv('../output/cn_table_status.tsv', sep='\t',dtype={'Chromosome':str}).drop(columns=['WGD_Status'])
good_samples =DataTools.get_good_samples()
wgd_status_calls = DataTools.get_wgd_calling_info()
wgd_status_calls = wgd_status_calls[wgd_status_calls['Sample_ID'].isin(good_samples)]
wgd_status_calls = wgd_status_calls[~((wgd_status_calls['Major_CN_Mode']==2) & (~wgd_status_calls['WGD_Status']))]
wgd_status_calls = wgd_status_calls[['Sample_ID','WGD_Status']]
cn_table = cn_table.merge(wgd_status_calls,how='inner')
autosomes = list(map(str,range(1,23)))
cn_table = cn_table[cn_table['Chromosome'].isin(autosomes)]


cn_table = cn_table[['Sample_ID','Chromosome','Segment_Start','Segment_End','Major_CN','Minor_CN','WGD_Status']].drop_duplicates()
cn_table['Segment_Start'] = cn_table['Segment_Start'].astype(int)
cn_table['Segment_End'] = cn_table['Segment_End'].astype(int)
cn_table['Segment_ID'] = cn_table['Chromosome'] + '-' + cn_table['Segment_Start'].astype(str) + '-' + cn_table['Segment_End'].astype(str)

cn_table_wgd = cn_table[cn_table['WGD_Status']].reset_index(drop=True).copy()
cn_table_no_wgd = cn_table[~cn_table['WGD_Status']].reset_index(drop=True).copy()

cn_table_wgd['Pre_Probability'] = cn_table_wgd['Minor_CN'] == 0
cn_table_wgd['Post_Probability'] = cn_table_wgd['Minor_CN'] == 1
cn_table_wgd['Non_Probability'] = 0

cn_table_no_wgd['Pre_Probability'] = 0
cn_table_no_wgd['Post_Probability'] = 0
cn_table_no_wgd['Non_Probability'] = cn_table_no_wgd['Minor_CN'] == 0

cn_table = pd.concat([cn_table_wgd,cn_table_no_wgd])
cn_table['N_Mutations'] = np.nan
out_table = cn_table[['Sample_ID','Segment_ID','Major_CN','Minor_CN','N_Mutations','Pre_Probability','Post_Probability','Non_Probability']]
out_table.to_csv('../output/pre_post_non/pre_post_non_loss.tsv',sep='\t',index=False)