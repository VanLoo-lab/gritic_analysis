import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
#font size 18
plt.rcParams.update({'font.size': 18})
#font type pdf 42
plt.rcParams['pdf.fonttype'] = 42
#nimbus sans font
plt.rcParams['font.sans-serif'] = "Nimbus Sans"

def load_metadata_table(input_dir,dataset):
    metadata_table = pd.read_csv(input_dir /f'complete_metadata_table_{dataset}.tsv',sep="\t").drop_duplicates()

    metadata_table['Segment_ID'] = metadata_table['Segment_ID'].str.replace('Segment_','')

    metadata_table['CN_State'] = metadata_table['Major_CN'].astype(str) + '_' + metadata_table['Minor_CN'].astype(str)
    metadata_table['N_Gains'] = np.maximum(0,metadata_table['Major_CN']-1) + np.maximum(0,metadata_table['Minor_CN']-1)
    #extract real gain timing from metadata table and convert to array
    metadata_table['Real_Gain_Timing'] = metadata_table['Gain_Timing'].apply(lambda x: np.array(list(eval(x).values())))
    metadata_table['Real_Gain_Timing'] = metadata_table.apply(lambda x: np.append(x['Real_Gain_Timing'],np.repeat(x['WGD_Timing'],x['N_Gains']-x['Real_Gain_Timing'].size)),axis=1)
    metadata_table['Real_Gain_Timing'] = metadata_table['Real_Gain_Timing'].apply(lambda x: np.sort(x))
    metadata_table = metadata_table[['Sample_ID','CN_State','Major_CN','Minor_CN','Segment_ID','Real_Gain_Timing','WGD_Timing','N_Mutations']].rename(columns={'Sample_ID':'Base_ID','WGD_Timing':'Real_WGD_Timing'})
    return metadata_table

def load_posterior_table(input_dir,dataset,apply_penalty):
    posterior_table = pd.read_csv(input_dir /f'complete_posterior_table_penalty_{apply_penalty}_{dataset}.tsv',sep="\t")
    posterior_table = posterior_table[(posterior_table['Major_CN']>2) & (posterior_table['N_Mutations']>=20)]
    posterior_table['WGD_Timing'] = posterior_table['WGD_Timing'].fillna(-1)
    

    posterior_table['WGD_Constraint'] = posterior_table['Sample_ID'].str.split('_').str[2]
    posterior_table = posterior_table[posterior_table['WGD_Constraint'].isin(['True','False'])]
    posterior_table['WGD_Constraint'] = posterior_table['WGD_Constraint'].apply(lambda x: x=='True')
    posterior_table['Base_ID'] = posterior_table['Sample_ID'].str.split('_').str[0]
    posterior_table['N_Gains'] = np.maximum(0,posterior_table['Major_CN']-1) + np.maximum(0,posterior_table['Minor_CN']-1)


    posterior_table_measured_timing = posterior_table.groupby(['Base_ID','Sample_ID','Segment_ID','Posterior_Sample_Index','WGD_Constraint','WGD_Timing']).apply(lambda x: np.append(x['Gain_Timing'].values,np.repeat(x['WGD_Timing'].values[0],x['N_Gains'].values[0]-x['Gain_Timing'].values.size)))

    posterior_table_measured_timing = posterior_table_measured_timing.apply(lambda x: np.sort(x))


    #flatten
    posterior_table_measured_timing = posterior_table_measured_timing.reset_index().rename(columns={0:'Measured_Gain_Timing'})
   
    return posterior_table_measured_timing

dataset = 'state_validation_wgd_constraint_test'
input_dir = Path(f'../../output/{dataset}/')
apply_penalty = False
posterior_table = load_posterior_table(input_dir,dataset,apply_penalty)
metadata_table = load_metadata_table(input_dir,dataset)

posterior_table = posterior_table.merge(metadata_table,how='inner')

#min_real_gain_timing
assert posterior_table['Real_Gain_Timing'].apply(lambda x: x.min()).min()>=-1e-9

posterior_table_average_error = posterior_table.groupby(['Base_ID','Sample_ID','WGD_Constraint','Segment_ID','Posterior_Sample_Index','WGD_Timing']).apply(lambda x: np.mean(np.mean(np.abs(x['Real_Gain_Timing'].values-x['Measured_Gain_Timing'].values)))).groupby(['Base_ID','WGD_Constraint','Segment_ID']).mean().reset_index().rename(columns={0:'Average_Error'})

posterior_table_no_wgd = posterior_table.copy()
posterior_table_no_wgd['WGD_Gain'] = posterior_table_no_wgd.apply(lambda x: np.isclose(x['Real_WGD_Timing'],x['Real_Gain_Timing']),axis=1)

posterior_table_no_wgd=  posterior_table_no_wgd[posterior_table_no_wgd['WGD_Gain'].apply(lambda x: x.mean())<0.999]


posterior_table_no_wgd['Real_Gain_Timing'] = posterior_table_no_wgd.apply(lambda x: x['Real_Gain_Timing'][~x['WGD_Gain']],axis=1)
posterior_table_no_wgd['Measured_Gain_Timing'] = posterior_table_no_wgd.apply(lambda x: x['Measured_Gain_Timing'][~x['WGD_Gain']],axis=1)
posterior_table_average_error_no_wgd = posterior_table_no_wgd.groupby(['Base_ID','WGD_Constraint','Segment_ID','Posterior_Sample_Index']).apply(lambda x: np.mean(np.mean(np.abs(x['Real_Gain_Timing'].values-x['Measured_Gain_Timing'].values)))).groupby(['Base_ID','WGD_Constraint','Segment_ID']).mean().reset_index().rename(columns={0:'Average_Error'})

#pivot on wgd constraint
posterior_table_average_error = posterior_table_average_error.pivot(index=['Base_ID','Segment_ID'],columns='WGD_Constraint',values='Average_Error').reset_index().dropna()
posterior_table_average_error_no_wgd = posterior_table_average_error_no_wgd.pivot(index=['Base_ID','Segment_ID'],columns='WGD_Constraint',values='Average_Error').reset_index().dropna()

posterior_table_average_error['Error_Diff'] = posterior_table_average_error[False]-posterior_table_average_error[True]
posterior_table_average_error_no_wgd['Error_Diff'] = posterior_table_average_error_no_wgd[False]-posterior_table_average_error_no_wgd[True]
posterior_table_average_error.to_csv('../output/wgd_constraint_analysis/posterior_table_average_error.tsv',sep='\t',index=False)
posterior_table_average_error_no_wgd.to_csv('../output/wgd_constraint_analysis/posterior_table_average_error_no_wgd.tsv',sep='\t',index=False)

posterior_table_average_error = pd.read_csv('../output/wgd_constraint_analysis/posterior_table_average_error.tsv',sep='\t')
posterior_table_average_error_no_wgd = pd.read_csv('../output/wgd_constraint_analysis/posterior_table_average_error_no_wgd.tsv',sep='\t')
posterior_table_average_error['Error_Diff'] = posterior_table_average_error['False']-posterior_table_average_error['True']
posterior_table_average_error_no_wgd['Error_Diff'] = posterior_table_average_error_no_wgd['False']-posterior_table_average_error_no_wgd['True']

reduction_prop = np.mean(posterior_table_average_error['Error_Diff']>0)*100
reduction_prop_no_wgd = np.mean(posterior_table_average_error_no_wgd['Error_Diff']>0)*100

bins = np.linspace(-0.1,0.2,20)
fig,axs = plt.subplots(1,2,figsize=(20,10))
axs[0].hist(posterior_table_average_error['Error_Diff'],bins=bins,color='#11219c',density=True)
axs[0].set_xlabel('Error Without WGD Constraint - Error With WGD Constraint')
axs[0].set_ylabel('Density')
axs[0].set_title(f'All Gains \n Reduced error in {reduction_prop:.1f}% of simulated segments')
axs[1].hist(posterior_table_average_error_no_wgd['Error_Diff'],bins=bins,color='#11219c',density=True)
axs[1].set_xlabel('Error Without WGD Constraint - Error With WGD Constraint')
axs[1].set_ylabel('Density')
axs[1].set_title(f'Non-WGD Gains \n Reduced Error in {reduction_prop_no_wgd:.1f}% of simulated segments')
plt.tight_layout()
plt.savefig('../plots/wgd_constraint_analysis/wgd_constraint_histograms.pdf',bbox_inches='tight')
