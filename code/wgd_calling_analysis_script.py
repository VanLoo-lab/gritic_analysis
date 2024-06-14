import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import  mannwhitneyu
#font type 42
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
#font size 18
plt.rcParams.update({'font.size': 32})
RNG = np.random.default_rng(42)
import DataTools

def get_ploidy(sample_cn_table):
    return np.average(sample_cn_table['Major_CN']+sample_cn_table['Minor_CN'],weights=sample_cn_table['Segment_Width'])

def get_loh(sample_cn_table):
    return np.sum(sample_cn_table[sample_cn_table['Minor_CN']==0]['Segment_Width'])/np.sum(sample_cn_table['Segment_Width'])
def get_ploidy_loh_analysis(cn_table):
    loh_table = {'Sample_ID':[],'Ploidy':[],'LOH':[]}
    for sample_id,sample_table in cn_table.groupby('Sample_ID'):
        ploidy = get_ploidy(sample_table)
        loh = get_loh(sample_table)
        loh_table['Sample_ID'].append(sample_id)
        loh_table['Ploidy'].append(ploidy)
        loh_table['LOH'].append(loh)
    return pd.DataFrame(loh_table)


def classify_wgd(row):
    if pd.isnull(row['GRITIC_WGD_Status']):
        return 'Major CN Mode > 2'
    elif row['Copy_Number_WGD_Status'] and row['GRITIC_WGD_Status']:
        return 'Both WGD'
    elif not row['Copy_Number_WGD_Status'] and not row['GRITIC_WGD_Status']:
        return 'Neither WGD'
    elif row['Copy_Number_WGD_Status'] and not row['GRITIC_WGD_Status']:
        return 'Dataset WGD Only'
    elif not row['Copy_Number_WGD_Status'] and row['GRITIC_WGD_Status']:
        return 'GRITIC WGD Only'
    else:
        return 'NA'


cn_table = pd.read_csv('../output/cn_table_status.tsv', sep='\t')
autosomes = list(map(str,range(1,23)))
cn_table = cn_table[cn_table['Chromosome'].isin(autosomes)]
cn_table.loc[:,'Segment_Width']+=1

ploidy_loh = get_ploidy_loh_analysis(cn_table)

copy_number_wgd_calls = cn_table[['Sample_ID','WGD_Status']].drop_duplicates().rename(columns={'WGD_Status':'Copy_Number_WGD_Status'})
assert len(copy_number_wgd_calls['Sample_ID'].unique())==len(copy_number_wgd_calls)

gritic_wgd_calls = DataTools.get_wgd_calling_info()[['Sample_ID','WGD_Status','Overlap_Proportion','Major_CN_Mode']].drop_duplicates().rename(columns={'WGD_Status':'GRITIC_WGD_Status'})
assert len(gritic_wgd_calls['Sample_ID'].unique())==len(gritic_wgd_calls)

ploidy_loh = ploidy_loh.merge(copy_number_wgd_calls,on='Sample_ID',how='left').merge(gritic_wgd_calls,on='Sample_ID',how='left')
good_samples = DataTools.get_good_samples()
ploidy_loh = ploidy_loh[ploidy_loh['Sample_ID'].isin(good_samples)]

ploidy_loh['GRITIC_WGD_Status'].value_counts()

ploidy_loh[(ploidy_loh['Overlap_Proportion']>0.6) & (ploidy_loh['Overlap_Proportion']<0.8)]

ploidy_loh['Major_CN_Mode'].value_counts()

wgd_calling_sim_table = pd.read_csv('../../output/state_validation_top_states_wgd_calling_test/complete_wgd_calling_info_state_validation_top_states_wgd_calling_test.tsv', sep='\t')
wgd_calling_sim_table = wgd_calling_sim_table[~pd.isnull(wgd_calling_sim_table['Best_Overlap_Timing'])]


sim_wgd_overlap_prop = wgd_calling_sim_table[~wgd_calling_sim_table['Sample_ID'].str.contains('_No_WGD')]['Overlap_Proportion'].values
sim_no_wgd_overlap_prop = wgd_calling_sim_table[wgd_calling_sim_table['Sample_ID'].str.contains('_No_WGD')]['Overlap_Proportion'].values
overlap_proportion = ploidy_loh[ (ploidy_loh['Major_CN_Mode']>1.5)]['Overlap_Proportion'].values
overlap_bins = np.linspace(0,1,21)
fig,axs = plt.subplots(2,1,figsize=(21,21),sharex=True)

axs[0].hist(sim_wgd_overlap_prop,bins=overlap_bins,density=True,alpha=0.5,label='Simulated WGD',color='#3DADB8')
axs[0].hist(sim_no_wgd_overlap_prop,bins=overlap_bins,density=True,alpha=0.5,label='Simulated Non-WGD',color='#C33C54')
axs[1].hist(overlap_proportion,bins=overlap_bins,density=True,alpha=0.5,label='True',color='#6331D8')

axs[0].set_title('Simulations')

axs[1].set_title('True Data')
for i in range(2):
    label = 'WGD Calling Threshold' if i == 0 else ''
    axs[i].axvline(0.6,color='red',linestyle='--',lw=5,label=label)
axs[0].legend()
axs[1].set_xlabel('Proportion of Major Copy Number 2 Segments With Overlapping Timing')
axs[0].set_ylabel('Density')
plt.savefig('../plots/wgd_overlap_analysis/wgd_overlap_vs_sim.pdf',bbox_inches='tight')




ploidy_loh['WGD_Status_Class'] = ploidy_loh.apply(classify_wgd,axis=1)
#add color map for WGD status class
wgd_status_class_color_map = {'Both WGD':'#FFC914','Neither WGD':'#E94F37','Dataset WGD Only':'#B47EAC','GRITIC WGD Only':'#3F88C5','Major CN Mode > 2':'#44BBA4'}
wgd_status_class_color_map['Dataset WGD\nOnly'] = wgd_status_class_color_map['Dataset WGD Only']
wgd_status_class_color_map['GRITIC WGD\nOnly'] = wgd_status_class_color_map['GRITIC WGD Only']
wgd_status_class_color_map['Major CN Mode\n> 2'] = wgd_status_class_color_map['Major CN Mode > 2']



ploidy_loh['WGD_Status_Class_Color'] = ploidy_loh['WGD_Status_Class'].apply(lambda x: wgd_status_class_color_map[x])

fig,ax = plt.subplots(figsize=(15,15))
ax.scatter(ploidy_loh['LOH'],ploidy_loh['Ploidy'],s=20,c=ploidy_loh['WGD_Status_Class_Color'],alpha=0.5)
ax.set_xlabel('LOH')
ax.set_ylabel('Ploidy')
#add custom legend
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label='Both WGD',markerfacecolor=wgd_status_class_color_map['Both WGD'], markersize=15),
                     plt.Line2D([0], [0], marker='o', color='w', label='Neither WGD',markerfacecolor=wgd_status_class_color_map['Neither WGD'], markersize=15),
                        plt.Line2D([0], [0], marker='o', color='w', label='Dataset WGD Only',markerfacecolor=wgd_status_class_color_map['Dataset WGD Only'], markersize=15),
                        plt.Line2D([0], [0], marker='o', color='w', label='Gritic WGD Only',markerfacecolor=wgd_status_class_color_map['GRITIC WGD Only'], markersize=15),
                        plt.Line2D([0], [0], marker='o', color='w', label='Major CN Mode > 2',markerfacecolor=wgd_status_class_color_map['Major CN Mode > 2'], markersize=15)]
ax.legend(handles=legend_elements,loc='upper left',bbox_to_anchor=(1.05, 1))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig('../plots/wgd_overlap_analysis/ploidy_loh_scatter.pdf',bbox_inches='tight')

wgd_class_counts = ploidy_loh['WGD_Status_Class'].value_counts()

wgd_class_counts.index = wgd_class_counts.index.str.replace('Copy Number WGD Only','Copy Number WGD\nOnly').str.replace('GRITIC WGD Only','GRITIC WGD\nOnly').str.replace('Major CN Mode > 2','Major CN Mode\n> 2')
fig,ax = plt.subplots(figsize=(25,15))
ax.bar(wgd_class_counts.index,wgd_class_counts.values,color=[wgd_status_class_color_map[x] for x in wgd_class_counts.index])
ax.set_ylabel('Number of Samples')
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig('../plots/wgd_overlap_analysis/wgd_class_counts.pdf',bbox_inches='tight')


wgd_calling_info = DataTools.get_wgd_calling_info()


classifications = DataTools.load_cancer_type_classifications()
wgd_calling_info = wgd_calling_info.merge(classifications,how='inner')

wgd_calling_info_major_cn_2 = wgd_calling_info[wgd_calling_info['Major_CN_Mode']==2]

wgd_calling_info_major_cn_2 = wgd_calling_info_major_cn_2[~pd.isnull(wgd_calling_info_major_cn_2['Best_Overlap_Timing'])]
wgd_calling_info_major_cn_2_wgd_false = wgd_calling_info_major_cn_2[wgd_calling_info_major_cn_2['WGD_Status']==False]
wgd_calling_info_major_cn_2_wgd_true = wgd_calling_info_major_cn_2[wgd_calling_info_major_cn_2['WGD_Status']==True]

wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts = wgd_calling_info_major_cn_2_wgd_false['Cancer_Type_Code'].value_counts()
wgd_calling_info_major_cn_2_wgd_true_cancer_type_counts = wgd_calling_info_major_cn_2_wgd_true['Cancer_Type_Code'].value_counts()
#set other if counts are less than 10 and sum
wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts['Other'] = np.sum(wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts[wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts<5])
wgd_calling_info_major_cn_2_wgd_true_cancer_type_counts['Other'] = np.sum(wgd_calling_info_major_cn_2_wgd_true_cancer_type_counts[wgd_calling_info_major_cn_2_wgd_true_cancer_type_counts<50])
#drop counts less than 10
wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts = wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts[wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts>=5]
wgd_calling_info_major_cn_2_wgd_true_cancer_type_counts = wgd_calling_info_major_cn_2_wgd_true_cancer_type_counts[wgd_calling_info_major_cn_2_wgd_true_cancer_type_counts>=50]
wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts

fig,axs= plt.subplots(2,1,figsize=(29,15))
axs[0].bar(wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts.index,wgd_calling_info_major_cn_2_wgd_false_cancer_type_counts.values,color='#E94F37')
axs[0].set_title('Major CN Mode 2, GRITIC WGD False')

axs[1].bar(wgd_calling_info_major_cn_2_wgd_true_cancer_type_counts.index,wgd_calling_info_major_cn_2_wgd_true_cancer_type_counts.values,color='#44BBA4')
axs[1].set_title('Major CN Mode 2, GRITIC WGD True')


for i in range(2):
    axs[i].spines['top'].set_visible(False)
    axs[i].spines['right'].set_visible(False)
    axs[i].set_ylabel('Number of Samples')

plt.savefig('../plots/wgd_overlap_analysis/major_cn_mode_2_cancer_type_counts.pdf',bbox_inches='tight')



from scipy.stats import mannwhitneyu
fig,ax = plt.subplots(figsize=(20,15))
ax.hist(wgd_calling_info_major_cn_2_wgd_false['Best_Overlap_Timing'],bins=np.linspace(0,1,11),color='#E94F37',label='GRITIC Insufficient WGD Overlap',density=True,alpha=0.5)
ax.hist(wgd_calling_info_major_cn_2_wgd_true['Best_Overlap_Timing'],bins=np.linspace(0,1,11),color='#37E9E9',label='GRITIC Sufficient WGD Overlap',density=True,alpha=0.5)

u,p = mannwhitneyu(wgd_calling_info_major_cn_2_wgd_false['Best_Overlap_Timing'],wgd_calling_info_major_cn_2_wgd_true['Best_Overlap_Timing'])

ax.legend()
ax.set_xlabel('Timing of Maximum Overlap')
ax.set_ylabel('Density')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#annotate with p-value
ax.text(0.5,0.5,'***',fontsize=30,transform=ax.transAxes)

plt.savefig('../plots/wgd_overlap_analysis/major_cn_mode_2_best_overlap_timing.pdf',bbox_inches='tight')
