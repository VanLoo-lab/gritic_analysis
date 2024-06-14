import pandas as pd
import numpy as np

import matplotlib
#pdf font type 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Sans'

import matplotlib.pyplot as plt
#font size 18
plt.rcParams.update({'font.size': 18})
import DataTools


def get_diff(sample_gains):
    true_wgd_frac = np.mean(sample_gains[sample_gains['WGD_Status']]['Post_Gain'] > 0)
    true_no_wgd_frac = np.mean(sample_gains[~sample_gains['WGD_Status']]['Post_Gain'] > 0)
    return np.abs(true_wgd_frac - true_no_wgd_frac)
def get_p_value(sample_gains,n_permutations=10000):
    true_diff = get_diff(sample_gains)
    sample_gains_permute = sample_gains.copy()
    diff_store = []
    for i in range(n_permutations):
        sample_gains_permute['WGD_Status'] = RNG.choice(sample_gains_permute['WGD_Status'],size=len(sample_gains_permute),replace=False)
        diff_store.append(get_diff(sample_gains_permute))
    return np.sum(np.array(diff_store) > true_diff)/n_permutations
def get_p_value_star(p_value):
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''


RNG = np.random.default_rng(42)
if __name__ == '__main__':
    wgd_calling_info = DataTools.get_wgd_calling_info()
    bad_samples = wgd_calling_info[(~wgd_calling_info['WGD_Status']) & (wgd_calling_info['Major_CN_Mode']==2)]['Sample_ID'].values
    permute = True
    for apply_penalty in [True,False]:
        min_mutations =0
        threshold = 0

        pre_post_close_path = f'../output/pre_post_close/{apply_penalty}_with_permutations_{permute}_pre_post.tsv'
    
        pre_post_close = pd.read_csv(pre_post_close_path,sep='\t')



        pre_post_close.loc[:,'Permuted_ID'] = pre_post_close['Sample_ID']
        pre_post_close.loc[:,'Sample_ID'] = pre_post_close['Sample_ID'].apply(lambda x: x.split('_permuted')[0])
        pre_post_close = pre_post_close[~pre_post_close['Sample_ID'].isin(bad_samples)]


        route_table = DataTools.load_route_table(apply_penalty=False,wgd_status=False)
        seg_info = route_table[['Sample_ID','Segment_ID','Chromosome','Segment_Start','Segment_End','Segment_Width','Cancer_Type','WGD_Status']].drop_duplicates()


        pre_post_close = pre_post_close.merge(seg_info,on=['Sample_ID','Segment_ID'],how='inner')
        pre_post_close.loc[:,'Pre_Gain'] = (pre_post_close['Pre_Probability'] >= 0.5).astype(int)
        pre_post_close.loc[:,'Post_Gain'] = (pre_post_close['Post_Probability'] >= 0.5).astype(int)

        sample_gains_pre = pre_post_close.groupby(['Cancer_Type','WGD_Status','Permuted_ID','WGD_Timing_Median']).apply(lambda x: np.sum(x['Pre_Gain']*x['Segment_Width'])).reset_index().rename(columns={0:'Pre_Gain'})
        sample_gains_post = pre_post_close.groupby(['Permuted_ID']).apply(lambda x: np.sum(x['Post_Gain']*x['Segment_Width'])).reset_index().rename(columns={0:'Post_Gain'})

        #merge
        sample_gains = sample_gains_pre.merge(sample_gains_post,on='Permuted_ID',how='inner')

        sample_gain_wgd_summary = sample_gains[['Permuted_ID']]
        sample_gains = sample_gains.merge(sample_gain_wgd_summary,on='Permuted_ID',how='inner')

        sample_gains_single_permutation = sample_gains[(sample_gains['Permuted_ID'].str.contains('permuted_0')) | (~(sample_gains['Permuted_ID'].str.contains('permuted')))]
        sample_gains_single_permutation_wgd = sample_gains_single_permutation[sample_gains_single_permutation['WGD_Status']]
        sample_gains_single_permutation_no_wgd = sample_gains_single_permutation[~sample_gains_single_permutation['WGD_Status']]

        true_wgd_frac = np.mean(sample_gains_single_permutation_wgd['Post_Gain'] > 0)
        true_no_wgd_frac = np.mean(sample_gains_single_permutation_no_wgd['Post_Gain'] > 0)

        bootstrap_data = {'WGD':[],'No_WGD':[]}
        
        n_bootstraps = 250
        for i in range(n_bootstraps):
            sample_gains_single_permutation_wgd_bootstrap = sample_gains_single_permutation_wgd.sample(frac=1,replace=True)
            sample_gains_single_permutation_no_wgd_bootstrap = sample_gains_single_permutation_no_wgd.sample(frac=1,replace=True)
            bootstrap_data['WGD'].append(np.mean(sample_gains_single_permutation_wgd_bootstrap['Post_Gain'] > 0))
            bootstrap_data['No_WGD'].append(np.mean(sample_gains_single_permutation_no_wgd_bootstrap['Post_Gain'] > 0))
        p_value = get_p_value(sample_gains_single_permutation)
        p_value_star = get_p_value_star(p_value)

        wgd_ci = np.percentile(bootstrap_data['WGD'],[5/2,100-5/2])
        no_wgd_ci = np.percentile(bootstrap_data['No_WGD'],[5/2,100-5/2])
        fig,ax = plt.subplots(figsize=(4,10))
        ax.bar(['WGD','Non-WGD Control'],[true_wgd_frac,true_no_wgd_frac],yerr=[[true_wgd_frac-wgd_ci[0],true_no_wgd_frac-no_wgd_ci[0]],[wgd_ci[1]-true_wgd_frac,no_wgd_ci[1]-true_no_wgd_frac]],capsize=10,color=['#D6102C','#808080'])
        ax.set_ylabel('Fraction of Samples with Gains Post-WGD')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #annotate p-value stard
        ax.annotate(p_value_star,xy=(0.5,0.5),xytext=(0.5,0.95),xycoords='axes fraction',textcoords='axes fraction',ha='center',va='center',fontsize=20)
        plt.savefig(f'../plots/fraction_gains_post_wgd/fraction_gains_post_wgd_prior_{apply_penalty}.pdf',bbox_inches='tight')

