import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Sans'
import matplotlib.pyplot as plt
#font size 18
matplotlib.rcParams.update({'font.size': 18})
import DataTools
import sys
import os

from scipy.stats import gaussian_kde

RNG = np.random.default_rng(int(sys.argv[1])*42)
def load_params(run_number):
    params = [(True,True),(True,False),(False,True),(False,False),(True,None),(False,None)]
    apply_penalty,wgd_status = params[run_number]
    return apply_penalty,wgd_status

def get_posterior(apply_penalty,wgd_status):
    posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty)

    if wgd_status is not None:
        if wgd_status:
            posterior_table = posterior_table[~posterior_table['WGD_Timing'].isnull()]
        else:
            posterior_table = posterior_table[posterior_table['WGD_Timing'].isnull()]
    return posterior_table


def get_kde_evaluations(posterior_table,points,major_cn_range):
    full_kde_store = {}
    min_kde_store = {}
    for major_cn in major_cn_range:
        full_kde = gaussian_kde(posterior_table[posterior_table['Major_CN']==major_cn]['Gain_Rank'],weights=posterior_table[posterior_table['Major_CN']==major_cn]['Segment_Width'])
        min_kde = gaussian_kde(posterior_table_min[posterior_table_min['Major_CN']==major_cn]['Gain_Rank'],weights=posterior_table_min[posterior_table_min['Major_CN']==major_cn]['Segment_Width'])
        full_kde_store[major_cn] = full_kde.evaluate(points)
        min_kde_store[major_cn] = min_kde.evaluate(points)
    return full_kde_store,min_kde_store

def get_bootstrapped_posteriors(posterior_table_sample_dict,posterior_table_min_sample_dict):
    sample_ids = list(posterior_table_sample_dict.keys())
    assert len(sample_ids) == len(posterior_table_min_sample_dict.keys())
    bootstrapped_sample_ids = RNG.choice(sample_ids,size=len(sample_ids),replace=True)
    posterior_table_bootstrapped = []
    posterior_table_min_bootstrapped = []
    for i,sample_id in enumerate(bootstrapped_sample_ids):
        print('doing bootstrap',i,'of',len(bootstrapped_sample_ids))
        posterior_sample_bootstrapped = posterior_table_sample_dict[sample_id].copy()
        posterior_sample_bootstrapped = posterior_sample_bootstrapped.rename(columns={'Sample_ID':f'Sample_ID_{i}'})
        posterior_sample_min_bootstrapped = posterior_table_min_sample_dict[sample_id].copy()
        posterior_sample_min_bootstrapped = posterior_sample_min_bootstrapped.rename(columns={'Sample_ID':f'Sample_ID_{i}'})

        posterior_table_bootstrapped.append(posterior_sample_bootstrapped)
        posterior_table_min_bootstrapped.append(posterior_sample_min_bootstrapped)
    print('concatenating')
    posterior_table_bootstrapped = pd.concat(posterior_table_bootstrapped,axis=0)
    posterior_table_min_bootstrapped = pd.concat(posterior_table_min_bootstrapped,axis=0)
    return posterior_table_bootstrapped,posterior_table_min_bootstrapped

def get_title_text(apply_penalty,wgd_status):
    if apply_penalty:
        title = 'With Non-Parsimony Penalty'
    else:
        title = 'Without Non-Parsimony Penalty'
    
    if wgd_status is not None:
        if wgd_status:
            wgd_text = 'WGD Samples Only'
        else:
            wgd_text = 'Non-WGD Samples Only'
        title = title + '-' + wgd_text
    return title

    

if __name__ == '__main__':
    run_number = int(sys.argv[1])
    apply_penalty,wgd_status = load_params(run_number)
    posterior_table = get_posterior(apply_penalty,wgd_status)
    posterior_table = posterior_table[['Sample_ID','Major_CN','Minor_CN','Segment_Width','Gain_Timing','Segment_ID','Posterior_Sample_Index']]

    wgd_calling_info = DataTools.get_wgd_calling_info()
    bad_samples = wgd_calling_info[(~wgd_calling_info['WGD_Status']) & (wgd_calling_info['Major_CN_Mode']==2)]['Sample_ID'].values
    posterior_table = posterior_table[~posterior_table['Sample_ID'].isin(bad_samples)]

    posterior_table.loc[:,'Gain_Rank'] = posterior_table.groupby(['Sample_ID'])['Gain_Timing'].rank(pct=True)
    #filter for min gain timing per posterior_sample_index
    posterior_table_min = posterior_table[posterior_table['Gain_Timing'] == posterior_table.groupby(['Sample_ID','Segment_ID','Posterior_Sample_Index'])['Gain_Timing'].transform('min')]



    full_kdes = {}
    min_kdes = {}

    for major_cn in [3,4,5,6,7,8]:
        full_kdes[major_cn] = gaussian_kde(posterior_table[posterior_table['Major_CN']==major_cn]['Gain_Rank'],weights=posterior_table[posterior_table['Major_CN']==major_cn]['Segment_Width'])
        min_kdes[major_cn] = gaussian_kde(posterior_table_min[posterior_table_min['Major_CN']==major_cn]['Gain_Rank'],weights=posterior_table_min[posterior_table_min['Major_CN']==major_cn]['Segment_Width'])

    points = np.linspace(0,1,100)

    fig, axs = plt.subplots(2,1,figsize=(14,10))
    for major_cn in [3,4,5,6,7,8]:
        axs[0].plot(points*100, min_kdes[major_cn](points), label=f'Major CN = {major_cn}', linewidth=3)
        axs[1].plot(points*100, full_kdes[major_cn](points), linewidth=3)
        


    axs[0].set_xlabel('Timing of First Gain (Percentile)')
    axs[1].set_xlabel('Timing of All Gains (Percentile)')

    axs[0].set_ylabel('Density')
    axs[1].set_ylabel('Density')
    axs[0].legend()

    fig.suptitle(get_title_text(apply_penalty,wgd_status))
    plt.tight_layout()


    plt.savefig(f'../plots/major_cn_timing/major_cn_timing_apply_penalty_{apply_penalty}_wgd_status_{wgd_status}.png', dpi=300)
    plt.savefig(f'../plots/major_cn_timing/major_cn_timing_apply_penalty_{apply_penalty}_wgd_status_{wgd_status}.pdf', dpi=300)





