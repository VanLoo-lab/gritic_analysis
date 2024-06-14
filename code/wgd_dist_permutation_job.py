import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#font type 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = 'Nimbus Sans'
plt.rcParams['ps.fonttype'] = 42
#font size 18
plt.rcParams['font.size'] = 18

from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu

import DataTools

def get_sample_distributions(posterior_table, column_name):
    sample_distributions = {}
    for sample_id,sample_table in posterior_table.groupby('Sample_ID'):
        if column_name =='WGD_Timing':
            sample_table = sample_table.drop_duplicates(subset=['Segment_ID','WGD_Timing'],keep='first')
        sample_distributions[sample_id] = sample_table[column_name].to_numpy()
    return sample_distributions

def get_wgd_distance_distribution(wgd_distributions, gain_distributions, permute=False,runs_per_sample=25):
    wgd_ids = sorted(wgd_distributions.keys())
    gain_ids = sorted(gain_distributions.keys())
    if permute:
        RNG.shuffle(gain_ids)
    distance_distribution = []
    for index in range(len(wgd_ids)):
        wgd_id = wgd_ids[index]
        gain_id = gain_ids[index]
        total_distance = 0
        for i in range(runs_per_sample):
            wgd_dist = RNG.choice(wgd_distributions[wgd_id],size=len(gain_distributions[gain_id]))
            distance = np.median(np.abs(wgd_dist-gain_distributions[gain_id]))
            total_distance += distance
        distance_distribution.append(total_distance/runs_per_sample)
    return np.array(distance_distribution)
def format_p_value(p):
    if p<0.01:
        #return f'{p:.1e}'
        #format as x10^-n using superscript
        print(p)
        p_str = f'{p:.1e}'
        p_str = f'{p_str.split("e")[0]}$\\times$10$^{{{int(p_str.split("e")[1])}}}$'
        print(p_str)
        print('---------')
        return p_str

    else:
        #2sf
        return f'{p:.2f}'



RNG = np.random.default_rng(42)

for apply_penalty in [True,False]:
    posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty)
    posterior_table_wgd = posterior_table[posterior_table['WGD_Status']]

    distribution_store = {}
    for cancer_type,cancer_type_posterior_table in posterior_table_wgd.groupby('Cancer_Type'):
        if len(cancer_type_posterior_table['Sample_ID'].unique())<=50:
            continue
        wgd_distributions = get_sample_distributions(cancer_type_posterior_table, 'WGD_Timing')
        gain_distributions = get_sample_distributions(cancer_type_posterior_table, 'Gain_Timing')

        true_distance_dist = get_wgd_distance_distribution(wgd_distributions, gain_distributions,permute=False)
        permuted_distance_dist = get_wgd_distance_distribution(wgd_distributions, gain_distributions,permute=True)

        distribution_store[cancer_type] = {'Permuted':permuted_distance_dist,'True':true_distance_dist}


    fig,axs = plt.subplots(3,5,figsize=(30,15))
    axs = axs.flatten()
    ax_index = 0
    for cancer_type,data in distribution_store.items():
        axs[ax_index].hist(data['Permuted'],bins=20,color='#59baff',alpha=0.5,density=True,label='WGD Permuted')
        axs[ax_index].hist(data['True'],bins=20,color='#e80f00',alpha=0.2,density=True,label='Truth')
        axs[ax_index].set_title(cancer_type)
        axs[ax_index].set_xlabel('Median Time To WGD')
        axs[ax_index].set_ylabel('Density')

        median_permuted = np.median(data['Permuted'])
        median_true = np.median(data['True'])
        axs[ax_index].axvline(median_permuted,color='#59baff',linestyle='--',lw=4)
        axs[ax_index].axvline(median_true,color='#e80f00',linestyle='--',lw=4)
        u,p = mannwhitneyu(data['Permuted'],data['True'])
        axs[ax_index].text(0.1,0.9,f'p={format_p_value(p)}',transform=axs[ax_index].transAxes)
        if ax_index==0:
            legend_elements = [Line2D([0], [0], color='#59baff', lw=5, label='WGD Permuted'),Line2D([0], [0], color='#e80f00', lw=5, label='Truth'),Line2D([0], [0], color='black', lw=3, linestyle='--', label='Median')]
            axs[ax_index].legend(handles=legend_elements, loc='upper right')
        ax_index += 1


    plt.tight_layout()
    plt.savefig(f'../plots/wgd_permutation_test/wgd_permutation_cancer_types_separate_penalty_{apply_penalty}.pdf',bbox_inches='tight')


    combined_permuted = np.concatenate([data['Permuted'] for data in distribution_store.values()])
    combined_true = np.concatenate([data['True'] for data in distribution_store.values()])
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax.hist(combined_permuted,bins=20,color='#59baff',alpha=0.5,density=True)
    ax.hist(combined_true,bins=20,color='#e80f00',alpha=0.2,density=True)
    median_permuted = np.median(combined_permuted)
    median_true = np.median(combined_true)
    ax.axvline(median_permuted,color='#59baff',linestyle='--',lw=4)
    ax.axvline(median_true,color='#e80f00',linestyle='--',lw=4)
    ax.set_xlabel('Median Time To WGD')
    ax.set_ylabel('Density')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    u,p = mannwhitneyu(combined_permuted,combined_true)
    ax.text(0.1,0.9,f'p={format_p_value(p)}',transform=ax.transAxes)
    legend_elements = [Line2D([0], [0], color='#59baff', lw=5, label='WGD Permuted'),Line2D([0], [0], color='#e80f00', lw=5, label='Truth'),Line2D([0], [0], color='black', lw=3, linestyle='--', label='Median')]
            

    ax.legend(handles=legend_elements, loc='upper right')
    plt.savefig(f'../plots/wgd_permutation_test/wgd_permutation_all_cancer_types_penalty_{apply_penalty}.pdf',bbox_inches='tight')