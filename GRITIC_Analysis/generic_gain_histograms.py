import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Sans'
import matplotlib.pyplot as plt
import matplotlib as mpl

import DataTools
import os

if __name__ == '__main__':
    apply_penalty = False
    posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty)

    posterior_table_wgd = posterior_table[posterior_table['WGD_Status']]
    posterior_table_non_wgd = posterior_table[~posterior_table['WGD_Status']]

    wgd_calling_info = DataTools.get_wgd_calling_info()
    bad_samples = wgd_calling_info[(~wgd_calling_info['WGD_Status']) & (wgd_calling_info['Major_CN_Mode']==2)]['Sample_ID'].values

    posterior_table_non_wgd_filtered = posterior_table_non_wgd[~posterior_table_non_wgd['Sample_ID'].isin(bad_samples)]

    fig,axs = plt.subplots(2,1,figsize=(8,8),sharex=True)
    bins = np.linspace(0,1,101)
    axs[0].hist(posterior_table_non_wgd_filtered['Gain_Timing'],bins=bins,color='#6EC7ED',density=True)
    axs[1].hist(posterior_table_wgd['Gain_Timing'],bins=bins,color='#DB7F67',density=True)
    axs[0].set_ylabel('Density')
    axs[1].set_ylabel('Density')
    axs[1].set_xlabel('Gain Timing')
    axs[0].set_title('Non-WGD Samples')
    axs[1].set_title('WGD Samples')
    plt.tight_layout()
    plt.savefig('../plots/generic_gain_histogram/generic_gain_histogram.pdf')

