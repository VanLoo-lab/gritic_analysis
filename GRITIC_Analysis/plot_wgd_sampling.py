import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Sans'
import DataTools
import os
plt.rcParams.update({'font.size': 21})



def plot_bootstrap_hists(cancer_type,tumor_type,apply_penalty,min_mutations,wgd_bins,bin_method,n_bins):

    out_dir = '../output/pan_wgd_sampling'
    base_out_path = f"{out_dir}/{cancer_type.replace(' ','_')}_{tumor_type.replace(' ','_')}_apply_penalty_{apply_penalty}_min_mutations_{min_mutations}_bin_method_{bin_method}_n_bins_{n_bins}"


    bootstrap_hists = np.load(f'{base_out_path}_bootstrap_hists.npy',allow_pickle=True)
    true_hists = np.load(f'{base_out_path}_true_hists.npy',allow_pickle=True)
    true_hists = [np.array(true_hists[i]).astype(float) for i in range(true_hists.shape[0])]
    true_bins = np.load(f'{base_out_path}_true_bins.npy',allow_pickle=True)
    true_bins = [np.array(true_bins[i]).astype(float) for i in range(true_bins.shape[0])]

    
    n_samples  =np.load(f'{base_out_path}_n_samples_true.npy')
    bootstrap_hists = [[bootstrap_hists[i,j] for i in range(bootstrap_hists.shape[0])] for j in range(bootstrap_hists.shape[1])]

    boostrap_hist_store_low_ci =[np.percentile(bootstrap_hists[i],5/2,axis=0).astype(float) for i in range(n_bins)]
    boostrap_hist_store_high_ci = [np.percentile(bootstrap_hists[i],100-5/2,axis=0).astype(float) for i in range(n_bins)]



    plt_colors = ['#EF476F','#FFB60A','#118AB2','#66008F']
    fill_colors = ['#F7A1B5','#FFD470','#8FDBF5','#E299FF']

    plt_title = cancer_type
    if cancer_type =='All':
        plt_title = 'All tumors'
    if tumor_type != 'All':
        plt_title += f' {tumor_type}'
    if apply_penalty:
        plt_title += ' with non-parsimony penalty'
    else:
        plt_title += ' without non-parsimony penalty'
    fig,ax = plt.subplots(1,1,figsize=(12,8))
    if bin_method == 'Even_Spacing':
        ax.axvspan(-0.05,0.05,color='#a3a3a3',alpha=0.5)


    for i in range(len(true_hists)):
        
        ax.stairs(np.array(boostrap_hist_store_high_ci[i]),np.array(true_bins[i]),baseline=np.array(boostrap_hist_store_low_ci[i]),color=fill_colors[i],fill=True,alpha=0.5)
        ax.stairs(np.array(true_hists[i]),np.array(true_bins[i]),label=f'{wgd_bins[i]:.2f} - {wgd_bins[i+1]:.2f} (n={n_samples[i]})',color=plt_colors[i],baseline=None,lw=2,alpha=0.6)
    ax.legend(title='WGD Timing')
    ax.set_title(plt_title)
    ax.set_xlabel('Gain Timing - WGD Timing')
    ax.set_ylabel('Relative Gain Rate')
    if n_bins ==4:
        ax.set_ylim(0,7.5)
    else:
        ax.set_ylim(0,5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)   
    
    out_dir = f'../plots/pan_wgd_gain_sampling/apply_penalty_{apply_penalty}_min_mutations_{min_mutations}_bin_method_{bin_method}_n_bins_{n_bins}'
    os.makedirs(out_dir,exist_ok=True)
    out_path  = f'{out_dir}/{cancer_type}_{tumor_type}_pan_wgd_sampling'.replace(' ','_')
    plt.savefig(out_path+'.png')
    plt.savefig(out_path+'.pdf')
    plt.close()
    print(f'Saved {out_path}')
