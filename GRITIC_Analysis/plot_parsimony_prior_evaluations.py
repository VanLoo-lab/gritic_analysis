import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#font size 18
plt.rcParams.update({'font.size': 18})
#pdf font type 42
plt.rcParams['pdf.fonttype'] = 42
#nimbus sans font
plt.rcParams['font.sans-serif'] = "Nimbus Sans"

if __name__ == '__main__':
    plot_out_dir = '../plots/parsimony_penalty_evaluations'

    agg_table = pd.read_csv('../output/evaluate_parsimony_penalty_out/agg_table.tsv',sep='\t').sort_values(by=['L_Values'])
    prop_table = pd.read_csv('../output/evaluate_parsimony_penalty_out/prop_table.tsv',sep='\t').sort_values(by=['L_Values'])

    fig,ax = plt.subplots(1,2,figsize=(20,5))
    ax[0].plot(agg_table['L_Values'],agg_table['Parsimony_Prob'],'#19138f',label='Parsimony',lw=5)
    ax[0].set_ylim(0,1.04)
    ax[0].set_title('Proportion of probability on parsimonious routes')
    ax[0].set_xlabel('Penalty term on additional events (L)')
    ax[0].axvline(x=2.7,ls='--',color='red',lw=5)
    ax[1].plot(agg_table['L_Values'],agg_table['Prop_Accuracy'],'#a60f2b',label='Accuracy',lw=5)
    ax[1].set_ylim(0,1.04)
    ax[1].set_title('Proportion of segments with maximum probability on true route')
    ax[1].set_xlabel('Penalty term on additional events (L)')
    ax[1].axvline(x=2.7,ls='--',color='red',lw=5)
    plt.tight_layout()
    plt.savefig(f'{plot_out_dir}/parsimony_penalty_evaluations_combined.pdf')


    prop_table['CN_State'] = prop_table['Major_CN'].astype(str) + '_' + prop_table['Minor_CN'].astype(str)

    prop_table = prop_table.sort_values(by=['Major_CN','Minor_CN','L_Values'])

    fig,axs = plt.subplots(6,5,figsize=(20,20),sharex=True,sharey=True)
    axs = axs.flatten()
    processed_count = 0
    for i,(cn_state,cn_state_table) in enumerate(prop_table.groupby('CN_State')):
        axs[i].plot(cn_state_table['L_Values'],cn_state_table['Parsimony_Prob']/cn_state_table['N_Segments'],color='#19138f',lw=5)
        
        axs[i].set_title(f'{cn_state.replace("_","+")}')
        axs[i].axvline(x=2.7,ls='--',color='red',lw=5)
        processed_count += 1
    #remove empty subplots
    for i in range(processed_count,len(axs)):
        axs[i].axis('off')


    fig.supxlabel('Penalty term on additional events (L)')
    fig.supylabel('Proportion of probability on parsimonious routes')
    plt.tight_layout()
    plt.savefig(f'{plot_out_dir}/parsimony_penalty_evaluations_probs_cn_state.pdf')

    fig,axs = plt.subplots(6,5,figsize=(20,20),sharex=True,sharey=True)
    axs = axs.flatten()
    processed_count = 0
    for i,(cn_state,cn_state_table) in enumerate(prop_table.groupby('CN_State')):
        axs[i].plot(cn_state_table['L_Values'],cn_state_table['N_Accurate']/cn_state_table['N_Segments'],color='#a60f2b',lw=5)
        
        axs[i].set_title(f'{cn_state.replace("_","+")}')
        axs[i].axvline(x=2.7,ls='--',color='red',lw=5)
        processed_count += 1
    #remove empty subplots
    for i in range(processed_count,len(axs)):
        axs[i].axis('off')

    fig.supxlabel('Penalty term on additional events (L)')
    fig.supylabel('Proportion of segments with maximum probability on true route')
    plt.tight_layout()
    plt.savefig(f'{plot_out_dir}/parsimony_penalty_evaluations_accuracy_cn_state.pdf')