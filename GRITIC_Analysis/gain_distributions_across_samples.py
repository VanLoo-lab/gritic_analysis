import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Sans'
import matplotlib.pyplot as plt
#font size 22
plt.rcParams.update({'font.size': 22})
import DataTools

def get_color_scheme(most_common_chromosomes):
    #a rainbow 10 color palette
    color_pallette = ['#e6194b','#3cb44b','#ffe119','#4363d8','#f58231','#911eb4','#46f0f0','#f032e6','#bcf60c','#fabebe']
    np.random.shuffle(color_pallette)
   
    color_scheme = {'NA':'#000000'}
    colors = []
    for chromosome in most_common_chromosomes:
        if chromosome not in color_scheme:
            color_scheme[chromosome] = color_pallette.pop()
        
        colors.append(color_scheme[chromosome])
    return colors

def get_sample_histograms(posterior_table):
    sample_histograms = {}
    sample_bins = np.linspace(0,1,11)
    for sample_id,sample_table in posterior_table.groupby('Sample_ID'):
        sample_table = sample_table.copy()
        hist = np.histogram(sample_table['Gain_Timing'],bins=sample_bins,weights=sample_table['Segment_Width'])[0]
        hist = hist/np.sum(hist)
        
        full_chr_width = sample_table[['Chromosome','Segment_ID','Segment_Width']].drop_duplicates().groupby(['Chromosome']).agg({'Segment_Width':'sum'}).reset_index().rename(columns={'Segment_Width':'Full_Chromosome_Width'})
        norm_chr_width = sample_table[['Chromosome','Segment_ID','Segment_Width']].drop_duplicates().merge(full_chr_width)
        norm_chr_width.loc[:,'Norm_Segment_Width'] = norm_chr_width['Segment_Width']/norm_chr_width['Full_Chromosome_Width']
        sample_table = sample_table.merge(norm_chr_width[['Segment_ID','Norm_Segment_Width']],on='Segment_ID')
        most_common_chromosomes = []
        for i in range(10):
            sample_table_bin = sample_table[(sample_table['Gain_Timing']>=sample_bins[i])&(sample_table['Gain_Timing']<sample_bins[i+1])]
            #most common chromosome
            x = sample_table_bin.groupby('Chromosome').agg({'Norm_Segment_Width':'sum'}).reset_index().sort_values(by='Norm_Segment_Width',ascending=False)
            if len(x)>0:
                most_common_chromosome = x.index[0]
            else:
                most_common_chromosome = 'NA'
            most_common_chromosomes.append(most_common_chromosome)
        color_scheme = get_color_scheme(most_common_chromosomes)
        sample_histograms[sample_id] = {'hist':hist,'colors':color_scheme}
    return sample_histograms

def get_sample_summary(posterior_table):
    sample_summary = posterior_table.groupby(['Sample_ID','Posterior_Sample_Index']).agg({'WGD_Timing':np.median,'Segment_Width':np.sum,'Gain_Timing':np.median}).reset_index()
    sample_summary = sample_summary.groupby('Sample_ID').agg({'WGD_Timing':np.mean,'Segment_Width':np.mean,'Gain_Timing':np.median}).reset_index()
    return sample_summary

def plot_histograms(sample_histograms,sample_summary,cancer_type,apply_penalty):
    row_size = np.floor(np.sqrt(len(sample_histograms))).astype(int)
    col_size = np.ceil(len(sample_histograms)/row_size).astype(int)
    sample_ids = list(sample_summary.sort_values(by='WGD_Timing')['Sample_ID'].values)
    fig,ax = plt.subplots(row_size,col_size,figsize=(20,20))
    i = 0
    for row in range(row_size):
        start_col =i
        row_sample_ids =sample_ids[i:i+col_size]
        row_table = sample_summary[sample_summary['Sample_ID'].isin(row_sample_ids)].sort_values(by='Gain_Timing')
        row_sample_ids = list(row_table['Sample_ID'].values)
        for col in range(col_size):
            ax_i = ax[row,col]
            plotting = i < len(sample_ids)
            if plotting:
                sample_id = row_sample_ids[i-start_col]
                #ax_i.bar(np.linspace(0,0.9,10),sample_histograms[sample_id]['hist'],width=0.1,align='edge',color=sample_histograms[sample_id]['colors'],alpha=0.5)
                ax_i.bar(np.linspace(0,0.9,10),sample_histograms[sample_id]['hist'],width=0.1,align='edge',color='#2d915f',edgecolor='black',linewidth=1)
                ax_i.axvline(sample_summary[sample_summary['Sample_ID']==sample_id]['WGD_Timing'].values[0],color='red',linestyle='--',lw=4)
            #no axes
            ax_i.set_xticks([])
            ax_i.set_yticks([])
            i+=1
    fig.suptitle(cancer_type)
    fig.supxlabel('Gain Timing')
    fig.supylabel('Frequency')
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.savefig(f'../plots/all_sample_hists/{cancer_type.replace(" ","_")}_all_WGD_timing_hists_prior_{apply_penalty}.pdf')
    plt.savefig(f'../plots/all_sample_hists/{cancer_type.replace(" ","_")}_all_WGD_timing_hists_prior_{apply_penalty}.png')
    
    plt.close(fig)

if __name__ == '__main__':
    for apply_penalty in [False,True]:
        posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty,min_mutations=0)
        posterior_table = posterior_table[posterior_table['WGD_Status']]



        for cancer_type,cancer_type_posterior_table in posterior_table.groupby('Cancer_Type'):
            if len(cancer_type_posterior_table['Sample_ID'].unique())<20:
                continue
            sample_summary = get_sample_summary(cancer_type_posterior_table)

            sample_histograms = get_sample_histograms(cancer_type_posterior_table)
            plot_histograms(sample_histograms,sample_summary,cancer_type,apply_penalty)
            print('Done with',cancer_type)

        RNG = np.random.default_rng(1234)
        random_ids = RNG.choice(posterior_table['Sample_ID'].unique(),size=36,replace=False)

        random_cancer_type_posterior_table = posterior_table[posterior_table['Sample_ID'].isin(random_ids)]
        sample_summary = get_sample_summary(random_cancer_type_posterior_table)
        plot_histograms(get_sample_histograms(random_cancer_type_posterior_table),sample_summary,'Random',apply_penalty)