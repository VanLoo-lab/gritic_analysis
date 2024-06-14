
import pandas as pd
import numpy as np
import os
import DataTools
#font type 42
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Sans'
#font size 18
matplotlib.rcParams.update({'font.size': 18})
import matplotlib.pyplot as plt

import plotly.graph_objs as go
RNG = np.random.default_rng(42)


def get_sample_parsimony_table(route_table):
    route_table = route_table.copy()
    #add minimum number of events per sample,segment_id
    route_table.loc[:,'Min_N_Events'] = route_table.groupby(['Sample_ID','Segment_ID'])['Average_N_Events'].transform('min')
    #get extra number of events from minumum, rounded and cast to int
    route_table.loc[:,'Extra_N_Events'] = route_table['Average_N_Events'] - route_table['Min_N_Events']
    route_table.loc[:,'Extra_N_Events'] = route_table['Extra_N_Events'].round(0).astype(int)
    #classify extra events as 0,1,2 or 3+
    route_table.loc[:,'Extra_N_Events_Class'] = route_table['Extra_N_Events'].apply(lambda x: '0' if x == 0 else '1' if x == 1 else '2' if x == 2 else '3+')
    #get product of segment_width and probability
    #route_table.loc[:,'Segment_Width_Probability'] = route_table['Segment_Width'] * route_table['Probability']
    route_table['Norm_Segment_Width'] = route_table['Segment_Width']/route_table.groupby(['Sample_ID','Segment_ID','Chromosome'])['Segment_Width'].transform('sum')
    #route_table.loc[:,'Segment_Width_Probability'] = route_table['Probability']*route_table['Norm_Segment_Width']
    route_table.loc[:,'Segment_Width_Probability'] = route_table['Probability']
    return route_table
    



def get_cancer_type_parsimony_summary_bootstrap(cancer_type_parsimony_table,chromosome,n_bootstraps=250):
    table_store = []
    sample_table_dict = DataTools.get_sample_dict(cancer_type_parsimony_table)
    #adding the one true samples to the test
    p_value_store = {3:[False],4:[True]}
    for i in range(n_bootstraps):
        if i % 5 == 0:
            print(i,'/',n_bootstraps)
        bootstrap_table = DataTools.get_bootstrap_table(sample_table_dict)
        parsimony_summary = get_cancer_type_parsimony_summary(bootstrap_table,chromosome)

        table_store.append(parsimony_summary)
    parsimony_summary_bootstrap = pd.concat(table_store)
    #get 95% confidence interval
    parsimony_summary_bootstrap = parsimony_summary_bootstrap.groupby(['Major_CN','Location'])['Segment_Width_Probability'].agg(Low_CI=lambda x: np.quantile(x,0.05/2),High_CI=lambda x: np.quantile(x,1-0.05/2)).reset_index()

    return parsimony_summary_bootstrap

def get_cancer_type_parsimony_summary(cancer_type_parsimony_table,chromosome):
    cancer_type_parsimony_table = cancer_type_parsimony_table.copy()
    cancer_type_parsimony_table.loc[:,'Location'] = cancer_type_parsimony_table['Chromosome'].apply(lambda x: chromosome if x == chromosome else 'Elsewhere')
    #get class if extra_n_events_class is greater than 0
    cancer_type_parsimony_table.loc[:,'Extra_N_Events_Class'] = cancer_type_parsimony_table['Extra_N_Events_Class'].apply(lambda x: '1+' if x != '0' else '0')
    #aggregate total segment_width_probability per major_cn,location and normalize within major cn and location
    cancer_type_parsimony_summary = cancer_type_parsimony_table.groupby(['Major_CN','Location','Extra_N_Events_Class'])['Segment_Width_Probability'].sum().reset_index()
    
    cancer_type_parsimony_summary.loc[:,'Segment_Width_Probability'] = cancer_type_parsimony_summary['Segment_Width_Probability']/cancer_type_parsimony_summary.groupby(['Major_CN','Location'])['Segment_Width_Probability'].transform('sum')
    cancer_type_parsimony_summary = cancer_type_parsimony_summary[cancer_type_parsimony_summary['Extra_N_Events_Class']=='1+'].drop(columns=['Extra_N_Events_Class'])
    return cancer_type_parsimony_summary



def get_p_value_star(p_value):
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''
    
def plot_cancer_type_parsimony_summary_matplotlib(cancer_type_parsimony_table_summary,apply_penalty,p_values,chromosome,cancer_type):
    # matplotlib bar chart with error bars
    chromosome_colors = {'5':'#2EB3FF','7':'#0FA35A'}
    colors = {'Elsewhere':'grey',chromosome:chromosome_colors[chromosome]}
    fig,ax = plt.subplots(1,1,figsize=(10,5))
    for location,location_table in cancer_type_parsimony_table_summary.groupby(['Location']):
        location_table.loc[:,'Major_CN_Location'] = location_table['Major_CN'] + np.maximum(location_table['Major_CN']-3.75,0.0)
        if location == 'Elsewhere':
            offset = 0
        else:
            offset = 0.5
        
        ax.bar(location_table['Major_CN_Location']-offset,location_table['Segment_Width_Probability'],yerr=(location_table['Segment_Width_Probability']-location_table['Low_CI'],location_table['High_CI']-location_table['Segment_Width_Probability']),color=colors[location],label=location,width=0.5,align='edge',capsize=10)
    ax.set_xlabel('Major Copy Number')
    ax.set_ylabel('Average Probability Non-Parsimonious')

    ax.annotate(f'{get_p_value_star(p_values[0])}',xy=(0.3,0.9),xycoords='axes fraction',fontsize=24)
    ax.annotate(f'{get_p_value_star(p_values[1])}',xy=(0.7,0.9),xycoords='axes fraction',fontsize=24)
    ax.set_xticks([3,4.25])
    ax.set_xticklabels(['3','4'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig(f'../plots/parsimony_plots_cancer_type/{cancer_type}_np_prior_{apply_penalty}.pdf',bbox_inches='tight')


def get_cancer_type_diffs(cancer_type_summary_table,chromosome):
    
    
    three_table = cancer_type_summary_table[cancer_type_summary_table['Major_CN']==3]
    four_table = cancer_type_summary_table[cancer_type_summary_table['Major_CN']==4]
    three_diff = np.abs(three_table[three_table['Location']=='Elsewhere']['Segment_Width_Probability'].to_numpy()[0]-three_table[three_table['Location']==chromosome]['Segment_Width_Probability'].to_numpy()[0])
    four_diff = np.abs(four_table[four_table['Location']=='Elsewhere']['Segment_Width_Probability'].to_numpy()[0]-four_table[four_table['Location']==chromosome]['Segment_Width_Probability'].to_numpy()[0])
    return three_diff,four_diff

def get_random_chr_table(chr_store):
    chr_table = {'Sample_ID':[],'Old_Chromosome':[],'Chromosome':[]}
    for sample_id in chr_store:
        chr_table['Sample_ID'].extend([sample_id]*len(chr_store[sample_id]))
        chr_table['Old_Chromosome'].extend(chr_store[sample_id])
        chr_table['Chromosome'].extend(RNG.choice(chr_store[sample_id],len(chr_store[sample_id]),replace=False))
    chr_table = pd.DataFrame(chr_table)
    return chr_table

def get_cancer_type_p_values(cancer_type_parsimony_table,chromosome,n_bootstraps=10000):
    cancer_type_parsimony_table = cancer_type_parsimony_table.copy()
    cancer_type_parsimony_table_summary = get_cancer_type_parsimony_summary(cancer_type_parsimony_table,chromosome)
    three_diff,four_diff = get_cancer_type_diffs(cancer_type_parsimony_table_summary,chromosome)
    cancer_type_parsimony_table = cancer_type_parsimony_table.rename(columns={'Chromosome':'Old_Chromosome'})
    cancer_type_chr_table = cancer_type_parsimony_table[['Sample_ID','Old_Chromosome']].drop_duplicates()
    chr_store = {sample_id:cancer_type_chr_table[cancer_type_chr_table['Sample_ID']==sample_id]['Old_Chromosome'].to_list() for sample_id in cancer_type_chr_table['Sample_ID'].unique()}
    three_diff_store,four_diff_store = [],[]
    for i in range(n_bootstraps):
        chr_table = get_random_chr_table(chr_store)
        bootstrap_summary = get_cancer_type_parsimony_summary(cancer_type_parsimony_table.merge(chr_table,on=['Sample_ID','Old_Chromosome']),chromosome)
        bootstrap_summary = bootstrap_summary.copy()
        
    
        bootstrap_three_diff,bootstrap_four_diff = get_cancer_type_diffs(bootstrap_summary,chromosome)
        three_diff_store.append(bootstrap_three_diff)
        four_diff_store.append(bootstrap_four_diff)
       
    three_p_value = np.mean(three_diff_store>=three_diff)
    four_p_value = np.mean(four_diff_store>=four_diff)
    return three_p_value,four_p_value

def get_good_routes(route_table):
    route_table = route_table[(route_table['Major_CN'].isin([3,4])) & (route_table['Minor_CN']<=2)].copy()
    good_routes = route_table[['Major_CN','Minor_CN','Route','Average_N_Events']].groupby('Route').mean().reset_index()
    good_routes = good_routes.sort_values(['Major_CN','Minor_CN','Average_N_Events'],ascending=True).reset_index()
    
    final_routes = []
    for (major_cn,minor_cn),cn_table in good_routes.groupby(['Major_CN','Minor_CN']):
        max_n_events = np.sort(cn_table['Average_N_Events'].to_numpy())[1]+1e-6
        cn_table = cn_table[cn_table['Average_N_Events']<=max_n_events]
        final_routes.extend(cn_table['Route'].to_list())
    return final_routes

def load_parsimony_table(cancer_type,chromosome,apply_penalty):
    cancer_type_table = DataTools.load_route_table(cancer_types=[cancer_type],apply_penalty=apply_penalty)
    print('Total number of samples')
    print(len(cancer_type_table[['Sample_ID']].drop_duplicates()))
    print('--')
    cancer_type_table = cancer_type_table[cancer_type_table['Major_CN'].isin([3,4]) & (cancer_type_table['Minor_CN']<=2)].copy()
    
    good_routes = get_good_routes(cancer_type_table)


    cancer_type_table = cancer_type_table[cancer_type_table['Route'].isin(good_routes)]

    cancer_type_table['Probability_Sum'] = cancer_type_table.groupby(['Sample_ID','Segment_ID'])['Probability'].transform('sum')
    print('Total number of segments, total segment width')
    print(len(cancer_type_table[['Sample_ID','Segment_ID']].drop_duplicates()))
    print(cancer_type_table[['Sample_ID','Segment_ID','Segment_Width']].drop_duplicates()['Segment_Width'].sum())
    print('filtered number of segments, total segment width')
    cancer_type_table = cancer_type_table[cancer_type_table['Probability_Sum']>=0.5]
    print(len(cancer_type_table[['Sample_ID','Segment_ID']].drop_duplicates()))
    print(cancer_type_table[['Sample_ID','Segment_ID','Segment_Width']].drop_duplicates()['Segment_Width'].sum())
    print('----------')
    cancer_type_table['Probability'] = cancer_type_table['Probability']/cancer_type_table['Probability_Sum']
    
    cancer_type_parsimony_table = get_sample_parsimony_table(cancer_type_table)
   
    return cancer_type_parsimony_table
def write_sample_parsimony_plot_table(cancer_type_parsimony_table,cancer_type,apply_penalty,chromosome,plot_dir):
    arm_gain_table_path  = f'../output/arm_pre_post_clean/gain_timing_{apply_penalty}_bin_counts.tsv'
    arm_gain_table = pd.read_csv(arm_gain_table_path,sep='\t')

    arm_loss_table_path  = f'../output/arm_pre_post_clean/loss_timing_bin_counts.tsv'
    arm_loss_table = pd.read_csv(arm_loss_table_path,sep='\t')
    arm_table = arm_gain_table.merge(arm_loss_table,how='inner')
    arm_table = arm_table[arm_table['Cancer_Type']==cancer_type]
    arm_table['In_Parsimony_Table'] = arm_table['Sample_ID'].apply(lambda x: x in cancer_type_parsimony_table['Sample_ID'].unique())
    arm_table_out_path = f'{plot_dir}/arm_timing_table.tsv'
    arm_table.to_csv(arm_table_out_path,sep='\t',index=False)
    
def get_sample_plots(cancer_type_parsimony_table,chromosome,cancer_type,apply_penalty):
    plot_dir = f'../plots/sample_parsimony_plots_cancer_type/{cancer_type.replace(" ","_")}_{chromosome}_prior_{apply_penalty}/'
    os.makedirs(plot_dir,exist_ok=True)
    for sample_id in cancer_type_parsimony_table['Sample_ID'].unique():
        sample_plot_path = f'../plots/sample_timing_plots/prior_{apply_penalty}/{sample_id}_sample_timing.pdf'
        if os.path.exists(sample_plot_path):
            os.system(f'cp {sample_plot_path} {plot_dir}/{sample_id}_sample_timing.pdf')
    write_sample_parsimony_plot_table(cancer_type_parsimony_table,cancer_type,apply_penalty,chromosome,plot_dir)

if __name__ == '__main__':
    for (cancer_type,chromosome) in [('Kidney renal clear cell carcinoma','5'),('Glioblastoma multiforme','7')]:
        for apply_penalty in [True,False]:
            if cancer_type == 'Glioblastoma multiforme' and chromosome == '7':
                continue
            cancer_type_parsimony_table = load_parsimony_table(cancer_type,chromosome,apply_penalty)
            print(cancer_type,chromosome,apply_penalty)
            print(len(cancer_type_parsimony_table['Sample_ID'].unique()))
            #get_sample_plots(cancer_type_parsimony_table,chromosome,cancer_type,apply_penalty)
            cancer_type_parsimony_table_summary = get_cancer_type_parsimony_summary(cancer_type_parsimony_table,chromosome)
            

            cancer_type_p_values = get_cancer_type_p_values(cancer_type_parsimony_table,chromosome)
            print('------')
            print(apply_penalty,cancer_type,chromosome)
            print(cancer_type_p_values)
            print('-----')
            #cancer_type_p_values.to_csv(f"../output/cancer_type_parsimony/{cancer_type.replace(' ','_')}_{chromosome}_prior_{apply_penalty}_p_values_prob_only.tsv",sep='\t',index=False)

            cancer_type_parsimony_table_bootstrap = get_cancer_type_parsimony_summary_bootstrap(cancer_type_parsimony_table,chromosome,n_bootstraps=250)
            cancer_type_parsimony_table_summary = cancer_type_parsimony_table_summary.merge(cancer_type_parsimony_table_bootstrap,on=['Major_CN','Location'])
            cancer_type_parsimony_table_summary.to_csv(f"../output/cancer_type_parsimony/{cancer_type.replace(' ','_')}_{chromosome}_prior_{apply_penalty}_bootstrap_values_prob_only.tsv",sep='\t',index=False)

            plot_cancer_type_parsimony_summary_matplotlib(cancer_type_parsimony_table_summary,apply_penalty,cancer_type_p_values,chromosome,cancer_type.replace(' ','_'))
            print('done',cancer_type,chromosome,apply_penalty)