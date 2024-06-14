import pandas as pd
import numpy as np

import DataTools
#font type 42
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Sans'
#font size 18
matplotlib.rcParams.update({'font.size': 26})
import matplotlib.pyplot as plt

import plotly.graph_objs as go
RNG = np.random.default_rng(42)
def get_sample_parsimony_table(route_table):
    route_table = route_table.copy()
    #add minimum number of events per sample,segment_id
    route_table['Min_N_Events'] = route_table.groupby(['Sample_ID','Segment_ID'])['Average_N_Events'].transform('min')
    #get  number of events from minumum, rounded and cast to int
    route_table['Extra_N_Events'] = route_table['Average_N_Events'] - route_table['Min_N_Events']
    route_table['Extra_N_Events'] = route_table['Extra_N_Events'].round(0).astype(int)
    #classify extra events as 0,1,2 or 3+
    route_table['Extra_N_Events_Class'] = route_table['Extra_N_Events'].apply(lambda x: '0' if x == 0 else '1' if x == 1 else '2' if x == 2 else '3+')
    #get product of segment_width and probability
    #route_table['Segment_Width_Probability'] = route_table['Segment_Width'] * route_table['Probability']
    route_table['Segment_Width_Probability'] = route_table['Probability']
    return route_table
    
def get_parsimony_summary(sample_parsimony_table):

    #aggregate total segment_width_probability per extra_n_events_class and normalize
    parsimony_summary = sample_parsimony_table.groupby(['Extra_N_Events_Class'])['Segment_Width_Probability'].sum().reset_index()
    parsimony_summary['Segment_Width_Probability'] = parsimony_summary['Segment_Width_Probability']/parsimony_summary['Segment_Width_Probability'].sum()
    return parsimony_summary

def get_parsimony_summary_bootstrap(sample_parsimony_table,n_bootstraps=250):
    table_store = []
    sample_table_dict = DataTools.get_sample_dict(sample_parsimony_table)
    for i in range(n_bootstraps):
        if i % 1 == 0:
            print(i,'/',n_bootstraps)
        bootstrap_table = DataTools.get_bootstrap_table(sample_table_dict)
        parsimony_summary = get_parsimony_summary(bootstrap_table)
        table_store.append(parsimony_summary)
    parsimony_summary_bootstrap = pd.concat(table_store)
    #get 95% confidence interval
    parsimony_summary_bootstrap = parsimony_summary_bootstrap.groupby(['Extra_N_Events_Class'])['Segment_Width_Probability'].agg(Low_CI=lambda x: np.quantile(x,0.05/2),High_CI=lambda x: np.quantile(x,1-0.05/2)).reset_index()
    return parsimony_summary_bootstrap
def get_kidney_parsimony_summary_bootstrap(kidney_parsimony_table,n_bootstraps=250):
    table_store = []
    sample_table_dict = DataTools.get_sample_dict(kidney_parsimony_table)
    #adding the one true samples to the test
    p_value_store = {3:[False],4:[True]}
    for i in range(n_bootstraps):
        if i % 5 == 0:
            print(i,'/',n_bootstraps)
        bootstrap_table = DataTools.get_bootstrap_table(sample_table_dict)
        parsimony_summary = get_kidney_parsimony_summary(bootstrap_table)
        
        #convert Segment_Width_Probability to dict with tuple index, major_cn
        parsimony_summary_dict = parsimony_summary.set_index(['Major_CN','Location',])['Segment_Width_Probability'].to_dict()

        for major_cn in [3,4]:
            p_value_store[major_cn].append(parsimony_summary_dict[(major_cn,'5')]<parsimony_summary_dict[(major_cn,'Elsewhere')])
        table_store.append(parsimony_summary)
    parsimony_summary_bootstrap = pd.concat(table_store)
   
    #get 95% confidence interval
    parsimony_summary_bootstrap = parsimony_summary_bootstrap.groupby(['Major_CN','Location'])['Segment_Width_Probability'].agg(Low_CI=lambda x: np.quantile(x,0.05/2),High_CI=lambda x: np.quantile(x,1-0.05/2)).reset_index()
    return parsimony_summary_bootstrap

def get_kidney_parsimony_summary(kidney_parsimony_table):
    #filter to major cn 3 and 4
    kidney_parsimony_table = kidney_parsimony_table[kidney_parsimony_table['Major_CN'].isin([3,4])].copy()

    #set class as to chromosome equals 5 or not
    kidney_parsimony_table['Location'] = kidney_parsimony_table['Chromosome'].apply(lambda x: '5' if x == '5' else 'Elsewhere')

    #get class if extra_n_events_class is greater than 0
    kidney_parsimony_table['Extra_N_Events_Class'] = kidney_parsimony_table['Extra_N_Events_Class'].apply(lambda x: '1+' if x != '0' else '0')
    #aggregate total segment_width_probability per major_cn,location and normalize within major cn and location
    kidney_parsimony_summary = kidney_parsimony_table.groupby(['Major_CN','Location','Extra_N_Events_Class'])['Segment_Width_Probability'].sum().reset_index()
    
    kidney_parsimony_summary['Segment_Width_Probability'] = kidney_parsimony_summary['Segment_Width_Probability']/kidney_parsimony_summary.groupby(['Major_CN','Location'])['Segment_Width_Probability'].transform('sum')
    kidney_parsimony_summary = kidney_parsimony_summary[kidney_parsimony_summary['Extra_N_Events_Class']=='1+'].drop(columns=['Extra_N_Events_Class'])
    return kidney_parsimony_summary


def plot_parsimony_summary_matplotlib(parsimony_summary,control_summary,apply_penalty,major_cn,sample_p_value_table,control_only=False):
    #a matplotlib bar chart with error bars
    fig,ax = plt.subplots(figsize=(14,7))
    width= 0.4
    plt_points = np.arange(len(parsimony_summary['Extra_N_Events_Class']))
    if control_only:
        ax.bar(plt_points-width/2,control_summary['Segment_Width_Probability'],color='grey',label='Parsimony Control',yerr=(control_summary['Segment_Width_Probability']-control_summary['Low_CI'],control_summary['High_CI']-control_summary['Segment_Width_Probability']),width=width,align='edge',capsize=10)
    else:
        ax.bar(plt_points-width,control_summary['Segment_Width_Probability'],color='grey',label='Parsimony Control',yerr=(control_summary['Segment_Width_Probability']-control_summary['Low_CI'],control_summary['High_CI']-control_summary['Segment_Width_Probability']),width=width,align='edge',capsize=10)
        ax.bar(plt_points,parsimony_summary['Segment_Width_Probability'],color='#d91437',label='WGD Tumor Cohort',yerr=(parsimony_summary['Segment_Width_Probability']-parsimony_summary['Low_CI'],parsimony_summary['High_CI']-parsimony_summary['Segment_Width_Probability']),width=width,align='edge',capsize=10)
    ax.set_xlabel('Events above Parsimony')
    ax.set_ylabel('Average Probability')
    max_height = np.maximum(parsimony_summary['High_CI'].max(),control_summary['High_CI'].max())
    ax.set_ylim(0,max_height+0.1)
    #ax.legend()
    ax.set_xticks(plt_points)
    ax.set_xticklabels(parsimony_summary['Extra_N_Events_Class'].astype(str))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #annotate with p-value stars
    if not control_only:
        for i in range(len(parsimony_summary)):
            p_value = sample_p_value_table.iloc[i]['P_Value']
            height = max(parsimony_summary.iloc[i]['High_CI'],control_summary.iloc[i]['High_CI'])
            ax.annotate(f'{get_p_value_star(p_value)}',(i,height+0.05),ha='center',va='center',fontsize=24)
    if control_only:
        out_path = f'../plots/parsimony_plots/parsimony_proportions_{apply_penalty}_control_only_major_cn_{major_cn}_matplotlib.pdf'
    else:
        out_path = f'../plots/parsimony_plots/parsimony_proportions_{apply_penalty}_major_cn_{major_cn}_matplotlib.pdf'
    plt.savefig(out_path,bbox_inches='tight')


def get_p_value_star(p_value):
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''


def get_parsimony_p_values(parsimony_table,control_parsimony_table):
    extra_events_class = ['0','1','2','3+']
    p_value_table ={'Extra_N_Events_Class':[],'P_Value':[]}
    for extra_events in extra_events_class:
        parsimony_events_subset = parsimony_table[parsimony_table['Extra_N_Events_Class']==extra_events]['Segment_Width_Probability'].to_numpy()
        control_parsimony_events_subset = control_parsimony_table[control_parsimony_table['Extra_N_Events_Class']==extra_events]['Segment_Width_Probability'].to_numpy()
        true_diff = np.abs(np.mean(parsimony_events_subset)-np.mean(control_parsimony_events_subset))
        combined_parsimony = np.concatenate((parsimony_events_subset,control_parsimony_events_subset))
        p_values = []
        for i in range(10000):

            RNG.shuffle(combined_parsimony)
            diff = np.abs(np.mean(combined_parsimony[:len(parsimony_events_subset)])-np.mean(combined_parsimony[len(parsimony_events_subset):]))
            p_values.append(diff>=true_diff)
        p_value_table['Extra_N_Events_Class'].append(extra_events)
        p_value_table['P_Value'].append(np.mean(p_values))
    return pd.DataFrame(p_value_table)

def get_kidney_diffs(kidney_summary_table):
    three_table = kidney_summary_table[kidney_summary_table['Major_CN']==3]
    four_table = kidney_summary_table[kidney_summary_table['Major_CN']==4]
    three_diff = np.abs(three_table[three_table['Location']=='Elsewhere']['Segment_Width_Probability'].to_numpy()[0]-three_table[three_table['Location']=='5']['Segment_Width_Probability'].to_numpy()[0])
    four_diff = np.abs(four_table[four_table['Location']=='Elsewhere']['Segment_Width_Probability'].to_numpy()[0]-four_table[four_table['Location']=='5']['Segment_Width_Probability'].to_numpy()[0])
    return three_diff,four_diff

def get_random_chr_table(chr_store):
    chr_table = {'Sample_ID':[],'Old_Chromosome':[],'Chromosome':[]}
    for sample_id in chr_store:
        chr_table['Sample_ID'].extend([sample_id]*len(chr_store[sample_id]))
        chr_table['Old_Chromosome'].extend(chr_store[sample_id])
        chr_table['Chromosome'].extend(RNG.choice(chr_store[sample_id],len(chr_store[sample_id]),replace=False))
    chr_table = pd.DataFrame(chr_table)
    return chr_table

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

for apply_penalty in [False,True]:
    for major_cn in ['All',3,4,5,6,7]:
        if major_cn != 'All':
            continue
        
        control_dataset = 'state_validation_all_routes'
        route_info = pd.read_csv('../../Analysis/output/all_routes_n_events.tsv',sep='\t')
        route_info = route_info[route_info['WGD_Status'] & (~route_info['Parsimonious'])]
        route_info['Route'] = route_info['Route'].str.slice(0,9)
        route_info = route_info[['Major_CN','Minor_CN','Route']]
        control_metadata_table = pd.read_csv(f'../../output/{control_dataset}/complete_metadata_table_{control_dataset}.tsv',sep='\t')
        control_metadata_table['Route'] = control_metadata_table['Route'].str.slice(0,9)
        
        control_metadata_table = control_metadata_table.merge(route_info,how='inner')
        control_metadata_table = control_metadata_table[['Sample_ID','Segment_ID']]
        control_metadata_table['Segment_ID'] = control_metadata_table['Segment_ID'].str.replace('Segment_', '')
        

        route_table = DataTools.load_route_table(apply_penalty=apply_penalty)

        route_table = route_table[route_table['N_Mutations']>=20]

        
        control_route_table = pd.read_csv(f'../../output/{control_dataset}/complete_route_table_{control_dataset}.tsv',sep='\t',dtype={'Sample_ID':str,'Segment_ID':str,'Chromosome':str})
        control_route_table = control_route_table[control_route_table['WGD_Status']]
        control_route_table = control_route_table[control_route_table['Major_CN']>2]
        control_route_table = control_route_table[control_route_table['N_Mutations']>=20]
        control_route_table = control_route_table.merge(control_metadata_table,how='inner')

        #non-parsimonious routes only
        
        
        control_route_table = control_route_table[['Sample_ID','Segment_ID','Segment_Start','Segment_End','Route','N_Mutations','Average_Pre_WGD_Losses', 'Average_Post_WGD_Losses','Average_N_Events','Major_CN','Minor_CN','Probability']].drop_duplicates()
        control_route_table = control_route_table[control_route_table['Major_CN']>2]
   
        if major_cn != 'All':
            route_table = route_table[route_table['Major_CN']==major_cn]
            control_route_table = control_route_table[control_route_table['Major_CN']==major_cn]

        control_route_table['Segment_Width'] = control_route_table['Segment_End']-control_route_table['Segment_Start']
        if apply_penalty:
            control_route_table['Probability'] = DataTools.apply_penalty_quick(control_route_table,2.7)


        sample_parsimony_table  =get_sample_parsimony_table(route_table)

        control_parsimony_table = get_sample_parsimony_table(control_route_table)
        

        #aggregate total segment_width_probability per sample,extra_n_events_class
        sample_parsimony_table = sample_parsimony_table.groupby(['Sample_ID','Extra_N_Events_Class'])['Segment_Width_Probability'].sum().reset_index()
        control_parsimony_table = control_parsimony_table.groupby(['Sample_ID','Extra_N_Events_Class'])['Segment_Width_Probability'].sum().reset_index()

        parsimony_summary = get_parsimony_summary(sample_parsimony_table)
        control_summary = get_parsimony_summary(control_parsimony_table)
       
        sample_p_value_table= get_parsimony_p_values(sample_parsimony_table,control_parsimony_table)
    

        
        parsimony_summary_bootstrap =get_parsimony_summary_bootstrap(sample_parsimony_table,n_bootstraps=250)
        control_summary_bootstrap = get_parsimony_summary_bootstrap(control_parsimony_table,n_bootstraps=250)
        parsimony_summary = parsimony_summary.merge(parsimony_summary_bootstrap,on='Extra_N_Events_Class')

        control_summary = control_summary.merge(control_summary_bootstrap,on='Extra_N_Events_Class')

        plot_parsimony_summary_matplotlib(parsimony_summary,control_summary,apply_penalty,major_cn,sample_p_value_table)
        plot_parsimony_summary_matplotlib(parsimony_summary,control_summary,apply_penalty,major_cn,sample_p_value_table,control_only=True)

