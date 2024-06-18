import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#font type 42
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
#font size 18
plt.rcParams.update({'font.size': 18})
#nimbus sans font
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Nimbus Sans'
RNG = np.random.default_rng(42)
import DataTools

from statsmodels.stats.proportion import proportions_ztest, proportion_confint

from scipy.stats import mannwhitneyu

primary_color = "#4AC1EB"
metastatic_color = "#334697" 

def get_p_value_star(p_value):
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''


def load_cn_table():
    good_samples = DataTools.get_good_samples()
    cn_table = pd.read_csv('../output/cn_table_status.tsv', sep='\t')
    cn_table['Segment_Width']+=1
    cn_table = cn_table[cn_table['Sample_ID'].isin(good_samples)]
    cn_table = cn_table[cn_table['Chromosome'].isin(list(map(str,range(1,23))))].copy()

    cn_table['Major_CN_Str'] = np.where(cn_table['Major_CN'] >=8, '8+', cn_table['Major_CN'].astype(str)) 
    return cn_table

def get_wgd_agg_values(wgd_status_cancer_type_table):
    wgd_status_cancer_type_agg = wgd_status_cancer_type_table.groupby(['Cancer_Type', 'Tumor_Type'])['WGD_Status'].mean().reset_index()
    wgd_status_cancer_type_agg = wgd_status_cancer_type_agg.groupby('Tumor_Type')['WGD_Status'].mean().reset_index()
    return wgd_status_cancer_type_agg

def get_bootstrap_dict(wgd_status_cancer_type_table):
    bootstrap_dict = {}
    for cancer_type,cancer_type_table in wgd_status_cancer_type_table.groupby('Cancer_Type'):
        bootstrap_dict[cancer_type] = DataTools.get_sample_dict(cancer_type_table)
    return bootstrap_dict

def get_bootstrap_agg(bootstrap_dict):
    bootstrap_table = []
    for cancer_type in bootstrap_dict:
        bootstrap_table.append(DataTools.get_bootstrap_table(bootstrap_dict[cancer_type]))
    bootstrap_table = pd.concat(bootstrap_table)
    bootstrap_table_agg = get_wgd_agg_values(bootstrap_table)
    return bootstrap_table_agg

def get_wgd_bootstrap_table(wgd_status_cancer_type_table, num_bootstrap=5):
    bootstrap_dict = get_bootstrap_dict(wgd_status_cancer_type_table)
    bootstrap_tables = []
    for i in range(num_bootstrap):
        bootstrap_table = get_bootstrap_agg(bootstrap_dict)
        bootstrap_table['bootstrap'] = i
        bootstrap_tables.append(bootstrap_table)
    bootstrap_tables = pd.concat(bootstrap_tables)
    bootstrap_ci = bootstrap_tables.groupby(['Tumor_Type'])['WGD_Status'].quantile([0.05/2, 1-0.05/2]).reset_index()
    bootstrap_ci = bootstrap_ci.pivot(index='Tumor_Type', columns='level_1', values='WGD_Status').reset_index()
    bootstrap_ci.columns = ['Tumor_Type', 'CI_Low', 'CI_High']
    bootstrap_ci = bootstrap_ci.sort_values(by=['Tumor_Type'], ascending=False)
    return bootstrap_ci

def get_diff(agg_table):
    return agg_table[agg_table['Tumor_Type']=='Metastatic']['WGD_Status'].values[0] - agg_table[agg_table['Tumor_Type']=='Primary']['WGD_Status'].values[0]

def get_permuted_table(cancer_type_dict):
    permutations = []
    for cancer_type in cancer_type_dict.keys():
        permuted_cancer_type_table = cancer_type_dict[cancer_type].copy()
        permuted_cancer_type_table['Tumor_Type'] = RNG.permutation(permuted_cancer_type_table['Tumor_Type'].values)
        permutations.append(permuted_cancer_type_table)
    return pd.concat(permutations)

def get_p_values(wgd_status_cancer_type_table,n_permutations=10000):
    true_diff = get_wgd_agg_values(wgd_status_cancer_type_table)
    true_diff = get_diff(true_diff)
    cancer_type_dict = {cancer_type:wgd_status_cancer_type_table[wgd_status_cancer_type_table['Cancer_Type']==cancer_type] for cancer_type in wgd_status_cancer_type_table['Cancer_Type'].unique()}

    permuted_diff_store = []
    for permutation in range(n_permutations):
        permuted_table = get_permuted_table(cancer_type_dict)
        permuted_agg_values = get_wgd_agg_values(permuted_table)
        permuted_diff = get_diff(permuted_agg_values)
        permuted_diff_store.append(permuted_diff)
    p_value = (np.sum(np.array(permuted_diff_store) >= true_diff) + 1) / (n_permutations + 1)
    return p_value

def plot_proportion_complex_genome(cn_table):
    proportion_complex_genome = cn_table.groupby(['Tumor_Type','WGD_Status','Sample_ID']).apply(lambda x: np.sum(x[x['Major_CN']>2]['Segment_Width'])/np.sum(x['Segment_Width'])).reset_index().rename(columns={0:'Proportion_Complex_Genome'})

    wgd_complex_genome = proportion_complex_genome[proportion_complex_genome['WGD_Status']]
    no_wgd_complex_genome = proportion_complex_genome[~proportion_complex_genome['WGD_Status']]
    #box plot split by WGD status and tumor type
    fig,ax = plt.subplots(figsize=(6,10))
    width=0.8
    ax.boxplot([no_wgd_complex_genome[no_wgd_complex_genome['Tumor_Type']=='Metastatic']['Proportion_Complex_Genome'],wgd_complex_genome[wgd_complex_genome['Tumor_Type']=='Metastatic']['Proportion_Complex_Genome']], positions=np.array([1,4]), widths=width, showfliers=False, patch_artist=True, boxprops=dict(facecolor=metastatic_color, color='black'),medianprops=dict(color='black'))
    ax.boxplot([no_wgd_complex_genome[no_wgd_complex_genome['Tumor_Type']=='Primary']['Proportion_Complex_Genome'],wgd_complex_genome[wgd_complex_genome['Tumor_Type']=='Primary']['Proportion_Complex_Genome'] ], positions=np.array([0,3]), widths=width, showfliers=False, patch_artist=True, boxprops=dict(facecolor=primary_color, color='black'),medianprops=dict(color='black'))
    ax.set_xticks([0.5,3.5])
    ax.set_xticklabels(['Non-WGD\nTumors','WGD Tumors'])
    ax.set_ylabel('Proportion Complex Genome')
    #do mann-whitney test
    u,p = mannwhitneyu(no_wgd_complex_genome[no_wgd_complex_genome['Tumor_Type']=='Primary']['Proportion_Complex_Genome'],no_wgd_complex_genome[no_wgd_complex_genome['Tumor_Type']=='Metastatic']['Proportion_Complex_Genome'])
    #annotate p-values with stars
    star = get_p_value_star(p)

    star_pos = np.percentile(no_wgd_complex_genome[no_wgd_complex_genome['Tumor_Type']=='Metastatic']['Proportion_Complex_Genome'],90)
    ax.annotate(star, xy=(0.2, star_pos), xycoords='data', fontsize=20)

    u,p = mannwhitneyu(wgd_complex_genome[wgd_complex_genome['Tumor_Type']=='Primary']['Proportion_Complex_Genome'],wgd_complex_genome[wgd_complex_genome['Tumor_Type']=='Metastatic']['Proportion_Complex_Genome'])
    #annotate p-values with stars
    star = get_p_value_star(p)
    star_pos = np.percentile(wgd_complex_genome[wgd_complex_genome['Tumor_Type']=='Metastatic']['Proportion_Complex_Genome'],75)

    ax.annotate(star, xy=(3.2, star_pos), xycoords='data', fontsize=20)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.savefig('../plots/cn_plots/complex_genome_boxplot.pdf',dpi=300, bbox_inches='tight')

def get_sufficient_metastatic_primary_counts(wgd_status_cancer_type_table,threshold=20):
    
    cancer_type_counts = wgd_status_cancer_type_table.groupby(['Cancer_Type','Tumor_Type']).size()
    #pivot on tumor type
    cancer_type_counts = cancer_type_counts.reset_index().pivot(index='Cancer_Type',columns='Tumor_Type',values=0).reset_index().fillna(0)
    cancer_type_counts['Metastatic'] = cancer_type_counts['Metastatic'].astype(int)
    cancer_type_counts['Primary'] = cancer_type_counts['Primary'].astype(int)
    sufficient_counts = cancer_type_counts[(cancer_type_counts['Metastatic'] >= threshold) & (cancer_type_counts['Primary'] >= threshold)]
    return sufficient_counts

def plot_wgd_diff(true_wgd_agg,bootstrap_wgd_table,p_value_cancer_type_weight):
    fig,ax = plt.subplots(figsize=(6,10))
    true_wgd_agg = true_wgd_agg.sort_values(by=['Tumor_Type'], ascending=False).reset_index(drop=True)
    bootstrap_wgd_table = bootstrap_wgd_table.sort_values(by=['Tumor_Type'], ascending=False).reset_index(drop=True)
    ax.bar(true_wgd_agg['Tumor_Type'], true_wgd_agg['WGD_Status'], color=[primary_color, metastatic_color])
    #add error bars
    ax.errorbar(true_wgd_agg['Tumor_Type'], true_wgd_agg['WGD_Status'], yerr=[true_wgd_agg['WGD_Status']-bootstrap_wgd_table['CI_Low'], bootstrap_wgd_table['CI_High']-true_wgd_agg['WGD_Status']], fmt='none', ecolor='black', capsize=5)
    #annotate p-value
    ax.text(0.5, 0.95, get_p_value_star(p_value_cancer_type_weight), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    #add labels
    ax.set_ylabel('Cancer Type Weighted Fraction of Samples with WGD')
    #ax.set_xlabel('Tumor Type')
    ax.set_xticklabels(['',''])
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig('../plots/cn_plots/cn_wgd_status_weighted.pdf', bbox_inches='tight')

def plot_wgd_diff_unweighted(wgd_status_table):
    


    wgd_count_table = wgd_status_table.groupby(['Tumor_Type']).agg({'WGD_Status': 'sum', 'Sample_ID': 'count'}).reset_index()
    wgd_count_table = wgd_count_table.sort_values(by=['Tumor_Type'], ascending=False).rename(columns={'WGD_Status': 'WGD_Count', 'Sample_ID': 'Total_Count'}).reset_index(drop=True)
    wgd_count_table['WGD_Count'] = wgd_count_table['WGD_Count'].astype(int)
    wgd_count_table['Total_Count'] = wgd_count_table['Total_Count'].astype(int)

    z,p_value = proportions_ztest(wgd_count_table['WGD_Count'], wgd_count_table['Total_Count'])
    conf_int = proportion_confint(wgd_count_table['WGD_Count'], wgd_count_table['Total_Count'],alpha=0.05)

    fig,ax = plt.subplots(figsize=(6,10))
    ax.bar(wgd_count_table['Tumor_Type'], wgd_count_table['WGD_Count']/wgd_count_table['Total_Count'], color=[primary_color, metastatic_color])
    #add error bars
    ax.errorbar(wgd_count_table['Tumor_Type'], wgd_count_table['WGD_Count']/wgd_count_table['Total_Count'], yerr=[(wgd_count_table['WGD_Count']/wgd_count_table['Total_Count'])-conf_int[0], conf_int[1]-(wgd_count_table['WGD_Count']/wgd_count_table['Total_Count'])], fmt='none', ecolor='black', capsize=5)
    #annotate p-value
    ax.text(0.5, 0.95, get_p_value_star(p_value), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    #add labels
    #ax.set_ylabel('Fraction of Samples with WGD')
    #ax.set_xlabel('Tumor Type')
    #ax.set_xticklabels(['',''])
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6])
    ax.set_xticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig('../plots/cn_plots/cn_wgd_status.pdf', bbox_inches='tight')


def get_diff_p_value(primary_sample_cn_table,metastatic_sample_cn_table):
    p_value_store =[]
    for major_cn_str in ['3','4','5','6','7','8+']:
        primary_major_cn_values = primary_sample_cn_table[primary_sample_cn_table['Major_CN_Str'] == major_cn_str]['Segment_Width_Proportion'].to_numpy()
        metastatic_major_cn_values = metastatic_sample_cn_table[metastatic_sample_cn_table['Major_CN_Str'] == major_cn_str]['Segment_Width_Proportion'].to_numpy()
        true_diff = np.mean(primary_major_cn_values) - np.mean(metastatic_major_cn_values)
        n_permutations = 50000
        diff_store = []
        combined_values = np.concatenate([primary_major_cn_values, metastatic_major_cn_values])
        for i in range(n_permutations):
            permuted_values = RNG.permutation(combined_values)
            permuted_primary_values = permuted_values[:len(primary_major_cn_values)]
            permuted_metastatic_values = permuted_values[len(primary_major_cn_values):]
            diff_store.append(np.mean(permuted_primary_values) - np.mean(permuted_metastatic_values))
        p_value = np.sum(np.abs(diff_store) > np.abs(true_diff))/n_permutations
        p_value_store.append(p_value)
    return np.array(p_value_store)

def plot_cn_proportions(cn_proportions,bootstrap_ci,diff_p_value_store,path,width=0.3,gained_genome=False):
    print('------')
    print(cn_proportions)
    print(gained_genome)
    print(path)
    print('-----')
    fig,ax = plt.subplots(figsize=(10,7))
    ax.bar(np.arange(cn_proportions['Primary'][3:].size)-width, cn_proportions['Primary'][3:], color=primary_color,width=width,align='edge')
    ax.errorbar(np.arange(cn_proportions['Primary'][3:].size)-width/2, cn_proportions['Primary'][3:], yerr=[cn_proportions['Primary'][3:]-bootstrap_ci['Primary'][0,3:], bootstrap_ci['Primary'][1,3:]-cn_proportions['Primary'][3:]], fmt='none', ecolor='black', capsize=5)
    ax.bar(np.arange(cn_proportions['Metastatic'][3:].size), cn_proportions['Metastatic'][3:], color=metastatic_color,width=width,align='edge')
    ax.errorbar(np.arange(cn_proportions['Metastatic'][3:].size)+width/2, cn_proportions['Metastatic'][3:], yerr=[cn_proportions['Metastatic'][3:]-bootstrap_ci['Metastatic'][0,3:], bootstrap_ci['Metastatic'][1,3:]-cn_proportions['Metastatic'][3:]], fmt='none', ecolor='black', capsize=5)
    ax.set_xticks(np.arange(cn_proportions['Primary'][3:].size))
    ax.set_xticklabels(['3','4','5','6','7','8+'])
    if gained_genome:
        ax.set_ylabel('Average Proportion of Gained Genome')
    else:
        ax.set_ylabel('Average Proportion of Genome')
    ax.set_xlabel('Major Copy Number')
    #annotate p-values with stars
    for i in range(cn_proportions['Primary'][3:].size):
        star = get_p_value_star(diff_p_value_store[i])
        ax.annotate(star, xy=(i-width+0.1,cn_proportions['Metastatic'][3+i] ), xycoords='data', fontsize=20)
    ax.set_ylim(0,np.max(cn_proportions['Metastatic'][3:])+0.02)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig(path,dpi=300,bbox_inches='tight')
    plt.close()

def plot_wgd_status_cancer_type(wgd_status_table,classifications):
    wgd_status_table_cancer_type = wgd_status_table.merge(classifications[['Sample_ID','Cancer_Type']], left_on='Sample_ID', right_on='Sample_ID', how='inner')
    cancer_type_counts_tumor_type = wgd_status_table_cancer_type.groupby(['Cancer_Type','Tumor_Type']).agg({'Sample_ID': 'nunique'}).reset_index().pivot(index='Cancer_Type',columns='Tumor_Type',values='Sample_ID').fillna(0).astype(int).reset_index()
    cancer_type_counts_tumor_type = cancer_type_counts_tumor_type[np.minimum(cancer_type_counts_tumor_type['Primary'],cancer_type_counts_tumor_type['Metastatic'])>=10]

    cancer_type_counts = wgd_status_table_cancer_type.groupby(['Cancer_Type']).agg({'Sample_ID': 'nunique'}).reset_index()
    cancer_type_counts = cancer_type_counts.rename(columns={'Sample_ID': 'Sample_Count'})

    cancer_type_counts = cancer_type_counts[cancer_type_counts['Cancer_Type'].isin(cancer_type_counts_tumor_type['Cancer_Type'])]

    cancer_type_counts['Count_Rank'] = cancer_type_counts['Sample_Count'].rank(ascending=False).astype(int)
    cancer_type_counts = cancer_type_counts[cancer_type_counts['Count_Rank'] <= 12]


    fig,axs = plt.subplots(3,4,figsize=(20,15))
    axs = axs.flatten()
    for i in range(len(cancer_type_counts)):
        cancer_type = cancer_type_counts['Cancer_Type'].iloc[i]
        wgd_status_table_cancer_type_sub = wgd_status_table_cancer_type[wgd_status_table_cancer_type['Cancer_Type'] == cancer_type]
        wgd_count_table_cancer_type_sub = wgd_status_table_cancer_type_sub.groupby(['Tumor_Type','WGD_Status']).agg({'Sample_ID': 'nunique'}).reset_index()
        wgd_count_table_cancer_type_sub = wgd_count_table_cancer_type_sub.rename(columns={'Sample_ID': 'Sample_Count'})
        #pivot wider to get WGD status as columns
        wgd_count_table_cancer_type_sub = wgd_count_table_cancer_type_sub.pivot(index='Tumor_Type', columns='WGD_Status', values='Sample_Count').reset_index()
        #add prefix to column names
        wgd_count_table_cancer_type_sub.columns = ['Tumor_Type'] + ['WGD_' + str(col) for col in wgd_count_table_cancer_type_sub.columns[1:]]
        wgd_count_table_cancer_type_sub['Total_Count'] = wgd_count_table_cancer_type_sub['WGD_False'] + wgd_count_table_cancer_type_sub['WGD_True']
        wgd_count_table_cancer_type_sub = wgd_count_table_cancer_type_sub.sort_values(by=['Tumor_Type'], ascending=False).reset_index(drop=True)
        wgd_count_table_cancer_type_sub = wgd_count_table_cancer_type_sub.rename(columns={'WGD_True': 'WGD_Count'})
        try:
            z,p_value = proportions_ztest(wgd_count_table_cancer_type_sub['WGD_Count'], wgd_count_table_cancer_type_sub['Total_Count'])
            conf_int = proportion_confint(wgd_count_table_cancer_type_sub['WGD_Count'], wgd_count_table_cancer_type_sub['Total_Count'],alpha=0.05)
        except:
            print(cancer_type)
            continue
        
        axs[i].bar(wgd_count_table_cancer_type_sub['Tumor_Type'], wgd_count_table_cancer_type_sub['WGD_Count']/wgd_count_table_cancer_type_sub['Total_Count'], color=[primary_color, metastatic_color])
        #add error bars
        axs[i].errorbar(wgd_count_table_cancer_type_sub['Tumor_Type'], wgd_count_table_cancer_type_sub['WGD_Count']/wgd_count_table_cancer_type_sub['Total_Count'], yerr=[(wgd_count_table_cancer_type_sub['WGD_Count']/wgd_count_table_cancer_type_sub['Total_Count'])-conf_int[0], conf_int[1]-(wgd_count_table_cancer_type_sub['WGD_Count']/wgd_count_table_cancer_type_sub['Total_Count'])], fmt='none', ecolor='black', capsize=5)
        #annotate p-value
        max_height =(wgd_count_table_cancer_type_sub['WGD_Count']/wgd_count_table_cancer_type_sub['Total_Count']).max()

        axs[i].text(0.475, max_height+0.05, get_p_value_star(p_value), horizontalalignment='center', verticalalignment='center', transform=axs[i].transAxes)
        #add labels
        axs[i].set_xticklabels(['',''])
        axs[i].set_ylim(0,1)
        axs[i].set_title(cancer_type,pad=15)
        
    plt.tight_layout()
    plt.savefig('../plots/cn_plots/cn_wgd_status_cancertype.pdf', bbox_inches='tight')



def get_cn_proportions_table(cn_table):
    cn_proportions = cn_table.groupby(['Major_CN_Str']).agg({'Segment_Width': 'sum'}).reset_index().rename(columns={'Segment_Width': 'Segment_Width_Sum'})
    cn_proportions['CN_Proportion'] = cn_proportions['Segment_Width_Sum']/cn_proportions['Segment_Width_Sum'].sum()
    return cn_proportions
def get_cn_proportions(cn_table):
    cn_proportions = cn_table.groupby(['Major_CN_Str']).agg({'Segment_Width': 'sum'}).reset_index().rename(columns={'Segment_Width': 'Segment_Width_Sum'})
    print(cn_proportions)
    cn_proportions['CN_Proportion'] = cn_proportions['Segment_Width_Sum']/cn_proportions['Segment_Width_Sum'].sum()
    return cn_proportions.sort_values(by=['Major_CN_Str'], ascending=True)['CN_Proportion'].to_numpy()

def get_cn_bootstrap_ci(cn_table,n_bootstraps=250):
    primary_sample_dict = DataTools.get_sample_dict(cn_table[cn_table['Tumor_Type'] == 'Primary'])
    metastatic_sample_dict = DataTools.get_sample_dict(cn_table[cn_table['Tumor_Type'] == 'Metastatic'])
    
    primary_bootstrap_store = []
    metastatic_bootstrap_store = []
    for i in range(n_bootstraps):
        primary_bootstrap_table = DataTools.get_bootstrap_table(primary_sample_dict)
        metastatic_bootstrap_table = DataTools.get_bootstrap_table(metastatic_sample_dict)
        primary_cn_proportions = get_cn_proportions(primary_bootstrap_table)
        metastatic_cn_proportions = get_cn_proportions(metastatic_bootstrap_table)
        primary_bootstrap_store.append(primary_cn_proportions)
        metastatic_bootstrap_store.append(metastatic_cn_proportions)

    primary_bootstrap_ci = np.percentile(np.array(primary_bootstrap_store), [5/2, 100-5/2], axis=0)
    metastatic_bootstrap_ci = np.percentile(np.array(metastatic_bootstrap_store), [5/2, 100-5/2], axis=0)
    return primary_bootstrap_ci,metastatic_bootstrap_ci

if __name__ == '__main__':
    RNG = np.random.RandomState(42)

    cn_table = load_cn_table()
    wgd_status_table = cn_table[['Sample_ID','WGD_Status','Tumor_Type']].drop_duplicates()

    sample_id_table = cn_table[['Sample_ID']].drop_duplicates()
    classifications = DataTools.load_cancer_type_classifications()

    plot_proportion_complex_genome(cn_table)
    wgd_status_cancer_type_table = cn_table.merge(classifications,how='inner')[['Sample_ID','Cancer_Type','Tumor_Type','WGD_Status']].drop_duplicates()
    
    
    sufficient_counts = get_sufficient_metastatic_primary_counts(wgd_status_cancer_type_table,threshold=20)

    wgd_status_cancer_type_table = wgd_status_cancer_type_table[wgd_status_cancer_type_table['Cancer_Type'].isin(sufficient_counts['Cancer_Type'])]
    
    p_value_cancer_type_weight = get_p_values(wgd_status_cancer_type_table)


    true_wgd_agg = get_wgd_agg_values(wgd_status_cancer_type_table)
    
    bootstrap_wgd_table = get_wgd_bootstrap_table(wgd_status_cancer_type_table, num_bootstrap=250)
    
    plot_wgd_diff(true_wgd_agg,bootstrap_wgd_table,p_value_cancer_type_weight)
    
    plot_wgd_diff_unweighted(wgd_status_table)




    #stats for data export
    cn_table_gained =cn_table[cn_table['Major_CN']>2].copy()
    cn_table_gained['Total_CN_9_or_less'] = (cn_table_gained['Total_CN']<=9) & (cn_table_gained['Major_CN']<=7)
    #segment width
    total_9 = cn_table_gained.groupby('Total_CN_9_or_less').agg({'Segment_Width':'sum'}).reset_index()
    total_9['Segment_Width'] = total_9['Segment_Width']/total_9['Segment_Width'].sum()





    plot_wgd_status_cancer_type(wgd_status_table,classifications)



    true_primary_cn_proportions = get_cn_proportions(cn_table[cn_table['Tumor_Type'] == 'Primary'])
    true_metastatic_cn_proportions = get_cn_proportions(cn_table[cn_table['Tumor_Type'] == 'Metastatic'])
    all_cn_proportions = get_cn_proportions(cn_table)

    true_primary_cn_proportions_table = get_cn_proportions_table(cn_table[cn_table['Tumor_Type'] == 'Primary'])
    true_metastatic_cn_proportions_table = get_cn_proportions_table(cn_table[cn_table['Tumor_Type'] == 'Metastatic'])
    all_cn_proportions_table = get_cn_proportions_table(cn_table)

    primary_bootstrap_ci,metastatic_bootstrap_ci = get_cn_bootstrap_ci(cn_table)

    primary_cn_table = cn_table[cn_table['Tumor_Type'] == 'Primary']
    metastatic_cn_table = cn_table[cn_table['Tumor_Type'] == 'Metastatic']
    primary_sample_cn_table = primary_cn_table.groupby(['Sample_ID','Major_CN_Str']).agg({'Segment_Width': 'sum'}).reset_index().rename(columns={'Segment_Width': 'Segment_Width_Sum'})
    primary_sample_cn_table['Segment_Width_Proportion'] = primary_sample_cn_table.groupby(['Sample_ID'])['Segment_Width_Sum'].apply(lambda x: x/x.sum())
    
    #primary_sample_cn_table = primary_sample_cn_table.merge(primary_sample_cn_table_segment_width_proportion,how='inner')

    metastatic_sample_cn_table = metastatic_cn_table.groupby(['Sample_ID','Major_CN_Str']).agg({'Segment_Width': 'sum'}).reset_index().rename(columns={'Segment_Width': 'Segment_Width_Sum'})
    metastatic_sample_cn_table['Segment_Width_Proportion'] = metastatic_sample_cn_table.groupby(['Sample_ID'])['Segment_Width_Sum'].apply(lambda x: x/x.sum())
    #fill missing values with 0
    primary_sample_cn_table = primary_sample_cn_table.set_index(['Sample_ID','Major_CN_Str']).unstack(level=-1).fillna(0).stack().reset_index()
    metastatic_sample_cn_table = metastatic_sample_cn_table.set_index(['Sample_ID','Major_CN_Str']).unstack(level=-1).fillna(0).stack().reset_index()
    

    print('true_primary_cn_proportions',true_primary_cn_proportions)
    print('true_metastatic_cn_proportions',true_metastatic_cn_proportions)

    diff_p_value_store = get_diff_p_value(primary_sample_cn_table, metastatic_sample_cn_table)
    true_primary_cn_proportions_gained = (true_primary_cn_proportions/sum(true_primary_cn_proportions[2:]))
    true_metastatic_cn_proportions_gained = (true_metastatic_cn_proportions/sum(true_metastatic_cn_proportions[2:]))
    primary_bootstrap_ci_gained = (primary_bootstrap_ci/np.sum(primary_bootstrap_ci[:,2:],axis=1)[:,None])
    metastatic_bootstrap_ci_gained = (metastatic_bootstrap_ci/np.sum(metastatic_bootstrap_ci[:,2:],axis=1)[:,None])

    print('true_primary_cn_proportions_gained',true_primary_cn_proportions_gained)
    print('true_metastatic_cn_proportions_gained',true_metastatic_cn_proportions_gained)

    cn_proportions_store = {'Primary':true_primary_cn_proportions,'Metastatic':true_metastatic_cn_proportions}
    bootstrap_ci_store = {'Primary':primary_bootstrap_ci,'Metastatic':metastatic_bootstrap_ci}

    cn_proportions_gained_store = {'Primary':true_primary_cn_proportions_gained,'Metastatic':true_metastatic_cn_proportions_gained}
    bootstrap_ci_gained_store = {'Primary':primary_bootstrap_ci_gained,'Metastatic':metastatic_bootstrap_ci_gained}
    

    plot_cn_proportions(cn_proportions_store,bootstrap_ci_store,diff_p_value_store,'../plots/cn_plots/major_cn_proportions.pdf')
    plot_cn_proportions(cn_proportions_gained_store,bootstrap_ci_gained_store,diff_p_value_store,'../plots/cn_plots/major_cn_proportions_cn_2_above.pdf',gained_genome=True)







    