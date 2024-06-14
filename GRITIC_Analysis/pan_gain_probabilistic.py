import pandas as pd
import numpy as np

import DataTools

import sys
import plot_wgd_sampling

import matplotlib.pyplot as plt
RNG = np.random.default_rng()

def get_bins_even_spacing(timing_samples,bin_width=0.1):
    bin_lower_target = np.percentile(timing_samples,0.5)
    bin_upper_target = np.percentile(timing_samples,99.5)
    potential_bins = np.linspace(-1.05,1.05,22)
    bin_lower_index = np.abs(potential_bins-bin_lower_target).argmin()
    bin_upper_index = np.abs(potential_bins-bin_upper_target).argmin()

    if potential_bins[bin_lower_index] < bin_lower_target:
        bin_lower_index += 1
    if potential_bins[bin_upper_index] > bin_upper_target:
        bin_upper_index -= 1
    bins = np.arange(potential_bins[bin_lower_index],potential_bins[bin_upper_index]+bin_width,bin_width)
    
    return bins

def get_bins(timing_samples,n_bins,bin_method,target_width=0.95):
    if bin_method =='Even_Spacing':
        return get_bins_even_spacing(timing_samples)
    bin_lower_target = np.percentile(timing_samples,0.5)
    bin_upper_target = np.percentile(timing_samples,99.5)

    width_offset = (target_width-(bin_upper_target-bin_lower_target))/2
    bin_lower_target -= width_offset
    bin_upper_target += width_offset
    bins = np.linspace(bin_lower_target,bin_upper_target,n_bins+1)
    if bin_method =='Original':
        bins -= bins[np.argmin(np.abs(bins))]
    
    return bins

def get_background_hist(posterior_table,bins,n_runs= 20):
    full_hist = np.zeros(bins.size-1)
        
    posterior_table_run = posterior_table.copy()
    posterior_table_sample_width = posterior_table_run.groupby(['Sample_ID']).agg({'Segment_Width':'sum'}).reset_index()
    posterior_table_sample_width = posterior_table_sample_width.rename(columns={'Segment_Width':'Sample_Width'})
    posterior_table_run = posterior_table_run.merge(posterior_table_sample_width,on='Sample_ID')
    posterior_table_run['Sample_Weighting'] = posterior_table_run['Segment_Width']/posterior_table_run['Sample_Width']
    
    wgd_timing = posterior_table_run['WGD_Timing'].values
    sample_weighting = posterior_table_run['Sample_Weighting'].values
    cancer_type_norm = posterior_table_run['Cancer_Type_Norm'].values
    timing_weights = np.multiply(sample_weighting,cancer_type_norm)
    for i in range(n_runs):
        random_timing = RNG.uniform(0,1,posterior_table.shape[0])
        random_hist = np.histogram(random_timing-wgd_timing,bins=bins,weights=timing_weights)[0]
        random_hist = random_hist/np.sum(random_hist)
        full_hist += random_hist
    full_hist = full_hist/n_runs
    return full_hist
def get_posterior_hist(posterior_table,bins=None,bin_method=None,n_bins=10):
    timing_samples = posterior_table['Gain_Timing']-posterior_table['WGD_Timing']

    if bins is None:
        bins = get_bins(timing_samples,n_bins,bin_method)
    
    background_hist = get_background_hist(posterior_table,bins)
   
    posterior_hist,bins = np.histogram(posterior_table['Gain_Timing']-posterior_table['WGD_Timing'],bins=bins,weights=posterior_table['Segment_Width'].values*posterior_table['Cancer_Type_Norm'].values)
    posterior_hist,bins = np.histogram(posterior_table['Gain_Timing']-posterior_table['WGD_Timing'],bins=bins)
    posterior_hist = posterior_hist/np.sum(posterior_hist)

    return posterior_hist/background_hist,bins


def run_wgd_timing_histograms(posterior_table_wgd,wgd_bins,hist_bins=None,bin_method=None):
    
    if bin_method is not None and hist_bins is not None:
        raise ValueError('Can only specify one of hist_bins and bin_method')
    sample_wgd_summary = posterior_table_wgd.groupby(['Sample_ID']).agg({'WGD_Timing':'mean'}).reset_index()

    hist_store = []
    bin_store = []
    n_samples_store = []

    for index in range(wgd_bins.size-1):
        wgd_bin_low, wgd_bin_high = wgd_bins[index], wgd_bins[index+1]
        wgd_bin_samples = sample_wgd_summary[(sample_wgd_summary['WGD_Timing']>=wgd_bin_low)&(sample_wgd_summary['WGD_Timing']<wgd_bin_high)]['Sample_ID'].values

        wgd_bin_posterior_table = posterior_table_wgd[posterior_table_wgd['Sample_ID'].isin(wgd_bin_samples)].copy()
   

        n_samples_store.append(len(wgd_bin_posterior_table['Sample_ID'].unique()))
        if hist_bins is not None:
            bins =  hist_bins[index]
        else:
            bins = None
        hist,bins = get_posterior_hist(wgd_bin_posterior_table,bins=bins,bin_method=bin_method)
        hist_store.append(hist)
        bin_store.append(bins)

    return np.array(hist_store,dtype=object),np.array(bin_store,dtype=object),np.array(n_samples_store)
def get_run_types():
    '''cancer_types = ['All','All_Balanced','Ovary','Kidney_ClearCell','Breast','Pancreas','Gastric','CNS_Glioma'
    ,'CNS_Medullo','Sarcoma_Leiomyo','Colorectum','Skin_Melanoma'
    ,'Kidney_Chromophobe','HeadAndNeck_Other','Lung_NSC','Liver','XXX'
    ,'Prostate','Uterus','Urothelial','Pancreas_NET','Cervix','Sarcoma_Lipo'
    ,'Biliary','Lymphoid','Kidney_Papillary','Sarcoma_Other','Lung_SC'
    ,'Sarcoma_GIST','Unknown','Gastrointestinal_NET','Thyroid','Vulva'
    ,'Mesothelium','Skin_Other','Anus','Adrenal','Testis','Penis'
    ,'HeadAndNeck_SG','Vagina','Thymus_NET','Kidney_Other','CNS_Other','Eye'
    ,'Thymus','Lung_NET']'''

    classifications = DataTools.load_classifications()
    cancer_type_counts = classifications['Cancer_Type'].value_counts().reset_index().rename(columns={'Cancer_Type':'Count'})
    cancer_type_counts = cancer_type_counts[cancer_type_counts['Count']>=100]
    cancer_types = list(cancer_type_counts['index'])
    cancer_types.extend(['All','All_Balanced','Random','Simulated'])
    min_mutations = [0,100,250]
    bin_methods = ['Original','No_Offset','Even_Spacing']

    tumor_types = ['All']
    run_types = []
    for cancer_type in cancer_types:
        for tumor_type in tumor_types:
            for apply_penalty in [True,False]:
                for min_mutation in min_mutations:
                    for bin_method in bin_methods:
                        for n_bins in [3,4]:
                            run_types.append((cancer_type,tumor_type,apply_penalty,min_mutation,bin_method,n_bins))
                            
    return run_types


def load_metadata_posterior_table(good_df):
    if good_df:
        metadata_path = '../../output/state_validation_all_routes/complete_metadata_table_state_validation_all_routes.tsv'
    else:
        metadata_path = '../../output/state_validation_uniform_sampling_new11/combined_metadata.tsv'
    metadata = pd.read_csv(metadata_path,sep='\t',dtype={'Chromosome':str})
    add_cols = ['Sample_ID','Segment_ID','Major_CN','Minor_CN','Chromosome','Segment_Start','Segment_End','WGD_Timing','N_Mutations','Route','Dataset_ID']
    new_posterior = {col:[] for col in add_cols}
    new_posterior['Gain_Timing'] = []
    new_posterior['Gain_Timing_True'] = []
    
    for index,row in metadata.iterrows():
        gain_timing = np.sort(np.array(list(eval(row['Gain_Timing']).values())))
        new_posterior['Gain_Timing'].extend(gain_timing)
        for col in add_cols:
            new_posterior[col].extend([row[col]]*len(gain_timing))
        new_posterior['Gain_Timing_True'].extend([str(row['Gain_Timing'])]*len(gain_timing))
    metadata = pd.DataFrame(new_posterior)
    metadata['Posterior_Sample_Index'] = 0
    metadata['Dataset'] = 'Simulated'
    metadata['Cancer_Type'] = 'Simulated'
    metadata['Segment_Width'] = metadata['Segment_End'] - metadata['Segment_Start']
   

    
    return metadata





def load_simulated_posterior_table(apply_penalty):

    posterior_path = f'../../output/state_validation_uniform_sampling/complete_posterior_table_penalty_{apply_penalty}_state_validation_uniform_sampling.tsv'
    
    posterior_table = pd.read_csv(posterior_path,sep='\t',dtype={'Chromosome':str})
    
    wgd_table = posterior_table.copy()
    wgd_table['WGD_Status'] =~wgd_table['WGD_Timing'].isnull()
    wgd_table = wgd_table[['Sample_ID','WGD_Status']].drop_duplicates()
    
    posterior_table = posterior_table[~posterior_table['WGD_Timing'].isnull()]
    posterior_table['Segment_Width'] = posterior_table['Segment_End'] - posterior_table['Segment_Start']

    return posterior_table
def get_cancer_type_norm(posterior_table):
    cancer_type_count = posterior_table.groupby(['Cancer_Type']).agg({'Sample_ID':'nunique'}).reset_index()
    cancer_type_count = cancer_type_count.rename(columns={'Sample_ID':'N_Samples'})
    cancer_type_count.loc[:,'Cancer_Type_Norm'] = 1/cancer_type_count['N_Samples']
    return cancer_type_count[['Cancer_Type','Cancer_Type_Norm']]



if __name__ == '__main__':
    run_number = int(sys.argv[1])
    run_types = get_run_types()
    
    
        
    cancer_type,tumor_type,apply_penalty,min_mutations,bin_method,n_bins = run_types[run_number]

    out_dir = '../output/pan_wgd_sampling'
    base_out_path = f"{out_dir}/{cancer_type.replace(' ','_')}_{tumor_type.replace(' ','_')}_apply_penalty_{apply_penalty}_min_mutations_{min_mutations}_bin_method_{bin_method}_n_bins_{n_bins}"
    print(cancer_type,tumor_type,apply_penalty,min_mutations,bin_method,n_bins)
    wgd_bins = np.linspace(0,1,n_bins+1)
    n_bootstraps = 250
    if cancer_type == 'Simulated':
        posterior_table_wgd = load_simulated_posterior_table(apply_penalty=apply_penalty)
        print(posterior_table_wgd)
    else:
        posterior_table = DataTools.load_posterior_table(apply_penalty=apply_penalty)
        posterior_table_wgd = posterior_table[~posterior_table['WGD_Timing'].isnull()].copy()
        if not cancer_type in ['All','All_Balanced','Random']:
            posterior_table_wgd = posterior_table_wgd[posterior_table_wgd['Cancer_Type']==cancer_type]
        if cancer_type == 'All_Balanced':
            cancer_type_counts = posterior_table_wgd[['Sample_ID','Cancer_Type']].drop_duplicates()['Cancer_Type'].value_counts()
            cancer_type_counts = cancer_type_counts[cancer_type_counts>=100]
            print(cancer_type_counts)
            posterior_table_wgd = posterior_table_wgd[posterior_table_wgd['Cancer_Type'].isin(cancer_type_counts.index)]
        if cancer_type == 'Random':
            posterior_table_wgd.loc[:,'Gain_Timing'] = np.random.uniform(0,1,len(posterior_table_wgd))
        if tumor_type =='Primary':
            posterior_table_wgd = posterior_table_wgd[posterior_table_wgd['Dataset']=='PCAWG']
        elif tumor_type =='Metastasis':
            posterior_table_wgd = posterior_table_wgd[posterior_table_wgd['Dataset']=='Hartwig']

        print(posterior_table_wgd)
        if len(posterior_table_wgd['Sample_ID'].unique())<10:
            print(f'Not enough samples for {cancer_type} {tumor_type}')
            exit()
    if cancer_type == 'All_Balanced':
        cancer_type_norm = get_cancer_type_norm(posterior_table_wgd)
        posterior_table_wgd = pd.merge(posterior_table_wgd,cancer_type_norm,on='Cancer_Type')
    else:
        posterior_table_wgd.loc[:,'Cancer_Type_Norm'] = 1
    

    posterior_table_wgd = posterior_table_wgd[posterior_table_wgd['N_Mutations']>=min_mutations]


    seg_table = posterior_table_wgd.groupby(['Sample_ID','Segment_ID'])['WGD_Timing'].median().reset_index()

    for i in range(len(wgd_bins)-1):
        wgd_bin_low = wgd_bins[i]
        wgd_bin_high = wgd_bins[i+1]
        posterior_table_wgd_mid = posterior_table_wgd[(posterior_table_wgd['WGD_Timing']>=wgd_bin_low) & (posterior_table_wgd['WGD_Timing']<wgd_bin_high)]
        if len(posterior_table_wgd_mid['Sample_ID'].unique())<10:
            print(f'Not enough samples for {cancer_type} {tumor_type} {wgd_bin_low} {wgd_bin_high}')
            exit()
    true_hists,true_bins,n_samples_true = run_wgd_timing_histograms(posterior_table_wgd,wgd_bins,bin_method=bin_method)
    
    sample_table_dict = DataTools.get_sample_dict(posterior_table_wgd)
    bootstrap_hist_store = []
    for i in range(n_bootstraps):
        bootstrap_posterior_table_wgd = DataTools.get_bootstrap_table(sample_table_dict)
        bootstrap_hists,bootstrap_bins,_= run_wgd_timing_histograms(bootstrap_posterior_table_wgd,wgd_bins,hist_bins=true_bins)

        bootstrap_hist_store.append(bootstrap_hists)

  
    np.save(f'{base_out_path}_true_hists.npy',true_hists)
    np.save(f'{base_out_path}_true_bins.npy',true_bins)
    np.save(f'{base_out_path}_n_samples_true.npy',n_samples_true)
    
    bootstrap_hist_store = np.array(bootstrap_hist_store)
    np.save(f'{base_out_path}_bootstrap_hists.npy',bootstrap_hist_store)
    
    plot_wgd_sampling.plot_bootstrap_hists(cancer_type,tumor_type,apply_penalty,min_mutations,wgd_bins,bin_method,n_bins)
    

        
        
    
