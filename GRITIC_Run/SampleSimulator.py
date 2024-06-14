import numpy as np
import pandas as pd

import string
from gritic import gritictimer


RNG = np.random.default_rng()

def generate_random_sample_id(length):
    id_list = list(string.ascii_letters+string.digits)
    sample_id = ''.join(RNG.choice(id_list,size=length,replace=True))
    return sample_id


def sample_mutation_multiplicities(multiplicity_proportions,n_mutations,subclone_table):
    multiplicity_states = np.arange(multiplicity_proportions.size)+1
    if subclone_table is not None:
        clonal_share = 1-subclone_table['Subclone_Fraction'].sum()
        multiplicity_proportions = np.concatenate([multiplicity_proportions*clonal_share,subclone_table['Subclone_Fraction'].to_numpy()])
        multiplicity_states = np.concatenate([multiplicity_states,subclone_table['Subclone_CCF'].to_numpy()])
    
    multiplicity_proportions = np.clip(multiplicity_proportions,0,1)
    mutation_samples = RNG.multinomial(n_mutations,multiplicity_proportions)
    multiplicities = []
    for index in range(multiplicity_proportions.size):
        multiplicity = multiplicity_states[index]
        multiplicities.extend([multiplicity]*mutation_samples[index])
    return np.array(multiplicities)

def sample_read_counts(mutation_multiplicities,average_coverage,sample_purity,total_cn,normal_cn):
    mult_one_vaf = sample_purity/(total_cn*sample_purity+normal_cn*(1-sample_purity))
    vafs = mutation_multiplicities*mult_one_vaf
    coverages =RNG.poisson(average_coverage,size=vafs.size)
    alt_counts = RNG.binomial(n=coverages,p=vafs)
    ref_counts = coverages - alt_counts
    return np.vstack([ref_counts,alt_counts]),vafs

def run_single_route_sample(classifier,cn_row):
    route = classifier.routes[cn_row['Route']]
    mults,timing = route.sample_mults(cn_row['WGD_Timing'],250)
    
    random_index = RNG.integers(0,mults.shape[0]-1)
    mult_sample = mults[random_index]
    
    timing_sample = timing[:,random_index]
    node_gain_timing = {}
    for index,node in enumerate(route.route_tree.non_phased_node_order):
        if node in route.route_tree.timeable_nodes:
            node_gain_timing[node] = timing_sample[index]
    return mult_sample,node_gain_timing

def run_route_sample_non_wgd(classifier,cn_row):
    route = classifier.routes[cn_row['Route']]
    mults,timing = route.sample_mults(None,250)
    
    random_index = RNG.integers(0,mults.shape[0]-1)
    mult_sample = mults[random_index]
    
    timing_sample = timing[:,random_index]
    node_gain_timing = {}
    for index,node in enumerate(route.route_tree.node_order):
        if node in route.route_tree.timeable_nodes:
            node_gain_timing[node] = timing_sample[index]
    return mult_sample,node_gain_timing

def get_cn_timing_distance(timing_target,timing_sim):
    #TODO implement a function that can handle timing_target being bigger than timing_sim

    assert timing_target.size <= timing_sim.size
    while timing_target.size >0:
        distances = np.abs(timing_sim-timing_target.reshape(-1,1))
        min_distance = np.min(distances)
        
        
        timing_target = np.delete(timing_target,np.argmin(distances,axis=0))
        timing_sim = np.delete(timing_sim,np.argmin(distances,axis=1))
    
    return min_distance
def run_route_time_sampling(classifier,cn_row,target=True,n_samples=500,distance_threshold=0.005):
    for i in range(n_samples):
        mults,timing = run_single_route_sample(classifier,cn_row)
        if not target:
            return mults,timing
        if cn_row['Gain_Timing_Target'].size ==0:
            return mults,timing
        else:
            distance= get_cn_timing_distance(cn_row['Gain_Timing_Target'],np.array(list(timing.values())))
            if distance < distance_threshold:
                return mults,timing
    return None,None


def get_major_mut_share(major_mult_sample,minor_mult_sample,major_cn,minor_cn,subclone_table):
    if minor_cn ==0:
        return 1
    major_array = np.arange(1,major_cn+1)
    minor_array= np.arange(1,minor_cn+1)
    if subclone_table is not None:
        n_subclones = len(subclone_table)
        subclone_fraction = subclone_table['Subclone_Fraction'].values
        clonal_fraction = 1 - np.sum(subclone_fraction)

        major_mult_sample = np.concatenate([subclone_fraction,major_mult_sample*clonal_fraction])
        minor_mult_sample = np.concatenate([subclone_fraction,minor_mult_sample*clonal_fraction])
        major_array = np.concatenate([np.ones(n_subclones),major_array])
        minor_array = np.concatenate([np.ones(n_subclones),minor_array])

 

    normed_major_mult_sample = major_cn/np.sum(np.multiply(major_array,major_mult_sample))
    normed_minor_mult_sample = minor_cn/np.sum(np.multiply(minor_array,minor_mult_sample))

    major_mut_share = normed_major_mult_sample/(normed_major_mult_sample+normed_minor_mult_sample)
    return major_mut_share
def simulate_sample(sample_cn_table,subclone_table,wgd_status,sex):
    sample_cn_table = sample_cn_table.copy()
    sample_cn_table['Index'] = np.arange(len(sample_cn_table.index))
    if wgd_status:
        wgd_trees_status = 'Default'
    else:
        wgd_trees_status = 'No_WGD'

    segment_frame_store = []
    segment_metadata_store = {'Index':[],'Segment_ID':[],'Gain_Timing':[]}
    
    segment_count = 0
    unable_count = 0
    for index,cn_row in sample_cn_table.iterrows():
       
        major_cn = cn_row['Major_CN']
        minor_cn = cn_row['Minor_CN']
        sample_id =cn_row['Sample_ID']
        segment_id = cn_row['Segment_ID']
        classifier = gritictimer.RouteClassifier(major_cn,minor_cn,wgd_status,wgd_trees_status,mult_store_dir=None)

        
        segment_id = f'Segment_{segment_id}'
        if (wgd_status and major_cn >2) or (not wgd_status and major_cn >1):
            mult_sample,node_gain_timing = run_route_time_sampling(classifier,cn_row)
            major_mult_sample = mult_sample[major_cn:major_cn*2]
            minor_mult_sample = mult_sample[major_cn*2:]
        
            if mult_sample is None:
                unable_count +=1
                continue
        elif (major_cn == 2 and wgd_status):
            
            major_mult_sample = np.array([1-cn_row['WGD_Timing']/(2-cn_row['WGD_Timing']),cn_row['WGD_Timing']/(2-cn_row['WGD_Timing'])])
            if minor_cn ==1:
                minor_mult_sample = np.array([1])
            elif minor_cn ==0:
                minor_mult_sample = np.array([])
            elif minor_cn ==2:
                minor_mult_sample = major_mult_sample
            node_gain_timing = {0:cn_row['WGD_Timing']}
        else:
            raise ValueError('Major CN must be greater than 1 for non-WGD samples')
        
        #this is clonal only
        major_mut_share = get_major_mut_share(major_mult_sample,minor_mult_sample,major_cn,minor_cn,subclone_table)
        
        n_major_mutations = np.maximum(RNG.binomial(cn_row['N_Mutations'],major_mut_share),1)
        n_minor_mutations = cn_row['N_Mutations']-n_major_mutations
        major_mutation_multiplicities = sample_mutation_multiplicities(major_mult_sample,n_major_mutations,subclone_table)
        
        if n_minor_mutations >0:
            minor_mutation_multiplicities = sample_mutation_multiplicities(minor_mult_sample,n_minor_mutations,subclone_table)
        else:
            minor_mutation_multiplicities = np.array([])
        
        mutation_multiplicities = np.concatenate([major_mutation_multiplicities,minor_mutation_multiplicities])

        mutation_phasing = ['Major']*n_major_mutations+['Minor']*n_minor_mutations
        normal_cn = 1 if (sex=='XY' and cn_row['Chromosome']=='X') else 2
        mutation_read_counts,vafs = sample_read_counts(mutation_multiplicities,cn_row['Adjusted_Coverage'],cn_row['Purity'],major_cn+minor_cn,normal_cn)
        
        segment_frame_data = {'Sample_ID':sample_id,'Segment_ID':segment_id,'Major_CN':major_cn,'Minor_CN':minor_cn,'Tumor_Ref_Count':list(mutation_read_counts[0,:]),'Tumor_Alt_Count':list(mutation_read_counts[1,:]),'True_VAF':vafs,'True_Multiplicity':mutation_multiplicities,'Phasing':np.nan,'Chromosome':cn_row['Chromosome'],'Segment_Start':cn_row['Segment_Start'],'Segment_End':cn_row['Segment_End']}
        

        segment_frame = pd.DataFrame(segment_frame_data)

        
        segment_frame['Total_CN'] = major_cn +minor_cn
        segment_frame['Position'] = np.arange(cn_row['Segment_Start']+1,cn_row['Segment_Start']+cn_row['N_Mutations']+1)
        segment_frame['Phasing'] = mutation_phasing
        
        
        #segment_frame = segment_frame[segment_frame['Tumor_Alt_Count']>=3].reset_index(drop=True)
        segment_frame_store.append(segment_frame)
        segment_count +=1
        if major_cn>2:
            segment_metadata_store['Segment_ID'].append(segment_id)
            segment_metadata_store['Gain_Timing'].append(str(node_gain_timing))
            segment_metadata_store['Index'].append(cn_row['Index'])


    mutation_table = pd.concat(segment_frame_store)

    #no phasing where major and minor CN are the same
    mutation_table['Phasing'] = np.where(mutation_table['Major_CN']==mutation_table['Minor_CN'],pd.NA,mutation_table['Phasing'])
    #only a third of SNVs are phased
    #mutation_table['Phasing'] = np.where(RNG.random(mutation_table.shape[0])<2/3,pd.NA,mutation_table['Phasing'])
    mutation_table = mutation_table.drop(columns=['Segment_Start','Segment_End'])
    segment_metadata = pd.DataFrame(segment_metadata_store)
    
    segment_metadata = segment_metadata.merge(sample_cn_table.drop(columns=['Segment_ID']),how='inner',on=['Index']).drop(columns=['Index'])
    return mutation_table,segment_metadata



