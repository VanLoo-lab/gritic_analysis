import numpy as np
import pandas as pd
from numba.typed import List
from numba import njit
import matplotlib.pyplot as plt
def process_cohort(cohort_routes):
    route_probabilities = []
    exponential_likelihoods = []
    true_route_indicies = []
    parsimony_indicies = []
    
    for (sample_id,segment_id),segment_table in cohort_routes.groupby(['Sample_ID','Segment_ID']):
        
        route_prob = segment_table['Probability'].to_numpy()

        route_probabilities.append(route_prob)
        exponential_likelihood = np.exp(-segment_table['Average_N_Events'].to_numpy())
        
        exponential_likelihood[np.isnan(exponential_likelihood)] = 0
        exponential_likelihoods.append(exponential_likelihood)
        true_route_indicies.append(np.where(segment_table['Route'].isin(segment_table['True_Route']))[0])
        parsimony_indicies.append(np.where(segment_table['Parsimonious'])[0])

    return List(route_probabilities),List(exponential_likelihoods),List(true_route_indicies),List(parsimony_indicies)
@njit
def evaluate_l(route_probabilities,exponential_likelihoods,true_route_indicies,parsimony_indicies,l):
    n_correct = 0
    parsimony_prob = 0
    for i in range(len(route_probabilities)):
        prob_state = np.multiply(route_probabilities[i],np.power(exponential_likelihoods[i],l))

        prob_state = prob_state/np.sum(prob_state)
        if np.isnan(prob_state).any():
            continue
       
        max_prob_index = np.argmax(prob_state)
        for index in true_route_indicies[i]:
            if index == max_prob_index:
                n_correct+=1
        for index in parsimony_indicies[i]:
            parsimony_prob += prob_state[index]

    return n_correct,parsimony_prob

if __name__ == '__main__':
    out_dir = '../output/evaluate_parsimony_penalty_out'
    dataset = 'state_validation_parsimony_only'
    cohort_metadata = pd.read_csv(f'../../output/{dataset}/complete_metadata_table_{dataset}.tsv',sep="\t").rename(columns={'Route':'True_Route'}).drop(columns=['Average_Pre_WGD_Losses','Average_Post_WGD_Losses','N_Mutations'])
    cohort_metadata['True_Route'] = cohort_metadata['True_Route'].str.slice(0,9)
    cohort_metadata['Segment_ID'] = cohort_metadata['Segment_ID'].str.replace('Segment_','')

    control_dataset = 'state_validation_parsimony_only'
    cohort_routes = pd.read_csv(f'../../output/{control_dataset}/complete_route_table_{control_dataset}.tsv',sep='\t',dtype={'Sample_ID':str,'Segment_ID':str,'Chromosome':str})
    cohort_routes = cohort_routes[cohort_routes['WGD_Status']]
    cohort_routes = cohort_routes[cohort_routes['Major_CN']>2]

    cohort_routes = cohort_routes[cohort_routes['N_Mutations']>=20]


    cohort_routes = cohort_routes[['Sample_ID','Segment_ID','Segment_Start','Segment_End','Route','N_Mutations','Average_Pre_WGD_Losses', 'Average_Post_WGD_Losses','Average_N_Events','Major_CN','Minor_CN','Probability']].drop_duplicates()


    cohort_routes['Average_Losses'] = cohort_routes['Average_Pre_WGD_Losses'] + cohort_routes['Average_Post_WGD_Losses']

    cohort_routes = cohort_routes.merge(cohort_metadata[['Sample_ID','Segment_ID','True_Route']],how='inner')

    parsimony_table = pd.read_csv('../../Analysis/output/all_routes_n_events.tsv',sep="\t")
    parsimony_table = parsimony_table[parsimony_table['WGD_Status']]
    parsimony_table['Average_Losses'] = parsimony_table['Average_Pre_WGD_Losses'] + parsimony_table['Average_Post_WGD_Losses']
    parsimony_table['Route'] = parsimony_table['Route'].str.slice(0,9)

    cohort_routes = cohort_routes.merge(parsimony_table[['Major_CN','Minor_CN','Route','Parsimonious']].drop_duplicates(),how='inner')

    l_values = np.linspace(0,5,501)
    prop_store = {}
    prop_tables = []
    for (major_cn,minor_cn),cohort_cn in cohort_routes.groupby(['Major_CN','Minor_CN']):

        prop_data_cn = {'Major_CN':major_cn,'Minor_CN':minor_cn,'N_Accurate':[],'Parsimony_Prob':[],'L_Values':l_values}
        if major_cn not in prop_store:
            prop_store[major_cn] = {}
        if minor_cn not in prop_store[major_cn]:
            prop_store[major_cn][minor_cn] = {}
        
        route_probabilities,exponential_likelihoods,true_route_indicies,parsimony_indicies = process_cohort(cohort_cn)
        prop_data_cn['N_Segments'] = len(route_probabilities)
        for l in l_values:
            n_accurate,parsimony_prob = evaluate_l(route_probabilities,exponential_likelihoods,true_route_indicies,parsimony_indicies,l)
            prop_data_cn['N_Accurate'].append(n_accurate)
            prop_data_cn['Parsimony_Prob'].append(parsimony_prob)
        prop_tables.append(pd.DataFrame(prop_data_cn))
    prop_table = pd.concat(prop_tables)

    prop_table.to_csv('../output/evaluate_parsimony_penalty_out/prop_table.tsv',sep='\t',index=False)


    agg_table = prop_table.groupby(['L_Values']).agg({'N_Segments':'sum','Parsimony_Prob':'sum','N_Accurate':'sum'}).reset_index()
    agg_table['Prop_Accuracy'] = agg_table['N_Accurate']/agg_table['N_Segments']
    agg_table['Parsimony_Prob'] = agg_table['Parsimony_Prob']/agg_table['N_Segments']
    agg_table['Non_Parsimony_Prob'] = 1 - agg_table['Parsimony_Prob']

    agg_table = agg_table.sort_values(by='Parsimony_Prob', key=lambda col: np.abs(col-0.95))

    agg_table.to_csv('../output/evaluate_parsimony_penalty_out/agg_table.tsv',sep='\t',index=False)

