import pandas as pd
import numpy as np
import os
import plotly.graph_objects as go

import DataTools
import copy
import traceback
import warnings
import sys
import time

def get_combined_absolute_segments(segment_table):
    absolute_segments = []
    for chromosome,chromosome_table in segment_table.groupby('Chromosome'):
        chromosome_break_data = {'Chromosome':chromosome}
        all_segment_breaks = np.concatenate([chromosome_table['Segment_Abs_Start'].to_numpy(),chromosome_table['Segment_Abs_End'].to_numpy()])
        all_segment_breaks = np.append(all_segment_breaks,[chromosome_table['Chrom_Abs_Start'].iloc[0],chromosome_table['Chrom_Abs_End'].iloc[0]])
        all_segment_breaks = np.sort(np.unique(all_segment_breaks))

        chromosome_break_data['Abs_Seg_Start'] = all_segment_breaks[:-1]
        chromosome_break_data['Abs_Seg_End'] = all_segment_breaks[1:]
        absolute_segments.append(pd.DataFrame(chromosome_break_data))

    absolute_segments = pd.concat(absolute_segments)
    absolute_segments = absolute_segments.sort_values(by=['Abs_Seg_Start']).reset_index()
    return absolute_segments
 
def get_absolute_positions(segment_table):
    segment_table = segment_table.copy()
    chromosome_table =load_chromosome_table()
    segment_table = segment_table.merge(chromosome_table)
    
    segment_table.loc[:,'Segment_End'] = np.minimum(segment_table['Segment_End'],segment_table['Chrom_End'])
    segment_table.loc[:,'Segment_Abs_Start'] = segment_table['Segment_Start']+segment_table['Chrom_Abs_Start']
    segment_table.loc[:,'Segment_Abs_End'] = segment_table['Segment_End']+segment_table['Chrom_Abs_Start'] 
    return segment_table

def get_proportion_line(segment_table,combined_absolute_segments,event_type=None,n_samples=None,weight_col=None):
    if n_samples is None:
        n_samples = len(segment_table['Sample_ID'].unique())
    if event_type is not None:
        event_table= segment_table[segment_table[event_type]]
    else:
        event_table = segment_table

    segment_positions = []
    segment_proportions = []
    event_table_chr = {chr:chr_table for chr,chr_table in event_table.groupby('Chromosome')}

    for index,segment_row in combined_absolute_segments.iterrows():
        matching_event_table = event_table_chr[segment_row['Chromosome']]
        intersecting_rows = matching_event_table[(matching_event_table['Segment_Abs_Start']<=segment_row['Abs_Seg_Start']) &(matching_event_table['Segment_Abs_End']>segment_row['Abs_Seg_End'])]
        #intersecting_rows = event_table[(event_table['Segment_Abs_Start']<=segment_row['Abs_Seg_End']) &(event_table['Segment_Abs_End']>segment_row['Abs_Seg_Start'])]
        n_intersecting_samples = len(intersecting_rows['Sample_ID'].unique())
        
        if weight_col is not None:
            proportion_intersecting = intersecting_rows[weight_col].sum()/n_samples
        else:
            proportion_intersecting = n_intersecting_samples/n_samples

        segment_positions.extend([segment_row['Abs_Seg_Start'],segment_row['Abs_Seg_End']])
        segment_proportions.extend([proportion_intersecting,proportion_intersecting])
    return segment_positions,segment_proportions

def get_event_proportions(cancer_type_table,event_type=None,n_samples=None,weight_col =None):
    cancer_type_table = get_absolute_positions(cancer_type_table)

    combined_absolute_segments = get_combined_absolute_segments(cancer_type_table)
    
    proportion_positions,proportions = get_proportion_line(cancer_type_table,combined_absolute_segments,event_type,n_samples=n_samples,weight_col=weight_col)
    return {'Positions':np.array(proportion_positions),'Proportions':np.array(proportions)}

def load_chromosome_arms():
    cytoband_path = '../resources/hg19_cytoBand.txt'
    cytoband_table = pd.read_csv(cytoband_path,sep="\t",names=['chrom','chromStart','chromEnd','name','gieStain'])
    cytoband_table = cytoband_table[~cytoband_table['chrom'].str.contains('_')]
    cytoband_table = cytoband_table.dropna()
    cytoband_table.loc[:,'Arm'] = np.nan
    cytoband_table.loc[cytoband_table['name'].str.contains('p'),'Arm']='p'
    cytoband_table.loc[cytoband_table['name'].str.contains('q'),'Arm']='q'
    cytoband_table = cytoband_table.dropna()

    chromosome_arm_data = {'Chromosome_Arm':[],'Chromosome':[],'Arm':[],'Arm_Start':[],"Arm_End":[]}
    
    for group_vars,group_data in cytoband_table.groupby(['chrom','Arm']):
        chromosome,arm = group_vars
        chromosome_arm = chromosome.replace("chr","")+arm
        chromosome_arm_data['Chromosome_Arm'].append(chromosome_arm)
        chromosome_arm_data['Chromosome'].append(chromosome.replace('chr',''))
        chromosome_arm_data['Arm'].append(arm)
        chromosome_arm_data['Arm_Start'].append(group_data['chromStart'].min())
        chromosome_arm_data['Arm_End'].append(group_data['chromEnd'].max())
    chromosome_arm_table = pd.DataFrame(chromosome_arm_data)
    return chromosome_arm_table

def load_chromosome_table():
    chromosome_order = list(map(str,range(1,23)))
    chromosome_order.extend(['X','Y'])
    chromosome_arms = load_chromosome_arms()
    chromosome_table =chromosome_arms.groupby('Chromosome').agg({'Arm_Start':'min','Arm_End':'max'}).reset_index()
    chromosome_table = chromosome_table.rename(columns={'Arm_Start':'Chrom_Start','Arm_End':'Chrom_End'})
    chromosome_table.loc[:,'Chromosome'] = pd.Categorical(chromosome_table['Chromosome'],categories=chromosome_order,ordered=True)
    chromosome_table = chromosome_table.sort_values(by=['Chromosome']).reset_index(drop=True)
    cum_abs = np.insert(np.cumsum(chromosome_table['Chrom_End'].to_numpy()),0,0)
 
    chromosome_table.loc[:,'Chrom_Abs_Start'] = cum_abs[:-1]
    chromosome_table.loc[:,'Chrom_Abs_End'] = cum_abs[1:]
    chromosome_table = chromosome_table[~chromosome_table['Chromosome'].isin(['X','Y'])]
    return chromosome_table


def norm_proportions(event_proportions):
    positions_diff = np.array(event_proportions['Positions'])[1::2]-np.array(event_proportions['Positions'])[::2]
    proportions = np.array(event_proportions['Proportions'])[::2]


    event_proportions['Proportions'] = np.array(event_proportions['Proportions'])/(np.sum(np.multiply(positions_diff,proportions))*1e-9)
    
def plot_events(event_proportions_store,cancer_type,events_to_plot,plot_type,apply_penalty,event_type,n_samples,norm=False,corrected=False,major_cn='All'):
    os.makedirs(f'../plots/wgd_diff_gain_proportion_plots_probabilistic_rewrite/{event_type}/',exist_ok=True)
    color_lookup = {'Pre_WGD':'#168a35','Close':'#8acfe6','Post_WGD':'#5d168a','Non_WGD':'#de4e26','WGD':'#08bdba','Post_WGD_Corrected':'#5d168a'}
    chromosome_table =load_chromosome_table()
    event_proportions_store = copy.deepcopy(event_proportions_store)
    if norm:
        for event_name in event_proportions_store:
            norm_proportions(event_proportions_store[event_name])
    fig = go.Figure()
    for i,event_name in enumerate(events_to_plot):
        fig.add_trace(go.Scatter(x=event_proportions_store[event_name]['Positions'], y=event_proportions_store[event_name]['Proportions'],
                            mode='lines',
                            marker={'color':color_lookup[event_name]},
                            name=event_name.replace('_','-')))

    fig.update_layout(
        xaxis = dict(
            tickmode = 'array',
            tickvals = chromosome_table['Chrom_Abs_Start'],
            ticktext = chromosome_table['Chromosome']
        ),

    )
    fig.update_layout(
    legend_title="",
    font=dict(
        family="Helvetica",
        size=18,
        color="Black"
    )
)
    plot_title = f"{cancer_type} (n={n_samples})" if major_cn == 'All' else f"{cancer_type} (n={n_samples}) Major CN = {major_cn}"
    y_axis_title = "Frequency" if not norm else "Normalized Frequency"
    fig.update_layout(
    title=plot_title,
    xaxis_title="Genome Position",
    yaxis_title=y_axis_title,
    legend_title="Relative Timing")

    fig.update_layout(template='plotly_white')
    fig.update_annotations(font_size=18)
    if event_type =='Gain':
        output_path = f"../plots/wgd_diff_gain_proportion_plots_probabilistic_rewrite/{event_type}/{cancer_type}_wgd_diff_proportion_{plot_type}_norm_{norm}_apply_penalty_{apply_penalty}_major_cn_{major_cn}.pdf".replace(' ','_')
    elif event_type == 'Loss':
        output_path = f"../plots/wgd_diff_gain_proportion_plots_probabilistic_rewrite/{event_type}/{cancer_type}_wgd_diff_proportion_{plot_type}_norm_{norm}_corrected_{corrected}_major_cn_{major_cn}.pdf".replace(' ','_')
    fig.write_image(output_path,width=1200)

def load_gain_timing_table(apply_penalty,major_cn,mode):
    in_path = f'../output/pre_post_non/pre_post_non_gains_apply_penalty_{apply_penalty}_major_cn_{major_cn}_mode_{mode}.tsv'
    timing_table = pd.read_csv(in_path,sep='\t')
    timing_table['Chromosome'] = timing_table['Segment_ID'].str.split('-').str[0]
    timing_table['Segment_Start'] = timing_table['Segment_ID'].str.split('-').str[1].astype(int)
    timing_table['Segment_End'] = timing_table['Segment_ID'].str.split('-').str[2].astype(int)

    timing_table = timing_table.rename(columns={'Pre_Probability':'Pre_WGD_Gain','Post_Probability':'Post_WGD_Gain','Non_Probability':'Non_WGD_Gain'})
    timing_table = timing_table[['Sample_ID','Segment_ID','Chromosome','Segment_Start','Segment_End','Pre_WGD_Gain','Post_WGD_Gain','Non_WGD_Gain']]
    return timing_table

def load_loss_timing_table():
    in_path = f'../output/pre_post_non/pre_post_non_loss.tsv'
    timing_table = pd.read_csv(in_path,sep='\t')
    timing_table['Chromosome'] = timing_table['Segment_ID'].str.split('-').str[0]
    timing_table['Segment_Start'] = timing_table['Segment_ID'].str.split('-').str[1].astype(int)
    timing_table['Segment_End'] = timing_table['Segment_ID'].str.split('-').str[2].astype(int)
    timing_table = timing_table.rename(columns={'Pre_Probability':'Pre_WGD_Loss','Post_Probability':'Post_WGD_Loss','Non_Probability':'Non_WGD_Loss'})
    timing_table = timing_table[['Sample_ID','Segment_ID','Chromosome','Segment_Start','Segment_End','Pre_WGD_Loss','Post_WGD_Loss','Non_WGD_Loss']]
    return timing_table

if __name__ == '__main__':
    RNG = np.random.default_rng(42)
    run_n = int(sys.argv[1])
    for apply_penalty in  [True,False]:
        for major_cn in ['All',3,4,5,6,7]:
            if major_cn !='All':
                continue
            for event_type in ['Gain','Loss']:
                if event_type == 'Loss':
                    if apply_penalty or major_cn != 'All':
                        continue
                    
                
                wgd_status_table = DataTools.load_classifications()[['Cancer_Type','Sample_ID','WGD_Status']]
            
                gain_timing_table = load_gain_timing_table(apply_penalty=apply_penalty,major_cn=major_cn,mode='absolute')
                gain_timing_table = gain_timing_table.merge(wgd_status_table)
                
                loss_timing_table = load_loss_timing_table()
                loss_timing_table = loss_timing_table.merge(wgd_status_table)

                if event_type == 'Gain':
                    timing_table = gain_timing_table
                elif event_type == 'Loss':
                    timing_table = loss_timing_table



                cancer_type_intersection = set(gain_timing_table['Cancer_Type'].unique()).intersection(set(loss_timing_table['Cancer_Type'].unique()))
                cancer_type_intersection = list(cancer_type_intersection)

                
                timing_table = timing_table[timing_table['Cancer_Type'].isin(cancer_type_intersection)]
        
                
                
                timing_table = timing_table.dropna(subset=['Segment_Start','Segment_End'])
                #round to nearest 1000bp
                timing_table.loc[:,'Segment_Start'] = np.round(timing_table['Segment_Start']/1000).astype(int)*1000
                timing_table.loc[:,'Segment_End'] = np.round(timing_table['Segment_End']/1000).astype(int)*1000
                timing_table = timing_table[timing_table['Segment_Start']!=timing_table['Segment_End']]


                timing_table = timing_table[~pd.isnull(timing_table[f'Pre_WGD_{event_type}'])]
                timing_table.loc[:,f'WGD_{event_type}'] = timing_table[f'Pre_WGD_{event_type}'].to_numpy()+timing_table[f'Post_WGD_{event_type}'].to_numpy()
                timing_table.loc[:,f'WGD_{event_type}'] = np.minimum(timing_table[f'WGD_{event_type}'],1)
                
                timing_table_all = timing_table.copy()
                timing_table_all.loc[:,'Cancer_Type'] = 'All'
                
                
                timing_table = pd.concat([timing_table,timing_table_all])
                cancer_types = list(timing_table['Cancer_Type'].astype(str).unique())
                cancer_types.sort()
        
            
            
                
                for cancer_type,cancer_type_table in timing_table.groupby('Cancer_Type'):

                    if cancer_type != cancer_types[run_n]:
                        continue


                    if cancer_type == 'All':
                        wgd_status_table_cancer_type = wgd_status_table
                    else:
                        wgd_status_table_cancer_type = wgd_status_table[wgd_status_table['Cancer_Type']==cancer_type]
                    n_wgd = len(wgd_status_table_cancer_type[wgd_status_table_cancer_type['WGD_Status']]['Sample_ID'].unique())
                    n_non_wgd = len(wgd_status_table_cancer_type[~wgd_status_table_cancer_type['WGD_Status']]['Sample_ID'].unique())
                    if n_wgd <10 or n_non_wgd <10:
                        continue

                    proportion_store = {}
                    proportion_store_normed = {}
                    proportion_line_pre= get_event_proportions(cancer_type_table,weight_col=f'Pre_WGD_{event_type}',n_samples=n_wgd)
                    proportion_line_post= get_event_proportions(cancer_type_table,weight_col=f'Post_WGD_{event_type}',n_samples=n_wgd)

                    proportion_line_non_wgd= get_event_proportions(cancer_type_table,weight_col=f'Non_WGD_{event_type}',n_samples=n_non_wgd)
                    proportion_line_wgd = get_event_proportions(cancer_type_table,weight_col=f'WGD_{event_type}',n_samples=n_wgd)
                    
                    
                    proportion_store['Pre_WGD'] = proportion_line_pre
                    proportion_store['Post_WGD'] = proportion_line_post
                    
                    proportion_store['Non_WGD'] = proportion_line_non_wgd
                    proportion_store['WGD'] = proportion_line_wgd

                    if event_type =='Loss':
                        corrected_post_proportions = {}
                        corrected_post_proportions['Positions'] = proportion_line_post['Positions']
                        corrected_post_proportions['Proportions'] = proportion_line_post['Proportions']/(1-np.clip(proportion_line_pre['Proportions'],0.00,0.95))
                        
                        
                        proportion_store['Post_WGD_Corrected'] = corrected_post_proportions
                    plot_events(proportion_store,cancer_type,['Pre_WGD','Post_WGD'],'WGD_Only',apply_penalty,event_type,n_wgd,major_cn=major_cn)
                    plot_events(proportion_store,cancer_type,['WGD','Non_WGD'],'WGD_Non_WGD',apply_penalty,event_type,n_wgd+n_non_wgd,major_cn=major_cn)

                    plot_events(proportion_store,cancer_type,['Pre_WGD','Post_WGD'],'WGD_Only',apply_penalty,event_type,n_wgd,norm=True,major_cn=major_cn)
                    plot_events(proportion_store,cancer_type,['WGD','Non_WGD'],'WGD_Non_WGD',apply_penalty,event_type,n_wgd+n_non_wgd,norm=True,major_cn=major_cn)

                    if event_type =='Loss':
                        plot_events(proportion_store,cancer_type,['Pre_WGD','Post_WGD_Corrected'],'WGD_Only',apply_penalty,event_type,n_wgd,corrected=True,major_cn=major_cn)
                        plot_events(proportion_store,cancer_type,['Pre_WGD','Post_WGD_Corrected'],'WGD_Only',apply_penalty,event_type,n_wgd,norm=True,corrected=True,major_cn=major_cn)
