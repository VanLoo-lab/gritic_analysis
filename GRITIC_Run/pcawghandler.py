import os
import sys
import numpy as np
import pandas as pd
from gritic import gritictimer,sampletools
import sys
import PCAWGDataLoader
import glob

import pickle


run_number = int(sys.argv[1])


sample_table = pd.read_csv('/camp/project/proj-vanloo-secure/PCAWG/ICGC_annotations/summary_table_combined_annotations_v2.txt',sep="\t")



sample_ids = list(sample_table['samplename'])


sample_id = sample_ids[run_number]


    
output_dir = f"../output/PCAWG/{sample_id}"
os.makedirs(output_dir,exist_ok=True)



nrpcc_table = {'Sample_ID':[],'Average_Coverage':[],'NRPCC':[],'Ploidy':[]}


sample_data_kataegis_filtered = PCAWGDataLoader.load_sample_data(sample_id,phased=False,filter_kataegis=True,cn_method='All')

if sample_data_kataegis_filtered is None:
    raise ValueError(f'No data for {sample_id}')


nrpcc_table['Sample_ID'].append(sample_id)
nrpcc_table['Average_Coverage'].append(sample_data_kataegis_filtered['Average_Coverage'])
nrpcc_table['NRPCC'].append(sample_data_kataegis_filtered['NRPCC'])
nrpcc_table['Ploidy'].append(sample_data_kataegis_filtered['Ploidy'])

nrpcc_table = pd.DataFrame(nrpcc_table,index=[0])
nrpcc_table.to_csv(f'{output_dir}/nrpcc_table.tsv',sep='\t',index=False)


sample_kataegis_filtered = sampletools.Sample(sample_data_kataegis_filtered['Mutation_Table'],sample_data_kataegis_filtered['CN_Table'],sample_data_kataegis_filtered['Subclone_Table'],sample_data_kataegis_filtered['Sample_ID'],sample_data_kataegis_filtered['Sample_Purity'],sex=sample_data_kataegis_filtered['Sex'])


try:
    gritictimer.process_sample(sample_kataegis_filtered,output_dir,plot_trees=False)
except ValueError as e:
    #print (e)
    print('error',sample_id,e)

sample_data_kataegis_filtered['Mutation_Table'].to_csv(f'{output_dir}/mutation_table.tsv',sep='\t',index=False)

nrpcc_table.to_csv(f'{output_dir}/nrpcc_table.tsv',sep='\t',index=False)

with open(f'{output_dir}/timing_complete.txt','w') as f:
    f.write('done')
