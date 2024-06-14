import os
import sys
import numpy as np
import pandas as pd
from gritic import gritictimer,sampletools
import HartwigDataLoader
import glob


hartwig_metadata_path = '/camp/project/proj-vanloo-secure/Hartwig/DR-163-update2/metadata/metadata.tsv'
hartwig_metadata_table = pd.read_csv(hartwig_metadata_path,sep="\t")


processed_count = 0

run_number = int(sys.argv[1])

sample_row = hartwig_metadata_table.iloc[run_number]
    
cn_method = "PURPLE"


sample_id = sample_row['sampleId']


output_dir = f"../output/Hartwig_mut_clean_revisions/{sample_id}"

set_name = sample_row['setName']

sample_data_kataegis_filtered = HartwigDataLoader.load_sample_data(sample_id,set_name,cn_method,filter_kataegis=True)

sample_data_kataegis_filtered['Mutation_Table'].to_csv(f'{output_dir}/mutation_table.tsv',sep='\t',index=False)

nrpcc_table = {'Sample_ID':[],'Average_Coverage':[],'NRPCC':[],'Ploidy':[]}
nrpcc_table['Sample_ID'].append(sample_id)
nrpcc_table['Average_Coverage'].append(sample_data_kataegis_filtered['Average_Coverage'])
nrpcc_table['NRPCC'].append(sample_data_kataegis_filtered['NRPCC'])
nrpcc_table['Ploidy'].append(sample_data_kataegis_filtered['Ploidy'])
nrpcc_table = pd.DataFrame(nrpcc_table)
nrpcc_table.to_csv(f'{output_dir}/nrpcc_table.tsv',sep='\t',index=False)

sample_kataegis_filtered = sampletools.Sample(sample_data_kataegis_filtered['Mutation_Table'],sample_data_kataegis_filtered['CN_Table'],sample_data_kataegis_filtered['Subclone_Table'],sample_data_kataegis_filtered['Sample_ID'],sample_data_kataegis_filtered['Sample_Purity'],sex=sample_data_kataegis_filtered['Sex'])

try:
    gritictimer.process_sample(sample_kataegis_filtered,output_dir)
except ValueError as e:
    print(e)

sample_data_kataegis_filtered['Mutation_Table'].to_csv(f'{output_dir}/mutation_table.tsv',sep='\t',index=False)

nrpcc_table.to_csv(f'{output_dir}/nrpcc_table.tsv',sep='\t',index=False)