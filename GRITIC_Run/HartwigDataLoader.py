import os
import pandas as pd
import numpy as np
from gritic import dataloader


def load_purity_wgd_status_sex(sample_id,set_name,hartwig_dir):

    purity_path = f"{hartwig_dir}/{sample_id}/purple/{sample_id}.purple.purity.tsv"
    
    purity_table = pd.read_csv(purity_path,sep="\t")
    sample_purity = purity_table['purity'].iloc[0]
    sample_wgd_status = str(purity_table['wholeGenomeDuplication'].iloc[0]).lower() =='true'
    gender = purity_table['gender'].iloc[0]
    sex = None
    if gender.lower() == 'male':
        sex = 'XY'
    elif gender.lower() =='female':
        sex = 'XX'
    return sample_purity,sample_wgd_status,sex


def load_cn_table(sample_id,set_name,hartwig_dir,method):
    if method == "PURPLE":

        cn_path = f'{hartwig_dir}/{sample_id}/purple/{sample_id}.purple.cnv.somatic.tsv'
        cn_table = pd.read_csv(cn_path,sep="\t",dtype={'chromosome':str})
        
        cn_table = cn_table[['chromosome','start','end','majorAlleleCopyNumber','minorAlleleCopyNumber']]
        cn_table = cn_table.dropna().copy()
        cn_table['majorAlleleCopyNumber'] = cn_table['majorAlleleCopyNumber'].round().astype(int)
        cn_table['minorAlleleCopyNumber'] = cn_table['minorAlleleCopyNumber'].round().astype(int)
        cn_table = cn_table.rename(columns={'chromosome':'Chromosome','start':'Segment_Start','end':'Segment_End','majorAlleleCopyNumber':'Major_CN','minorAlleleCopyNumber':'Minor_CN'})

    return cn_table.reset_index(drop=True)

def load_vcf(sample_id,hartwig_dir):
    
    vcf_path = f'{hartwig_dir}/{sample_id}/purple/{sample_id}.purple.somatic.vcf.gz'
    
    vcf = pd.read_csv(vcf_path,sep='\t',comment='#',header=None,names=['Chromosome','Position','ID','Reference_Allele','Tumor_Seq_Allele2','QUAL','FILTER','INFO','FORMAT','Ref','Tumor'],dtype={'Chromosome':str},index_col=False)

    vcf = vcf[vcf['Reference_Allele'].isin(['A','C','G','T'])]
    vcf = vcf[vcf['Tumor_Seq_Allele2'].isin(['A','C','G','T'])]
    
    vcf = vcf[vcf['FILTER']=='PASS']
    ad_index = vcf['FORMAT'].str.split(':').apply(lambda x: x.index('AD'))
    gt_index = vcf['FORMAT'].str.split(':').apply(lambda x: x.index('GT'))
    assert len(ad_index.unique())==1
    assert len(gt_index.unique())==1
    ad_index = ad_index.iloc[0]
    gt_index = gt_index.iloc[0]
    vcf['Tumor_AD'] = vcf['Tumor'].str.split(':').apply(lambda x: x[ad_index])
    vcf['Tumor_GT'] = vcf['Tumor'].str.split(':').apply(lambda x: x[gt_index])
    vcf['Ref_GT'] = vcf['Ref'].str.split(':').apply(lambda x: x[gt_index])
    vcf['Ref_AD'] = vcf['Ref'].str.split(':').apply(lambda x: x[ad_index])
    #vcf = vcf[vcf['Tumor_GT'] == '0/1']
    vcf = vcf[vcf['Ref_GT'] == '0/0']
    #assert len(vcf['Tumor_GT'].unique()) ==1
    #assert vcf['Tumor_GT'].iloc[0] == '0/1'
    vcf['Tumor_Ref_Count'] = vcf['Tumor_AD'].str.split(',').apply(lambda x: x[0])
    vcf['Tumor_Alt_Count'] = vcf['Tumor_AD'].str.split(',').apply(lambda x: x[1])
    
    
    
    vcf = vcf[['Chromosome','Position','Tumor_Ref_Count','Tumor_Alt_Count']]
    
    vcf['Tumor_Alt_Count'] = vcf['Tumor_Alt_Count'].astype(int)
    vcf['Tumor_Ref_Count'] = vcf['Tumor_Ref_Count'].astype(int)
    
    vcf = vcf[vcf['Tumor_Alt_Count']>=3]
    
    return vcf
def load_snv_table(sample_id,hartwig_dir):
    return load_vcf(sample_id,hartwig_dir)
def run_kataegis_filtering(mutation_table,sample_id):
    mutation_table = mutation_table.copy()
    kataegis_calls_path = f'/KATAEGIS_RESULTS/{sample_id}_kataegis_annotated.txt'

    if not os.path.exists(kataegis_calls_path):
        return mutation_table
    kataegis_calls = pd.read_csv(kataegis_calls_path,sep="\t")
    kataegis_calls = kataegis_calls.rename(columns={'chr':'Chromosome','pos':'Position'})
    kataegis_calls['Position'] = kataegis_calls['Position']
    mutation_table['Pos_ID'] = mutation_table['Chromosome'].astype(str)+"_"+mutation_table['Position'].astype(str)
    kataegis_calls['Pos_ID'] = kataegis_calls['Chromosome'].astype(str)+"_"+kataegis_calls['Position'].astype(str)
    mutation_table = mutation_table[~mutation_table['Pos_ID'].isin(kataegis_calls['Pos_ID'])].copy()

    return mutation_table.drop(columns=['Pos_ID']).reset_index(drop=True)

def load_subclone_table(sample_id,set_name):
    subclone_path = f'/DPCLUST_RESULTS/DPClust/OUT/{set_name}/DPClust_output_V2/{sample_id}_optimaInfo.txt'

    subclone_table = pd.read_csv(subclone_path,sep="\t")

    subclone_table = subclone_table[['cluster.no','location','no.of.mutations']]
    subclone_table = subclone_table.rename(columns={'cluster.no':'Cluster','location':'Subclone_CCF'})

    subclone_table['Subclone_Fraction'] = subclone_table['no.of.mutations']/subclone_table['no.of.mutations'].sum()
    
    subclone_table = subclone_table[subclone_table['Subclone_CCF']<0.9].copy()
    subclone_table = subclone_table[subclone_table['Subclone_Fraction']>0.001].copy()
    subclone_table = subclone_table.drop(columns=['no.of.mutations'])
    if len(subclone_table.index) ==0:
        return None
    return subclone_table


def load_sample_data(sample_id,set_name,cn_method,filter_kataegis=False,wgd_mode_error=True):
    hartwig_dir = 'HARTWIG_DIR'
    sample_purity,sample_wgd_status,sex = load_purity_wgd_status_sex(sample_id,set_name,hartwig_dir)

    cn_table = load_cn_table(sample_id,set_name,hartwig_dir,cn_method)
    cn_table = dataloader.merge_segments(cn_table)

    snv_table = load_snv_table(sample_id,set_name,hartwig_dir)
    if filter_kataegis:
        snv_table = run_kataegis_filtering(snv_table,sample_id)
    mutation_table = dataloader.assign_cn_to_snv(snv_table,cn_table)

    
    mutation_table = mutation_table[mutation_table['Major_CN']>0]
    
    subclone_table = load_subclone_table(sample_id,set_name)

    subclone_table = dataloader.filter_excess_subclones(subclone_table)
    nrpcc,average_coverage,ploidy = dataloader.calculate_nrpcc(cn_table,mutation_table,sample_purity,sex)
    sample_data = {
        "Sample_ID": sample_id,
        "Sample_Purity": sample_purity,
        "Mutation_Table": mutation_table,
        "Subclone_Table": subclone_table,
        'CN_Table':cn_table,
        'NRPCC':nrpcc,
        'Ploidy':ploidy,
        'Average_Coverage':average_coverage,
        'Sex':sex
    }
    return sample_data
