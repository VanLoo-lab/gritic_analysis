import os
import pandas as pd
import numpy as np
from gritic import dataloader
import PCAWGPhaser
def load_purity_wgd_status_sex(sample_dir,donor_id):
    #use f string instead

    purity_path= f'{sample_dir}/{donor_id}T.purple.purity.tsv'
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


def load_cn_table(sample_dir,donor_id):
    cn_path = f'{sample_dir}/{donor_id}T.purple.cnv.somatic.tsv'
    cn_table = pd.read_csv(cn_path,sep="\t",dtype={'chromosome':str})
    
    cn_table = cn_table[['chromosome','start','end','majorAlleleCopyNumber','minorAlleleCopyNumber']]
    cn_table = cn_table.dropna().copy()
    cn_table['majorAlleleCopyNumber'] = cn_table['majorAlleleCopyNumber'].round().astype(int)
    cn_table['minorAlleleCopyNumber'] = cn_table['minorAlleleCopyNumber'].round().astype(int)
    cn_table = cn_table.rename(columns={'chromosome':'Chromosome','start':'Segment_Start','end':'Segment_End','majorAlleleCopyNumber':'Major_CN','minorAlleleCopyNumber':'Minor_CN'})

    return cn_table.reset_index(drop=True)


def load_vcf(sample_dir,donor_id):
    
    vcf_path = f'{sample_dir}/{donor_id}T.purple.somatic.vcf.gz'
    
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
   

def load_snv_table(sample_dir,donor_id):
    return load_vcf(sample_dir,donor_id)


def load_subclone_table(sample_dir,donor_id):

    subclone_dir = '/PCAWG_DIR/'
    subclone_path = f'{subclone_dir}/DPClust_output/{donor_id}T_optimaInfo.txt'
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

def filter_mutation_table_for_phasing(mutation_table,wgd_status):
    filtered_table = []
    for segment_id,segment_table in mutation_table.groupby('Segment_ID'):
        if  not segment_table['Phasing'].isnull().all():
            filtered_table.append(segment_table)
        if segment_table['Major_CN'].iloc[0]==2 and wgd_status:
            filtered_table.append(segment_table)
    if len(filtered_table) ==0:
        return None
    return pd.concat(filtered_table).reset_index(drop=True)

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
def load_sex(sample_id):
    sample_metadata_path = '/PCAWG_DIR/summary_table_combined_annotations_v2.txt'
    sample_metadata = pd.read_csv(sample_metadata_path,sep="\t")
    sample_metadata = sample_metadata[['samplename','inferred_sex']]
    sample_row = sample_metadata[sample_metadata['samplename']==sample_id].iloc[0]
    if sample_row['inferred_sex'] == 'male':
        return 'XY'
    if sample_row['inferred_sex'] == 'female':
        return 'XX'
    return None
def load_associated_codes(sample_id):
    sample_metadata_path = '/PCAWG_DIR/summary_table_combined_annotations_v2.txt'
    sample_metadata = pd.read_csv(sample_metadata_path,sep="\t")
    sample_row = sample_metadata[sample_metadata['samplename']==sample_id].iloc[0]
    donor_id = sample_row['icgc_donor_id']
    project_code = sample_row['projectcode']
    is_tcga = not pd.isnull(sample_row['tcga_donor_uuid'])
    return donor_id,project_code,is_tcga

def get_sample_dir(is_tcga, project_code, donor_id):
    snv_dir = '/PCAWG_DIR/'
    if is_tcga:
        sample_dir = f"{snv_dir}/TCGA/{project_code}/{donor_id}/purple/"
    else:
        sample_dir = f"{snv_dir}/ICGC/{project_code}/{donor_id}/purple/"
    
    return sample_dir
def load_sample_data(sample_id,phased=False,filter_kataegis=False,cn_method='All'):

    donor_id,project_code,tcga_status = load_associated_codes(sample_id)
    sample_dir = get_sample_dir(tcga_status,project_code,donor_id)
    if not os.path.exists(sample_dir):
        return None
    sample_purity,sample_wgd_status,sample_sex = load_purity_wgd_status_sex(sample_dir,donor_id)

    cn_table = load_cn_table(sample_dir,donor_id)
    cn_table = dataloader.merge_segments(cn_table)


    snv_table = load_snv_table(sample_dir,donor_id)

    mutation_table = dataloader.assign_cn_to_snv(snv_table,cn_table)

    mutation_table['Phasing'] = np.nan
    if filter_kataegis:
        mutation_table = run_kataegis_filtering(mutation_table,sample_id)
        
    subclone_table = load_subclone_table(sample_dir,donor_id)
    subclone_table = dataloader.filter_excess_subclones(subclone_table)
    mutation_table = mutation_table[mutation_table['Major_CN']>0]
    

    nrpcc,average_coverage,ploidy = dataloader.calculate_nrpcc(cn_table,mutation_table,sample_purity,sample_sex)

    sample_data = {
        "Sample_ID": sample_id,
        "Sample_Purity": sample_purity,
        "Mutation_Table": mutation_table,
        "Subclone_Table": subclone_table,
        'CN_Table':cn_table,
        'NRPCC':nrpcc,
        'Ploidy':ploidy,
        'Average_Coverage':average_coverage,
        'Sex':sample_sex
    }
    return sample_data
