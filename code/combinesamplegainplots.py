import os
import PyPDF2
import numpy as np

import DataTools

import matplotlib.pyplot as plt

RNG = np.random.default_rng(1234)

n_samples = 5
wgd_samples = DataTools.load_route_table(apply_penalty=False)[['Sample_ID','Cancer_Type','WGD_Timing']].drop_duplicates()
wgd_samples = wgd_samples[~np.isnan(wgd_samples['WGD_Timing'])]

wgd_bins = np.linspace(0,1,n_samples+1)

def check_paths(sample_id):
    for apply_penalty in [False,True]:
        test_path = f'../plots/sample_timing_plots/prior_{apply_penalty}/{sample_id}_sample_timing.pdf'
        if not os.path.exists(test_path):
            return False
    return True
def get_sample(wgd_samples,bins,cancer_types):
    sample = wgd_samples[(wgd_samples['WGD_Timing'] >= bins[0]) & (wgd_samples['WGD_Timing'] < bins[1])]
    while True:
        sample_id = RNG.choice(sample['Sample_ID'])
        if sample_id not in sample_ids:
            if not check_paths(sample_id):
                continue
            cancer_type = sample[sample['Sample_ID'] == sample_id]['Cancer_Type'].iloc[0]
            if cancer_type in cancer_types:
                continue
            return sample_id,cancer_type

def combine_pdfs(pdf_files,output_path):
    pdfs = []
    for pdf_file in pdf_files:
        pdf = open(pdf_file, 'rb')
        pdfs.append(PyPDF2.PdfFileReader(pdf))
    
    #reverse the order of the pdfs
    pdfs = pdfs[::-1]
    # Get the first page of each PDF and determine their size
    pages = [pdf.getPage(0) for pdf in pdfs]
    widths, heights = [], []
    for page in pages:
        llx, lly, urx, ury = page.mediaBox.lowerLeft + page.mediaBox.upperRight
        widths.append(urx - llx)
        heights.append(ury - lly)

    output_pdf = PyPDF2.PdfFileWriter()

    # Create a new page for each row
    page = output_pdf.addBlankPage(max(widths), sum(heights))
    y = 0
    for i in range(0, len(pages)):
        # Add each PDF to the corresponding row
        page.mergeTranslatedPage(pages[i], 0, y)
        y += heights[i]
    output_pdf.addPage(page)
    with open(output_path, 'wb') as f:
        output_pdf.write(f)

    # Close the PDF files
    for pdf in pdfs:
        pdf.stream.close()

sample_ids = []
cancer_types = []
pdf_paths = {True:[],False:[]}
for i in range(n_samples):
    new_sample,new_cancer_type= get_sample(wgd_samples,(wgd_bins[i],wgd_bins[i+1]),cancer_types)
    sample_ids.append(new_sample)
    cancer_types.append(new_cancer_type)
    for apply_penalty in [False,True]:
        pdf_paths[apply_penalty].append(f'../plots/sample_timing_plots/prior_{apply_penalty}/{new_sample}_sample_timing.pdf')

#combine pdfs to one page
for apply_penalty in [False,True]:
    output_path = f'../plots/combined_sample_timing_plots/combined_sample_timing_prior_{apply_penalty}.pdf'
    combine_pdfs(pdf_paths[apply_penalty],output_path)