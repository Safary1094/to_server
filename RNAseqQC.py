


""" MultiQC module to parse output from bcbioRNASeq Quality control """

from multiqc.modules.base_module import BaseMultiqcModule
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import scipy.stats as st
from multiqc.plots import scatter, heatmap
from math import log2, log10
import logging
from os.path import join


# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        mod_name = 'RNAseqQC'

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='bcbiornaseqqc', anchor=mod_name)
        file_names, roots = [], []
        for f in self.find_log_files(mod_name, filecontents=False):
            print(f)
            dirpath, fname = f['root'], f['fn']
            if f['s_name'] == 'rawCounts':
                raw_counts = pd.read_csv(join(dirpath, fname))
                # print(raw_counts.head())
            if f['s_name'] == 'normalizedCounts':
                norm_counts = pd.read_csv(join(dirpath, fname))
                # print(norm_counts.head())
            if f['s_name'] == 'corMatrix':
                raw_data = pd.read_csv(join(dirpath, fname))
                # print(raw_data.head())
            if f['s_name'] == 'pca':
                pca_data = pd.read_csv(join(dirpath, fname))


        col_names = list(raw_counts)[1:]
        group_num = len(col_names)

        raw_counts['sum'] = raw_counts.sum(axis=1)

        self.plot_correlation_heatmap(raw_counts, norm_counts, col_names, group_num)
        #self.plot_mean_sd(raw_counts, norm_counts, col_names, group_num, vst, rlog, combined_counts)
        #self.plot_disp_ests(combined_counts, genes_est,genes_final,genes_fitted)
        self.plot_covariates(raw_data)
        self.plot_pca(pca_data)

    def plot_pca(self, pca_data):
        standard_colors = [
            '#0000FF', '#008000', '#FFA500', '#FF00FF', '#CCCC00', '#800000',
            '#00CCCC', '#808080', '#800080', '#808000', '#000080', '#008080', '#00FF00',]

        color_by_sam = {}
        for i in pca_data.index.tolist():
            pca_data.at[i, 'color'] = standard_colors[i % len(standard_colors)]
            color_by_sam[pca_data.at[i, 'name']] = standard_colors[i % len(standard_colors)]

        pca_data = pca_data.set_index('rowname')
        pca_data = pca_data[['pc1', 'pc2', 'name', 'color']]
        pca_data = pca_data.rename(columns = {'pc1':'x', 'pc2':'y'})

        # color_by_sam = pd.DataFrame(pca_data['color'])
        # color_by_sam = color_by_sam.to_dict('index')

        pca_data = pca_data.to_dict('index')

        pconfig = {'colors': color_by_sam}

        self.add_section (
            name = 'PCA plot',
            anchor = 'pca',
            description = 'PCA is a popular method that is based on the principles of dimensional reduction. Below is a PCA plot of the samples within the space of the first two principal components that explain the most variation in the data. These were calculated using the read counts of the top 1000 most variable genes within the dataset.',

            plot = scatter.plot(pca_data, pconfig)
        )


    def plot_covariates(self, raw_data):

        # find number of PCs
        pc_num = raw_data['covar'].value_counts().tolist()
        pc_num = pc_num[0]

        comp_names = []
        for iter, row in raw_data.iterrows():
            comp_names.append(row['compare'])
            if iter % pc_num == 2:
                break

        cor_names = []
        for _, row in raw_data.iterrows():
            if not (row['covar'] in cor_names):
                cor_names.append(row['covar'])

        hmdata = [[None for _ in range(len(cor_names))] for _ in range(3)]
        for iter,row in raw_data.iterrows():
            j = int(iter / pc_num)
            i = iter % pc_num
            if row['fdr'] < 0.1:
                hmdata[i][j] = row['r']

        pconfig = {
            'title': "bcbioRNASeq Quality Control: PCA Covariates",                 # Plot title - should be in format "Module Name: Plot Title"
            'square': False,                # Force the plot to stay square? (Maintain aspect ratio)
            'reverseColors': False,        # Reverse the order of the colour axis
            'decimalPlaces': 2,            # Number of decimal places for tooltip
            'legend': True,                # Colour axis key enabled or not
            'borderWidth': 0,              # Border width between cells
            'datalabels': True,            # Show values in each cell. Defaults True when less than 20 samples.
            'datalabel_colour': '<auto>',  # Colour of text for values. Defaults to auto contrast.
        }

        hm_html = heatmap.plot(hmdata, cor_names, comp_names, pconfig)

        self.add_section (
            name = 'PCA Covariates',
            anchor = 'covar_section',
            description = 'When multiple factors may influence the results of a given experiment, it is useful to assess which of them is responsible for the most variance as determined by PCA. '
                          'We adapted the method described by Daily et al. where they integrated a method to correlate covariates with principal components values to determine the importance of each factor. '
                          'Here we are showing the correlational analysis of the rlog transformed count data’s principal components with the metadata covariates of interest. Significant correlations (FDR < 0.1) are shaded from blue (anti-correlated) to red (correlated), with non-significant correlations set to zero and shaded in gray.',
            #helptext = 'Inter-correlation analysis (ICA) is another way to look at how well samples cluster by plotting the correlation between the expression profiles of the samples. Pearson`s correlation coefficient is a measure of how well your data would be fitted by a linear regression.',
            plot = hm_html
        )


    def plot_correlation_heatmap(self, raw_counts, norm_counts, col_names, group_num):

        norm_counts = norm_counts[raw_counts['sum'] > 0]
        hmdata = [[st.pearsonr(norm_counts[col_names[i]], norm_counts[col_names[j]])[0]
                    for i in range(group_num)] for j in range(group_num)]

        pconfig = {
            'title': "bcbioRNASeq Quality Control: Correlation Heatmap",
            'xTitle': None,                # X-axis title
            'yTitle': None,                # Y-axis title
            'min': None,                   # Minimum value (default: auto)
            'max': None,                   # Maximum value (default: auto)
            'square': False,               # Force the plot to stay square? (Maintain aspect ratio)
            'reverseColors': False,        # Reverse the order of the colour axis
            'decimalPlaces': 2,            # Number of decimal places for tooltip
            'legend': True,                # Colour axis key enabled or not
            'borderWidth': 0,              # Border width between cells
            'datalabels': True,            # Show values in each cell. Defaults True when less than 20 samples.
            'datalabel_colour': '<auto>',  # Colour of text for values. Defaults to auto contrast.
        }

        hm_html = heatmap.plot(hmdata, col_names, col_names, pconfig)

        self.add_section (
            name = 'Correlation Heatmap',
            anchor = 'heatmap_section',
            description = 'This heatmap shows Pearson`s correlation values between groups.',
            helptext = 'Inter-correlation analysis (ICA) is another way to look at how well samples cluster by plotting the correlation between the expression profiles of the samples. Pearson`s correlation coefficient is a measure of how well your data would be fitted by a linear regression.',
            plot = hm_html
        )
