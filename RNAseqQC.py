


""" MultiQC module to parse output from bcbioRNASeq Quality control """

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import scatter
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from multiqc.plots import heatmap
from multiqc.plots import scatter
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
            if f['s_name'] == 'rlog':
                rlog = pd.read_csv(join(dirpath, fname))
                # print(rlog.head())
            if f['s_name'] == 'vst':
                vst = pd.read_csv(join(dirpath, fname))
                # print(vst.head())
            if f['s_name'] == 'combined':
                combined_counts = pd.read_csv(join(dirpath, fname), sep='\t')
                # print(combined_counts.head())
            if f['s_name'] == 'gene.est':
                genes_est = pd.read_csv(join(dirpath, fname))
                # print(genes_est.head())
            if f['s_name'] == 'gene.final':
                genes_final = pd.read_csv(join(dirpath, fname))
                # print(genes_final.head())
            if f['s_name'] == 'gene.fitted':
                genes_fitted = pd.read_csv(join(dirpath, fname))
                # print(genes_fitted.head())
            if f['s_name'] == 'corMatrix':
                raw_data = pd.read_csv(join(dirpath, fname))
                # print(raw_data.head())

        col_names = list(raw_counts)[1:]
        group_num = len(col_names)

        raw_counts['sum'] = raw_counts.sum(axis=1)

        self.plot_correlation_heatmap(raw_counts, norm_counts, col_names, group_num)
        self.plot_mean_sd(raw_counts, norm_counts, col_names, group_num, vst, rlog, combined_counts)
        self.plot_disp_ests(combined_counts, genes_est,genes_final,genes_fitted)
        self.plot_covariates(raw_data)


    def plot_covariates(self, raw_data):

        comp_names = []
        for iter, row in raw_data.iterrows():
            comp_names.append(row['compare'])
            if iter % 3 == 2:
                break

        cor_names = []
        for _, row in raw_data.iterrows():
            if not (row['covar'] in cor_names):
                cor_names.append(row['covar'])

        hmdata = [[None for _ in range(len(cor_names))] for _ in range(3)]
        for iter,row in raw_data.iterrows():
            j = int(iter / 3)
            i = iter % 3
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



    def plot_disp_ests(self, combined_counts, genes_est,genes_final,genes_fitted):

        pconfig = {
            'title': 'bcbioRNASeq Quality Control: Dispersion Estimates plot',
            'xlab': 'mean of normalized counts',
            'ylab': 'dispersion',
            'xLog': True,
            'yLog': True,
            'showInLegend': True,
            'tt_label': 'mean: {point.x}<br/> dispersion: {point.y}'
        }

        genes_est = set_gene_names(genes_est, combined_counts)
        genes_final = set_gene_names(genes_final, combined_counts)
        genes_fitted = set_gene_names(genes_fitted, combined_counts)

        data = dict()
        colors = ['#fd0000','#1803ff','#000000']
        labels = ['fitted', 'final', 'est']
        for _, row in genes_final.iterrows():
            data[labels[0] + '_' + row['HUGO']] = [{
                'x': row['px'],
                'y': row['py3'],
                'color': colors[0]
            }]
        for _, row in genes_fitted.iterrows():
            data[labels[1] + '_' + row['HUGO']] = [{
                    'x': row['px'],
                    'y': row['py2'],
                    'color': colors[1]
                }]
        for _, row in genes_est.iterrows():
            data[labels[2] + '_' + row['HUGO']] = [{
                'x': row['px'],
                'y': row['V2'],
                'color': colors[2]
            }]

        legend = ''
        label_style = 'font-family: \'Lucida Grande\', \'Lucida Sans Unicode\', Arial, Helvetica, sans-serif; ' \
                      'font-size: 12px; ' \
                      'font-weight: bold; '
        legend += '<center><div>'
        legend += '<span style="' + label_style + ' margin-right: 10px;">Legend: </span>'
        for i in range(3):
            legend += '<span style="white-space: nowrap;">'
            legend += '<span style="display: inline-block; width: 16px; height: 12px; ' + \
                      '            margin-bottom: -1px; margin-right: 1px; background-color: ' + colors[i] + '"></span>'
            legend += '<span style="' + label_style + ' margin-right: 20px; white-space: normal;"> ' + labels[i] + '</span>'
            legend += '</span>'
        legend += '</div></center>'
        self.add_section (
            name = 'Dispersion',
            anchor = 'dispests_section',
            description = 'The following plot shows the dispersion by mean of normalized counts. '
                          'We expect the dispersion to decrease as the mean of normalized counts increases.',
            plot = legend + scatter.plot(data, pconfig)
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

    def plot_mean_sd(self, raw_counts, norm_counts, col_names, group_num, vst, rlog, combined_counts):

        pconfig1 = get_config('log2')
        pconfig2 = get_config('rlog')
        pconfig3 = get_config('vst')

        desc = 'The plots below show the standard deviation of normalized counts (normalized_counts) ' \
               'using log2(), rlog(), and variance stabilizing (vst()) transformations by rank (mean). ' \
               'The transformations greatly reduce the standard deviation, with rlog() stabilizing the variance best across the mean.'

        norm_counts = norm_counts[raw_counts['sum'] > 0]
        non_zero = norm_counts

        non_zero['mean'] = non_zero[col_names].mean(axis=1)
        non_zero['mean_rank'] = non_zero['mean'].rank()

        rlog['mean'] = rlog[col_names].mean(axis=1)
        rlog['mean_rank'] = rlog['mean'].rank()

        # print(col_names)
        # print(vst.head())

        vst['mean'] = vst[col_names].mean(axis=1)
        vst['mean_rank'] = vst['mean'].rank()

        non_zero = non_zero.loc[non_zero['rowname'].isin(combined_counts['id'])]
        d = combined_counts.set_index('id')['HUGO'].to_dict()
        non_zero['HUGO'] = non_zero['rowname'].map(d)
        data1 = dict()
        for _, row in non_zero.iterrows():
            dlog2 = [log2(row[col]) for col in col_names if row[col]]
            if any(dlog2):
                data1[row['HUGO']] = [{
                        'x': row['mean_rank'],
                        'y': np.std(dlog2)
                    }]

        rlog = rlog.loc[rlog['rowname'].isin(combined_counts['id'])]
        d = combined_counts.set_index('id')['HUGO'].to_dict()
        rlog['HUGO'] = rlog['rowname'].map(d)
        data2 = dict()
        for _, row in rlog.iterrows():
            sd = row[col_names].std()
            data2[row['HUGO']] = [{
                    'x': row['mean_rank'],
                    'y': sd
                }]

        vst = vst.loc[vst['rowname'].isin(combined_counts['id'])]
        d = combined_counts.set_index('id')['HUGO'].to_dict()
        vst['HUGO'] = vst['rowname'].map(d)
        data3 = dict()
        for _, row in vst.iterrows():
            sd = row[col_names].std()
            data3[row['HUGO']] = [{
                    'x': row['mean_rank'],
                    'y': sd
                }]

        self.add_section (
            name = 'Variance stabilization',
            anchor = 'meansd_section',
            description = desc,
            plot = scatter.plot([data1, data2, data3], pconfig1)
                    )

def get_config(plot_id):
        return({
            'title': 'bcbioRNASeq Quality Control: MeanSD plot',
            'xlab': 'rank (mean)',
            'data_labels':[
                {'name': 'log2', 'ylab': 'log2'},
                {'name': 'rlog', 'ylab': 'rlog'},
                {'name': 'vst', 'ylab': 'vst'}]})
                #'tt_label': 'rank(mean): {point.x}<br/> sd('+ plot_id + '): {point.y}'
           # )

def set_gene_names(data_frame, combined_counts):
    # print(data_frame.head())
    # print(combined_counts.head())

    data_frame = data_frame.loc[data_frame['gene.names'].isin(combined_counts['id'])]

    d = combined_counts.set_index('id')['HUGO'].to_dict()
    data_frame['HUGO'] = data_frame['gene.names'].map(d)
    data = dict()
    return data_frame

