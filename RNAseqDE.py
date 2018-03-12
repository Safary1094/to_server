#!/usr/bin/env python

""" MultiQC module to add link to Bcl2fastq reports """

import logging
from os.path import join
import pandas as pd
import numpy as np


from multiqc.modules.base_module import BaseMultiqcModule
from multiqc import config
from multiqc.plots import table, scatter, heatmap
from collections import OrderedDict
from ngs_utils.bcbio import BcbioProject
from ngs_reporting.reference_data import get_key_genes_file
import xml.etree.ElementTree as ET



class MultiqcModule(BaseMultiqcModule):


    def __init__(self):

        mod_name = 'RNAseqDE'

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='RNA Differential Expression', anchor=mod_name)
        for f in self.find_log_files(mod_name, filecontents=False):
            print(f)
            dirpath, fname = f['root'], f['fn']
            if f['s_name'] == 'RNA_DE':
                de_path = join(dirpath, fname)
            if f['s_name'] == 'RNA_HM':
                hm_path = join(dirpath, fname)

        data_parsed = pd.read_csv(de_path)
        data_parsed = data_parsed.dropna()

        key_genes_ind = data_parsed['is_key'] == True
        data_key = data_parsed[key_genes_ind]

        def addVolcanoPlot(data_key):

            data = {}

            for d in data_key.iterrows():

                point = {}
                point['x'] = d[1]['lfc']
                point['y'] = d[1]['p']

                if abs(d[1]['p']) > 1 and abs(d[1]['lfc']) > 1:
                    point['color'] = '#b2df8a'
                else:
                    point['color'] = '#1f78b4'

                data[d[1]['HUGO']] = point

            config = {
                'xlab': 'log2 FoldChange',
                'ylab': '-log10 p-value adjusted',

            }

            html_content = scatter.plot(data, config)

            self.add_section(
                name='Volcano',
                anchor='Volcano',
                content=html_content,
                description='Each dot at volcano plot copares corresponding gene fold change between tested groups of samples and evidence (p-value) in differential expression between these groups. Dots marked with green have both: p-value less than 0.1, and absolute log2 fold change greater than 1.'
            )

        def addBaseMeanPlot(data_key):
            data = {}

            for d in data_key.iterrows():
                point = {}

                point['x'] = d[1]['baseMean']
                point['y'] = d[1]['lfc']

                if d[1]['p'] < 0.1:
                    point['color'] = '#1f78b4'
                else:
                    point['color'] = '#b2df8a'

                data[d[1]['HUGO']] = point

            config = {
                'xLog': True,
                'xlab': 'log10 MeanCounts',
                'ylab': 'log2 FoldChange',

            }

            html_content = scatter.plot(data, config)
            self.add_section(
                name='MeanAverage',
                anchor='MeanAverage',
                content=html_content,
                description='Each dot at meanAverage plot copares corresponding gene mean expression across all samples and fold change between tested groups. Dots marked with green corresponds to genes with high evidence in differential expression between tested groups (-log10 p-value greater than 1)'
            )

        def addTopGenes(data_key):

            data_key = data_key.drop(columns=['Unnamed: 0', 'baseMean', 'gene_names', 'is_key'])

            data_key.sort_values(by='p', inplace=True, ascending=False)
            data_key['lfc'] = abs(data_key['lfc'])
            part = data_key[0:20]
            part = part.set_index('HUGO')

            headers = OrderedDict()
            headers['p'] = {'title': 'log10 p-value adjusted'}
            headers['lfc'] = {'title': 'log2 FoldChange'}

            table_config = {'col1_header': 'Gene Name'}

            data = part.to_dict(orient='index')


            html_content = table.plot(data, headers, table_config)

            self.add_section(
                name='Top genes',
                anchor='Top_genes',
                content=html_content,
                description='Top genes table shows top 20 genes with most weighted evidence in differential expression of these genes between tested groups'
            )

        def addHeatMap(HM_data):

            gene_names = HM_data['aux_sorted$dataHUGO'].tolist()
            HM_data.drop(columns=['Unnamed: 0', 'aux_sorted$dataHUGO'], inplace=True)

            val = HM_data.values

            lis = val.tolist()
            names = list(HM_data.columns.values)

            pconfig = {
                'xTitle': 'Sample Name',
                'yTitle': 'Gene id'
            }

            hm_html = heatmap.plot(lis, names, gene_names, pconfig)
            self.add_section(
                name='Heatmap',
                anchor='Heatmap',
                content=hm_html,
                description='This plot shows counts of differentially expressed genes on a per-sample basis. We have scaled the data by row'
            )

        addVolcanoPlot(data_key)
        addBaseMeanPlot(data_key)
        addTopGenes(data_key)

        # HM_data = pd.read_csv(hm_path)
        # addHeatMap(HM_data)