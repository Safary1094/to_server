#!/usr/bin/env python

""" MultiQC module to add link to Bcl2fastq reports """

import logging
from os.path import join, dirname, abspath
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

    def makeVolcanoData(self, data_key):

        data = {}

        for d in data_key.iterrows():

            point = {}
            point['x'] = d[1]['lfc']
            point['y'] = d[1]['p']

            if abs(d[1]['p']) > 1 and abs(d[1]['lfc']) > 1:
                point['color'] = '#b2df8a'
            else:
                point['color'] = '#1f78b4'

            data[d[1]['gene']] = point

        return data

    def addVolcano(self, de):
        data = []
        contrast = []
        for cont in de:
            data.append(self.makeVolcanoData(de[cont]))
            contrast.append({'name': cont})

        config = {'data_labels': contrast}

        self.add_section(
            name='Volcano',
            anchor='Volcano',
            content=scatter.plot(data, config)
        )

    def addNumDE_perContrast(self, de):
        de_num = {}
        p_val = 0.1
        logFold = 2
        for cont in de:
            curr_de = de[cont]
            num = len(curr_de.loc[(curr_de['p'] > 1) & (abs(curr_de['lfc']) > 0.5)])
            de_num[cont] = {'numbef_of_de_genes': num}


        print(de_num)

        headers = OrderedDict()
        headers['Sample Name'] = {'title': 'Contrast'}
        headers['numbef_of_de_genes'] = {'title': 'number of DE genes'}

        self.add_section(
            name='DEs per contrast',
            anchor='DEs percontrast',

            content=table.plot(de_num, headers, {'col1_header': 'Contrast'})
        )

        return de_num

    def addDE_overlap(self, de):

        de_overlap = {}

        for i in de:
            de_overlap[i] = {}
            for j in de:
                de_i = de[i].loc[(de[i]['p'] > 1) & (abs(de[i]['lfc']) > 0.5)]
                de_j = de[j].loc[(de[j]['p'] > 1) & (abs(de[j]['lfc']) > 0.5)]
                i_de_genes = set(de_i.index.tolist())
                j_de_genes = set(de_j.index.tolist())

                de_overlap[i][j] = len(i_de_genes & j_de_genes)

        de_overlap = pd.DataFrame(de_overlap)

        hmdata = de_overlap.values.tolist()
        names = de_overlap.index.tolist()

        hm_html = heatmap.plot(hmdata, names)

        self.add_section(
            name='DE overlap',
            anchor='DE overlap',
            content=hm_html,
            description='Table of numbers of overlapping genes across contrasts'
        )

    def addTopGenes(self, de):

        tab_header_volcano = '<ul>'
        tab_content_volcano = '<div>'

        for c in de:
            de[c] = de[c].drop('gene_id', axis=1)

            de[c].sort_values(by='p', inplace=True, ascending=False)
            top = de[c][0:20]
            top = top.set_index('gene')

            headers = OrderedDict()
            headers['p'] = {'title': 'log10 p-value adjusted'}
            headers['lfc'] = {'title': 'log2 FoldChange'}

            tab_header_volcano += '<li>' + c + '</li>'

            table_config = {'col1_header': 'Gene Name'}
            mqc_table = table.plot(top.to_dict(orient='index'), headers, table_config)
            tab_content_volcano += '<div>' + mqc_table + '</div>'

        tab_header_volcano += '</ul>'
        tab_content_volcano += '</div>'

        html_table = '<div class="tabs">' + tab_header_volcano + tab_content_volcano + '</div>'

        style_file = open(join(dirname(abspath(__file__)), "style.txt"), 'r')
        style = style_file.read()

        script_file = open(join(dirname(abspath(__file__)), "js.txt"), 'r')
        script = script_file.read()

        self.add_section(
            name='Top DE genes',
            anchor='Top DE genes',
            content=style +'\n'+ script+'\n' + html_table+'\n' + '<script> $(document).ready(function(){ $(".tabs").lightTabs(); }); </script>',
        )

    def addBaseMeanPlot(self, ata_key):
        data1 = {
            'sample 1': {
                'x': 1,
                'y': 1
            }}
        data2 = {'sample 2': {
            'x': 2,
            'y': 3
        }
        }

        config = {
            'data_labels': [
                {'name': 'rlog', 'ylab': 'rlog'},
                {'name': 'vst', 'ylab': 'vst'}]}

        html_content = scatter.plot([data1, data2], config)

        self.add_section(
            name='MeanAverage',
            anchor='MeanAverage',
            content=html_content,
            description='Each dot at meanAverage plot copares corresponding gene mean expression across all samples and fold change between tested groups. Dots marked with green corresponds to genes with high evidence in differential expression between tested groups (-log10 p-value greater than 1)'
        )

    def __init__(self):
        mod_name = 'RNAseqDE'
        super(MultiqcModule, self).__init__(name='RNA Differential Expression', anchor=mod_name)
        # make dict of de-tables per contrast
        de = {}
        for f in self.find_log_files(mod_name, filecontents=False):
            print(f)
            dirpath, fname = f['root'], f['fn']
            if f['s_name'] == 'de_gene_key':
                de_path = join(dirpath, fname)
                contrast = f['root'].split('/')
                data = pd.read_csv(de_path)
                de[contrast[-1]] = data

        if len(de)==0:
            print('no DE files')
        else:
            self.addVolcano(de)
            self.addNumDE_perContrast(de)
            self.addDE_overlap(de)
            self.addTopGenes(de)



        #
        #
        #
        # def addHeatMap(HM_data):
        #
        #     gene_names = HM_data['aux_sorted$dataHUGO'].tolist()
        #     HM_data.drop(columns=['Unnamed: 0', 'aux_sorted$dataHUGO'], inplace=True)
        #
        #     val = HM_data.values
        #
        #     lis = val.tolist()
        #     names = list(HM_data.columns.values)
        #
        #     pconfig = {
        #         'xTitle': 'Sample Name',
        #         'yTitle': 'Gene id'
        #     }
        #
        #     hm_html = heatmap.plot(lis, names, gene_names, pconfig)
        #     self.add_section(
        #         name='Heatmap',
        #         anchor='Heatmap',
        #         content=hm_html,
        #         description='This plot shows counts of differentially expressed genes on a per-sample basis. We have scaled the data by row'
        #     )
        #
        # # addVolcanoPlot(data)

        # addBaseMeanPlot(data)
        # addTopGenes(data)
        #
        # # HM_data = pd.read_csv(hm_path)
        # # addHeatMap(HM_data)