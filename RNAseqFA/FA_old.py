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
from ngs_reporting.rnaseq.table_css import table_css_string
import xml.etree.ElementTree as ET



class MultiqcModule(BaseMultiqcModule):


    def __init__(self):
        mod_name = 'DE_RNAseq'
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='RNA Differential Expression', anchor=mod_name)

        rnaseq_pca_files = self.find_log_files(mod_name, filecontents=False)
        for f in rnaseq_pca_files:
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
            print(html_content)

            self.add_section(
                name='Volcano',
                anchor='Volcano',
                content=html_content,
                description='Each dot at volcano plot copares corresponding gene fold change between tested groups of samples and evidence (p-value) in differential expression between these groups. Dots marked with green have both: -log10 p-value greater than 1, and absolute log2 fold change greater than 1.'
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
                description='This plot shows only differentially expressed genes on a per-sample basis. We have scaled the data by row'
            )

        addVolcanoPlot(data_key)
        addBaseMeanPlot(data_key)
        addTopGenes(data_key)

        HM_data = pd.read_csv(hm_path)
        addHeatMap(HM_data)


        # make full external table
        data_parsed= data_parsed[['HUGO', 'gene_names', 'p', 'lfc']]
        data_parsed = data_parsed.rename(index=str, columns={'HUGO':'HUGO name', 'gene_names':'Gene index', 'p':'Log10 p value adjusted', 'lfc':'Log2 fold change'})
        html_table_code = data_parsed.to_html(border=0, justify='left', index_names=False, index=False)

        # add table id
        table_id = 'diff_exp'
        html_table_code = html_table_code.replace('<table border="0" class="dataframe">',
                                                  '<table id="' + table_id + '" class="display">')

        title = '<h3>Differentially expressed genes</h3><p>Shown only annotated genes</p>'

        # jquery scripts
        style = table_css_string
        script1 = '<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jQuery/jquery-1.8.2.min.js"></script>'
        script2 = '<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/jquery.dataTables.min.js"></script>'
        script3 = '<script> $(function(){$("#' + table_id + '").dataTable({"iDisplayLength": 50}); })</script>'

        # write combined html code
        file_out = open('/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/2017-09-21_bcbio_rnaseq/DiffExp/diff.html', 'w')

        file_out.write(title + style + html_table_code + script1 + script2 + script3)

        # pathviewxperimental

        xml_path = '/home/alexey/hsa04015.xml'

        tree = ET.parse(xml_path)

        root = tree.getroot()

        elem = 'elements: ['

        for el in root:
            if el.attrib['type'] == 'gene':
                el_cont = iter(el).__next__()

                node_id = el.attrib['id']
                x = el_cont.attrib['x']
                y = el_cont.attrib['y']

                elem += '{ ' + \
                        "classes: 'gene'," + \
                        'data: {id:' + node_id + '},' + \
                        'position:{' + 'x:' + x + ',y:' + y + '}' + \
                                                              '},'

            if el.attrib['type'] == 'compound':
                el_cont = iter(el).__next__()

                node_id = el.attrib['id']
                x = el_cont.attrib['x']
                y = el_cont.attrib['y']

                elem += '{ '  + \
                        'classes: "compound",' + \
                        'data: {id:' + node_id + '},' + \
                        'position:{' + 'x:' + x + ',y:' + y + '}' + \
                                                              '},'

            if el.tag == 'relation':
                b = el.attrib['entry1']
                e = el.attrib['entry2']
                elem += '{ ' + \
                        'classes: "relation",' + \
                        'data: {id: ' + b + e + ', source: ' + b + ',' + 'target: ' + e + '}' + \
                        '},'

        elem = elem[:-1]
        elem += ']'
        cyto_scape_path = "<script src='https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.2.7/cytoscape.js'></script>"
        div_style = "<style> #cy { width: 100%; height: 900px; position: relative; top: 0px; left: 0px; } </style>"
        layout = "layout: {    name: 'preset'  }"
        cy_style = "style: [" \
                   "{ "+\
                   "selector: '.gene',"+\
                   "style: {"+\
                   "'content': 'data(id)',"+ \
                   "'width': 46," + \
                   "'height': 17," + \
                   "'border-width': 0.5," + \
                   "'border-style': 'solid'," + \
                   "'border-color': '#333'," + \
                   "'content': 'data(id)'," + \
                   "'shape': 'rectangle',"+ \
                   "'text-valign': 'center'," + \
                   "'background-color': '#bfffbf'" + \
                   "}}," \
                   "{ " + \
                   "selector: '.compound'," + \
                   "style: {" + \
                   "'height': 25," + \
                   "'width': 25," + \
                   "'border-width': 0.5," + \
                   "'border-style': 'solid'," + \
                   "'border-color': '#333'," + \
                   "'content': 'data(id)'," + \
                   "'text-valign': 'center'," + \
                   "'background-color': '#ddd'" + \
                   "}},"+ \
                   "{" +\
                   "selector: '.relation',"+\
                   "style: {"+\
                   "'width': 1,"+ \
                   "'curve-style': 'bezier',"+\
                   "'line-color': '#aaa',"+\
                   "'target-arrow-color': '#aaa',"+\
                   "'target-arrow-shape': 'triangle'"+\
                   "}}" +\
                   "]"

        self.add_section(
            name='Pathview',
            anchor='Pathview',
            content="<div id='table_qonl_container' class='mqc_table_container'>"+
                    cyto_scape_path + div_style +
                    " <div id='cy'></div> <script> var cy = cytoscape({ container: document.getElementById('cy'), " + elem + "," + layout +","+cy_style+" }); </script>",
            description='Metabolic path: hsa4015'
        ) 
