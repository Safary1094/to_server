import os
import pandas as pd
from ngs_reporting.rnaseq import gene_expression
from ngs_utils.file_utils import which, can_reuse, safe_mkdir
from ngs_reporting.rnaseq.table_css import table_css_string
from os.path import join, dirname
from ngs_utils.logger import debug, info
import numpy as np

def annotateGeneCounts(proj, key_genes):
    genes = pd.read_csv(proj.raw_expression_dir + '/combined.counts', sep = '\t')
    # set gene_by_transcriptid and gene_by_geneid
    gene_by_transcriptid, gene_by_geneid = gene_expression._parse_gene_transcripts_id(proj.genome_build, '')

    for i in range(len(genes)):
        gene_id = genes.at[i,'id']
        if gene_id in gene_by_transcriptid:
            genes.at[i, 'HUGO'] = gene_by_transcriptid[gene_id]
        elif gene_id in gene_by_geneid:
            genes.at[i, 'HUGO'] = gene_by_geneid[gene_id]

    genes.dropna(inplace=True)
    genes['index'] = range(len(genes))
    genes.set_index('index', inplace=True)

    for i in range(len(genes)):
        gene_name= genes.at[i,'HUGO']
        if gene_name in key_genes:
            genes.at[i, 'is_key'] = True
        else:
            genes.at[i, 'is_key'] = False

    genes.sort_values(by='id', inplace=True)
    # for postproc
    genes.to_csv(proj.expression_dir + '/combined.counts', index = False, sep = '\t')
    # for RNAseq
    genes.to_csv(proj.date_dir + '/combined.counts', index = False, sep = '\t')


def make_full_expreesion_table(de_path, full_table_html_path):
    if not can_reuse(full_table_html_path, de_path):
        data = pd.read_csv(de_path)
        data = data.dropna()
        # make full external table
        data = data[['HUGO', 'gene_names', 'p', 'lfc']]
        data = data.rename(index=str, columns={'HUGO': 'HUGO name', 'gene_names': 'Gene index', 'p': 'Log10 p value adjusted', 'lfc': 'Log2 fold change'})
        html_table_code = data.to_html(border=0, justify='left', index_names=False, index=False)
        # add table id
        table_id = 'diff_exp'
        html_table_code = html_table_code.replace('<table border="0" class="dataframe">', '<table id="' + table_id + '" class="display">')
        title = '<h3>Differentially expressed genes</h3><p>Shown only annotated genes</p>'
        # jquery scripts
        style = table_css_string
        script1 = '<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jQuery/jquery-1.8.2.min.js"></script>'
        script2 = '<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/jquery.dataTables.min.js"></script>'
        script3 = '<script> $(function(){$("#' + table_id + '").dataTable({"iDisplayLength": 50}); })</script>'
        # write combined html code
        file_out = open(full_table_html_path, 'w')
        file_out.write(title + style + html_table_code + script1 + script2 + script3)
        file_out.close()

    return full_table_html_path


def transcript_level_html(tpm, TxLvl_html_path):

    def fl_for(x):
        return "%.3f" % x

    tpm['Transcript index'] = tpm.index
    tpm = tpm.drop(['is_key'], axis = 1)
    tpm = tpm.rename(index=str,columns={'gene': 'Gene name'})
    # reorder columns
    gene_tx = tpm[['Gene name', 'Transcript index']]
    tpm = tpm.drop(['Gene name', 'Transcript index'], axis = 1)
    tpm = pd.concat([gene_tx, tpm], axis = 1)

    html_table_code = tpm.to_html(float_format=fl_for, border=0, justify='left', index_names=False, index=False)
    table_id = 'TxLevel'
    html_table_code = html_table_code.replace('<table border="0" class="dataframe">','<table id="' + table_id + '" class="display">')
    title = '<h3>Isoform level</h3><p>Shown only key genes</p>'
    # jquery scripts
    style = table_css_string
    script1 = '<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jQuery/jquery-1.8.2.min.js"></script>'
    script2 = '<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/jquery.dataTables.min.js"></script>'
    script3 = '<script> $(function(){$("#' + table_id + '").dataTable({"iDisplayLength": 50}); })</script>'
    # write combined html code
    file_out = open(TxLvl_html_path, 'w')
    file_out.write(title + style + html_table_code + script1 + script2 + script3)
    file_out.close()


def transcript_summary(proj, key_gene_names):

    for sam in proj.samples:

        quant_path = join(sam.dirpath, 'salmon', 'quant.sf')
        sam_qua = pd.read_table(quant_path)

        sam_qua = sam_qua.set_index('Name')

        if 'tpm' not in locals():
            tpm = pd.DataFrame(sam_qua['TPM'])
            sam_names = [sam.name]

        else:
            tpm = pd.concat([tpm, sam_qua['TPM']], axis=1)
            sam_names.append(sam.name)

    tpm.columns = sam_names

    tx2gene = pd.read_csv(proj.date_dir + '/tx2gene_name.csv')
    tx2gene = tx2gene.set_index('ensTx')
    tx2gene = tx2gene.to_dict()
    tx2gene = tx2gene['Gene']

    tpm['gene'] = 'NA'

    for i in tpm.index.tolist():
        tpm.at[i, 'gene'] = tx2gene[i]

    tpm['is_key'] = 'False'

    for i in tpm.index.tolist():
        if tpm.at[i,'gene'] in key_gene_names:
            tpm.at[i, 'is_key'] = 'True'

    tpm_key = tpm.loc[tpm['is_key'] == 'True']

    tx_lvl_html_path = proj.expression_dir + '/html/isoform.sf.tpm.html'
    transcript_level_html(tpm_key, tx_lvl_html_path)

    1


def run_DE(proj, de_out, hm_out):
    pathRScript = which('Rscript')
    pathDE_R = dirname(__file__) + '/DE.R'
    comp_file_for_reuse = proj.date_dir + '/combined.counts'
    if not can_reuse(de_out, comp_file_for_reuse):
        cmd = ' '.join([pathRScript, pathDE_R, proj.final_dir, de_out, hm_out])
        print('Running:')
        print(cmd)
        # os.system(cmd)


def run_QC(proj, out_file):
    pathRScript = which('Rscript')
    pathQC_R = dirname(__file__) + '/QC.R'
    comp_file_for_reuse = proj.date_dir + '/combined.counts'
    if not can_reuse(out_file, comp_file_for_reuse):
        cmd = ' '.join([pathRScript, pathQC_R, proj.final_dir, out_file])
        print('Running:')
        print(cmd)
        # os.system(cmd)

    return [join(out_file, f) for f in os.listdir(out_file)[1:]]



def run_FA(fa_in, fa_out):
    pathRScript = which('Rscript')
    pathFA_R = dirname(__file__) + '/FA.R'
    if not can_reuse(fa_out, fa_in):
        cmd = ' '.join([pathRScript, pathFA_R, fa_in, fa_out])
        print('Running:')
        print(cmd)
        # os.system(cmd)

    return join(fa_out, 'RNA_PW.csv')

def run(proj, key_gene_names):
    info('*' * 70)
    info('running RNA analysis')

    transcript_summary(proj, key_gene_names)

    # if not os.path.isfile(proj.expression_dir + '/combined.counts'):
    #     annotateGeneCounts(proj, key_gene_names)

    # rna_files_list = []
    #
    # safe_mkdir(proj.dir + '/work/postproc')
    # de_out = proj.work_dir + '/RNA_DE.csv'
    # hm_out = proj.work_dir + '/RNA_HM.csv'
    # run_DE(proj, de_out, hm_out)
    # rna_files_list.extend([de_out, hm_out])
    #
    # full_table_html_path = proj.expression_dir + '/html/diff_exp.html'
    # make_full_expreesion_table(de_out, full_table_html_path)
    # proj.full_expression_dir = full_table_html_path
    #
    # qc_out_dir = safe_mkdir(proj.work_dir + '/RNA_QC')
    # qc_out_files = run_QC(proj, qc_out_dir)
    #
    # rna_files_list.extend(qc_out_files)
    #
    # fa_in = de_out
    # fa_out = proj.work_dir + '/'
    # fa_file = run_FA(fa_in, fa_out)
    # rna_files_list.extend([fa_file])

    gene_expression.make_heatmaps(proj, key_gene_names)

    # for p in rna_files_list:
    #     proj.postproc_mqc_files.append(p)
