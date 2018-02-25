import os
import az
import subprocess
import pandas as pd
from ngs_reporting.rnaseq import gene_expression
from ngs_utils.file_utils import which, can_reuse, safe_mkdir, verify_file
from ngs_reporting.rnaseq.table_css import table_css_string
from os.path import join, dirname, isfile, abspath
from ngs_utils.logger import debug, info, critical, err
from ngs_utils.call_process import run
import numpy as np
from shutil import copyfile


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
    file_out.write(title + html_table_code + script1 + script2 + script3)
    file_out.close()


def transcript_summary(proj, key_gene_names):

    local_quant_path = ['sailfish/quant', 'salmon']

    for sam in proj.samples:

        for p in local_quant_path:
            quant_path = join(sam.dirpath, p, 'quant.sf')
            if isfile(quant_path):
                sam_qua = pd.read_table(quant_path)

        sam_qua = sam_qua.set_index('Name')

        if 'tpm' not in locals():
            tpm = pd.DataFrame(sam_qua['TPM'])
            sam_names = [sam.name]

        else:
            tpm = pd.concat([tpm, sam_qua['TPM']], axis=1)
            sam_names.append(sam.name)

    tpm.columns = sam_names

    tx2gene = pd.read_csv(proj.date_dir + '/name_by_transcript_id.csv')
    tx2gene = tx2gene.set_index('ensTx')
    tx2gene = tx2gene.to_dict()
    tx2gene = tx2gene['Gene']

    tpm['gene'] = 'NA'

    for i in tpm.index.tolist():
        if i in tx2gene:
            tpm.at[i, 'gene'] = tx2gene[i]

    tpm['is_key'] = 'False'

    for i in tpm.index.tolist():
        if tpm.at[i,'gene'] in key_gene_names:
            tpm.at[i, 'is_key'] = 'True'

    tpm_key = tpm.loc[tpm['is_key'] == 'True']

    tx_lvl_html_path = proj.expression_dir + '/html/isoform.html'
    transcript_level_html(tpm_key, tx_lvl_html_path)

    tpm_key.to_csv(join(proj.final_dir, 'isoforms.tpm'), sep='\t')


def exon_level_html(cnt, exLvl_html_path):

    def fl_for(x):
        return "%.3f" % x

    cnt = cnt.drop(['is_key', 'GeneID'], axis = 1)
    # cnt = cnt.rename(index=str,columns={'gene': 'Gene name'})
    # reorder columns
    ex_cnt = cnt[['gene', 'start', 'end']]
    cnt = cnt.drop(['gene', 'start', 'end'], axis = 1)
    ex_cnt = pd.concat([ex_cnt, cnt], axis = 1)

    html_table_code = ex_cnt.to_html(float_format=fl_for, border=0, justify='left', index_names=False, index=False)
    table_id = 'TxLevel'
    html_table_code = html_table_code.replace('<table border="0" class="dataframe">','<table id="' + table_id + '" class="display">')
    title = '<h3>Exon level</h3><p>Shown only key genes</p>'
    # jquery scripts
    style = table_css_string
    script1 = '<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jQuery/jquery-1.8.2.min.js"></script>'
    script2 = '<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/jquery.dataTables.min.js"></script>'
    script3 = '<script> $(function(){$("#' + table_id + '").dataTable({"iDisplayLength": 50}); })</script>'
    # write combined html code
    file_out = open(exLvl_html_path, 'w')
    file_out.write(title + style + html_table_code + script1 + script2 + script3)
    file_out.close()


def create_exon_counts_file(proj, key_gene_names):

    transcripts_file = az.get_refdata(proj.genome_build)['all_transcripts']

    for sam in proj.samples:
        bam_name = sam.name + '-ready.bam'
        bam_path = join(sam.dirpath, bam_name)

        safe_mkdir(proj.dir + '/work/postproc')
        out_file = join(proj.work_dir, sam.name + '_exon_counts.csv')

        if not can_reuse(out_file, bam_name):
            pathEXC_R = dirname(__file__) + '/exon_counts.R'
            cmdl = ' '.join(['Rscript', pathEXC_R, bam_path, out_file, transcripts_file])
            print(cmdl)
            os.system(cmdl)

        ex_cnt = pd.read_csv(out_file, skiprows=0, header=0)
        # ex_cnt = ex_cnt.set_index('GeneID')

        if 'cnt' not in locals():
            cnt = ex_cnt
            cnt = cnt.drop('Unnamed: 0', axis=1)
        else:
            cnt = pd.concat([cnt, ex_cnt['counts']], axis=1)

        cnt = cnt.rename(columns={'counts': sam.name})

    # add gene name
    id2name = pd.read_csv(proj.date_dir + '/name_by_gene_id.csv')
    id2name = id2name.set_index('ensId')
    id2name = id2name.to_dict()
    id2name = id2name['Name']

    cnt['gene'] = 'NA'

    for i in cnt.index.tolist():
        cnt.at[i, 'gene'] = id2name[cnt.at[i, 'GeneID']]

    # add is_key
    cnt['is_key'] = 'False'

    for i in cnt.index.tolist():
        if cnt.at[i, 'gene'] in key_gene_names:
            cnt.at[i, 'is_key'] = 'True'

    cnt_key = cnt.loc[cnt['is_key'] == 'True']

    # make html-table
    ex_lvl_html_path = proj.expression_dir + '/html/exon.html'
    exon_level_html(cnt_key, ex_lvl_html_path)


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

    # prepare file-links

    fnames = ['rawCounts.csv', 'normalizedCounts.csv', 'rlog.csv', 'vst.csv', \
              'gene.est.csv', 'gene.final.csv', 'gene.fitted.csv', 'corMatrix.csv']
    fnames = [join(out_file, f) for f in fnames]

    fnames.append(join(proj.date_dir, 'combined.counts'))

    return fnames


def run_FA(fa_in, fa_out):
    pathRScript = which('Rscript')
    pathFA_R = dirname(__file__) + '/FA.R'
    if not can_reuse(fa_out, fa_in):
        cmd = ' '.join([pathRScript, pathFA_R, fa_in, fa_out])
        print('Running:')
        print(cmd)
        # os.system(cmd)

    return join(fa_out, 'RNA_PW.csv')


def run_analysis(proj, key_gene_names):
    info('*' * 70)
    info('running RNA analysis')

    # gene_expression.make_heatmaps(proj, key_gene_names)

    # transcript_summary(proj, key_gene_names)

    # create_exon_counts_file(proj, key_gene_names)

    if not os.path.isfile(proj.expression_dir + '/combined.counts'):
        annotateGeneCounts(proj, key_gene_names)

    rna_files_list = []

    if not isfile(join(proj.date_dir, 'project-summary.yaml')):
        copyfile(join(proj.log_dir, 'project-summary.yaml'), join(proj.date_dir, 'project-summary.yaml'))
    if not isfile(join(proj.date_dir, 'combined.counts')):
        copyfile(join(proj.expression_dir, 'combined.counts'), join(proj.date_dir, 'combined.counts'))

    safe_mkdir(proj.dir + '/work/postproc')
    de_out = proj.work_dir + '/RNA_DE.csv'
    hm_out = proj.work_dir + '/RNA_HM.csv'
    run_DE(proj, de_out, hm_out)
    rna_files_list.extend([de_out, hm_out])

    full_table_html_path = proj.expression_dir + '/html/diff_exp.html'
    make_full_expreesion_table(de_out, full_table_html_path)
    proj.full_expression_dir = full_table_html_path

    qc_out_dir = safe_mkdir(proj.work_dir + '/RNA_QC')
    qc_out_files = run_QC(proj, qc_out_dir)

    rna_files_list.extend(qc_out_files)

    fa_in = de_out
    fa_out = proj.work_dir + '/'
    fa_file = run_FA(fa_in, fa_out)
    rna_files_list.extend([fa_file])



    for p in rna_files_list:
        proj.postproc_mqc_files.append(p)
