import os
import pandas as pd
from ngs_reporting.rnaseq import gene_expression
from ngs_utils.file_utils import which, can_reuse, safe_mkdir
from ngs_reporting.rnaseq.table_css import table_css_string
from os.path import join, dirname
from ngs_utils.logger import debug, info

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
    if not os.path.isfile(proj.expression_dir + '/combined.counts'):
        annotateGeneCounts(proj, key_gene_names)

    rna_files_list = []

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
    #gene_expression.make_heatmaps(proj, key_gene_names)

    for p in rna_files_list:
        proj.postproc_mqc_files.append(p)
