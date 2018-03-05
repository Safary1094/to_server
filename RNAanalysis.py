import os
import az
import subprocess
import pandas as pd
from ngs_utils.file_utils import which, can_reuse, safe_mkdir, verify_file
from ngs_reporting.rnaseq.table_to_html import table_to_html
from os.path import join, dirname, isfile, abspath
from ngs_utils.logger import debug, info, critical, err
from ngs_utils.call_process import run
import numpy as np
from shutil import copyfile


def load_id2name(proj):
    id2name = pd.read_csv(join(proj.local_ref_data, proj.genome_build, 'name_by_gene_id.csv'))
    id2name = id2name.set_index('ensId')
    id2name = id2name.to_dict()
    return id2name['Name']


def load_tx2name(proj):
    tx2gene = pd.read_csv(join(proj.local_ref_data, proj.genome_build, 'name_by_transcript_id.csv'))
    tx2gene = tx2gene.set_index('ensTx')
    tx2gene = tx2gene.to_dict()
    return tx2gene['Gene']


def run_tximport(proj):
    # computes, gene_counts, gene_tpms, isoform tpms from quant.sf
    pathEXC_R = dirname(__file__) + '/run_tximport.R'
    cmdl = ['Rscript', pathEXC_R]

    local_quant_path = ['sailfish/quant', 'salmon']
    for sam in proj.samples:
        for p in local_quant_path:
            quant_path = join(sam.dirpath, p, 'quant.sf')
            if isfile(quant_path):
                cmdl.append(quant_path)

    cmdl.append(proj.expression_dir)
    cmdl = ' '.join(cmdl)
    print(cmdl)
    os.system(cmdl)


def run_feature_counts(proj):
    # computes, exon_counts from .bam
    transcripts_file = az.get_refdata(proj.genome_build)['exon_annotation']
    pathEXC_R = dirname(__file__) + '/exon_counts.R'
    cmdl = ['Rscript', pathEXC_R]
    for sam in proj.samples:
        cmdl.append(join(sam.dirpath, sam.name + '-ready.bam'))
    cmdl.append(proj.expression_dir)
    cmdl.append(transcripts_file)

    cmdl = ' '.join(cmdl)
    print(cmdl)
    os.system(cmdl)


def calculate_expression_levels(proj):

    run_tximport(proj)
    run_feature_counts(proj)


def exon_level_html(proj, key_gene_names):

    exon = pd.read_csv(join(proj.expression_dir, 'exon_counts.csv'), index_col=0)

    # add gene names
    id2gene = load_id2name(proj)
    exon['gene'] = 'NA'
    for i in exon.index.tolist():
        if i in id2gene:
            exon.at[i, 'gene'] = id2gene[i]
        else:
            print('Warning! gene index ' + str(i) + ' not found')

    # rewrite annotated isoforms
    full_link = join(proj.expression_dir, 'exon_counts.csv')
    exon.to_csv(full_link)

    # select key genes
    exon_key = exon.loc[exon['gene'].isin(key_gene_names)]

    # save html for key-genes
    gradient_cols = [sam.name for sam in proj.samples]
    title = '<h3>Exon level</h3>Showing only key genes. \
    The full results can be downloaded from here: \
    <a href="'+ full_link + '">exon_counts.csv</a>'

    html_path = join(proj.expression_dir, 'html', 'exon_counts.html')
    table_to_html(exon_key, gradient_cols, title, html_path)
    proj.counts_names.append(html_path)


def isoform_level_html(proj, key_gene_names):

    tpm = pd.read_csv(join(proj.expression_dir, 'isoform_tpm.csv'), index_col=0)

    # add gene names
    tx2gene = load_tx2name(proj)
    tpm['gene'] = 'NA'
    for i in tpm.index.tolist():
        if i in tx2gene:
            tpm.at[i, 'gene'] = tx2gene[i]
        else:
            print('Warning! transcript index ' + str(i) + ' not found')

    # rewrite annotated isoforms
    full_link = join(proj.expression_dir, 'isoform_tpm.csv')
    tpm.to_csv(full_link)

    # select key genes
    tpm_key = tpm.loc[tpm['gene'].isin(key_gene_names)]

    # save html for key-genes
    gradient_cols = [sam.name for sam in proj.samples]
    title = '<h3>Isoform level</h3>Showing only key genes. \
    The full results can be downloaded from here: \
    <a href="'+ full_link + '">isoform_tpm.csv</a>'

    html_path = join(proj.expression_dir, 'html', 'isoforms_tpm.html')
    table_to_html(tpm_key, gradient_cols, title, html_path)
    proj.counts_names.append(html_path)


def gene_counts_html(proj, key_gene_names):

    gcounts = pd.read_csv(join(proj.expression_dir, 'gene_counts.csv'), index_col=0)

    # add gene names
    id2gene = load_id2name(proj)
    gcounts['gene'] = 'NA'
    for i in gcounts.index.tolist():
        if i in id2gene:
            gcounts.at[i, 'gene'] = id2gene[i]
        else:
            print('Warning! gene index ' + str(i) + ' not found')

    # rewrite annotated isoforms
    full_link = join(proj.expression_dir, 'gene_counts.csv')
    gcounts.to_csv(full_link)

    # select key genes
    gcounts_key = gcounts.loc[gcounts['gene'].isin(key_gene_names)]

    # save html for key-genes
    gradient_cols = [sam.name for sam in proj.samples]
    title = '<h3>Gene counts</h3>Showing only key genes. \
    The full results can be downloaded from here: \
    <a href="'+ full_link + '">gene_counts.csv</a>'

    html_path = join(proj.expression_dir, 'html', 'gene_counts.html')
    table_to_html(gcounts_key, gradient_cols, title, html_path)
    proj.counts_names.append(html_path)


def gene_tpm_html(proj, key_gene_names):

    gene_tpm = pd.read_csv(join(proj.expression_dir, 'gene_tpm.csv'), index_col=0)

    # add gene names
    id2gene = load_id2name(proj)
    gene_tpm['gene'] = 'NA'
    for i in gene_tpm.index.tolist():
        if i in id2gene:
            gene_tpm.at[i, 'gene'] = id2gene[i]
        else:
            print('Warning! gene index ' + str(i) + ' not found')

    # rewrite annotated isoforms
    full_link = join(proj.expression_dir, 'gene_tpm.csv')
    gene_tpm.to_csv(full_link)

    # select key genes
    gene_tpm_key = gene_tpm.loc[gene_tpm['gene'].isin(key_gene_names)]

    # save html for key-genes
    gradient_cols = [sam.name for sam in proj.samples]
    title = '<h3>Gene TPM</h3>Showing only key genes. \
    The full results can be downloaded from here: \
    <a href="'+ full_link + '">gene_tpm.csv</a>'

    html_path = join(proj.expression_dir, 'html', 'gene_tpm.html')
    table_to_html(gene_tpm_key, gradient_cols, title, html_path)
    proj.counts_names.append(html_path)


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
        os.system(cmd)

    # prepare file-links

    fnames = ['rawCounts.csv', 'normalizedCounts.csv', 'rlog.csv', 'vst.csv', \
              'gene.est.csv', 'gene.final.csv', 'gene.fitted.csv', 'corMatrix.csv', 'pca.csv']
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


def diff_exp_genes_html(de_out, proj, key_gene_names):
    de_gene = pd.read_csv(de_out, index_col=0)

    # add gene names
    id2gene = load_id2name(proj)
    de_gene['gene'] = 'NA'
    for i in de_gene.index.tolist():
        if i in id2gene:
            de_gene.at[i, 'gene'] = id2gene[i]
        else:
            print('Warning! gene index ' + str(i) + ' not found')

    # rewrite annotated isoforms
    full_link = join(proj.expression_dir, 'de_gene.csv')
    de_gene.to_csv(full_link)

    # select key genes
    de_gene_key = de_gene.loc[de_gene['gene'].isin(key_gene_names)]

    # save html for key-genes
    gradient_cols = [sam.name for sam in proj.samples]
    title = '<h3>DE genes</h3>Showing only key genes. \
    The full results can be downloaded from here: \
    <a href="'+ full_link + '">de_genes.csv</a>'

    html_path = join(proj.expression_dir, 'html', 'de_genes.html')
    table_to_html(de_gene_key, gradient_cols, title, html_path)
    proj.counts_names.append(html_path)


def run_analysis(proj, key_gene_names):
    info('*' * 70)
    info('running RNA analysis')

    # expression levels
    #calculate_expression_levels(proj)

    isoform_level_html(proj, key_gene_names)
    # exon_level_html(proj, key_gene_names)
    gene_counts_html(proj, key_gene_names)
    gene_tpm_html(proj, key_gene_names)

    # DE analysis
    # rna_files_list = []
    #
    # if not isfile(join(proj.date_dir, 'project-summary.yaml')):
    #     copyfile(join(proj.log_dir, 'project-summary.yaml'), join(proj.date_dir, 'project-summary.yaml'))
    # if not isfile(join(proj.date_dir, 'combined.counts')):
    #     copyfile(join(proj.expression_dir, 'combined.counts'), join(proj.date_dir, 'combined.counts'))
    #
    # safe_mkdir(proj.dir + '/work/postproc')
    # de_out = proj.work_dir + '/RNA_DE.csv'
    # hm_out = proj.work_dir + '/RNA_HM.csv'
    # run_DE(proj, de_out, hm_out)
    # rna_files_list.extend([de_out, hm_out])
    #
    # full_table_html_path = proj.expression_dir + '/html/diff_exp.html'
    # diff_exp_genes_html(de_out, proj, key_gene_names)
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
    #
    # for p in rna_files_list:
    #     proj.postproc_mqc_files.append(p)

