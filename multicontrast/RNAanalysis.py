import os
import az
import subprocess
import pandas as pd
from ngs_utils.file_utils import which, can_reuse, safe_mkdir, verify_file
from ngs_reporting.rnaseq.table_to_html import table_to_html
from os.path import join, dirname, isfile, abspath
from ngs_utils.logger import debug, info, critical, err
from shutil import copyfile
import yaml


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
    cmdl.append(join(proj.date_dir, 'tx2gene.csv'))
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
    print('exon_level_html')
    exon = pd.read_csv(join(proj.expression_dir, 'exon_counts.csv'), index_col=0)

    # add gene names
    if 'gene' not in list(exon):
        id2gene = load_id2name(proj)
        exon_index = exon.index.tolist()
        gene_index = []
        for i in exon_index:
            if i in id2gene:
                gene_index.append(id2gene[i])
            else:
                print('Warning! transcript index ' + str(i) + ' not found')
                gene_index.append('NA')
        se = pd.Series(gene_index)
        exon.insert(0, 'gene', se.values)

    # rewrite annotated isoforms
    full_link = join(proj.expression_dir, 'exon_counts.csv')
    exon.to_csv(full_link)

    # select key genes
    exon_key = exon.loc[exon['gene'].isin(key_gene_names)]
    # drop index, to make it unique
    exon_key = exon_key.reset_index()

    # save html for key-genes
    gradient_cols = [sam.name for sam in proj.samples]
    title = '<h3>Exon level</h3>Showing only key genes. \
    The full results can be downloaded from here: \
    <a href="'+ full_link + '">exon_counts.csv</a>'

    html_path = join(proj.expression_dir, 'html', 'exon_counts.html')
    table_to_html(exon_key, gradient_cols, title, html_path)
    proj.counts_names.append(html_path)


def isoform_level_html(proj, key_gene_names):
    print('isoform_level_html')

    tpm = pd.read_csv(join(proj.expression_dir, 'isoform_tpm.csv'), index_col=0)

    # add gene names
    if 'gene' not in list(tpm):
        tx2gene = load_tx2name(proj)

        isof_index = tpm.index.tolist()
        gene_index = []
        for i in isof_index:
            if i in tx2gene:
                gene_index.append(tx2gene[i])
            else:
                print('Warning! transcript index ' + str(i) + ' not found')
                gene_index.append('NA')
        se = pd.Series(gene_index)
        tpm.insert(0, 'gene', se.values)

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
    print('gene_counts_html')
    gcounts = pd.read_csv(join(proj.expression_dir, 'gene_counts.csv'), index_col=0)

    # add gene names
    if 'gene' not in list(gcounts):
        id2gene = load_id2name(proj)

        isof_index = gcounts.index.tolist()
        gene_index = []
        for i in isof_index:
            if i in id2gene:
                gene_index.append(id2gene[i])
            else:
                print('Warning! gene id ' + str(i) + ' not found')
                gene_index.append('NA')
        se = pd.Series(gene_index)
        gcounts.insert(0, 'gene', se.values)

    # rewrite annotated isoforms
    full_link = join(proj.expression_dir, 'gene_counts.csv')
    gcounts.to_csv(full_link)

    # select key genes
    gcounts_key = gcounts.loc[gcounts['gene'].isin(key_gene_names)]

    # save html for key-genes
    gradient_cols = [sam.name for sam in proj.samples]
    title = '<h3>Gene counts</h3>Showing only key genes. \
    The full results can be downloaded from here: \
    <a href="' + full_link + '">gene_counts.csv</a>'

    html_path = join(proj.expression_dir, 'html', 'gene_counts.html')
    table_to_html(gcounts_key, gradient_cols, title, html_path)
    proj.counts_names.append(html_path)


def gene_tpm_html(proj, key_gene_names):
    print('gene_tpm_html')
    gene_tpm = pd.read_csv(join(proj.expression_dir, 'gene_tpm.csv'), index_col=0)

    # add gene names
    if 'gene' not in list(gene_tpm):
        id2gene = load_id2name(proj)

        isof_index = gene_tpm.index.tolist()
        gene_index = []
        for i in isof_index:
            if i in id2gene:
                gene_index.append(id2gene[i])
            else:
                print('Warning! gene id ' + str(i) + ' not found')
                gene_index.append('NA')
        se = pd.Series(gene_index)
        gene_tpm.insert(0, 'gene', se.values)

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


def diff_exp_genes_html(out_dir, proj, key_gene_names):
    de_gene = pd.read_csv(join(out_dir, 'RNA_DE.csv'), index_col=0)
    id2gene = load_id2name(proj)

    isof_index = de_gene.index.tolist()
    gene_index = []
    for i in isof_index:
        if i in id2gene:
            gene_index.append(id2gene[i])
        else:
            print('Warning! gene id ' + str(i) + ' not found')
            gene_index.append('NA')
    se = pd.Series(gene_index)
    de_gene.insert(0, 'gene', se.values)

    # rewrite annotated isoforms
    full_link = join(out_dir, 'de_gene_all.csv')
    de_gene.to_csv(full_link)

    # select key genes
    de_gene_key = de_gene.loc[de_gene['gene'].isin(key_gene_names)]
    de_gene_key = de_gene_key.dropna()

    # save csv for key-genes
    key_link = join(out_dir, 'de_gene_key.csv')
    de_gene_key.to_csv(key_link)

    # save html for key-genes
    gradient_cols = ['p','lfc', 'baseMean']
    title = '<h3>DE genes</h3>Showing only key genes. \
    The full results can be downloaded from here: \
    <a href="'+ full_link + '">de_genes.csv</a>'

    html_path = join(proj.expression_dir, 'html', 'de_genes.html')
    table_to_html(de_gene_key, gradient_cols, title, html_path)
    proj.counts_names.append(html_path)


def run_QC(proj, key_gene_names):
    out_dir = safe_mkdir(proj.work_dir + '/RNAanalysis')
    pathRScript = which('Rscript')
    pathQC_R = dirname(__file__) + '/QC.R'
    comp_file_for_reuse = proj.date_dir + '/combined.counts'

    if not can_reuse(out_dir, comp_file_for_reuse):
        cmd = ' '.join([pathRScript, pathQC_R, proj.final_dir, out_dir])
        print('Running:')
        print(cmd)
        # os.system(cmd)

    fnames = ['rawCounts.csv', 'normalizedCounts.csv', 'rlog.csv', 'vst.csv', \
              'gene.est.csv', 'gene.final.csv', 'gene.fitted.csv', 'corMatrix.csv', 'pca.csv']
    fnames = [join(out_dir, f) for f in fnames]

    fnames.append(join(proj.date_dir, 'combined.counts'))

    return fnames

def run_DE(proj, key_gene_names):

    pathRScript = which('Rscript')
    pathQC_R = dirname(__file__) + '/DE.R'

    contrast = pd.read_csv(join(proj.date_dir, 'contrasts.csv'), header=None)
    contrast = [tuple(x) for x in contrast.to_records(index=False)]
    for c in contrast:
        print(c)
        de_out = safe_mkdir(join(proj.work_dir, 'RNAanalysis', c[0]+'_'+c[1]))
        cmd = ' '.join([pathRScript, pathQC_R, proj.final_dir, de_out, c[0], c[1]])
        print('Running:')
        print(cmd)
        os.system(cmd)

        # diff_exp_genes_html(de_out, proj, key_gene_names)

    # prepare file-links

    fnames = ['rawCounts.csv', 'normalizedCounts.csv', 'rlog.csv', 'vst.csv', \
              'gene.est.csv', 'gene.final.csv', 'gene.fitted.csv', 'corMatrix.csv', 'pca.csv', \
              'de_gene_key.csv']
    fnames = [join(out_dir, f) for f in fnames]

    fnames.append(join(proj.date_dir, 'combined.counts'))

    return fnames


def run_FA(proj):

    fa_in = join(proj.work_dir, 'RNAanalysis', 'de_gene_key.csv')
    fa_out = join(proj.work_dir, 'RNAanalysis')

    pathRScript = which('Rscript')
    pathFA_R1 = dirname(__file__) + '/FA1.R'

    contr = pd.read_csv(join(proj.date_dir, 'contrasts.csv'), header=None)
    contr = [tuple(x) for x in contr.to_records(index=False)]
    for c in contr:
        fa_out = safe_mkdir(join(proj.work_dir, 'RNAanalysis', c[0] + '_' + c[1]))
        cmd = ' '.join([pathRScript, pathFA_R1, fa_in, fa_out])
        print('Running:')
        print(cmd)
        # os.system(cmd)

        pathRScript = which('Rscript')
        g_obj = join(fa_out, 'genes_obj.csv')
        pw = join(fa_out, 'RNA_PW.csv')
        pathFA_R2 = dirname(__file__) + '/FA2.R'

        cmd = ' '.join([pathRScript, pathFA_R2, g_obj, pw, fa_out])
        print('Running:')
        print(cmd)
        # os.system(cmd)

    return join(fa_out, 'RNA_PW.csv')


def prepare_project_summary(proj):

    path = join(proj.date_dir, 'project-summary.yaml')

    if not isfile(path):
        copyfile(join(proj.log_dir, 'project-summary.yaml'), path)

    gr = pd.read_csv(join(proj.date_dir, 'sample_group.csv'), index_col=0)

    stream = open(path, 'r')
    data = yaml.load(stream)

    for s in data['samples']:
            s_name = s['summary']['metrics']['Name']
            s['metadata']['group'] = gr.at[s_name, 'condition']

    stream2 = open(path, 'w')
    yaml.dump(data, stream2)

def run_analysis(proj, key_gene_names):
    info('*' * 70)
    info('running RNA analysis')

    # expression levels
    # calculate_expression_levels(proj)
    # safe_mkdir(join(proj.expression_dir, 'html'))

    # isoform_level_html(proj, key_gene_names)
    # exon_level_html(proj, key_gene_names)
    # gene_counts_html(proj, key_gene_names)
    # gene_tpm_html(proj, key_gene_names)

    rna_files_list = []
    prepare_project_summary(proj)

    bcbioRNASeq_out_files = run_QC(proj, key_gene_names)
    rna_files_list.extend(bcbioRNASeq_out_files)

    bcbioRNASeq_out_files = run_DE(proj, key_gene_names)
    rna_files_list.extend(bcbioRNASeq_out_files)

    # FA analysis

    fa_file = run_FA(proj)
    rna_files_list.extend([fa_file])

    for p in rna_files_list:
        proj.postproc_mqc_files.append(p)

