import os
import az
import pandas as pd
from os.path import join, dirname, isfile, basename
from shutil import copyfile
import yaml

from ngs_utils.call_process import run
from ngs_utils.file_utils import which, can_reuse, safe_mkdir, verify_file
from ngs_utils.logger import debug, info, critical, err, warn
from ngs_utils.parallel import parallel_view
from ngs_reporting.bcbio.multiqc import _make_link
from ngs_reporting.rnaseq.table_to_html import table_to_html
from ngs_reporting.reference_data import get_name_by_gene_id, get_name_by_transcript_id
from ngs_utils.utils import is_us

def load_id2name(proj):
    id2name = pd.read_csv(get_name_by_gene_id(proj.genome_build))
    id2name = id2name.set_index('ensId')
    id2name = id2name.to_dict()
    return id2name['Name']


def load_tx2name(proj):
    tx2gene = pd.read_csv(get_name_by_transcript_id(proj.genome_build))
    tx2gene = tx2gene.set_index('ensTx')
    tx2gene = tx2gene.to_dict()
    return tx2gene['Gene']


def correct_tx2gene(tx2gene_fpath, output_dir):
    # tx2gene can contain transcript indexes in wrong format (ENST00000346219.5 instead of ENST00000346219)
    tx2gene = pd.read_csv(tx2gene_fpath, header=None)
    tx2gene[0] = tx2gene[0].str.replace(r'\.\d+', '')
    corrected_tx2gene_fpath = join(output_dir, basename(tx2gene_fpath))
    tx2gene.to_csv(corrected_tx2gene_fpath, header=None, index=None)
    return corrected_tx2gene_fpath


def get_sf_file(sample):
    local_quant_path = ['salmon', 'sailfish/quant']
    for p in local_quant_path:
        quant_fpath = join(sample.dirpath, p, 'quant.sf')
        if verify_file(quant_fpath):
            return quant_fpath


def run_tximport(proj, tx2gene_fpath):
    # computes gene_counts, gene_tpms, isoform tpms from quant.sf
    pathEXC_R = dirname(__file__) + '/run_tximport.R'
    cmdl = ['Rscript', pathEXC_R]

    local_quant_path = ['sailfish/quant', 'salmon']
    for s in proj.samples:
        for p in local_quant_path:
            quant_path = join(s.dirpath, p, 'quant.sf')
            if isfile(quant_path):
                cmdl.append(quant_path)

    cmdl.append(proj.expression_dir)
    cmdl.append(tx2gene_fpath)
    run(' '.join(cmdl))


def run_feature_counts(sample, output_dir, transcripts_fpath):
    # computes, exon_counts from .bam
    pathEXC_R = dirname(__file__) + '/exon_counts.R'
    cmdl = ['Rscript', pathEXC_R, sample.bam, output_dir, transcripts_fpath]
    run(cmdl)


def calculate_expression_levels(proj, parallel_cfg):
    tx2gene_fpath = correct_tx2gene(join(proj.date_dir, 'tx2gene.csv'), proj.expression_dir)
    work_dir = safe_mkdir(join(proj.work_dir, 'expression'))
    tximport_output_fpaths = [join(proj.expression_dir, fname) for fname in ['gene_counts.csv', 'gene_tpm.csv', 'isoform_tpm.csv']]
    if not all(can_reuse(fp, [get_sf_file(s) for s in proj.samples]) for fp in tximport_output_fpaths):
            info()
            info('Running txImport...')
            with parallel_view(1, parallel_cfg, work_dir) as view:
                view.run(run_tximport, [[proj, tx2gene_fpath]])

    else:
        info('Reusing txImport output...')

    if not can_reuse(join(proj.expression_dir, 'exon_counts.csv'), [s.bam for s in proj.samples]):
            info()
            info('Running featureCounts...')
            transcripts_fpath = az.get_refdata(proj.genome_build)['exon_annotation']
            with parallel_view(len(proj.samples), parallel_cfg, work_dir) as view:
                view.run(run_feature_counts, [[sample, work_dir, transcripts_fpath] for sample in proj.samples])
    else:
        info('Reusing featureCount output...')


def merge_csv(input_dir, samples, suffix, merged_fpath):
    csv_fpaths = [join(input_dir, s.name + suffix) for s in samples if verify_file(join(input_dir, s.name + suffix), silent=True)]
    if not can_reuse(merged_fpath, csv_fpaths):
        dfs = [pd.read_csv(csv_fpath, index_col=0) for csv_fpath in csv_fpaths]
        if dfs:
            merged_df = pd.concat(dfs, axis=1)
            merged_df.to_csv(merged_fpath)


def build_report_html(proj, fname, plot_title, key_gene_names):
    csv_fpath = join(proj.expression_dir, fname + '.csv')
    if not verify_file(csv_fpath):
        warn('Could not calculate ' + fname.replace('_', ' '))
        return

    counts = pd.read_csv(csv_fpath, index_col=0)
    # add gene names
    if 'gene' not in list(counts):
        id2gene = load_tx2name(proj) if 'isoform' in fname else load_id2name(proj)

        feature_index = counts.index.tolist()
        gene_index = []
        for i in feature_index:
            if i in id2gene:
                gene_index.append(id2gene[i])
            else:
                # debug('Warning! transcript index ' + str(i) + ' not found')
                gene_index.append('NA')
        se = pd.Series(gene_index)
        counts.insert(0, 'gene', se.values)

    # rewrite annotated counts
    counts.to_csv(csv_fpath)
    # select key genes
    key_counts = counts.loc[counts['gene'].isin(key_gene_names)]
    if 'exon' in fname:
        key_counts = key_counts.reset_index()

    # save html for key genes
    html_dirpath = join(proj.expression_dir, 'html')
    gradient_cols = [s.name for s in proj.samples]
    title = '<h3>' + plot_title + '</h3>Showing only key genes. \
    The full results can be downloaded from here: ' + \
            _make_link(csv_fpath, html_dirpath, text=None, blank=False)

    html_fpath = join(proj.expression_dir, 'html', fname + '.html')
    table_to_html(key_counts, gradient_cols, title, html_fpath)
    return html_fpath


def diff_exp_genes_html(out_dir, proj, key_gene_names):
    csv_fpath = join(out_dir, 'RNA_DE.csv')

    de_gene = pd.read_csv(csv_fpath, index_col=0)
    id2gene = load_id2name(proj)
    isof_index = de_gene.index.tolist()
    gene_index = []
    for i in isof_index:
        if i in id2gene:
            gene_index.append(id2gene[i])
        else:
            # debug('Warning! gene id ' + str(i) + ' not found')
            gene_index.append('NA')
    se = pd.Series(gene_index)
    de_gene.insert(0, 'gene', se.values)

    # rewrite annotated isoforms
    html_dirpath = join(proj.expression_dir, 'html')
    out_csv_fpath = join(out_dir, 'de_gene_all.csv')
    de_gene.to_csv(out_csv_fpath)

    # select key genes
    de_gene_key = de_gene.loc[de_gene['gene'].isin(key_gene_names)]
    de_gene_key = de_gene_key.dropna()

    # save csv for key-genes
    key_link = join(out_dir, 'de_gene_key.csv')
    de_gene_key.to_csv(key_link)

    # save html for key-genes
    gradient_cols = ['p','lfc', 'baseMean']
    title = '<h3>DE genes</h3>Showing only key genes. \
    The full results can be downloaded from here: ' + \
            _make_link(out_csv_fpath, html_dirpath, text=None, blank=False)

    html_path = join(proj.expression_dir, 'html', 'de_genes.html')
    table_to_html(de_gene_key, gradient_cols, title, html_path)
    proj.counts_names.append(html_path)


def run_bcbioRNASeq(proj, key_gene_names):
    out_dir = safe_mkdir(proj.work_dir + '/RNAanalysis')
    pathRScript = which('Rscript')
    path_R = dirname(__file__) + '/bcbioRNAseq.R'
    comp_file_for_reuse = proj.date_dir + '/combined.counts'
    contrasts_path = join(proj.date_dir, 'contrasts.csv')

    if not can_reuse(join(out_dir, 'corMatrix.csv'), comp_file_for_reuse):
        cmdl = ' '.join([pathRScript, path_R, proj.final_dir, out_dir])
        if is_us():
            cmdl = ' '.join([cmdl, 'is_us'])
        else:
            cmdl = ' '.join([cmdl, 'no_us'])
        if isfile(contrasts_path):
            cmdl = ' '.join([cmdl, contrasts_path])
        else:
            cmdl = ' '.join([cmdl, 'no_contrast'])

        print(cmdl)
        run(cmdl)

    fnames_key, fnames_all = [], []

    if isfile(contrasts_path):

        contrasts = pd.read_csv(contrasts_path, header=None)
        for ind, c in contrasts.iterrows():
            contrast_dir_path = join(proj.work_dir, 'RNAanalysis', 'group_' + c[0] + '_vs_' + c[1])
            diff_exp_genes_html(contrast_dir_path, proj, key_gene_names)

            # prepare file-links
            fnames_key.append(join(contrast_dir_path, 'de_gene_key.csv'))
            fnames_all.append(join(contrast_dir_path, 'de_gene_all.csv'))

    # prepare file-links

    qc = ['rawCounts.csv', 'normalizedCounts.csv', 'rlog.csv', 'vst.csv',
              'gene.est.csv', 'gene.final.csv', 'gene.fitted.csv', 'corMatrix.csv', 'pca.csv']
    qc_names = [join(out_dir, f) for f in qc]
    qc.append(join(proj.date_dir, 'combined.counts'))


    return qc_names, fnames_key, fnames_all


def run_FA(de_files_list):
    pathRScript = which('Rscript')
    pathFA_R1 = join(dirname(__file__), 'FA1.R')
    pathFA_R2 = join(dirname(__file__), 'FA2.R')

    fnames = []

    # loop over all contrasts and determine all pathways
    for de_csv in de_files_list:
        dir_path = dirname(de_csv)
        cmd = ' '.join([pathRScript, pathFA_R1, dir_path, de_csv])
        print('Running:')
        print(cmd)
        # os.system(cmd)

    # summarize pathways
    pathways = set()
    for de_csv in de_files_list:
        pathway_table_path = join(dirname(de_csv), 'pathway_table.csv')
        d = pd.read_csv(pathway_table_path, index_col=[0])
        pathways = pathways | set(d.index.tolist())

    # complete each pathway table with missing pathways
    for de_csv in de_files_list:
        pathway_table_path = join(dirname(de_csv), 'pathway_table.csv', )
        d = pd.read_csv(pathway_table_path, index_col=[0])

        for pw in pathways:
            if pw not in d.index.tolist():
                d.loc[pw] = 'NA'

        d.sort_index(inplace=True)

        d.to_csv(pathway_table_path)

    for de_csv in de_files_list:
        dir_path = dirname(de_csv)

        cmd = ' '.join([pathRScript, pathFA_R2, dir_path])
        print('Running:')
        print(cmd)
        # os.system(cmd)

        # write pathways for multiqc
        fnames.append(join(dir_path, 'pathway_table.csv'))

    return fnames


def prepare_project_summary(proj):

    path = join(proj.date_dir, 'project-summary.yaml')

    if not verify_file(path):
        copyfile(join(proj.log_dir, 'project-summary.yaml'), path)

    #gr = pd.read_table(join(proj.date_dir, 'sample_group'), index_col=0)

    stream = open(path, 'r')
    data = yaml.load(stream)

    # stupid group assignment
    gr_id = 1
    for s in data['samples']:
        if 'group' not in s['metadata']:
            gr_id *= -1
            #s_name = s['summary']['metrics']['Name']
            s['metadata']['group'] = 'g' + str(gr_id)
            #s['metadata']['group'] = gr.at[s_name, 'condition']

    stream2 = open(path, 'w')
    yaml.dump(data, stream2)


def run_analysis(parallel_cfg, proj, key_gene_names):
    info('*' * 70)
    info('Running RNA analysis...')

    # expression levels
    calculate_expression_levels(proj, parallel_cfg)
    merge_csv(join(proj.work_dir, 'expression'), proj.samples, '_exon_counts.csv', join(proj.expression_dir, 'exon_counts.csv'))
    html_dir = safe_mkdir(join(proj.expression_dir, 'html'))

    info()
    info('Building RNA reports...')
    report_titles = [('isoform_tpm', 'Isoform TPM'), ('gene_counts', 'Gene counts'), ('gene_tpm', 'Gene TPM'),
                     ('exon_counts', 'Exon counts')]
    report_fpaths = [join(html_dir, fname + '.html') for fname, title in report_titles]
    if not all(can_reuse(fp, [get_sf_file(s) for s in proj.samples]) for fp in report_fpaths):
        for fname, plot_title in report_titles:
            build_report_html(proj, fname, plot_title, key_gene_names)
    proj.counts_names = [fp for fp in report_fpaths if verify_file(fp)]
    info('Done')

    # bcbioRNAseq
    prepare_project_summary(proj)

    rna_files_list = []
    qc_files, de_files_key, de_files_all = run_bcbioRNASeq(proj, key_gene_names)
    rna_files_list.extend(qc_files)
    rna_files_list.extend(de_files_key)

    # FA
    if isfile(join(proj.date_dir, 'contrasts.csv')):
        fa_files_list = run_FA(de_files_all)
        rna_files_list.extend(fa_files_list)

    # Prepare files for multiQC

    for p in rna_files_list:
        proj.postproc_mqc_files.append(p)

