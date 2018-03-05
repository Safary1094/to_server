import csv
import traceback
import subprocess

from collections import defaultdict
from os import listdir
from os.path import join, dirname, abspath

from ngs_utils.call_process import run
from ngs_utils.file_utils import verify_file, which, can_reuse
from ngs_utils.logger import info, critical, err


def create_csv_file(bcbio_proj):
    csv_fpath = join(bcbio_proj.work_dir, 'pca_data.csv')
    csv_files_in_config_dir = [
        join(bcbio_proj.config_dir, fname)
        for fname in listdir(bcbio_proj.config_dir)
        if fname.endswith('.csv')]

    csv_data = defaultdict(lambda : defaultdict(dict))
    if csv_files_in_config_dir:
        with open(csv_files_in_config_dir[0]) as f:
            csv_reader = csv.DictReader(f)
            name_col = csv_reader.fieldnames[0]
            for r in csv_reader:
                sample_name = r[name_col]
                csv_data[sample_name] = {'description': r['description'] if 'description' in r else None,
                                         'condition': r['condition'] if 'condition' in r else None}
    else:
        info('No CSV file found in config dir ' + bcbio_proj.config_dir)
    with open(csv_fpath, 'w') as f:
        f.write(','.join(['samplename', 'description', 'condition']) + '\n')
        for sample in bcbio_proj.samples:
            name, description, condition = sample.name, sample.name, sample.name
            old_name = sample.raw_name
            if csv_data[old_name] and csv_data[old_name]['description']:
                description = csv_data[old_name]['description']
            if csv_data[old_name] and csv_data[old_name]['condition']:
                condition = csv_data[old_name]['condition']
            f.write(','.join([name, description, condition]) + '\n')
    return csv_fpath


def create_rnaseq_pca_plot(bcbio_proj, gene_counts_fpath):
    info('Making RNASeq PCA plot')
    output_fpath = join(bcbio_proj.expression_dir, 'pca_data.txt')
    if can_reuse(output_fpath, gene_counts_fpath):
        info('PCA plot ' + output_fpath + ' exists, reusing...')
        return output_fpath

    rscript = which('Rscript')
    if not rscript: critical('Rscript not found in PATH')

    pca_r = verify_file(join(dirname(abspath(__file__)), 'pca.R'), is_critical=True)
    gene_counts_fpath = verify_file(gene_counts_fpath, is_critical=True)
    csv_fpath = create_csv_file(bcbio_proj)

    cmdl = '{rscript} {pca_r} {csv_fpath} {output_fpath} {gene_counts_fpath}'.format(**locals())
    try:
        run(cmdl, output_fpath=output_fpath, stdout_to_outputfile=False)
    except subprocess.CalledProcessError:
        err('Error running pca.R:\n' + traceback.format_exc())
        return None
    else:
        return output_fpath
