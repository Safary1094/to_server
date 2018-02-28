import os
import shutil
import subprocess
from inspect import getsourcefile
from os.path import join, dirname, abspath, pardir, isfile, isdir, islink, basename, relpath, getmtime, exists
from time import localtime, strftime
import yaml

from ngs_reporting import config
from ngs_reporting.rnaseq.pca import create_rnaseq_pca_plot
from ngs_utils import logger, sambamba
from ngs_utils.bcbio import BcbioProject
from ngs_utils.bed_utils import get_total_bed_size
from ngs_utils.call_process import run
from ngs_utils.file_utils import verify_file, safe_mkdir, can_reuse, file_transaction
from ngs_utils.logger import info, warn, err, critical, timestamp, debug
from ngs_utils.utils import mean, is_us


def make_report_metadata(bcbio_proj, base_dirpath,
                         oncoprints_url=None, call_vis_html_fpath=None,
                         combined_ngs_rep_html_fpath=None, analysis_dir=None):
    conf = dict()
    conf['az'] = dict()
    additional_files = []

    conf['az']['is_rnaseq'] = bcbio_proj.is_rnaseq

    # Paired samples relations
    normal_samples = [s for s in bcbio_proj.samples if s.phenotype == 'normal']
    if normal_samples:
        sample_match_on_hover_js = '<script type="text/javascript">\n'
        for s in bcbio_proj.samples:
            if s.phenotype != 'normal' and s.normal_match:
                sample_match_on_hover_js += ('' +
                    '\tdocument.getElementById("' + s.name + '_match").onmouseover = function() { document.getElementById("' + s.normal_match.name + '").style.backgroundColor = "#EEE"; };\n' +
                    '\tdocument.getElementById("' + s.name + '_match").onmouseleave = function() { document.getElementById("' + s.normal_match.name + '").style.backgroundColor = "white"; };\n'
                )
        sample_match_on_hover_js += '</script>\n'
        conf['az']['sample_match_on_hover_js'] = sample_match_on_hover_js

    # General links
    conf['title'] = bcbio_proj.project_name
    conf['az']['run_section'] = get_run_info(bcbio_proj, base_dirpath, analysis_dir=analysis_dir)

    if bcbio_proj.is_rnaseq:
        conf['az']['expression_links'] = _rna_general_links(bcbio_proj, base_dirpath)
    else:
        mutations_links = _dna_general_links(bcbio_proj, base_dirpath)
        conf['az']['mutations_links'] = mutations_links
        if oncoprints_url:
            mutations_links.append(('<a href="{oncoprints_url}" target="_blank">oncoprints</a> ' +
                                    '(loading may take 5-10 seconds)').format(**locals()))
        if call_vis_html_fpath:
            call_vis_link = _make_link(call_vis_html_fpath, base_dirpath, 'call visualisation', blank=True)
            mutations_links.append(call_vis_link)
        if combined_ngs_rep_html_fpath:
            link = _make_link(combined_ngs_rep_html_fpath, base_dirpath, 'known mutations')
            mutations_links.append(link)

    # Sample-level details
    gender_by_sample = dict()
    for s in bcbio_proj.samples:
        gender_fpath = join(s.dirpath, BcbioProject.ngs_report_name, 'gender.txt')
        if isfile(gender_fpath):
            gender = open(gender_fpath).read()
            gender_by_sample[s.name] = gender
    if gender_by_sample:
        conf['az']['gender_by_sample'] = gender_by_sample
    
    ngs_report_by_sample = {s.name:
        relpath(s.find_ngs_report(silent=True), base_dirpath)
        for s in bcbio_proj.samples if s.find_ngs_report(silent=True)}
    if ngs_report_by_sample:
        conf['az']['ngs_report_by_sample'] = ngs_report_by_sample
    
    # Coverage thresholds
    avg_depths = [s.get_avg_depth() for s in bcbio_proj.samples]
    cov_threshs = config.coverage_thresholds
    cov_threshs_hidden = []
    for i, thres in enumerate(cov_threshs):
        if thres > 30:
            cov_threshs_hidden.append(thres)
            if any(cov_threshs[i - 1] <= cov < cov_threshs[i + 1] for cov in avg_depths
                   if cov and i - 1 >= 0 and i + 1 < len(cov_threshs)):
                cov_threshs_hidden.remove(thres)
    conf['qualimap_config'] = {
        'general_stats_coverage': [str(t) for t in cov_threshs],
        'general_stats_coverage_hidden': [str(t) for t in cov_threshs_hidden]}

    # Preseq
    preseq_files, preseq_config = _add_preseq_data(bcbio_proj)
    if preseq_files:
        additional_files.extend(preseq_files)
    if preseq_config:
        conf['preseq'] = preseq_config

    # PCA plot
    if bcbio_proj.is_rnaseq:
        gene_counts_fpath = join(bcbio_proj.expression_dir, 'counts.tsv')
        if not isfile(gene_counts_fpath):
            gene_counts_fpath = join(bcbio_proj.expression_dir, 'combined.counts')
        if gene_counts_fpath:
            pca_plot_fpath = create_rnaseq_pca_plot(bcbio_proj, gene_counts_fpath)
            debug()
            if pca_plot_fpath and verify_file(pca_plot_fpath):
                additional_files.append(pca_plot_fpath)

    # QC DE FA
    for l in bcbio_proj.postproc_mqc_files:
        additional_files.append(l)
        print('QC DE FA ' + l)

    return conf, additional_files


def _rna_general_links(bcbio_proj, base_dirpath):
    expression_dir = join(bcbio_proj.date_dir, BcbioProject.expression_dir)
    links = []
    for fname, text in zip(BcbioProject.counts_names, ['Gene counts', 'Exon counts', 'Gene TPM', 'Isoform TPM']):
        html_path = join(expression_dir, 'html', fname.replace('combined.', '').replace('.tsv', '') + '.html')
        if isfile(html_path):
            links.append(_make_link(html_path, base_dirpath, text))
    return links


def _dna_general_links(bcbio_proj, base_dirpath):
    links = []

    for (caller, is_germline), samples in bcbio_proj.samples_by_caller.items():
        for fpath in bcbio_proj.find_mutation_files(caller=caller, is_germline=is_germline):
            link = _make_link(fpath, base_dirpath) + (' (germline)' if is_germline else '')
            if link not in links:
                links.append(link)

    seq2c_file = bcbio_proj.find_seq2c_file()
    if seq2c_file:
        links.append(_make_link(seq2c_file, base_dirpath))

    return links


def _make_url(fpath, base_dirpath):
    return relpath(fpath, base_dirpath) if verify_file(fpath) else None


def _make_link(fpath, base_dirpath, text=None, blank=False):
    url = _make_url(fpath, base_dirpath)
    if url:
        return '<a href="' + url + '"' + (' target="_blank"' if blank else '') + '>' + (text or basename(fpath)) + '</a>'
    else:
        return '<span>' + (text or basename(fpath)) + '</span>'


def get_run_info(bcbio_proj, base_dirpath, analysis_dir=None):
    info('Getting run and codebase information...')
    run_info_dict = dict()
    cur_fpath = abspath(getsourcefile(lambda: 0))
    reporting_suite_dirpath = dirname(dirname(dirname(cur_fpath)))

    run_info_dict["run_date"] = strftime('%d %b %Y, %H:%M (GMT%z)', localtime())

    version = None
    try:
        from ngs_reporting.version import __version__
    except ImportError:
        err('Cannot import __version__ from ngs_reporting.version and get version')
    else:
        version = __version__

    last_modified_datestamp = ''
    try:
        py_fpaths = set()
        for rootpath, dirnames, fnames in os.walk(reporting_suite_dirpath):
            for fn in fnames:
                if fn.endswith('.py'):
                    fpath = abspath(join(rootpath, fn))
                    if isfile(fpath):
                        py_fpaths.add(fpath)
        last_modification_time = max(getmtime(fpath) for fpath in py_fpaths)
    except Exception:
        warn('Cannot get codebase files datestamps')
    else:
        last_modified_datestamp = strftime('%d %b %Y, %H:%M (GMT%z)', localtime(last_modification_time))

    if last_modified_datestamp or version:
        version_text = 'Reporting Suite '
        if version:
            version_text += 'v.' + version
        if version and last_modified_datestamp:
            version_text += ', '
        if last_modified_datestamp:
            version_text += 'last modified ' + last_modified_datestamp
        run_info_dict['suite_version'] = version_text

    # if bcbio_structure.is_rnaseq:
    program_versions_fpath = bcbio_proj.find_in_log('programs.txt')
    program_versions = dict(l.strip().split(',') for l in open(program_versions_fpath).readlines())
    program_versions['reporting-suite'] = version
    try:
        with open(program_versions_fpath, 'w') as f:
            for p, v in sorted(program_versions.items(), key=lambda kv: kv[0]):
                f.write(p + ',' + v + '\n')
    except OSError as e:
        err(e)
    programs_url = relpath(program_versions_fpath, base_dirpath) \
        if verify_file(program_versions_fpath) else None

    data_versions_fpath = bcbio_proj.find_in_log('data_versions.csv')
    datas_url = relpath(data_versions_fpath, base_dirpath) if verify_file(data_versions_fpath) else None

    run_info_dict['program_versions'] = '<a href="{programs_url}">program versions</a>'.format(**locals())
    run_info_dict['data_versions'] = '<a href="{datas_url}">data versions</a>'.format(**locals())
    run_info_dict['analysis_dir'] = analysis_dir or bcbio_proj.final_dir
    return run_info_dict


def find_multiqc_file_list(proj):
    multiqc_bcbio_dir = join(proj.date_dir, 'multiqc')
    new_multiqc_bcbio_dir = join(proj.date_dir, 'log', 'multiqc_bcbio')
    
    postproc_input_list_file = join(proj.date_dir, 'log', 'multiqc_list_files.txt')
    if isfile(postproc_input_list_file):
        return postproc_input_list_file

    bcbio_list_files = abspath(join(new_multiqc_bcbio_dir, 'list_files_final.txt'))
    if isfile(bcbio_list_files):
        return bcbio_list_files

    bcbio_list_files = abspath(join(multiqc_bcbio_dir, 'list_files_final.txt'))
    if isfile(bcbio_list_files):
        return bcbio_list_files


def make_multiqc_report(proj, az_multiqc_config, additional_files, multiqc_fpath,
                        work_dir=None, move_files=True, extra_params=''):
    work_dir = work_dir or proj.work_dir
    multiqc_bcbio_dir = join(proj.date_dir, 'multiqc')
    if move_files:
        new_multiqc_bcbio_dir = join(proj.date_dir, 'log', 'multiqc_bcbio')
        safe_mkdir(dirname(new_multiqc_bcbio_dir))
        if isdir(multiqc_bcbio_dir):
            if isdir(new_multiqc_bcbio_dir):
                try:
                    debug('Remove ' + new_multiqc_bcbio_dir)
                    shutil.rmtree(new_multiqc_bcbio_dir)
                except OSError:
                    debug('Cannot remove ' + new_multiqc_bcbio_dir + ', backing up')
                    os.rename(new_multiqc_bcbio_dir, new_multiqc_bcbio_dir + '.' + timestamp().replace(':', '_').replace(' ', '_'))
            try:
                debug('Rename ' + multiqc_bcbio_dir + ' -> ' + new_multiqc_bcbio_dir)
                os.rename(multiqc_bcbio_dir, new_multiqc_bcbio_dir)
            except OSError:
                err('Cannot rename ' + multiqc_bcbio_dir + ' to ' + new_multiqc_bcbio_dir)
                raise
    else:
        new_multiqc_bcbio_dir = multiqc_bcbio_dir

    multiqc_postproc_dir = safe_mkdir(dirname(multiqc_fpath))

    cmdl = 'multiqc -f' + (' -v' if logger.is_debug else '') + ' -o ' + multiqc_postproc_dir

    if az_multiqc_config:
        az_multiqc_yaml = join(work_dir, 'az_multiqc_config.yaml')
        with file_transaction(None, az_multiqc_yaml) as tx:
            with open(tx, 'w') as f:
                yaml.dump(az_multiqc_config, f, default_flow_style=False)
        cmdl += ' -c ' + az_multiqc_yaml
    bcbio_multiqc_yaml = join(new_multiqc_bcbio_dir, 'multiqc_config.yaml')
    if isfile(bcbio_multiqc_yaml):
        cmdl += ' -c ' + bcbio_multiqc_yaml

    if is_us():
        clarity_config_fpath = '/users/klpf990/.genologicsrc'
        cmdl += ' --clarity_config ' + clarity_config_fpath
        fastq_dir = dirname(proj.samples[0].sample_info['files'][0])
        samplesheet_fpath = join(fastq_dir.split('/Unalign')[0], 'SampleSheet.csv')
        if verify_file(samplesheet_fpath):
            cmdl += ' --samplesheet ' + samplesheet_fpath
        csv_files_in_config_dir = [
            join(proj.config_dir, fname)
            for fname in os.listdir(proj.config_dir)
            if fname.endswith('.csv')]
        if csv_files_in_config_dir:
            cmdl += ' --bcbio_csv ' + abspath(csv_files_in_config_dir[0])

    bcbio_list_files_final = find_multiqc_file_list(proj)
    postproc_list_file = join(work_dir or join(proj.date_dir, 'log'), 'multiqc_list_files.txt')
    qc_files_not_found = []
    if not can_reuse(postproc_list_file, bcbio_list_files_final):
        with file_transaction(None, postproc_list_file) as tx:
            try:
                with open(bcbio_list_files_final) as inp, open(tx, 'w') as out:
                    for l in inp:
                        fpath = join(proj.final_dir, l.strip()).replace(multiqc_bcbio_dir, new_multiqc_bcbio_dir)
                        if not verify_file(fpath):
                            qc_files_not_found.append(fpath)
                            continue
                        out.write(fpath + '\n')

                    if qc_files_not_found:
                        warn('-')
                        warn('Some QC files from list ' + bcbio_list_files_final + ' were not found:' +
                            ''.join('\n  ' + fpath for fpath in qc_files_not_found))
                    for fp in additional_files:
                        if fp:
                            out.write(fp + '\n')
            except OSError as e:
                err(e)



    if verify_file(postproc_list_file , silent=True):
        cmdl += ' -l ' + postproc_list_file
    else:
        if not isfile(postproc_list_file ):
            critical('Critical: MultiQC files list was not found in ' + postproc_list_file )

    if extra_params:
        cmdl += ' ' + extra_params

    if isfile(multiqc_fpath):
        try:
            os.remove(multiqc_fpath)
        except OSError as e:
            err(e)
    try:
        run(cmdl, env_vars={'https_proxy': None})
    except subprocess.CalledProcessError:
        pass
    if not verify_file(multiqc_fpath, silent=True):
        critical('Error: MultiQC has failed. Please, make sure bcbio_postproc is properly loaded into PATH. '
                 'Try starting from a clean session.')
    return multiqc_fpath


def _add_preseq_data(proj):
    preseq_files = []
    samples = []
    for s in proj.samples:
        fp = join(join(s.dirpath, 'preseq'), s.name + '.txt')
        if verify_file(fp, silent=True):
            preseq_files.append(fp)
            samples.append(s)
    if not samples:
        return [], None

    work_dir = safe_mkdir(proj.work_dir, 'preseq')
    actual_counts_file = join(work_dir, 'preseq_real_counts.txt')
    
    preseq_config = dict()

    if proj.coverage_bed:
        target_size = get_total_bed_size(proj.coverage_bed)
        debug('target_size: ' + str(target_size))

        ontarget_counts = []
        ontarget_unique_counts = []
        ontarget_unique_depths = []
        ontarget_alignments_avg_lengths = []
        for s in samples:
            u_dp = s.get_avg_depth()
            ontarget_unique_depths.append(u_dp)

            cnt = sambamba.number_mapped_reads_on_target(work_dir, bed=proj.coverage_bed, bam=s.bam, dedup=False)
            ontarget_counts.append(cnt)
            if s.is_dedupped():
                u_cnt = sambamba.number_mapped_reads_on_target(work_dir, bed=proj.coverage_bed, bam=s.bam, dedup=True)
                ontarget_unique_counts.append(u_cnt)
            else:
                u_cnt = cnt
                ontarget_unique_counts.append(None)

            # avg depth ~~ unique mapped reads on target (usable reads) * on target aln length / target_size
            read_len = u_dp * target_size // u_cnt
            ontarget_alignments_avg_lengths.append(read_len)

        avg_ontarget_alignments_avg_length = mean(ontarget_alignments_avg_lengths)

        debug('ontarget_counts: ' + str(ontarget_counts))
        debug('ontarget_unique_counts: ' + str(ontarget_counts))
        debug('ontarget_unique_depths: ' + str(ontarget_unique_depths))
        debug('ontarget_alignments_avg_lengths: ' + str(ontarget_alignments_avg_lengths))
        debug('avg_ontarget_alignments_avg_length: ' + str(avg_ontarget_alignments_avg_length))

        with open(actual_counts_file, 'w') as f:
            for i, s in enumerate(samples):
                line = s.name + '\t' + str(ontarget_counts[i])
                if s.is_dedupped():
                    line += '\t' + str(ontarget_unique_counts[i])
                line += '\n'
                f.write(line)
        
        preseq_config['genome_size'] = target_size
        preseq_config['read_length'] = avg_ontarget_alignments_avg_length

    else:  # WGS
        preseq_config['genome_size'] = proj.genome_build
        try:
            read_lengths = [int(s.get_metric('Sequence length').split('-')[-1]) for s in samples]
        except:
            pass
        else:
            preseq_config['read_length'] = mean(read_lengths)
            debug('read_lengths: ' + str(read_lengths))
            debug('mean(read_lengths): ' + str(mean(read_lengths)))

        with open(actual_counts_file, 'w') as f:
            for s in samples:
                count = s.get_reads_count()
                line = s.name + '\t' + str(count)
                if s.is_dedupped():
                    dup_count = s.get_metric('Duplicates')
                    unique_count = s.get_reads_count() - dup_count
                    line += '\t' + str(unique_count)
                line += '\n'
                f.write(line)

    preseq_files = [fp for fp in preseq_files if fp] + [actual_counts_file]
    return preseq_files, preseq_config
