import datetime
import math
import os
import shutil
import sys
import traceback
from collections import defaultdict, OrderedDict
from os import listdir
from os.path import join, dirname, isdir, islink, basename, samefile, isfile, relpath
from traceback import format_exc

from ngs_utils.bcbio import BcbioProject, MAIN_CALLER, CALLER_PRIORITY
from ngs_utils.bed_utils import clean_bed, cut, sort_bed, count_bed_cols, annotate_target
from ngs_utils.call_process import run
from ngs_utils.file_utils import verify_file, add_suffix, safe_mkdir, verify_dir, symlink_plus, file_transaction, \
    can_reuse, iterate_file, which, open_gzipsafe, intermediate_fname
from ngs_utils.logger import debug, info, warn, err, critical, send_email, CriticalError
from ngs_utils.proc_args import set_up_log
from ngs_utils.utils import is_us, is_local, is_uk, is_sweden
import ngs_utils.variant_filtering as vf

import az
from az.ngb import add_bcbio_project_to_ngb
from az.webserver.exposing import convert_gpfs_path_to_url, sync_with_ngs_server
from ngs_reporting import is_small_target
from ngs_reporting.bcbio.bcbio_filtering import finish_filtering_for_bcbio
from ngs_reporting.bcbio.multiqc import make_report_metadata, make_multiqc_report
from ngs_reporting.clinical_report.clinical_parser import NgsReportSample, ClinicalExperimentInfo
from ngs_reporting.clinical_report.clinical_reporting import make_clinical_report, GeneralMutationReporting
from ngs_reporting.cnvkit.cnvkit import merge_cnvkit_files
from ngs_reporting.coverage import extract_features_for_genes
from ngs_reporting.oncoprints import create_oncoprints_link
from ngs_reporting.panel_evaluation import evaluate_panel, threshs_by_cov_interval
from ngs_reporting.reference_data import get_key_genes_bed, get_key_genes_file
from ngs_reporting.rnaseq import gene_expression, RNAanalysis
from ngs_reporting.seq2c.seq2c import seq2c_summarize_bcbio_proj
from ngs_reporting.visualisation import make_circos_and_linear_plots

from variant_filtering.anno import finialize_annotate_file, run_vcfanno
from variant_filtering.config import get_dbsnp_multi_mafs, detect_run_info_in_config_dir, save_filt_cnf
from variant_filtering.filtering import run_filtering
from variant_filtering.vcf import clean_up_vcf


def run_postproc(proj, parallel_cfg, filt_cnf, jira_url=None, re_summarize=False, seq2c_controls=None,
                 run_eval_panel=False, skip_ngb=False, preseq=False):
    safe_mkdir(proj.work_dir)
    
    # If variants present, save filtering config
    if proj.samples_by_caller:
        _save_run_cnf(proj, filt_cnf)

    error_msg = []
    try:
        clean_up_bcbio_dir(proj)

        if proj.sv_regions_bed:
            proj.sv_regions_bed = _prep_seq2c_bed(proj)
        if proj.coverage_bed:
            proj.coverage_bed = _prep_bed(safe_mkdir(join(proj.work_dir, 'coverage_bed')), proj.coverage_bed, proj.genome_build)
        for s in proj.samples:
            s.sv_regions_bed = proj.sv_regions_bed
            s.coverage_bed = proj.coverage_bed
        key_gene_names = proj.get_target_genes(get_key_genes_file)

        oncoprints_url = None
        avg_depths = None
        call_vis_html_fpath = None
        combined_ngs_rep_html_fpath = None

        from ngs_utils.parallel import parallel_view

        genome_cfg = az.get_refdata(proj.genome_build)

        if proj.is_rnaseq:
            info('*' * 70)
            info('Analysing expression')
            # gene_expression.make_heatmaps(proj, key_gene_names)
            RNAanalysis.run_analysis(proj, key_gene_names)

        with parallel_view(len(proj.samples), parallel_cfg, safe_mkdir(proj.work_dir)) as view:
            if preseq:
                info()
                info('Analysing comlexity using Preseq')
                run_preseq(view, proj)

            if filt_cnf or not proj.is_rnaseq:
                info()
                info('*' * 70)
                info('Analysing genomic variation')
                info('*' * 70)
                info('Processing CNV...')
                seq2c_files = _symlink_cnv(proj)
                if seq2c_files:
                    seq2c_summarize_bcbio_proj(view, proj, re_summarize=re_summarize, seq2c_controls=seq2c_controls)

                if filt_cnf:
                    for (caller, is_germline), samples in proj.samples_by_caller.items():
                        info('Processing variants for variant caller ' + str(caller) + (', germline' if is_germline else ''))
                        raw_vcf_files = [_f for _f in [s.find_raw_vcf(silent=True, caller=caller) for s in samples] if _f]
                        proj_mut_files = proj.find_mutation_files(passed=False, caller=caller, is_germline=is_germline)
                        if not proj_mut_files or not all(can_reuse(fp, raw_vcf_files) for fp in proj_mut_files):
                            debug()
                            info('Annotating variants...')
                            # We have to use len(proj.samples) number of subprocesses. So for samples not in this caller,
                            #    we pass the `ignore` flag that equals `s in samples` and make it skip this sample.
                            view.run(annotate_fn, [[proj, s, filt_cnf, genome_cfg, caller, view.cores_per_job, s not in samples]
                                                   for s in proj.samples])
                            debug()
                            info('Filtering variants...')
                            dbsnp_multi_mafs = get_dbsnp_multi_mafs(genome_cfg)
                            view.run(filtering_fn, [[proj, s, filt_cnf, dbsnp_multi_mafs, caller, view.cores_per_job, s not in samples]
                                                    for s in proj.samples])
                            finish_filtering_for_bcbio(parallel_cfg, proj, filt_cnf, caller=caller, is_germline=is_germline)

                            info()
                            info('*' * 70)

            if not skip_ngb and (is_us() or is_uk() or is_sweden()):
                debug()
                info('Exposing to NGB...')
                try:
                    add_bcbio_project_to_ngb(proj.work_dir, proj, view)
                except Exception:
                    err(traceback.format_exc())
                    err('Error: cannot export to NGB')
                info('*' * 70)

            if filt_cnf:
                info('Reporting')
                info('*' * 70)
                combined_ngs_rep_html_fpath, call_vis_html_fpath = \
                    _visualization(view, proj, filt_cnf, key_gene_names, skip_ngb=skip_ngb)
            
        if run_eval_panel:
            debug()
            info('Evaluating capture panels...')
            work_dir = safe_mkdir(join(proj.work_dir, BcbioProject.evaluate_panel_dir))
            output_dir = safe_mkdir(join(proj.date_dir, 'qc', BcbioProject.evaluate_panel_dir))
            for depth in threshs_by_cov_interval.get(proj.coverage_interval or 'regional'):
                evaluate_panel([proj], output_dir, work_dir, min_depth=depth, annotate_with_tricky=False)

        if is_us():
            debug()
            info('Exporting to oncoprints...')
            sample_names = [s.name for s in proj.samples]
            work_dir = safe_mkdir(join(proj.work_dir, BcbioProject.oncoprints_dir))
            oncoprints_url = create_oncoprints_link(
                work_dir, proj.genome_build, sample_names, proj.final_dir, proj.project_name,
                proj.find_mutation_files(), proj.find_cnv_filt_file(),
                bed_fpath=proj.coverage_bed or proj.sv_regions_bed)
        
        debug()
        info('Building MultiQC report...')
        if proj.coverage_bed:
            shutil.copyfile(proj.coverage_bed, join(proj.date_dir, basename(proj.coverage_bed)))
        multiqc_fpath = join(proj.date_dir, BcbioProject.multiqc_report_name)
        az_multiqc_conf, additional_files = make_report_metadata(proj, dirname(multiqc_fpath),
            oncoprints_url=oncoprints_url, call_vis_html_fpath=call_vis_html_fpath,
            combined_ngs_rep_html_fpath=combined_ngs_rep_html_fpath)
        make_multiqc_report(proj, az_multiqc_conf, additional_files, multiqc_fpath)

        multiqc_url = None
        if multiqc_fpath:
            info()
            info('Exposing to NGS webserver...')
            try:
                multiqc_url = sync_with_ngs_server(
                    work_dir=proj.work_dir,
                    jira_url=jira_url,
                    project_name=proj.project_name,
                    sample_names=[s.name for s in proj.samples],
                    bcbio_final_dirpath=proj.final_dir,
                    summary_report_fpath=multiqc_fpath)
            except:
                err('Could not expose reports:')
                warn(traceback.format_exc())

        final_email_notification(multiqc_url, multiqc_fpath, jira_url, proj.final_dir, proj.project_name)
        info()
        info('Report path: ' + multiqc_fpath)
        if multiqc_url:
            info('Report URL: ' + multiqc_url)

    except KeyboardInterrupt:
        warn('Interrupted.')
    except SystemExit:
        warn('Interrupted.')
    except CriticalError as e:
        error_msg = e.args[0]
    except Exception:
        err(traceback.format_exc())
    finally:
        if error_msg:
            warn('')
            warn('-' * 70)
            err('Error: ' + error_msg)
            warn('')
            warn('Post-processing finished with error. Try rerunning in single-threaded debug mode ' +
                 'to get more detailed error information (add -d -t 1 to the command line). '
                 'Send the project location on the HPC to Vlad Saveliev to get support.')
            sys.exit(1)
        else:
            _move_work_to_scratch(proj.work_dir)
            warn('')
            warn('-' * 70)
            info('Done post-processing.')


def run_preseq(view, proj):
    read_counts = [s.get_reads_count() for s in proj.samples]
    max_read_count = max(read_counts)
    info('max_read_count = ' + str(max_read_count))
    
    unrounded__max_extrap_read_count = max_read_count * 3
    total_steps = 300
    unrounded__step = unrounded__max_extrap_read_count // total_steps
    power_of_10 = 10**math.floor(math.log(unrounded__step, 10))
    rounded__step = int(math.floor(unrounded__step // power_of_10) * power_of_10)
    rounded__max_extrap_read_count = int(rounded__step) * total_steps
    info('Total ' + str(total_steps) + ' steps of size ' + str(rounded__step) +
         ', X limit = ' + str(rounded__max_extrap_read_count))

    view.run(preseq_fn, [[proj, s, rounded__step, rounded__max_extrap_read_count]
                         for s in proj.samples])


def preseq_fn(proj, sample, step, max_extrap_read_count):
    if not sample.bam:
        return None

    output_file = join(safe_mkdir(join(sample.dirpath, 'preseq')), sample.name + '.txt')
    if can_reuse(output_file, sample.bam):
        return output_file

    if is_us():
        preseq = 'LD_LIBRARY_PATH="/usr/lib64:$LD_LIBRARY_PATH" /users/klpf990/bin/preseq'
    elif is_local():
        preseq = 'preseq'
    else:
        preseq = which('preseq')
    if not preseq:
        warn('Preseq is not found in PATH')
        return

    # sorted_bam = intermediate_fname(safe_mkdir(join(proj.work_dir, sample.name)), sample.bam, 'sorted')
    # run('sambamba sort ')

    cmdl = '{preseq} lc_extrap -bam -pe {sample.bam} -o {output_file} ' \
           '-s {step} -e {max_extrap_read_count} -l 100000'.format(**locals())
    run(cmdl, output_fpath=output_file, stdout_to_outputfile=False)


def _prep_bed(work_dir, bed_file, genome_build):
    if not bed_file:
        return None

    debug('Remove incorrect lines in BED...')
    bed_file = clean_bed(bed_file, work_dir)
    
    debug('Sorting BED...')
    bed_file = sort_bed(bed_file, work_dir=work_dir, genome=genome_build)
    
    cols = count_bed_cols(bed_file)
    if cols < 4:
        info()
        info('Number columns in SV bed is ' + str(cols) + '. Annotating regions with gene names...')
        bed_file = annotate_target(work_dir, bed_file, genome_build)
    if 8 > cols > 4:
        bed_file = cut(bed_file, 4)
    elif cols > 8:
        bed_file = cut(bed_file, 8)
    debug('Done cleaning BED: ' + bed_file)
    return bed_file


def _prep_seq2c_bed(proj):
    if not proj.sv_regions_bed:
        return None
    if proj.coverage_bed == proj.sv_regions_bed:
        sv_regions_dir = 'coverage_bed'
    else:
        sv_regions_dir = 'sv_regions_bed'
    work_dir = safe_mkdir(join(proj.work_dir, sv_regions_dir))

    sv_regions_bed = _prep_bed(work_dir, proj.sv_regions_bed, proj.genome_build)
    if not sv_regions_bed:
        return None
    
    debug('Seq2C bed: removing regions with no gene annotation')
    def f(l, i):
        if l.split('\t')[3].strip() == '.':
            return None
        else:
            return l
    return iterate_file(work_dir, sv_regions_bed, f, suffix='filt')


def _move_work_to_scratch(work_dir):
    bcbio_work_dirpath = dirname(work_dir).replace('/Analysis/', '/analysis/')

    scratch_root_dirpath = None
    analysis_root_dirpath = None
    if is_us() and '/ngs/oncology/analysis/' in bcbio_work_dirpath:
        scratch_root_dirpath = '/ngs/scratch/'
        analysis_root_dirpath = '/ngs/oncology/analysis/'
    elif is_local() and '/Users/vlad/googledrive/az/analysis/' in bcbio_work_dirpath:
        scratch_root_dirpath = '/Users/vlad/scratch/'
        analysis_root_dirpath = '/Users/vlad/googledrive/az/analysis/'
    if scratch_root_dirpath and isdir(bcbio_work_dirpath) and not islink(bcbio_work_dirpath):
        info()
        work_scratch_dirpath = bcbio_work_dirpath.replace(analysis_root_dirpath, scratch_root_dirpath)
        assert work_scratch_dirpath != bcbio_work_dirpath, (work_scratch_dirpath, bcbio_work_dirpath)
        safe_mkdir(dirname(work_scratch_dirpath))
        if isdir(work_scratch_dirpath):
            move_existing_work_dir_to = work_scratch_dirpath + '-' + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
            info('Work dir ' + work_scratch_dirpath + ' in scratch already exists, backing it up into ' + move_existing_work_dir_to)
            try:
                os.rename(work_scratch_dirpath, move_existing_work_dir_to)
            except OSError:
                err('Cannot move "work" into scratch: it exists and cannot backup the existing one')
                raise
        assert work_scratch_dirpath != bcbio_work_dirpath, (work_scratch_dirpath, bcbio_work_dirpath)
        safe_mkdir(dirname(work_scratch_dirpath))
        info('Moving work directory to scratch: ' + bcbio_work_dirpath + ' -> ' + work_scratch_dirpath)
        shutil.move(bcbio_work_dirpath, work_scratch_dirpath)
        os.symlink(work_scratch_dirpath, bcbio_work_dirpath)
        info('Symlinked work directory ' + bcbio_work_dirpath + ' -> ' + work_scratch_dirpath)


def clean_up_bcbio_dir(bcbio_proj):
    set_up_log(safe_mkdir(bcbio_proj.postproc_log_dir), 'postproc.log')

    def _move(src_fpath, dst_fpath):
        safe_mkdir(dirname(dst_fpath))
        info('Moving ' + src_fpath + ' to ' + dirname(dst_fpath))
        try:
            os.rename(src_fpath, dst_fpath)
        except OSError:
            critical('Cannot move ' + src_fpath + ' to ' + dst_fpath + ': dst exists, and permissions do now allow to remove it.')

    # Moving raw expression files in the date dir to expression/raw
    for fname in os.listdir(bcbio_proj.date_dir):
        if fname.startswith('combined.') or fname.startswith('annotated_combined.'):
            safe_mkdir(bcbio_proj.raw_expression_dir)
            _move(join(bcbio_proj.date_dir, fname), join(bcbio_proj.raw_expression_dir, fname))

    # Moving raw variants in the date dir to var/raw
    for fname in os.listdir(bcbio_proj.date_dir):
        if '.vcf' in fname and '.anno.filt' not in fname:
            _move(join(bcbio_proj.date_dir, fname), join(bcbio_proj.raw_var_dir, fname))

    # cleaning date dir
    if bcbio_proj.log_dir:
        for fname in listdir(bcbio_proj.date_dir):
            if fname.endswith('.log') or fname in ['project-summary.yaml', 'programs.txt', 'data_versions.csv']:
                os.rename(join(bcbio_proj.date_dir, fname), join(bcbio_proj.log_dir, fname))


def annotate_fn(proj, s, filt_cnf, genome_cfg, caller, cores_per_job=None, ignore=False):
    if ignore:
        return

    raw_vcf = s.find_raw_vcf(caller=caller)
    if not raw_vcf:
        f = info if s.phenotype == 'normal' else err
        f('No VCF for ' + s.phenotype + ' sample, skipping variant annotation')
        return

    output_dir = safe_mkdir(join(s.dirpath, BcbioProject.varannotate_dir))

    # Checking if bcbio-generated VCF already annotated with vcfanno with our config
    need_running = True
    with open_gzipsafe(raw_vcf) as f:
        for l in f:
            if not l.startswith('##'):
                break
            if 'cosmic_id' in l:
                need_running = False

    work_dir = safe_mkdir(join(proj.work_dir, BcbioProject.varannotate_dir, s.name))

    if can_reuse(s.find_annotated_vcf(caller=caller), raw_vcf):
        info('Annotated VCF exists, reusing')
        return

    debug('Removing rejected records...')
    pass_vcf_fpath = clean_up_vcf(work_dir, filt_cnf, raw_vcf, s.name, caller)
    passed_records = sum(1 for l in open_gzipsafe(pass_vcf_fpath) if not l.startswith('#'))
    if passed_records == 0:
        err('No passed variants found in ' + raw_vcf + '. Skipping annotation and filtering.')
        return

    if need_running:
        info('Needs annotation, running vcfanno')
        work_vcf = run_vcfanno(work_dir, pass_vcf_fpath, genome_cfg, threads=cores_per_job)
    else:
        work_vcf = pass_vcf_fpath
        info('File annotated in bcbio, skipping')

    finialize_annotate_file(work_vcf, output_dir, s.name, caller=caller)


def filtering_fn(proj, s, filt_cnf, dbsnp_multi_mafs, caller, cores_per_job=None, ignore=False):
    if ignore:
        return

    if not s.find_annotated_vcf(caller=caller):
        f = info if s.phenotype == 'normal' else err
        f('No VCF for ' + s.phenotype + ' sample, skipping variant annotation')
        return
    work_dir = safe_mkdir(join(proj.work_dir, BcbioProject.varfilter_dir, s.name))
    output_dir = safe_mkdir(join(s.dirpath, BcbioProject.varfilter_dir))
    run_filtering(
        work_dir, output_dir, s.find_annotated_vcf(caller=caller), s.name,
        filt_cnf, proj.project_name, s.genome_build, dbsnp_multi_mafs, caller, is_wgs=proj.is_wgs)


def _visualization(view, proj, filt_cnf, key_gene_names, skip_ngb=False):
    call_vis_html_fpath = None
    combined_ngs_rep_html_fpath = None

    if not proj.is_rnaseq:
        info('Creating call visualizations (linear and circos plots) for all samples...')
        call_vis_html_fpath = make_circos_and_linear_plots(proj, key_gene_names, parallel_view=view)
        info()

        info('Getting coordinates for key/target genes')
        cds_bed_file = join(proj.work_dir, 'target_genes_cds.bed')
        if is_small_target(proj.coverage_bed):
            if not can_reuse(cds_bed_file, proj.coverage_bed):
                extract_features_for_genes(proj.genome_build, proj.get_target_genes(), cds_bed_file)
        else:
            if verify_file(get_key_genes_bed(proj.genome_build)):
                run('gunzip -c ' + get_key_genes_bed(proj.genome_build), output_fpath=cds_bed_file)
            else:
                extract_features_for_genes(proj.genome_build, proj.get_target_genes(get_key_genes_file), cds_bed_file)

        info('Building NGS reports...')
        view.run(_ngs_report_fn, [[
            proj, s, filt_cnf, view.cores_per_job, call_vis_html_fpath, cds_bed_file, skip_ngb
        ] for s in proj.samples])

        info('Building combined mutations NGS report...')
        combined_ngs_rep_html_fpath = _combined_mut_ngs_report(proj, view, filt_cnf, skip_ngb)

    return combined_ngs_rep_html_fpath, call_vis_html_fpath


def _combined_mut_ngs_report(proj, view, filt_cnf, skip_ngb=False):
    work_dir = join(proj.work_dir, BcbioProject.ngs_report_name + '_combined')

    paths = proj.find_mutation_files(passed=True, caller=MAIN_CALLER)
    if not paths:
        err('No mutation files found, cannot generate the combined NGS report')
        return None

    mutations_fpath = paths[0]
    key_or_target_gene_names = proj.get_target_genes(get_key_genes_file)

    full_clin_info = ClinicalExperimentInfo(
        genome=proj.genome_build,
        key_or_target_gene_names=proj.get_target_genes(get_key_genes_file),
        bed_fpath=proj.coverage_bed,
        is_small_target=proj.is_small_target(),
        coverage_interval=proj.coverage_interval,
        mutations_fpath=mutations_fpath,
        project_name=proj.project_name,
        work_dir=safe_mkdir(work_dir)
    )
    ngs_rep_by_sample = {
        s.name: relpath(s.find_ngs_report(silent=True), join(proj.date_dir, BcbioProject.reports_dir))
            if s.find_ngs_report(silent=True) else None
        for s in proj.samples
    }
    html_fpath = join(proj.date_dir, BcbioProject.reports_dir, 'combined_mutations.html')
    GeneralMutationReporting(
            work_dir, proj.samples, proj.genome_build, full_clin_info, key_or_target_gene_names,
            filt_cnf['min_freq'], filt_cnf['act_min_freq'], ngs_rep_by_sample, skip_ngb=skip_ngb)\
        .write_report(html_fpath)
    return html_fpath


def _ngs_report_fn(proj, s, filt_cnf, cores_per_job, call_vis_html_fpath, cds_bed_file, skip_ngb=False):
    work_dir = join(proj.work_dir, BcbioProject.ngs_report_name, s.name)
    output_dir = join(proj.date_dir, BcbioProject.reports_dir)

    caller = next((c for c in CALLER_PRIORITY if c in s.variantcallers), None)
    var_txt_fpaths = s.find_mutation_files(passed=False, caller=caller) if caller else None
    if not var_txt_fpaths:
        var_txt_fpath = None
    elif len(var_txt_fpaths) == 1:
        var_txt_fpath = var_txt_fpaths[0]
    else:
        single_var_fpath, paired_var_fpath = var_txt_fpaths
        if s.normal_match:
            var_txt_fpath = paired_var_fpath
        else:
            var_txt_fpath = single_var_fpath
    if var_txt_fpath:
        mutations_fpath = verify_file(add_suffix(var_txt_fpath, vf.mut_pass_suffix))
    else:
        mutations_fpath = None

    key_or_target_gene_names = proj.get_target_genes(get_key_genes_file)

    ngs_report_path = join(proj.date_dir, BcbioProject.reports_dir, s.name + '.html')
    if can_reuse(ngs_report_path, [proj.find_seq2c_filt_file(), s.find_sv_tsv(), mutations_fpath]):
        return

    clin_info = ClinicalExperimentInfo(
        genome=s.genome_build,
        threads=cores_per_job,
        sample=NgsReportSample(
            s.name,
            output_dir,
            phenotype=s.phenotype,
            work_dir=work_dir,
            targqc_dirpath=output_dir,
            clinical_report_dirpath=output_dir,
            normal_match=s.normal_match
        ),
        key_or_target_gene_names=key_or_target_gene_names,
        cds_bed_file=cds_bed_file,
        bam_fpath=s.bam,
        bed_fpath=proj.coverage_bed,
        is_small_target=proj.is_small_target(),
        coverage_interval=proj.coverage_interval,
        mutations_fpath=mutations_fpath,
        vcf_fpath=s.find_raw_vcf(silent=True, caller=caller) if caller else None,
        filt_vcf_fpath=s.find_filt_vcf(caller=caller) if caller else None,
        sv_fpath=s.find_sv_tsv(),
        sv_vcf_fpath=s.find_sv_vcf(),
        seq2c_tsv_fpath=proj.find_seq2c_filt_file(),
        cnvkit_tsv_fpath=proj.find_cnvkit_filt_file(),
        project_name=proj.project_name,
        project_report_path=join(proj.date_dir, BcbioProject.multiqc_report_name),
        avg_depth=s.get_avg_depth(),
        bcbio_coverage_file=verify_file(join(s.dirpath, 'qc', 'coverage', s.name + '_coverage.bed')),
        circos_data_fpath=verify_file(join(proj.work_dir, s.name, BcbioProject.ngs_report_name, 'circos_data.json')),
        work_dir=safe_mkdir(join(proj.work_dir, s.name, BcbioProject.ngs_report_name)),
        caller=caller,
    )
    make_clinical_report(work_dir, s.genome_build,
        clin_info, key_or_target_gene_names, join(output_dir, s.name + '.html'),
        filt_cnf['min_freq'], filt_cnf['act_min_freq'], call_vis_html_fpath, skip_ngb=skip_ngb)


def _symlink_cnv(proj):
    """ Moves cnv files around, and sets fields
    """
    cnv_summary_dir = join(proj.date_dir, BcbioProject.cnv_dir)

    for sample in proj.samples:
        sample_cnv_dirpath = join(sample.dirpath, BcbioProject.cnv_dir)

        for cnv_caller in ['-seq2c', '-seq2c-coverage.tsv', '-cnvkit', '-cn_mops',
                           '-lumpy', '-manta', '-wham', '-sv-prioritize', '-metasv']:
            for fname in os.listdir(sample.dirpath):
                if cnv_caller in fname:
                    # Copy to <sample>/cnv
                    safe_mkdir(sample_cnv_dirpath)
                    try:
                        os.rename(join(sample.dirpath, fname), join(sample_cnv_dirpath, fname))
                    except OSError:
                        verify_file(join(sample.dirpath, fname))
                        verify_dir(sample_cnv_dirpath)
                        err(format_exc())
                        info()

            if isdir(sample_cnv_dirpath):
                for fname in os.listdir(sample_cnv_dirpath):
                    if cnv_caller in fname:
                        # Symlink to <datestamp>/cnv/<cnvcaller>
                        dst_dirpath = cnv_summary_dir
                        dst_fname = fname
                        if sample.name not in fname:
                            dst_fname = sample.name + '.' + dst_fname
                        dst_fpath = join(dst_dirpath, dst_fname)
                        try:
                            safe_mkdir(dst_dirpath)
                            if islink(dst_fpath):
                                os.unlink(dst_fpath)
                            symlink_plus(join(sample_cnv_dirpath, fname), dst_fpath)
                        except OSError:
                            err(format_exc())
                            info()

    # Merging all Seq2C into one
    seq2c_files = [_f for _f in (s.find_seq2c_file() for s in proj.samples) if _f]
    if not seq2c_files:
        debug('No Seq2C results found.')
    else:
        info('Merging Seq2C files')
        merged_seq2c_fpath = join(safe_mkdir(cnv_summary_dir), BcbioProject.seq2c_fname)
        if can_reuse(merged_seq2c_fpath, cmp_f=seq2c_files):
            with file_transaction(proj.work_dir, merged_seq2c_fpath) as tx:
                with open(tx, 'w') as out:
                    for file_index, fpath in enumerate(seq2c_files):
                        with open(fpath) as inp:
                            for line_index, l in enumerate(inp):
                                if file_index == 0 or line_index > 0:
                                    out.write(l)
    # Merging CNVkit into one
    merge_cnvkit_files(proj)
    return seq2c_files


def final_email_notification(html_report_url, html_report_fpath, jira_url, analysis_dir, project_name):
    subj = project_name
    txt = 'Post-processing finished for ' + project_name + '\n'
    txt += '\n'
    txt += 'Path: ' + analysis_dir + '\n'
    if convert_gpfs_path_to_url(analysis_dir):
        txt += 'URL: ' + convert_gpfs_path_to_url(analysis_dir) + '\n'
    if html_report_url:
        txt += 'Report: ' + html_report_url + '\n'
    txt += 'Report path: ' + html_report_fpath + '\n'
    if jira_url:
        txt += 'Jira: ' + jira_url
    send_email(txt, subj)


def change_permissions(path):
    try:
        run('chmod -R g+w ' + path)
    except Exception:
        debug('Cannot change permissions to ' + path)
        debug(traceback.format_exc())


def _save_run_cnf(proj, filt_cnf):
    if filt_cnf is None:
        return
    filt_cnf_fpath = filt_cnf.get('filt_cnf_fpath')
    pre_existing_fpath = detect_run_info_in_config_dir(proj.config_dir)
    if filt_cnf_fpath:
        if pre_existing_fpath:
            if samefile(filt_cnf_fpath, pre_existing_fpath):
                return
            else:
                pre_existing_fpath_bak = pre_existing_fpath + '.bak'
                if isfile(pre_existing_fpath_bak):
                    os.remove(pre_existing_fpath_bak)
                os.rename(pre_existing_fpath, pre_existing_fpath + '.bak')
        shutil.copyfile(filt_cnf_fpath, join(proj.config_dir, basename(filt_cnf_fpath)))
    if not filt_cnf_fpath:
        if pre_existing_fpath:
            pre_existing_fpath_bak = pre_existing_fpath + '.bak'
            if isfile(pre_existing_fpath_bak):
                os.remove(pre_existing_fpath_bak)
            os.rename(pre_existing_fpath, pre_existing_fpath + '.bak')
        filt_cnf['filt_cnf_fpath'] = save_filt_cnf(filt_cnf, proj.config_dir)
