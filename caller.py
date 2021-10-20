#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
from multiprocessing import Pool, TimeoutError

from utils import log, check_files, index_bam, remove_temp_files, check_args
from runners import run, run_HC
from merger import merge_all_sub_files

# Pipeline

def bam_call(args) :
    """Pipeline to obtain the g.VCF files per sample with HaplotypeCaller"""

    # 0. parse arguments & check files
    # Get reference, dict and intervals
    ref = check_files(args.Reference)[0]
    ref_dict = check_files([os.path.splitext(ref)[0] + ".dict"])[0]
    intervals = check_files([ref + ".intervals.bed"])[0]

    # Check if valid directory and if exist raise a warning
    sample = args.Sample[0]
    sample_name = sample
    sample_dir, bam, outdir = find_files_and_sample_directories(sample, sample_name)

    # Get other arguments
    dc_args = {
        "nproc":args.processes[0], "pcr":args.pcr_indel_model[0],
        "het":args.heterozygosity[0], "java":args.java_options[0],
        "hmm":args.hmm_threads[0], "mrpas":args.max_reads_per_alignment_start[0],
        "bam_out":args.bam_out, "indel_het":args.indel_heterozygosity[0],
        "fid":args.founder_id[0], "ERC":args.emit_ref_confidence[0],
        "om":args.output_mode[0], "pjo":args.picard_java_options[0],
        "output_logs":args.log, "chunk_size":args.max_merge_number[0],
    }
    dc_args = check_args(dc_args)

    keep = args.keep_temp

    print("# runGATK.py call")
    print("Reference genome:\t{}\n".format(ref))
    print("Reference dictionary:\t{}\n".format(os.path.splitext(ref)[0] + ".dict"))
    print("Sample name:\t\t{}\n".format(sample_name))
    print("Sample alignment:\t{}\n".format(bam))
    print("Output stored in:\t{}\n".format(outdir))
    print("Keeping intermediate files: {}".format(keep))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    if args.dry_run :
        sys.exit(0)

    # Final file path
    merge_out = os.path.join(outdir, "merged.g.vcf")
    if os.path.isfile(merge_out) :
        log("SKIP: found final file.")
        return

    # 1. Make a queue of jobs to run in parallel
    # Create a list of jobs
    jobs, sub_files, interval_files = make_call_jobs(ref, bam, sample_name, intervals, outdir, merge_out, **dc_args)
    # Run jobs
    if len(jobs) > 0 :
        log("Created {} jobs for GATK4 HaplotypeCaller!".format(len(jobs)))
        p = Pool(dc_args["nproc"])
        p.map(run_HC, jobs)
        p.close() # Required so that pool stops correctly and program does not hang
        p.terminate()
    else :
        log("SKIP: No HaplotypeCaller job created.")

    # 2. Gather all g.VCFs files back to one main sample G.VCF file
    chunks_to_remove = []
    if len(sub_files) != 1 :
        log("Merging HaplotypeCaller results.")
        chunks_to_remove = merge_all_sub_files(sub_files, merge_out, outdir, dc_args["chunk_size"], ref_dict, dc_args["pjo"])
    elif len(sub_files) == 1 :
        log("SKIP: No merging required.")
        shutil.copyfile(sub_files[0], merge_out)
    else :
        # SHOULD NEVER COME HERE
        log("SKIP: found final file.")

    log("Finished calling with GATK4 HaplotypeCaller!")

    # 3. Cleanup intermediate files
    if not keep :
        remove_temp_files(sub_files + interval_files + chunks_to_remove)


def find_files_and_sample_directories(sample, sample_name) :
    """Gather relevant files and exceptions if not found"""
    # Check if valid directory and if exist raise a warning
    if not os.path.isdir(sample) :
        raise Exception("Could not find sample directory")
    else :
        sample = os.path.abspath(sample)

    alignment = os.path.join(sample, sample_name) + ".sorted.CALL.bam"
    if not os.path.isfile(alignment) :
        #print(alignment)
        raise Exception("Could not find sample .CALL.bam alignment file")

    index = alignment + ".bai"
    if not os.path.isfile(index) :
        index_bam(alignment, 1) # only one thread here

    outdir = os.path.join(sample, "call")
    if os.path.isdir(outdir) :
        log("WARNING: Output directory already exists: {}.".format(outdir))
    else :
        os.makedirs(outdir) # Create directory following path
    outdir = os.path.abspath(outdir)

    return sample, alignment, outdir


def make_call_jobs(ref, bam, sample_name, intervals, out, merge_out, **kwargs) :
    """Create calling jobs for running them"""

    cmd = "gatk HaplotypeCaller --java-options \"{java}\" -R {ref} -I {bam} -O {sub_out} -L {sub_interval} --emit-ref-confidence {ERC} --output-mode {om} --pcr-indel-model {pcr} --heterozygosity {het} --indel-heterozygosity {indel_het} --max-reads-per-alignment-start {mrpas} --sample-name {sample} --native-pair-hmm-threads {hmm}"
    dc_haplo = {
        "java":kwargs["java"], "pcr":kwargs["pcr"], "het":kwargs["het"],
        "hmm":kwargs["hmm"], "ERC":kwargs["ERC"], "om":kwargs["om"],
        "indel_het":kwargs["indel_het"], "mrpas":kwargs["mrpas"],
        "fid":kwargs["fid"],
        "ref":ref, "bam":bam, "sample": sample_name,
        }
    if kwargs["bam_out"] :
        cmd += " --bam-output {bam_out}"
    if kwargs["fid"] :
        cmd += " --founder-id {fid}"

    jobs = []
    sub_files = []
    interval_files = []

    for n, line in enumerate(open(intervals, "r")) :
        sub_interval = os.path.join(out, "sub_" + str(n) + ".bed")
        sub_out = os.path.join(out, "sub_" + str(n) + ".g.vcf")

        # Add a sub bam file
        if kwargs["bam_out"] :
            sub_bam = os.path.join(out, "sub_" + str(n) + ".realigned.bam")
            dc_haplo["bam_out"] = sub_bam

        if kwargs["output_logs"] :
            log_file = os.path.join(out, "sub_" + str(n) + ".log")
        else :
            log_file = None

        f = open(sub_interval, "w")
        f.write(line)
        f.close()

        # Add paths to command dictionary
        dc_haplo["sub_interval"] = sub_interval
        dc_haplo["sub_out"] = sub_out

        # If g.vcf of subfile does not exist -> run the job
        if not os.path.isfile(sub_out) :
            job = cmd.format(**dc_haplo)
            if kwargs["output_logs"] :
                jobs.append([job, log_file, n])
            else :
                jobs.append([job, False, n])

        sub_files.append(sub_out)
        interval_files.append(sub_interval)
        #else : # in case subfile found
        #    print("SKIP: g.VCF file already found: {}".format(subout))

    return jobs, sub_files, interval_files
