#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
from multiprocessing import Pool

from utils import log, check_files, read_dict, remove_temp_files, check_args
from runners import run, run_GG
from merger import merge_all_sub_files

#     _____   ______   _   _    ____    _______  __     __  _____    ______
#    / ____| |  ____| | \ | |  / __ \  |__   __| \ \   / / |  __ \  |  ____|
#   | |  __  | |__    |  \| | | |  | |    | |     \ \_/ /  | |__) | | |__
#   | | |_ | |  __|   | . ` | | |  | |    | |      \   /   |  ___/  |  __|
#   | |__| | | |____  | |\  | | |__| |    | |       | |    | |      | |____
#    \_____| |______| |_| \_|  \____/     |_|       |_|    |_|      |______|
#
#

def genotype_samples(args) :
    """Pipeline to obtain the g.VCF files per sample with HaplotypeCaller"""

    # Get reference, dict and intervals
    ref = check_files(args.Reference)[0]
    refdict = check_files([os.path.splitext(ref)[0] + ".dict"])[0]
    intervals = check_files([ref + ".intervals.bed"])[0]
    chromosomes_intervals = check_files([ref + ".chromosomes.bed"])[0]

    # Get sample name
    all_samples = args.Samples[0]
    out = args.Output[0]

    # Create output directory
    dc_samples, out, gvcfs = find_files_and_sample_directories(all_samples, out)

    # Get other arguments
    dc_args = {
        "nproc":args.processes[0], "java":args.java_options[0],
        "het":args.heterozygosity[0], "indel_het":args.indel_heterozygosity[0],
        "batch_size":args.batch_size[0], "sproc":args.sample_processes[0],
        "chunk_size":args.max_merge_number[0], "output_logs":args.log,
        "pjo":args.picard_java_options[0], "dbjo":args.db_java_options[0],
    }
    dc_args = check_args(dc_args)

    keep = args.keep_temp
    do_all_sites = args.all_sites

    print("# runGATK.py call")
    print("Reference genome:\t{}\n".format(ref))
    print("Reference dictionary:\t{}\n".format(os.path.splitext(ref)[0] + ".dict"))
    print("Samples:\n{}\n".format("\n".join("- " + sm + ": " + os.path.abspath(dir["sampledir"]) for sm, dir in dc_samples.items())))
    print("Output stored in:\t{}\n".format(out))
    print("Outputting all sites: {}".format(do_all_sites))
    print("Keeping intermediate files: {}".format(keep))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    # 0. check for final file
    if do_all_sites :
        merge_out = os.path.join(out, "merged_allsites.vcf")
    else :
        merge_out = os.path.join(out, "merged_raw.vcf")

    if os.path.isfile(merge_out) :
        log("SKIP: found merged final file.")
        return

    # 1. CombineGVCFS
    combined_gvcf = os.path.join(out, "combined.g.vcf")
    database = os.path.join(out, "combined.db")
    sample_map = os.path.join(out, "sample.map.txt")

    # Make and run command
    if not os.path.isdir(database) :
        genomics_db_import(gvcfs, ref, refdict, out, database, sample_map, **dc_args)
    else  :
        log("SKIP: found database.")

    # 2. GenotypeGVCFS
    sub_files, interval_files = genotype_gvcfs(ref, database, chromosomes_intervals, out, all=do_all_sites, **dc_args)

    # 2. Gather all raw .VCFs files back to one main raw .VCF file
    chunks_to_remove = []
    if len(sub_files) != 1 :
        log("Merging GenotypeGVCFs sub files.")
        chunks_to_remove = merge_all_sub_files(sub_files, merge_out, out, dc_args["chunk_size"], refdict, dc_args["pjo"])
    elif len(sub_files) == 1 :
        log("SKIP: No merging required.")
        shutil.copyfile(sub_files[0], merge_out)
    else :
        # SHOULD NEVER COME HERE
        log("SKIP: found final .g.vcf.")

    log("Finished calling with GATK4 GenotypeGVCFs.")

    # 3. Cleanup intermediate files
    if not keep :
        remove_temp_files(sub_files + chunks_to_remove + interval_files)


def find_files_and_sample_directories(samples, out) :
    """Find path of {sample}.merged.g.vcf of each sample to genotype"""
    dc_samples = {}
    gvcfs = {}
    for sample in samples.split(",") :

        # Find sample and call directories
        if not os.path.isdir(sample) :
            raise Exception("ERROR: Sample directory does not exist! {}".format(sample))
        sample_dir = os.path.abspath(sample)

        # Check if valid directory and if exist raise a warning
        call_dir = os.path.join(sample_dir, "call")
        if not os.path.isdir(call_dir) :
            log("WARNING: sample directory does not contain the \"/call\" subdirectory! Did you run 'runGATK.py call'?")

        dc_samples[sample] = {"sampledir":sample_dir, "calldir":call_dir}

        gvcf = os.path.join(call_dir, "merged.g.vcf")
        if not os.path.isfile(gvcf) :
            raise Exception("Could not find .g.vcf file for sample {}.".format(sample))
        else :
            gvcfs[sample] = gvcf

    if os.path.isdir(out) :
        log("WARNING: Output directory already exists! {}".format(out))
    else :
        os.makedirs(out) # Create directory following path
    out = os.path.abspath(out)

    return dc_samples, out, gvcfs

def genomics_db_import(gvcfs, ref, refdict, outdir, database, sample_map, **kwargs) :
    """Run gatk GenomicsDBImport"""

    log("Importing samples to combined.db.")

    chromosomes = read_dict(refdict)
    intervals = os.path.join(outdir, "intervals.list")
    f = open(intervals, "w")
    for chrom in chromosomes.keys() :
        f.write(chrom + "\n")
    f.close()

    f = open(sample_map, 'w')
    for sample, gvcf in gvcfs.items() :
        f.write("{}\t{}\n".format(sample, gvcf))
    f.close()

    # Prepare the command and run it
    cmd = "gatk GenomicsDBImport --java-options \"{dbjo}\" -R {ref} -L {intervals} --sample-name-map {sample_map} --genomicsdb-workspace-path {database} --batch-size {batch_size} --max-num-intervals-to-import-in-parallel {sproc}"
    dc_combine = {"dbjo":kwargs["dbjo"], "ref":ref, "sample_map":sample_map,
                  "database":database, "batch_size":kwargs["batch_size"],
                  "intervals":intervals, "sproc":kwargs["sproc"]}
    cmd = cmd.format(**dc_combine)
    print(cmd + "\n")
    if kwargs["output_logs"] :
        log_file = os.path.join(outdir, "GenomicsDBImport.log")
        f = open(log_file, "w")
        run(cmd, ERR=f)
        f.close()
    else :
        run(cmd)

def genotype_gvcfs(ref, database, chr_intervals, out, all=False, **kwargs) :
    """Run gatk GenotypeGVCFs"""
    # Make a queue of jobs to run in parallel
    # Create a list of jobs

    log("Joint-genotyping samples.")

    cmd = "gatk GenotypeGVCFs --java-options \"{java}\" -R {ref} -O {sub} -L {sub_interval} -V {input} --heterozygosity {het} --indel-heterozygosity {indel_het}" # RAW
    if all :
        cmd += " --all-sites" # All sites

    dc_gg = {
        "ref":ref, "input": "gendb://" + database, "java":kwargs["java"],
        "het":kwargs["het"], "indel_het":kwargs["indel_het"],
    }

    jobs = []
    sub_files = []
    interval_files = []
    for n, line in enumerate(open(chr_intervals, "r")) :

        chr_name = line.strip().split("\t")[0]
        dc_gg["sub_interval"] = os.path.join(out, "sub_" + str(n) + ".bed")
        f = open(dc_gg["sub_interval"], "w")
        f.write(line)
        f.close()

        if all :
            dc_gg["sub"] = os.path.join(out, "sub_" + str(n) + ".all.vcf")
            dc_gg["log"] = os.path.join(out, "sub_" + str(n) + ".all.log")
        else :
            dc_gg["sub"] = os.path.join(out, "sub_" + str(n) + ".raw.vcf")
            dc_gg["log"] = os.path.join(out, "sub_" + str(n) + ".raw.log")

        if not kwargs["output_logs"] :
            dc_gg["log"] = False

        if not os.path.isfile(dc_gg["sub"]) :
            job = cmd.format(**dc_gg)
            jobs.append([job, dc_gg["log"], n])

        sub_files.append(dc_gg["sub"]) # For later merging & removing
        interval_files.append(dc_gg["sub_interval"]) # For later removing

    if len(jobs) > 0 :
        log("Created {} jobs for GATK4 GenotypeGVCFs.".format(len(jobs)))

        p = Pool(kwargs["nproc"])
        p.map(run_GG, jobs)
        p.close() # Required so that pool stops correctly and program does not hang
        p.terminate()
    else :
        log("SKIP: found interval files.")

    return sub_files, interval_files
