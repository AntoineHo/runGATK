#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Small GATK alignment and variant calling pipeline using python"""

#import errno
import os
import sys
import shutil
import datetime
from time import localtime, strftime
import argparse
#import binascii
import gzip

import subprocess
from multiprocessing import Pool, TimeoutError

#import pandas as pd
import numpy as np
from Bio import SeqIO # Need BIOPYTHON SEQ/IO









#               _        _____    _____   _   _
#       /\     | |      |_   _|  / ____| | \ | |
#      /  \    | |        | |   | |  __  |  \| |
#     / /\ \   | |        | |   | | |_ | | . ` |
#    / ____ \  | |____   _| |_  | |__| | | |\  |
#   /_/    \_\ |______| |_____|  \_____| |_| \_|
#
#

def align(args) :
    """Pipeline to align reads to reference with automatic reads group"""
    """Uses BWA MEM / SAMBAMBA / """
    # 0. parse arguments & check files
    #print(args.R1[0])
    R1_reads = check_files(args.R1[0]) # Get R1 reads per lane
    R2_reads = check_files(args.R2[0]) # Get R2 reads per lane
    ref = check_files(args.Reference)[0]
    out = args.Sample[0]
    if os.path.isdir(out) :
        print("WARNING: Output directory already exists! {}".format(out))
    else :
        os.makedirs(out) # Create directory following path
    out = os.path.abspath(out)

    dc_args = {"threads":args.threads[0], "ram":args.ram[0], "mapq":args.mapq[0], "sample":args.Sample[0]} # Dict to make commands

    keep = args.keep_temp

    print("# runGATK.py align")
    print("Reference genome:\t{}\n".format(ref))
    print("R1 reads:\t{}".format(" | ".join(x for x in R1_reads)))
    print("R2 reads:\t{}\n".format(" | ".join(x for x in R2_reads)))
    print("Output stored in:\t{}\n".format(out))
    print("Keeping intermediate files:\t{}".format(keep))
    print("Other arguments:\t" + str(dc_args))
    print("===============================================================================\n")

    # 1. ALIGN READS

    # Index the reference
    indexed_ref = ref + ".bwt"
    cmd = "bwa index {ref}".format(**{"ref":ref})
    if not os.path.isfile(indexed_ref) :
        print(cmd + "\n")
        run(cmd)
    else :
        print("SKIP: bwa index file found: {}".format(indexed_ref))

    # Loop through lanes
    #print("Starting alignment...")
    i = 0
    bam_to_merge = []
    intermediate_files = []
    for r1, r2 in zip(R1_reads, R2_reads) :
        #print("Aligning paired: {} & {}".format(r1, r2))
        # Get read group for alignment
        header, id = get_read_group(r1)
        readgroup = "@RG\\tID:{id}\\tSM:{sample}\\tLB:{id}_{sample}\\tPL:{sequencer}"
        dc_rg = {"id":id, "sample":dc_args["sample"], "sequencer":"ILLUMINA"}
        readgroup = readgroup.format(**dc_rg)

        # Preparing outputs
        out_bwa_sam = os.path.join(out, dc_args["sample"] + "_" + os.path.basename(r1) + "_" + os.path.basename(r1) + ".sam")
        out_view = os.path.join(out, dc_args["sample"] + "_" + os.path.basename(r1) + "_" + os.path.basename(r1) + ".view.bam")
        out_sort_first = os.path.join(out, dc_args["sample"] + "_" + os.path.basename(r1) + "_" + os.path.basename(r1) + ".sorted.bam")

        # 1. BWA mem command
        dc_bwa = {"threads":dc_args["threads"], "readgroup":readgroup, "output":out_bwa_sam, "ref":ref, "r1":r1, "r2":r2}
        cmd = "bwa mem -K 100000000 -t {threads} -R \"{readgroup}\" -o {output} {ref} {r1} {r2}"
        cmd = cmd.format(**dc_bwa)
        if not os.path.isfile(out_bwa_sam) :
            print("\n{}: Aligning reads".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
            print(cmd + "\n")
            run(cmd)
        else :
            print("SKIP: alignement sam file found: {}".format(out_bwa_sam))
        intermediate_files.append(out_bwa_sam)

        # 2. Sambamba view command
        dc_view = {"threads":dc_args["threads"], "mapq":dc_args["mapq"], "output":out_view, "input_sam":out_bwa_sam}
        cmd = "sambamba view -h -S -f bam -l 7 -t {threads} -F \"mapping_quality >= {mapq}\" -o {output} {input_sam}"
        cmd = cmd.format(**dc_view)
        if not os.path.isfile(out_view) :
            print("\n{}: Converting to BAM file".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
            print(cmd + "\n")
            run(cmd)
        else :
            print("SKIP: bam file found: {}".format(out_view))
        intermediate_files.append(out_view)

        # 3. sort intermetdiate
        dc_sort = {"threads":dc_args["threads"], "output":out_sort_first, "input":out_view}
        cmd = "sambamba sort -l 7 -t {threads} -o {output} {input}"
        cmd = cmd.format(**dc_sort)
        if not os.path.isfile(out_sort_first) :
            print("\n{}: Sorting bam file".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
            print(cmd + "\n")
            run(cmd)
        else :
            print("SKIP: bam file found: {}".format(out_sort_first))
        intermediate_files.append(out_sort_first)
        bam_to_merge.append(out_sort_first)

    # Preparing outputs
    out_merge = os.path.join(out, dc_args["sample"] + ".merged.bam")
    out_markdup = os.path.join(out, dc_args["sample"] + ".markdup.bam")
    out_sort = os.path.join(out, dc_args["sample"] + ".sorted.CALL.bam")

    # Merging all alignments
    if len(bam_to_merge) > 1 :
        input_bams = " ".join(x for x in bam_to_merge)
        dc_merge = {"threads":dc_args["threads"], "output":out_merge, "input_bams":input_bams}
        cmd = "sambamba merge -l 7 -t {threads} {output} {input_bams}"
        cmd = cmd.format(**dc_merge)
        if not os.path.isfile(out_merge) :
            print("\n{}: Merging alignments".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
            print(cmd + "\n")
            run(cmd)
        else :
            print("SKIP: merged bam file found: {}".format(out_merge))
        intermediate_files.append(out_merge)
    else :
        print("SKIP: Merging is not necessary, only one sample!")
        # In case only one sample change out_merge for the rest of the pipeline
        out_merge = bam_to_merge[0]

    # Marking duplicates
    dc_markdup = {"threads":dc_args["threads"], "output":out_markdup, "input":out_merge}
    cmd = "sambamba markdup --overflow-list-size 600000 --hash-table-size 500000 -l 7 -t {threads} {input} {output}"
    cmd = cmd.format(**dc_markdup)
    if not os.path.isfile(out_markdup) :
        print("\n{}: Marking duplicates".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        print(cmd + "\n")
        run(cmd)
    else :
        print("SKIP: markdup bam file found: {}".format(out_markdup))
    intermediate_files.append(out_markdup)

    # Sorting final file
    dc_sort = {"threads":dc_args["threads"], "output":out_sort, "input":out_markdup}
    cmd = "sambamba sort -l 7 -t {threads} -o {output} {input}"
    cmd = cmd.format(**dc_sort)
    if not os.path.isfile(out_sort) :
        print("\n{}: Sorting bam file".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        print(cmd + "\n")
        run(cmd)
    else :
        print("SKIP: final file found: {}".format(out_sort))

    # Indexing final file
    dc_index = {"threads":dc_args["threads"], "input":out_sort}
    cmd = "sambamba index -t {threads} {input}"
    cmd = cmd.format(**dc_index)
    if not os.path.isfile(out_sort + ".bai") :
        print("\n{}: Indexing final bam file".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        print(cmd + "\n")
        run(cmd)
    else :
        print("SKIP: bam index file found: {}".format(out_sort + ".bai"))

    # IF REMOVE TEMP : remove sam file, filtered bam file, merged file, marked file
    # Keep only filtered, merged, marked and sorted final file
    if not keep :
        print("\n{} Removing:\n".format(strftime("%Y-%m-%d %H:%M:%S", localtime()) + "\n".join(["- " + str(x) for x in intermediate_files])))
        for file in intermediate_files :
            if os.path.isfile(file) :
                os.remove(file)

    sys.exit(0)
















#    _____    _____    ______   _____               _____    ______
#   |  __ \  |  __ \  |  ____| |  __ \      /\     |  __ \  |  ____|
#   | |__) | | |__) | | |__    | |__) |    /  \    | |__) | | |__
#   |  ___/  |  _  /  |  __|   |  ___/    / /\ \   |  _  /  |  __|
#   | |      | | \ \  | |____  | |       / ____ \  | | \ \  | |____
#   |_|      |_|  \_\ |______| |_|      /_/    \_\ |_|  \_\ |______|
#
#

def prepare_genome(args) :
    """Pipeline to prepare reference genome files for variant calling"""

    # 0. parse arguments & check files
    ref = check_files(args.Reference)[0]

    dc_args = {"window_size":args.window_size[0]*1000}

    print("# runGATK.py prepare")
    print("Reference genome:\t{}\n".format(ref))
    print("Other arguments:\t" + str(dc_args))
    print("===============================================================================\n")

    # 1.1 Create a bed file with intervals <- haplotype caller intervals
    out_intervals = ref + ".intervals.bed"
    if not os.path.isfile(out_intervals) :
        print("\n{}: Making intervals .bed file".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        windows = []
        for record in SeqIO.parse(ref, "fasta"):
            sid = record.id
            step = dc_args["window_size"]
            if len(record) < step :
                windows.append("{}\t{}\t{}\n".format(sid, 0, len(record)))
                continue
            else :
                all_steps = np.arange(0, len(record), step)
                for x in range(len(all_steps)-1) :
                    windows.append("{}\t{}\t{}\n".format(sid, all_steps[x], all_steps[x+1]))
                if len(record) > all_steps[-1] :
                    windows.append("{}\t{}\t{}\n".format(sid, all_steps[-1], len(record)))

        f = open(out_intervals, "w")
        for line in windows :
            f.write(line)
        f.close()
    else :
        print("SKIP: Intervals file found: {}".format(out_intervals))

    # 1.2 Create a bed file with chromosomes <- GenotypeGVCFS intervals
    out_chromosomes = ref + ".chromosomes.bed"
    if not os.path.isfile(out_chromosomes) :
        print("\n{}: Making chromosomes .bed file".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        f = open(out_chromosomes, "w")
        for record in SeqIO.parse(ref, "fasta"):
            f.write("{}\t0\t{}\n".format(record.id, len(record.seq)))
        f.close()
    else :
        print("SKIP: Chromosomes .bed file found: {}".format(out_chromosomes))

    # 2. Create sequence dictionnary with picard tools
    out_dict = os.path.splitext(ref)[0] + ".dict"
    dc_dict = {"input":ref, "output":out_dict}
    cmd = "picard CreateSequenceDictionary R={input} O={output}"
    cmd = cmd.format(**dc_dict)
    if not os.path.isfile(out_dict) :
        print("\n{}: Making reference dictionary file".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        print(cmd + "\n")
        run(cmd)
    else :
        print("SKIP: Reference dictionary file found: {}".format(out_dict))

    # 3. Indexing fasta file
    out_index = ref + ".fai"
    dc_index = {"input":ref}
    cmd = "samtools faidx {input}"
    cmd = cmd.format(**dc_index)
    if not os.path.isfile(out_index) :
        print("\n{}: Indexing reference".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        print(cmd + "\n")
        run(cmd)
    else :
        print("SKIP: Reference index found: {}".format(out_index))

    sys.exit(0)















#     _____              _        _
#    / ____|     /\     | |      | |
#   | |         /  \    | |      | |
#   | |        / /\ \   | |      | |
#   | |____   / ____ \  | |____  | |____
#    \_____| /_/    \_\ |______| |______|
#
#

def call(args) :
    """Pipeline to obtain the g.VCF files per sample with HaplotypeCaller"""
    # 0. parse arguments & check files
    # Get reference, dict and intervals
    ref = check_files(args.Reference)[0]
    refdict = check_files([os.path.splitext(ref)[0] + ".dict"])[0]
    intervals = check_files([ref + ".intervals.bed"])[0]

    # Check if valid directory and if exist raise a warning
    out = os.path.join(args.Sample[0], "call")
    if os.path.isdir(out) :
        print("WARNING: Output directory already exists! {}".format(out))
    else :
        os.makedirs(out) # Create directory following path
    out = os.path.abspath(out)

    # Get sample name
    sample_name = args.Sample[0]

    # Find sample dir and bam alignment
    sample_dir = args.Sample[0]
    if not os.path.isdir(sample_dir) :
        raise Exception("ERROR: Sample directory does not exist! {}".format(out))
    sample_dir = os.path.abspath(sample_dir)

    sample_files = os.listdir(sample_dir)
    for file in sample_files :
        if file.find(".sorted.CALL.bam") != -1 and file.find(".bai") == -1 :
            bam = os.path.abspath(os.path.join(sample_dir, file))
        else :
            raise Exception("ERROR: Could not find alignment file: {}!".format(file))
    bam = check_files([bam])[0]

    # Get other arguments
    dc_args = {"nproc":args.processes[0], "pcr":args.pcr_indel_model[0],
               "het":args.heterozygosity[0], "java":args.java_options[0],
               "hmm":args.hmm_threads[0], "mrpas":args.max_reads_per_alignment_start[0],
               "bamout":args.bam_out, "indel_het":args.indel_heterozygosity[0],
               "fid":args.founder_id[0], "ERC":args.emit_ref_confidence[0],
               "om":args.output_mode[0]}

    # Check "" in the java options and remove them if necessary
    if dc_args["java"][0] == '"' :
        dc_args["java"] = dc_args["java"][1:]
    if dc_args["java"][-1] == '"' :
        dc_args["java"] = dc_args["java"][:-1]

    keep = args.keep_temp

    print("# runGATK.py call")
    print("Reference genome:\t{}\n".format(ref))
    print("Reference dictionary:\t{}\n".format(os.path.splitext(ref)[0] + ".dict"))
    print("Sample name:\t\t{}\n".format(sample_name))
    print("Sample alignment:\t{}\n".format(bam))
    print("Output stored in:\t{}\n".format(out))
    print("Keeping intermediate files: {}".format(keep))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    if args.dry_run :
        sys.exit(0)

    # Final file path
    merge_out = os.path.join(out, "merged.g.vcf")

    # 1. Make a queue of jobs to run in parallel
    # Create a list of jobs
    cmd = "gatk HaplotypeCaller --java-options \"{java}\" -R {ref} -I {bam} -O {subout} -L {subinterval} --emit-ref-confidence {ERC} --output-mode {om} --pcr-indel-model {pcr} --heterozygosity {het} --indel-heterozygosity {indel_het} --max-reads-per-alignment-start {mrpas} --sample-name {sample} --native-pair-hmm-threads {hmm}"
    dc_haplo = {"java":dc_args["java"], "pcr":dc_args["pcr"], "het":dc_args["het"], "hmm":dc_args["hmm"], "ERC":dc_args["ERC"], "om":dc_args["om"],
                "indel_het":dc_args["indel_het"], "mrpas":dc_args["mrpas"], "fid":dc_args["fid"], "ref":ref, "bam":bam, "sample": sample_name,}

    # Add bamout to cmd if True
    if dc_args["bamout"] :
        cmd +=  " --bam-output {bamout}"

    # Add founder id if available
    if dc_args["fid"] != "" :
        cmd +=  " --founder-id {fid}"

    jobs = []
    subfiles_out = []
    interval_files = []
    for n, line in enumerate(open(intervals, "r")) :
        # Set filepaths
        subinterval = os.path.join(out, "sub_" + str(n) + ".bed")
        subout = os.path.join(out, "sub_" + str(n) + ".g.vcf")

        # Add bamout to cmd if True
        if dc_args["bamout"] :
            subbam = os.path.join(out, "sub_" + str(n) + ".realigned.bam")
            dc_haplo["bamout"] = subbam

        logfile = os.path.join(out, "sub_" + str(n) + ".log")
        # Write interval file
        f = open(subinterval, "w")
        f.write(line)
        f.close()
        # Add paths to command dictionary
        dc_haplo["subinterval"] = subinterval
        dc_haplo["subout"] = subout
        # If g.vcf of subfile does not exist and final file does not exist
        if not os.path.isfile(subout) and not os.path.isfile(merge_out):
            job = cmd.format(**dc_haplo)
            jobs.append([job, logfile, n])
        elif os.path.isfile(merge_out) : # in case final file found
            continue
        else : # in case subfile found
            print("SKIP: g.VCF file already found: {}".format(subout))

        # Print if final file is found in dir
        if os.path.isfile(merge_out) :
            print("SKIP: Final file was found in directory!")

        subfiles_out.append(subout) # For later merging & removing
        interval_files.append(subinterval) # For later removing

    print("\n{}: Created {} jobs for GATK4 HaplotypeCaller!".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), len(jobs)))
    p = Pool(dc_args["nproc"])
    p.map(run_HC, jobs)

    # 2. Gather all g.VCFs files back to one main sample G.VCF file
    if len(subfiles_out) != 1 and not os.path.isfile(merge_out):
        print("\n{}: Merging HaplotypeCaller subfiles results!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        # Write g.vcf file list to "all_g_vcf.list" file
        list_file = os.path.join(out, "all_g_vcf.list")
        f = open(list_file, "w")
        for file in subfiles_out :
            f.write(file + "\n")
        f.close()

        # Making and running merge command
        dc_merge = {"java":dc_args["java"], "listfile":list_file, "dict":refdict, "out":merge_out}
        cmd = "picard MergeVcfs {java} I={listfile} D={dict} O={out}"
        cmd = cmd.format(**dc_merge)
        print(cmd + "\n")
        run(cmd)

    elif len(subfiles_out) == 1 and not os.path.isfile(merge_out) :
        print("SKIP: No merging required! Copying output file...")
        shutil.copyfile(subfiles_out[0], merge_out)

    else :
        print("SKIP: Final file was found in directory!")

    print("\n{}: Finished calling with GATK4 HaplotypeCaller!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))

    # 3. Cleanup intermediate files
    if not keep :
        print("\n{}: Removing temporary files!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        for file in subfiles_out :
            # Remove file
            if os.path.isfile(file) :
                os.remove(file)
            # Remove index file
            if os.path.isfile(file + ".idx") :
                os.remove(file + ".idx")

        for file in interval_files :
            # Remove file
            if os.path.isfile(file) :
                os.remove(file)
            else :
                continue

    p.close() # Required so that pool stops correctly and program does not hang
    p.terminate()
    sys.exit(0)









#     _____   ______   _   _    ____    _______  __     __  _____    ______
#    / ____| |  ____| | \ | |  / __ \  |__   __| \ \   / / |  __ \  |  ____|
#   | |  __  | |__    |  \| | | |  | |    | |     \ \_/ /  | |__) | | |__
#   | | |_ | |  __|   | . ` | | |  | |    | |      \   /   |  ___/  |  __|
#   | |__| | | |____  | |\  | | |__| |    | |       | |    | |      | |____
#    \_____| |______| |_| \_|  \____/     |_|       |_|    |_|      |______|
#
#

def genotype(args) :
    """Pipeline to obtain the g.VCF files per sample with HaplotypeCaller"""
    # 0. parse arguments & check files

    # Get reference, dict and intervals
    ref = check_files(args.Reference)[0]
    refdict = check_files([os.path.splitext(ref)[0] + ".dict"])[0]
    intervals = check_files([ref + ".intervals.bed"])[0]
    chromosomes_intervals = check_files([ref + ".chromosomes.bed"])[0]

    # Get sample name
    all_samples = args.Samples[0]

    dc_samples = {}
    for sample in all_samples.split(",") :

        # Find sample and call directories
        if not os.path.isdir(sample) :
            raise Exception("ERROR: Sample directory does not exist! {}".format(sample))
        sample_dir = os.path.abspath(sample)

        # Check if valid directory and if exist raise a warning
        call_dir = os.path.join(sample_dir, "call")
        if not os.path.isdir(call_dir) :
            print("WARNING: sample directory does not contain the \"/call\" subdirectory! Did you run 'runGATK.py call'?")

        dc_samples[sample] = {"sampledir":sample_dir, "calldir":call_dir}

    # Create output directory
    out = args.Output[0]
    if os.path.isdir(out) :
        print("WARNING: Output directory already exists! {}".format(out))
    else :
        os.makedirs(out) # Create directory following path
    out = os.path.abspath(out)

    # Get other arguments
    dc_args = {"nproc":args.processes[0], "java":args.java_options[0],
               "het":args.heterozygosity[0], "indel_het":args.indel_heterozygosity[0],
               "use_combine":args.use_combine, "batch_size":args.batch_size[0],
               "sproc":args.sample_processes[0],
              }

    # Check "" in the java options and remove them if necessary
    if dc_args["java"][0] == '"' :
        dc_args["java"] = dc_args["java"][1:]
    if dc_args["java"][-1] == '"' :
        dc_args["java"] = dc_args["java"][:-1]

    keep = args.keep_temp
    not_all_sites = args.not_all

    print("# runGATK.py call")
    print("Reference genome:\t{}\n".format(ref))
    print("Reference dictionary:\t{}\n".format(os.path.splitext(ref)[0] + ".dict"))
    print("Samples:\n{}\n".format("\n".join("- " + sm + ": " + os.path.abspath(dir["sampledir"]) for sm, dir in dc_samples.items())))
    print("Output stored in:\t{}\n".format(out))
    print("Do not output all sites vcf: {}".format(not_all_sites))
    print("Keeping intermediate files: {}".format(keep))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    # 1. CombineGVCFS
    gvcfs = { sm : os.path.join(dc["calldir"], "merged.g.vcf") for sm, dc in dc_samples.items() }

    # Prepare IO files
    combined_gvcf = os.path.join(out, "combined.g.vcf")
    database = os.path.join(out, "combined.db")
    sample_map = os.path.join(out, "sample.map.txt")

    if dc_args["use_combine"] : # In case CombineGVCFS is preferred to GenomicsDBImport
        # Make and run command
        if not os.path.isfile(combined_gvcf) :
            # Prepare merged.g.vcf list to a file
            gvcfs_list = os.path.join(out, "gvcf.list")
            f = open(gvcfs_list, "w")
            for k, v in gvcfs.items() :
                f.write(v + "\n")
            f.close()

            # Prepare the command and run it
            cmd = "gatk CombineGVCFs --java-options \"{java}\" -R {ref} -V {gvcflist} -O {output}"
            dc_combine = {"java":dc_args["java"], "ref":ref, "gvcflist":gvcfs_list, "output":combined_gvcf}
            cmd = cmd.format(**dc_combine)
            print(cmd + "\n")
            run(cmd)
        else :
            print("SKIP: Combined g.vcf was found in directory!")
    else : # DEFAULT: use GenomicsDBImport
        # Make and run command
        if not os.path.isdir(database) :
            # Prepare interval & sample map file
            chromosomes = read_dict(os.path.splitext(ref)[0] + ".dict")
            intervals = os.path.join(out, "intervals.list")
            f = open(intervals, "w")
            for chrom in chromosomes.keys() :
                f.write(chrom + "\n")
            f.close()

            f = open(sample_map, 'w')
            for sample, gvcf in gvcfs.items() :
                f.write("{}\t{}\n".format(sample, gvcf))
            f.close()

            # Prepare the command and run it
            cmd = "gatk GenomicsDBImport --java-options \"{java}\" -R {ref} -L {intervals} --sample-name-map {sample_map} --genomicsdb-workspace-path {database} --batch-size {batch_size} --max-num-intervals-to-import-in-parallel {sproc}"
            dc_combine = {"java":dc_args["java"], "ref":ref, "sample_map":sample_map,
                          "database":database, "batch_size":dc_args["batch_size"],
                          "intervals":intervals, "sproc":dc_args["sproc"]}
            cmd = cmd.format(**dc_combine)
            print(cmd + "\n")
            run(cmd)
        else  :
            print("SKIP: Database was found!")

    # 2. GenotypeGVCFS
    # Make a queue of jobs to run in parallel
    # Create a list of jobs
    merge_out1 = os.path.join(out, "merged_raw.vcf")
    merge_out2 = os.path.join(out, "merged_allsites.vcf")
    cmd1 = "gatk GenotypeGVCFs --java-options \"{java}\" -R {ref} -O {subout1} -L {subinterval} -V {input} --heterozygosity {het} --indel-heterozygosity {indel_het}" # RAW
    cmd2 = "gatk GenotypeGVCFs --java-options \"{java}\" -R {ref} -O {subout2} -L {subinterval} -V {input} --heterozygosity {het} --indel-heterozygosity {indel_het} --all-sites" # All sites
    dc_gg = {"java":dc_args["java"], "ref":ref,
             "het":dc_args["het"], "indel_het":dc_args["indel_het"], }

    if dc_args["use_combine"] :
        dc_gg["input"] = combined_gvcf
    else :
        dc_gg["input"] = "gendb://" + database

    jobs = []
    subfiles_out1 = []
    subfiles_out2 = []
    interval_files = []
    for n, line in enumerate(open(chromosomes_intervals, "r")) :
        chrom_name = line.strip().split("\t")[0]
        # Set filepaths
        subinterval = os.path.join(out, "sub_" + str(chrom_name) + ".bed")
        subout1 = os.path.join(out, "sub_" + str(chrom_name) + ".raw.vcf")
        subout2 = os.path.join(out, "sub_" + str(chrom_name) + ".allsites.vcf")
        logfile1 = os.path.join(out, "sub_" + str(chrom_name) + ".raw.log")
        logfile2 = os.path.join(out, "sub_" + str(chrom_name) + ".allsites.log")

        # Write interval file
        f = open(subinterval, "w")
        f.write(line)
        f.close()

        # Add paths to command dictionary
        dc_gg["subinterval"] = subinterval
        dc_gg["subout1"] = subout1
        dc_gg["subout2"] = subout2

        # If .vcf of subfile does not exist and final file does not exist
        if not os.path.isfile(subout1) and not os.path.isfile(merge_out1):
            job1 = cmd1.format(**dc_gg)
            jobs.append([job1, logfile1, n, False])
        #elif os.path.isfile(merge_out1) : # in case final file found
        #    pass
        else : # in case subfile found
            print("SKIP: .vcf file already found: {}".format(subout1))

        # If .vcf of subfile does not exist and final file does not exist
        if not os.path.isfile(subout2) and not os.path.isfile(merge_out2):
            job2 = cmd2.format(**dc_gg)
            jobs.append([job2, logfile2, n, True])
        #elif os.path.isfile(merge_out2) : # in case final file found
        #    continue
        else : # in case subfile found
            print("SKIP: .vcf file already found: {}".format(subout2))

        # Print if final file is found in dir
        if os.path.isfile(merge_out1) :
            print("SKIP: Final raw .vcf file was found in directory!")
        if os.path.isfile(merge_out2) :
            print("SKIP: Final all sites .vcf file was found in directory!")

        subfiles_out1.append(subout1) # For later merging & removing
        subfiles_out2.append(subout2) # For later merging & removing
        interval_files.append(subinterval) # For later removing

    print("\n{}: Created {} jobs for GATK4 GenotypeGVCFs!".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), len(jobs)))
    p = Pool(dc_args["nproc"])
    p.map(run_GG, jobs)

    # 2. Gather all raw .VCFs files back to one main raw .VCF file
    if len(subfiles_out1) != 1 and not os.path.isfile(merge_out1):
        print("\n{}: Merging GenotypeGVCFs (raw) subfiles results!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        # Write g.vcf file list to "all_g_vcf.list" file
        list_file = os.path.join(out, "raw_vcf.list")
        f = open(list_file, "w")
        for file in subfiles_out1 :
            f.write(file + "\n")
        f.close()

        # Making and running merge command
        dc_merge = {"java":dc_args["java"], "listfile":list_file, "dict":refdict, "out":merge_out1}
        cmd = "picard MergeVcfs {java} I={listfile} D={dict} O={out}"
        cmd = cmd.format(**dc_merge)
        print(cmd + "\n")
        run(cmd)
    elif len(subfiles_out1) == 1 and not os.path.isfile(merge_out1) :
        print("SKIP: No merging required (raw)! Copying output file...")
        shutil.copyfile(subfiles_out1[0], merge_out1)
    else :
        print("SKIP: Final (raw) file was found in directory!")

    # 3. Gather all allsites .VCFs files back to one main allsites .VCF file
    if not not_all_sites : # In case all sites is required
        if len(subfiles_out2) != 1 and not os.path.isfile(merge_out2):
            print("\n{}: Merging GenotypeGVCFs (allsites) subfiles results!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
            # Write g.vcf file list to "all_g_vcf.list" file
            list_file = os.path.join(out, "allsites_vcf.list")
            f = open(list_file, "w")
            for file in subfiles_out2 :
                f.write(file + "\n")
            f.close()

            # Making and running merge command
            dc_merge = {"java":dc_args["java"], "listfile":list_file, "dict":refdict, "out":merge_out2}
            cmd = "picard MergeVcfs {java} I={listfile} D={dict} O={out}"
            cmd = cmd.format(**dc_merge)
            print(cmd + "\n")
            run(cmd)
        elif len(subfiles_out2) == 1 and not os.path.isfile(merge_out2) :
            print("SKIP: No merging required (allsites)! Copying output file...")
            shutil.copyfile(subfiles_out2[0], merge_out2)
        else :
            print("SKIP: Final (allsites) file was found in directory!")

    print("\n{}: Finished calling with GATK4 GenotypeGVCFs!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))

    # 3. Cleanup intermediate files
    if not keep :
        print("\n{}: Removing temporary files!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        for file in subfiles_out1 :
            # Remove file
            if os.path.isfile(file) :
                os.remove(file)
            # Remove index file
            if os.path.isfile(file + ".idx") :
                os.remove(file + ".idx")

        for file in subfiles_out2 :
            # Remove file
            if os.path.isfile(file) :
                os.remove(file)
            # Remove index file
            if os.path.isfile(file + ".idx") :
                os.remove(file + ".idx")

        for file in interval_files :
            # Remove file
            if os.path.isfile(file) :
                os.remove(file)
            else :
                continue

    p.close() # Required so that pool stops correctly and program does not hang
    p.terminate()
    sys.exit(0)


#   __      __   ____     _____   _____
#   \ \    / /  / __ \   / ____| |  __ \
#    \ \  / /  | |  | | | (___   | |__) |
#     \ \/ /   | |  | |  \___ \  |  _  /
#      \  /    | |__| |  ____) | | | \ \
#       \/      \___\_\ |_____/  |_|  \_\
#
#

def VQSR(args) :
    """Pipeline to obtain the g.VCF files per sample with HaplotypeCaller"""
    # 0. parse arguments & check files

    # Get reference, dict and intervals
    ref = check_files(args.Reference)[0]
    refdict = check_files([os.path.splitext(ref)[0] + ".dict"])[0]

    # Check if valid directory and if exist raise a warning
    out = args.Output[0]
    if os.path.isdir(out) :
        print("WARNING: Output directory already exists! {}".format(out))
    else :
        os.makedirs(out) # Create directory following path
    out = os.path.abspath(out)

    # Get VCF to recalibrate
    vcf = check_files(args.VCF)[0]

    # Get a list of resources for recalibration
    resources = []
    for res in args.resources :
        if len(res) != 4 : # Check length
            raise Exception("ERROR: Invalid resource format! {}".format(res))
        elif res[1].lower() not in ["training", "truth", "known"] : # Check type
            raise Exception("ERROR: Resource must be one of the following \"training\", \"truth\" or \"known\"")
        else :
            try : # Check prior
                prior = int(res[2])
            except :
                raise Exception("ERROR: Prior must be a valid integer <INT>! {}".format(res[2]))
            if not os.path.isfile(res[3]) :
                raise Exception("ERROR: resource VCF must be a valid file path <STRING>! {}".format(res[3]))
            else :
                istraining = "true" if res[1].lower() == "training" else "false"
                istruth = "true" if res[1].lower() == "truth" else "false"
                isknown = "true" if res[1].lower() == "known" else "false"
                resources.append({"name":res[0], "training":istraining, "truth":istruth, "known":isknown, "prior":int(res[2]), "vcf":os.path.abspath(res[3])})

    # Get the list of annotations to use
    annotations = []
    if not args.annotation :
        annotations = ["MQ", "QD", "FS", "SOR", "ReadPosRankSum", "MQRankSum"]
    else :
        annotations = args.annotation

    # Get other arguments
    dc_args = {"java":args.java_options[0]}

    # Check "" in the java options and remove them if necessary
    if dc_args["java"][0] == '"' :
        dc_args["java"] = dc_args["java"][1:]
    if dc_args["java"][-1] == '"' :
        dc_args["java"] = dc_args["java"][:-1]

    keep = args.keep_temp

    print("# runGATK.py call")
    print("Reference genome:\t{}\n".format(ref))
    print("Reference dictionary:\t{}\n".format(refdict))
    print("Output stored in:\t{}\n".format(out))
    print("Annotations used: " + ", ".join(an for an in annotations))
    for n, res in enumerate(resources) :
        print("Resource {}: Name({}), Training({}), Truth({}), Known({}), Prior({})".format(n, res["name"], res["training"], res["truth"], res["known"], res["prior"]))
    print("Keeping intermediate files: {}".format(keep))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    # 1. VariantRecalibrator
    output_recal = os.path.join(out, "vqsr.recal")
    output_tranches = os.path.join(out, "vqsr.tranches")
    output_rscript = os.path.join(out, "vqsr.rscript")
    if not os.path.isfile(output_recal) :
        dc_vq = {"output_recal":output_recal, "output_tranches":output_tranches,
                 "output_rscript":output_rscript, "ref":ref, "vcf":vcf,
                 "java":dc_args["java"]}
        cmd = "gatk VariantRecalibrator --java-options \"{java}\" -mode SNP -R {ref} -V {vcf} -O {output_recal} --tranches-file {output_tranches} --rscript-file {output_rscript}".format(**dc_vq)
        for an in annotations :
            cmd += " -an {}".format(an)
        for res in resources :
            cmd += " --resource:{name},known={known},training={training},truth={truth},prior={prior}.0 {vcf}".format(**res)
        print(cmd)
    sys.exit(0) # Exit no errors





#    _____    _   _          _    _
#   |  __ \  (_) | |        | |  | |
#   | |__) |  _  | |   ___  | |  | |  _ __
#   |  ___/  | | | |  / _ \ | |  | | | '_ \
#   | |      | | | | |  __/ | |__| | | |_) |
#   |_|      |_| |_|  \___|  \____/  | .__/
#                                    | |
#                                    |_|


def AlleleCount(args) :
    """Pipeline to obtain the g.VCF files per sample with HaplotypeCaller"""
    # 0. parse arguments & check files
    # Get reference, fai, dict and intervals
    ref = check_files(args.Reference)[0]
    refindex = check_files([ref + ".fai"])[0]
    refdict = check_files([os.path.splitext(ref)[0] + ".dict"])[0]
    intervals = check_files([ref + ".intervals.bed"])[0]

    # Check if valid directory and if exist raise a warning
    out = os.path.join(os.getcwd(), args.out)
    outdir = os.path.dirname(out)
    if os.path.isdir(outdir) :
        print("WARNING: Output directory already exists! {}".format(out))
    else :
        os.makedirs(outdir) # Create directory following path

    # Input alignment bam file
    bam = os.path.abspath(args.Bam[0])
    bam = check_files([bam])[0]

    # Get other arguments
    dc_args = {"nproc":args.processes[0], "java":args.java_options[0],
               "mbq":args.minimum_base_quality}

    # Check "" in the java options and remove them if necessary
    if dc_args["java"][0] == '"' :
        dc_args["java"] = dc_args["java"][1:]
    if dc_args["java"][-1] == '"' :
        dc_args["java"] = dc_args["java"][:-1]

    print("# runGATK.py call")
    print("Reference genome:\t{}\n".format(ref))
    print("Reference dictionary:\t{}\n".format(os.path.splitext(ref)[0] + ".dict"))
    print("Alignment:\t{}\n".format(bam))
    print("Output:\t{}\n".format(out))
    print("Other arguments: " + str(dc_args))
    print("===============================================================================\n")

    # 1. Make a queue of jobs to run in parallel
    # Create a list of jobs
    cmd = "gatk CollectAllelicCounts --java-options \"{java}\" -R {ref} -I {bam} -O {subout} -L {subinterval} --minimum-base-quality {mbq} --sequence-dictionary {dict}"
    dc_haplo = {"java":dc_args["java"], "mbq":dc_args["mbq"], "ref":ref, "dict":refdict, "bam":bam}

    jobs = []
    subfiles_out = []
    interval_files = []
    for n, line in enumerate(open(intervals, "r")) :
        # Set filepaths
        subinterval = os.path.join(outdir, "sub_" + str(n) + ".bed")
        subout = os.path.join(outdir, "sub_" + str(n) + ".sub.tsv")
        logfile = os.path.join(outdir, "sub_" + str(n) + ".log")
        # Write interval file
        f = open(subinterval, "w")
        f.write(line)
        f.close()
        # Add paths to command dictionary
        dc_haplo["subinterval"] = subinterval
        dc_haplo["subout"] = subout
        # If g.vcf of subfile does not exist and final file does not exist
        if not os.path.isfile(subout) and not os.path.isfile(out):
            job = cmd.format(**dc_haplo)
            jobs.append([job, logfile, n])
        elif os.path.isfile(out) : # in case final file found
            print("STOP: output file already found: {}".format(out))
            break
        else : # in case subfile found
            print("SKIP: output interval file already found: {}".format(subout))

        # Print if final file is found in dir
        if os.path.isfile(out) :
            print("SKIP: Final file was found in directory!")

        subfiles_out.append(subout) # For later merging & removing
        interval_files.append(subinterval) # For later removing

    print("\n{}: Created {} jobs for GATK4 CollectAllelicCounts!".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), len(jobs)))
    p = Pool(dc_args["nproc"])
    p.map(run_CAC, jobs)

    # 2. Gather all output files back to one main out file
    if len(subfiles_out) != 1 and not os.path.isfile(out):
        print("\n{}: Concatenating subfiles results!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        # Making and running merge command
        concatenate_output = {"files":" ".join(file for file in sorted(subfiles_out)), "out":out}
        cmd = "cat {files} > {out}"
        cmd = cmd.format(**concatenate_output)
        print(cmd + "\n")
        run(cmd)

    elif len(subfiles_out) == 1 and not os.path.isfile(out) :
        print("SKIP: No merging required! Copying output file...")
        shutil.copyfile(subfiles_out[0], out)

    else :
        print("SKIP: Final file was found in directory!")

    print("\n{}: Finished allele counting with GATK4 CollectAllelicCount!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))

    # 3. Cleanup intermediate files
    print("\n{}: Removing temporary files!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
    for file in subfiles_out :
        # Remove file
        if os.path.isfile(file) :
            os.remove(file)
        # Remove index file
        if os.path.isfile(file + ".idx") :
            os.remove(file + ".idx")

    for file in interval_files :
        # Remove file
        if os.path.isfile(file) :
            os.remove(file)
        else :
            continue

    p.close() # Required so that pool stops correctly and program does not hang
    p.terminate()
    print("\n{}: Done!".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
    sys.exit(0)


#     ____    _______   _    _   ______   _____     _____
#    / __ \  |__   __| | |  | | |  ____| |  __ \   / ____|
#   | |  | |    | |    | |__| | | |__    | |__) | | (___
#   | |  | |    | |    |  __  | |  __|   |  _  /   \___ \
#   | |__| |    | |    | |  | | | |____  | | \ \   ____) |
#    \____/     |_|    |_|  |_| |______| |_|  \_\ |_____/
#
#


def read_dict(dict_file) :
    chromosomes = {}
    f = open(dict_file, 'r')
    for line in f :
        if "@SQ" in line :
            s = line.strip().split("\t")
            name, length = None, None
            for d in s :
                data = d.split(":")
                if data[0] == "SN" :
                    name = data[1]
                elif data[0] == "LN" :
                    length = int(data[1])
                else :
                    continue
            if name is not None :
                chromosomes[name] = length
            else :
                pass
        else :
            continue
    return chromosomes

def run(cmd) :
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    proc.communicate()

def run_GG(job) :
    cmd = job[0]
    number = job[2]
    is_allsite = job[3]
    if is_allsite :
        print("{}: Running GenotypeGVCFs (--all-sites) job number {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), number))
    else :
        print("{}: Running GenotypeGVCFs (default) job number {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), number))
    logfile = open(job[1], "w")
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=logfile)
    proc.communicate()
    logfile.close()
    print("{}: Finished GenotypeGVCFs job number {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), number))

def run_HC(job) :
    cmd = job[0]
    number = job[2]
    print("{}: Running HaplotypeCaller job number {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), number))
    logfile = open(job[1], "w")
    logfile.write(cmd + "\n")
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=logfile)
    proc.communicate()
    logfile.close()
    print("{}: Finished HaplotypeCaller job number {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), number))

def run_CAC(job) :
    cmd = job[0]
    number = job[2]
    print("{}: Running CollectAllelicCount job number {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), number))
    logfile = open(job[1], "w")
    logfile.write(cmd + "\n")
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=logfile)
    proc.communicate()
    logfile.close()
    print("{}: Finished CollectAllelicCount job number {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), number))

def get_read_group(path) :
    try :
        f = gzip.GzipFile(path, 'rb')
        head = str(f.readline().decode('UTF-8').strip()) # read first line
        f.close()
        id = "_".join(x for x in head.split(':')[0:4]).replace("@", '')
        return head, id
    except :
        f = open(path, "r")
        head = str(f.readline().strip()) # read first line
        f.close()
        id = "_".join(x for x in head.split(':')[0:4]).replace("@", '')
        return head, id

def check_dirs(dirs) :
    """Returns absolute paths and raise exception if dir does not exist"""
    absdirs = []
    for d in dirs :
        if not os.path.isdir(d) :
            raise Exception("ERROR: {} is not found!".format(d))
        else :
            absdirs.append(os.path.abspath(d))
    return absdirs

def check_files(files) :
    """Returns absolute file paths and raise exception if file does not exist"""
    absfiles = []
    for file in files :
        if not os.path.isfile(file) :
            raise Exception("ERROR: {} is not found!".format(file))
        else :
            absfiles.append(os.path.abspath(file))
    return absfiles

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return False

def list_str(v) :
    return v.split(',')

def str_to_bool(v) :
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')



#               _____     _____   _    _   __  __   ______   _   _   _______    _____
#       /\     |  __ \   / ____| | |  | | |  \/  | |  ____| | \ | | |__   __|  / ____|
#      /  \    | |__) | | |  __  | |  | | | \  / | | |__    |  \| |    | |    | (___
#     / /\ \   |  _  /  | | |_ | | |  | | | |\/| | |  __|   | . ` |    | |     \___ \
#    / ____ \  | | \ \  | |__| | | |__| | | |  | | | |____  | |\  |    | |     ____) |
#   /_/    \_\ |_|  \_\  \_____|  \____/  |_|  |_| |______| |_| \_|    |_|    |_____/
#
#

def main() :
    """Argument parser"""
    parser = argparse.ArgumentParser(description='Runs modules of a variant calling pipeline with bwa-mem and GATK.')
    subparsers = parser.add_subparsers(required=True, dest="align || prepare || call || genotype || vqsr || allelecount")

    # Align reads to the reference
    aln = subparsers.add_parser('align', help="Aligns reads from ONE sample on a reference genome\nReads are considered Illumina Paired-End\nReads list correspond to lanes and must be given in the same order or the alignment will fail")
    aln.add_argument('Sample',nargs=1,type=str,default=['sample1'],help="<STRING> sample name for output and adding to read group. Default: 'sample1'")
    aln.add_argument('Reference',nargs=1,type=str,help="<STRING> A path to the reference fenome fasta file.")
    aln.add_argument('R1',nargs=1,type=list_str,help="<STRING> A comma separated list of FastQ files FORWARD reads sequences separated by LANES (can be compressed).")
    aln.add_argument('R2',nargs=1,type=list_str,help="<STRING> A comma separated list of FastQ files REVERSE reads sequences separated by LANES (can be compressed).")
    aln.add_argument('-t','--threads',nargs=1,type=int,default=[4],required=False,help="<INT> Maximum threads to use. Default: %(default)s")
    aln.add_argument('-r','--ram',nargs=1,type=int,default=[4],required=False,help="<INT> Number of GB of RAM per thread to use. Default: %(default)s")
    aln.add_argument('-m','--mapq',nargs=1,type=int,default=[20],required=False,help="<INT> Sambamba view filtering to keep only alignements with MAPQ >= <INT>. Default: %(default)s")
    aln.add_argument('-kp', '--keep-temp',type=str_to_bool, nargs='?', const=True, default=False, help="Do not remove intermediate steps alignments (cannot be before positional argument). Default: %(default)s")
    #aln.add_argument('Output',nargs=1,type=str,help="<STRING> An output directory path.")
    #aln.add_argument('Output',nargs=1,type=str,help="<STRING> The technology reference for read group. Default: ILLUMINA")
    #aln.add_argument('--bwa-options', '-bo',nargs=1,type=str,help="<STRING> supplementary options to pass to bwa mem.")
    #aln.add_argument('--sambamba-options', '-so',nargs=1,type=str,help="<STRING> supplementary options to pass to sambamba.")
    aln.set_defaults(func=align)

    # Prepare genome assembly for the HaplotypeCaller
    # WARNING WINDOW SIZE IN KB
    prep = subparsers.add_parser('prepare', help="Create several files from a reference to prepare for the variant calling")
    prep.add_argument('Reference',nargs=1,type=str,help="<STRING> A fasta file containing the reference genome.")
    prep.add_argument('-ws','--window-size',nargs=1,type=int,default=[100], required=False,help="<INT> Window size (in kb!!!) to split genome. Default: %(default)s")
    prep.set_defaults(func=prepare_genome)

    # Parallel HaplotypeCaller
    cal = subparsers.add_parser('call', help="Runs GATK HaplotypeCaller with a sample and the prepared reference genome assembly.")
    cal.add_argument('Sample',nargs=1,type=str,help="<STRING> sample name for output and SAME AS \"runGATK.py align\"")
    cal.add_argument('Reference',nargs=1,type=str,help="<STRING> A path to the reference genome fasta file")
    cal.add_argument('-pim', '--pcr-indel-model',nargs=1,type=str,default=['NONE'],choices=["NONE","CONSERVATIVE","HOSTILE","AGGRESSIVE"],help="<STRING> Argument to pass to HaplotypeCaller --pcr-indel-model. Default: %(default)s")
    cal.add_argument('-jo', '--java-options',nargs=1,type=str,default=['-Xmx4G'],help="<STRING> Java Virtual Machine options (Ram Per Process is defined here). Default: %(default)s")
    cal.add_argument('-he','--heterozygosity',nargs=1,type=float,default=[0.01], required=False,help="<FLOAT> Heterozygosity value to pass to HaplotypeCaller. Default: %(default)s")
    cal.add_argument('-p','--processes',nargs=1,type=int,default=[4], required=False,help="<INT> Maximum processes in pool to use. Default: %(default)s")
    cal.add_argument('-ht','--hmm-threads',nargs=1,type=int,default=[4], required=False,help="<INT> Number of threads to use for HMM. Default: %(default)s (per process).")
    cal.add_argument('-mrpas','--max-reads-per-alignment-start',nargs=1,type=int,default=[50], required=False,help="<INT> See GATK doc. Default: %(default)s (set to 0 to deactivate).")
    cal.add_argument('-ihe','--indel-heterozygosity',nargs=1,type=float,default=[0.0001], required=False,help="<FLOAT> Indel heterozygosity value to pass to HaplotypeCaller. Default: %(default)s")
    cal.add_argument('-fi', '--founder-id',nargs=1,type=str,default=[''],help="<STRING> Founder population sample ID. Default: empty")
    cal.add_argument('-erc', '--emit-ref-confidence',nargs=1,type=str,default=['BP_RESOLUTION'],choices=['BP_RESOLUTION','GVCF'],help="<STRING> ERC (see HaplotypeCaller documentation). Default: %(default)s")
    cal.add_argument('-om', '--output-mode',nargs=1,type=str,default=['EMIT_ALL_ACTIVE_SITES'],choices=['EMIT_VARIANTS_ONLY','EMIT_ALL_CONFIDENT_SITES','EMIT_ALL_ACTIVE_SITES'],help="<STRING> Output mode (see HaplotypeCaller documentation). Default: %(default)s")
    cal.add_argument('-kp', '--keep-temp',type=str_to_bool, nargs='?', const=True, default=False, help="Do not remove intermediate steps files (cannot be before positional argument). Default: %(default)s")
    cal.add_argument('-dr', '--dry-run',type=str_to_bool, nargs='?', const=True, default=False, help="Only parse arguments for testing and debugging commands. Default: %(default)s")
    cal.add_argument('-bo','--bam-out',type=str_to_bool, nargs='?', const=True, default=False, help="<INT> Output realigned reads bam. Default: %(default)s")
    cal.set_defaults(func=call)
    #cal.add_argument('Output',nargs=1,type=str,help="<STRING> An output directory path storing the result of HaplotypeCaller (one sample).")


    # AlleleCount from bam
    """
    alc = subparsers.add_parser('allelecount', help="Runs GATK CollectAllelicCounts with the prepared reference genome assembly and a bam file.")
    alc.add_argument('Bam',nargs=1,type=str,help="<STRING> A path to an alignment bam file i.e.: from \"runGATK.py align\"")
    alc.add_argument('Reference',nargs=1,type=str,help="<STRING> A path to the reference genome fasta file")
    alc.add_argument('-jo', '--java-options',nargs=1,type=str,default=['-Xmx4G'],help="<STRING> Java Virtual Machine options (Ram Per Process is defined here). Default: %(default)s")
    alc.add_argument('-p','--processes',nargs=1,type=int,default=[4], required=False,help="<INT> Maximum processes in pool to use. Default: %(default)s")
    alc.add_argument('-mbq','--minimum-base-quality',nargs=1,type=int,default=20, required=False,help="<INT> Minimum base quality. Base calls with lower quality will be filtered out of pileups. Default: %(default)s")
    alc.add_argument('-o','--out',nargs=1,type=str,default="allelecount.tsv", required=False, help="<STRING> Output filename. Default: %(default)s")
    alc.set_defaults(func=AlleleCount)
    """

    # GenotypeGVCFS
    gen = subparsers.add_parser('genotype',             help="Runs GATK joint calling using outputs of different samples calls and the prepared reference output directory.")
    gen.add_argument('Reference',                       nargs=1,             type=str,           help="<STRING> A path to the output directory of \"runGATK.py prepare\"")
    gen.add_argument('Samples',                         nargs=1,             type=str,           help="<STRING> A comma separated list of samples (directory must contain call/merged.g.vcf)")
    gen.add_argument('Output',                          nargs=1,             type=str,           help="<STRING> An output directory path for the joint calling.")
    gen.add_argument('-p','--processes',                nargs=1,type=int,    default=[4],        required=False, help="<INT> Maximum threads to use. Default: %(default)s")
    gen.add_argument('-sp','--sample-processes',        nargs=1,type=int,    default=[4],        required=False, help="<INT> Maximum threads to read samples simultaneously (GenomicsDBImport only). Default: %(default)s")
    gen.add_argument('-bs','--batch-size',              nargs=1,type=int,    default=[50],       required=False, help="<INT> Default batch size for GenomicsDBImport. Default: %(default)s")
    gen.add_argument('-jo', '--java-options',           nargs=1,type=str,    default=['-Xmx4G'], help="<STRING> Java Virtual Machine options (Ram Per Process is defined here). Default: %(default)s")
    gen.add_argument('-fi', '--founder-id',             nargs=1,type=str,    default=[''],       help="<STRING> Founder population sample ID. Default: \"\" (empty)")
    gen.add_argument('-he','--heterozygosity',          nargs=1,type=float,  default=[0.01],     required=False, help="<FLOAT> Heterozygosity value to pass to HaplotypeCaller. Default: %(default)s")
    gen.add_argument('-ihe','--indel-heterozygosity',   nargs=1,type=float,  default=[0.0001],   required=False, help="<FLOAT> Indel heterozygosity value to pass to HaplotypeCaller. Default: %(default)s")
    gen.add_argument('-na', '--not-all',                action="store_true", help="Do not use allsites option in GenotypeGVCFs (in addition to normal vcf output). Default: will output allsites vcf and normal vcf output.")
    gen.add_argument('-kp', '--keep-temp',              action="store_true", help="Do not remove intermediate steps files (cannot be before positional argument). Default: %(default)s")
    gen.add_argument("--use-combine",                   action="store_true", help="Use CombineGVCFS instead of GenomicsDBImport")
    gen.set_defaults(func=genotype)
    #gen.add_argument('-w','--windows',nargs=1,type=int,default=[100000], required=False,help="<INT> Window size (in kb!!!) to split genome. Default: 100.")


    # VariantScoreRecalibration
    """
    vqs = subparsers.add_parser('vqsr', help="Runs GATK VariantRecalibrator using truth, training and known variant sets.")
    vqs.add_argument('Reference',nargs=1,type=str,help="<STRING> A path to the output directory of \"runGATK.py prepare\"")
    vqs.add_argument('VCF',nargs=1,type=str,help="<STRING> A VCF file containing variants to apply recalibration to.")
    vqs.add_argument('Output',nargs=1,type=str,help="<STRING> An output directory path for the recalibration files.")
    vqs.add_argument('-r', '--resources',nargs=4, metavar=('NAME','{KNOWN,TRAINING,TRUE}','PRIOR','VCF'),action="append",type=str,help="<STRINGS> A resource for recalibration. Format: NAME KNOWN|TRAINING|TRUE PRIOR VCF. e.g.: hapmap KNOWN 15 ./input/hapmap.vcf")
    vqs.add_argument('-an', '--annotation',nargs="?", type=str, action="append", help="<STRINGS> Annotations to use for recalibration calculation. Default: MQ, QD, FS, SOR, ReadPosRankSum, MQRankSum")
    #vqs.add_argument('-ihe','--indel-heterozygosity',nargs=1,type=float,default=[0.0001], required=False,help="<FLOAT> Indel heterozygosity value to pass to HaplotypeCaller. Default: %(default)s")
    vqs.add_argument('-kp', '--keep-temp', type=str_to_bool, metavar="", nargs='?', const=True, default=False, help="Do not remove intermediate steps files (cannot be before positional argument). Default: %(default)s")
    vqs.add_argument('-jo', '--java-options',nargs=1,type=str,default=['-Xmx4G'],help="<STRING> Java Virtual Machine options (Ram Per Process is defined here). Default: %(default)s")
    vqs.set_defaults(func=VQSR)
    """
    #prep.add_argument('-ts','--tm-size',nargs=1,type=int,default=[15], required=False,help="<INT> Minimum aligned length to perform thermodynamics analysis. Default: 15.")
    #test.add_argument('-ma','--min-align',nargs=1,type=int,default=[15], required=False,help="<INT> Minimum BLAST aligned length to report. Default: 15.")

    args = parser.parse_args()

    args.func(args) # test

    #print("Done")
    sys.exit(0)




if __name__ == '__main__':
	main()
