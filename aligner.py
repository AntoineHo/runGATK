#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

from runners import run, run_log
from utils import log, check_files, get_read_group, index_bam

# Pipeline

def align_reads(args) :
    """Pipeline to align reads to reference with automatic reads group"""
    """Uses BWA MEM / SAMBAMBA / """
    # 0. parse arguments & check files
    #print(args.R1[0])
    R1_reads = check_files(args.R1[0]) # Get R1 reads per lane
    R2_reads = check_files(args.R2[0]) # Get R2 reads per lane
    ref = check_files(args.Reference)[0]
    out = args.Sample[0]

    dc_args = {"threads":args.threads[0], "ram":args.ram[0], "mapq":args.mapq[0], "sample":args.Sample[0], "output_logs":args.log} # Dict to make commands
    keep = args.keep_temp

    print("# runGATK.py align")
    print("Reference genome:\t{}\n".format(ref))
    print("R1 reads:\t{}".format(" | ".join(x for x in R1_reads)))
    print("R2 reads:\t{}\n".format(" | ".join(x for x in R2_reads)))
    print("Output stored in:\t{}\n".format(out))
    print("Keeping intermediate files:\t{}".format(keep))
    print("Other arguments:\t" + str(dc_args))
    print("===============================================================================\n")

    # 0. Make output
    if os.path.isdir(out) :
        log("WARNING: Output directory already exists: {}.".format(out))
    else :
        os.makedirs(out) # Create directory following path
    out = os.path.abspath(out)

    # 1. Index the reference
    indexed_ref = ref + ".bwt"
    cmd = "bwa index {ref}".format(**{"ref":ref})
    if not os.path.isfile(indexed_ref) :
        print(cmd + "\n")
        run(cmd)
    else :
        log("SKIP: found bwa index: {}".format(indexed_ref))

    # 2. Preparing outputs and required variables
    out_merge = os.path.join(out, dc_args["sample"] + ".merged.bam")
    out_merge_exists = os.path.isfile(out_merge)

    out_markdup = os.path.join(out, dc_args["sample"] + ".markdup.bam")
    out_markdup_exists = os.path.isfile(out_markdup)

    out_sort = os.path.join(out, dc_args["sample"] + ".sorted.CALL.bam")
    out_sort_exists = os.path.isfile(out_sort)

    intermediate_files = []

    # 3. Run alignment jobs if output files are not found
    if out_merge_exists or out_markdup_exists or out_sort_exists :
        log("SKIP: found output .bam file.")
        bam_to_merge = []
        for r1, r2 in zip(R1_reads, R2_reads) :
            out_sort_first = os.path.join(out, dc_args["sample"] + "_" + os.path.basename(r1) + "_" + os.path.basename(r1) + ".sorted.bam")
            bam_to_merge.append(out_sort_first)
    else :
        # Loop through given libraries
        i = 0
        bam_to_merge = []
        for r1, r2 in zip(R1_reads, R2_reads) :

            # 3.1 Get read group for alignment
            header, id = get_read_group(r1)
            readgroup = "@RG\\tID:{id}\\tSM:{sample}\\tLB:{id}_{sample}\\tPL:{sequencer}"
            dc_rg = {"id":id, "sample":dc_args["sample"], "sequencer":"ILLUMINA"}
            readgroup = readgroup.format(**dc_rg)

            # 3.2 Preparing outputs
            out_bwa_sam = os.path.join(out, dc_args["sample"] + "_" + os.path.basename(r1) + "_" + os.path.basename(r1) + ".sam")
            out_bwa_sam_exists = os.path.isfile(out_bwa_sam)

            out_view = os.path.join(out, dc_args["sample"] + "_" + os.path.basename(r1) + "_" + os.path.basename(r1) + ".view.bam")
            out_view_exists = os.path.isfile(out_view)

            out_sort_first = os.path.join(out, dc_args["sample"] + "_" + os.path.basename(r1) + "_" + os.path.basename(r1) + ".sorted.bam")
            out_sort_first_exists = os.path.isfile(out_sort_first)

            # 3.4 make a BWA mem command
            if out_bwa_sam_exists :
                log("SKIP: found .sam file.")
            else :
                bwa_align_reads(ref, r1, r2, out_bwa_sam, readgroup, dc_args["threads"], dc_args["output_logs"])
                intermediate_files.append(out_bwa_sam)

            # 3.5. make a sambamba view command: filtering on mapq
            if out_view_exists :
                log("SKIP: found filtered .bam file.")
            else :
                filter_mapq_bam(out_bwa_sam, out_view, dc_args["threads"], dc_args["mapq"], dc_args["output_logs"])
                intermediate_files.append(out_view)
                intermediate_files.append(out_view+".bai")

            # 3.6 sort the created alignment file
            if out_sort_first_exists :
                log("SKIP: found sorted .bam file.")
            else :
                sort_bam(out_view, out_sort_first, dc_args["threads"], dc_args["output_logs"])
                intermediate_files.append(out_sort_first)
                intermediate_files.append(out_sort_first+".bai")

            # 3.7 Add a merged file to merge into a final file
            if os.path.isfile(out_sort_first) : # cannot call out_sort_first_exists because it might have been created in the meantime
                bam_to_merge.append(out_sort_first)

    # 4. Merging all sorted-filtered-per library alignment files
    if out_markdup_exists or out_sort_exists :
        log("SKIP: found merged .bam file.")
    else :
        merge_bam_list(bam_to_merge, out_merge, dc_args["threads"], dc_args["output_logs"])
        intermediate_files.append(out_merge)
        intermediate_files.append(out_merge+".bai")

    # 5. Marking duplicates in the merged-filtered alignment file
    if out_sort_exists :
        log("SKIP: found sorted .bam file.")
    else :
        markdup_bam(out_merge, out_markdup, dc_args["threads"], dc_args["output_logs"])
        intermediate_files.append(out_markdup)
        intermediate_files.append(out_markdup+".bai")

    # 6. Sorting the markedup-merged-filtered final file
    sort_bam(out_markdup, out_sort, dc_args["threads"], dc_args["output_logs"])

    # 7. Indexing final file (may be useless, not sure because sambamba sort might already index)
    index_bam(out_sort, dc_args["threads"], dc_args["output_logs"])

    # IF REMOVE TEMP : remove sam file, filtered bam file, merged file, marked file
    # Keep only filtered, merged, marked and sorted final file
    if not keep :
        log("Removing:\n")
        for file in intermediate_files :
            if os.path.isfile(file) :
                print("- " + str(file))
                os.remove(file)


# Functions
def bwa_align_reads(ref, r1, r2, output, readgroup, threads, output_logs) :
    """Align a library with BWA-MEM"""
    if os.path.isfile(output) :
        log("SKIP: found alignement: {}.".format(output))
    else :
        dc_bwa = {"threads":threads, "readgroup":readgroup, "output":output, "ref":ref, "r1":r1, "r2":r2}
        cmd = "bwa mem -K 100000000 -t {threads} -R \"{readgroup}\" -o {output} {ref} {r1} {r2}"
        cmd = cmd.format(**dc_bwa)
        log("Aligning reads.")
        print(cmd + "\n")
        log_file = output + ".log" if output_logs else None
        run_log(cmd, log_file)

def filter_mapq_bam(bamin, bamout, threads, mapq, output_logs) :
    """Filter a bam file"""
    if os.path.isfile(bamout) :
        log("SKIP: found bam file: {}.".format(bamout))
    else :
        dc_view = {"threads":threads, "mapq":mapq, "output":bamout, "input_sam":bamin}
        cmd = "sambamba view -h -S -f bam -l 7 -t {threads} -F \"mapping_quality >= {mapq}\" -o {output} {input_sam}"
        cmd = cmd.format(**dc_view)
        log("Converting to BAM file.")
        print(cmd + "\n")
        #run(cmd)
        log_file = bamout + ".log" if output_logs else None
        run_log(cmd, log_file)

def merge_bam_list(bamlist, bamout, threads, output_logs) :
    """Merge a list of .bam files"""
    if os.path.isfile(bamout) :
        log("SKIP: found merged .bam file.")
    else :
        if len(bamlist) > 1 : # If more than one file is found
            # 4.1 create and run sambamba merge command
            input_bams = " ".join(x for x in bamlist)
            dc_merge = {"threads":threads, "output":bamout, "input_bams":input_bams}
            cmd = "sambamba merge -l 7 -t {threads} {output} {input_bams}"
            cmd = cmd.format(**dc_merge)
            log("Merging alignments.")
            print(cmd + "\n")
            log_file = bamout + ".log" if output_logs else None
            run_log(cmd, log_file)
            #run(cmd)
        else :
            log("SKIP: Merging is not necessary.")
            # In case only one sample change rename file
            os.rename(bamlist[0], bamout)

def sort_bam(bamin, bamout, threads, output_logs) :
    """Sort a .bam files"""
    if os.path.isfile(bamout) :
        log("SKIP: found sorted .bam file.")
    else :
        dc_sort = {"threads":threads, "output":bamout, "input":bamin}
        cmd = "sambamba sort -l 7 -t {threads} -o {output} {input}"
        cmd = cmd.format(**dc_sort)
        log("Sorting bam file.")
        print(cmd + "\n")
        log_file = bamout + ".log" if output_logs else None
        run_log(cmd, log_file)
        #run(cmd)

def markdup_bam(bamin, bamout, threads, output_logs) :
    """Mark duplicates in a .bam file"""
    if os.path.isfile(bamout) :
        log("SKIP: found marked duplicates .bam file.")
    else :
        dc_markdup = {"threads":threads, "output":bamout, "input":bamin}
        cmd = "sambamba markdup --overflow-list-size 600000 --hash-table-size 500000 -l 7 -t {threads} {input} {output}"
        cmd = cmd.format(**dc_markdup)
        log("Marking duplicates.")
        print(cmd + "\n")
        log_file = bamout + ".log" if output_logs else None
        run_log(cmd, log_file)
        #run(cmd)
