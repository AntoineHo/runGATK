#!/usr/bin/python
# -*- coding: utf-8 -*-

from utils import log, check_files

import os
from multiprocessing import Pool
from runners import run_FP

# Pipeline
def trim_reads(args) :
    """Pipeline to prepare reference genome files for variant calling"""

    # 0. parse arguments & check files
    freads = check_files(args.FOFN)[0]

    dc_args = {"threads":args.threads[0], "nproc":args.processes[0],
               "fo":args.fastp_options[0], "no_report":args.no_report,
              }

    print("# runGATK.py trim")
    print("Reads filepath list:\t{}\n".format(freads))
    print("Other arguments:\t" + str(dc_args))
    print("===============================================================================\n")

    # Read file and create jobs
    jobs = []
    number = 0 # job number
    f = open(freads, "r")
    for line in f :
        if len(line.strip()) == 0 :
            continue
        elif line[0] == "#" :
            continue
        dc = {k:v for k,v in dc_args.items()}
        pairs = line.strip().split("\t")
        if len(pairs) != 2 :
            raise Exception("Reads file appear unpaired!\nPlease use format: R1\tR2.\nExample: /path/to/R1.fq(.gz)\t/path/to/R2.fq(.gz)")

        dc["R1"] = check_files([pairs[0]])[0]
        dc["R2"] = check_files([pairs[1]])[0]
        dc["O1"] = dc["R1"] + ".trimmed.gz"
        dc["O2"] = dc["R2"] + ".trimmed.gz"

        if dc_args["no_report"] : # skip fastp report files
            dc["json"] = "/dev/null"
            dc["html"] = "/dev/null"
        else :
            dc["json"] = dc["O1"] + ".fastp_report_R1_R2.json" # default
            dc["html"] = dc["O1"] + ".fastp_report_R1_R2.html"

        if os.path.isfile(dc["O1"]) and os.path.isfile(dc["O2"]) :
            log("SKIP: found trimmed file: {}.".format(dc["O1"]))
            continue
        else :
            cmd = "fastp --json {json} --html {html} -w {threads} -i {R1} -I {R2} -o {O1} -O {O2} {fo}"
            cmd = cmd.format(**dc)
            jobs.append((cmd, number))
            number += 1
    f.close()

    # Run jobs in parallel
    p = Pool(dc_args["nproc"])
    p.map(run_FP, jobs)
    p.close() # Required so that pool stops correctly and program does not hang
    p.terminate()
