#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import yaml # REQUIRE pyyaml <- not in python package

from utils import log, check_files

#     _____                         _                      _                  _   _
#    / ____|                       | |                    (_)                | | (_)
#   | |       _ __    ___    __ _  | |_    ___     _ __    _   _ __     ___  | |  _   _ __     ___
#   | |      | '__|  / _ \  / _` | | __|  / _ \   | '_ \  | | | '_ \   / _ \ | | | | | '_ \   / _ \
#   | |____  | |    |  __/ | (_| | | |_  |  __/   | |_) | | | | |_) | |  __/ | | | | | | | | |  __/
#    \_____| |_|     \___|  \__,_|  \__|  \___|   | .__/  |_| | .__/   \___| |_| |_| |_| |_|  \___|
#                                                 | |         | |
#                                                 |_|         |_|

def create_pipeline(args) :
    """Pipeline to prepare reference genome files for variant calling"""

    # Parse arguments & check files
    config_file = args.config[0]
    config_file = check_files([config_file])[0]

    print("# runGATK.py pipeline")
    print("Configuration file:\t{}\n".format(config_file))
    print("===============================================================================\n")

    # Read config file and create a .bash file
    f = open(config_file, "r")
    cfg = yaml.load(f, Loader=yaml.FullLoader)
    f.close()

    # Check output directory
    if not os.path.isdir(cfg["output_directory"]) :
        os.makedirs(cfg["output_directory"])
    else :
        raise Exception("ERROR: Output directory already exists!")

    # Check runGATK, reference and input reads
    runGATK = check_files([cfg["rungatk"]])[0]
    ref = check_files([cfg["reference"]])[0]

    for sample, libraries in cfg["samples"].items() :
        for lib, reads in libraries.items() :
            reads["R1"] = check_files([reads["R1"]])[0]
            reads["R2"] = check_files([reads["R2"]])[0]

    # Write pipeline
    f = open(os.path.join(cfg["output_directory"], "pipeline.sh"), 'w')
    f.write("#!/bin/bash\n\n# Global variables\nrungatk={}\nref={}\n\n".format(runGATK, ref))

    # Step preparing reference
    f.write("\n## Preparing reference\n")
    f.write("echo 'Preparing reference files:'\n")
    f.write("python ${{rungatk}} prepare -ws {} ${{ref}}\n".format(cfg["prepare_options"]["window_size"]))

    # Step Trimming
    if cfg["trim_options"]["skip_trimming"] : # Skip trimmming
        pass
    else : # Do not skip trimmming
        f.write("\n## Trimming reads\n")
        f.write("echo 'Trimming:'\n")

        # Create a reads file
        readslist = os.path.join(os.path.abspath(cfg["output_directory"]), "sample_reads.list")
        fr = open(readslist, "w")
        for sample, libraries in cfg["samples"].items() : # for each library
            for lib, reads in libraries.items() : # for reads in library
                fr.write("{}\t{}\n".format(reads["R1"], reads["R2"]))
        fr.close()

        dc = {
              "threads":cfg["trim_options"]["threads_per_job"],
              "jobs":cfg["trim_options"]["njobs"],
              "options":cfg["trim_options"]["fastp_options"],
              "reads":readslist
              }
        f.write("python ${{rungatk}} trim -t {threads} -p {jobs} -fo=\"{options}\" {reads}\n".format(**dc))

    # Step aligning
    f.write("\n## Aligning trimmed reads\n")
    f.write("echo 'Aligning:'\n")

    # Per sample alignment
    trimreads = {sample:{"R1":[], "R2":[]} for sample in cfg["samples"].keys()}
    # Store libraries with trimled.gz suffix
    for sample, libraries in cfg["samples"].items() :
        for lib, reads in libraries.items() :
            trimreads[sample]["R1"].append(reads["R1"] + ".trimmed.gz")
            trimreads[sample]["R2"].append(reads["R2"] + ".trimmed.gz")

    for sample in cfg["samples"].keys() :
        f.write("# Aligning {}\n".format(sample))

        dc = {
              "sample":sample,
              "threads":cfg["align_options"]["threads"],
              "ram":cfg["align_options"]["ram_per_thread"],
              "mapq":cfg["align_options"]["mapq"],
              "R1":",".join(r for r in trimreads[sample]["R1"]),
              "R2":",".join(r for r in trimreads[sample]["R2"]),
              }

        cmd = "python ${{rungatk}} align -t {threads} -r {ram} -m {mapq} {sample} ${{ref}} {R1} {R2}".format(**dc)
        if cfg["align_options"]["keep_temp"] == "true" :
            cmd += " --keep-temp\n"
        else :
            cmd += "\n"
        f.write(cmd)

    # Step calling samples
    f.write("\n## Calling samples\n")
    f.write("echo 'Calling:'\n")

    for sample in cfg["samples"].keys() :
        f.write("# Calling {}\n".format(sample))
        fi = cfg["call_options"]["fi"] if cfg["call_options"]["fi"] != '' else "\'\'"
        dc = {
              "sample":sample,
              "njobs":cfg["call_options"]["njobs"],
              "ht":cfg["call_options"]["ht"],
              "pim":cfg["call_options"]["pim"],
              "he":cfg["call_options"]["he"],
              "ihe":cfg["call_options"]["ihe"],
              "mrpas":cfg["call_options"]["mrpas"],
              "fi":fi,
              "erc":cfg["call_options"]["erc"],
              "om":cfg["call_options"]["om"],
              "jo":cfg["call_options"]["java_options"],
              "pjo":cfg["call_options"]["picard_java_options"],
              "mmn":cfg["call_options"]["mmn"],
        }
        cmd = "python ${{rungatk}} call {sample} ${{ref}} -pjo=\"{pjo}\" -jo=\"{jo}\" "
        cmd += "-p {njobs} -ht {ht} -he {he} -ihe {ihe} -mrpas {mrpas} -erc {erc} "
        cmd += "-om {om} -fi {fi} -mmn {mmn}"
        if cfg["call_options"]["bam_output"] == "true" :
            cmd += " --bam-out"
        if cfg["call_options"]["log"] == "true" :
            cmd += " --log"
        if cfg["call_options"]["keep_temp"] == "true" :
            cmd += " --keep-temp"
        cmd = cmd.format(**dc) + "\n"
        f.write(cmd)

    # Genotyping
    f.write("\n## Joint-genotyping\n")
    f.write("echo 'Genotyping:'\n")

    fi = cfg["genotype_options"]["fi"] if cfg["genotype_options"]["fi"] != '' else "\'\'"
    dc = {
          "sample":sample,
          "njobs":cfg["genotype_options"]["njobs"],
          "sp":cfg["genotype_options"]["sp"],
          "fi":fi,
          "he":cfg["genotype_options"]["he"],
          "ihe":cfg["genotype_options"]["ihe"],
          "jo":cfg["genotype_options"]["java_options"],
          "pjo":cfg["genotype_options"]["picard_java_options"],
          "dbjo":cfg["genotype_options"]["db_java_options"],
          "mmn":cfg["genotype_options"]["mmn"],
          "bs":cfg["genotype_options"]["bs"],
          "samples":",".join(sample for sample in cfg["samples"].keys())
    }

    cmd = "python ${{rungatk}} genotype ${{ref}} {samples} jointgenotyping "
    cmd += "-jo=\"{jo}\" -pjo=\"{pjo}\" -dbjo=\"{dbjo}\" -p {njobs} "
    cmd += "-he {he} -ihe {ihe} -fi {fi} -mmn {mmn} -bs {bs}"

    if cfg["genotype_options"]["all_sites"] == "true" :
        cmd += " --all-sites"
    if cfg["genotype_options"]["keep_temp"] == "true" :
        cmd += " --keep-temp"
    if cfg["genotype_options"]["log"] == "true" :
        cmd += " --log"

    cmd = cmd.format(**dc)
    cmd += "\n"
    f.write(cmd)

    f.write("\necho \"Done!\"\n")
    f.close()

    print("Done!")
    print("A file pipeline.sh was created in the ./output directory.")
    print("You can make it executable with `chmod +x pipeline.sh`")
    print("A final joint-genotyped .vcf file will be stored in \"output/jointegenotyping\".")
