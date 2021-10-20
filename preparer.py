#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os
from time import localtime, strftime
from utils import parse_fasta, check_files, faidx, fadict

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
        seqs = parse_fasta(ref)
        for sid, record in seqs.items() :
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
        seqs = parse_fasta(ref)
        for sid, record in seqs.items() :
            f.write("{}\t0\t{}\n".format(sid, len(record)))
        f.close()
    else :
        print("SKIP: Chromosomes .bed file found: {}".format(out_chromosomes))

    # 2. Create sequence dictionnary with picard tools
    out_dict = os.path.splitext(ref)[0] + ".dict"
    if not os.path.isfile(out_dict) :
        fadict(ref, out_dict)
    else :
        print("SKIP: Reference dictionary file found: {}".format(out_dict))

    # 3. Indexing fasta file
    out_index = ref + ".fai"
    if not os.path.isfile(out_index) :
        print("\n{}: Indexing reference".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        faidx(ref)
    else :
        print("SKIP: Reference index found: {}".format(out_index))
