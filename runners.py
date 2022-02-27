#!/usr/bin/python
# -*- coding: utf-8 -*-

import subprocess
from utils import log
#from multiprocessing import Pool, TimeoutError

def run(cmd, shell=True, OUT=subprocess.PIPE, ERR=subprocess.PIPE) :
    if ERR == "devnull" :
        ERR = subprocess.DEVNULL
    proc = subprocess.Popen(cmd, shell=shell, stdout=OUT, stderr=ERR)
    proc.communicate()

def run_log(cmd, log_file) :
    if not log_file : # log_file == None
        run(cmd, ERR=subprocess.DEVNULL)
    else :
        lf = open(log_file, "w")
        lf.write(cmd + "\n")
        run(cmd, OUT=lf, ERR=lf)
        lf.close()

def run_FP(job) :
    cmd = job[0]
    number = job[1]
    log("Running FastP job number {}".format(number), newline=False)
    run(cmd)
    log("Finished FastP job number {}".format(number), newline=False)

def run_picard_(job) :
    cmd = job[0]
    log_file = job[1]
    number = job[2]
    log("Running GenotypeGVCFs job number {}".format(number), newline=False)
    if not log_file :
        run(cmd, ERR=subprocess.DEVNULL)
    else :
        lf = open(log_file, "w")
        lf.write(cmd + "\n")
        run(cmd, ERR=lf)
        lf.close()
    log("Finished GenotypeGVCFs job number {}".format(number), newline=False)

def run_GG(job) :
    cmd = job[0]
    log_file = job[1]
    number = job[2]
    log("Running GenotypeGVCFs job number {}".format(number), newline=False)
    if not log_file :
        run(cmd, ERR=subprocess.DEVNULL)
    else :
        lf = open(log_file, "w")
        lf.write(cmd + "\n")
        run(cmd, ERR=lf)
        lf.close()
    log("Finished GenotypeGVCFs job number {}".format(number), newline=False)

def run_HC(job) :
    cmd = job[0]
    log_file = job[1]
    number = job[2]
    log("Running HaplotypeCaller job number {}".format(number), newline=False)
    if not log_file : # if do not log do not keep stderr
        run(cmd, ERR=subprocess.DEVNULL)
    else :
        lf = open(log_file, "w")
        lf.write(cmd + "\n")
        run(cmd, ERR=lf)  # if output logs -> go to log_file
        lf.close()
    log("Finished HaplotypeCaller job number {}".format(number), newline=False)
