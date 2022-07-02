#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from time import localtime, strftime

def log(string, newline=True) :
    if newline :
        print("\n{}: {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), string))
    else :
        print("{}: {}".format(strftime("%Y-%m-%d %H:%M:%S", localtime()), string))

def get_read_group(path) :
    """Obtains the read group from fq header"""
    try :
        f = open(path, "r")
        head = str(f.readline().strip()) # read first line
        f.close()
        id = "_".join(x for x in head.split(':')[0:4]).replace("@", '')
        return head, id
    except :
        import gzip
        f = gzip.GzipFile(path, 'rb')
        head = str(f.readline().decode('UTF-8').strip()) # read first line
        f.close()
        id = "_".join(x for x in head.split(':')[0:4]).replace("@", '')
        return head, id


def check_files(files) :
    """Returns absolute file paths and raise exception if file does not exist"""
    absfiles = []
    for file in files :
        if not os.path.isfile(file) :
            raise Exception("ERROR: {} is not found!".format(file))
        else :
            absfiles.append(os.path.abspath(file))
    return absfiles

def check_dirs(dirs) :
    """Returns absolute paths and raise exception if dir does not exist"""
    absdirs = []
    for d in dirs :
        if not os.path.isdir(d) :
            raise Exception("ERROR: {} is not found!".format(d))
        else :
            absdirs.append(os.path.abspath(d))
    return absdirs

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

def parse_fasta(file) :
    """parse a fasta file"""
    seq = {}
    f = open(file, "r")
    for line in f :
        if line[0] == ">" :
            name = line[1:].strip().split()[0] # remove the end at first space
            seq[name] = ""
        else :
            seq[name] += line.strip()
    f.close()
    return seq

def fadict(file, output) :
    """substitute for picard CreateSequenceDictionary"""

    from hashlib import md5

    if os.path.isfile(file) :
        fpath = os.path.abspath(file)
    else :
        raise Exception("ERROR: file you are trying to index does not exist")

    seqs = parse_fasta(fpath)
    odc = open(output, "w")
    odc.write("@HD\tVN:1.6\n")
    for sid, seq in seqs.items() :
        string = "@SQ\tSN:{0}\tLN:{1}\tM5:{2}\tUR:file:{3}\n"
        md5sum = md5(seq.encode('utf-8')).hexdigest()
        #print(sid, len(seq), md5sum, fpath)
        string = string.format(sid, len(seq), md5sum, fpath)
        odc.write(string)

    odc.close()

def faidx(file) :
    """Fasta file must be uncompressed and have perfect format because no error checks"""

    if os.path.isfile(file) :
        fpath = os.path.abspath(file)
    else :
        raise Exception("ERROR: file you are trying to index does not exist")

    fid = fpath + ".fai"

    fst = open(fpath, "rb")
    idx = open(fid, "w")

    # From pyfaidx
    ref_name = None  # reference sequence name
    offset = 0  # binary offset of end of current line
    p_offset = offset
    ref_len = 0  # reference character length
    bin_len = 0  # binary line length (includes newline)
    char_len = 0  # character line length
    bad_lines = []  # lines > || < than blen
    valid_entry = False
    last_line = None

    for i, line in enumerate(fst) :
        line_bin_len = len(line)
        line = line.decode()
        line_char_len = len(line.rstrip('\n\r'))
        last_line = i
        if line[0] == '>' :
            if i > 0 :
                string = "{0}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\n"
                string = string.format(ref_name, ref_len, p_offset, char_len, bin_len)
                idx.write(string)
            bin_len = 0
            ref_len = 0
            char_len = 0
            ref_name = line.rstrip('\n\r')[1:].split()[0]
            offset += line_bin_len
            p_offset = offset
        else :
            if not bin_len:
                bin_len = line_bin_len
            if not char_len:
                char_len = line_char_len
            offset += line_bin_len
            ref_len += line_char_len

    # Write last line
    if last_line is not None :
        string = "{0}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\n"
        string = string.format(ref_name, ref_len, p_offset, char_len, bin_len)
        idx.write(string)

    fst.close()
    idx.close()

def read_dict(dict_file) :
    chromosomes = {}
    f = open(dict_file, 'r')
    for line in f :
        if "@SQ" in line :
            s = line.strip().split("\t")
            name, length = None, None
            for d in s :
                if d.startswith("SN:") :
                    name = d[3:]
                elif d.startswith("LN:") :
                    length = int(d[3:])
                else :
                    continue
            if name is not None :
                chromosomes[name] = length
            else :
                pass
        else :
            continue
    return chromosomes

def index_bam(bam, threads, output_logs) :
    """Index a .bam file"""

    import subprocess
    def run(cmd, shell=True, OUT=subprocess.PIPE, ERR=subprocess.PIPE) :
        if ERR == "devnull" :
            ERR = subprocess.DEVNULL
        proc = subprocess.Popen(cmd, shell=shell, stdout=OUT, stderr=ERR)
        proc.communicate()

    bamidx = bam + ".bai"
    if not os.path.isfile(bamidx) :
        dc_index = {"input":bam, "threads":threads}
        cmd = "sambamba index -t {threads} {input}"
        cmd = cmd.format(**dc_index)
        print("\n{}: Indexing final bam file".format(strftime("%Y-%m-%d %H:%M:%S", localtime())))
        print(cmd + "\n")
        log_file = bamidx + ".log" if output_logs else None
        if not log_file :
            run(cmd, ERR="devnull")
        else :
            lf = open(log_file, "w")
            lf.write(cmd + "\n")
            run(cmd, OUT=lf, ERR=lf)
            lf.close()
    else :
        print("SKIP: bam index file found: {}".format(bamidx))

def remove_temp_files(interm_files) :
    """Clean up intermediate files"""
    log("Removing temporary files:")
    for file in interm_files :
        # Remove file
        if os.path.isfile(file) :
            os.remove(file)
        # Remove index file
        if os.path.isfile(file + ".idx") :
            os.remove(file + ".idx")

def check_args(dc_args) :
    """Check that arguments are valid"""
    # MMN cannot be 1
    if dc_args["chunk_size"] < 2 :
        log("WARNING: Max merge number must be >= 2: setting to 2")
        dc_args["chunk_size"] = 2

    """ # Probably useless
    for s in ["java", "pjo"] :
        # Check "" in the java options and remove them if necessary
        if dc_args[s][0] == '"' :
            dc_args[s] = dc_args[s][1:]
        if dc_args[s][-1] == '"' :
            dc_args[s] = dc_args[s][:-1]
    """

    return dc_args
