#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Small GATK alignment and variant calling pipeline using python"""

#import errno
#import os
#import shutil
#import datetime

import argparse
import sys

#import subprocess
#from multiprocessing import Pool, TimeoutError

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

def pipeline(args) :
    """Runs the pipeline module"""
    from pipeliner import create_pipeline
    create_pipeline(args)

def prepare(args) :
    """Runs the prepare module"""
    from preparer import prepare_genome
    prepare_genome(args)

def align(args) :
    """Runs the align module"""
    from aligner import align_reads
    align_reads(args)

def trim(args) :
    """Runs the trim module"""
    from trimmer import trim_reads
    trim_reads(args)

def call(args) :
    """Runs the call module"""
    from caller import bam_call
    bam_call(args)

def genotype(args) :
    """Runs the genotype module"""
    from genotyper import genotype_samples
    genotype_samples(args)

def main() :
    """Argument parser"""
    parser = argparse.ArgumentParser(description='Runs modules of a variant calling pipeline with bwa-mem and GATK.')
    subparsers = parser.add_subparsers(required=True, dest="align || prepare || call || genotype || vqsr || allelecount")

    # Trim Illumina short reads with fastp
    pip = subparsers.add_parser('pipeline',       help="Create a pipeline bash file to run variant calling")
    pip.add_argument('config', nargs=1, type=str, help="<STRING> A YAML file containing the samples, reads and options.")
    pip.set_defaults(func=pipeline)

    # Align reads to the reference
    aln = subparsers.add_parser('align', help="Aligns reads from ONE sample on a reference genome\nReads are considered Illumina Paired-End\nReads list correspond to lanes and must be given in the same order or the alignment will fail")
    aln.add_argument('Sample',        nargs=1, type=str, default=['sample1'],help="<STRING> sample name for output and adding to read group. Default: 'sample1'")
    aln.add_argument('Reference',     nargs=1, type=str,               help="<STRING> A path to the reference fenome fasta file.")
    aln.add_argument('R1',            nargs=1, type=list_str,          help="<STRING> A comma separated list of FastQ files FORWARD reads sequences separated by LANES (can be compressed).")
    aln.add_argument('R2',            nargs=1, type=list_str,          help="<STRING> A comma separated list of FastQ files REVERSE reads sequences separated by LANES (can be compressed).")
    aln.add_argument('-t','--threads',nargs=1, type=int, default=[4],  required=False,help="<INT> Maximum threads to use. Default: %(default)s")
    aln.add_argument('-r','--ram',    nargs=1, type=int, default=[4],  required=False,help="<INT> Number of GB of RAM per thread to use. Default: %(default)s")
    aln.add_argument('-m','--mapq',   nargs=1, type=int, default=[20], required=False,help="<INT> Sambamba view filtering to keep only alignements with MAPQ >= <INT>. Default: %(default)s")
    aln.add_argument('-l','--log',             action="store_true",    help="Output log files of aligner steps, useful for debugging. Default: %(default)s")
    aln.add_argument('-kp', '--keep-temp',     action="store_true",    help="Do not remove intermediate steps alignments (cannot be before positional argument). Default: %(default)s")
    aln.set_defaults(func=align)

    # Prepare genome assembly for the HaplotypeCaller
    # WARNING WINDOW SIZE IN KB
    prep = subparsers.add_parser('prepare',     help="Create several files from a reference to prepare for the variant calling")
    prep.add_argument('Reference',              nargs=1, type=str, help="<STRING> A fasta file containing the reference genome.")
    prep.add_argument('-ws','--window-size',    nargs=1, type=int, default=[100], required=False, help="<INT> Window size (in kb!!!) to split genome. Default: %(default)s")
    prep.set_defaults(func=prepare)

    # Trim Illumina short reads with fastp
    trm = subparsers.add_parser('trim',       help="Trim input sequencing reads with fastp")
    trm.add_argument('FOFN',                  nargs=1, type=str, help="<STRING> A file of filenames containing paths to the reads to trim. Format: R1\\tR2.")
    trm.add_argument('-t','--threads',        nargs=1, type=int, default=[4], required=False, help="<INT> Maximum threads for fastp process. Default: %(default)s")
    trm.add_argument('-p','--processes',      nargs=1, type=int, default=[4], required=False, help="<INT> Maximum parallel fastp processes in pool. Each process uses the number of threads defined by -t. Default: %(default)s")
    trm.add_argument('-fo','--fastp-options', nargs=1, type=str, default=['--detect_adapter_for_pe --trim_poly_g --compression 9'], help="<STRING> fastp filtering options. Default: %(default)s")
    trm.add_argument('-nr','--no-report',     action="store_true",    help="Redirect fastp reports to /dev/null. Default: %(default)s")
    trm.set_defaults(func=trim)

    # Parallel HaplotypeCaller
    cal = subparsers.add_parser('call', help="Runs GATK HaplotypeCaller with a sample and the prepared reference genome assembly.")
    cal.add_argument('Sample',                          nargs=1, type=str, help="<STRING> sample name for output and SAME AS \"runGATK.py align\"")
    cal.add_argument('Reference',                       nargs=1, type=str, help="<STRING> A path to the reference genome fasta file")
    cal.add_argument('-p','--processes',                nargs=1, type=int,   default=[4],        required=False, help="<INT> Maximum processes in pool to use. Default: %(default)s")
    cal.add_argument('-ht','--hmm-threads',             nargs=1, type=int,   default=[4],        required=False, help="<INT> Number of threads to use for HMM. Default: %(default)s (per process).")
    cal.add_argument('-he','--heterozygosity',          nargs=1, type=float, default=[0.01],     required=False, help="<FLOAT> Heterozygosity value to pass to HaplotypeCaller. Default: %(default)s")
    cal.add_argument('-mrpas','--max-reads-per-alignment-start', nargs=1,type=int,default=[50],  required=False, help="<INT> See GATK doc. Default: %(default)s (set to 0 to deactivate).")
    cal.add_argument('-ihe','--indel-heterozygosity',   nargs=1, type=float, default=[0.0001],   required=False, help="<FLOAT> Indel heterozygosity value to pass to HaplotypeCaller. Default: %(default)s")
    cal.add_argument('-fi', '--founder-id',             nargs=1, type=str,   default=[''],                       help="<STRING> Founder population sample ID. Default: empty")
    cal.add_argument('-pim', '--pcr-indel-model',       nargs=1, type=str,   default=['NONE'],                   choices=["NONE","CONSERVATIVE","HOSTILE","AGGRESSIVE"], help="<STRING> Argument to pass to HaplotypeCaller --pcr-indel-model. Default: %(default)s")
    cal.add_argument('-erc', '--emit-ref-confidence',   nargs=1, type=str,   default=['BP_RESOLUTION'],          choices=['BP_RESOLUTION','GVCF'], help="<STRING> ERC (see HaplotypeCaller documentation). Default: %(default)s")
    cal.add_argument('-om', '--output-mode',            nargs=1, type=str,   default=['EMIT_ALL_ACTIVE_SITES'],  choices=['EMIT_VARIANTS_ONLY','EMIT_ALL_CONFIDENT_SITES','EMIT_ALL_ACTIVE_SITES'], help="<STRING> Output mode (see HaplotypeCaller documentation). Default: %(default)s")
    cal.add_argument('-jo', '--java-options',           nargs=1, type=str,   default=['-Xmx4G'],                 help="<STRING> Java Virtual Machine options (Ram Per Process is defined here). Default: %(default)s")
    cal.add_argument('-pjo', '--picard-java-options',   nargs=1, type=str,   default=['-Xmx4G'],                 help="<STRING> Java Virtual Machine options. Default: %(default)s")
    cal.add_argument('-mmn','--max-merge-number',       nargs=1, type=int,   default=[500],      required=False, help="<INT> Maximum number of files to merge with picard MergeVcfs (Higher number, higher RAM usage). Default: %(default)s")
    cal.add_argument('-kp', '--keep-temp',  action="store_true", help="Do not remove intermediate steps files (cannot be before positional argument). Default: %(default)s")
    cal.add_argument('-dr', '--dry-run',    action="store_true", help="Only parse arguments for testing and debugging commands. Default: %(default)s")
    cal.add_argument('-bo','--bam-out',     action="store_true", help="Output realigned reads bam. Default: %(default)s")
    cal.add_argument('-l','--log',          action="store_true", help="Output sub files log of HaplotypeCaller, useful for debugging. Default: %(default)s")
    cal.set_defaults(func=call)

    # GenotypeGVCFS
    gen = subparsers.add_parser('genotype',             help="Runs GATK joint calling using outputs of different samples calls and the prepared reference output directory.")
    gen.add_argument('Reference',                       nargs=1,             type=str,           help="<STRING> A path to the output directory of \"runGATK.py prepare\"")
    gen.add_argument('Samples',                         nargs=1,             type=str,           help="<STRING> A comma separated list of samples (directory must contain call/merged.g.vcf)")
    gen.add_argument('Output',                          nargs=1,             type=str,           help="<STRING> An output directory path for the joint calling.")
    gen.add_argument('-p','--processes',                nargs=1, type=int,   default=[4],        required=False, help="<INT> Maximum threads to use. Default: %(default)s")
    gen.add_argument('-sp','--sample-processes',        nargs=1, type=int,   default=[4],        required=False, help="<INT> Maximum threads to read samples simultaneously (GenomicsDBImport only). Default: %(default)s")
    gen.add_argument('-bs','--batch-size',              nargs=1, type=int,   default=[50],       required=False, help="<INT> Default batch size for GenomicsDBImport. Default: %(default)s")
    gen.add_argument('-mmn','--max-merge-number',       nargs=1, type=int,   default=[50],       required=False, help="<INT> Maximum number of files to merge with picard MergeVcfs (Higher number, higher RAM usage). Default: %(default)s")
    gen.add_argument('-jo', '--java-options',           nargs=1, type=str,   default=['-Xmx4G'], help="<STRING> Java Virtual Machine options (Ram Per Process is defined here). Default: %(default)s")
    gen.add_argument('-dbjo', '--db-java-options',      nargs=1, type=str,   default=['-Xmx56G'],help="<STRING> Java Virtual Machine options for GenomicsDBImport (Ram Per Process is defined here). Default: %(default)s")
    gen.add_argument('-pjo', '--picard-java-options',   nargs=1, type=str,   default=['-Xmx4G'], help="<STRING> Java Virtual Machine options for picard commands. Default: %(default)s")
    gen.add_argument('-fi', '--founder-id',             nargs=1, type=str,   default=[''],       help="<STRING> Founder population sample ID. Default: \"\" (empty)")
    gen.add_argument('-he','--heterozygosity',          nargs=1, type=float, default=[0.01],     required=False, help="<FLOAT> Heterozygosity value to pass to HaplotypeCaller. Default: %(default)s")
    gen.add_argument('-ihe','--indel-heterozygosity',   nargs=1, type=float, default=[0.0001],   required=False, help="<FLOAT> Indel heterozygosity value to pass to HaplotypeCaller. Default: %(default)s")
    gen.add_argument('-a', '--all-sites',               action="store_true", help="Output all sites in VCFs (in addition to normal vcf output).")
    gen.add_argument('-kp', '--keep-temp',              action="store_true", help="Do not remove intermediate steps files (cannot be before positional argument).")
    gen.add_argument('-l','--log',                      action="store_true", help="Output sub files log of HaplotypeCaller, useful for debugging.")
    gen.set_defaults(func=genotype)

    args = parser.parse_args()
    args.func(args)
    sys.exit(0)

if __name__ == '__main__':
	main()
